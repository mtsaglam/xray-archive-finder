"""Download HEASARC tables and write them as static sky-cell JSON files for the browser UI.

The browser cannot query HEASARC directly (no CORS), so this runs in GitHub Actions and the
static UI reads only the cells around the requested position.

  --observations  mission master tables (changes daily)  -> OUT/data/<mission>/<cell>.json
  --catalogs      X-ray / radio source catalogs (static)  -> OUT/data/cat/<catalog>/<cell>.json

Usage:  python build_index.py OUT_DIR [--observations] [--catalogs]   (default: both)
"""
import json
import math
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone

import archive

CELL = 2.0  # degrees
N_DEC = int(180 / CELL)
N_RA = int(360 / CELL)
MAX_ROWS = 90000  # HEASARC TAP truncates at 100000 rows

# flux errors: "err" symmetric, "eneg"/"epos" asymmetric, or "lo"/"hi" confidence limits.
# match: cross-match radius (arcsec); params: extra columns shown next to the flux.
CATALOGS = [
    dict(key="csc", label="Chandra Source Catalog 2.2", kind="x", table="csc", band="0.5–7 keV",
         unit="erg/s/cm²", match=5, flux="b_flux_ap", lo="b_flux_ap_lo", hi="b_flux_ap_hi",
         poserr="error_ellipse_r0", poserr_note="95%",
         params=[("HR h-m", "hardness_ratio_hm"), ("HR m-s", "hardness_ratio_ms"),
                 ("var", "var_flag"), ("ext", "extent_flag")]),
    dict(key="xmm", label="5XMM-DR15 (XMM-Newton)", kind="x", table="xmmssc", band="0.2–12 keV",
         unit="erg/s/cm²", match=10, flux="ep_flux", err="ep_flux_error", poserr="error_radius",
         poserr_note="1σ",
         params=[("HR2", "ep_hr2"), ("HR3", "ep_hr3"), ("var", "var_flag"), ("ext″", "extent")]),
    dict(key="sxps", label="2SXPS (Swift-XRT)", kind="x", table="swift2sxps", band="0.3–10 keV",
         unit="erg/s/cm²", match=15, flux="pow_flux_a", eneg="pow_flux_a_neg_err",
         epos="pow_flux_a_pos_err", poserr="error_radius", poserr_note="90%",
         params=[("rate c/s", "count_rate_fb"), ("HR1", "hardness_ratio_1"), ("HR2", "hardness_ratio_2")]),
    dict(key="erass1", label="eRASS1 (eROSITA)", kind="x", table="erass1main", band="0.2–2.3 keV",
         unit="erg/s/cm²", match=20, flux="b1_flux", err="b1_flux_error", poserr="error_radius",
         poserr_note="1σ",
         params=[("rate c/s", "b1_count_rate"), ("ext″", "source_extent"), ("DET_LIKE", "b0_detect_likelihood")]),
    dict(key="rxs2", label="2RXS (ROSAT)", kind="x", table="rass2rxs", band="0.1–2.4 keV",
         unit="erg/s/cm²", match=40, flux="plaw_flux", rate_err=("count_rate", "count_rate_error"),
         params=[("rate c/s", "count_rate"), ("HR1", "hardness_ratio_1"), ("HR2", "hardness_ratio_2")]),
    dict(key="nustar", label="NuSTAR Serendipitous Survey", kind="x", table="nustarssc", band="3–24 keV",
         unit="erg/s/cm²", match=30, flux="fb_flux", err="fb_flux_error",
         params=[("8–24 keV", "hb_flux")]),
    dict(key="bat", label="Swift-BAT 157-month", kind="x", table="swbat157m", band="14–195 keV",
         unit="erg/s/cm²", match=420, flux="flux", lo="flux_lower", hi="flux_upper",
         params=[("S/N", "snr"), ("counterpart", "ctrpart_name"), ("type", "source_type")]),
    dict(key="ibis", label="INTEGRAL/IBIS (5th cat.)", kind="x", table="ibiscat", band="20–40 keV",
         unit="mCrab", match=300, flux="sb_flux", err="sb_flux_error", poserr="error_radius",
         params=[("40–100 keV", "hb_flux"), ("type", "source_type"), ("var", "variability_flag")]),
    dict(key="nvss", label="NVSS", kind="radio", table="nvss", band="1.4 GHz", unit="mJy", match=30,
         flux="flux_20_cm", err="flux_20_cm_error", noname=True,
         params=[("maj″", "major_axis"), ("min″", "minor_axis")]),
    dict(key="sumss", label="SUMSS", kind="radio", table="sumss", band="843 MHz", unit="mJy", match=30,
         flux="int_flux_36_cm", err="int_flux_36_cm_error", noname=True, params=[]),
    dict(key="wenss", label="WENSS", kind="radio", table="wenss", band="325 MHz", unit="mJy", match=45,
         flux="int_flux_92_cm", err="flux_92_cm_error", noname=True, params=[]),
]


def cell_index(ra, dec):
    i = min(N_DEC - 1, max(0, int(math.floor((dec + 90.0) / CELL))))
    j = int(math.floor((ra % 360.0) / CELL)) % N_RA
    return i * N_RA + j


def rnd(v, nd):
    return None if v is None else round(float(v), nd)


def sig(v, n=4):
    """Round numbers to n significant digits; pass other values through."""
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        return v
    if v == 0 or v != v:
        return 0 if v == 0 else None
    return float("{:.{}g}".format(v, n))


def fetch_rows(key, table, cols, lo, hi):
    """Rows with lo <= dec < hi, splitting the band until each query is below the row cap."""
    where = "dec >= {} AND dec < {}".format(lo, hi) if hi < 90 else "dec >= {}".format(lo)
    n = archive.tap_query("SELECT COUNT(*) AS n FROM {} WHERE {}".format(table, where))[0]["n"]
    if n == 0:
        return []
    if n > MAX_ROWS:
        mid = (lo + hi) / 2.0
        return fetch_rows(key, table, cols, lo, mid) + fetch_rows(key, table, cols, mid, hi)
    sel = ", ".join("{} AS {}".format(expr, alias) for alias, expr in cols)
    rows = archive.tap_query("SELECT {}, ra, dec FROM {} WHERE {}".format(sel, table, where))
    if len(rows) != n:
        raise RuntimeError("{} band [{}, {}): expected {} rows, got {}".format(key, lo, hi, n, len(rows)))
    return rows


def fetch_all(key, table, cols):
    rows = []
    for lo in range(-90, 90, 10):
        rows += fetch_rows(key, table, cols, lo, lo + 10)
    return rows


def write_cells(directory, cells):
    os.makedirs(directory, exist_ok=True)
    for idx, recs in cells.items():
        with open(os.path.join(directory, "{}.json".format(idx)), "w") as f:
            json.dump(recs, f, separators=(",", ":"), ensure_ascii=False)


# ---------------------------------------------------------------- observations

def mission_columns(m):
    if m.get("grouped"):
        return [("oid", "obs_id"), ("t", "start_date"), ("expo", "(end_date - start_date) * 86400")]
    cols = []
    for alias, key in (("oid", "id"), ("tname", "name"), ("t", "time"), ("expo", "expo"),
                       ("st", "status"), ("pub", "public")):
        if m.get(key):
            cols.append((alias, m[key]))
    for alias, col in m.get("extra", {}).items():
        cols.append(("x_" + alias, col))
    return cols


def build_mission(m, out_dir):
    t0 = time.time()
    rows = fetch_all(m["key"], m["table"], mission_columns(m))
    extra_keys = list(m.get("extra", {}).keys())
    statuses = sorted({(r.get("st") or "").lower() for r in rows})
    st_index = {s: k for k, s in enumerate(statuses)}
    cells = {}
    for r in rows:
        if r.get("ra") is None or r.get("dec") is None:
            continue
        oid = r.get("oid")
        rec = [rnd(r["ra"], 5), rnd(r["dec"], 5), "" if oid is None else str(oid), rnd(r.get("t"), 5),
               None if r.get("expo") is None else int(round(float(r["expo"])))]
        if not m.get("grouped"):
            rec += [r.get("tname") or "", st_index[(r.get("st") or "").lower()], rnd(r.get("pub"), 3)]
            rec += [r.get("x_" + k) for k in extra_keys]
        cells.setdefault(cell_index(r["ra"], r["dec"]), []).append(rec)
    write_cells(os.path.join(out_dir, "data", m["key"]), cells)
    info = dict(key=m["key"], label=m["label"], radius=m["radius"], grouped=bool(m.get("grouped")),
                statuses=statuses, extra_keys=extra_keys, expo_note=m.get("expo_note"),
                link=m.get("link"), rows=len(rows), cells=sorted(cells))
    print("{:<28} {:>8} rows {:>6} cells {:6.1f}s".format(m["label"], len(rows), len(cells), time.time() - t0),
          flush=True)
    return info


def build_surveys(out_dir):
    erass = archive.tap_query("SELECT obsid AS tile, ra_min, ra_max, dec_min, dec_max FROM erassmastr")
    rass = archive.tap_query("SELECT seq_id AS fid, ra, dec FROM rassmaster")
    if len(erass) < 1000 or len(rass) < 1000:
        raise RuntimeError("survey tables look incomplete")
    data = dict(
        erass1=[[str(r["tile"]), r["ra_min"], r["ra_max"], r["dec_min"], r["dec_max"]] for r in erass],
        rass=[[str(r["fid"]), rnd(r["ra"], 4), rnd(r["dec"], 4)] for r in rass],
    )
    with open(os.path.join(out_dir, "data", "surveys.json"), "w") as f:
        json.dump(data, f, separators=(",", ":"))
    print("surveys: {} eRASS1 tiles, {} RASS fields".format(len(erass), len(rass)), flush=True)


def build_observations(out_dir):
    os.makedirs(os.path.join(out_dir, "data"), exist_ok=True)
    with ThreadPoolExecutor(max_workers=4) as ex:
        missions = list(ex.map(lambda m: build_mission(m, out_dir), archive.MISSIONS))
    build_surveys(out_dir)
    manifest = dict(
        updated=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        cell_deg=CELL, n_ra=N_RA, n_dec=N_DEC,
        done_statuses=sorted(archive.DONE_STATUSES),
        missions=missions,
        surveys=[dict(key=s["key"], label=s["label"]) for s in archive.SURVEYS],
    )
    with open(os.path.join(out_dir, "data", "manifest.json"), "w") as f:
        json.dump(manifest, f, separators=(",", ":"))
    print("observations: total rows", sum(m["rows"] for m in missions), flush=True)


# ---------------------------------------------------------------- catalogs

def catalog_columns(c):
    cols = [] if c.get("noname") else [("cname", "name")]
    cols.append(("flux", c["flux"]))
    for k in ("err", "eneg", "epos", "lo", "hi", "poserr"):
        if c.get(k):
            cols.append((k, c[k]))
    if c.get("rate_err"):
        cols += [("rate", c["rate_err"][0]), ("rate_e", c["rate_err"][1])]
    cols += [("p{}".format(i), col) for i, (_, col) in enumerate(c["params"])]
    return cols


def flux_errors(c, r):
    f = r.get("flux")
    if f is None:
        return None, None
    if c.get("err"):
        return r.get("err"), r.get("err")
    if c.get("eneg"):
        neg = r.get("eneg")
        return (abs(neg) if neg is not None else None), r.get("epos")
    if c.get("lo"):
        lo, hi = r.get("lo"), r.get("hi")
        return (f - lo if lo is not None else None), (hi - f if hi is not None else None)
    if c.get("rate_err"):
        rate, rate_e = r.get("rate"), r.get("rate_e")
        if rate and rate_e is not None:
            e = abs(f * rate_e / rate)
            return e, e
    return None, None


def build_catalog(c, out_dir):
    t0 = time.time()
    rows = fetch_all(c["key"], c["table"], catalog_columns(c))
    cells = {}
    for r in rows:
        if r.get("ra") is None or r.get("dec") is None:
            continue
        eneg, epos = flux_errors(c, r)
        rec = [rnd(r["ra"], 5), rnd(r["dec"], 5), r.get("cname") or "",
               sig(r.get("poserr"), 3), sig(r.get("flux")), sig(eneg, 3), sig(epos, 3)]
        rec += [sig(r.get("p{}".format(i))) for i in range(len(c["params"]))]
        cells.setdefault(cell_index(r["ra"], r["dec"]), []).append(rec)
    write_cells(os.path.join(out_dir, "data", "cat", c["key"]), cells)
    info = {k: c.get(k) for k in ("key", "label", "kind", "band", "unit", "match", "poserr_note")}
    info.update(params=[p[0] for p in c["params"]], rows=len(rows), cells=sorted(cells))
    print("{:<28} {:>8} rows {:>6} cells {:6.1f}s".format(c["label"], len(rows), len(cells), time.time() - t0),
          flush=True)
    return info


def build_catalogs(out_dir):
    os.makedirs(os.path.join(out_dir, "data", "cat"), exist_ok=True)
    with ThreadPoolExecutor(max_workers=4) as ex:
        catalogs = list(ex.map(lambda c: build_catalog(c, out_dir), CATALOGS))
    manifest = dict(updated=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
                    cell_deg=CELL, n_ra=N_RA, n_dec=N_DEC, catalogs=catalogs)
    with open(os.path.join(out_dir, "data", "cat", "manifest.json"), "w") as f:
        json.dump(manifest, f, separators=(",", ":"))
    print("catalogs: total rows", sum(c["rows"] for c in catalogs), flush=True)


def main(argv):
    args = [a for a in argv[1:] if not a.startswith("--")]
    flags = {a for a in argv[1:] if a.startswith("--")}
    unknown = flags - {"--observations", "--catalogs"}
    if len(args) != 1 or unknown:
        print(__doc__)
        return 1
    if not flags:
        flags = {"--observations", "--catalogs"}
    if "--observations" in flags:
        build_observations(args[0])
    if "--catalogs" in flags:
        build_catalogs(args[0])
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
