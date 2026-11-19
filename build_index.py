"""Download the HEASARC X-ray master tables and write them as static sky-cell JSON files.

The browser cannot query HEASARC directly (no CORS), so this runs daily in GitHub Actions and
the static UI reads only the cells around the requested position.

Usage:  python build_index.py OUT_DIR
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


def cell_index(ra, dec):
    i = min(N_DEC - 1, max(0, int(math.floor((dec + 90.0) / CELL))))
    j = int(math.floor((ra % 360.0) / CELL)) % N_RA
    return i * N_RA + j


def columns(m):
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


def fetch_band(m, lo, hi):
    """Rows with lo <= dec < hi, splitting the band until each query is below the row cap."""
    where = "dec >= {} AND dec < {}".format(lo, hi) if hi < 90 else "dec >= {}".format(lo)
    n = archive.tap_query("SELECT COUNT(*) AS n FROM {} WHERE {}".format(m["table"], where))[0]["n"]
    if n == 0:
        return []
    if n > MAX_ROWS:
        mid = (lo + hi) / 2.0
        return fetch_band(m, lo, mid) + fetch_band(m, mid, hi)
    sel = ", ".join("{} AS {}".format(expr, alias) for alias, expr in columns(m))
    rows = archive.tap_query("SELECT {}, ra, dec FROM {} WHERE {}".format(sel, m["table"], where))
    if len(rows) != n:
        raise RuntimeError("{} band [{}, {}): expected {} rows, got {}".format(m["key"], lo, hi, n, len(rows)))
    return rows


def rnd(v, nd):
    return None if v is None else round(float(v), nd)


def build_mission(m, out_dir):
    t0 = time.time()
    rows = []
    for lo in range(-90, 90, 10):
        rows += fetch_band(m, lo, lo + 10)
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
    mdir = os.path.join(out_dir, "data", m["key"])
    os.makedirs(mdir, exist_ok=True)
    for idx, recs in cells.items():
        with open(os.path.join(mdir, "{}.json".format(idx)), "w") as f:
            json.dump(recs, f, separators=(",", ":"), ensure_ascii=False)
    info = dict(key=m["key"], label=m["label"], radius=m["radius"], grouped=bool(m.get("grouped")),
                statuses=statuses, extra_keys=extra_keys, expo_note=m.get("expo_note"),
                link=m.get("link"), rows=len(rows), cells=sorted(cells))
    print("{:<20} {:>7} rows {:>5} cells {:6.1f}s".format(m["label"], len(rows), len(cells), time.time() - t0),
          flush=True)
    return info


def build_surveys(out_dir):
    erass = archive.tap_query("SELECT obsid AS tile, ra_min, ra_max, dec_min, dec_max FROM erassmastr")
    rass = archive.tap_query("SELECT seq_id AS fid, ra, dec FROM rassmaster")
    data = dict(
        erass1=[[str(r["tile"]), r["ra_min"], r["ra_max"], r["dec_min"], r["dec_max"]] for r in erass],
        rass=[[str(r["fid"]), rnd(r["ra"], 4), rnd(r["dec"], 4)] for r in rass],
    )
    with open(os.path.join(out_dir, "data", "surveys.json"), "w") as f:
        json.dump(data, f, separators=(",", ":"))
    if len(erass) < 1000 or len(rass) < 1000:
        raise RuntimeError("survey tables look incomplete")
    print("surveys: {} eRASS1 tiles, {} RASS fields".format(len(erass), len(rass)), flush=True)


def main(out_dir):
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
    print("total rows:", sum(m["rows"] for m in missions))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1])
