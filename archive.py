"""SIMBAD name resolution + HEASARC X-ray archive search.

CLI usage:  python archive.py "Cyg X-1"   |   python archive.py "299.59 35.20"
"""
import io
import re
import math
import sys
import time
import warnings
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone
from functools import lru_cache

warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402
import requests  # noqa: E402
from astropy import units as u  # noqa: E402
from astropy.coordinates import SkyCoord  # noqa: E402
from astropy.io.votable import parse  # noqa: E402

HEASARC_TAP = "https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync"
SIMBAD_MIRRORS = [
    "https://simbad.cds.unistra.fr/simbad",
    "https://simbad.cfa.harvard.edu/simbad",
]
SESAME = "https://cds.unistra.fr/cgi-bin/nph-sesame/-oxp/SNV?"
SDSS_SQL = "https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch"
SDSS_BASE = "https://skyserver.sdss.org/dr17"
UA = {"User-Agent": "xray-archive-finder/1.0"}
TIMEOUT = 90

# "yes" is BeppoSAX obs_perf; anything else (scheduled, accepted, unobserved, "no"...) = not observed.
DONE_STATUSES = {"archived", "processed", "observed", "preprocessed", "public", "yes", ""}

# radius: default search radius in arcmin, roughly the instrument field of view.
MISSIONS = [
    dict(key="chandra", label="Chandra", table="chanmaster", radius=10, id="obsid", name="name",
         time="time", expo="exposure", status="status", public="public_date",
         extra={"detector": "detector", "grating": "grating"},
         link="https://cda.harvard.edu/chaser/startViewer.do?menuItem=details&obsid={id}"),
    dict(key="xmm", label="XMM-Newton", table="xmmmaster", radius=15, id="obsid", name="name",
         time="time", expo="duration", status="status", public="public_date",
         link="https://nxsa.esac.esa.int/nxsa-web/#obsid={id}"),
    dict(key="swift", label="Swift", table="swiftmastr", radius=12, id="obsid", name="name",
         time="start_time", expo="xrt_exposure", expo_note="XRT"),
    dict(key="nustar", label="NuSTAR", table="numaster", radius=6, id="obsid", name="name",
         time="time", expo="exposure_a", status="status", public="public_date"),
    dict(key="nicer", label="NICER", table="nicermastr", radius=3, id="obsid", name="name",
         time="time", expo="exposure", public="public_date"),
    dict(key="xrism", label="XRISM", table="xrismmastr", radius=10, id="obsid", name="name",
         time="time", expo="exposure", status="status", public="public_date"),
    dict(key="ixpe", label="IXPE", table="ixmaster", radius=6, id="obsid", name="name",
         time="time", expo="exposure_1", status="status", public="public_date"),
    dict(key="erosita", label="eROSITA (pointed)", table="erosmaster", radius=30, id="obsid",
         name="name", extra={"type": "obs_type"}),
    dict(key="integral", label="INTEGRAL", table="intscwpub", radius=300, grouped=True),
    dict(key="suzaku", label="Suzaku", table="suzamaster", radius=9, id="obsid", name="name",
         time="time", expo="exposure", public="public_date"),
    dict(key="hitomi", label="Hitomi", table="hitomaster", radius=5, id="obsid", name="name",
         time="time", expo="exposure", public="public_date"),
    dict(key="rxte", label="RXTE", table="xtemaster", radius=30, id="obsid", name="target_name",
         time="time", expo="exposure", status="status"),
    dict(key="rosat", label="ROSAT (pointed)", table="rosmaster", radius=30, id="seq_id",
         name="name", time="start_time", expo="exposure", extra={"instrument": "instrument"}),
    dict(key="asca", label="ASCA", table="ascamaster", radius=20, id="sequence_number",
         name="target_name", time="time", expo="gis_exposure", status="status",
         public="public_date", expo_note="GIS"),
    dict(key="bepposax", label="BeppoSAX NFI", table="saxnfilog", radius=15, id="seq_id",
         name="name", time="start_date", expo="mecs_exposure", status="obs_perf",
         expo_note="MECS"),
    dict(key="exosat", label="EXOSAT", table="exomaster", radius=30, id="sequence_number",
         name="name", time="time", expo="(end_time - time) * 86400", expo_note="duration",
         extra={"instrument": "instrument"}),
    dict(key="ginga", label="Ginga LAC", table="gingalog", radius=60, name="name",
         time="time", expo="exposure"),
]
MISSION_BY_KEY = {m["key"]: m for m in MISSIONS}

SURVEYS = [
    dict(key="erass1", label="eROSITA eRASS1 (DR1)", table="erassmastr"),
    dict(key="rass", label="ROSAT All-Sky Survey", table="rassmaster"),
]
SURVEY_BY_KEY = {s["key"]: s for s in SURVEYS}


class ArchiveError(Exception):
    pass


# ---------------------------------------------------------------- helpers

def mjd_to_date(mjd):
    if mjd is None:
        return None
    return (datetime(1858, 11, 17) + timedelta(days=float(mjd))).strftime("%Y-%m-%d")


def today_mjd():
    now = datetime.now(timezone.utc).replace(tzinfo=None)
    return (now - datetime(1858, 11, 17)).total_seconds() / 86400.0


def _clean(v):
    """Convert numpy/masked VOTable values into JSON-friendly Python values."""
    if v is None or v is np.ma.masked:
        return None
    if hasattr(v, "item"):
        v = v.item()
    if isinstance(v, bytes):
        v = v.decode("utf-8", "replace")
    if isinstance(v, float) and v != v:
        return None
    if isinstance(v, str):
        v = v.strip()
    return v


def _post_with_retry(url, data, timeout, attempts=3):
    last = None
    for i in range(attempts):
        try:
            r = requests.post(url, data=data, headers=UA, timeout=timeout)
            if r.status_code < 500:
                return r
            last = "HTTP {}".format(r.status_code)
        except requests.RequestException as e:
            last = str(e)
        if i < attempts - 1:
            time.sleep(1.5 * (i + 1))
    raise ArchiveError("{} -> {}".format(url, last))


def tap_query(adql):
    r = _post_with_retry(HEASARC_TAP, {"REQUEST": "doQuery", "LANG": "ADQL", "QUERY": adql}, TIMEOUT)
    r.raise_for_status()
    vot = parse(io.BytesIO(r.content))
    for res in vot.resources:
        for info in res.infos:
            if info.name == "QUERY_STATUS" and info.value == "ERROR":
                raise ArchiveError("HEASARC: " + (info.content or "query error"))
    tables = list(vot.iter_tables())
    if not tables:
        return []
    arr = tables[0].array
    names = arr.dtype.names
    return [{n: _clean(row[n]) for n in names} for row in arr]


def cone(ra, dec, r_deg):
    return ("CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{:.7f},{:.7f},{:.6f}))=1"
            .format(ra, dec, r_deg))


def offset_expr(ra, dec):
    return "DISTANCE(POINT('ICRS',ra,dec),POINT('ICRS',{:.7f},{:.7f}))*60".format(ra, dec)


# ---------------------------------------------------------------- resolving

COORD_RE = re.compile(r"^[\s0-9+\-.:,hmsd°'\"]+$", re.I)


def parse_coords(text):
    """Return SkyCoord if text looks like coordinates, else None."""
    t = text.strip()
    if not COORD_RE.match(t) or not re.search(r"\d", t):
        return None
    t = t.replace(",", " ")
    parts = t.split()
    try:
        if len(parts) == 2 and all(re.match(r"^[+-]?\d+(\.\d*)?$", p) for p in parts):
            c = SkyCoord(float(parts[0]) * u.deg, float(parts[1]) * u.deg, frame="icrs")
        else:
            c = SkyCoord(t, unit=(u.hourangle, u.deg), frame="icrs")
    except Exception:
        return None
    return c


def _simbad_script(script):
    last_err = None
    for base in SIMBAD_MIRRORS:
        try:
            r = _post_with_retry(base + "/sim-script", {"script": script}, 40, attempts=2)
        except ArchiveError as e:
            last_err = str(e)
            continue
        if r.status_code != 200 or "<html" in r.text[:500].lower():
            last_err = "{} -> HTTP {}".format(base, r.status_code)
            continue
        return base, r.text
    raise ArchiveError("SIMBAD unreachable ({})".format(last_err))


SIMBAD_FMT = 'format object "%MAIN_ID|%COO(d;A D)|%OTYPE(V)|%IDLIST[%*;]"'


def _parse_simbad_line(line):
    parts = line.split("|")
    if len(parts) < 4:
        return None
    radec = parts[1].split()
    if len(radec) != 2:
        return None
    aliases = [" ".join(a.split()) for a in parts[3].split(";") if a.strip()]
    return dict(main_id=" ".join(parts[0].split()), ra=float(radec[0]), dec=float(radec[1]),
                otype=parts[2].strip(), aliases=aliases)


def simbad_by_name(name):
    script = "output console=off script=off\n{}\nquery id {}".format(SIMBAD_FMT, name)
    base, text = _simbad_script(script)
    if "::error::" in text:
        return base, None
    for line in text.splitlines():
        obj = _parse_simbad_line(line)
        if obj:
            return base, obj
    return base, None


def simbad_nearest(ra, dec, radius_arcsec=10):
    script = ("output console=off script=off\n{}\nquery coo {:.7f} {:+.7f} radius={}s"
              .format(SIMBAD_FMT, ra, dec, radius_arcsec))
    base, text = _simbad_script(script)
    if "::error::" in text:
        return None
    best = None
    for line in text.splitlines():
        obj = _parse_simbad_line(line)
        if obj:
            sep = SkyCoord(ra * u.deg, dec * u.deg).separation(
                SkyCoord(obj["ra"] * u.deg, obj["dec"] * u.deg)).arcsec
            obj["sep_arcsec"] = round(sep, 2)
            if best is None or sep < best["sep_arcsec"]:
                best = obj
    return best


def sesame(name):
    """Fallback resolver (SIMBAD/NED/VizieR via CDS Sesame); coordinates only."""
    r = requests.get(SESAME + requests.utils.quote(name), headers=UA, timeout=30)
    m_ra = re.search(r"<jradeg>([-\d.]+)</jradeg>", r.text)
    m_de = re.search(r"<jdedeg>([-\d.]+)</jdedeg>", r.text)
    if not (m_ra and m_de):
        return None
    m_res = re.search(r'<Resolver name="([^"]+)"', r.text)
    m_on = re.search(r"<oname>([^<]+)</oname>", r.text)
    m_ot = re.search(r"<otype>([^<]+)</otype>", r.text)
    return dict(main_id=m_on.group(1) if m_on else name, ra=float(m_ra.group(1)),
                dec=float(m_de.group(1)), otype=m_ot.group(1) if m_ot else "",
                aliases=[], resolver=(m_res.group(1) if m_res else "Sesame"))


@lru_cache(maxsize=256)
def resolve(text):
    text = " ".join(text.split())
    if not text:
        raise ArchiveError("Empty query")
    c = parse_coords(text)
    out = dict(query=text)
    if c is not None:
        out.update(input_type="coord", ra=c.ra.deg, dec=c.dec.deg, main_id=None, otype=None,
                   aliases=[], resolver="coordinates")
        try:
            near = simbad_nearest(c.ra.deg, c.dec.deg)
        except ArchiveError:
            near = None
        out["nearest"] = near
    else:
        obj, base = None, None
        try:
            base, obj = simbad_by_name(text)
        except ArchiveError:
            obj = None
        if obj:
            obj["resolver"] = "SIMBAD"
            out["simbad_url"] = base + "/sim-id?Ident=" + requests.utils.quote(obj["main_id"])
        else:
            obj = sesame(text)
        if not obj:
            raise ArchiveError("Could not resolve '{}' (SIMBAD/NED/VizieR)".format(text))
        out.update(input_type="name", **obj)
    sc = SkyCoord(out["ra"] * u.deg, out["dec"] * u.deg, frame="icrs")
    out["ra_hms"] = sc.ra.to_string(u.hourangle, sep=":", precision=2, pad=True)
    out["dec_dms"] = sc.dec.to_string(u.deg, sep=":", precision=1, alwayssign=True, pad=True)
    out["l"] = round(sc.galactic.l.deg, 4)
    out["b"] = round(sc.galactic.b.deg, 4)
    return out


# ---------------------------------------------------------------- missions

def _build_query(m, ra, dec, r_deg):
    if m.get("grouped"):  # INTEGRAL: science windows grouped by observation id
        return ("SELECT obs_id AS oid, MIN(start_date) AS t, "
                "SUM((end_date - start_date) * 86400) AS expo, COUNT(*) AS nscw, "
                "MIN({off}) AS off, AVG(ra) AS ra, AVG(dec) AS dec "
                "FROM intscwpub WHERE {cone} GROUP BY obs_id"
                .format(off=offset_expr(ra, dec), cone=cone(ra, dec, r_deg)))
    cols = []
    for alias, key in (("oid", "id"), ("tname", "name"), ("t", "time"), ("expo", "expo"),
                       ("st", "status"), ("pub", "public")):
        if m.get(key):
            cols.append("{} AS {}".format(m[key], alias))
    for alias, col in m.get("extra", {}).items():
        cols.append("{} AS x_{}".format(col, alias))
    cols += ["ra", "dec", "{} AS off".format(offset_expr(ra, dec))]
    return "SELECT {} FROM {} WHERE {}".format(", ".join(cols), m["table"], cone(ra, dec, r_deg))


@lru_cache(maxsize=1024)
def _mission_cached(key, ra, dec, radius_arcmin):
    m = MISSION_BY_KEY[key]
    rows = tap_query(_build_query(m, ra, dec, radius_arcmin / 60.0))
    now = today_mjd()
    obs = []
    for i, r in enumerate(rows):
        status = (r.get("st") or "").lower()
        done = status in DONE_STATUSES
        pub = r.get("pub")
        proprietary = bool(done and pub is not None and pub > now)
        oid = r.get("oid")
        if oid is None:
            oid = "#{}".format(i + 1)
        expo = r.get("expo")
        o = dict(id=str(oid), name=r.get("tname") or "", mjd=r.get("t"), date=mjd_to_date(r.get("t")),
                 exposure_ks=round(expo / 1000.0, 2) if expo is not None else None,
                 offset_arcmin=round(r["off"], 2) if r.get("off") is not None else None,
                 status=status or ("archived" if done else ""), done=done,
                 proprietary=proprietary, public_date=mjd_to_date(pub) if pub else None,
                 ra=r.get("ra"), dec=r.get("dec"),
                 extra={k[2:]: v for k, v in r.items() if k.startswith("x_") and v not in (None, "")})
        if r.get("nscw") is not None:
            o["extra"]["ScW"] = r["nscw"]
        if m.get("link") and done:
            o["link"] = m["link"].format(id=oid)
        obs.append(o)
    obs.sort(key=lambda o: (o["mjd"] is None, o["mjd"] or 0))
    done_obs = [o for o in obs if o["done"]]
    dates = [o["date"] for o in done_obs if o["date"]]
    return dict(
        key=key, label=m["label"], table=m["table"], radius_arcmin=radius_arcmin,
        n_done=len(done_obs), n_planned=len(obs) - len(done_obs),
        n_proprietary=sum(o["proprietary"] for o in obs),
        n_scw=sum(o["extra"].get("ScW", 0) for o in obs) if m.get("grouped") else None,
        exposure_ks=round(sum(o["exposure_ks"] or 0 for o in done_obs), 1),
        expo_note=m.get("expo_note"),
        first=min(dates) if dates else None, last=max(dates) if dates else None,
        observations=obs,
    )


def query_mission(key, ra, dec, radius_arcmin=None):
    if key not in MISSION_BY_KEY:
        raise ArchiveError("Unknown mission: " + key)
    if radius_arcmin is None:
        radius_arcmin = MISSION_BY_KEY[key]["radius"]
    return _mission_cached(key, round(ra, 6), round(dec, 6), float(radius_arcmin))


def _in_box(ra, dec, ra_min, ra_max, dec_min, dec_max):
    width = (ra_max - ra_min) % 360
    return dec_min <= dec <= dec_max and (ra - ra_min) % 360 <= width


@lru_cache(maxsize=512)
def query_survey(key, ra, dec):
    s = SURVEY_BY_KEY[key]
    if key == "erass1":
        rows = tap_query("SELECT obsid AS tile, ra_min, ra_max, dec_min, dec_max FROM erassmastr WHERE "
                         + cone(ra, dec, 3.0))
        # DR1 is the German half of the sky: 179.944 < l < 359.944 deg.
        l = SkyCoord(ra * u.deg, dec * u.deg).galactic.l.deg
        in_half = 179.9442 < l < 359.9442
        tiles = [str(r["tile"]) for r in rows
                 if _in_box(ra, dec, r["ra_min"], r["ra_max"], r["dec_min"], r["dec_max"])] if in_half else []
        note = None if tiles else "outside_dr1"
    else:
        rows = tap_query("SELECT seq_id AS fid, ra, dec FROM rassmaster WHERE " + cone(ra, dec, 4.6))
        tiles = []
        for r in rows:
            dra = ((ra - r["ra"] + 180) % 360 - 180) * math.cos(math.radians(dec))
            if abs(dra) <= 3.2 and abs(dec - r["dec"]) <= 3.2:
                tiles.append(str(r["fid"]))
        note = None
    return dict(key=key, label=s["label"], covered=bool(tiles), fields=tiles, note=note)


# ---------------------------------------------------------------- SDSS

@lru_cache(maxsize=512)
def sdss_spectra(ra, dec, radius_arcmin=0.05):
    cmd = ("SELECT TOP 10 s.specobjid, s.ra, s.dec, s.class, s.subclass, s.z, s.zerr, "
           "s.plate, s.mjd, s.fiberid, s.run2d, s.survey, n.distance "
           "FROM dbo.fGetNearbySpecObjEq({:.7f},{:.7f},{}) n "
           "JOIN specobj s ON s.specobjid = n.specobjid ORDER BY n.distance"
           .format(ra, dec, radius_arcmin))
    r = requests.get(SDSS_SQL, params={"cmd": cmd, "format": "json"}, headers=UA, timeout=60)
    r.raise_for_status()
    rows = r.json()[0].get("Rows", [])
    out = []
    for s in rows:
        sid = str(s["specobjid"])
        out.append(dict(
            specobjid=sid, cls=s.get("class"), subclass=s.get("subclass"),
            z=s.get("z"), zerr=s.get("zerr"), plate=s.get("plate"), mjd=s.get("mjd"),
            fiberid=s.get("fiberid"), survey=s.get("survey"),
            distance_arcsec=round(s["distance"] * 60, 2),
            image=SDSS_BASE + "/en/get/SpecById.ashx?id=" + sid,
            explore=SDSS_BASE + "/VisualTools/explore/summary?sid=" + sid,
        ))
    return out


# ---------------------------------------------------------------- CLI

def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 1
    try:
        src = resolve(" ".join(argv[1:]))
    except ArchiveError as e:
        print("HATA:", e)
        return 2
    print("Kaynak : {}  ({})".format(src.get("main_id") or src["query"], src.get("otype") or "-"))
    print("Konum  : RA={:.5f} Dec={:+.5f}  [{} {}]  l={} b={}".format(
        src["ra"], src["dec"], src["ra_hms"], src["dec_dms"], src["l"], src["b"]))
    if src.get("aliases"):
        print("Adlar  : {} ad ({}...)".format(len(src["aliases"]), ", ".join(src["aliases"][:6])))
    print()
    print("{:<22}{:>8}{:>10}{:>12}{:>9}   {}".format("Uydu", "Gözlem", "Planlı", "Pozlama ks", "R(')", "Tarih aralığı"))

    def run(m):
        try:
            return query_mission(m["key"], src["ra"], src["dec"])
        except Exception as e:  # noqa: BLE001
            return dict(label=m["label"], error=str(e))

    with ThreadPoolExecutor(max_workers=8) as ex:
        results = list(ex.map(run, MISSIONS))
    total = 0
    for res in results:
        if "error" in res:
            print("{:<22}  HATA: {}".format(res["label"], res["error"][:80]))
            continue
        total += res["n_done"]
        span = "{} – {}".format(res["first"], res["last"]) if res["first"] else ""
        print("{:<22}{:>8}{:>10}{:>12}{:>9.0f}   {}".format(
            res["label"], res["n_done"], res["n_planned"] or "", res["exposure_ks"],
            res["radius_arcmin"], span))
    print("\nToplam gözlem: {}".format(total))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
