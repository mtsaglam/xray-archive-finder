import os
import threading
import time
from collections import defaultdict, deque

from flask import Flask, jsonify, request, send_from_directory

import archive

HERE = os.path.dirname(os.path.abspath(__file__))
app = Flask(__name__, static_folder=os.path.join(HERE, "static"), static_url_path="/static")
app.json.ensure_ascii = False

ALLOWED_ORIGINS = {o.strip().rstrip("/") for o in os.environ.get(
    "ALLOWED_ORIGINS",
    "https://mustafaturansaglam.com,https://www.mustafaturansaglam.com,https://mtsaglam.github.io,"
    "http://127.0.0.1:5050,http://localhost:5050,http://127.0.0.1:4000,http://localhost:4000",
).split(",") if o.strip()}

# A full search makes ~21 API calls, so this allows a few searches per minute per client.
RATE_LIMIT = int(os.environ.get("RATE_LIMIT_PER_MIN", "150"))
_hits = defaultdict(deque)
_hits_lock = threading.Lock()


def _client_ip():
    fwd = request.headers.get("X-Forwarded-For", "")
    return fwd.split(",")[0].strip() if fwd else (request.remote_addr or "?")


@app.before_request
def _rate_limit():
    if not request.path.startswith("/api/") or request.path == "/api/health":
        return None
    now = time.time()
    with _hits_lock:
        q = _hits[_client_ip()]
        while q and now - q[0] > 60:
            q.popleft()
        if len(q) >= RATE_LIMIT:
            return _err("Too many requests, please wait a minute", 429)
        q.append(now)
        if len(_hits) > 5000:
            for ip in [ip for ip, d in _hits.items() if not d or now - d[-1] > 60]:
                del _hits[ip]
    return None


@app.after_request
def _cors(resp):
    origin = (request.headers.get("Origin") or "").rstrip("/")
    if request.path.startswith("/api/") and origin in ALLOWED_ORIGINS:
        resp.headers["Access-Control-Allow-Origin"] = origin
        resp.headers["Vary"] = "Origin"
    return resp


def _err(msg, code=400):
    return jsonify(error=msg), code


def _coords():
    try:
        ra = float(request.args["ra"])
        dec = float(request.args["dec"])
    except (KeyError, ValueError):
        raise ValueError("ra/dec parameters (degrees) are required")
    if not (0 <= ra < 360 and -90 <= dec <= 90):
        raise ValueError("ra/dec out of range")
    return ra, dec


@app.get("/")
def index():
    return send_from_directory(app.static_folder, "index.html")


@app.get("/api/health")
def health():
    return jsonify(ok=True)


@app.get("/api/missions")
def missions():
    return jsonify(
        missions=[dict(key=m["key"], label=m["label"], radius=m["radius"]) for m in archive.MISSIONS],
        surveys=[dict(key=s["key"], label=s["label"]) for s in archive.SURVEYS],
    )


@app.get("/api/resolve")
def resolve():
    q = request.args.get("q", "").strip()
    if not q:
        return _err("Enter a source name or coordinates")
    if len(q) > 200:
        return _err("Query too long")
    try:
        return jsonify(archive.resolve(q))
    except archive.ArchiveError as e:
        return _err(str(e), 404)
    except Exception as e:  # noqa: BLE001
        return _err("Resolver error: {}".format(e), 502)


@app.get("/api/mission/<key>")
def mission(key):
    if key not in archive.MISSION_BY_KEY:
        return _err("Unknown mission", 404)
    try:
        ra, dec = _coords()
        radius = request.args.get("radius")
        radius = float(radius) if radius else None
        if radius is not None and not (0 < radius <= 600):
            return _err("Radius must be 0–600 arcmin")
    except ValueError as e:
        return _err(str(e))
    try:
        return jsonify(archive.query_mission(key, ra, dec, radius))
    except Exception as e:  # noqa: BLE001
        return _err(str(e), 502)


@app.get("/api/survey/<key>")
def survey(key):
    if key not in archive.SURVEY_BY_KEY:
        return _err("Unknown survey", 404)
    try:
        ra, dec = _coords()
    except ValueError as e:
        return _err(str(e))
    try:
        return jsonify(archive.query_survey(key, round(ra, 6), round(dec, 6)))
    except Exception as e:  # noqa: BLE001
        return _err(str(e), 502)


@app.get("/api/sdss")
def sdss():
    try:
        ra, dec = _coords()
    except ValueError as e:
        return _err(str(e))
    try:
        return jsonify(spectra=archive.sdss_spectra(round(ra, 6), round(dec, 6)))
    except Exception as e:  # noqa: BLE001
        return _err("SDSS: {}".format(e), 502)


if __name__ == "__main__":
    port = int(os.environ.get("PORT", "5050"))
    print("\n  X-ray Archive Finder:  http://127.0.0.1:{}\n".format(port))
    app.run(host="127.0.0.1", port=port, threaded=True, debug=False)
