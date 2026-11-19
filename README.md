# X-ray Archive Finder

Enter a source name (any SIMBAD identifier) or coordinates and see how many observations each X-ray mission has of it in the NASA HEASARC archive, with dates, exposures and observation IDs. It also shows the field in optical, X-ray and radio surveys, ROSAT/eROSITA all-sky survey coverage, and the SDSS spectrum when one exists.

Live: https://mustafaturansaglam.com/tools/xray-archive-finder/

HEASARC does not allow browser (CORS) queries, so a nightly GitHub Actions job (`build_index.py`) downloads the HEASARC master tables, splits them into 2°×2° sky cells and publishes them with the static page (`static/index.html`) on GitHub Pages. The page reads only the cells around the requested position; name resolution (SIMBAD) and SDSS spectra are queried live from the browser.

`archive.py` queries HEASARC directly from the command line:

```
pip install -r requirements.txt
python archive.py "Cyg X-1"
```
