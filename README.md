# X-ray Archive Finder

Enter a source name (any SIMBAD identifier) or coordinates and see how many observations each X-ray mission has of it in the NASA HEASARC archive, with dates, exposures and observation IDs. It also shows synced optical / X-ray / radio image panels, X-ray and radio catalogue fluxes (CSC 2.2, 5XMM-DR15, 2SXPS, eRASS1, 2RXS, NuSTAR, Swift-BAT, INTEGRAL/IBIS, NVSS, SUMSS, WENSS from HEASARC; VLASS, FIRST, RACS, TGSS, LoTSS, GLEAM live from VizieR), ROSAT/eROSITA all-sky survey coverage, and the SDSS spectrum when one exists. A batch mode runs a list of up to 500 names or coordinates from a file and exports a CSV. The observation lists include PI and proposal titles and an interactive timeline, and a generated `wget` script downloads the public data of every listed observation from the HEASARC archive. A long-term X-ray light curve combines per-observation fluxes from the Chandra Source Catalog 2.1 (queried live), the 5XMM-DR15 stacked catalogue, the XMM slew survey, eRASS1 and 2RXS.

Live: https://mustafaturansaglam.com/tools/xray-archive-finder/

HEASARC does not allow browser (CORS) queries, so a nightly GitHub Actions job (`build_index.py`) downloads the HEASARC master tables (and, weekly, the source catalogues), splits them into 2°×2° sky cells and publishes them with the static page (`static/index.html`) on GitHub Pages. The page reads only the cells around the requested position; name resolution (SIMBAD) and SDSS spectra are queried live from the browser.

`archive.py` queries HEASARC directly from the command line:

```
pip install -r requirements.txt
python archive.py "Cyg X-1"
```
