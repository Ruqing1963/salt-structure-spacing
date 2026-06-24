# Regular spacing of salt structures as a Rayleigh–Taylor signature

Edge-corrected point-pattern and structure-factor analysis of salt-structure spacing in the
**U.S. Gulf Coast** (USGS *gcdiapirg*) and **North German Basin** (BGR *InSpEE*) salt provinces.

This repository contains the data, code, figures and paper for the study. Within coherent
loading-driven provinces the salt-structure populations are significantly **regularly spaced**
(Clark–Evans *R* = 1.10–1.47, *p* ≈ 0.002, robust to a density-gradient null), with hard-core
exclusion and a 12–15 km characteristic spacing — the map-view signature of Rayleigh–Taylor
wavelength selection plus salt-withdrawal exclusion. Regularity strengthens with salt thickness
(Glückstadt Graben > central basin > Gulf shelf). The order is robust at short and intermediate
range; long-range hyperuniformity cannot be confirmed because the intrinsic basin density gradient
masks the long-wavelength limit.

## Repository layout

```
salt-structure-spacing/
├── code/                       analysis code (pure NumPy/SciPy/Matplotlib)
│   ├── stats_common.py         shared point-pattern & structure-factor estimators
│   ├── parse_shapefiles.py     regenerate the derived CSVs from raw shapefiles
│   ├── analyze_gulf.py         Gulf: Clark–Evans table + Figs 1, 3
│   ├── analyze_ngb.py          N. Germany: type map, NNI scan, Figs 2, 4, 5
│   └── make_summary_table.py   Table 1 -> results/clark_evans_summary.csv
├── data/                       derived coordinate datasets (see data/README.md)
│   ├── gulf_salt_diapirs.csv   1253 diapir centroids (lon, lat, polygon area)
│   └── ngb_salt_structures.csv 693 structure centroids (UTM-32N km, type, name)
├── figures/                    generated figures (.pdf and .png)
├── results/                    numeric output (results.txt, clark_evans_summary.csv)
├── paper/                      manuscript (PDF + LaTeX source + figures)
├── requirements.txt
├── LICENSE                     MIT (code)
└── DATA_LICENSE.md             data provenance and licences
```

## Reproducing the analysis

Requires Python 3.9+ with NumPy, SciPy and Matplotlib (`pip install -r requirements.txt`).

```bash
cd code
python analyze_gulf.py          # Gulf table + gulf_overview, gulf_pair_correlation
python analyze_ngb.py           # N. German scan + 3 figures
python make_summary_table.py    # results/clark_evans_summary.csv  (Table 1)
```

Figures are written to `../figures/` (both `.pdf` and `.png`); numeric results are appended to
`../results/results.txt`. The analysis scripts read the derived CSVs in `data/` and do **not**
require the raw shapefiles.

To regenerate the derived CSVs from the original open shapefiles (download separately — see
`data/README.md`):

```bash
cd code
python parse_shapefiles.py gulf /path/to/gcdiapirg.shp
python parse_shapefiles.py ngb  /path/to/InSpEE_unzipped_directory
```

## Method summary

| statistic | scale probed | what it measures |
|---|---|---|
| Clark–Evans *R* (MC-CSR null) | short range | nearest-neighbour dispersion vs. clustering |
| pair correlation *g(r)* | short range | hard-core exclusion + characteristic spacing |
| Ripley *L(r) − r* | short–intermediate | cumulative regularity (negative = regular) |
| structure factor *S(k)* | intermediate–long | suppression of density fluctuations |
| number variance *σ²(R)* | long range | fluctuation scaling (hyperuniform ~ *R*¹) |

All tests use Monte-Carlo **CSR** and **inhomogeneous-Poisson** (density-gradient) null models;
global significance uses the studentised maximum-deviation envelope test. Because a smooth basin
density gradient biases the statistics toward clustering, the detected regularity is conservative.

## Citation

If you use this code or the derived data, please cite the paper (`paper/salt_paper.pdf`) and the
original data sources listed in `DATA_LICENSE.md`.

## Licence

Code: MIT (see `LICENSE`). Data: see `DATA_LICENSE.md` (USGS public domain; BGR open geodata
licence, attribution required).
