# Data provenance and licences

This repository redistributes only **derived** coordinate data (structure centroids computed from
the original polygon outlines). The raw shapefiles are **not** redistributed here; download them
from the official sources below. The `code/parse_shapefiles.py` script regenerates the derived CSVs
from those raw files.

## U.S. Gulf Coast — `data/gulf_salt_diapirs.csv`
- **Source:** U.S. Geological Survey, *Salt Diapirs in the Gulf Coast* (`gcdiapirg`).
- **DOI:** https://doi.org/10.5066/P9KH3WB4
- **Original compilation:** Martin, R.G. (1980), *Distribution of salt structures, Gulf of Mexico*,
  USGS Misc. Field Studies Map MF-1213.
- **Licence:** U.S. Government work / public domain.
- **Derived field:** area-weighted centroid (lon, lat) and polygon area of each diapir outline.

## North German Basin — `data/ngb_salt_structures.csv`
- **Source:** Bundesanstalt für Geowissenschaften und Rohstoffe (BGR),
  *InSpEE-Salzstrukturen* (salt-structure outlines, 1:500 000).
- **Download:** https://www.bgr.bund.de/ (InSpEE / Salzstrukturen)
- **Base mapping:** Reinhold, K., Krull, P., Kockel, F. (2008), *Salzstrukturen Norddeutschlands
  1:500 000*, BGR, Hannover.
- **Licence:** open geodata licence Germany (Datenlizenz Deutschland / geonutzv).
  **Required attribution:** "InSpEE-Salzstrukturen, © BGR, Hannover, 2015".
- **Derived fields:** area-weighted centroid (UTM-32N, km), structure class, internal-structure
  type, and structure name.

Please retain the attributions above when using the derived datasets.
