# Derived data

These CSVs contain **structure centroids** derived from the original open polygon datasets.
See `../DATA_LICENSE.md` for sources, DOIs and required attributions.

| file | rows | columns |
|---|---|---|
| `gulf_salt_diapirs.csv` | 1253 | `lon, lat, poly_area_deg2` |
| `ngb_salt_structures.csv` | 693 | `easting_km, northing_km, structure_class, internal_type, name` |

`gulf_salt_diapirs.csv` is in geographic coordinates (NAD83); the analysis projects it to a local
kilometre frame. `ngb_salt_structures.csv` is already in ETRS89 / UTM zone 32N, expressed in
kilometres. `internal_type` values are transliterated to ASCII
(e.g. `Salzdiapir kompressiv ueberpraegt`).

To regenerate from the raw shapefiles, see `../code/parse_shapefiles.py`.
