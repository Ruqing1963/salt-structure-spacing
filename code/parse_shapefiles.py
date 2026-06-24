"""
parse_shapefiles.py
Pure-Python (no GIS libraries) parsers that regenerate the derived CSVs in
../data/ from the raw open shapefiles. Provided for full reproducibility from
source; the analysis scripts read the CSVs and do not require this step.

Raw data (download separately, see ../data/README.md):
  USGS gcdiapirg  (Gulf Coast salt diapirs)      -> gulf_salt_diapirs.csv
  BGR InSpEE      (North German salt structures) -> ngb_salt_structures.csv

Usage:
  python parse_shapefiles.py gulf  /path/to/gcdiapirg.shp
  python parse_shapefiles.py ngb   /path/to/InSpEE_unzipped_dir
"""
import sys, struct, csv, os
import numpy as np


def _read_polygons(shp):
    """Return list of (cx, cy, area) area-weighted centroids of the exterior
    ring of every polygon record (shapefile type 5)."""
    d = open(shp, "rb").read(); pos = 100; out = []
    while pos < len(d):
        _, clen = struct.unpack(">ii", d[pos:pos + 8]); pos += 8
        rec = d[pos:pos + clen * 2]; pos += clen * 2
        if len(rec) < 44 or struct.unpack("<i", rec[0:4])[0] != 5:
            out.append(None); continue
        nP, nPt = struct.unpack("<ii", rec[36:44])
        parts = struct.unpack("<%di" % nP, rec[44:44 + 4 * nP])
        xy = np.frombuffer(rec[44 + 4 * nP:44 + 4 * nP + 16 * nPt],
                           dtype="<f8").reshape(-1, 2)
        end = parts[1] if nP > 1 else nPt
        ring = xy[parts[0]:end]; x, y = ring[:, 0], ring[:, 1]
        A = 0.5 * np.sum(x[:-1] * y[1:] - x[1:] * y[:-1])
        if abs(A) > 1e-9:
            cx = np.sum((x[:-1] + x[1:]) * (x[:-1] * y[1:] - x[1:] * y[:-1])) / (6 * A)
            cy = np.sum((y[:-1] + y[1:]) * (x[:-1] * y[1:] - x[1:] * y[:-1])) / (6 * A)
        else:
            cx, cy = x.mean(), y.mean()
        out.append((cx, cy, abs(A)))
    return out


def _read_dbf(dbf):
    d = open(dbf, "rb").read()
    nrec = struct.unpack("<i", d[4:8])[0]
    hs = struct.unpack("<h", d[8:10])[0]; rs = struct.unpack("<h", d[10:12])[0]
    fields = []; p = 32
    while d[p] != 0x0D:
        fields.append((d[p:p + 11].split(b"\x00")[0].decode("latin1"), d[p + 16]))
        p += 32
    recs = []
    for i in range(nrec):
        r = d[hs + i * rs: hs + (i + 1) * rs]; o = 1; row = {}
        for nm, fl in fields:
            row[nm] = r[o:o + fl].decode("latin1").strip(); o += fl
        recs.append(row)
    return recs


def parse_gulf(shp, out="../data/gulf_salt_diapirs.csv"):
    cent = _read_polygons(shp)
    with open(out, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["lon", "lat", "poly_area_deg2"])
        for c in cent:
            if c is None:
                continue
            w.writerow([f"{c[0]:.6f}", f"{c[1]:.6f}", f"{c[2]:.8g}"])
    print("wrote", out)


def parse_ngb(folder, out="../data/ngb_salt_structures.csv"):
    base = os.path.join(folder, "Salzstrukturen_Inspee__v1_poly")
    cent = _read_polygons(base + ".shp")
    recs = _read_dbf(base + ".dbf")
    itype = _read_dbf(os.path.join(folder, "Internbautyp___v1_poly.dbf"))
    name2type = {r["Bezeichnun"]: r["Internbaut"] for r in itype}
    agg = {}
    for c, r in zip(cent, recs):
        if c is None:
            continue
        nm = r["Bezeichnun"]; cx, cy, A = c
        a = agg.setdefault(nm, [0.0, 0.0, 0.0, r["Salzstrukt"]])
        a[0] += cx * A; a[1] += cy * A; a[2] += A
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["easting_km", "northing_km", "structure_class",
                    "internal_type", "name"])
        for nm, (sx, sy, sA, sk) in agg.items():
            it = name2type.get(nm, "NA").replace("\xc3\xbc", "ue")
            it = it.encode("latin1", "ignore").decode("latin1").replace("Ã¼", "ue")
            w.writerow([f"{sx / sA / 1000:.4f}", f"{sy / sA / 1000:.4f}",
                        sk, it, nm])
    print("wrote", out)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__); sys.exit(1)
    if sys.argv[1] == "gulf":
        parse_gulf(sys.argv[2])
    elif sys.argv[1] == "ngb":
        parse_ngb(sys.argv[2])
    else:
        print(__doc__)
