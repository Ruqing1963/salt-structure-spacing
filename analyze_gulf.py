"""
analyze_gulf.py  --  U.S. Gulf Coast salt diapirs (USGS gcdiapirg)

Reads ../data/gulf_salt_diapirs.csv, computes the Clark-Evans table for the
whole basin and for individual provinces, and regenerates:
  ../figures/gulf_overview.{pdf,png}     (Fig. 1)
  ../figures/gulf_pair_correlation.{pdf,png}  (Fig. 3)
Numeric output is appended to ../results/results.txt.
"""
import os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import stats_common as sc

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data", "gulf_salt_diapirs.csv")
FIG = os.path.join(HERE, "..", "figures")
RES = os.path.join(HERE, "..", "results")
PURPLE = "#534AB7"

lon, lat = [], []
with open(DATA) as f:
    for row in csv.DictReader(f):
        lon.append(float(row["lon"])); lat.append(float(row["lat"]))
lon, lat = np.array(lon), np.array(lat)


def subset(a, b, c, d):
    m = (lon >= a) & (lon <= b) & (lat >= c) & (lat <= d)
    return sc.project_lonlat_km(lon[m], lat[m])


def main():
    rng = np.random.default_rng(101)
    lines = ["", "=" * 60, "GULF OF MEXICO  (USGS gcdiapirg, N=%d)" % len(lon), "=" * 60,
             f"{'province / subset':38} {'N':>4} {'NND':>5} {'R':>5} {'p':>6} verdict"]
    provinces = [
        ("Whole basin (all pooled)", (-100, -84, 15.5, 33)),
        ("N. Gulf shelf/slope + coast", (-97, -88, 26.5, 30)),
        ("Texas-Louisiana shelf core", (-96, -90, 27, 29)),
        ("Bay of Campeche", (-96, -90, 16, 22)),
        ("Onshore interior basins", (-96, -88, 30, 33)),
    ]
    for name, reg in provinces:
        P = subset(*reg)
        R, p, _, med, _ = sc.clark_evans(P, S=499, rng=rng)
        v = "clustered" if R < 0.95 else ("regular" if R > 1.05 else "random")
        lines.append(f"{name:38} {len(P):>4} {med:5.1f} {R:5.2f} {p:6.3f} {v}")
    print("\n".join(lines))
    with open(os.path.join(RES, "results.txt"), "a") as f:
        f.write("\n".join(lines) + "\n")

    # ---- Fig 1: overview map + density ----
    lat0 = lat.mean()
    fig, ax = plt.subplots(1, 2, figsize=(11, 5))
    ax[0].scatter(lon, lat, s=4, c=PURPLE, alpha=0.6, lw=0)
    ax[0].set_xlabel("lon"); ax[0].set_ylabel("lat")
    ax[0].set_title(f"Gulf Coast salt diapirs (N={len(lon)})")
    ax[0].set_aspect(1 / np.cos(np.radians(lat0)))
    ax[1].hist2d(lon, lat, bins=60, cmap="viridis")
    ax[1].set_xlabel("lon"); ax[1].set_ylabel("lat")
    ax[1].set_title("binned counts (sub-basins / gradient)")
    ax[1].set_aspect(1 / np.cos(np.radians(lat0)))
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"gulf_overview.{ext}"), dpi=130, bbox_inches="tight")
    plt.close(fig)

    # ---- Fig 3: Texas-Louisiana pair correlation with envelopes ----
    P = subset(-96, -90, 27, 29); V, A, path = sc.hull_window(P); n = len(P)
    nnd = np.median(np.sort(np.sqrt(((P[:, None] - P[None]) ** 2).sum(2)), 1)[:, 1])
    redges = np.arange(0, 40, 2.0)
    rc, g = sc.pair_correlation(P, redges)
    rng = np.random.default_rng(11)
    gC = np.array([sc.pair_correlation(sc.csr_in_hull(n, V, path, rng), redges)[1]
                   for _ in range(999)])
    isamp = sc.inhomogeneous_sampler(P, path, 2 * nnd, rng)
    gI = np.array([sc.pair_correlation(isamp(n), redges)[1] for _ in range(399)])
    mask = rc <= 16
    p_hom = sc.global_envelope_p(g, gC, mask)
    p_inh = sc.global_envelope_p(g, gI, mask)
    glo, ghi = np.nanpercentile(gC, 2.5, 0), np.nanpercentile(gC, 97.5, 0)
    ilo, ihi = np.nanpercentile(gI, 2.5, 0), np.nanpercentile(gI, 97.5, 0)
    fig, ax = plt.subplots(figsize=(6, 4.2))
    ax.fill_between(rc, glo, ghi, color="#bbb", alpha=0.5, label="CSR 95%")
    ax.fill_between(rc, ilo, ihi, color="#D9943A", alpha=0.25, label="inhomogeneous 95%")
    ax.plot(rc, g, "-o", color=PURPLE, ms=3, lw=2, label="observed g(r)")
    ax.axhline(1, ls="--", color="k", lw=0.8); ax.set_xlim(0, 30)
    ax.set_xlabel("r (km)"); ax.set_ylabel("g(r)"); ax.legend(fontsize=8, frameon=False)
    ax.set_title(f"Texas-Louisiana shelf diapirs: hard-core + characteristic spacing\n"
                 f"global p(CSR)={p_hom:.3f}, p(inhom)={p_inh:.3f}")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"gulf_pair_correlation.{ext}"), dpi=130, bbox_inches="tight")
    plt.close(fig)
    msg = f"\nGulf TX-LA g(r) global test: p(CSR)={p_hom:.4f}, p(inhom)={p_inh:.4f}"
    print(msg)
    with open(os.path.join(RES, "results.txt"), "a") as f:
        f.write(msg + "\n")
    print("Gulf figures written to ../figures/")


if __name__ == "__main__":
    main()
