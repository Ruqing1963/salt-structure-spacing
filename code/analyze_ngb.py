"""
analyze_ngb.py  --  North German Basin salt structures (BGR InSpEE)

Reads ../data/ngb_salt_structures.csv (coordinates already in km, UTM 32N) and
regenerates:
  ../figures/ngb_type_map.{pdf,png}          (Fig. 2)
  ../figures/ngb_pair_correlation_ripley.{pdf,png}  (Fig. 4)
  ../figures/ngb_structure_factor.{pdf,png}  (Fig. 5)
plus the Clark-Evans scan. Numeric output is appended to ../results/results.txt.
"""
import os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
import stats_common as sc

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data", "ngb_salt_structures.csv")
FIG = os.path.join(HERE, "..", "figures")
RES = os.path.join(HERE, "..", "results")
PURPLE = "#534AB7"

X, Y, ITYPE = [], [], []
with open(DATA) as f:
    for row in csv.DictReader(f):
        X.append(float(row["easting_km"])); Y.append(float(row["northing_km"]))
        ITYPE.append(row["internal_type"])
X, Y, ITYPE = np.array(X), np.array(Y), np.array(ITYPE)
P_all = np.c_[X, Y]
COMP = "Salzdiapir kompressiv ueberpraegt"
HALOK = np.isin(ITYPE, ["Salzdiapir", "Salzkissen", "Doppelsalinar", COMP])


def box(P, a, b, c, d):
    m = (P[:, 0] >= a) & (P[:, 0] <= b) & (P[:, 1] >= c) & (P[:, 1] <= d)
    return P[m]


def scan():
    rng = np.random.default_rng(101)
    lines = ["", "=" * 60, "NORTH GERMAN BASIN  (BGR InSpEE, %d structures)" % len(X), "=" * 60,
             f"{'subset':38} {'N':>4} {'NND':>5} {'R':>5} {'p':>6} verdict"]
    subsets = [
        ("All structures", P_all),
        ("All halokinetic", P_all[HALOK]),
        ("Central basin halokinetic", box(P_all[HALOK], 350, 650, 5850, 6090)),
        ("Glueckstadt Graben halokinetic", box(P_all[HALOK], 460, 600, 5900, 6070)),
        ("Salt diapirs only (mature)", P_all[ITYPE == "Salzdiapir"]),
        ("Salt pillows only", P_all[ITYPE == "Salzkissen"]),
    ]
    for name, P in subsets:
        R, p, _, med, _ = sc.clark_evans(P, S=499, rng=rng)
        v = "clustered" if R < 0.95 else ("regular" if R > 1.05 else "random")
        lines.append(f"{name:38} {len(P):>4} {med:5.1f} {R:5.2f} {p:6.3f} {v}")
    print("\n".join(lines))
    with open(os.path.join(RES, "results.txt"), "a") as f:
        f.write("\n".join(lines) + "\n")


def type_map():
    col = {"Salzdiapir": PURPLE, COMP: "#D9543A", "Salzkissen": "#D9A23A",
           "Doppelsalinar": "#3AA0D9", "tektonische Strukturen": "#888",
           "Keine Zuordnung": "#ccc", "NA": "#ccc"}
    label = {"Salzdiapir": "salt diapir", COMP: "diapir, compr. overprinted",
             "Salzkissen": "salt pillow", "Doppelsalinar": "double-salinar",
             "tektonische Strukturen": "tectonic", "Keine Zuordnung": "unassigned",
             "NA": "unassigned"}
    fig, ax = plt.subplots(figsize=(9, 6))
    for t in ["Salzdiapir", COMP, "Salzkissen", "Doppelsalinar",
              "tektonische Strukturen", "Keine Zuordnung", "NA"]:
        m = ITYPE == t
        if m.sum():
            ax.scatter(X[m], Y[m], s=14, c=col[t], alpha=0.85, lw=0,
                       label=f"{label[t]} ({m.sum()})")
    ax.set_aspect("equal"); ax.set_xlabel("UTM-E (km)"); ax.set_ylabel("UTM-N (km)")
    ax.set_title("North German salt structures by internal type")
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"ngb_type_map.{ext}"), dpi=130, bbox_inches="tight")
    plt.close(fig)


def hardened():
    P = box(P_all[HALOK], 350, 650, 5850, 6090)
    V, A, path = sc.hull_window(P); n = len(P)
    nnd = np.median(cKDTree(P).query(P, k=2)[0][:, 1])
    rng = np.random.default_rng(21)
    redges = np.arange(0, 46, 2.0); rc, g = sc.pair_correlation(P, redges)
    redL = np.arange(2, 40, 2.0); Lobs = sc.ripley_L_minus_r(P, redL)
    gC = np.array([sc.pair_correlation(sc.csr_in_hull(n, V, path, rng), redges)[1]
                   for _ in range(999)])
    LC = np.array([sc.ripley_L_minus_r(sc.csr_in_hull(n, V, path, rng), redL)
                   for _ in range(999)])
    isamp = sc.inhomogeneous_sampler(P, path, 2 * nnd, rng)
    gI = np.array([sc.pair_correlation(isamp(n), redges)[1] for _ in range(399)])
    R, p_ce, _, _, _ = sc.clark_evans(P, S=999, rng=np.random.default_rng(7))
    mask = rc <= 14
    p_hom = sc.global_envelope_p(g, gC, mask); p_inh = sc.global_envelope_p(g, gI, mask)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    glo, ghi = np.nanpercentile(gC, 2.5, 0), np.nanpercentile(gC, 97.5, 0)
    ilo, ihi = np.nanpercentile(gI, 2.5, 0), np.nanpercentile(gI, 97.5, 0)
    ax[0].fill_between(rc, glo, ghi, color="#bbb", alpha=0.5, label="CSR 95%")
    ax[0].fill_between(rc, ilo, ihi, color="#D9943A", alpha=0.25, label="inhomogeneous 95%")
    ax[0].plot(rc, g, "-o", color=PURPLE, ms=3, lw=2, label="observed")
    ax[0].axhline(1, ls="--", color="k", lw=0.8); ax[0].set_xlim(0, 40)
    ax[0].set_xlabel("r (km)"); ax[0].set_ylabel("g(r)"); ax[0].legend(fontsize=8, frameon=False)
    ax[0].set_title(f"g(r): hard-core + characteristic spacing\n"
                    f"global p(CSR)={p_hom:.3f}, p(inhom)={p_inh:.3f}")
    Llo, Lhi = np.nanpercentile(LC, 2.5, 0), np.nanpercentile(LC, 97.5, 0)
    ax[1].fill_between(redL, Llo, Lhi, color="#bbb", alpha=0.5, label="CSR 95%")
    ax[1].plot(redL, Lobs, "-o", color=PURPLE, ms=3, lw=2, label="observed")
    ax[1].axhline(0, ls="--", color="k", lw=0.8)
    ax[1].set_xlabel("r (km)"); ax[1].set_ylabel("L(r) - r (km)"); ax[1].legend(fontsize=8, frameon=False)
    ax[1].set_title(f"Ripley L: negative = regular (R={R:.2f})")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"ngb_pair_correlation_ripley.{ext}"), dpi=130, bbox_inches="tight")
    plt.close(fig)
    msg = (f"\nNGB central-basin halokinetic (N={n}): R={R:.3f}, p(disp)={p_ce:.4f}; "
           f"g(r) global p(CSR)={p_hom:.4f}, p(inhom)={p_inh:.4f}")
    print(msg)
    with open(os.path.join(RES, "results.txt"), "a") as f:
        f.write(msg + "\n")


def struct_factor():
    P = box(P_all[HALOK], 350, 650, 5850, 6090)
    V, A, path = sc.hull_window(P); n = len(P)
    lo, hi = P.min(0), P.max(0); Lx, Ly = hi - lo
    nnd = np.median(cKDTree(P).query(P, k=2)[0][:, 1])
    rng = np.random.default_rng(31)
    kc, Sobs = sc.structure_factor(P, kmax=1.3)
    isamp = sc.inhomogeneous_sampler(P, path, 2 * nnd, rng)
    SC = np.array([sc.structure_factor(sc.csr_in_hull(n, V, path, rng), 1.3)[1] for _ in range(120)])
    SI = np.array([sc.structure_factor(isamp(n), 1.3)[1] for _ in range(120)])
    clo, chi = np.nanpercentile(SC, 2.5, 0), np.nanpercentile(SC, 97.5, 0)
    ilo, ihi = np.nanpercentile(SI, 2.5, 0), np.nanpercentile(SI, 97.5, 0)
    dipsel = (kc > 0.05) & (kc < 0.35)
    Sdip = np.nanmin(Sobs[dipsel]); kdip = kc[dipsel][np.nanargmin(Sobs[dipsel])]
    Rs = np.array([6, 9, 12, 16, 20, 26, 32])
    v_obs = sc.number_variance(P, Rs, path, rng)
    v_csr = np.mean([sc.number_variance(sc.csr_in_hull(n, V, path, rng), Rs, path, rng) for _ in range(6)], 0)
    v_inh = np.mean([sc.number_variance(isamp(n), Rs, path, rng) for _ in range(6)], 0)
    b1 = np.polyfit(np.log(Rs[2:]), np.log(v_obs[2:]), 1)[0]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].fill_between(kc, clo, chi, color="#bbb", alpha=0.5, label="CSR 95%")
    ax[0].fill_between(kc, ilo, ihi, color="#D9943A", alpha=0.3, label="inhomog. 95% (gradient)")
    ax[0].plot(kc, Sobs, "-o", color=PURPLE, ms=3, lw=2, label="observed S(k)")
    ax[0].axhline(1, ls="--", color="k", lw=0.8); ax[0].set_ylim(0, 4)
    ax[0].annotate("gradient spike\n(lowest mode)", xy=(kc[0], Sobs[0]), xytext=(0.25, 3.0),
                   fontsize=7, arrowprops=dict(arrowstyle="->", lw=0.7))
    ax[0].annotate(f"suppression dip\nS\u2248{Sdip:.2f} (~{2*np.pi/kdip:.0f} km)",
                   xy=(kdip, Sdip), xytext=(0.45, 0.3), fontsize=7,
                   arrowprops=dict(arrowstyle="->", lw=0.7))
    ax[0].set_xlabel("k (1/km)"); ax[0].set_ylabel("S(k)"); ax[0].legend(fontsize=7, frameon=False)
    ax[0].set_title("Direct S(k): intermediate-range suppression,\nk\u21920 masked by basin gradient")
    ax[1].loglog(Rs, v_obs, "-o", color=PURPLE, label=f"observed ~R^{b1:.2f}")
    ax[1].loglog(Rs, v_csr, "--", color="#888", label="CSR ~R^2")
    ax[1].loglog(Rs, v_inh, ":", color="#D9943A", label="inhom Poisson (gradient)")
    ax[1].loglog(Rs, v_obs[0] * (Rs / Rs[0]), "-.", color="#3a7", alpha=0.7, label="hyperuniform ~R^1")
    ax[1].set_xlabel("R (km)"); ax[1].set_ylabel(r"$\sigma^2(R)$"); ax[1].legend(fontsize=7, frameon=False)
    ax[1].set_title("Number variance: suppressed vs gradient-matched\nnull, scaling gradient-dominated")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"ngb_structure_factor.{ext}"), dpi=130, bbox_inches="tight")
    plt.close(fig)
    msg = (f"\nNGB S(k): gradient spike S={Sobs[0]:.2f} @ k={kc[0]:.3f}; "
           f"dip S={Sdip:.2f} @ k={kdip:.2f} (~{2*np.pi/kdip:.0f} km); "
           f"number-variance scaling R^{b1:.2f}")
    print(msg)
    with open(os.path.join(RES, "results.txt"), "a") as f:
        f.write(msg + "\n")


if __name__ == "__main__":
    scan(); type_map(); hardened(); struct_factor()
    print("NGB figures written to ../figures/")
