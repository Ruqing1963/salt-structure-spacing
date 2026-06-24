"""
stats_common.py
Shared spatial point-pattern and structure-factor statistics for the
salt-structure-spacing analysis.

Pure NumPy/SciPy/Matplotlib; no GIS dependencies.

Estimators
----------
- Clark-Evans nearest-neighbour index with Monte-Carlo CSR null (edge-bias removed)
- pair correlation g(r) with the border / reduced-sample edge correction
- Ripley's L(r) - r with the border correction
- direct structure factor S(k) at allowed wavevectors, radially averaged
- number variance sigma^2(R)
- CSR and inhomogeneous-Poisson null samplers
- global studentised maximum-deviation envelope test
"""
import numpy as np
from scipy.spatial import cKDTree, ConvexHull
from matplotlib.path import Path as MplPath

EARTH_R_KM = 6371.0


# ---------------------------------------------------------------- projection
def project_lonlat_km(lon, lat):
    """Equirectangular projection (km) about the mean latitude. For regional
    extents this preserves distances to well within a per-cent."""
    lon = np.asarray(lon, float); lat = np.asarray(lat, float)
    lon0, lat0 = lon.mean(), lat.mean()
    x = np.radians(lon - lon0) * EARTH_R_KM * np.cos(np.radians(lat0))
    y = np.radians(lat - lat0) * EARTH_R_KM
    return np.c_[x, y]


# ------------------------------------------------------------------- windows
def hull_window(P):
    h = ConvexHull(P)
    return P[h.vertices], h.volume, MplPath(P[h.vertices])


def dist_to_boundary(P, V):
    """Minimum distance from each point to the polygon boundary defined by
    hull vertices V (used by the border edge correction)."""
    n = len(V); b = np.full(len(P), np.inf)
    for k in range(n):
        a = V[k]; c = V[(k + 1) % n]; ac = c - a; L2 = ac @ ac
        t = np.clip(((P - a) @ ac) / L2, 0, 1)
        proj = a + t[:, None] * ac
        b = np.minimum(b, np.sqrt(((P - proj) ** 2).sum(1)))
    return b


# ------------------------------------------------------------- null samplers
def csr_in_hull(n, V, path, rng):
    lo, hi = V.min(0), V.max(0); out = []
    while len(out) < n:
        q = rng.uniform(lo, hi, (n * 3, 2)); q = q[path.contains_points(q)]
        out.extend(q.tolist())
    return np.array(out[:n])


def inhomogeneous_sampler(P, path, bw, rng):
    """Sampler for an inhomogeneous-Poisson null whose intensity is a Gaussian
    kernel estimate of the observed intensity (bandwidth bw). Carries the
    basin density gradient into the null."""
    def sample(n):
        out = []
        while len(out) < n:
            q = P[rng.integers(0, len(P), n * 3)] + rng.normal(0, bw, (n * 3, 2))
            q = q[path.contains_points(q)]; out.extend(q.tolist())
        return np.array(out[:n])
    return sample


# --------------------------------------------------------------- estimators
def clark_evans(P, S=499, rng=None):
    """Clark-Evans R = mean NND / (MC-CSR mean NND). R>1 dispersed, <1 clustered.
    Returns (R, p_one_sided_dispersion, mean_nnd, median_nnd, hull_area)."""
    if rng is None:
        rng = np.random.default_rng(0)
    nn = cKDTree(P).query(P, k=2)[0][:, 1]; obs = nn.mean()
    V, A, path = hull_window(P)
    nulls = np.empty(S)
    for i in range(S):
        Q = csr_in_hull(len(P), V, path, rng)
        nulls[i] = cKDTree(Q).query(Q, k=2)[0][:, 1].mean()
    R = obs / nulls.mean()
    p = (1 + np.sum(nulls >= obs)) / (S + 1)
    return R, p, obs, np.median(nn), A


def pair_correlation(P, redges, interior_only=True):
    """g(r) with the border (reduced-sample) edge correction."""
    V, A, _ = hull_window(P)
    D = np.sqrt(((P[:, None] - P[None]) ** 2).sum(2)); np.fill_diagonal(D, np.inf)
    b = dist_to_boundary(P, V); rho = len(P) / A
    rc = 0.5 * (redges[1:] + redges[:-1]); g = np.full(len(rc), np.nan)
    for i in range(len(rc)):
        keep = b > redges[i + 1]
        if keep.sum() < 3:
            continue
        ring = ((D[keep] > redges[i]) & (D[keep] <= redges[i + 1])).sum()
        area = np.pi * (redges[i + 1] ** 2 - redges[i] ** 2)
        g[i] = ring / (keep.sum() * rho * area)
    return rc, g


def ripley_L_minus_r(P, redges):
    """L(r) - r with the border correction. Negative => regular/dispersed."""
    V, A, _ = hull_window(P)
    D = np.sqrt(((P[:, None] - P[None]) ** 2).sum(2)); np.fill_diagonal(D, np.inf)
    b = dist_to_boundary(P, V); n = len(P); K = np.full(len(redges), np.nan)
    for ii, r in enumerate(redges):
        keep = b > r
        if keep.sum() < 3:
            continue
        K[ii] = A * (D[keep] <= r).sum() / (n * keep.sum())
    return np.sqrt(np.clip(K, 0, None) / np.pi) - redges


def structure_factor(P, kmax, mmax=60, nbin=28):
    """Direct radially-averaged S(k) at the wavevectors allowed by the
    rectangular bounding window."""
    q = P - P.mean(0)
    Lx, Ly = (P.max(0) - P.min(0))
    ms = np.arange(-mmax, mmax + 1)
    KX, KY = np.meshgrid(2 * np.pi * ms / Lx, 2 * np.pi * ms / Ly)
    K = np.c_[KX.ravel(), KY.ravel()]; kk = np.sqrt((K ** 2).sum(1))
    sel = (kk > 0) & (kk <= kmax); K = K[sel]; kk = kk[sel]
    S = (np.abs(np.exp(-1j * (K @ q.T)).sum(1)) ** 2) / len(q)
    bins = np.linspace(0, kmax, nbin); idx = np.digitize(kk, bins)
    kc = 0.5 * (bins[1:] + bins[:-1])
    Sm = np.array([S[idx == i].mean() if np.any(idx == i) else np.nan
                   for i in range(1, len(bins))])
    return kc, Sm


def number_variance(P, Rs, path, rng, M=2500, ncenters=1500):
    """Number variance sigma^2(R) over interior discs of radius R."""
    lo, hi = P.min(0), P.max(0); t = cKDTree(P); out = []
    for R in Rs:
        c = rng.uniform(lo + R, hi - R, (M, 2)); c = c[path.contains_points(c)]
        cnt = np.array([len(t.query_ball_point(cc, R)) for cc in c[:ncenters]])
        out.append(cnt.var())
    return np.array(out)


def global_envelope_p(obs, sims, mask=None):
    """Global studentised maximum-deviation (MAD) envelope test.
    obs: (nr,), sims: (S, nr). Returns p-value."""
    if mask is None:
        mask = np.ones(obs.shape, bool)
    m = np.nanmean(sims, 0); s = np.nanstd(sims, 0) + 1e-9
    Dobs = np.nanmax(np.abs(obs - m)[mask] / s[mask])
    Ds = np.nanmax(np.abs(sims - m)[:, mask] / s[mask], axis=1)
    return (1 + np.sum(Ds >= Dobs)) / (len(sims) + 1)
