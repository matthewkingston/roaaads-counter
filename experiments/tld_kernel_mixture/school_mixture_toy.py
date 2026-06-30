"""
School mixture toy — primary/secondary/tertiary WITH real clustering.

Unlike the retail toy (uniform scatter, so all TLD spread came from the kernel),
schools have real spatial structure: primaries dense & local, secondaries sparser,
tertiary RARE and CLUSTERED (a town may have a small local college and/or a far
university in a regional city). This adds back the geometry term n(t) the retail
toy suppressed, so we can test the honest question: for a *splittable* component,
how much of the multi-scale TLD is GEOMETRY (already in the real model via real
attractor locations) vs genuine KERNEL heterogeneity?

We compare three reconstructions of the SAME school flows on the SAME real layout:
  M1  POOLED generation + single kernel  — "lump all school trips, one budget,
      one kernel over all schools."
  M2  SPLIT generation + per-level kernels — the full model (the 'truth').
  M3  SPLIT generation + ONE shared (moderate) kernel — does geometry carry it
      once generation is split by level?

Production constraint: within each generation group g, flow to attractor j is
  T_j = budget_g * f(d_j) / Σ_{k in g} f(d_k)   (each group its OWN denominator).

Kernel: Tanner f(t)=u^B exp(B(1-u)), u=t/P, with P FIXED; the per-level scale is
carried by the tail BETA (1/gamma = P/BETA). Time proportional to distance.
"""

import numpy as np
from scipy.optimize import least_squares
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rng = np.random.default_rng(20260629)

P0 = 2.0           # FIXED peak/rise scale (min)
R = 45.0           # square half-side (min) — edge-free interior
BMAX = 45.0        # TLD analysis range (min)

# trip-rate mix (= generation budgets) and per-level tail BETA
LEVELS = [
    # name        budget  BETA   (1/gamma = P0/BETA)
    ("primary",   0.45,   4.0),   # dense, very local, sharp  (1/g 0.50)
    ("secondary", 0.40,   2.0),   # sparser, short            (1/g 1.00)
    ("tertiary",  0.15,   0.8),   # rare, clustered, fat tail (1/g 2.50)
]


def tanner(t, B, P=P0):
    u = np.asarray(t, float) / P
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.power(u, B) * np.exp(B * (1.0 - u))
    return np.nan_to_num(out, nan=0.0)


def scatter_dist(n):
    """Uniform 2D scatter in the square -> distance from the central house."""
    return np.hypot(rng.uniform(-R, R, n), rng.uniform(-R, R, n))


def cluster_dist(centre, n, sigma=2.0):
    """A tight blob of n attractors around a point at radius `centre` (min)."""
    cx, cy = centre, 0.0
    return np.hypot(cx + rng.normal(0, sigma, n), cy + rng.normal(0, sigma, n))


# ---- real attractor layout -------------------------------------------------
prim_d = scatter_dist(600)               # primaries everywhere, dense
sec_d = scatter_dist(150)                # secondaries sparser
tert_local = cluster_dist(3.0, 2)        # small local FE college (~3 min)
tert_far = cluster_dist(33.0, 4)         # university in a regional city (~33 min)

bins = np.linspace(0, BMAX, 120)
centres = 0.5 * (bins[:-1] + bins[1:])
binw = bins[1] - bins[0]


def flow_tld(groups):
    """groups = list of (distances, budget, beta). Production-constrained flow
    per group (own denominator), binned to a normalised TLD density."""
    hist = np.zeros(len(centres))
    for d, budget, beta in groups:
        w = tanner(d, beta)
        s = w.sum()
        if s <= 0:
            continue
        T = budget * w / s                      # production constraint
        h, _ = np.histogram(d, bins=bins, weights=T)
        hist += h
    dens = hist / (hist.sum() * binw)           # normalise to a density
    return dens


def stats(dens):
    cdf = np.cumsum(dens) * binw
    cdf /= cdf[-1]
    med = np.interp(0.5, cdf, centres)
    p90 = np.interp(0.9, cdf, centres)
    far = dens[centres > 20].sum() * binw       # share of trips > 20 min
    return med, p90, far


def rms(a, b):
    return np.sqrt(np.mean((a - b) ** 2)) / max(b.max(), 1e-9) * 100


def run_scenario(tert_d, label):
    # M2 truth: split generation + per-level kernels
    truth = flow_tld([
        (prim_d, 0.45, 4.0),
        (sec_d, 0.40, 2.0),
        (tert_d, 0.15, 0.8),
    ])

    # M1 pooled: one budget over ALL schools, single kernel (fit BETA, P fixed)
    all_d = np.concatenate([prim_d, sec_d, tert_d])

    def resid_pooled(b):
        return flow_tld([(all_d, 1.0, b[0])]) - truth
    bp = least_squares(resid_pooled, [2.0], bounds=([0.1], [8.0])).x[0]
    m1 = flow_tld([(all_d, 1.0, bp)])

    # M3 split + shared kernel: per-level budgets, ONE shared BETA (fit)
    def resid_shared(b):
        return flow_tld([(prim_d, 0.45, b[0]), (sec_d, 0.40, b[0]),
                         (tert_d, 0.15, b[0])]) - truth
    bs = least_squares(resid_shared, [2.0], bounds=([0.1], [8.0])).x[0]
    m3 = flow_tld([(prim_d, 0.45, bs), (sec_d, 0.40, bs), (tert_d, 0.15, bs)])

    print(f"\n=== scenario: {label} ===")
    print(f"  tertiary attractors at distances (min): "
          f"{np.sort(tert_d).round(1).tolist()}")
    for name, d, tag in [("TRUTH (split+per-level)", truth, ""),
                         (f"M1 pooled+single (β={bp:.2f})", m1, "pool"),
                         (f"M3 split+shared (β={bs:.2f})", m3, "shared")]:
        med, p90, far = stats(d)
        extra = f"   RMS vs truth {rms(d, truth):5.2f}%" if tag else ""
        print(f"  {name:<30} median {med:5.2f}  p90 {p90:5.2f}  "
              f">20min {far*100:4.1f}%{extra}")
    return dict(label=label, tert_d=tert_d, truth=truth, m1=m1, m3=m3,
                bp=bp, bs=bs)


print("=" * 76)
print(f"SCHOOL MIXTURE TOY  —  real clustering, FIXED P = {P0:.1f} min")
print("=" * 76)
print("Primaries dense+local, secondaries sparser, tertiary clustered.")
print("M1 pooled generation | M2 split+per-level (truth) | M3 split+shared kernel.")

s_far = run_scenario(tert_far, "tertiary = far university only")
s_both = run_scenario(np.concatenate([tert_local, tert_far]),
                      "tertiary = local college + far university")


# ---- plot ------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(14, 5.2))
for ax, s in zip(axes, [s_far, s_both]):
    ax.plot(centres, s["truth"], "k-", lw=2.6, label="truth (split + per-level)")
    ax.plot(centres, s["m1"], "r-", lw=2.0,
            label=f"M1 pooled + single (β={s['bp']:.2f})")
    ax.plot(centres, s["m3"], "b--", lw=2.0,
            label=f"M3 split + shared (β={s['bs']:.2f})")
    ax.set_xlabel("trip time (min)")
    ax.set_ylabel("school-trip density")
    ax.set_title(s["label"])
    ax.legend(fontsize=8)
    ax.set_xlim(0, BMAX)
fig.suptitle("Geometry vs kernel for a splittable component: "
             "split generation carries the tail; pooling loses it", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig("school_mixture.png", dpi=120)
print("\nplot -> school_mixture.png")
