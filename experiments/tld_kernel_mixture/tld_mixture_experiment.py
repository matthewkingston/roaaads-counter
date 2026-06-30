"""
TLD mixture thought experiment — fixed-P kernels, tail (BETA) carries the scale.

Context
-------
We are moving toward a kernel that is NOT pure Tanner but still has an
exponential tail, with the peak/rise scale P held FIXED as a structural
constant. So here P is fixed everywhere and the per-trip-type differences are
carried entirely by the tail shape BETA (characteristic decay 1/gamma = P/BETA):
short frequent trips = sharp tail (high BETA), long infrequent trips = fat tail
(low BETA), all sharing the same fixed peak P.

Question
--------
"Retail" is a bucket of sub-trip-types with different tail scales. If we tune a
kernel so the model aggregate trip-length distribution (TLD) matches the observed
aggregate TLD, with P fixed:
 (1) what effective single-tail kernel do we get, and how bad is the fit?
 (2) does splitting into just TWO regimes (short + long) — a 2-component fit —
     essentially recover the mixture? i.e. is short/long enough?

Setup (uniform random scatter; time proportional to distance)
-------------------------------------------------------------
 - House at the centre of a large square measured in minutes of travel time.
 - Businesses scattered uniformly; 4 types at map fractions 0.50/0.30/0.15/0.05.
 - Each type: Tanner kernel f_c(t)=u^B exp(B(1-u)), u=t/P, with the SAME fixed P
   and its own BETA. Production-constrained: each type gets a fixed trip share.

Identity: observed TLD(t) ∝ n(t) f(t); uniform 2D => n(t) ∝ t, so each type's
TLD is t*f_c(t) and the aggregate is the rate-weighted sum. Large square => the
n(t) ∝ t geometry holds over the kernel support (edge effects isolated out).
"""

import numpy as np
from scipy.optimize import least_squares
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rng = np.random.default_rng(20260629)

# ---- fixed peak + per-type tails ------------------------------------------
P0 = 3.0   # FIXED peak/rise scale (min) — structural constant of the kernel family
# frac = trip rate (= map frequency); BETA = tail shape at the shared fixed P.
# 1/gamma = P0/BETA is the characteristic decay; spread of scales via BETA only.
#               name           frac   BETA
TYPES = [
    ("convenience", 0.50, 4.0),   # 1/g=0.75 min : short, sharp
    ("weekly_shop", 0.30, 2.0),   # 1/g=1.50
    ("comparison",  0.15, 1.0),   # 1/g=3.00
    ("specialist",  0.05, 0.5),   # 1/g=6.00 : long reach, fat (but real) tail
]

R = 45.0            # square half-side, minutes (large => edge-free interior)
T_ANALYSIS = 30.0   # restrict all fitting/plotting to the edge-free region
N_BUS = 250_000
N_TRIPS = 400_000


def tanner(t, B, P=P0):
    """Tanner deterrence kernel at fixed peak P: peak f(P)=1, f(0)=0, tail ~exp(-B t/P)."""
    u = np.asarray(t, float) / P
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.power(u, B) * np.exp(B * (1.0 - u))
    return np.nan_to_num(out, nan=0.0)


def tld_shape(t, B, P=P0):
    """Unnormalised continuum TLD shape t*f(t) for uniform 2D scatter."""
    return t * tanner(t, B, P)


t_grid = np.linspace(1e-3, T_ANALYSIS, 4000)
dt = t_grid[1] - t_grid[0]


def norm_density(y, x=t_grid):
    area = np.trapz(y, x)
    return y / area if area > 0 else y


# ---- continuum aggregate TLD (clean target) -------------------------------
agg = np.zeros_like(t_grid)
comp_densities = []
for name, frac, B in TYPES:
    d_c = norm_density(tld_shape(t_grid, B))
    comp_densities.append(d_c)
    agg += frac * d_c
agg = norm_density(agg)


def stats_of(density, x=t_grid):
    cdf = np.cumsum(density) * (x[1] - x[0])
    cdf /= cdf[-1]
    median = np.interp(0.5, cdf, x)
    p90 = np.interp(0.9, cdf, x)
    mean = np.trapz(density * x, x) / np.trapz(density, x)
    return mean, median, p90


# ---- (1) single fixed-P Tanner fit (only BETA free) -----------------------
def resid_1(theta):
    (B,) = theta
    return norm_density(tld_shape(t_grid, B)) - agg


fit1 = least_squares(resid_1, x0=[2.0], bounds=([0.2], [8.0]))
B1 = fit1.x[0]
fit1_density = norm_density(tld_shape(t_grid, B1))
rms1 = np.sqrt(np.mean((fit1_density - agg) ** 2))


# ---- (2) two fixed-P Tanners: short + long regime -------------------------
# density = w*norm(t f(B_s)) + (1-w)*norm(t f(B_l)); P fixed for both.
def two_density(w, Bs, Bl):
    return w * norm_density(tld_shape(t_grid, Bs)) + \
           (1.0 - w) * norm_density(tld_shape(t_grid, Bl))


def resid_2(theta):
    w, Bs, Bl = theta
    return two_density(w, Bs, Bl) - agg


# multi-start to dodge local minima / label-swap
best = None
for w0 in (0.5, 0.7, 0.85):
    for Bs0, Bl0 in ((4.0, 0.8), (3.0, 1.0), (5.0, 0.5)):
        f = least_squares(resid_2, x0=[w0, Bs0, Bl0],
                          bounds=([0.0, 0.5, 0.1], [1.0, 8.0, 8.0]))
        if best is None or f.cost < best.cost:
            best = f
w2, Bs2, Bl2 = best.x
if Bs2 < Bl2:                      # order short (sharp) first
    w2, Bs2, Bl2 = 1 - w2, Bl2, Bs2
fit2_density = two_density(w2, Bs2, Bl2)
rms2 = np.sqrt(np.mean((fit2_density - agg) ** 2))


# ---- Monte-Carlo realisation (faithful random scatter) --------------------
def scatter_radii(n):
    return np.hypot(rng.uniform(-R, R, n), rng.uniform(-R, R, n))


mc_t = []
for name, frac, B in TYPES:
    n_bus_c = int(N_BUS * frac)
    n_trip_c = int(N_TRIPS * frac)
    t_bus = scatter_radii(n_bus_c)
    w = tanner(t_bus, B)
    w /= w.sum()
    mc_t.append(t_bus[rng.choice(n_bus_c, size=n_trip_c, p=w)])
mc_t = np.concatenate(mc_t)
mc_t = mc_t[mc_t < T_ANALYSIS]
bins = np.linspace(0, T_ANALYSIS, 80)
mc_hist, edges = np.histogram(mc_t, bins=bins, density=True)
centres = 0.5 * (edges[:-1] + edges[1:])


# ---- report ----------------------------------------------------------------
def g(B):
    return P0 / B


print("=" * 76)
print(f"TLD MIXTURE  —  FIXED P = {P0:.1f} min;  tail BETA carries the scale")
print("=" * 76)
print(f"\nWorld {2*R:.0f}x{2*R:.0f} min, {N_BUS:,} businesses, {N_TRIPS:,} trips; "
      f"analysis t<{T_ANALYSIS:.0f} min.\n")
print("Ground-truth components (mixed reality):")
print(f"  {'type':<12} {'rate':>5} {'BETA':>5} {'1/gamma':>8} {'median':>7} {'p90':>6}")
for (name, frac, B), d_c in zip(TYPES, comp_densities):
    m, med, p90 = stats_of(d_c)
    print(f"  {name:<12} {frac:>5.2f} {B:>5.2f} {g(B):>8.2f} {med:>7.2f} {p90:>6.2f}")
am, amed, ap90 = stats_of(agg)
print(f"  {'AGGREGATE':<12} {'':>5} {'':>5} {'':>8} {amed:>7.2f} {ap90:>6.2f}")

print("\n" + "-" * 76)
print("(1) SINGLE fixed-P Tanner fit (BETA only):")
m1, med1, p901 = stats_of(fit1_density)
print(f"    BETA_fit = {B1:.3f}   1/gamma = {g(B1):.2f} min")
print(f"    median {med1:.2f}  p90 {p901:.2f} min   (truth p90 {ap90:.2f})")
print(f"    residual RMS = {rms1/agg.max()*100:.2f}% of peak")
minB = min(b for _, _, b in TYPES)
print(f"    -> effective BETA {B1:.2f} vs flattest component {minB:.2f}: "
      f"{'FLATTER than any component' if B1 < minB else 'within component range'}")

print("\n(2) TWO fixed-P Tanners (short + long regime):")
ms, meds, p90s = stats_of(norm_density(tld_shape(t_grid, Bs2)))
ml, medl, p90l = stats_of(norm_density(tld_shape(t_grid, Bl2)))
print(f"    short: weight {w2:5.2f}  BETA {Bs2:5.2f}  1/gamma {g(Bs2):5.2f}  "
      f"median {meds:5.2f}  p90 {p90s:5.2f}")
print(f"    long : weight {1-w2:5.2f}  BETA {Bl2:5.2f}  1/gamma {g(Bl2):5.2f}  "
      f"median {medl:5.2f}  p90 {p90l:5.2f}")
m2, med2, p902 = stats_of(fit2_density)
print(f"    combined median {med2:.2f}  p90 {p902:.2f} min   (truth p90 {ap90:.2f})")
print(f"    residual RMS = {rms2/agg.max()*100:.2f}% of peak")
print(f"\n    residual drop 1->2 Tanner: {rms1/agg.max()*100:.2f}%  ->  "
      f"{rms2/agg.max()*100:.2f}%   ({rms1/rms2:.1f}x better)")
print("=" * 76)


# ---- plot ------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

ax = axes[0]
ax.bar(centres, mc_hist, width=(centres[1] - centres[0]) * 0.9, color="0.85",
       label="Monte-Carlo scatter (observed)")
ax.plot(t_grid, agg, "k-", lw=2.6, label="aggregate TLD (truth)")
ax.plot(t_grid, fit1_density, "r-", lw=2.0, label=f"1 Tanner (β={B1:.2f})")
ax.plot(t_grid, fit2_density, "b-", lw=2.0,
        label=f"2 Tanner (β={Bs2:.2f}/{Bl2:.2f})")
ax.plot(t_grid, w2 * norm_density(tld_shape(t_grid, Bs2)), "b--", lw=1.0, alpha=0.7)
ax.plot(t_grid, (1 - w2) * norm_density(tld_shape(t_grid, Bl2)), "b:", lw=1.2, alpha=0.7)
ax.set_xlabel("trip time (min)"); ax.set_ylabel("trip-length density")
ax.set_title(f"Fixed P={P0:.0f} min: 1- vs 2-Tanner fit")
ax.legend(fontsize=8); ax.set_xlim(0, T_ANALYSIS)

ax = axes[1]
ax.semilogy(t_grid, agg, "k-", lw=2.6, label="truth")
ax.semilogy(t_grid, fit1_density, "r-", lw=2.0, label="1 Tanner")
ax.semilogy(t_grid, fit2_density, "b-", lw=2.0, label="2 Tanner")
for (name, frac, B), d_c in zip(TYPES, comp_densities):
    ax.semilogy(t_grid, frac * d_c, lw=0.9, ls="--", alpha=0.5)
ax.set_xlabel("trip time (min)"); ax.set_ylabel("density (log)")
ax.set_title("Log scale — the tail")
ax.set_ylim(agg.max() * 1e-4, agg.max() * 1.6); ax.set_xlim(0, T_ANALYSIS)
ax.legend(fontsize=8)

fig.tight_layout()
fig.savefig("tld_mixture.png", dpi=120)
print("\nplot -> tld_mixture.png")
