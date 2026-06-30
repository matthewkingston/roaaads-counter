"""
Mix-sweep transportability test for the 2-regime (short+long) retail kernel.

Fixed P, fixed 4 ground-truth component tails (BETA 4/2/1/0.5). We tilt the
trip-rate mix from local-heavy to comparison-heavy and ask:

  (A) Free 2-Tanner fit: do the two regime tails (BETA_short, BETA_long) stay
      stable while only the split weight w moves?  -> if so, the tails are the
      transportable primitives and w is the per-town knob.
  (B) TRANSPORT TEST: freeze BETA_short/BETA_long at their baseline-mix values
      and re-fit ONLY w for each new mix. How well does that fit the aggregate
      TLD across the whole sweep?  -> this is exactly "carry the two kernels
      between towns, recalibrate only the local/comparison split fraction."
  Compared against the single-Tanner (one BETA) baseline.

x-axis = the long fraction (comparison+specialist share of trips).
"""

import numpy as np
from scipy.optimize import least_squares
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

P0 = 3.0
BASE_RATES = np.array([0.50, 0.30, 0.15, 0.05])   # conv, weekly, comparison, specialist
BETAS = np.array([4.0, 2.0, 1.0, 0.5])            # fixed component tails
T_ANALYSIS = 30.0
t_grid = np.linspace(1e-3, T_ANALYSIS, 4000)


def tanner(t, B, P=P0):
    u = np.asarray(t, float) / P
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.power(u, B) * np.exp(B * (1.0 - u))
    return np.nan_to_num(out, nan=0.0)


def tld_shape(t, B):
    return t * tanner(t, B)


def norm_density(y):
    a = np.trapz(y, t_grid)
    return y / a if a > 0 else y


COMP = [norm_density(tld_shape(t_grid, B)) for B in BETAS]   # normalised type TLDs


def aggregate(rates):
    agg = np.zeros_like(t_grid)
    for r, d in zip(rates, COMP):
        agg += r * d
    return norm_density(agg)


def p90_of(density):
    cdf = np.cumsum(density) * (t_grid[1] - t_grid[0])
    cdf /= cdf[-1]
    return np.interp(0.9, cdf, t_grid)


def tilt_rates(theta):
    """Geometric tilt: rate_k ∝ base_k * exp(theta*k). theta>0 -> longer mix."""
    k = np.arange(4)
    w = BASE_RATES * np.exp(theta * k)
    return w / w.sum()


def two_density(w, Bs, Bl):
    return w * norm_density(tld_shape(t_grid, Bs)) + \
           (1 - w) * norm_density(tld_shape(t_grid, Bl))


def fit_two_free(agg):
    best = None
    for w0 in (0.5, 0.7, 0.85):
        for Bs0, Bl0 in ((4.0, 0.8), (3.0, 1.0), (5.0, 0.5)):
            f = least_squares(lambda th: two_density(*th) - agg,
                              x0=[w0, Bs0, Bl0],
                              bounds=([0.0, 0.5, 0.1], [1.0, 8.0, 8.0]))
            if best is None or f.cost < best.cost:
                best = f
    w, Bs, Bl = best.x
    if Bs < Bl:
        w, Bs, Bl = 1 - w, Bl, Bs
    return w, Bs, Bl


def fit_one(agg):
    f = least_squares(lambda th: norm_density(tld_shape(t_grid, th[0])) - agg,
                      x0=[2.0], bounds=([0.2], [8.0]))
    return f.x[0]


def rms(density, agg):
    return np.sqrt(np.mean((density - agg) ** 2)) / agg.max() * 100


# ---- baseline (theta=0) frozen tails for the transport test ---------------
agg0 = aggregate(tilt_rates(0.0))
w0, Bs_fix, Bl_fix = fit_two_free(agg0)
print(f"baseline frozen tails:  BETA_short={Bs_fix:.3f}  BETA_long={Bl_fix:.3f}  "
      f"(w={w0:.3f})\n")


# ---- sweep -----------------------------------------------------------------
thetas = np.linspace(-0.7, 1.2, 20)
rows = []
for th in thetas:
    rates = tilt_rates(th)
    long_frac = rates[2] + rates[3]
    agg = aggregate(rates)

    w_f, Bs_f, Bl_f = fit_two_free(agg)
    rms2_free = rms(two_density(w_f, Bs_f, Bl_f), agg)

    # transport test: tails frozen at baseline, refit only w
    fw = least_squares(lambda th_: two_density(th_[0], Bs_fix, Bl_fix) - agg,
                       x0=[0.6], bounds=([0.0], [1.0]))
    w_t = fw.x[0]
    rms2_frozen = rms(two_density(w_t, Bs_fix, Bl_fix), agg)

    B1 = fit_one(agg)
    rms1 = rms(norm_density(tld_shape(t_grid, B1)), agg)

    rows.append(dict(theta=th, long_frac=long_frac, p90=p90_of(agg),
                     w_free=w_f, Bs=Bs_f, Bl=Bl_f, rms2_free=rms2_free,
                     w_frozen=w_t, rms2_frozen=rms2_frozen, B1=B1, rms1=rms1))

lf = np.array([r["long_frac"] for r in rows])

print(f"{'long%':>6} {'p90':>5} | {'free: w':>7} {'Bs':>5} {'Bl':>5} {'rms%':>5} "
      f"| {'froz w':>6} {'rms%':>5} | {'1T B':>5} {'rms%':>5}")
for r in rows:
    print(f"{r['long_frac']*100:6.1f} {r['p90']:5.1f} | "
          f"{r['w_free']:7.3f} {r['Bs']:5.2f} {r['Bl']:5.2f} {r['rms2_free']:5.2f} | "
          f"{r['w_frozen']:6.3f} {r['rms2_frozen']:5.2f} | "
          f"{r['B1']:5.2f} {r['rms1']:5.2f}")

Bs_arr = np.array([r["Bs"] for r in rows])
Bl_arr = np.array([r["Bl"] for r in rows])
print(f"\nfree tail stability over the sweep:")
print(f"  BETA_short: {Bs_arr.min():.2f}–{Bs_arr.max():.2f} "
      f"(spread {Bs_arr.max()-Bs_arr.min():.2f})")
print(f"  BETA_long : {Bl_arr.min():.2f}–{Bl_arr.max():.2f} "
      f"(spread {Bl_arr.max()-Bl_arr.min():.2f})")
print(f"  w_free    : {min(r['w_free'] for r in rows):.2f}–"
      f"{max(r['w_free'] for r in rows):.2f}")
print(f"\ntransport-test (frozen tails, refit w only) RMS: "
      f"{min(r['rms2_frozen'] for r in rows):.2f}–"
      f"{max(r['rms2_frozen'] for r in rows):.2f}%")
print(f"single-Tanner RMS over sweep: "
      f"{min(r['rms1'] for r in rows):.2f}–{max(r['rms1'] for r in rows):.2f}%")


# ---- plot ------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))

ax = axes[0]
ax.plot(lf * 100, Bs_arr, "o-", color="C0", label="β_short (free)")
ax.plot(lf * 100, Bl_arr, "s-", color="C3", label="β_long (free)")
ax.plot(lf * 100, [r["B1"] for r in rows], "^--", color="0.5",
        label="β single-Tanner")
ax.axhline(Bs_fix, color="C0", ls=":", alpha=0.6)
ax.axhline(Bl_fix, color="C3", ls=":", alpha=0.6)
ax.set_xlabel("long fraction (comparison+specialist) %")
ax.set_ylabel("fitted BETA"); ax.set_title("Tail stability vs mix")
ax.legend(fontsize=8)

ax = axes[1]
ax.plot(lf * 100, [r["w_free"] for r in rows], "o-", color="C0",
        label="w_short (free fit)")
ax.plot(lf * 100, [r["w_frozen"] for r in rows], "x--", color="C2",
        label="w_short (frozen tails)")
ax.set_xlabel("long fraction %"); ax.set_ylabel("short-regime weight w")
ax.set_title("The split weight is the moving knob"); ax.legend(fontsize=8)

ax = axes[2]
ax.plot(lf * 100, [r["rms1"] for r in rows], "^--", color="0.5",
        label="1 Tanner")
ax.plot(lf * 100, [r["rms2_frozen"] for r in rows], "x-", color="C2",
        label="2 Tanner, frozen tails (transport)")
ax.plot(lf * 100, [r["rms2_free"] for r in rows], "o-", color="C0",
        label="2 Tanner, free")
ax.set_xlabel("long fraction %"); ax.set_ylabel("residual RMS (% of peak)")
ax.set_title("Fit quality across the sweep"); ax.legend(fontsize=8)
ax.set_ylim(bottom=0)

fig.tight_layout()
fig.savefig("tld_mix_sweep.png", dpi=120)
print("\nplot -> tld_mix_sweep.png")
