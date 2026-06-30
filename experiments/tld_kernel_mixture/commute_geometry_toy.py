"""
Commute geometry-vs-kernel check — the structural opposite of the school toy.

Commute differs from school in two decisive ways:
  (1) you CANNOT split the generation budget by sub-purpose — there is one
      "resident workers" budget per origin; there are no level-like sub-budgets to
      hand a far attractor its own production constraint. So the school lever
      (split generation + shared kernel) is UNAVAILABLE.
  (2) the attractors (jobs) are clustered: abundant LOCAL jobs (town centre +
      estates) plus a big FAR employment centre (Belfast). A commuter genuinely
      chooses between near and far, so the local/far split is set by the kernel.

This probes: for the long commute tail (home -> far centre), how much is carried by
attractor MASS/geometry (Belfast's job count — which the model already has from
workplace census) vs the KERNEL tail (willingness to travel) — and how the two
trade off / are identified.

Production constraint: one worker budget distributed over all jobs,
  T_j = budget * a_j f(d_j) / Σ_k a_k f(d_k).
Kernel: Tanner f(t)=u^B exp(B(1-u)), u=t/P, P FIXED (per-component constant; commute's
fixed peak is larger than school's — longer modal commute). Tail carried by BETA.
Time proportional to distance.
"""

import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rng = np.random.default_rng(20260630)

P0 = 8.0     # FIXED commute peak (min) — larger than school's (longer modal commute)


def tanner(t, B, P=P0):
    u = np.asarray(t, float) / P
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.power(u, B) * np.exp(B * (1.0 - u))
    return np.nan_to_num(out, nan=0.0)


def kern(t, spec):
    """spec = ('single', B) or ('mix', w, Bs, Bl) — a 2-exp/2-Tanner sum."""
    if spec[0] == "single":
        return tanner(t, spec[1])
    _, w, Bs, Bl = spec
    return w * tanner(t, Bs) + (1 - w) * tanner(t, Bl)


# ---- clustered job geography (distances from a central home, minutes) ------
def clust(centre, n, sd):
    return np.clip(rng.normal(centre, sd, n), 0.4, None)


local_d = clust(4.0, 300, 2.5)     # town centre + local estates: many near jobs
mid_d = clust(13.0, 70, 2.5)       # Bangor/Comber/Dundonald: moderate
far_d = clust(27.0, 200, 4.0)      # Belfast: big far employment centre
ALL_D = np.concatenate([local_d, mid_d, far_d])
FAR_MASK = np.arange(len(ALL_D)) >= len(local_d) + len(mid_d)
MID_MASK = (np.arange(len(ALL_D)) >= len(local_d)) & ~FAR_MASK


def shares(spec, far_mult=1.0):
    """Production-constrained flow over all jobs; far_mult scales Belfast job MASS."""
    a = np.ones(len(ALL_D))
    a[FAR_MASK] = far_mult                       # geometry knob: far attractor mass
    w = a * kern(ALL_D, spec)
    T = w / w.sum()
    far = T[FAR_MASK].sum()
    mid = T[MID_MASK].sum()
    loc = 1 - far - mid
    med = np.median(np.repeat(ALL_D, (T * 1e6).astype(int)))  # crude weighted median
    return loc, mid, far, med


def far_share(beta, far_mult=1.0):
    return shares(("single", beta), far_mult)[2]


print("=" * 74)
print(f"COMMUTE geometry-vs-kernel  —  FIXED commute peak P = {P0:.0f} min")
print("=" * 74)
print("Jobs: ~300 local (~4min) + 70 mid (~13min) + 200 far/Belfast (~27min).")
print("ONE worker budget (no sub-purpose split possible).\n")

print(f"  {'kernel':<26} {'1/gamma':>8} {'local%':>7} {'mid%':>6} {'far%':>6} {'med':>5}")
KERNELS = [
    ("single sharp  B=2.5", ("single", 2.5), P0 / 2.5),
    ("single moderate B=1.3", ("single", 1.3), P0 / 1.3),
    ("single fat   B=0.6", ("single", 0.6), P0 / 0.6),
    ("2-exp (0.7@2.5 + 0.3@0.4)", ("mix", 0.7, 2.5, 0.4), None),
]
for name, spec, g in KERNELS:
    loc, mid, far, med = shares(spec)
    gs = f"{g:8.1f}" if g else f"{'mix':>8}"
    print(f"  {name:<26} {gs} {loc*100:6.1f} {mid*100:5.1f} {far*100:5.1f} {med:5.1f}")

print("\n-- GEOMETRY lever: far(Belfast) job MASS at fixed moderate kernel (B=1.3) --")
for fm in [0.25, 0.5, 1.0, 1.5, 2.0, 3.0]:
    print(f"   far_mass x{fm:<4} -> far share {far_share(1.3, fm)*100:5.1f}%")

print("\n-- KERNEL lever: BETA at fixed (real) far mass --")
for b in [2.5, 1.8, 1.3, 0.9, 0.6, 0.4]:
    print(f"   BETA {b:<4} (1/g {P0/b:4.1f}) -> far share {far_share(b)*100:5.1f}%")

# ---- mass/kernel degeneracy: which (mass, BETA) give a target far share -----
TARGET = 0.25   # illustrative far-commute share — NOT a measured Newtownards
                # figure; sets where reality sits on the curves. (Note the hard
                # CEILING = the far job-MASS fraction ~35%, reached only as the
                # kernel flattens; a single Tanner tops out ~28%, so a notably
                # higher far-share would itself force a 2-exp/heavier tail.)
print(f"\n-- DEGENERACY: (far_mass, BETA) pairs giving far share = {TARGET:.0%} --")
betas = np.linspace(0.35, 2.2, 40)
iso_mass = []
for b in betas:
    try:
        m = brentq(lambda mm: far_share(b, mm) - TARGET, 0.05, 50.0)
    except ValueError:
        m = np.nan
    iso_mass.append(m)
iso_mass = np.array(iso_mass)
for b, m in list(zip(betas, iso_mass))[::8]:
    print(f"   BETA {b:4.2f} needs far_mass x{m:5.2f}  (same {TARGET:.0%} tail)")
print("  -> a sharp kernel can mimic a fat one by inflating far mass: DEGENERATE")
print("     from commute counts alone. The model PINS far mass from workplace")
print("     census (Belfast jobs are known) -> that vertical cut identifies BETA.")

CENSUS_MASS = 1.0   # the data-known far mass (workplace census) breaks the tie
try:
    beta_id = brentq(lambda b: far_share(b, CENSUS_MASS) - TARGET, 0.3, 4.0)
    print(f"  At census-known far_mass x{CENSUS_MASS:.1f}: identified BETA = "
          f"{beta_id:.2f} (1/gamma = {P0/beta_id:.1f} min)")
except ValueError:
    beta_id = np.nan
    print(f"  Target {TARGET:.0%} not reachable at census mass with a single Tanner "
          f"(needs a fatter/2-exp tail).")

# ---- plot ------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
bins = np.linspace(0, 45, 90)
ctr = 0.5 * (bins[:-1] + bins[1:])

ax = axes[0]
for name, spec, g in KERNELS:
    w = kern(ALL_D, spec)
    T = w / w.sum()
    h, _ = np.histogram(ALL_D, bins=bins, weights=T)
    ax.plot(ctr, h / (h.sum() * (bins[1] - bins[0])), lw=1.8, label=name)
ax.set_xlabel("commute time (min)"); ax.set_ylabel("commute-trip density")
ax.set_title("Commute TLD vs kernel (real job clustering)")
ax.legend(fontsize=7); ax.set_xlim(0, 45)

ax = axes[1]
fm = np.linspace(0.2, 3.0, 30)
ax.plot(fm, [far_share(1.3, m) * 100 for m in fm], "o-", color="C2",
        label="vary far MASS (geometry), B=1.3")
bb = np.linspace(0.4, 2.5, 30)
ax2 = ax.twiny()
ax2.plot(bb, [far_share(b) * 100 for b in bb], "s-", color="C3",
         label="vary BETA (kernel), real mass")
ax.set_xlabel("far job-mass multiplier", color="C2")
ax2.set_xlabel("kernel BETA (sharp→fat ←)", color="C3")
ax.set_ylabel("far(Belfast) commute share %")
ax.set_title("Both levers move the tail")
ax.legend(loc="upper left", fontsize=8); ax2.legend(loc="lower right", fontsize=8)

ax = axes[2]
ax.plot(betas, iso_mass, "k-", lw=2)
ax.axvline(beta_id if np.isfinite(beta_id) else np.nan, color="C0", ls="--",
           label=f"census-pinned mass → BETA≈{beta_id:.2f}" if np.isfinite(beta_id) else "")
ax.axhline(CENSUS_MASS, color="0.6", ls=":")
ax.set_xlabel("kernel BETA"); ax.set_ylabel("far job-mass needed for 30% tail")
ax.set_title("Mass↔kernel degeneracy (broken by workplace census)")
ax.legend(fontsize=8); ax.set_ylim(0, min(8, np.nanmax(iso_mass) * 1.1))

fig.tight_layout()
fig.savefig("commute_geometry.png", dpi=120)
print("\nplot -> commute_geometry.png")
