# TLD kernel mixture — thought experiments (NOT part of the pipeline)

Standalone, self-contained exploratory scripts. They read **nothing** from the
repo and touch no model state — pure numpy/scipy/matplotlib. They explore how a
**deterrence kernel** (Tanner / exponential-tail) behaves when the demand it
serves is really a *mixture* of trip-length scales (the "retail bucket" problem),
and what that means for **transferring** a calibrated kernel to a new town.

Run from this directory:

```
python3 tld_mixture_experiment.py   # -> tld_mixture.png
python3 tld_mix_sweep.py            # -> tld_mix_sweep.png
```

## Key identity used throughout

Observed trip-length distribution `TLD(t) ∝ n(t)·f(t)`, where `f` is the kernel
and `n(t)` is the opportunity density at travel-time `t`. Uniform 2D scatter ⇒
`n(t) ∝ t`. So **matching a TLD is not setting `f = TLD`** — you deconvolve the
geometry. A large square is used so `n(t) ∝ t` holds over the kernel support and
the *mixture* effect is isolated from finite-domain edge clipping.

`P` (the Tanner peak) is held **fixed** as a structural constant; the tail shape
`BETA` (decay `1/γ = P/BETA`) carries the per-trip-type scale — simulating a move
toward a non-Tanner kernel that keeps an exponential tail.

## `tld_mixture_experiment.py`

A house at the centre of a square; businesses scattered uniformly; 4 retail
sub-types at rates 0.50/0.30/0.15/0.05, each with its own tail `BETA`
(production-constrained: each type gets a fixed trip share). Fits a **single**
fixed-P Tanner and a **2-Tanner (short+long)** to the aggregate TLD.

Finding: a single Tanner can't hold the multi-scale tail (under-predicts long
trips); a 2-regime fit drops the residual ~13× and matches the aggregate p90
almost exactly. Short/long is essentially enough.

## `tld_mix_sweep.py`

Transportability test. Tilts the trip-rate mix from local-heavy to
comparison-heavy (long fraction 6%→70%) and asks whether the 2-regime kernel
transports: (A) do the two tail `BETA`s stay stable while only the split weight
`w` moves; (B) freeze both tails at the baseline-mix values and re-fit **only**
`w` per mix — how well does that hold.

Finding: 2 regimes fit to <0.7% RMS at *every* mix. Carrying fixed tails and
re-fitting only `w` holds <1% RMS within ~±15pp of the calibration mix. `w` is
the smooth, interpretable knob — and under a production constraint it is set by
**generation**, not the kernel. ⇒ split retail into 2 generation sub-components,
each with a transported tail, and the mix-transportability problem dissolves into
generation.

See the project memory note `project_tld_retail_mixture` for the full conclusion.
