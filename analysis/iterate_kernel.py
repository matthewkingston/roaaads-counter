"""Production-constraint (1/D_i) iteration → constrained double-exp willingness kernels.

The v1 kernels (analysis/kernel_fit.json) were recovered from the UNCONSTRAINED n_Ire(t), which
inflates the tail tau_l: the deployed model is production-constrained (T_ij = P_i A_j f / D_i,
D_i = Σ_k A_k f(c_ik)), so the correct geometry is Ñ(t) = Σ_{i,j} (P_i/D_i) A_j δ(c_ij−t). D_i
depends on f — a fixed point. This iterates it and writes a new set of iterated kernel parameters.

    f^(0) = kernel_fit.json double-exps  →  D_i[f]  →  Ñ[f]  →  f = fit_double([TLD/Ñ]/driveshare)  →  repeat

Six independent kernels (the three school levels each their own — NOT shared). **Artifact only** —
writes analysis/kernel_fit_constrained.json + reports/kernel_fit_constrained.png; no model wiring,
no re-tune. Remaining caveats (documented, not fixed): still n_Ire not n_Eng; finite-island truncation.

Route-once-iterate-cheap: the routed times don't depend on f, so Phase A samples+routes+caches once
(resumable to data/_kernel_iter_cache_<purpose>.npz); Phase B only re-weights cached times per iter.

Run:  python3 analysis/iterate_kernel.py            (all six; needs OSRM up on :5000)
      python3 analysis/iterate_kernel.py --purpose res
"""
import argparse
import json
import os
import sys
import time
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# --- reuse the n(t) sampler geometry + OSRM plumbing ---
import build_n_of_t as B
from build_n_of_t import (load_area_masses, load_poi_layers, build_point_cache, _geo, _band_masses,
                          _dpt, osrm_table, _check_osrm, BANDS_KM, PURPOSES, SNAP_TOL_M)
# --- reuse the willingness divide + double-exp fit ---
from fit_kernel import (_tld_interp, equiv_miles_deriv, fit_double, fit_single, _wrms,
                        DRIVESHARE_MIN_FRAC, EFFN_MIN, TLD_FILE, NT_FILE)
from driveshare import driveshare, CURVES
from equiv_miles import equiv_miles

KERNEL_FIT = "analysis/kernel_fit.json"          # f^(0)
OUT_JSON = "analysis/kernel_fit_constrained.json"
PLOT_PATH = "reports/kernel_fit_constrained.png"
CACHE_TMPL = "data/_kernel_iter_cache_{p}.npz"    # gitignored, resumable

SEED = 20260703
M_ACC = 64          # dests per origin for D_i
B_ACC = 35          # origins per accessibility /table batch (B_ACC + M_ACC <= max-table-size)
NEAR_BUDGET = 150_000
FAR_BUDGET = 1_000_000
DESTS_PER_CALL = 49
B_FAR = 45
MAX_ITERS = 20
CONV_TOL = 0.01     # max relative change in (w, tau_s, tau_l) to stop
DAMP = 0.7          # relaxation on the param update (log-space for tau); 1.0 = pure Picard
T_FLOOR = 30.0      # clamp c before equiv_miles (avoid the log-quadratic blow-up at t→0)


# ── kernel evaluation ─────────────────────────────────────────────────────────
def kernel_f(comp, p, t):
    """Full kernel f(c) = driveshare(equiv_miles(c), comp) · W(c); W = w·e^{-c/τs}+(1-w)·e^{-c/τl}.
    Amplitude cancels in 1/D_i (global per purpose), so W is unit-amplitude."""
    tc = np.clip(t, T_FLOOR, None)
    W = p["w"] * np.exp(-tc / p["tau_s_s"]) + (1.0 - p["w"]) * np.exp(-tc / p["tau_l_s"])
    return driveshare(equiv_miles(tc), comp) * W


def _p0(comp):
    d = json.load(open(KERNEL_FIT))["components"][comp]["double"]
    return {"logA": float(np.log(max(d.get("A", 1.0), 1e-300))),
            "w": d["w"], "tau_s_s": d["tau_s_s"], "tau_l_s": d["tau_l_s"]}


def fit_double_warm(t, W, wts, x0):
    """Global double-exp fit via MULTI-START (keep the lowest-cost basin).  fit_kernel.fit_double's
    single fixed start lets the optimizer collapse to the w→1 degenerate basin (τl free, HIGHER
    residual) on the short-domain school fits, and flip basins across iterations. The good basin has
    the lower residual, so trying several starts (incl. the warm x0) and keeping the best finds it
    deterministically — no w cap needed, so the genuinely w≈1 non-school fits are untouched."""
    from scipy.optimize import least_squares
    logW = np.log(W); sw = np.sqrt(wts)
    A0 = float(np.average(logW, weights=wts))
    bounds = ([-np.inf, 0.0, 30.0, 100.0], [np.inf, 1.0, 1e5, 1e6])

    def resid(pp):
        logA, wgt, ts, tl = pp
        model = np.exp(logA) * (wgt * np.exp(-t / ts) + (1 - wgt) * np.exp(-t / tl))
        return (np.log(np.clip(model, 1e-300, None)) - logW) * sw

    starts = [x0, [A0, 0.6, 300.0, 4000.0], [A0, 0.9, 200.0, 1500.0], [A0, 0.97, 150.0, 800.0]]
    best, best_cost = None, np.inf
    for s in starts:
        s = [s[0], min(max(s[1], 1e-3), 1 - 1e-3), max(s[2], 30.0), max(s[3], 100.0)]
        try:
            r = least_squares(resid, s, bounds=bounds)
        except Exception:
            continue
        if r.cost < best_cost:
            best_cost, best = r.cost, r.x
    logA, wgt, ts, tl = best
    if ts > tl:
        ts, tl, wgt = tl, ts, 1.0 - wgt
    return float(logA), float(wgt), float(ts), float(tl)


# ── Phase A: sample + route + cache (once per purpose) ─────────────────────────
def build_accessibility(kind, ocent, prod, dloc, dw, dll, cache, rng):
    """Per-origin M dest-times (∝ A, global) → acc[N, M] (nan on route/snap fail). For D_i."""
    N = len(prod); dwn = dw / dw.sum()
    acc = np.full((N, M_ACC), np.nan, np.float32)
    t0 = time.time()
    for s in range(0, N, B_ACC):
        oidx = np.arange(s, min(s + B_ACC, N))
        di = rng.choice(len(dwn), size=M_ACC, p=dwn)          # M dests shared within this batch
        src = [cache[int(i)][rng.integers(len(cache[int(i)]))] for i in oidx]
        dst = [_dpt(kind, int(k), cache, dll, rng) for k in di]
        res = osrm_table(src, dst)
        if res is None:
            continue
        dur, ts, td = res
        good = (ts[:, None] < SNAP_TOL_M) & (td[None, :] < SNAP_TOL_M) & np.isfinite(dur)
        acc[oidx, :] = np.where(good, dur, np.nan).astype(np.float32)
        if (s // B_ACC) % 200 == 0:
            print(f"    acc {s:,}/{N:,} ({time.time()-t0:.0f}s)", flush=True)
    return acc


def build_density(kind, ocent, prod, dloc, dw, dll, cache, rng, edges_m):
    """Stratified pair sample, tagged by origin/band and cached (for Ñ reconstruction)."""
    tree, S, _, _ = _band_masses(ocent, prod, dloc, dw, edges_m)
    Sfar = np.maximum(0.0, dw.sum() - S.sum(axis=1))
    nb = len(edges_m) - 1
    n_oi, n_b, n_t = [], [], []
    for b in range(nb):
        ow = prod * S[:, b]
        if ow.sum() <= 0:
            continue
        ow = ow / ow.sum(); n = 0
        while n < NEAR_BUDGET:
            i = int(rng.choice(len(ow), p=ow))
            idx = np.asarray(tree.query_ball_point(ocent[i], edges_m[b + 1]))
            if idx.size:
                d = np.hypot(dloc[idx, 0] - ocent[i, 0], dloc[idx, 1] - ocent[i, 1])
                idx = idx[(d >= edges_m[b]) & (d < edges_m[b + 1])]
            if idx.size == 0:
                n += DESTS_PER_CALL; continue
            w = dw[idx]; w = w / w.sum()
            pick = idx[rng.choice(len(idx), size=DESTS_PER_CALL, p=w)]
            src = [cache[i][rng.integers(len(cache[i]))]]
            dst = [_dpt(kind, int(k), cache, dll, rng) for k in pick]
            n += DESTS_PER_CALL
            res = osrm_table(src, dst)
            if res is None:
                continue
            dur, ts, td = res
            ok = np.isfinite(dur[0]) & (td < SNAP_TOL_M) & (ts[0] < SNAP_TOL_M)
            tt = dur[0][ok]
            if tt.size:
                n_oi.append(np.full(tt.size, i, np.int32)); n_b.append(np.full(tt.size, b, np.int8))
                n_t.append(tt.astype(np.float32))
        print(f"    near band {BANDS_KM[b]:g}-{BANDS_KM[b+1]:g}km cached "
              f"{sum(x.size for x in n_t):,}", flush=True)
    # far tail (outer product, far-masked)
    f_oi, f_t = [], []
    owf = prod / prod.sum(); dwn = dw / dw.sum(); n = 0
    while n < FAR_BUDGET:
        oi = rng.choice(len(owf), size=B_FAR, p=owf); di = rng.choice(len(dwn), size=B_FAR, p=dwn)
        src = [cache[int(a)][rng.integers(len(cache[int(a)]))] for a in oi]
        dst = [_dpt(kind, int(k), cache, dll, rng) for k in di]
        n += B_FAR * B_FAR
        res = osrm_table(src, dst)
        if res is None:
            continue
        dur, ts, td = res
        dd = np.hypot(ocent[oi][:, None, 0] - dloc[di][None, :, 0],
                      ocent[oi][:, None, 1] - dloc[di][None, :, 1])
        ok = np.isfinite(dur) & (ts[:, None] < SNAP_TOL_M) & (td[None, :] < SNAP_TOL_M) \
            & (dd >= edges_m[-1])
        rows, _cols = np.where(ok)
        f_oi.append(oi[rows].astype(np.int32)); f_t.append(dur[ok].astype(np.float32))
    print(f"    far cached {sum(x.size for x in f_t):,}", flush=True)
    return (S.astype(np.float32), Sfar.astype(np.float32),
            np.concatenate(n_oi), np.concatenate(n_b), np.concatenate(n_t),
            np.concatenate(f_oi), np.concatenate(f_t))


def phase_a(name, df, pois, cache, rng):
    """Route + cache (resumable). Returns the geometry + cached samples for one purpose."""
    _, _, kind, ocent, prod, dloc, dw, dll = _geo(name, df, pois)
    cf = CACHE_TMPL.format(p=name)
    if os.path.exists(cf):
        z = np.load(cf)
        print(f"  [{name}] loaded cached sample from {cf}")
        return (kind, ocent, prod, dw, z["acc"], z["S"], z["Sfar"],
                z["n_oi"], z["n_b"], z["n_t"], z["f_oi"], z["f_t"])
    edges_m = np.array(BANDS_KM) * 1000.0
    print(f"  [{name}] accessibility pass (D_i) …", flush=True)
    acc = build_accessibility(kind, ocent, prod, dloc, dw, dll, cache, rng)
    print(f"  [{name}] stratified density pass (Ñ) …", flush=True)
    S, Sfar, n_oi, n_b, n_t, f_oi, f_t = build_density(kind, ocent, prod, dloc, dw, dll,
                                                       cache, rng, edges_m)
    np.savez(cf, acc=acc, S=S, Sfar=Sfar, n_oi=n_oi, n_b=n_b, n_t=n_t, f_oi=f_oi, f_t=f_t)
    print(f"  [{name}] cached → {cf}")
    return kind, ocent, prod, dw, acc, S, Sfar, n_oi, n_b, n_t, f_oi, f_t


# ── Phase B: fixed-point iteration (cheap) ────────────────────────────────────
def reconstruct_Ntilde(inv, S, Sfar, n_oi, n_b, n_t, f_oi, f_t, edges):
    """Ñ(t) per-bin mass = Σ_b M̃_b·ŝ_b, with origins reweighted by inv=P_i/D_i."""
    nb = S.shape[1]
    Mtil = (inv[:, None] * S).sum(axis=0)             # per near band
    Mtil_far = float((inv * Sfar).sum())
    Nt = np.zeros(len(edges) - 1)
    for b in range(nb):
        m = n_b == b
        if not m.any():
            continue
        h, _ = np.histogram(n_t[m], bins=edges, weights=inv[n_oi[m]])
        s = h.sum()
        if s > 0:
            Nt += Mtil[b] * h / s
    hf, _ = np.histogram(f_t, bins=edges, weights=inv[f_oi]); sf = hf.sum()
    if sf > 0:
        Nt += Mtil_far * hf / sf
    return Nt


def willingness_from_ndens(tld, comp, tc, n_dens_s):
    """W(t)=[TLD/Ñ]/driveshare on the seconds grid (fit_kernel.willingness with Ñ injected)."""
    dens_at, effn_at = _tld_interp(tld, comp)
    d = equiv_miles(tc); jac = equiv_miles_deriv(tc)
    tld_dens_s = dens_at(d) * jac
    ds = driveshare(d, comp); plateau = CURVES[comp][0]
    with np.errstate(divide="ignore", invalid="ignore"):
        W = (tld_dens_s / n_dens_s) / ds
    effn = effn_at(d)
    dom = (n_dens_s > 0) & (tld_dens_s > 0) & np.isfinite(W) & (W > 0) \
        & (ds >= DRIVESHARE_MIN_FRAC * plateau) & (effn >= EFFN_MIN)
    return W, effn, dom


def _rel_change(a, b):
    return max(abs(b[k] - a[k]) / max(abs(a[k]), 1e-9) for k in ("w", "tau_s_s", "tau_l_s"))


def _damp(old, new, alpha):
    """Relaxed update: w linear, tau geometric (log-space); logA carried from the latest fit."""
    return {"logA": new["logA"],
            "w": (1 - alpha) * old["w"] + alpha * new["w"],
            "tau_s_s": float(np.exp((1 - alpha) * np.log(old["tau_s_s"]) + alpha * np.log(new["tau_s_s"]))),
            "tau_l_s": float(np.exp((1 - alpha) * np.log(old["tau_l_s"]) + alpha * np.log(new["tau_l_s"])))}


def iterate(name, tld, edges, tc, width, geom, single=False):
    kind, ocent, prod, dw, acc, S, Sfar, n_oi, n_b, n_t, f_oi, f_t = geom
    if single:                                                        # single-exp: W=exp(-t/τ), w≡1
        p = {"logA": 0.0, "w": 1.0,
             "tau_s_s": json.load(open(KERNEL_FIT))["components"][name]["tau_single_s"],
             "tau_l_s": 1e9}
    else:
        p = _p0(name)
    p0 = {k: p[k] for k in ("w", "tau_s_s", "tau_l_s")}
    trace = []
    for k in range(MAX_ITERS):
        with warnings.catch_warnings():                                # some origins are unroutable
            warnings.simplefilter("ignore", RuntimeWarning)            # (all-nan rows → median below)
            D = np.nanmean(kernel_f(name, p, acc), axis=1)             # per-origin accessibility
        D = np.where(np.isfinite(D) & (D > 0), D, np.nanmedian(D[np.isfinite(D)]))
        inv = prod / D                                                 # P_i / D_i
        Nt = reconstruct_Ntilde(inv, S, Sfar, n_oi, n_b, n_t, f_oi, f_t, edges)
        n_dens_s = Nt / width
        W, effn, dom = willingness_from_ndens(tld, name, tc, n_dens_s)
        if single:
            tau, logA = fit_single(tc[dom], W[dom], effn[dom])        # identifiable (1 param)
            new = {"logA": logA, "w": 1.0, "tau_s_s": tau, "tau_l_s": 1e9}
            model = lambda tt, _a=logA, _t=tau: np.exp(_a) * np.exp(-tt / _t)     # noqa: E731
        else:
            x0 = [p["logA"], p["w"], p["tau_s_s"], p["tau_l_s"]]       # warm-start (stay in-basin)
            logA, w, ts, tl = fit_double_warm(tc[dom], W[dom], effn[dom], x0)
            new = {"logA": logA, "w": w, "tau_s_s": ts, "tau_l_s": tl}
            model = lambda tt, _a=logA, _w=w, _s=ts, _l=tl: \
                np.exp(_a) * (_w * np.exp(-tt / _s) + (1 - _w) * np.exp(-tt / _l))   # noqa: E731
        wr = _wrms(tc[dom], W[dom], effn[dom], model)
        rc = _rel_change(p, new)
        trace.append({"iter": k, "w": new["w"], "tau_s_s": new["tau_s_s"], "tau_l_s": new["tau_l_s"],
                      "wrms": wr, "rel_change": rc, "n_points": int(dom.sum())})
        print(f"    it{k}: w={new['w']:.3f} tau_s={new['tau_s_s']:6.0f}s tau_l={new['tau_l_s']:7.0f}s  "
              f"wrms={wr:.3f}  Δ={rc:.3f}", flush=True)
        p = _damp(p, new, DAMP)
        if rc < CONV_TOL:
            break
    converged = trace[-1]["rel_change"] < CONV_TOL
    if converged:
        di = {"w": p["w"], "tau_s_s": p["tau_s_s"], "tau_l_s": p["tau_l_s"]}
    else:
        # weakly-identified tail (school postprimary/tertiary): the fit alternates between an
        # identifiable basin and a τl→∞ runaway. The runaways are a minority, so the per-iteration
        # MEDIAN is a robust identifiable-basin estimate (τs stays stable regardless).
        med = lambda k: float(np.median([t[k] for t in trace]))            # noqa: E731
        di = {"w": med("w"), "tau_s_s": med("tau_s_s"), "tau_l_s": med("tau_l_s")}
    return {"p0": p0, "converged": converged, "tail_weakly_identified": not converged,
            "n_iter": len(trace), "double_iterated": di, "trace": trace}


# ── driver ────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--purpose", default=None, choices=list(PURPOSES))
    ap.add_argument("--single", action="store_true",
                    help="prototype: single-exp willingness (default school levels); print-only, no save")
    args = ap.parse_args()

    _check_osrm()
    tld = json.load(open(TLD_FILE))
    nt = json.load(open(NT_FILE))
    edges = np.array(nt["bin_edges_s"], float)
    fin = np.isfinite(edges[1:])
    tc = 0.5 * (edges[:-1][fin] + edges[1:][fin])
    width = np.diff(edges)[fin]
    edges_fin = np.append(edges[:-1][fin], edges[1:][fin][-1])         # finite bin edges for histograms

    df, geoms = load_area_masses()
    pois = load_poi_layers()
    rng = np.random.default_rng(SEED)
    cache = build_point_cache(df, geoms, rng)

    if args.single:                                                   # prototype: single-exp, print-only
        names = [args.purpose] if args.purpose else \
            ["school_primary", "school_postprimary", "school_tertiary"]
        kf = json.load(open(KERNEL_FIT))["components"]
        print("\nSINGLE-EXP prototype (print-only; compare fit vs double-exp):")
        print(f"  {'purpose':20s}{'conv':>6}{'it':>4}  {'tau uncon->iter':>18}  "
              f"{'wrms single':>12}  {'wrms double(ref)':>16}")
        for name in names:
            geom = phase_a(name, df, pois, cache, rng)
            res = iterate(name, tld, edges_fin, tc, width, geom, single=True)
            it = res["double_iterated"]; last = res["trace"][-1]
            u_single = kf[name]["tau_single_s"]; u_dbl_wrms = kf[name]["wrms_double"]
            print(f"  {name:20s}{str(res['converged']):>6}{res['n_iter']:>4}  "
                  f"{u_single:7.0f}->{it['tau_s_s']:<7.0f}s  {last['wrms']:>12.3f}  {u_dbl_wrms:>16.3f}")
        return

    names = [args.purpose] if args.purpose else list(PURPOSES)
    results = {}
    for name in names:
        print(f"\n=== {name} ===")
        geom = phase_a(name, df, pois, cache, rng)
        # reconstruct on finite bins only (drop the [last,inf) overflow, matching tc/width)
        res = iterate(name, tld, edges_fin, tc, width, geom)
        u = json.load(open(KERNEL_FIT))["components"][name]["double"]
        it = res["double_iterated"]
        print(f"  {name}: tau_l {u['tau_l_s']:.0f}s (uncon) -> {it['tau_s_s']:.0f}/{it['tau_l_s']:.0f}s "
              f"(iterated tau_s/tau_l); converged={res['converged']}")
        results[name] = res
    _write(results, tld)
    _plot(results, tld, tc, width, edges_fin)


def _write(results, tld):
    out = {"_meta": {
        "purpose": "constrained (1/D_i) iterated double-exp willingness kernels",
        "method": "fixed-point f -> D_i -> Ñ=Σ(P_i/D_i)A_j δ -> fit_double([TLD/Ñ]/driveshare)",
        "input_kernel": KERNEL_FIT, "tld_file": TLD_FILE, "n_of_t_file": NT_FILE,
        "schools": "per-level (one kernel each, NOT shared)",
        "damp": DAMP, "conv_tol": CONV_TOL, "seed": SEED,
        "removes": "the 1/D_i production-constraint mis-attribution that inflated tau_l in the "
                   "unconstrained fit",
        "still_caveated": "n_Ire not n_Eng (source-region geometry ratio); finite-island truncation; "
                          "n(t) sampling from the same v1 machinery.",
        "converged_note": "res/commute/retail/school_primary converge cleanly; tau_l shortens (the "
                          "1/D_i tail de-inflation).",
        "weak_tail_note": "school_postprimary + school_tertiary have tail_weakly_identified=true: "
                          "short fit domains (school trips are short-range, thin long-tail TLD) make "
                          "tau_l bimodal (identifiable basin vs a τl→∞ runaway), so it does not "
                          "converge. double_iterated.tau_l is the robust per-iteration MEDIAN "
                          "(identifiable basin; tau_s is stable). Single-exp was worse — it "
                          "destabilises the iteration (τ collapses, D_i underflows).",
    }, "components": {}}
    for name, r in results.items():
        out["components"][name] = {k: r[k] for k in
                                   ("p0", "double_iterated", "converged", "tail_weakly_identified",
                                    "n_iter", "trace")}
    json.dump(out, open(OUT_JSON, "w"), indent=2)
    print(f"\nSaved -> {OUT_JSON}")


def _plot(results, tld, tc, width, edges):
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    except Exception as e:
        print(f"(plot skipped: {e})"); return
    comps = list(results)
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    for ax, name in zip(axes.flat, comps):
        r = results[name]
        for lbl, key, style in (("uncon (f0)", "p0", "--"), ("iterated", "double_iterated", "-")):
            p = r[key]
            W = p["w"] * np.exp(-tc / p["tau_s_s"]) + (1 - p["w"]) * np.exp(-tc / p["tau_l_s"])
            W = W / W[np.searchsorted(tc, 600)]                       # norm @10 min
            ax.plot(tc / 60, W, style, lw=1.8,
                    label=f"{lbl}: τs={p['tau_s_s']:.0f} τl={p['tau_l_s']:.0f}")
        ax.set_yscale("log"); ax.set_title(name, fontsize=10)
        ax.set_xlabel("min"); ax.set_ylabel("W (norm@10min)"); ax.legend(fontsize=7.5)
        ax.grid(alpha=0.3); ax.set_xlim(0, 120)
    fig.suptitle("Constrained (1/D_i) iterated willingness vs unconstrained")
    fig.tight_layout(); os.makedirs(os.path.dirname(PLOT_PATH), exist_ok=True)
    fig.savefig(PLOT_PATH, dpi=115); print(f"saved {PLOT_PATH}")


if __name__ == "__main__":
    main()
