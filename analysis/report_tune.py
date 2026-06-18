"""
Generate a structured report from a tuning history entry.

Usage:
  python3 analysis/report_tune.py             # last history entry
  python3 analysis/report_tune.py <run_id>    # specific entry by ID prefix

Outputs (in reports/):
  tune_report_{id}.txt   — text tables: summary, chi² by measurement, chi² by link,
                            gravity params, external zones, slot fractions
  slot_pulls_{id}.png    — horizontal pull plot for slot fraction constraints
"""

import json, math, os, sys

HISTORY_FILE  = "simulation/tuning_history.jsonl"
TUNER_CONFIG  = "simulation/tuner_config.json"
REPORTS_DIR   = "reports"

_DT_NAMES = {0: "Wkday", 1: "Sat", 2: "Sun"}


# ── Load history entry ────────────────────────────────────────────────────────

def _load_entry(run_id_prefix=None):
    with open(HISTORY_FILE) as f:
        lines = [l.strip() for l in f if l.strip()]
    entries = [json.loads(l) for l in lines]
    if not entries:
        raise SystemExit("tuning_history.jsonl is empty")
    if run_id_prefix is None:
        return entries[-1], entries[-2] if len(entries) >= 2 else None
    matches = [e for e in entries if e["id"].startswith(run_id_prefix)]
    if not matches:
        raise SystemExit(f"No history entry matching prefix '{run_id_prefix}'")
    idx = entries.index(matches[-1])
    return entries[idx], entries[idx - 1] if idx > 0 else None


# ── Helpers ───────────────────────────────────────────────────────────────────

def _rule(title, width=72):
    return f"── {title} {'─' * max(0, width - len(title) - 4)}"


def _fmt_pct(val):
    return f"{val:+.1f}%"


# ── Section: SUMMARY ──────────────────────────────────────────────────────────

def _section_summary(e):
    lines = [_rule("SUMMARY"), ""]
    ts = e.get("timestamp", "")[:16].replace("T", "  ")
    lines.append(
        f"  Run {e['id']}  {ts}  stage={e['stage']}  "
        f"N={e['n_obs']}  N_eff={e['n_eff']}  N_slots={e['n_slots']}"
    )
    chi2_n_eff = e["chi2"] / e["n_eff"] if e["n_eff"] else float("nan")
    lines.append(
        f"  χ²/N={e['chi2_per_n']:.4f}  χ²/N_eff={chi2_n_eff:.4f}  "
        f"N_params={e['n_params']}  n_evals={e['n_evals']}"
    )
    note = e.get("note")
    if note:
        lines.append(f"  note: {note}")
    lines.append("")
    return lines


# ── Section: CHI² BY MEASUREMENT ─────────────────────────────────────────────

def _section_by_measurement(e):
    obs = e["observations"]
    lines = [_rule("CHI² TABLE: BY MEASUREMENT"), ""]
    LW = 54
    header = (
        f"  {'':1s}  {'Kind':<8}  {'Label':<{LW}}  "
        f"{'Obs':>9}  {'σ':>8}  {'Model':>9}  {'z':>7}"
    )
    lines.append(header)
    lines.append("  " + "─" * (len(header) - 2))

    rows = sorted(obs, key=lambda o: abs(o["z"]), reverse=True)
    for o in rows:
        marker = "*" if abs(o["z"]) > 2 else " "
        lbl = o.get("label", "?")
        if len(lbl) > LW:
            lbl = lbl[:LW - 1] + "…"
        lines.append(
            f"  {marker}  {o['kind']:<8}  {lbl:<{LW}}  "
            f"{o['observed']:>9,.0f}  {o['sigma']:>8,.0f}  "
            f"{o['model']:>9,.0f}  {o['z']:>+7.2f}"
        )

    chi2 = e["chi2"]
    n = e["n_obs"]
    n_eff = e["n_eff"]
    zs = [abs(o["z"]) for o in obs]
    n_out2 = sum(1 for z in zs if z > 2)
    n_out3 = sum(1 for z in zs if z > 3)
    lines.append("")
    lines.append(
        f"  n={n}  χ²={chi2:.2f}  χ²/N={chi2/n:.4f}  "
        f"mean|z|={sum(zs)/len(zs):.2f}  |z|>2: {n_out2}  |z|>3: {n_out3}"
    )
    lines.append(f"  N_eff={n_eff}  χ²/N_eff={chi2/n_eff:.4f}")
    lines.append("")
    return lines


# ── Section: CHI² BY LINK ─────────────────────────────────────────────────────

def _section_by_link(e):
    obs = e["observations"]
    lines = [_rule("CHI² TABLE: BY LINK"), ""]
    LW = 54

    # Group observations by site/link.
    # Walking obs: group by label (unique per directed link, same across sessions).
    # Official hourly obs: group by site (strip the " hXX" hour suffix from label).
    link_groups = {}
    for o in obs:
        kind = o["kind"]
        lbl  = o.get("label", "?")
        if kind in ("official", "official_hourly"):
            # "731 h08" → "731";  "site 507, A21 ... h08" → "site 507, A21 ..."
            site_key = lbl.rsplit(" h", 1)[0] if " h" in lbl else lbl
            key = ("official_hourly", site_key)
        else:
            key = ("walking", lbl)
        link_groups.setdefault(key, []).append(o)

    # Build rows: (Σz², max|z|, z_min, z_max, N_sess, lbl, kind)
    rows = []
    for (kind, grp_key), group in link_groups.items():
        zvals  = [o["z"] for o in group]
        sum_z2 = sum(z * z for z in zvals)
        z_min  = min(zvals)
        z_max  = max(zvals)
        max_az = max(abs(z) for z in zvals)
        lbl    = group[0].get("label", str(grp_key))
        rows.append((sum_z2, max_az, z_min, z_max, len(group), lbl, kind))

    rows.sort(key=lambda r: r[0], reverse=True)

    header = (
        f"  {'':1s}  {'Kind':<8}  {'Label':<{LW}}  "
        f"{'N':>3}  {'Σz²':>7}  {'z_min':>7}  {'z_max':>7}  {'|z|_max':>7}"
    )
    lines.append(header)
    lines.append("  " + "─" * (len(header) - 2))

    for (sum_z2, max_az, z_min, z_max, n_sess, lbl, kind) in rows:
        marker = "*" if max_az > 2 else " "
        if len(lbl) > LW:
            lbl = lbl[:LW - 1] + "…"
        lines.append(
            f"  {marker}  {kind:<8}  {lbl:<{LW}}  "
            f"{n_sess:>3}  {sum_z2:>7.2f}  {z_min:>+7.2f}  {z_max:>+7.2f}  {max_az:>7.2f}"
        )

    lines.append("")
    return lines


# ── Section: GRAVITY PARAMETERS ───────────────────────────────────────────────

def _section_gravity(e, prev_entry, config):
    params      = e["params"]
    grav_ref    = config.get("gravity_ref", {})
    grav_lam    = config.get("gravity_lambda", 0.0)

    # σ in log-space = 1/sqrt(grav_lam); pull = log(final/ref) * sqrt(grav_lam)
    def _log_pull(val, ref):
        if grav_lam <= 0 or ref <= 0 or val <= 0:
            return float("nan")
        return math.log(val / ref) * math.sqrt(grav_lam)

    # Determine initial gravity: prefer entry's own field, else previous entry, else grav_ref
    init_src = "run start"
    if "initial_gravity" in e:
        init_g = e["initial_gravity"]
    elif prev_entry is not None:
        init_g = {k: prev_entry["params"].get(k) for k in ("W_BIZ", "P", "ALPHA", "BETA", "THETA")}
        init_g = {k: v for k, v in init_g.items() if v is not None}
        init_src = f"prev run ({prev_entry['id'][:8]})"
    else:
        init_g = {k: grav_ref.get(k) for k in ("W_BIZ", "P", "ALPHA", "BETA")}
        init_g = {k: v for k, v in init_g.items() if v is not None}
        init_src = "gravity_ref (no prior run)"

    lines = [_rule("GRAVITY PARAMETERS"), ""]
    lines.append(f"  Initial from: {init_src}")
    lines.append("")

    grav_keys = [k for k in ("W_BIZ", "P", "ALPHA", "BETA", "THETA") if k in params]
    header = (
        f"  {'Param':<8}  {'Initial':>12}  {'Final':>12}  "
        f"{'Δ%':>8}  {'ref':>12}  {'Δ%_ref':>8}  {'pull_ref':>9}"
    )
    lines.append(header)
    lines.append("  " + "─" * (len(header) - 2))

    for k in grav_keys:
        final = params[k]
        init  = init_g.get(k)
        ref   = grav_ref.get(k)
        d_pct  = 100.0 * (final - init) / init if init else float("nan")
        dr_pct = 100.0 * (final - ref)  / ref  if ref  else float("nan")
        pull   = _log_pull(final, ref)
        init_s = f"{init:>12.4f}" if init is not None else f"{'—':>12}"
        ref_s  = f"{ref:>12.4f}"  if ref  is not None else f"{'—':>12}"
        d_s    = _fmt_pct(d_pct)  if not math.isnan(d_pct)  else "    —"
        dr_s   = _fmt_pct(dr_pct) if not math.isnan(dr_pct) else "    —"
        pl_s   = f"{pull:>+9.3f}" if not math.isnan(pull) else f"{'—':>9}"
        lines.append(
            f"  {k:<8}  {init_s}  {final:>12.4f}  "
            f"{d_s:>8}  {ref_s}  {dr_s:>8}  {pl_s}"
        )

    # K components (analytical, not tuned)
    K     = params.get("K")
    K_res = params.get("K_res")
    K_biz = params.get("K_biz")
    if K is not None:
        lines.append(
            f"  {'K':<8}  {'(analytical)':>12}  {K:>12.4e}  "
            f"{'':>8}  {'—':>12}  {'':>8}  {'':>9}"
        )
    if K_res is not None and K_biz is not None and K:
        phi = K_biz / K
        lines.append(f"  {'K_res':<8}  {'':>12}  {K_res:>12.4e}  "
                     f"{'':>8}  {'':>12}  {'':>8}  {'':>9}  ({100*(1-phi):.1f}% of K)")
        lines.append(f"  {'K_biz':<8}  {'':>12}  {K_biz:>12.4e}  "
                     f"{'':>8}  {'':>12}  {'':>8}  {'':>9}  ({100*phi:.1f}% of K, phi)")
        phi_p = e.get("tuner_hyperparams", {}).get("phi_prior", 0.35)
        lines.append(f"  {'phi':<8}  {'':>12}  {phi:>12.4f}  "
                     f"{'':>8}  {phi_p:>12.4f}  {'':>8}  {'':>9}  (phi_prior)")

    lines.append("")
    lines.append(
        f"  gravity_ref: " +
        "  ".join(f"{k}={v}" for k, v in grav_ref.items())
    )
    if grav_lam > 0:
        sigma_log = 1.0 / math.sqrt(grav_lam)
        lines.append(
            f"  gravity_lambda={grav_lam}  "
            f"→ L2 σ_log={sigma_log:.3f} (pull = log(final/ref)/σ_log)"
        )
    lines.append("")
    return lines


# ── Section: EXTERNAL ZONES ───────────────────────────────────────────────────

def _section_external(e, config):
    if e["stage"] != "full":
        return [_rule("EXTERNAL ZONES"), "", "  (gravity-only run — no external zone params)", ""]

    params    = e["params"]
    city_pop  = params.get("external_city_pop", {})
    city_wp   = params.get("external_city_wp", {})
    dampings  = params.get("external_dampings", {})
    cities    = config.get("cities", {})

    lines = [_rule("EXTERNAL ZONES"), ""]

    def _flag(dpct):
        a = abs(dpct)
        if a > 100: return "***"
        if a > 50:  return " **"
        if a > 25:  return "  *"
        return "   "

    header = (
        f"  {'City':<12}  {'ref_pop':>10}  {'tuned_pop':>10}  {'Δpop%':>8}  {'':3}  "
        f"{'ref_wp':>9}  {'tuned_wp':>9}  {'Δwp%':>8}  {'':3}"
    )
    lines.append(header)
    lines.append("  " + "─" * (len(header) - 2))

    for city_name, city_cfg in cities.items():
        rp = city_cfg["ref_pop"]
        rw = city_cfg["ref_wp"]
        tp = city_pop.get(city_name, rp)
        tw = city_wp.get(city_name, rw)
        dp = 100.0 * (tp - rp) / rp if rp else float("nan")
        dw = 100.0 * (tw - rw) / rw if rw else float("nan")
        lines.append(
            f"  {city_name:<12}  {rp:>10,.0f}  {tp:>10,.0f}  {_fmt_pct(dp):>8}  {_flag(dp)}  "
            f"{rw:>9,.0f}  {tw:>9,.0f}  {_fmt_pct(dw):>8}  {_flag(dw)}"
        )

    if dampings:
        lines.append("")
        lines.append(f"  Tunable dampings (ref → tuned):")
        for node_str, tuned_d in dampings.items():
            # Find city and ref damping
            city_label = "?"
            ref_d = None
            for cn, cc in cities.items():
                d_map = cc.get("dampings", {})
                if node_str in d_map:
                    city_label = cn
                    ref_d = d_map[node_str]
                    break
            if ref_d is not None:
                dd = 100.0 * (tuned_d - ref_d) / ref_d
                lines.append(
                    f"    node {node_str:<6}  ({city_label})  "
                    f"ref={ref_d:.4f}  →  tuned={tuned_d:.4f}  {_fmt_pct(dd)}"
                )

    lines.append("")
    lines.append("  *** = |Δ|>100%,  ** = |Δ|>50%,  * = |Δ|>25%")
    lines.append("")
    return lines


# ── Section: SLOT FRACTIONS ───────────────────────────────────────────────────

def _section_slots(e):
    params         = e["params"]
    slot_fracs_res = params.get("slot_fracs_res", {})
    slot_fracs_biz = params.get("slot_fracs_biz", {})
    slot_prior     = e.get("slot_prior", {})

    if not slot_fracs_res and not slot_fracs_biz:
        return [_rule("SLOT FRACTIONS"), "", "  (no slot fraction data)", ""]

    lines = [_rule("SLOT FRACTIONS"), ""]
    # Columns: type, hour, prior_agg, f_res, pull_res, f_biz, pull_biz, coupling
    header = (
        f"  {'Type':<5}  {'Hr':>2}  {'Prior_agg':>9}  "
        f"{'f_res':>8}  {'Δ/σ_res':>7}  "
        f"{'f_biz':>8}  {'Δ/σ_biz':>7}  "
        f"{'f_r+f_b':>8}  {'2·prior':>8}"
    )
    lines.append(header)
    lines.append("  " + "─" * (len(header) - 2))

    all_keys = sorted(set(slot_fracs_res) | set(slot_fracs_biz),
                      key=lambda s: (int(s.split(",")[0]), int(s.split(",")[1])))

    for sk_str in all_keys:
        dt, h = int(sk_str.split(",")[0]), int(sk_str.split(",")[1])
        f_r = slot_fracs_res.get(sk_str, float("nan"))
        f_b = slot_fracs_biz.get(sk_str, float("nan"))

        if sk_str in slot_prior:
            mfa, std_f, mfr, mfb = slot_prior[sk_str]
        else:
            mfa = std_f = mfr = mfb = float("nan")

        def _pull(f, prior):
            if std_f > 0 and not any(math.isnan(x) for x in (f, prior, std_f)):
                return (f - prior) / std_f
            return float("nan")

        pr = _pull(f_r, mfr)
        pb = _pull(f_b, mfb)
        f_sum  = f_r + f_b if not (math.isnan(f_r) or math.isnan(f_b)) else float("nan")
        f_agg2 = mfa if not math.isnan(mfa) else float("nan")

        def _fs(v): return f"{v:>8.5f}" if not math.isnan(v) else f"{'—':>8}"
        def _ps(v): return f"{v:>+7.2f}" if not math.isnan(v) else f"{'—':>7}"

        lines.append(
            f"  {_DT_NAMES[dt]:<5}  {h:>2d}  {mfa:>9.5f}  "
            f"{_fs(f_r)}  {_ps(pr)}  "
            f"{_fs(f_b)}  {_ps(pb)}  "
            f"{_fs(f_sum)}  {_fs(f_agg2)}"
        )

    lines.append("")
    return lines


# ── Pull plot ─────────────────────────────────────────────────────────────────

def _make_pull_plot(e, out_path):
    """
    Two side-by-side heatmaps: residential pull and business pull.
    Each heatmap is 24 rows (hour 0–23) × 3 columns (Wkday / Sat / Sun).
    Cell colour = (f_inferred − f_prior_component) / σ_agg  (diverging: blue < 0 < red).
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("matplotlib not available — skipping pull plot")
        return

    params         = e["params"]
    slot_fracs_res = params.get("slot_fracs_res", {})
    slot_fracs_biz = params.get("slot_fracs_biz", {})
    slot_prior     = e.get("slot_prior", {})

    if (not slot_fracs_res and not slot_fracs_biz) or not slot_prior:
        print("  No slot prior data — skipping pull plot")
        return

    # Build pull matrices: shape (24 hours, 3 day-types)
    # Day-type columns: 0=Wkday, 1=Sat, 2=Sun
    pull_res = np.full((24, 3), np.nan)
    pull_biz = np.full((24, 3), np.nan)

    for sk_str, prior_vals in slot_prior.items():
        dt, h = int(sk_str.split(",")[0]), int(sk_str.split(",")[1])
        if len(prior_vals) < 4:
            continue
        mfa, std_f, mfr, mfb = prior_vals
        if std_f <= 0 or math.isnan(std_f):
            continue
        f_r = slot_fracs_res.get(sk_str)
        f_b = slot_fracs_biz.get(sk_str)
        if f_r is not None:
            pull_res[h, dt] = (f_r - mfr) / std_f
        if f_b is not None:
            pull_biz[h, dt] = (f_b - mfb) / std_f

    vmax = max(
        np.nanmax(np.abs(pull_res)) if not np.all(np.isnan(pull_res)) else 1.0,
        np.nanmax(np.abs(pull_biz)) if not np.all(np.isnan(pull_biz)) else 1.0,
        1.0,
    )
    vmax = math.ceil(vmax * 10) / 10   # round up to 1 d.p.

    col_labels = ["Wkday", "Sat", "Sun"]
    hour_labels = [f"{h:02d}h" for h in range(24)]

    fig, axes = plt.subplots(1, 2, figsize=(10, 9), sharey=True)
    fig.subplots_adjust(wspace=0.08)

    for ax, pull_mat, comp_label, prior_key in [
        (axes[0], pull_res, "Residential (f_res)", "mean_f_res"),
        (axes[1], pull_biz, "Business (f_biz)",    "mean_f_biz"),
    ]:
        im = ax.imshow(pull_mat, aspect="auto", cmap="RdBu_r",
                       vmin=-vmax, vmax=vmax, origin="upper")
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(col_labels, fontsize=9)
        ax.set_title(f"{comp_label}\npull = (f − prior) / σ_agg", fontsize=9)
        # Annotate cells with pull value
        for h in range(24):
            for dt in range(3):
                v = pull_mat[h, dt]
                if not math.isnan(v):
                    ax.text(dt, h, f"{v:+.1f}", ha="center", va="center",
                            fontsize=6.5, color="white" if abs(v) > vmax * 0.6 else "black")

    axes[0].set_yticks(range(24))
    axes[0].set_yticklabels(hour_labels, fontsize=8)
    axes[0].set_ylabel("Hour of day")

    cbar = fig.colorbar(im, ax=axes, orientation="vertical", fraction=0.025, pad=0.02)
    cbar.set_label("Pull (σ)", fontsize=9)

    fig.suptitle(
        f"Temporal profile pulls — run {e['id'][:8]}  "
        f"χ²/N={e['chi2_per_n']:.4f}  stage={e['stage']}",
        fontsize=10, y=0.995
    )
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {out_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    run_prefix = sys.argv[1] if len(sys.argv) > 1 else None
    entry, prev_entry = _load_entry(run_prefix)

    with open(TUNER_CONFIG) as f:
        config = json.load(f)

    run_id = entry["id"]
    os.makedirs(REPORTS_DIR, exist_ok=True)

    # ── Text report ───────────────────────────────────────────────────────────
    report_path = os.path.join(REPORTS_DIR, f"tune_report_{run_id[:8]}.txt")
    sections = (
        _section_summary(entry)
        + _section_by_measurement(entry)
        + _section_by_link(entry)
        + _section_gravity(entry, prev_entry, config)
        + _section_external(entry, config)
        + _section_slots(entry)
    )
    text = "\n".join(sections)
    with open(report_path, "w") as f:
        f.write(text)
    print(f"  Saved → {report_path}")

    # ── Pull plot ─────────────────────────────────────────────────────────────
    plot_path = os.path.join(REPORTS_DIR, f"slot_pulls_{run_id[:8]}.png")
    _make_pull_plot(entry, plot_path)

    print(f"\nDone.  Reports in {REPORTS_DIR}/")


if __name__ == "__main__":
    main()
