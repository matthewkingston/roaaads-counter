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
        lbl = o.get("label", str(o["target"]))
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

    # Group walking obs by canonical target key; official sites kept as-is
    link_groups = {}  # key → list of obs
    for o in obs:
        if o["kind"] in ("official", "official_hourly"):
            key = ("official_hourly", o["label"])
        else:
            key = ("walking", tuple(o["target"]))
        link_groups.setdefault(key, []).append(o)

    # Build rows: (Σz², max|z|, z_min, z_max, N_sess, lbl, kind)
    rows = []
    for (kind, target), group in link_groups.items():
        zvals  = [o["z"] for o in group]
        sum_z2 = sum(z * z for z in zvals)
        z_min  = min(zvals)
        z_max  = max(zvals)
        max_az = max(abs(z) for z in zvals)
        lbl    = group[0].get("label", str(target))
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
        init_g = {k: prev_entry["params"].get(k) for k in ("W_BIZ", "P", "ALPHA", "THETA")}
        init_g = {k: v for k, v in init_g.items() if v is not None}
        init_src = f"prev run ({prev_entry['id'][:8]})"
    else:
        init_g = {k: grav_ref.get(k) for k in ("W_BIZ", "P", "ALPHA")}
        init_g = {k: v for k, v in init_g.items() if v is not None}
        init_src = "gravity_ref (no prior run)"

    lines = [_rule("GRAVITY PARAMETERS"), ""]
    lines.append(f"  Initial from: {init_src}")
    lines.append("")

    grav_keys = [k for k in ("W_BIZ", "P", "ALPHA", "THETA") if k in params]
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

    # K (analytical, not tuned)
    K = params.get("K")
    if K is not None:
        lines.append(
            f"  {'K':<8}  {'(analytical)':>12}  {K:>12.4e}  "
            f"{'':>8}  {'—':>12}  {'':>8}  {'':>9}"
        )

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
    params     = e["params"]
    slot_fracs = params.get("slot_fracs", {})
    slot_prior = e.get("slot_prior", {})

    if not slot_fracs:
        return [_rule("SLOT FRACTIONS"), "", "  (no slot fraction data)", ""]

    lines = [_rule("SLOT FRACTIONS"), ""]
    header = (
        f"  {'Type':<5}  {'Hr':>2}  {'Prior':>10}  {'Inferred':>10}  "
        f"{'Δ':>10}  {'Δ%':>7}  {'|Δ|/σ':>7}  {'N':>3}"
    )
    lines.append(header)
    lines.append("  " + "─" * (len(header) - 2))

    # Build rows sorted by |pull| descending
    rows = []
    for sk_str, f_s in slot_fracs.items():
        dt_str, h_str = sk_str.split(",")
        dt, h = int(dt_str), int(h_str)
        if sk_str in slot_prior:
            mean_f, std_f = slot_prior[sk_str][:2]
        else:
            mean_f, std_f = f_s, float("nan")
        pull    = abs(f_s - mean_f) / std_f if std_f > 0 and not math.isnan(std_f) else float("nan")
        delta   = f_s - mean_f
        dpct    = 100.0 * delta / mean_f if mean_f != 0 else float("nan")
        rows.append((pull, dt, h, mean_f, std_f, f_s, delta, dpct))

    rows.sort(key=lambda r: r[0] if not math.isnan(r[0]) else -1, reverse=True)

    for (pull, dt, h, mean_f, std_f, f_s, delta, dpct) in rows:
        pull_s = f"{pull:>7.2f}σ" if not math.isnan(pull) else f"{'—':>7s} "
        dpct_s = _fmt_pct(dpct) if not math.isnan(dpct) else "    —"
        # N: not stored in history; show "?" to indicate unknown count
        lines.append(
            f"  {_DT_NAMES[dt]:<5}  {h:>2d}  {mean_f:>10.6f}  {f_s:>10.6f}  "
            f"{delta:>+10.6f}  {dpct_s:>7}  {pull_s}  {'?':>3}"
        )

    lines.append("")
    return lines


# ── Pull plot ─────────────────────────────────────────────────────────────────

def _make_pull_plot(e, out_path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("matplotlib not available — skipping pull plot")
        return

    params     = e["params"]
    slot_fracs = params.get("slot_fracs", {})
    slot_prior = e.get("slot_prior", {})

    if not slot_fracs or not slot_prior:
        print("  No slot prior data — skipping pull plot")
        return

    # Build (pull, label) pairs sorted by |pull| descending
    rows = []
    for sk_str, f_s in slot_fracs.items():
        if sk_str not in slot_prior:
            continue
        mean_f, std_f = slot_prior[sk_str]
        if std_f <= 0 or math.isnan(std_f):
            continue
        pull  = (f_s - mean_f) / std_f
        dt, h = int(sk_str.split(",")[0]), int(sk_str.split(",")[1])
        label = f"{_DT_NAMES[dt]} {h:02d}h"
        rows.append((abs(pull), pull, label))

    rows.sort(key=lambda r: r[0], reverse=True)
    if not rows:
        print("  No usable slot fraction data — skipping pull plot")
        return

    _, pulls, labels = zip(*rows)
    pulls  = np.array(pulls)
    y_pos  = np.arange(len(pulls))

    fig_h  = max(3.5, len(pulls) * 0.55 + 1.0)
    fig, ax = plt.subplots(figsize=(8, fig_h))

    # Shaded ±1 band
    ax.axvspan(-1, 1, color="lightgray", alpha=0.35, zorder=0)

    # Reference lines
    for xv, ls in [(-2, ":"), (-1, "--"), (0, "-"), (1, "--"), (2, ":")]:
        ax.axvline(xv, color="gray", linewidth=0.8, linestyle=ls, zorder=1)

    # Error bars: dot at pull, bar from (pull−1) to (pull+1)
    xerr_lo = np.ones_like(pulls)
    xerr_hi = np.ones_like(pulls)
    ax.errorbar(
        pulls, y_pos,
        xerr=[xerr_lo, xerr_hi],
        fmt="o", color="steelblue", ecolor="steelblue",
        markersize=6, capsize=4, linewidth=1.4, zorder=3
    )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Pull  (f_inferred − f_prior) / σ_prior")
    ax.set_title(
        f"Slot fraction pulls — run {e['id'][:8]}  "
        f"χ²/N={e['chi2_per_n']:.4f}  stage={e['stage']}",
        fontsize=10
    )
    ax.set_xlim(min(-3.0, pulls.min() - 1.2), max(3.0, pulls.max() + 1.2))
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.2)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
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
