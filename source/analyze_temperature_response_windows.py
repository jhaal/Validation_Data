#!/usr/bin/env python3
"""Distinguish sampled, fixed-density fuel slopes from operating coefficients.

No transport is rerun. Uses the 1312 fuel-sweep states already in the benchmark.
Four unweighted regression windows are evaluated independently at every burnup.
Uncertainty includes both the existing independent-state propagation and a
first-order upper bound for arbitrary cross-state covariance. Neither measures
uncertainty in converting a finite-window fit into an operating derivative.
Usage: python source/analyze_temperature_response_windows.py [--plot]
The optional plot requires matplotlib; data generation uses the standard library.
"""
import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WINDOWS = (("cold", 293.15, 900.), ("full", 293.15, 2593.15),
           ("within_fuel_tabulation", 293.15, 1993.15),
           ("near_nominal", 793.15, 993.15))


def rho(k):
    return 1e5 * (1. - 1. / k)


def fit(rows, code):
    t = [r["tfuel_K"] for r in rows]
    average = sum(t) / len(t)
    denominator = sum((v - average) ** 2 for v in t)
    weights = [(v - average) / denominator for v in t]
    kkey = "openmc_keff" if code == "openmc" else "k_dragon"
    slope = sum(w * rho(r[kkey]) for w, r in zip(weights, rows))
    errors = [1e5 * r["openmc_keff_sd"] / r[kkey] ** 2 if code == "openmc" else 0.
              for r in rows]
    independent = math.sqrt(sum((w * s) ** 2 for w, s in zip(weights, errors)))
    bound = sum(abs(w) * s for w, s in zip(weights, errors))
    return slope, independent, bound


def load():
    dragon = {r["case_id"]: float(r["k_corrected"]) for r in csv.DictReader(
        (ROOT / "data/corrected_matrix_dragon_states.csv").open(newline=""))
        if r["variant"] == "corrected_refined_mesh_J2_matched_inventory"}
    sweeps = defaultdict(list)
    for raw in csv.DictReader((ROOT / "data/openmc_reference_states.csv").open(newline="")):
        if float(raw["tcool_K"]) == float(raw["tmod_K"]) == 343.15:
            row = {k: float(raw[k]) for k in ("burnup_GWd_tU", "tfuel_K", "openmc_keff", "openmc_keff_sd")}
            row["k_dragon"] = dragon[raw["case_id"]]
            sweeps[row["burnup_GWd_tU"]].append(row)
    for rows in sweeps.values():
        rows.sort(key=lambda r: r["tfuel_K"])
    assert len(sweeps) == 41 and all(len(r) == 32 for r in sweeps.values())
    return sweeps


def crossovers(rows, key):
    out = []
    for a, b in zip(rows, rows[1:]):
        if a[key] * b[key] < 0:
            out.append(a["burnup_GWd_tU"] - a[key] * (b["burnup_GWd_tU"] - a["burnup_GWd_tU"]) / (b[key] - a[key]))
    return out


def plot(sweeps, results):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({"font.size": 9, "font.family": "serif", "pdf.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.65), gridspec_kw={"width_ratios": [1.1, 1.]})
    colors = {"openmc": "#c4450a", "dragon": "#1f5fa8"}
    states = sweeps[10.]
    nominal = next(r for r in states if r["tfuel_K"] == 900.)
    ax = axes[0]
    ax.axvspan(793.15, 993.15, color="#d6eae6", alpha=.8)
    ax.axvspan(2000, 2650, color="#f5dfcf", alpha=.7)
    for code, kkey, marker, ls in (("openmc", "openmc_keff", "o", "none"),
                                   ("dragon", "k_dragon", None, "-")):
        x = [r["tfuel_K"] for r in states]
        y = [rho(r[kkey]) - rho(nominal[kkey]) for r in states]
        ax.plot(x, y, marker=marker, linestyle=ls, markersize=3, color=colors[code])
    ax.axhline(0, color=".5", lw=.6)
    ax.axvline(900, color=".35", lw=.7, ls=":")
    ax.text(.03, .97, "Near-nominal window: green\nDRAGLIB extrapolation: orange",
            transform=ax.transAxes, va="top", fontsize=7.5)
    ax.set(xlabel="Uniform fuel temperature (K)", ylabel=r"$\rho(T_f)-\rho(900\,\mathrm{K})$ (pcm)",
           title="(a) Same supplied inventory, 10 GWd/tU", xlim=(250, 2650))
    ax = axes[1]
    order = ["cold", "full", "within_fuel_tabulation", "near_nominal"]
    labels = ["293-900 K", "293-2593 K", "293-1993 K", "793-993 K"]
    for i, window in enumerate(order):
        row = next(r for r in results if r["burnup_GWd_tU"] == 10. and r["window"] == window)
        ax.errorbar(row["openmc_slope_pcm_K"], i-.10, xerr=row["openmc_sd_independent_pcm_K"],
                    fmt="o", color=colors["openmc"], ms=4, capsize=3)
        ax.plot(row["dragon_slope_pcm_K"], i+.10, "s", color=colors["dragon"], ms=4)
    ax.axvline(0, color=".4", lw=.7)
    ax.set_yticks(range(4), labels)
    ax.invert_yaxis()
    ax.set(xlabel="Fitted fixed-density slope (pcm/K)", title="(b) Different fits to the same response")
    ax.grid(axis="x", alpha=.2)
    fig.legend(handles=[Line2D([], [], marker="o", color=colors["openmc"], ls="none",
                    label=r"OpenMC ($\pm1\sigma$, independent-state propagation)"),
                    Line2D([], [], marker="s", color=colors["dragon"], ls="none", label="DRAGON")],
               loc="lower center", ncol=2, fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, .11, 1, 1), w_pad=1.7)
    for extension in ("pdf", "png"):
        fig.savefig(ROOT / f"figures/fig_response_windows.{extension}", dpi=220)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    sweeps = load()
    rows = []
    for b, states in sorted(sweeps.items()):
        for name, low, high in WINDOWS:
            chosen = [r for r in states if low <= r["tfuel_K"] <= high]
            omc, sd, upper = fit(chosen, "openmc")
            drg, _, _ = fit(chosen, "dragon")
            rows.append({"burnup_GWd_tU": b, "window": name, "Tmin_K": low, "Tmax_K": high,
                         "points": len(chosen), "openmc_slope_pcm_K": omc,
                         "openmc_sd_independent_pcm_K": sd,
                         "openmc_sd_upper_bound_pcm_K": upper,
                         "dragon_slope_pcm_K": drg,
                         "openmc_sign_at_2sd_independent": "positive" if omc > 2*sd else
                         ("negative" if omc < -2*sd else "unresolved"),
                         "fuel_points_above_2000K": sum(r["tfuel_K"] > 2000 for r in chosen)})
    with (ROOT / "data/temperature_response_windows.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    summary = {"scope": "Fixed-density fitted lattice responses; no operating coefficients calculated",
               "uncertainty": "Independent-state standard errors and arbitrary-covariance first-order upper bounds; excludes model-form uncertainty",
               "crossovers_are_point_estimates": True,
               "windows": {name: {
                   "openmc_point_estimate_crossovers_GWd_tU": crossovers([r for r in rows if r["window"] == name], "openmc_slope_pcm_K"),
                   "dragon_point_estimate_crossovers_GWd_tU": crossovers([r for r in rows if r["window"] == name], "dragon_slope_pcm_K")}
                   for name, _, _ in WINDOWS},
               "example_at_10_GWd_tU": [r for r in rows if r["burnup_GWd_tU"] == 10.]}
    (ROOT / "qa/temperature_response_windows.json").write_text(json.dumps(summary, indent=2) + "\n")
    if args.plot:
        plot(sweeps, rows)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
