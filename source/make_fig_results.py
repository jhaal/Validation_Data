"""Result figures 2-3: benchmarked DRAGON library versus OpenMC.

DRAGON eigenvalues: variant corrected_refined_mesh_J2_matched_inventory of
../data/corrected_matrix_dragon_states.csv (JEF-2.2 DRAGLIB, create-mode branches, 20-ring moderator,
OpenMC inventories). OpenMC eigenvalues and uncertainties: ../data/openmc_reference_states.csv.

fig_delta_map     : k_OpenMC - k_DRAGON over the fuel-temperature sweep.
fig_coefficients  : fitted fixed-density fuel (two ranges), moderator and coolant slopes.

Coefficients are least-squares slopes of rho = (k-1)/k (pcm) against temperature at fixed burnup,
with the other two temperatures nominal (Tf = 900 K, Tc = Tm = 343.15 K). The error model is Eq. (2)
with full-range FTCs; model reactivity is converted back to an OpenMC eigenvalue, k = 1/(1 - rho).
Numbers agree with source/analyze_benchmark_library.py.
Usage: python make_fig_results.py
"""
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data" / "openmc_reference_states.csv"
REGEN = HERE.parent / "data" / "corrected_matrix_dragon_states.csv"
REGEN_VARIANT = "corrected_refined_mesh_J2_matched_inventory"
OUT = HERE.parent / "figures"
TF_NOM, TC_NOM = 900.0, 343.15
SCREEN = 300.0
DISCHARGE = 7.5
C_OMC, C_DRG = "#c4450a", "#1f5fa8"

plt.rcParams.update({"font.size": 8, "font.family": "serif", "mathtext.fontset": "dejavuserif",
                     "pdf.fonttype": 42, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6})


def load():
    regen = {r["case_id"]: float(r["k_corrected"]) for r in csv.DictReader(REGEN.open(newline=""))
             if r["variant"] == REGEN_VARIANT}
    rows = []
    for r in csv.DictReader(DATA.open(newline="")):
        d = {k: float(r[k]) for k in ("burnup_GWd_tU", "tfuel_K", "tcool_K", "tmod_K",
                                      "openmc_keff", "openmc_keff_sd")}
        d["dragon_k_eff"] = regen[r["case_id"]]
        rows.append(d)
    assert len(rows) == len(regen) == 2132
    return rows


def near(a, b):
    return abs(a - b) < 1e-6


def rho(k):
    return (k - 1.0) / k * 1e5


def slope(x, y, sy=None):
    x, y = np.asarray(x), np.asarray(y)
    dx = x - x.mean()
    s = float(np.sum(dx * (y - y.mean())) / np.sum(dx * dx))
    if sy is None:
        return s, 0.0
    return s, float(math.sqrt(np.sum(dx * dx * np.asarray(sy) ** 2)) / np.sum(dx * dx))


def sweeps(rows):
    sel = {
        "fuel": lambda r: near(r["tcool_K"], TC_NOM) and near(r["tmod_K"], TC_NOM),
        "cool": lambda r: near(r["tfuel_K"], TF_NOM) and near(r["tmod_K"], TC_NOM),
        "mod": lambda r: near(r["tfuel_K"], TF_NOM) and near(r["tcool_K"], TC_NOM),
    }
    out = {}
    for name, f in sel.items():
        by = defaultdict(list)
        for r in rows:
            if f(r):
                by[r["burnup_GWd_tU"]].append(r)
        out[name] = by
    return out


def coefficients(by, tkey, tmax=1e9):
    bu, drg, omc, somc = [], [], [], []
    for b in sorted(by):
        v = [r for r in by[b] if r[tkey] < tmax]
        t = [r[tkey] for r in v]
        sd, _ = slope(t, [rho(r["dragon_k_eff"]) for r in v])
        so, ss = slope(t, [rho(r["openmc_keff"]) for r in v],
                       [r["openmc_keff_sd"] / r["openmc_keff"] ** 2 * 1e5 for r in v])
        bu.append(b), drg.append(sd), omc.append(so), somc.append(ss)
    return np.array(bu), np.array(drg), np.array(omc), np.array(somc)


def error_model(fuel):
    recs = []
    for b, v in fuel.items():
        t = [r["tfuel_K"] for r in v]
        so, _ = slope(t, [rho(r["openmc_keff"]) for r in v])
        sd, _ = slope(t, [rho(r["dragon_k_eff"]) for r in v])
        nom = next(r for r in v if near(r["tfuel_K"], TF_NOM))
        d0 = rho(nom["openmc_keff"]) - rho(nom["dragon_k_eff"])
        for r in v:
            dr = (d0 + (so - sd) * (r["tfuel_K"] - TF_NOM)) / 1e5
            k_pred = 1.0 / (1.0 - (rho(r["dragon_k_eff"]) / 1e5 + dr))
            recs.append((b, r["tfuel_K"], (r["openmc_keff"] - r["dragon_k_eff"]) * 1e5,
                         (k_pred - r["dragon_k_eff"]) * 1e5))
    return np.array(recs)


def fig_delta_map(em):
    bu = np.unique(em[:, 0])
    tf = np.unique(em[:, 1])
    obs = np.full((tf.size, bu.size), np.nan)
    for b, t, o, _ in em:
        obs[np.searchsorted(tf, t), np.searchsorted(bu, b)] = o
    fig, ax = plt.subplots(figsize=(4.6, 3.0))
    vmax = 500
    B, T = np.meshgrid(bu, tf)
    levels = np.arange(-vmax, vmax + 1, 50)
    box = dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.8)
    pc = ax.contourf(B, T, obs, levels=levels, cmap="RdBu_r", norm=TwoSlopeNorm(0, -vmax, vmax), extend="both")
    ax.contour(B, T, obs, levels=levels, colors="white", linewidths=0.25, alpha=0.6)
    ax.contour(B, T, np.abs(obs), levels=[SCREEN], colors="black", linewidths=1.0)
    ax.axhline(TF_NOM, color="black", linewidth=0.6, linestyle=":")
    ax.axvline(DISCHARGE, color="0.35", linewidth=0.6, linestyle=(0, (1, 2)))
    ax.text(DISCHARGE + 0.2, 330, "discharge", fontsize=6.5, color="0.2", va="bottom", bbox=box)
    ax.text(19.7, TF_NOM + 40, "nominal $T_f$", fontsize=6.5, ha="right", va="bottom", bbox=box)
    ax.set_xticks([0, 5, 10, 15, 20])
    ax.set_xlim(bu[0], bu[-1])
    ax.set_ylim(tf[0], tf[-1])
    ax.set_xlabel("Burnup (GWd/tU)")
    ax.set_ylabel("Fuel temperature $T_f$ (K)")
    cb = fig.colorbar(pc, ax=ax, pad=0.03, fraction=0.05, ticks=np.arange(-vmax, vmax + 1, 100))
    cb.set_label("$k_{\\mathrm{OpenMC}}-k_{\\mathrm{DRAGON}}$ (pcm)")
    for y in (-SCREEN, SCREEN):
        cb.ax.axhline(y, color="black", linewidth=1.0)
    cb.outline.set_linewidth(0.5)
    fig.legend(handles=[Line2D([], [], color="black", linewidth=1.0, label="$|\\Delta k|=300$ pcm")],
               loc="lower center", frameon=False, fontsize=7, bbox_to_anchor=(0.45, -0.04))
    fig.subplots_adjust(bottom=0.24)
    save(fig, "fig_delta_map")


def fig_coefficients(sw):
    fig, axs = plt.subplots(1, 3, figsize=(6.6, 2.4))
    panels = [
        (r"(a) Fuel ($\widehat\beta_f$)", [(sw["fuel"], "tfuel_K", 1e9, "-", True), (sw["fuel"], "tfuel_K", 901, "--", False)]),
        (r"(b) Moderator ($\widehat\beta_m$)", [(sw["mod"], "tmod_K", 1e9, "-", True)]),
        (r"(c) Coolant ($\widehat\beta_c$)", [(sw["cool"], "tcool_K", 1e9, "-", True)]),
    ]
    for ax, (title, series) in zip(axs, panels):
        for by, key, tmax, ls, filled in series:
            bu, drg, omc, s = coefficients(by, key, tmax)
            ax.fill_between(bu, omc - s, omc + s, color=C_OMC, alpha=0.18, linewidth=0)
            ax.plot(bu, omc, linestyle="none", marker="o", markersize=2.8, markeredgewidth=0.6,
                    color=C_OMC, markerfacecolor=C_OMC if filled else "white")
            ax.plot(bu, drg, linestyle=ls, color=C_DRG, linewidth=1.1)
        ax.axhline(0, color="black", linewidth=0.5)
        ax.axvline(DISCHARGE, color="0.45", linewidth=0.6, linestyle=(0, (1, 2)))
        ax.set_title(title, fontsize=8, pad=3)
        ax.set_xlim(0, 20)
        ax.set_xticks([0, 5, 10, 15, 20])
        ax.set_xlabel("Burnup (GWd/tU)")
        ax.set_ylabel("pcm/K")
    axs[0].text(DISCHARGE + 0.3, axs[0].get_ylim()[0] + 0.05, "discharge", fontsize=6.5, color="0.3",
                va="bottom")
    handles = [
        Line2D([], [], color=C_OMC, marker="o", linestyle="none", markersize=3.5,
               label="OpenMC (band: $\\pm1\\sigma$)"),
        Line2D([], [], color=C_DRG, linewidth=1.1, label="DRAGON"),
        Line2D([], [], color="0.2", marker="o", markerfacecolor="0.2", linestyle="-", linewidth=0.8,
               markersize=3, label="Fuel fit 293-2593 K"),
        Line2D([], [], color="0.2", marker="o", markerfacecolor="white", linestyle="--", linewidth=0.8,
               markersize=3, label="Fuel fit 293-900 K"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False, fontsize=6.5,
               bbox_to_anchor=(0.5, 1.06))
    fig.tight_layout(w_pad=1.2)
    fig.subplots_adjust(top=0.82)
    save(fig, "fig_coefficients")


def save(fig, name):
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight", dpi=400)
    fig.savefig(OUT / f"{name}.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)
    rows = load()
    sw = sweeps(rows)
    em = error_model(sw["fuel"])
    fig_delta_map(em)
    fig_coefficients(sw)


if __name__ == "__main__":
    main()
