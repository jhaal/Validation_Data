#!/usr/bin/env python3
"""Plot completed operating central-difference estimates only (no partial data)."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]


def main():
    summary = json.loads((ROOT / "qa/operating_analysis.json").read_text())
    if not summary["complete"]: raise SystemExit("Operating matrix is incomplete")
    rows = list(csv.DictReader((ROOT / "data/operating_coefficients.csv").open(newline="")))
    assert len(rows) == 30
    plt.rcParams.update({"font.family":"serif", "font.size":9, "pdf.fonttype":42})
    fig, axes = plt.subplots(1, 3, figsize=(8.2, 3.35))
    colors = {"dragon":"#1f5fa8", "openmc":"#c4450a"}
    burnups = [0., 7.5, 20.]
    for ax, variable, title, step in zip(axes, ["fuel","coolant","moderator"],
                                        ["(a) Fuel", "(b) Coolant", "(c) Moderator"], [100.,10.,10.]):
        for i, b in enumerate(burnups):
            for mode in (["fixed"] if variable=="fuel" else ["fixed","pressure"]):
                row = next(r for r in rows if float(r["burnup_GWd_tU"])==b and r["variable"]==variable
                           and r["density_mode"]==mode and float(r["step_K"])==step)
                offset = 0 if variable=="fuel" else (-.12 if mode=="fixed" else .12)
                face = "white" if mode=="fixed" and variable!="fuel" else None
                ax.plot(i+offset-.035, float(row["dragon_pcm_K"]), marker="s", color=colors["dragon"],
                        markerfacecolor=face or colors["dragon"], markersize=4.5, linestyle="none")
                ax.errorbar(i+offset+.035, float(row["openmc_pcm_K"]), yerr=float(row["openmc_sd_pcm_K"]),
                            marker="o", color=colors["openmc"], markerfacecolor=face or colors["openmc"],
                            markersize=4.5, linestyle="none", capsize=3, elinewidth=.8)
        ax.axhline(0, color=".4", lw=.7)
        ax.set_xticks(range(3), ["0", "7.5", "20"])
        ax.set(xlim=(-.4, 2.4), xlabel="Burnup (GWd/tU)", ylabel="Central difference (pcm/K)", title=title)
        ax.grid(axis="y", alpha=.15)
    handles = [Line2D([],[],color=colors["dragon"],marker="s",ls="none",label="DRAGON"),
               Line2D([],[],color=colors["openmc"],marker="o",ls="none",label=r"OpenMC ($\pm1\sigma$)"),
               Line2D([],[],color=".3",marker="o",mfc="white",ls="none",label="Fluid: fixed density"),
               Line2D([],[],color=".3",marker="o",ls="none",label="Fluid: constant pressure")]
    fig.legend(handles=handles, ncol=4, loc="lower center", frameon=False, fontsize=8)
    fig.tight_layout(rect=(0,.12,1,1), w_pad=1.2)
    for ext in ("pdf","png"):
        fig.savefig(ROOT / f"figures/fig_operating_coefficients.{ext}", dpi=250)
    plt.close(fig)


if __name__ == "__main__": main()
