"""Cross-section of the 37-element PHWR lattice cell (Fig. geometry).

Dimensions are those of the production DRAGON (implement_phase2_dragon_37el.py)
and OpenMC (phase_f_openmc_candu_benchmark.py) models; see Table "tab:geometry".
Usage: python make_fig_geometry.py  -> ../figures/fig_geometry.{pdf,png}
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Patch, Rectangle

PITCH = 28.575
R_COOL, R_PT, R_GAP, R_CT = 5.1689, 5.6032, 6.4478, 6.5875
R_ZONES = (0.3061, 0.4898, 0.5816, 0.6122)
R_SHEATH = 0.6540
RINGS = ((0.0, 1), (1.4880, 6), (2.8750, 12), (4.3300, 18))

C = {
    "moderator": "#cfe3f3",
    "coolant": "#9cc3e4",
    "pt": "#5a5a5a",
    "gap": "#ffffff",
    "ct": "#9a9a9a",
    "sheath": "#3a3a3a",
    "fuel": ("#7f2704", "#c4450a", "#ef7d2d", "#fdb36b"),
}
EDGE = dict(edgecolor="black", linewidth=0.25)


def element_centres():
    for radius, count in RINGS:
        for k in range(count):
            a = 2.0 * np.pi * k / count
            yield radius * np.cos(a), radius * np.sin(a)


def draw_element(ax, x, y, lw=0.0):
    ax.add_patch(Circle((x, y), R_SHEATH, facecolor=C["sheath"], linewidth=lw, edgecolor="black"))
    for r, col in zip(reversed(R_ZONES), reversed(C["fuel"])):
        ax.add_patch(Circle((x, y), r, facecolor=col, linewidth=lw, edgecolor="black"))


def draw_cell(ax, element_lw=0.0):
    h = PITCH / 2.0
    ax.add_patch(Rectangle((-h, -h), PITCH, PITCH, facecolor=C["moderator"], edgecolor="black", linewidth=1.0))
    for r, key in ((R_CT, "ct"), (R_GAP, "gap"), (R_PT, "pt"), (R_COOL, "coolant")):
        ax.add_patch(Circle((0, 0), r, facecolor=C[key], **EDGE))
    for x, y in element_centres():
        draw_element(ax, x, y, lw=element_lw)


def main():
    plt.rcParams.update({"font.size": 8, "font.family": "serif", "mathtext.fontset": "dejavuserif",
                         "pdf.fonttype": 42})  # embed TrueType, not Type 3 bitmaps
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    draw_cell(ax)
    h = PITCH / 2.0
    ax.set_xlim(-h, h)
    ax.set_ylim(-h, h)
    ax.set_aspect("equal")
    ax.set_xlabel("x (cm)")
    ax.set_ylabel("y (cm)")
    labels = [f"−{h:.3f}", "−7", "0", "7", f"{h:.3f}"]
    ax.set_xticks([-h, -7, 0, 7, h], labels)
    ax.set_yticks([-h, -7, 0, 7, h], labels)
    ax.text(0, -h + 0.8, "reflective boundary, pitch 28.575 cm", ha="center", va="bottom", fontsize=7)

    # Inset 1: channel and cluster
    ins1 = ax.inset_axes([1.04, 0.50, 0.48, 0.48])
    draw_cell(ins1, element_lw=0.15)
    ins1.set_xlim(-6.9, 6.9)
    ins1.set_ylim(-6.9, 6.9)
    ins1.set_aspect("equal")
    ins1.set_xticks([])
    ins1.set_yticks([])
    for r in (1.4880, 2.8750, 4.3300):
        ins1.add_patch(Circle((0, 0), r, fill=False, linestyle=(0, (2, 2)), linewidth=0.5, edgecolor="white"))
    ins1.text(0.5, -0.03, "(b) fuel channel: 1 + 6 + 12 + 18 elements", transform=ins1.transAxes,
              ha="center", va="top", fontsize=7)
    ax.indicate_inset_zoom(ins1, edgecolor="black", linewidth=0.6)

    # Inset 2: single outer-ring element, zones labelled with leader lines
    xe, ye = 4.3300, 0.0
    ins2 = ax.inset_axes([1.04, 0.00, 0.48, 0.40])
    ins2.add_patch(Rectangle((xe - 1, ye - 1), 3, 2, facecolor=C["coolant"], linewidth=0))
    draw_element(ins2, xe, ye, lw=0.3)
    ins2.set_xlim(xe - 0.72, xe + 1.28)
    ins2.set_ylim(ye - 0.72, ye + 0.72)
    ins2.set_aspect("equal")
    ins2.set_xticks([])
    ins2.set_yticks([])
    mids = [R_ZONES[0] / 2] + [(R_ZONES[i - 1] + R_ZONES[i]) / 2 for i in range(1, 4)] + [(R_ZONES[3] + R_SHEATH) / 2]
    names = ["zone 1", "zone 2", "zone 3", "zone 4", "sheath"]
    for i, (rm, name) in enumerate(zip(mids, names)):
        ang = np.deg2rad(50 - 25 * i)
        ytxt = ye + 0.55 - 0.27 * i
        ins2.annotate(name, xy=(xe + rm * np.cos(ang), ye + rm * np.sin(ang)), xytext=(xe + 0.86, ytxt),
                      fontsize=6, va="center",
                      arrowprops=dict(arrowstyle="-", linewidth=0.5, color="black", shrinkA=0, shrinkB=0))
    ins2.text(0.5, -0.03, "(c) fuel element", transform=ins2.transAxes, ha="center", va="top", fontsize=7)
    ind = ins1.indicate_inset_zoom(ins2, edgecolor="black", linewidth=0.6)
    for conn in ind.connectors:  # box only; connectors would cross the (b) label
        conn.set_visible(False)
    ax.set_title("(a) lattice cell", fontsize=7, pad=3)

    handles = [
        Patch(facecolor=C["fuel"][1], **EDGE, label="Natural UO$_2$ (zones 1–4)"),
        Patch(facecolor=C["sheath"], **EDGE, label="Zircaloy sheath"),
        Patch(facecolor=C["coolant"], **EDGE, label="D$_2$O coolant"),
        Patch(facecolor=C["pt"], **EDGE, label="Zr–Nb pressure tube"),
        Patch(facecolor=C["gap"], **EDGE, label="He annulus gap"),
        Patch(facecolor=C["ct"], **EDGE, label="Zircaloy calandria tube"),
        Patch(facecolor=C["moderator"], **EDGE, label="D$_2$O moderator"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=7, bbox_to_anchor=(0.5, -0.01))
    fig.subplots_adjust(left=0.10, right=0.62, top=0.97, bottom=0.20)

    out = Path(__file__).resolve().parent.parent / "figures"
    fig.savefig(out / "fig_geometry.pdf", bbox_inches="tight")
    fig.savefig(out / "fig_geometry.png", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
