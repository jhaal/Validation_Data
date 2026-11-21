"""Derive graphical-abstract plots for the benchmarked DRAGON library.

DRAGON: variant corrected_refined_mesh_J2_matched_inventory of data/corrected_matrix_dragon_states.csv;
OpenMC: data/openmc_reference_states.csv. Run this script before compiling graphical_abstract.tex.
Fuel coefficients are least-squares slopes of reactivity over 293--2593 K (coolant and moderator at
343.15 K). Map cells use midpoint boundaries around the 1312 sampled fuel-sweep states and three bands
of |Delta k|: <= 100, 100--300 and > 300 pcm; the drawn boundary is |Delta k| = 300 pcm.
"""
import csv
import math
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / "non_upload" / "graphical_abstract_build"
VARIANT = "corrected_refined_mesh_J2_matched_inventory"
INNER, LIMIT = 100.0, 300.0


def rho(k):
    return (1.0 - 1.0 / k) * 1e5


def main():
    dragon = {r["case_id"]: float(r["k_corrected"])
              for r in csv.DictReader((ROOT / "data" / "corrected_matrix_dragon_states.csv").open(newline=""))
              if r["variant"] == VARIANT}
    rows = list(csv.DictReader((ROOT / "data" / "openmc_reference_states.csv").open(newline="")))
    assert len(rows) == len(dragon) == 2132
    for r in rows:
        r["k_dragon"] = dragon[r["case_id"]]
        r["k_openmc"] = float(r["openmc_keff"])
        r["dk"] = (r["k_openmc"] - r["k_dragon"]) * 1e5
    rms = math.sqrt(sum(r["dk"] ** 2 for r in rows) / len(rows))
    inside = sum(abs(r["dk"]) <= LIMIT for r in rows)
    assert round(rms) == 117 and inside == 2070
    fuel = defaultdict(list)
    for r in rows:
        if abs(float(r["tcool_K"]) - 343.15) < 1e-6 and abs(float(r["tmod_K"]) - 343.15) < 1e-6:
            fuel[float(r["burnup_GWd_tU"])].append(r)
    burnups = sorted(fuel)
    temperatures = sorted(float(r["tfuel_K"]) for r in fuel[burnups[0]])
    assert len(burnups) == 41 and len(temperatures) == 32
    BUILD.mkdir(parents=True, exist_ok=True)

    coefficients = ["burnup dragon openmc"]
    for b in burnups:
        states = fuel[b]
        mean_t = sum(float(r["tfuel_K"]) for r in states) / len(states)
        den = sum((float(r["tfuel_K"]) - mean_t) ** 2 for r in states)
        s = [sum((float(r["tfuel_K"]) - mean_t) * rho(r[k]) for r in states) / den for k in ("k_dragon", "k_openmc")]
        coefficients.append(f"{b:.6f} {s[0]:.9f} {s[1]:.9f}")
    (BUILD / "fuel_coefficients.dat").write_text("\n".join(coefficients) + "\n", encoding="ascii")

    operating = list(csv.DictReader((ROOT / "data/operating_coefficients.csv").open(newline="")))
    chosen = sorted((r for r in operating if r["variable"]=="fuel" and float(r["step_K"])==100), key=lambda r:float(r["burnup_GWd_tU"]))
    assert len(chosen)==3
    lines=["burnup dragon openmc sd"]
    for r in chosen:
        lines.append(" ".join(r[k] for k in ("burnup_GWd_tU","dragon_pcm_K","openmc_pcm_K","openmc_sd_pcm_K")))
    (BUILD / "operating_fuel.dat").write_text("\n".join(lines)+"\n")

    # Physical coordinates in cm, inside the map's TikZ scope.
    x = lambda b: b / 20.0 * 7.25
    y = lambda t: (t - 293.15) / 2300.0 * 3.5
    b_edges = [burnups[0]] + [(a + c) / 2 for a, c in zip(burnups, burnups[1:])] + [burnups[-1]]
    t_edges = [temperatures[0]] + [(a + c) / 2 for a, c in zip(temperatures, temperatures[1:])] + [temperatures[-1]]
    outside, commands = {}, []
    for i, b in enumerate(burnups):
        by_t = {float(r["tfuel_K"]): r for r in fuel[b]}
        for j, t in enumerate(temperatures):
            d = abs(by_t[t]["dk"])
            outside[i, j] = d > LIMIT
            color = "bandlow" if d <= INNER else ("bandmid" if d <= LIMIT else "bandhigh")
            commands.append(f"\\fill[{color}] ({x(b_edges[i]):.5f},{y(t_edges[j]):.5f}) "
                            f"rectangle ({x(b_edges[i + 1]):.5f},{y(t_edges[j + 1]):.5f});")
    for i in range(len(burnups)):
        for j in range(len(temperatures)):
            if i < len(burnups) - 1 and outside[i, j] != outside[i + 1, j]:
                commands.append(f"\\draw[boundary] ({x(b_edges[i + 1]):.5f},{y(t_edges[j]):.5f}) "
                                f"-- ({x(b_edges[i + 1]):.5f},{y(t_edges[j + 1]):.5f});")
            if j < len(temperatures) - 1 and outside[i, j] != outside[i, j + 1]:
                commands.append(f"\\draw[boundary] ({x(b_edges[i]):.5f},{y(t_edges[j + 1]):.5f}) "
                                f"-- ({x(b_edges[i + 1]):.5f},{y(t_edges[j + 1]):.5f});")
    (BUILD / "validity_cells.tex").write_text("\n".join(commands) + "\n", encoding="ascii")
    print(f"Verified 2132 states: RMS {rms:.0f} pcm, {inside} within {LIMIT:.0f} pcm; "
          f"wrote {len(burnups) * len(temperatures)} map cells and {len(burnups)} coefficient pairs.")


if __name__ == "__main__":
    main()
