#!/usr/bin/env python3
"""Compare the DRAGON branch matrices with the OpenMC reference (all 2132 states).

Reads data/corrected_matrix_dragon_states.csv (all DRAGON variants) and data/openmc_reference_states.csv.
Variant labels: J2_mesh_matched = benchmarked library (JEF-2.2, OpenMC inventories, 20-ring moderator);
J2 = same on the single-region moderator; E81, E81_mesh = ENDF/B-VIII.1 with its own depletion on the
single-region and 20-ring meshes; J2_mesh = JEF-2.2 depleted on the 20-ring mesh.
Delta k = k_OpenMC - k_DRAGON (pcm).
Writes data/corrected_matrix_summary.csv, data/corrected_matrix_coefficients.csv,
data/corrected_matrix_two_group.csv and data/corrected_matrix_states.csv.
Usage: python analyze_corrected_matrix.py
Standard library only.
"""
import csv
import math
import pathlib
import statistics as st
from collections import defaultdict

DATA = pathlib.Path(__file__).resolve().parent.parent / "data"
NOM = 343.15


def rho(k):
    return (1.0 - 1.0 / k) * 1e5


def lsq(points):
    """Unweighted slope of rho vs T with its standard error from OpenMC sd (points: T, k, sd)."""
    mt = sum(p[0] for p in points) / len(points)
    den = sum((p[0] - mt) ** 2 for p in points)
    b = sum((p[0] - mt) * rho(p[1]) for p in points) / den
    var = sum((p[0] - mt) ** 2 * (p[2] / p[1] ** 2 * 1e5) ** 2 for p in points) / den ** 2
    return b, math.sqrt(var)


def write(name, rows):
    with (DATA / name).open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    ref = {r["case_id"]: r for r in csv.DictReader((DATA / "openmc_reference_states.csv").open())}
    tal = {r["case_id"]: r for r in csv.DictReader((DATA / "openmc_two_group_rate_tallies.csv").open())}
    labels = {"corrected_J2": "J2", "corrected_E81": "E81", "corrected_refined_mesh_J2_own_depletion": "J2_mesh",
              "corrected_refined_mesh_E81": "E81_mesh",
              "corrected_refined_mesh_J2_matched_inventory": "J2_mesh_matched"}
    variants = {label: {} for label in labels.values()}
    for r in csv.DictReader((DATA / "corrected_matrix_dragon_states.csv").open()):
        variants[labels[r["variant"]]][r["case_id"]] = {
            "k": float(r["k_corrected"]), **{g: float(r[g]) for g in ("Sa1", "Sa2", "nuSf1", "nuSf2")}}
    for label, v in variants.items():
        assert set(v) == set(ref), f"{label}: state set differs from the OpenMC matrix"

    summary, coeffs, twog, states = [], [], [], []
    for label, v in variants.items():
        d, inside95, outside95 = [], 0, 0
        for c, r in ref.items():
            ko, sd = float(r["openmc_keff"]), float(r["openmc_keff_sd"]) * 1e5
            x = (ko - v[c]["k"]) * 1e5
            d.append((c, x, r))
            lo, hi = x - 1.96 * sd, x + 1.96 * sd
            inside95 += lo >= -1000 and hi <= 1000
            outside95 += lo > 1000 or hi < -1000
        a = [abs(x) for _, x, _ in d]
        nominal = [x for _, x, r in d if float(r["tfuel_K"]) == 900 and float(r["tcool_K"]) == float(r["tmod_K"]) == NOM]
        row = {"variant": label, "states": len(d), "mean_pcm": st.mean(x for _, x, _ in d),
               "mean_abs_pcm": st.mean(a), "median_abs_pcm": st.median(a),
               "rms_pcm": math.sqrt(st.mean(x * x for x in a)), "max_abs_pcm": max(a),
               "nominal_min_pcm": min(nominal), "nominal_max_pcm": max(nominal)}
        for lim in (200, 300, 500, 700, 1000):
            row[f"pass_{lim}"] = sum(x <= lim for x in a)
        row["interval_inside_1000"], row["interval_outside_1000"] = inside95, outside95
        row["borderline_1000"] = len(d) - inside95 - outside95
        # Temperature-dependent part: Delta k minus the nominal-state Delta k at the same burnup
        nom_by_bu = {r["burnup_GWd_tU"]: x for _, x, r in d
                     if float(r["tfuel_K"]) == 900 and float(r["tcool_K"]) == float(r["tmod_K"]) == NOM}
        dev = [x - nom_by_bu[r["burnup_GWd_tU"]] for _, x, r in d]
        row["temp_dependent_rms_pcm"] = math.sqrt(st.mean(e * e for e in dev))
        row["temp_dependent_max_abs_pcm"] = max(map(abs, dev))
        summary.append(row)
        for c, x, r in d:
            states.append({"variant": label, "case_id": c, "k_dragon": v[c]["k"], "delta_k_pcm": x})

        # Range-fitted coefficients at each burnup (fuel: other temperatures nominal; coolant, moderator likewise)
        groups = defaultdict(lambda: defaultdict(list))
        for c, r in ref.items():
            bu, tf, tc, tm = (float(r[k]) for k in ("burnup_GWd_tU", "tfuel_K", "tcool_K", "tmod_K"))
            ko, sd = float(r["openmc_keff"]), float(r["openmc_keff_sd"])
            if tc == tm == NOM:
                groups[bu]["fuel"].append((tf, v[c]["k"], ko, sd))
            if tf == 900 and tm == NOM:
                groups[bu]["cool"].append((tc, v[c]["k"], ko, sd))
            if tf == 900 and tc == NOM:
                groups[bu]["mod"].append((tm, v[c]["k"], ko, sd))
        for bu in sorted(groups):
            g = groups[bu]
            for name, kind, lo, hi in (("FTC_293_900", "fuel", 293, 900), ("FTC_293_2593", "fuel", 293, 2600),
                                       ("FTC_793_993", "fuel", 793, 994), ("CTC_293_363", "cool", 293, 364),
                                       ("MTC_293_363", "mod", 293, 364)):
                p = [q for q in g[kind] if lo <= q[0] <= hi]
                if len(p) < 3:
                    continue
                bd = lsq([(t, kd, 0.0) for t, kd, _, _ in p])[0]
                bo, so = lsq([(t, ko, sd) for t, _, ko, sd in p])
                coeffs.append({"variant": label, "burnup_GWd_tU": bu, "coefficient": name, "points": len(p),
                               "dragon_pcm_K": bd, "openmc_pcm_K": bo, "openmc_sd_pcm_K": so,
                               "z": (bd - bo) / so})
        # Two-group constants vs OpenMC whole-cell tallies
        for m in ("Sa1", "Sa2", "nuSf1", "nuSf2", "nuSf2_over_Sa2"):
            rel = []
            for c in ref:
                t = tal[c]
                o = (float(t["nuSf2_openmc_rate"]) / float(t["Sa2_openmc_rate"]) if m == "nuSf2_over_Sa2"
                     else float(t[f"{m}_openmc_rate"]))
                x = v[c]["nuSf2"] / v[c]["Sa2"] if m == "nuSf2_over_Sa2" else v[c][m]
                rel.append(o / x - 1.0)
            twog.append({"variant": label, "metric": m, "rows": len(rel),
                         "mean_rel_pct": 100 * st.mean(rel), "min_rel_pct": 100 * min(rel),
                         "max_rel_pct": 100 * max(rel), "pass_5pct": sum(abs(e) <= 0.05 for e in rel),
                         "pass_2pct": sum(abs(e) <= 0.02 for e in rel)})
    for row in summary:
        for name in ("FTC_293_900", "FTC_293_2593", "FTC_793_993", "MTC_293_363", "CTC_293_363"):
            c = [r for r in coeffs if r["variant"] == row["variant"] and r["coefficient"] == name]
            row[f"{name}_max_abs_diff_pcm_K"] = max(abs(r["dragon_pcm_K"] - r["openmc_pcm_K"]) for r in c)
            row[f"{name}_rms_z"] = math.sqrt(st.mean(r["z"] ** 2 for r in c))
    write("corrected_matrix_summary.csv", summary)
    write("corrected_matrix_coefficients.csv", coeffs)
    write("corrected_matrix_states.csv", states)
    if twog:
        write("corrected_matrix_two_group.csv", twog)
    for r in summary:
        print({k: (round(x, 1) if isinstance(x, float) else x) for k, x in r.items()})
    for r in twog:
        print({k: (round(x, 2) if isinstance(x, float) else x) for k, x in r.items()})


if __name__ == "__main__":
    main()
