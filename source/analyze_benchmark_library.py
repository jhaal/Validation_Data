#!/usr/bin/env python3
"""Results of the benchmarked DRAGON library against the OpenMC reference (Section 3 of the manuscript).

Library: variant corrected_refined_mesh_J2_matched_inventory of data/corrected_matrix_dragon_states.csv
(JEF-2.2 DRAGLIB, create-mode branches, 20-ring moderator, OpenMC inventories). Reference: OpenMC
eigenvalues in data/openmc_reference_states.csv and whole-cell two-group tallies in
data/openmc_two_group_rate_tallies.csv. Delta k = 1e5 (k_OpenMC - k_DRAGON).

Writes data/benchmark_library_states.csv, _families.csv, _thresholds.csv, _nominal.csv,
_coefficients.csv, _error_model.csv, _two_group.csv and qa/benchmark_library_analysis.json.
Standard library only.
"""
import csv
import json
import math
import statistics as st
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
VARIANT = "corrected_refined_mesh_J2_matched_inventory"
NOM_TF, NOM_T = 900.0, 343.15
SCREEN = 300.0
LIMITS = (100, 200, 300, 500)


def write(name, rows):
    with (DATA / name).open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def rho(k):
    return (1.0 - 1.0 / k) * 1e5


def slope(pts):
    """pts: (T, k, sd_k). Unweighted LSQ slope of rho vs T and its propagated sd."""
    mt = sum(p[0] for p in pts) / len(pts)
    den = sum((p[0] - mt) ** 2 for p in pts)
    b = sum((p[0] - mt) * rho(p[1]) for p in pts) / den
    var = sum((p[0] - mt) ** 2 * (p[2] / p[1] ** 2 * 1e5) ** 2 for p in pts) / den ** 2
    return b, math.sqrt(var)


def stats(xs):
    a = [abs(x) for x in xs]
    return {"states": len(xs), "mean_pcm": st.mean(xs), "mean_abs_pcm": st.mean(a),
            "rms_pcm": math.sqrt(st.mean(x * x for x in xs)), "max_abs_pcm": max(a),
            f"inside_{int(SCREEN)}": sum(x <= SCREEN for x in a)}


def family(r):
    tf, tc, tm = r["tfuel_K"], r["tcool_K"], r["tmod_K"]
    nf, nc, nm = tf == NOM_TF, tc == NOM_T, tm == NOM_T
    if nf and nc and nm:
        return "nominal"
    if nc and nm:
        return "fuel"
    if nf and nm:
        return "coolant"
    if nf and nc:
        return "moderator"
    return "corner"


def tf_region(tf):
    if tf <= 900:
        return "Tf<=900"
    if tf <= 1593.15:
        return "900<Tf<=1593"
    if tf <= 1993.15:
        return "1593<Tf<=1993"
    return "Tf>1993"


def main():
    rows = []
    for r in csv.DictReader((DATA / "openmc_reference_states.csv").open()):
        rows.append({"case_id": r["case_id"], **{k: float(r[k]) for k in
                     ("burnup_GWd_tU", "tfuel_K", "tcool_K", "tmod_K", "openmc_keff", "openmc_keff_sd")}})
    lib = {r["case_id"]: r for r in csv.DictReader((DATA / "corrected_matrix_dragon_states.csv").open())
           if r["variant"] == VARIANT}
    tal = {r["case_id"]: r for r in csv.DictReader((DATA / "openmc_two_group_rate_tallies.csv").open())}
    assert len(rows) == len(lib) == 2132
    for r in rows:
        r["k_dragon"] = float(lib[r["case_id"]]["k_corrected"])
        r["dk"] = (r["openmc_keff"] - r["k_dragon"]) * 1e5
        r["sigma"] = r["openmc_keff_sd"] * 1e5
        r["family"] = family(r)
    nominal = {r["burnup_GWd_tU"]: r for r in rows if r["family"] == "nominal"}
    for r in rows:
        r["dk_temp"] = r["dk"] - nominal[r["burnup_GWd_tU"]]["dk"]
        lo, hi = r["dk"] - 1.96 * r["sigma"], r["dk"] + 1.96 * r["sigma"]
        r["ci95"] = ("inside" if lo >= -SCREEN and hi <= SCREEN else
                     "outside" if lo > SCREEN or hi < -SCREEN else "borderline")
    write("benchmark_library_states.csv",
          [{k: r[k] for k in ("case_id", "burnup_GWd_tU", "tfuel_K", "tcool_K", "tmod_K", "k_dragon",
                              "openmc_keff", "openmc_keff_sd", "dk", "dk_temp", "family", "ci95")} for r in rows])

    fam = []
    for name in ("nominal", "coolant", "moderator", "fuel", "corner"):
        fam.append({"group": name, **stats([r["dk"] for r in rows if r["family"] == name])})
    for reg in ("Tf<=900", "900<Tf<=1593", "1593<Tf<=1993", "Tf>1993"):
        fam.append({"group": reg, **stats([r["dk"] for r in rows if tf_region(r["tfuel_K"]) == reg])})
    for name, f in (("BU<=2.5", lambda b: b <= 2.5), ("BU<=7.5", lambda b: b <= 7.5),
                    ("0.5<=BU<=7.5", lambda b: 0 < b <= 7.5),
                    ("BU>=2.5", lambda b: b >= 2.5), ("BU>7.5", lambda b: b > 7.5), ("BU=0", lambda b: b == 0)):
        fam.append({"group": name, **stats([r["dk"] for r in rows if f(r["burnup_GWd_tU"])])})
    fam.append({"group": "all", **stats([r["dk"] for r in rows])})
    fam.append({"group": "temperature_dependent_all", **stats([r["dk_temp"] for r in rows])})
    write("benchmark_library_families.csv", fam)

    thr = [{"limit_pcm": L, "nominal_inside": sum(abs(r["dk"]) <= L for r in nominal.values()),
            "all_inside": sum(abs(r["dk"]) <= L for r in rows),
            "all_inside_BU_ge_2p5": sum(abs(r["dk"]) <= L for r in rows if r["burnup_GWd_tU"] >= 2.5)}
           for L in LIMITS]
    write("benchmark_library_thresholds.csv", thr)
    write("benchmark_library_nominal.csv",
          [{"burnup_GWd_tU": b, "k_openmc": r["openmc_keff"], "k_dragon": r["k_dragon"], "dk_pcm": r["dk"]}
           for b, r in sorted(nominal.items())])

    # Range-fitted coefficients
    sweeps = defaultdict(lambda: defaultdict(list))
    for r in rows:
        b = r["burnup_GWd_tU"]
        if r["tcool_K"] == r["tmod_K"] == NOM_T:
            sweeps[b]["fuel"].append((r["tfuel_K"], r))
        if r["tfuel_K"] == NOM_TF and r["tmod_K"] == NOM_T:
            sweeps[b]["cool"].append((r["tcool_K"], r))
        if r["tfuel_K"] == NOM_TF and r["tcool_K"] == NOM_T:
            sweeps[b]["mod"].append((r["tmod_K"], r))
    coef = []
    for b in sorted(sweeps):
        for name, kind, lo, hi in (("FTC_293_900", "fuel", 293, 900), ("FTC_293_2593", "fuel", 293, 2600),
                                   ("FTC_793_993", "fuel", 793, 994), ("MTC_293_363", "mod", 293, 364),
                                   ("CTC_293_363", "cool", 293, 364)):
            p = [(t, r) for t, r in sweeps[b][kind] if lo <= t <= hi]
            bd = slope([(t, r["k_dragon"], 0.0) for t, r in p])[0]
            bo, so = slope([(t, r["openmc_keff"], r["openmc_keff_sd"]) for t, r in p])
            coef.append({"burnup_GWd_tU": b, "coefficient": name, "points": len(p), "dragon_pcm_K": bd,
                         "openmc_pcm_K": bo, "openmc_sd_pcm_K": so, "diff_pcm_K": bd - bo, "z": (bd - bo) / so})
    write("benchmark_library_coefficients.csv", coef)
    main_coef = [c for c in coef if c["coefficient"] != "FTC_793_993"]
    opposite = [c for c in main_coef if c["dragon_pcm_K"] * c["openmc_pcm_K"] < 0]
    coef_summary = {}
    for name in ("FTC_293_900", "FTC_293_2593", "FTC_793_993", "MTC_293_363", "CTC_293_363"):
        c = [x for x in coef if x["coefficient"] == name]
        coef_summary[name] = {
            "max_abs_diff": max(abs(x["diff_pcm_K"]) for x in c),
            "rms_z": math.sqrt(st.mean(x["z"] ** 2 for x in c)),
            "within_2sigma": sum(abs(x["z"]) <= 2 for x in c),
            "openmc_range": [min(x["openmc_pcm_K"] for x in c), max(x["openmc_pcm_K"] for x in c)],
            "dragon_range": [min(x["dragon_pcm_K"] for x in c), max(x["dragon_pcm_K"] for x in c)],
            "openmc_sd_mean": st.mean(x["openmc_sd_pcm_K"] for x in c)}

    # Linear discrepancy representation (Eq. 2) on the fuel sweep
    fuel = {b: sorted((r for _, r in sweeps[b]["fuel"]), key=lambda r: r["tfuel_K"]) for b in sweeps}

    def params(states):
        nom = next(r for r in states if r["tfuel_K"] == NOM_TF)
        bo = slope([(r["tfuel_K"], r["openmc_keff"], 0.0) for r in states])[0]
        bd = slope([(r["tfuel_K"], r["k_dragon"], 0.0) for r in states])[0]
        return rho(nom["openmc_keff"]) - rho(nom["k_dragon"]), bo - bd

    def assess(states, off, db):
        out = []
        for r in states:
            obs = rho(r["openmc_keff"]) - rho(r["k_dragon"])
            mod = off + db * (r["tfuel_K"] - NOM_TF)
            kp = 1.0 / (1.0 - (rho(r["k_dragon"]) + mod) / 1e5)
            out.append((obs, mod, (kp - r["openmc_keff"]) * 1e5))
        return out

    def metric(name, res, ncal):
        obs = [o for o, _, _ in res]
        mo = st.mean(obs)
        return {"assessment": name, "calibration_states": ncal, "assessment_states": len(res),
                "R2_delta_rho": 1 - sum((m - o) ** 2 for o, m, _ in res) / sum((o - mo) ** 2 for o in obs),
                "RMS_delta_k_residual_pcm": math.sqrt(st.mean(e * e for _, _, e in res)),
                "max_abs_delta_k_residual_pcm": max(abs(e) for _, _, e in res)}

    ins, ht, hb, ncal_t = [], [], [], 0
    for b, states in fuel.items():
        ins += assess(states, *params(states))
        other = [r for r in states if r["tfuel_K"] != NOM_TF]
        train = other[::2] + [r for r in states if r["tfuel_K"] == NOM_TF]
        ncal_t += len(train)
        ht += assess(other[1::2], *params(train))
        if b % 1.0:
            lp, rp = params(fuel[b - 0.5]), params(fuel[b + 0.5])
            hb += assess(states, (lp[0] + rp[0]) / 2, (lp[1] + rp[1]) / 2)
    em = [metric("in_sample", ins, 1312), metric("withheld_temperatures", ht, ncal_t),
          metric("withheld_burnups", hb, 21 * 32)]
    fuel_rows = [r for s in fuel.values() for r in s]
    em_extra = {"fuel_sweep_rms_dk_pcm": math.sqrt(st.mean(r["dk"] ** 2 for r in fuel_rows)),
                "fuel_sweep_temp_dependent_rms_pcm": math.sqrt(st.mean(r["dk_temp"] ** 2 for r in fuel_rows))}
    write("benchmark_library_error_model.csv", em)

    # Two-group constants vs OpenMC whole-cell tallies, all 2132 states
    twog = []
    for m in ("Sa1", "Sa2", "nuSf1", "nuSf2", "nuSf2_over_Sa2"):
        rel = []
        for r in rows:
            t, x = tal[r["case_id"]], lib[r["case_id"]]
            o = (float(t["nuSf2_openmc_rate"]) / float(t["Sa2_openmc_rate"]) if m == "nuSf2_over_Sa2"
                 else float(t[f"{m}_openmc_rate"]))
            d = float(x["nuSf2"]) / float(x["Sa2"]) if m == "nuSf2_over_Sa2" else float(x[m])
            rel.append(100 * (o / d - 1.0))
        a = [abs(e) for e in rel]
        twog.append({"metric": m, "rows": len(rel), "mean_rel_pct": st.mean(rel), "min_rel_pct": min(rel),
                     "max_rel_pct": max(rel), "max_abs_pct": max(a), "inside_1pct": sum(e <= 1 for e in a),
                     "inside_2pct": sum(e <= 2 for e in a), "inside_5pct": sum(e <= 5 for e in a)})
    write("benchmark_library_two_group.csv", twog)

    worst = max(rows, key=lambda r: abs(r["dk"]))
    worst_ge2p5 = max((r for r in rows if r["burnup_GWd_tU"] >= 2.5), key=lambda r: abs(r["dk"]))
    summary = {
        "variant": VARIANT, "families": fam, "thresholds": thr,
        "ci95_vs_%d" % SCREEN: {k: sum(r["ci95"] == k for r in rows) for k in ("inside", "outside", "borderline")},
        "sigma_pcm": {"mean": st.mean(r["sigma"] for r in rows), "max": max(r["sigma"] for r in rows)},
        "worst_state": {k: worst[k] for k in ("case_id", "dk", "sigma")},
        "worst_state_BU_ge_2p5": {k: worst_ge2p5[k] for k in ("case_id", "dk", "sigma")},
        "nominal_dk_range": [min(r["dk"] for r in nominal.values()), max(r["dk"] for r in nominal.values())],
        "nominal_dk_BU_ge_2p5_range": [min(r["dk"] for b, r in nominal.items() if b >= 2.5),
                                       max(r["dk"] for b, r in nominal.items() if b >= 2.5)],
        "nominal_selected": {b: nominal[b]["dk"] for b in (0, 0.5, 1, 2.5, 5, 7.5, 10, 15, 20)},
        "coefficients": coef_summary,
        "coefficient_pairs": len(main_coef),
        "opposite_sign_pairs": [{k: c[k] for k in ("burnup_GWd_tU", "coefficient", "dragon_pcm_K",
                                                   "openmc_pcm_K", "openmc_sd_pcm_K")} for c in opposite],
        "coefficients_selected": [c for c in coef if c["burnup_GWd_tU"] in (0, 5, 7.5, 10, 15, 20)],
        "error_model": em, "error_model_reference": em_extra, "two_group": twog}
    (ROOT / "qa").mkdir(exist_ok=True)
    (ROOT / "qa" / "benchmark_library_analysis.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
