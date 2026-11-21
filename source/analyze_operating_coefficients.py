#!/usr/bin/env python3
"""Summarize local central differences at the illustrative operating reference.

Use --allow-partial only for progress diagnostics. Submission results require
all 63 DRAGON and all 63 OpenMC states. Independent-state errors are justified
by distinct run seeds; systematic/model-form errors are not included.
"""
import argparse
import csv
import json
import math
import re
import statistics
import xml.etree.ElementTree as ET
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "non_upload/operating_coefficients"
DRAGON = Path("/Users/alokjha/DRAGON/5.1/Dragon/Darwin_arm64")
NUM = r"([-+\d.]+[eE][-+]?\d+)"


def rho(k): return (1. - 1. / k)*1e5


def write(path, rows):
    if not rows: return
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--allow-partial", action="store_true")
    args = p.parse_args()
    states = list(csv.DictReader((ROOT / "data/operating_coefficient_states.csv").open(newline="")))
    dragon, checks, complete_decks = {}, {}, []
    for bu in sorted({float(r["burnup_GWd_tU"]) for r in states}):
        name = f"OP_CONV_J2_BU{bu:04.1f}".replace(".", "p") + ".result"
        file = DRAGON / name
        if file.is_file():
            text = file.read_text(errors="replace")
            if "CONVERGENCE NOT REACHED" in text:
                # Keep diagnostic results separate; an end marker is insufficient.
                continue
            for case, k in re.findall(r"^>\|OP_RESULT\s+(\S+)\s+"+NUM, text, re.M):
                dragon[case] = float(k)
            for case, dc, dm in re.findall(r"^>\|OP_DENSITY\s+(\S+)\s+"+NUM+r"(?:\s+"+NUM+r")?", text, re.M):
                checks[case] = (float(dc), float(dm) if dm else None)
            if ">|VRF_COMPLETED" in text: complete_decks.append(bu)
    omc, convergence, tally_rows = {}, [], []
    for row in states:
        case = row["case_id"]
        if case in checks:
            dc, dm = checks[case]
            assert math.isclose(dc, .064886*float(row["coolant_scale"]), rel_tol=2e-6), (case, dc)
            if dm is not None:
                assert math.isclose(dm, .064886*float(row["moderator_scale"]), rel_tol=2e-6), (case, dm)
        folder = WORK / "openmc_cases" / case
        batches = int(ET.parse(folder / "settings.xml").getroot().findtext("batches"))
        sp = folder / f"statepoint.{batches}.h5"
        if not sp.is_file(): continue
        with h5py.File(sp) as f:
            k, sd = map(float, f["k_combined"][()])
            assert int(f["seed"][()]) == int(row["seed"])
            omc[case] = (k, sd)
            generations = np.asarray(f["k_generation"])
            n_in = int(f["n_inactive"][()])
            active, tail = generations[n_in:], generations[n_in//2:n_in]
            s = active.std(ddof=1)
            a, b = np.array_split(active, 2)
            z_half = (a.mean()-b.mean())/(s*math.sqrt(1/len(a)+1/len(b)))
            z_tail = (tail.mean()-active.mean())/(s*math.sqrt(1/len(tail)+1/len(active)))
            entropy = np.asarray(f["entropy"]) if "entropy" in f else np.array([])
            entropy_sd = entropy[n_in:].std(ddof=1) if len(entropy) else math.nan
            entropy_z = (entropy[n_in//2:n_in].mean()-entropy[n_in:].mean())/(entropy_sd*math.sqrt(1/len(tail)+1/len(active))) if len(entropy) and entropy_sd else math.nan
            convergence.append({"case_id": case, "z_active_halves": z_half, "z_inactive_tail": z_tail,
                                "entropy_recorded": bool(len(entropy)), "entropy_z_inactive_tail": entropy_z})
    output = []
    for row in states:
        case = row["case_id"]
        if case in dragon or case in omc:
            pair = omc.get(case, (math.nan, math.nan))
            output.append({**row, "k_dragon": dragon.get(case, math.nan),
                           "k_openmc": pair[0], "k_openmc_sd": pair[1],
                           "delta_k_pcm": (pair[0]-dragon[case])*1e5 if case in dragon else math.nan,
                           "dragon_coolant_density_checked": case in checks,
                           "dragon_moderator_density_checked": case in checks and checks[case][1] is not None})
    pairs = {}
    for row in states:
        if row["variable"] == "none": continue
        key = (float(row["burnup_GWd_tU"]), row["variable"], row["density_mode"], float(row["step_K"]))
        pairs.setdefault(key, {})[int(float(row["sign"]))] = row["case_id"]
    coefficients = []
    for (bu, var, mode, h), endpoints in sorted(pairs.items()):
        minus, plus = endpoints[-1], endpoints[1]
        drg = (rho(dragon[plus])-rho(dragon[minus]))/(2*h) if plus in dragon and minus in dragon else math.nan
        if plus in omc and minus in omc:
            kp, sp = omc[plus]; km, sm = omc[minus]
            mc = (rho(kp)-rho(km))/(2*h)
            error = math.hypot(sp/kp**2, sm/km**2)*1e5/(2*h)
        else: mc, error = math.nan, math.nan
        coefficients.append({"burnup_GWd_tU": bu, "variable": var, "density_mode": mode, "step_K": h,
                             "dragon_pcm_K": drg, "openmc_pcm_K": mc, "openmc_sd_pcm_K": error,
                             "dragon_minus_openmc_pcm_K": drg-mc,
                             "difference_z_independent": (drg-mc)/error if error>0 else math.nan,
                             "openmc_sign_2sd": "unresolved" if abs(mc)<=2*error else
                             ("positive" if mc>0 else "negative") if math.isfinite(mc) else "pending"})
    step_tests = []
    for bu in sorted({float(r["burnup_GWd_tU"]) for r in states}):
        for var, modes in (("fuel", ("fixed",)), ("coolant", ("fixed", "pressure")), ("moderator", ("fixed", "pressure"))):
            records = sorted([c for c in coefficients if c["burnup_GWd_tU"]==bu and c["variable"]==var and c["density_mode"] in modes], key=lambda c:(c["density_mode"],c["step_K"]))
            for mode in modes:
                small, large = [c for c in records if c["density_mode"] == mode]
                step_tests.append({"burnup_GWd_tU": bu, "variable": var, "density_mode": mode,
                                   "small_step_K": small["step_K"], "large_step_K": large["step_K"],
                                   "dragon_step_change_pcm_K": small["dragon_pcm_K"]-large["dragon_pcm_K"],
                                   "openmc_step_change_pcm_K": small["openmc_pcm_K"]-large["openmc_pcm_K"],
                                   "openmc_step_change_sd_pcm_K": math.hypot(small["openmc_sd_pcm_K"],large["openmc_sd_pcm_K"])})
    density_contributions = []
    for bu in sorted({float(r["burnup_GWd_tU"]) for r in states}):
        for variable in ("coolant", "moderator"):
            for h in (5., 10.):
                fixed = next(c for c in coefficients if c["burnup_GWd_tU"]==bu and c["variable"]==variable
                             and c["step_K"]==h and c["density_mode"]=="fixed")
                path = next(c for c in coefficients if c["burnup_GWd_tU"]==bu and c["variable"]==variable
                            and c["step_K"]==h and c["density_mode"]=="pressure")
                density_contributions.append({"burnup_GWd_tU":bu, "variable":variable, "step_K":h,
                    "dragon_density_contribution_pcm_K":path["dragon_pcm_K"]-fixed["dragon_pcm_K"],
                    "openmc_density_contribution_pcm_K":path["openmc_pcm_K"]-fixed["openmc_pcm_K"],
                    "openmc_sd_pcm_K":math.hypot(path["openmc_sd_pcm_K"], fixed["openmc_sd_pcm_K"]),
                    "interpretation":"Pressure-path central difference minus fixed-density central difference; includes O(h^2) finite-step terms"})
    done = len(dragon)==len(omc)==63 and len(complete_decks)==3
    summary = {"complete": done, "required_states_per_code":63, "dragon_states":len(dragon), "openmc_states":len(omc),
               "complete_dragon_burnups":complete_decks,
               "scope":"Local finite-difference coefficients at an illustrative operating reference; fixed shared inventories, structure and geometry; uncorrected resonance scattering",
               "fuel_sign_resolution":"Reported separately at each step; unresolved signs must not be assigned a crossover",
               "coefficients":coefficients,"step_size_comparisons":step_tests}
    # Partial files have explicit names and are never mistaken for final tables.
    prefix = "operating_partial" if not done else "operating"
    write(ROOT / f"data/{prefix}_states.csv", output)
    write(ROOT / f"data/{prefix}_coefficients.csv", coefficients)
    write(ROOT / f"data/{prefix}_step_convergence.csv", step_tests)
    write(ROOT / f"data/{prefix}_density_contributions.csv", density_contributions)
    write(ROOT / f"data/{prefix}_source_convergence.csv", convergence)
    (ROOT / f"qa/{prefix}_analysis.json").write_text(json.dumps(summary, indent=2, allow_nan=not done)+"\n")
    print(json.dumps({k:v for k,v in summary.items() if k not in ("coefficients","step_size_comparisons")}, indent=2))
    if not done and not args.allow_partial: raise SystemExit("Operating matrix is not complete; partial diagnostic files only")


if __name__ == "__main__": main()
