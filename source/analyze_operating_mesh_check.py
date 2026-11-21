#!/usr/bin/env python3
"""Compare 20/40-ring operating eigenvalues and larger-step coefficients."""
import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LISTINGS = Path("/Users/alokjha/DRAGON/5.1/Dragon/Darwin_arm64")
NUM = r"([-+\d.]+[eE][-+]?\d+)"


def load(prefix):
    results = {}
    complete = []
    for bu in (0., 7.5, 20.):
        name = f"{prefix}_BU{bu:04.1f}".replace(".", "p")+".result"
        p = LISTINGS/name
        if not p.is_file(): continue
        text = p.read_text(errors="replace")
        if "CONVERGENCE NOT REACHED" in text: continue
        for case, value in re.findall(r"^>\|OP_RESULT\s+(\S+)\s+"+NUM, text, re.M):
            results[case] = float(value)
        if ">|VRF_COMPLETED" in text: complete.append(bu)
    return results, complete


def write(name, rows):
    with (ROOT / "data" / name).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def main():
    states = list(csv.DictReader((ROOT / "data/operating_mesh_check_states.csv").open(newline="")))
    k20, done20 = load("OP_CONV_J2")
    k40, done40 = load("OP_M40_J2")
    expected = {r["case_id"] for r in states}
    if not expected <= k20.keys() or not expected <= k40.keys() or len(done40)!=3:
        raise SystemExit(f"Operating mesh check incomplete: {len(k40)}/33 fine-mesh states")
    rows = [{"case_id":r["case_id"], "burnup_GWd_tU":r["burnup_GWd_tU"],
             "k_20_rings":k20[r["case_id"]], "k_40_rings":k40[r["case_id"]],
             "change_40_minus_20_pcm":1e5*(k40[r["case_id"]]-k20[r["case_id"]])} for r in states]
    pairs = {}
    for r in states:
        if r["variable"]=="none": continue
        key = (float(r["burnup_GWd_tU"]),r["variable"],r["density_mode"],float(r["step_K"]))
        pairs.setdefault(key,{})[int(float(r["sign"]))] = r["case_id"]
    coefficients = []
    for (bu,var,mode,h), points in sorted(pairs.items()):
        minus, plus = points[-1], points[1]
        alpha20 = 1e5*(1/k20[minus]-1/k20[plus])/(2*h)
        alpha40 = 1e5*(1/k40[minus]-1/k40[plus])/(2*h)
        coefficients.append({"burnup_GWd_tU":bu,"variable":var,"density_mode":mode,"step_K":h,
                             "alpha_20_rings_pcm_K":alpha20,"alpha_40_rings_pcm_K":alpha40,
                             "change_40_minus_20_pcm_K":alpha40-alpha20})
    assert len(rows)==33 and len(coefficients)==15
    write("operating_mesh_states.csv",rows)
    write("operating_mesh_coefficients.csv",coefficients)
    summary = {"complete":True,"states_compared":33,"coefficients_compared":15,
               "moderator_rings":[20,40],"tracking":"8 angles, 10 lines/cm in both meshes",
               "max_abs_eigenvalue_change_pcm":max(abs(r["change_40_minus_20_pcm"]) for r in rows),
               "max_abs_coefficient_change_pcm_K":max(abs(r["change_40_minus_20_pcm_K"]) for r in coefficients),
               "scope":"Three selected operating reference states and larger central-difference endpoints; not full-domain convergence or an application acceptance criterion"}
    (ROOT / "qa/operating_mesh_analysis.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))


if __name__ == "__main__": main()
