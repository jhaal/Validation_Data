#!/usr/bin/env python3
"""Independently check completed central differences and their uncertainties."""
import csv
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def close(a,b):
    assert math.isclose(float(a),float(b),rel_tol=1e-10,abs_tol=1e-10),(a,b)


def main():
    summary = json.loads((ROOT / "qa/operating_analysis.json").read_text())
    assert summary["complete"] and summary["dragon_states"]==summary["openmc_states"]==63
    states = {r["case_id"]:r for r in csv.DictReader((ROOT/"data/operating_states.csv").open(newline=""))}
    coefficients = list(csv.DictReader((ROOT/"data/operating_coefficients.csv").open(newline="")))
    assert len(states)==63 and len(coefficients)==30
    assert len({r["seed"] for r in states.values()})==63
    for r in states.values():
        assert r["dragon_coolant_density_checked"]==r["dragon_moderator_density_checked"]=="True"
        assert float(r["k_dragon"])>0 and float(r["k_openmc"])>0 and float(r["k_openmc_sd"])>0
    for row in coefficients:
        matches = [r for r in states.values() if all(r[k]==row[k] for k in
                   ("burnup_GWd_tU","variable","density_mode","step_K"))]
        assert len(matches)==2
        minus = next(r for r in matches if float(r["sign"])<0)
        plus = next(r for r in matches if float(r["sign"])>0)
        h=float(row["step_K"])
        for code,key in (("dragon","k_dragon"),("openmc","k_openmc")):
            km,kp=float(minus[key]),float(plus[key])
            close(row[code+"_pcm_K"],1e5*(1/km-1/kp)/(2*h))
        km,kp=float(minus["k_openmc"]),float(plus["k_openmc"])
        sm,sp=float(minus["k_openmc_sd"]),float(plus["k_openmc_sd"])
        close(row["openmc_sd_pcm_K"],1e5*math.sqrt(sm*sm/km**4+sp*sp/kp**4)/(2*h))
    mesh=json.loads((ROOT/"qa/operating_mesh_analysis.json").read_text())
    assert mesh["complete"] and mesh["states_compared"]==33 and mesh["coefficients_compared"]==15
    solver=json.loads((ROOT/"qa/operating_solver_convergence.json").read_text())
    assert solver["complete"] and solver["states"]==3
    for row in json.loads((ROOT/"qa/operating_input_hashes.json").read_text()):
        assert hashlib.sha256((ROOT/row["path"]).read_bytes()).hexdigest()==row["sha256"],row["path"]
    convergence=list(csv.DictReader((ROOT/"data/operating_source_convergence.csv").open(newline="")))
    assert len(convergence)==63
    flags=[{"case_id":r["case_id"],"metric":key,"z":float(r[key])} for r in convergence
           for key in ("z_active_halves","z_inactive_tail","entropy_z_inactive_tail") if abs(float(r[key]))>3]
    (ROOT/"qa/operating_source_diagnostic_flags.json").write_text(json.dumps({
        "cases":63,"three_sigma_flags":flags,
        "interpretation":"Flags require review; these are pointwise diagnostics with multiple comparisons, not automatic proof of nonconvergence"},indent=2)+"\n")
    print(f"PASS 63 states, 30 independently checked central differences/errors, density records, 33 mesh checks, solver pilot and input hashes; source flags={len(flags)}")


if __name__=="__main__":main()
