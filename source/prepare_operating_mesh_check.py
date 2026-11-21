#!/usr/bin/env python3
"""Prepare 40-ring checks of the nominal and larger-step operating states."""
import csv
from prepare_operating_coefficient_cases import ROOT, install_dragon, write_csv


def main():
    rows = list(csv.DictReader((ROOT / "data/operating_coefficient_states.csv").open(newline="")))
    for row in rows:
        for key in ("burnup_GWd_tU", "tfuel_K", "tcool_K", "tmod_K", "coolant_scale", "moderator_scale", "step_K"):
            row[key] = float(row[key])
    selected = [r for r in rows if r["variable"]=="none" or
                r["step_K"]==(100. if r["variable"]=="fuel" else 10.)]
    assert len(selected)==33
    write_csv(ROOT / "data/operating_mesh_check_states.csv", selected)
    install_dragon(selected, rings=40, prefix="OP_M40_J2")


if __name__ == "__main__": main()
