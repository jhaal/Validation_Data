#!/usr/bin/env python3
"""Check the tally-moment correction without external statepoints or OpenMC."""
import csv
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "source"))
from extract_openmc_tally_statistics import GROUPS, SCORES, mean_sd, ratio_bound


def close(a, b):
    assert math.isclose(float(a), float(b), rel_tol=1e-12, abs_tol=1e-15), (a, b)


def main():
    # Independent known batches [1, 2, 3] distinguish SD from standard error.
    m, s = mean_sd(6., 14., 3)
    close(m, 2.)
    close(s, math.sqrt(1. / 3.))
    close(mean_sd(6., 12., 3)[1], 0.)
    for args in ((6., 10., 3), (6., 14., 1)):
        try:
            mean_sd(*args)
        except ValueError:
            pass
        else:
            raise AssertionError("Invalid moments/count accepted")
    value, bound = ratio_bound(2., .2, 4., .8)
    close(value, .5)
    close(bound, .15)  # negative perfect correlation gives the largest linearized SD
    raw = {r["case_id"]: r for r in csv.DictReader(
        (ROOT / "data/openmc_two_group_tally_moments.csv").open(newline=""))}
    rows = list(csv.DictReader((ROOT / "data/openmc_two_group_rate_tallies.csv").open(newline="")))
    assert len(raw) == len(rows) == len({r["case_id"] for r in rows}) == 2132
    for row in rows:
        original = raw[row["case_id"]]
        n = int(original["n_realizations"])
        assert n == int(row["n_realizations"]) == 150
        for group in GROUPS:
            for score in SCORES:
                key = f"{score}_{group}"
                m, s = mean_sd(float(original[key + "_sum"]), float(original[key + "_sum_sq"]), n)
                close(row[key + "_mean"], m)
                close(row[key + "_sd"], s)
            for score, name in (("absorption", "Sa"), ("fission", "fission"), ("nu_fission", "nuSf")):
                suffix = "1" if group == "g1_fast" else "2"
                key = f"{name}{suffix}_openmc_rate"
                v, b = ratio_bound(float(row[f"{score}_{group}_mean"]), float(row[f"{score}_{group}_sd"]),
                                   float(row[f"flux_{group}_mean"]), float(row[f"flux_{group}_sd"]))
                close(row[key], v)
                close(row[key + "_sd_upper_bound_delta"], b)
                # Legacy ratios used raw sums: normalize both or neither.
                close(v, float(original[f"{score}_{group}_sum"]) / float(original[f"flux_{group}_sum"]))
        for key, numerator, denominator in (
            ("thermal_to_fast_flux_ratio_from_tally", "flux_g2_thermal", "flux_g1_fast"),
            ("nuSf2_over_Sa2_openmc_rate", "nu_fission_g2_thermal", "absorption_g2_thermal")):
            v, b = ratio_bound(float(row[numerator + "_mean"]), float(row[numerator + "_sd"]),
                               float(row[denominator + "_mean"]), float(row[denominator + "_sd"]))
            close(row[key], v)
            close(row[key + "_sd_upper_bound_delta"], b)
    # Independently printed first-case tally values, six significant figures.
    first = next(r for r in rows if r["case_id"] == "BU00p0_TFUEL0293p15_TCOOL0293p15_TMOD0293p15")
    assert math.isclose(float(first["flux_g2_thermal_mean"]), 232.609, rel_tol=6e-6)
    assert math.isclose(float(first["flux_g2_thermal_sd"]), .102399, rel_tol=6e-6)
    summary = json.loads((ROOT / "qa/openmc_tally_statistics_summary.json").read_text())
    for key, reported in summary["max_relative_sd_upper_bound_pct"].items():
        close(reported, max(100 * float(r[key + "_sd_upper_bound_delta"]) / float(r[key]) for r in rows))
    print("PASS known-batch statistics, invalid-input rejection, all 2132 means/errors/ratios/bounds, printed reference and summary")


if __name__ == "__main__":
    main()
