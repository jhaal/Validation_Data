#!/usr/bin/env python3
"""Rebuild tally statistics and ratio uncertainty bounds from saved raw moments.

Usage: python source/extract_openmc_tally_statistics.py
Requires only Python's standard library. The input moments are OpenMC sums and
sums of squares over n_realizations, NOT means and standard deviations.
See https://docs.openmc.org/en/stable/io_formats/statepoint.html and OpenMC
v0.15.3 openmc/tallies.py (Tally.mean and Tally.std_dev).

Cross-score products are absent from ordinary statepoints. Consequently ratio
uncertainties are reported as first-order upper bounds for arbitrary covariance,
not estimates of the actual standard error and not confidence intervals.
"""
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GROUPS = ("g2_thermal", "g1_fast")
SCORES = ("flux", "absorption", "fission", "nu_fission")


def mean_sd(total, sum_sq, n):
    """OpenMC standard error of the mean, for n active batch realizations."""
    if n < 2 or not all(math.isfinite(x) for x in (total, sum_sq)):
        raise ValueError("Invalid tally moments or realization count")
    mean = total / n
    centered = sum_sq / n - mean * mean
    if centered < -1e-12 * max(abs(sum_sq / n), mean * mean, 1e-300):
        raise ValueError("Negative tally variance: inconsistent raw moments")
    return mean, math.sqrt(max(0.0, centered) / (n - 1))


def ratio_bound(a, sa, b, sb):
    """Delta-method SD upper bound using |Cov(a,b)| <= sa*sb."""
    if a <= 0 or b <= 0 or sa < 0 or sb < 0:
        raise ValueError("Expected positive tally means and nonnegative errors")
    ratio = a / b
    return ratio, ratio * (sa / a + sb / b)


def main():
    raw_path = ROOT / "data/openmc_two_group_tally_moments.csv"
    source = list(csv.DictReader(raw_path.open(newline="")))
    if len(source) != 2132 or len({r["case_id"] for r in source}) != 2132:
        raise ValueError("Expected exactly 2132 unique cases")
    rows = []
    meta = ("global_index", "case_id", "burnup_GWd_tU", "tfuel_K", "tcool_K",
            "tmod_K", "case_dir", "statepoint_file", "tally_status", "message",
            "n_realizations", "n_realizations_basis")
    for raw in source:
        r = {k: raw[k] for k in meta}
        n = int(raw["n_realizations"])
        for group in GROUPS:
            for score in SCORES:
                key = f"{score}_{group}"
                r[f"{key}_mean"], r[f"{key}_sd"] = mean_sd(
                    float(raw[f"{key}_sum"]), float(raw[f"{key}_sum_sq"]), n)
        for group, suffix in (("g1_fast", "1"), ("g2_thermal", "2")):
            for score, name in (("absorption", "Sa"), ("fission", "fission"),
                                ("nu_fission", "nuSf")):
                key = f"{name}{suffix}_openmc_rate"
                r[key], r[f"{key}_sd_upper_bound_delta"] = ratio_bound(
                    r[f"{score}_{group}_mean"], r[f"{score}_{group}_sd"],
                    r[f"flux_{group}_mean"], r[f"flux_{group}_sd"])
        key = "thermal_to_fast_flux_ratio_from_tally"
        r[key], r[f"{key}_sd_upper_bound_delta"] = ratio_bound(
            r["flux_g2_thermal_mean"], r["flux_g2_thermal_sd"],
            r["flux_g1_fast_mean"], r["flux_g1_fast_sd"])
        # Direct production/absorption ratio: the shared flux cancels algebraically.
        key = "nuSf2_over_Sa2_openmc_rate"
        r[key], r[f"{key}_sd_upper_bound_delta"] = ratio_bound(
            r["nu_fission_g2_thermal_mean"], r["nu_fission_g2_thermal_sd"],
            r["absorption_g2_thermal_mean"], r["absorption_g2_thermal_sd"])
        rows.append(r)
    out = ROOT / "data/openmc_two_group_rate_tallies.csv"
    with out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    keys = [k for k in rows[0] if k.endswith("_sd_upper_bound_delta")]
    summary = {
        "cases": len(rows), "n_realizations": sorted({int(r["n_realizations"]) for r in rows}),
        "method": "First-order ratio standard-error upper bound, arbitrary covariance; not a confidence interval",
        "covariance_available": False,
        "max_relative_sd_upper_bound_pct": {
            k.removesuffix("_sd_upper_bound_delta"): max(
                100 * r[k] / r[k.removesuffix("_sd_upper_bound_delta")] for r in rows)
            for k in keys},
    }
    (ROOT / "qa/openmc_tally_statistics_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
