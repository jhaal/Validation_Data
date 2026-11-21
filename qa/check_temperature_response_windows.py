#!/usr/bin/env python3
"""Verify window results against existing fits and independent regression."""
import csv
import json
import math
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "source"))
from analyze_temperature_response_windows import load, rho


def main():
    rows = list(csv.DictReader((ROOT / "data/temperature_response_windows.csv").open(newline="")))
    old = {(float(r["burnup_GWd_tU"]), r["coefficient"]): r for r in csv.DictReader(
        (ROOT / "data/benchmark_library_coefficients.csv").open(newline=""))}
    sweeps = load()
    assert len(rows) == len({(r["burnup_GWd_tU"], r["window"]) for r in rows}) == 164
    mapping = {"cold": "FTC_293_900", "full": "FTC_293_2593", "near_nominal": "FTC_793_993"}
    for row in rows:
        b = float(row["burnup_GWd_tU"])
        states = [r for r in sweeps[b] if float(row["Tmin_K"]) <= r["tfuel_K"] <= float(row["Tmax_K"])]
        t = [r["tfuel_K"] for r in states]
        assert len(states) == int(row["points"])
        for code, key in (("openmc", "openmc_keff"), ("dragon", "k_dragon")):
            result = statistics.linear_regression(t, [rho(r[key]) for r in states]).slope
            assert math.isclose(result, float(row[f"{code}_slope_pcm_K"]), abs_tol=1e-12)
        assert float(row["openmc_sd_upper_bound_pcm_K"]) >= float(row["openmc_sd_independent_pcm_K"])
        if row["window"] in mapping:
            previous = old[b, mapping[row["window"]]]
            for newkey, oldkey in (("openmc_slope_pcm_K", "openmc_pcm_K"),
                                    ("dragon_slope_pcm_K", "dragon_pcm_K"),
                                    ("openmc_sd_independent_pcm_K", "openmc_sd_pcm_K")):
                assert math.isclose(float(row[newkey]), float(previous[oldkey]), abs_tol=1e-12)
        if row["window"] == "within_fuel_tabulation":
            assert len(states) == 26 and int(row["fuel_points_above_2000K"]) == 0
    selected = {r["window"]: r for r in rows if float(r["burnup_GWd_tU"]) == 10.}
    assert selected["cold"]["openmc_sign_at_2sd_independent"] == "negative"
    assert selected["full"]["openmc_sign_at_2sd_independent"] == "positive"
    assert selected["near_nominal"]["openmc_sign_at_2sd_independent"] == "unresolved"
    summary = json.loads((ROOT / "qa/temperature_response_windows.json").read_text())
    assert round(summary["windows"]["full"]["openmc_point_estimate_crossovers_GWd_tU"][0], 2) == 5.04
    assert round(summary["windows"]["within_fuel_tabulation"]["openmc_point_estimate_crossovers_GWd_tU"][0], 2) == 6.67
    assert round(summary["windows"]["cold"]["openmc_point_estimate_crossovers_GWd_tU"][0], 2) == 14.03
    assert len(summary["windows"]["near_nominal"]["openmc_point_estimate_crossovers_GWd_tU"]) == 5
    print("PASS all 164 window fits, independent regression, agreement with existing results, tabulation exclusion and crossover example")


if __name__ == "__main__":
    main()
