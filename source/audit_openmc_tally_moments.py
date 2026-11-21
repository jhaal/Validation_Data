#!/usr/bin/env python3
"""Compare saved moments with locally available statepoints and printed tallies.

Usage: python source/audit_openmc_tally_moments.py
Requires h5py/numpy. Missing external statepoints are recorded, never counted as
verified. Counts for such cases are checked against the per-case run-settings
audit (one realization per active batch in these runs). No transport is rerun.
"""
import csv
import hashlib
import json
import math
import re
from collections import Counter
from pathlib import Path

import h5py
import numpy as np

from extract_openmc_tally_statistics import GROUPS, SCORES, ROOT, mean_sd


def decode(x):
    return x.decode() if isinstance(x, bytes) else str(x)


def main():
    moment_file = ROOT / "data/openmc_two_group_tally_moments.csv"
    moments = list(csv.DictReader(moment_file.open(newline="")))
    settings = {r["case_id"]: r for r in csv.DictReader(
        (ROOT / "data/openmc_case_run_settings.csv").open(newline=""))}
    checked, missing, printed, seeds = [], [], 0, Counter()
    for row in moments:
        case = row["case_id"]
        n = int(row["n_realizations"])
        setting = settings[case]
        assert n == int(setting["batches"]) - int(setting["inactive"]), case
        sp = Path(row["statepoint_file"])
        if not sp.is_file():
            missing.append(case)
            continue
        with h5py.File(sp, "r") as f:
            tallies = [g for k, g in f["tallies"].items() if k.startswith("tally ")
                       and decode(g["name"][()]) == "whole_cell_two_group_rates"]
            assert len(tallies) == 1, case
            g = tallies[0]
            assert n == int(g["n_realizations"][()]), case
            assert int(f["generations_per_batch"][()]) == 1, case
            assert int(f["n_batches"][()]) == int(setting["batches"]), case
            assert int(f["n_inactive"][()]) == int(setting["inactive"]), case
            scores = [decode(x).replace("-", "_") for x in g["score_bins"][()]]
            assert set(scores) == set(SCORES), case
            filters = g["filters"][()]
            assert len(filters) == 1, case
            bins = f[f"tallies/filters/filter {int(filters[0])}/bins"][()]
            assert np.array_equal(bins, [1e-5, 4., 2e7]), case
            raw = g["results"][()]
            assert raw.shape == (2, 4, 2), case
            seeds[str(int(f["seed"][()]))] += 1
            for i, group in enumerate(GROUPS):
                for j, score in enumerate(scores):
                    key = f"{score}_{group}"
                    assert float(row[key + "_sum"]) == float(raw[i, j, 0]), (case, key)
                    assert float(row[key + "_sum_sq"]) == float(raw[i, j, 1]), (case, key)
        text_path = sp.parent / "tallies.out"
        if text_path.is_file():
            pairs = re.findall(r"(?:Flux|Absorption Rate|Fission Rate|Nu-Fission Rate)\s+"
                               r"([-+\d.eE]+)\s+\+/-\s+([-+\d.eE]+)", text_path.read_text())
            assert len(pairs) == 8, case
            for (group, score), (m, s) in zip(
                    [(g, s) for g in GROUPS for s in SCORES], pairs):
                key = f"{score}_{group}"
                mean, sd = mean_sd(float(row[key + "_sum"]), float(row[key + "_sum_sq"]), n)
                # tallies.out prints six significant figures.
                assert math.isclose(mean, float(m), rel_tol=6e-6), (case, key, "mean")
                assert math.isclose(sd, float(s), rel_tol=6e-6, abs_tol=1e-12), (case, key, "sd")
            printed += 1
        checked.append(case)
    summary = {
        "cases": len(moments), "statepoints_verified": len(checked),
        "printed_tallies_verified": printed,
        "statepoints_unavailable": len(missing),
        "unavailable_case_ids": missing,
        "realization_count_basis_for_unavailable": "Per-case audited batches minus inactive; raw moment origin retained",
        "raw_moments_sha256": hashlib.sha256(moment_file.read_bytes()).hexdigest(),
        "seeds_in_available_statepoints": dict(seeds),
        "note": "Seed equality alone does not quantify cross-state covariance.",
    }
    (ROOT / "qa/openmc_tally_moments_audit.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "unavailable_case_ids"}, indent=2))


if __name__ == "__main__":
    main()
