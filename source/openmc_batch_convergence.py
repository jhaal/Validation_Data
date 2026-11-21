#!/usr/bin/env python3
"""Source-convergence diagnostics from the OpenMC matrix statepoints (batch k_generation).

For each case with statepoint.200.h5 (200 batches, 50 inactive):
  z_inactive_tail : (mean k of inactive batches 26-50 - mean k of active batches) / SE
  z_active_halves : (mean k of active batches 51-125 - mean of 126-200) / SE
  z_active_trend  : least-squares slope of k over active batches / its standard error
SE uses the sample variance of active batch k. Under a converged, stationary source each
z is approximately standard normal.
Usage: python openmc_batch_convergence.py <cases_dir> <out.csv>   (requires h5py, numpy)
"""
import csv
import pathlib
import sys

import h5py
import numpy as np


def main() -> None:
    cases, out = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
    rows = []
    for sp in sorted(cases.glob("*/statepoint.200.h5")):
        with h5py.File(sp, "r") as f:
            k = np.array(f["k_generation"])
            n_in = int(f["n_inactive"][()])
        act = k[n_in:]
        s = act.std(ddof=1)
        tail = k[n_in // 2:n_in]
        z_tail = (tail.mean() - act.mean()) / (s * np.sqrt(1 / len(tail) + 1 / len(act)))
        h1, h2 = act[: len(act) // 2], act[len(act) // 2:]
        z_half = (h1.mean() - h2.mean()) / (s * np.sqrt(1 / len(h1) + 1 / len(h2)))
        x = np.arange(len(act)) - (len(act) - 1) / 2
        b = (x * (act - act.mean())).sum() / (x * x).sum()
        resid = act - act.mean() - b * x
        se_b = np.sqrt((resid ** 2).sum() / (len(act) - 2) / (x * x).sum())
        rows.append({"case_id": sp.parent.name, "n_inactive": n_in, "batch_k_sd_pcm": 1e5 * s,
                     "z_inactive_tail": z_tail, "z_active_halves": z_half, "z_active_trend": b / se_b})
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    for key in ("z_inactive_tail", "z_active_halves", "z_active_trend"):
        z = np.array([r[key] for r in rows])
        print(f"{key:16s} n={len(z)} mean={z.mean():+.3f} sd={z.std(ddof=1):.3f} "
              f"|z|>2: {np.mean(abs(z) > 2):.3f} (normal 0.046)  |z|>3: {np.mean(abs(z) > 3):.4f} (normal 0.0027)  max|z|={abs(z).max():.2f}")


if __name__ == "__main__":
    main()
