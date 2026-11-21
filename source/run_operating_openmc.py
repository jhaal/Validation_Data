#!/usr/bin/env python3
"""Resumable operating-state OpenMC queue, with a progress record per case.

Usage: python source/run_operating_openmc.py --binary /path/to/openmc --threads 8
Requires h5py. Run under caffeinate on macOS for a long calculation.
Existing completed statepoints are checked and reused; failed cases stop the
queue. No manuscript result is claimed until the manifest is complete.
"""
import argparse
import csv
import hashlib
import json
import os
import subprocess
import time
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

import h5py

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "non_upload/operating_coefficients"


def read_state(row):
    folder = WORK / "openmc_cases" / row["case_id"]
    settings = ET.parse(folder / "settings.xml").getroot()
    batches = int(settings.findtext("batches"))
    sp = folder / f"statepoint.{batches}.h5"
    if not sp.is_file(): return None
    with h5py.File(sp) as f:
        assert int(f["seed"][()]) == int(row["seed"])
        assert int(f["current_batch"][()]) == batches
        k, sd = f["k_combined"][()]
        return {"case_id": row["case_id"], "openmc_keff": float(k), "openmc_keff_sd": float(sd),
                "active_histories": int(f["n_particles"][()]) * (batches-int(f["n_inactive"][()])),
                "seed": int(f["seed"][()]), "statepoint_file": str(sp)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--schedule", default="static", help="OpenMP work schedule; recorded per run")
    parser.add_argument("--exclude", nargs="*", default=[])
    args = parser.parse_args()
    rows = list(csv.DictReader((ROOT / "data/operating_coefficient_states.csv").open(newline="")))
    # Obtain interpretable paired results early, before the full matrix completes.
    def priority(r):
        variable = {"fuel": 0, "moderator": 1, "coolant": 2, "none": 3}[r["variable"]]
        mode = 0 if r["density_mode"] == "pressure" else 1
        return float(r["burnup_GWd_tU"]), variable, mode, -float(r["step_K"]), float(r["sign"])
    rows.sort(key=priority)
    binary_hash = hashlib.sha256(args.binary.read_bytes()).hexdigest()
    (WORK / "openmc_queue_build.json").write_text(json.dumps({
        "binary": str(args.binary), "binary_sha256": binary_hash,
        "version": subprocess.check_output([str(args.binary), "--version"], text=True),
        "threads": args.threads, "schedule": args.schedule,
        "started_utc": datetime.now(timezone.utc).isoformat()}, indent=2)+"\n")
    for row in rows:
        if row["case_id"] in args.exclude: continue
        result = read_state(row)
        if result is None:
            folder = WORK / "openmc_cases" / row["case_id"]
            start = time.monotonic()
            status = {"case_id": row["case_id"], "status": "running", "completed": sum(read_state(r) is not None for r in rows),
                      "total": len(rows), "started_utc": datetime.now(timezone.utc).isoformat()}
            (WORK / "progress.json").write_text(json.dumps(status, indent=2)+"\n")
            print(json.dumps(status), flush=True)
            with (folder / "openmc_stdout.log").open("w") as output:
                env = os.environ.copy()
                env["OMP_SCHEDULE"] = args.schedule
                process = subprocess.run([str(args.binary), "-s", str(args.threads)], cwd=folder,
                                         env=env, stdout=output, stderr=subprocess.STDOUT)
            if process.returncode != 0 or (result := read_state(row)) is None:
                status["status"] = "failed"
                (WORK / "progress.json").write_text(json.dumps(status, indent=2)+"\n")
                raise SystemExit(f"OpenMC failed: {row['case_id']}; see {folder / 'openmc_stdout.log'}")
            result["elapsed_seconds"] = time.monotonic()-start
            (folder / "run_provenance.json").write_text(json.dumps({**result, "binary_sha256": binary_hash,
                "binary": str(args.binary), "threads": args.threads, "schedule": args.schedule}, indent=2)+"\n")
        print(json.dumps(result), flush=True)
        done = [x for r in rows if (x := read_state(r)) is not None]
        # Fixed fields because reused early cross-check runs lack elapsed_seconds.
        keys = [k for k in done[0] if k != "elapsed_seconds"]
        with (ROOT / "data/operating_openmc_states.csv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore"); writer.writeheader(); writer.writerows(done)
    missing = [r["case_id"] for r in rows if read_state(r) is None]
    (WORK / "progress.json").write_text(json.dumps({"status": "complete" if not missing else "partial",
         "completed": len(rows)-len(missing), "total": len(rows), "missing": missing}, indent=2)+"\n")


if __name__ == "__main__": main()
