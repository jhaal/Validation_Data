#!/usr/bin/env python3
"""Run the Phase F full-matrix OpenMC benchmark with resumable results."""

from __future__ import annotations

import argparse
import csv
import math
import os
import pathlib
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone

import phase_f_openmc_candu_benchmark as bench


DEFAULT_OUT_DIR = bench.OUT_DIR / "full_matrix_openmc"
RESULTS_CSV = "phaseF_full_matrix_openmc_results.csv"
REPORT_MD = "PHASE_F_OPENMC_FULL_MATRIX_REPORT.md"
MAC_1399_SPLIT_DIR = pathlib.Path(__file__).resolve().parent / "phaseF_split_1399"
MAC_1399_SENTINEL = MAC_1399_SPLIT_DIR / "ENABLE_MAC_0001_1399_DEFAULT"
MAC_1399_CASE_LIST = MAC_1399_SPLIT_DIR / "mac_case_list_0001_1399.csv"


FIELDNAMES = [
    "case_id",
    "burnup_GWd_tU",
    "tfuel_K",
    "tcool_K",
    "tmod_K",
    "dragon_k_eff",
    "openmc_keff",
    "openmc_keff_sd",
    "delta_pcm",
    "thermal_to_fast_flux_ratio",
    "run_status",
    "acceptance_status",
    "case_dir",
    "message",
    "completed_utc",
]


def read_existing(path: pathlib.Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    with path.open(newline="", encoding="utf-8") as handle:
        return {row["case_id"]: row for row in csv.DictReader(handle)}


def write_results(path: pathlib.Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(1, 11):
        try:
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
                writer.writeheader()
                writer.writerows(rows)
            return
        except PermissionError:
            if attempt == 10:
                raise
            time.sleep(2.0)


def state_to_row(
    state: bench.BranchState,
    case_dir: pathlib.Path,
    run_status: str,
    acceptance_status: str,
    metrics: dict[str, float],
    message: str,
) -> dict[str, object]:
    dragon_k = state.k_eff if state.k_eff is not None else float("nan")
    openmc_keff = metrics.get("openmc_keff", float("nan"))
    delta_pcm = (openmc_keff - dragon_k) * 1.0e5 if math.isfinite(openmc_keff) else float("nan")
    return {
        "case_id": state.case_id,
        "burnup_GWd_tU": state.burnup,
        "tfuel_K": state.tfuel,
        "tcool_K": state.tcool,
        "tmod_K": state.tmod,
        "dragon_k_eff": dragon_k,
        "openmc_keff": openmc_keff,
        "openmc_keff_sd": metrics.get("openmc_keff_sd", float("nan")),
        "delta_pcm": delta_pcm,
        "thermal_to_fast_flux_ratio": metrics.get("thermal_to_fast_flux_ratio", float("nan")),
        "run_status": run_status,
        "acceptance_status": acceptance_status,
        "case_dir": str(case_dir),
        "message": message,
        "completed_utc": datetime.now(timezone.utc).isoformat(),
    }


def run_one(payload: tuple[dict[str, object], str, int, int, int, int, float, float]) -> dict[str, object]:
    state_data, case_root_text, particles, batches, inactive, threads, max_sigma, max_delta_pcm = payload
    state = bench.BranchState(
        case_id=str(state_data["case_id"]),
        burnup=float(state_data["burnup"]),
        tfuel=float(state_data["tfuel"]),
        tcool=float(state_data["tcool"]),
        tmod=float(state_data["tmod"]),
        k_eff=float(state_data["k_eff"]) if state_data["k_eff"] is not None else None,
        run_status=str(state_data["run_status"]),
    )
    case_dir = pathlib.Path(case_root_text) / state.case_id
    try:
        inventory_densities = bench.read_inventory_densities(bench.INVENTORY_CSV)
        bench.export_case(case_dir, state, particles, batches, inactive, inventory_densities)
        proc = bench.run_openmc(case_dir, threads)
        (case_dir / "openmc_stdout.log").write_text(proc.stdout, encoding="utf-8")
        if proc.returncode != 0:
            return state_to_row(state, case_dir, "FAIL", "FAIL", {}, proc.stdout[-1000:])
        metrics = bench.parse_statepoint(case_dir)
        sigma = metrics.get("openmc_keff_sd", float("inf"))
        delta_pcm = (metrics["openmc_keff"] - (state.k_eff or float("nan"))) * 1.0e5
        accepted = sigma <= max_sigma and abs(delta_pcm) <= max_delta_pcm
        message = f"sigma={sigma:.6g}; abs_delta_pcm={abs(delta_pcm):.3f}"
        return state_to_row(state, case_dir, "PASS", "PASS" if accepted else "FAIL", metrics, message)
    except Exception as exc:  # noqa: BLE001 - exact failure text is needed for audit.
        return state_to_row(state, case_dir, "FAIL", "FAIL", {}, repr(exc))


def write_report(
    out_dir: pathlib.Path,
    results: list[dict[str, object]],
    total_expected: int,
    particles: int,
    batches: int,
    inactive: int,
    max_sigma: float,
    max_delta_pcm: float,
) -> None:
    completed = len(results)
    run_pass = sum(1 for row in results if row["run_status"] == "PASS")
    accept_pass = sum(1 for row in results if row["acceptance_status"] == "PASS")
    run_fail = completed - run_pass
    accept_fail = completed - accept_pass
    remaining = total_expected - completed
    if remaining == 0 and accept_fail == 0:
        disposition = "PHASE_F_FULL_MATRIX_ACCEPTANCE_PASS"
    elif remaining == 0:
        disposition = "PHASE_F_FULL_MATRIX_COMPLETE_WITH_ACCEPTANCE_FAILURES"
    else:
        disposition = "PHASE_F_FULL_MATRIX_IN_PROGRESS"

    worst = sorted(
        (
            row for row in results
            if str(row.get("delta_pcm", "nan")).lower() != "nan"
        ),
        key=lambda row: abs(float(row["delta_pcm"])),
        reverse=True,
    )[:20]

    lines = [
        "# Phase F OpenMC Full-Matrix Benchmark Report",
        "",
        f"Generated UTC: `{datetime.now(timezone.utc).isoformat()}`",
        f"Disposition: `{disposition}`",
        "",
        "## Scope",
        "",
        "This report tracks full-matrix independent OpenMC execution for all",
        "completed DRAGON burnup-temperature branch states with complete CE",
        "inventory coverage. It is resumable and is updated after each batch.",
        "",
        "## Acceptance Criteria",
        "",
        f"- OpenMC run must complete with return code zero.",
        f"- OpenMC k-effective one-sigma uncertainty must be <= `{max_sigma:.6g}`.",
        f"- Absolute OpenMC minus DRAGON k-effective difference must be <= `{max_delta_pcm:.1f}` pcm.",
        "",
        "## Run Settings",
        "",
        f"- particles: `{particles}`",
        f"- batches: `{batches}`",
        f"- inactive: `{inactive}`",
        f"- Candidate branch states: `{total_expected}`",
        "",
        "## Status",
        "",
        f"- Completed OpenMC cases: `{completed}`",
        f"- Remaining OpenMC cases: `{remaining}`",
        f"- Run PASS: `{run_pass}`",
        f"- Run FAIL: `{run_fail}`",
        f"- Acceptance PASS: `{accept_pass}`",
        f"- Acceptance FAIL: `{accept_fail}`",
        "",
        "## Largest Absolute Delta-pcm Cases",
        "",
        "| case_id | burnup | Tfuel | Tcool | Tmod | DRAGON k | OpenMC k | sigma | delta pcm | acceptance |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in worst:
        lines.append(
            "| {case_id} | {burnup_GWd_tU:.1f} | {tfuel_K:.2f} | {tcool_K:.2f} | {tmod_K:.2f} | "
            "{dragon_k_eff:.6f} | {openmc_keff:.6f} | {openmc_keff_sd:.6f} | {delta_pcm:.1f} | {acceptance_status} |".format(
                **{
                    key: (float(value) if key not in {"case_id", "acceptance_status"} else value)
                    for key, value in row.items()
                    if key in {
                        "case_id", "burnup_GWd_tU", "tfuel_K", "tcool_K", "tmod_K",
                        "dragon_k_eff", "openmc_keff", "openmc_keff_sd", "delta_pcm",
                        "acceptance_status",
                    }
                }
            )
        )
    lines += [
        "",
        "## Artifacts",
        "",
        f"- Results CSV: `{out_dir / RESULTS_CSV}`",
        f"- Case directories: `{out_dir / 'openmc_cases'}`",
        "",
    ]
    report_path = out_dir / REPORT_MD
    for attempt in range(1, 11):
        try:
            report_path.write_text("\n".join(lines), encoding="utf-8")
            return
        except PermissionError:
            if attempt == 10:
                raise
            time.sleep(2.0)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=pathlib.Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--particles", type=int, default=20000)
    parser.add_argument("--batches", type=int, default=200)
    parser.add_argument("--inactive", type=int, default=50)
    parser.add_argument("--threads-per-case", type=int, default=4)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--max-sigma", type=float, default=5.0e-4)
    parser.add_argument("--max-delta-pcm", type=float, default=1000.0)
    parser.add_argument("--max-cases", type=int, default=0)
    parser.add_argument(
        "--case-list",
        type=pathlib.Path,
        default=None,
        help="Optional CSV/text file containing case_id values to run. Used for disjoint distributed workers.",
    )
    parser.add_argument("--rerun-failed", action="store_true")
    parser.add_argument("--rerun-all", action="store_true")
    args = parser.parse_args()
    mac_1399_default = False
    if args.case_list is None and MAC_1399_SENTINEL.exists():
        args.case_list = MAC_1399_CASE_LIST
        mac_1399_default = True

    args.out_dir.mkdir(parents=True, exist_ok=True)
    results_path = args.out_dir / RESULTS_CSV
    states = bench.read_branch_states(bench.BRANCH_STATUS)
    available = bench.openmc_materials_from_xml(bench.OPENMC_DATA / "cross_sections.xml")
    coverage = bench.read_inventory_coverage(bench.INVENTORY_CSV, available)
    runnable = [state for state in states if bench.state_openmc_supported(state, coverage)]
    runnable.sort(key=lambda state: (state.burnup, state.tfuel, state.tcool, state.tmod, state.case_id))
    allowed_case_ids: set[str] | None = None
    if args.case_list is not None:
        allowed_case_ids = set()
        with args.case_list.open(newline="", encoding="utf-8") as handle:
            sample = handle.read(4096)
            handle.seek(0)
            if "," in sample.splitlines()[0]:
                for row in csv.DictReader(handle):
                    case_id = (row.get("case_id") or "").strip()
                    if case_id:
                        allowed_case_ids.add(case_id)
            else:
                for line in handle:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        allowed_case_ids.add(line.split(",", 1)[0].strip())
        runnable = [state for state in runnable if state.case_id in allowed_case_ids]

    existing = read_existing(results_path)
    done_rows = list(existing.values())
    todo = []
    for state in runnable:
        prior = existing.get(state.case_id)
        if args.rerun_all:
            todo.append(state)
        elif prior is None:
            todo.append(state)
        elif args.rerun_failed and prior.get("acceptance_status") != "PASS":
            todo.append(state)
    if args.max_cases > 0:
        todo = todo[: args.max_cases]

    def payload_for(state: bench.BranchState) -> tuple[dict[str, object], str, int, int, int, int, float, float]:
        return (
            {
                "case_id": state.case_id,
                "burnup": state.burnup,
                "tfuel": state.tfuel,
                "tcool": state.tcool,
                "tmod": state.tmod,
                "k_eff": state.k_eff,
                "run_status": state.run_status,
            },
            str(args.out_dir / "openmc_cases"),
            args.particles,
            args.batches,
            args.inactive,
            args.threads_per_case,
            args.max_sigma,
            args.max_delta_pcm,
        )

    completed_by_id = {row["case_id"]: row for row in done_rows}
    case_index_by_id = {state.case_id: index for index, state in enumerate(runnable, start=1)}
    if todo:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = [pool.submit(run_one, payload_for(state)) for state in todo]
            for future in as_completed(futures):
                row = future.result()
                completed_by_id[str(row["case_id"])] = row
                ordered = [
                    completed_by_id[state.case_id]
                    for state in runnable
                    if state.case_id in completed_by_id
                ]
                write_results(results_path, ordered)
                write_report(
                    args.out_dir,
                    ordered,
                    len(runnable),
                    args.particles,
                    args.batches,
                    args.inactive,
                    args.max_sigma,
                    args.max_delta_pcm,
                )
                print(
                    f"{row['case_id']} {row['run_status']} {row['acceptance_status']} "
                    f"completed={len(ordered)}/{len(runnable)}",
                    flush=True,
                )
                case_index = case_index_by_id.get(str(row["case_id"]), "?")
                print(
                    f"CASE#{case_index} completed. Files written to {row['case_dir']}",
                    flush=True,
                )
    else:
        ordered = [
            completed_by_id[state.case_id]
            for state in runnable
            if state.case_id in completed_by_id
        ]
        write_report(
            args.out_dir,
            ordered,
            len(runnable),
            args.particles,
            args.batches,
            args.inactive,
            args.max_sigma,
            args.max_delta_pcm,
        )

    final_rows = read_existing(results_path)
    complete = len(final_rows) == len(runnable)
    all_accept = complete and all(row.get("acceptance_status") == "PASS" for row in final_rows.values())
    print(args.out_dir / REPORT_MD)
    print(f"completed={len(final_rows)} total={len(runnable)} all_accept={all_accept}")
    if mac_1399_default and complete:
        return 0
    return 0 if all_accept else 2


if __name__ == "__main__":
    raise SystemExit(main())
