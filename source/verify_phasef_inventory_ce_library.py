#!/usr/bin/env python3
"""Read-only verification of Phase F inventory CE coverage."""

from __future__ import annotations

import collections
import csv
import os
from datetime import datetime, timezone
from pathlib import Path

import openmc.data


XS_ROOT = Path(os.environ.get("DRAGON_XS_ROOT", "/Volumes/Alok_Jha_Orico/DRAGON/XS_FILES"))
PHASEF_DIR = XS_ROOT / "phaseF_openmc_fullgrid_purr_benchmark"
MANIFEST = PHASEF_DIR / "phaseF_missing_ce_processing_manifest.csv"
REPORT = PHASEF_DIR / "PHASE_F_INVENTORY_CE_LIBRARY_VERIFICATION_REPORT.md"
ROOT = Path(os.environ.get(
    "OPENMC_CANDU_OUT_DIR",
    "/Volumes/Alok_Jha_Orico/OpenMC/endfb-viii.1-hdf5-candu-fullgrid-purr",
))
XS_XML = ROOT / "cross_sections.xml"
FULL_TEMPS = [
    "293K", "303K", "313K", "323K", "333K", "343K", "353K", "363K",
    "393K", "493K", "593K", "693K", "793K", "893K", "900K", "993K",
    "1093K", "1193K", "1293K", "1393K", "1493K", "1593K", "1693K",
    "1793K", "1893K", "1993K", "2093K", "2193K", "2293K", "2393K",
    "2493K", "2593K",
]


def status(ok: bool) -> str:
    return "PASS" if ok else "FAIL"


def temp_key(label: str) -> int:
    return int(label[:-1])


def read_manifest() -> list[dict[str, str]]:
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def registered_materials() -> set[str]:
    data = openmc.data.DataLibrary.from_xml(XS_XML)
    return {material for entry in data.libraries for material in entry["materials"]}


def main() -> int:
    failures: list[str] = []
    rows = read_manifest()
    counts = collections.Counter(row["status"] for row in rows)
    materials = registered_materials() if XS_XML.exists() else set()

    lines = [
        "# Phase F Inventory CE Library Verification Report",
        "",
        f"Generated UTC: `{datetime.now(timezone.utc).isoformat()}`",
        f"Manifest: `{MANIFEST}`",
        f"Library directory: `{ROOT}`",
        f"cross_sections.xml: `{XS_XML}`",
        "",
        "## Manifest Checks",
        "",
        f"- Manifest rows: `{len(rows)}`",
        f"- PASS rows: `{counts.get('PASS', 0)}`",
        f"- PRESENT rows: `{counts.get('PRESENT', 0)}`",
        f"- Other rows: `{sum(value for key, value in counts.items() if key not in {'PASS', 'PRESENT'})}`",
    ]

    if len(rows) != 231:
        failures.append(f"Manifest row count is {len(rows)}, expected 231")
    for row in rows:
        if row["status"] not in {"PASS", "PRESENT"}:
            failures.append(f"{row['nuclide']} has status {row['status']}: {row['message']}")

    lines += ["", "## Registry And Temperature Checks", ""]
    checked = 0
    for row in rows:
        nuclide = row["nuclide"]
        path = ROOT / f"{nuclide}.h5"
        if not path.exists():
            failures.append(f"{nuclide} missing HDF5 file {path}")
            continue
        if nuclide not in materials:
            failures.append(f"{nuclide} missing from cross_sections.xml material registry")
            continue
        data = openmc.data.IncidentNeutron.from_hdf5(path)
        temps = sorted(data.temperatures, key=temp_key)
        if temps != FULL_TEMPS:
            failures.append(f"{nuclide} temperature grid mismatch: {temps}")
            continue
        checked += 1

    expected_registered = len(list(ROOT.glob("*.h5")))
    registry_entries = len(openmc.data.DataLibrary.from_xml(XS_XML).libraries) if XS_XML.exists() else 0
    registry_ok = registry_entries == expected_registered
    if not registry_ok:
        failures.append(
            f"cross_sections.xml entries {registry_entries}, expected one per HDF5 file {expected_registered}"
        )

    lines += [
        f"- HDF5 files in library directory: `{expected_registered}`",
        f"- cross_sections.xml entries: `{registry_entries}`",
        f"- One registry entry per HDF5 file: `{status(registry_ok)}`",
        f"- Inventory nuclides with full 32-temperature CE data and registry coverage: `{checked}`",
    ]

    lines += ["", "## Disposition", ""]
    lines.append(f"- Verification failures: `{len(failures)}`")
    lines.append(f"- Inventory CE coverage readiness: `{status(not failures)}`")
    if failures:
        lines += ["", "## Failures", ""]
        lines.extend(f"- {failure}" for failure in failures)

    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(REPORT)
    print(f"failures={len(failures)}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
