#!/usr/bin/env python3
"""Verify the OpenMC full-grid PURR CANDU data library."""

from __future__ import annotations

import os
from pathlib import Path

import h5py
import openmc.data


def first_existing(*paths: Path) -> Path:
    for path in paths:
        if path.exists():
            return path
    return paths[0]


ROOT = Path(os.environ.get(
    "OPENMC_CANDU_OUT_DIR",
    str(first_existing(
        Path("/Volumes/Alok_Jha_Orico/OpenMC/endfb-viii.1-hdf5-candu-fullgrid-purr"),
        Path("/Users/alokjha/OpenMC/endfb-viii.1-hdf5-candu-fullgrid-purr"),
    )),
))
REPORT = ROOT / "OPENMC_CANDU_FULLGRID_PURR_VERIFICATION_REPORT.md"
XS_XML = ROOT / "cross_sections.xml"

NEUTRONS = [
    "U235", "U238", "H1", "H2", "He4", "O16",
    "Cr50", "Cr52", "Cr53", "Cr54",
    "Fe54", "Fe56", "Fe57", "Fe58",
    "Ni58", "Ni60", "Ni61", "Ni62", "Ni64",
    "Zr90", "Zr91", "Zr92", "Zr93", "Zr94", "Zr95", "Zr96",
    "Nb93",
]
THERMALS = ["c_D_in_D2O", "c_O_in_D2O"]
FULL_TEMPS = [
    "293K", "303K", "313K", "323K", "333K", "343K", "353K", "363K",
    "393K", "493K", "593K", "693K", "793K", "893K", "900K", "993K",
    "1093K", "1193K", "1293K", "1393K", "1493K", "1593K", "1693K",
    "1793K", "1893K", "1993K", "2093K", "2193K", "2293K", "2393K",
    "2493K", "2593K",
]
THERMAL_TEMPS = [
    "284K", "294K", "300K", "324K", "350K", "374K", "400K", "424K",
    "450K", "474K", "500K", "524K", "550K", "574K", "600K", "624K",
    "650K",
]


def status(pass_condition: bool) -> str:
    return "PASS" if pass_condition else "FAIL"


def main() -> int:
    failures: list[str] = []
    lines = [
        "# OpenMC Full-Grid PURR Library Verification Report",
        "",
        f"- Library directory: `{ROOT}`",
        f"- cross_sections.xml: `{XS_XML}`",
        "",
        "## File And Registry Checks",
        "",
    ]

    library = openmc.data.DataLibrary()
    for material in NEUTRONS + THERMALS:
        path = ROOT / f"{material}.h5"
        ok = path.exists()
        lines.append(f"- `{material}.h5`: `{status(ok)}`")
        if ok:
            library.register_file(path)
        else:
            failures.append(f"Missing {path}")

    if not failures:
        library.export_to_xml(XS_XML)

    if XS_XML.exists():
        parsed = openmc.data.DataLibrary.from_xml(XS_XML)
        materials = sorted(
            item for entry in parsed.libraries for item in entry["materials"]
        )
        expected_materials = sorted(NEUTRONS + THERMALS)
        ok = materials == expected_materials
        lines.append(f"- `cross_sections.xml material set`: `{status(ok)}`")
        lines.append(f"- `cross_sections.xml entries`: `{len(parsed.libraries)}`")
        if not ok:
            failures.append(
                f"Registry materials mismatch: got {materials}, expected {expected_materials}"
            )
    else:
        lines.append("- `cross_sections.xml`: `FAIL`")
        failures.append("Missing cross_sections.xml")

    lines += ["", "## Temperature Coverage", ""]
    for material in NEUTRONS:
        path = ROOT / f"{material}.h5"
        if not path.exists():
            continue
        data = openmc.data.IncidentNeutron.from_hdf5(path)
        temps = sorted(data.temperatures, key=lambda item: int(item[:-1]))
        ok = temps == FULL_TEMPS
        lines.append(f"- `{material}` neutron temperatures: `{status(ok)}`")
        if not ok:
            failures.append(f"{material} temperatures mismatch: {temps}")

    for material in THERMALS:
        path = ROOT / f"{material}.h5"
        if not path.exists():
            continue
        data = openmc.data.ThermalScattering.from_hdf5(path)
        temps = sorted(data.temperatures, key=lambda item: int(item[:-1]))
        ok = temps == THERMAL_TEMPS
        lines.append(f"- `{material}` thermal S(a,b) temperatures: `{status(ok)}`")
        if not ok:
            failures.append(f"{material} thermal temperatures mismatch: {temps}")

    lines += ["", "## PURR / URR Checks", ""]
    for material in ("U235", "U238"):
        path = ROOT / f"{material}.h5"
        if not path.exists():
            continue
        with h5py.File(path, "r") as handle:
            urr = handle.get(f"{material}/urr")
            if urr is None:
                ok = False
                temps: list[str] = []
            else:
                temps = sorted(urr.keys(), key=lambda item: int(item[:-1]))
                ok = temps == FULL_TEMPS and all(
                    {"energy", "table"} <= set(urr[temp].keys()) for temp in temps
                )
        lines.append(f"- `{material}` unresolved probability tables: `{status(ok)}`")
        if not ok:
            failures.append(f"{material} missing/incomplete URR probability tables: {temps}")

    lines += [
        "",
        "## Disposition",
        "",
        f"- Verification failures: `{len(failures)}`",
        f"- Phase-F OpenMC data readiness: `{status(not failures)}`",
    ]
    if failures:
        lines += ["", "## Failures", ""]
        lines.extend(f"- {failure}" for failure in failures)

    REPORT.write_text("\n".join(lines) + "\n")
    print(REPORT)
    print(f"failures={len(failures)}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
