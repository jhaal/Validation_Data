#!/usr/bin/env python3
"""Generate a focused OpenMC CE library for CANDU PHWR validation.

This script processes the ENDF-B-VIII.1 files needed for the OpenMC Phase F
benchmark path:
  - U-235, U-238
  - D, O-16
  - Zr-90 through Zr-96 isotopes present in the ENDF tree
  - H-1, Nb-93, Fe, Cr, Ni, and He-4 structural/impurity nuclides
  - D in D2O and O in D2O thermal scattering tables, where available

The output is an OpenMC HDF5 data directory with a cross_sections.xml file.
"""

from __future__ import annotations

import csv
import os
import shutil
import traceback
from datetime import datetime
from pathlib import Path

import openmc.data


def first_existing(*paths: Path) -> Path:
    for path in paths:
        if path.exists():
            return path
    return paths[0]


ROOT = Path(os.environ.get(
    "OPENMC_CANDU_ROOT",
    str(first_existing(
        Path("/Volumes/Alok_Jha_Orico/OpenMC"),
        Path("/Users/alokjha/OpenMC"),
    )),
))
DATA_ROOT = Path(os.environ.get(
    "OPENMC_CANDU_DATA_ROOT",
    str(first_existing(
        Path("/Volumes/Alok_Jha_Orico/OpenMC"),
        ROOT,
        Path("/Users/alokjha/OpenMC"),
    )),
))
ENDF_ROOT = DATA_ROOT / "ENDF-B-VIII.1"
NEUTRON_DIR = ENDF_ROOT / "neutrons-version.VIII.1"
THERMAL_DIR = ENDF_ROOT / "thermal_scatt-version.VIII.1"
OUT_DIR = Path(os.environ.get(
    "OPENMC_CANDU_OUT_DIR",
    str(DATA_ROOT / "endfb-viii.1-hdf5-candu"),
))
LOG_DIR = OUT_DIR / "logs"
XS_XML = OUT_DIR / "cross_sections.xml"
MANIFEST = OUT_DIR / "openmc_candu_data_manifest.csv"
REPORT = OUT_DIR / "OPENMC_CANDU_DATA_GENERATION_REPORT.md"
NJOY_EXEC = "/Users/alokjha/NJOY21/build/njoy21"


NEUTRON_FILES = [
    "n-092_U_235.endf",
    "n-092_U_238.endf",
    "n-001_H_002.endf",
    "n-001_H_001.endf",
    "n-002_He_004.endf",
    "n-008_O_016.endf",
    "n-024_Cr_050.endf",
    "n-024_Cr_052.endf",
    "n-024_Cr_053.endf",
    "n-024_Cr_054.endf",
    "n-026_Fe_054.endf",
    "n-026_Fe_056.endf",
    "n-026_Fe_057.endf",
    "n-026_Fe_058.endf",
    "n-028_Ni_058.endf",
    "n-028_Ni_060.endf",
    "n-028_Ni_061.endf",
    "n-028_Ni_062.endf",
    "n-028_Ni_064.endf",
    "n-040_Zr_090.endf",
    "n-040_Zr_091.endf",
    "n-040_Zr_092.endf",
    "n-040_Zr_093.endf",
    "n-040_Zr_094.endf",
    "n-040_Zr_095.endf",
    "n-040_Zr_096.endf",
    "n-041_Nb_093.endf",
]

THERMAL_FILES = [
    ("n-001_H_002.endf", "tsl-DinD2O.endf", "c_D_in_D2O"),
    ("n-008_O_016.endf", "tsl-OinD2O.endf", "c_O_in_D2O"),
]


def branch_temperatures() -> list[float]:
    if os.environ.get("OPENMC_CANDU_FULL_TEMPS") == "1":
        return full_branch_temperatures()
    # Bracket the production branch grid for OpenMC interpolation. The full
    # branch grid can be requested with OPENMC_CANDU_FULL_TEMPS=1, but it is a
    # long production run for actinides.
    return [293.15, 343.15, 593.15, 900.0, 1193.15, 1593.15, 1993.15, 2393.15, 2593.15]


def full_branch_temperatures() -> list[float]:
    status = Path(
        "/Users/alokjha/DRAGON/XS_FILES/"
        "phase4_dragon_37el_candu_temperature_branches/branch_status.csv"
    )
    if not status.exists():
        return [293.15]
    temps: set[float] = set()
    with status.open(newline="") as handle:
        for row in csv.DictReader(handle):
            for key in ("tfuel_K", "tcool_K", "tmod_K"):
                value = row.get(key, "")
                if value:
                    temps.add(round(float(value), 2))
    return sorted(temps)


def reset_output() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    for path in OUT_DIR.glob("*.h5"):
        path.unlink()
    for path in LOG_DIR.glob("*"):
        if path.is_file():
            path.unlink()
        else:
            shutil.rmtree(path)
    if XS_XML.exists():
        XS_XML.unlink()
    if MANIFEST.exists():
        MANIFEST.unlink()
    if REPORT.exists():
        REPORT.unlink()


def selected_neutron_files() -> list[str]:
    only = os.environ.get("OPENMC_CANDU_ONLY_NEUTRONS")
    if not only:
        return NEUTRON_FILES
    requested = {item.strip() for item in only.split(",") if item.strip()}
    selected = [item for item in NEUTRON_FILES if item in requested or Path(item).stem in requested]
    missing = requested - set(selected) - {Path(item).stem for item in selected}
    if missing:
        raise ValueError(f"Unknown neutron selections: {sorted(missing)}")
    return selected


def selected_thermal_files() -> list[tuple[str, str, str]]:
    if os.environ.get("OPENMC_CANDU_SKIP_THERMAL") == "1":
        return []
    return THERMAL_FILES


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(path)


def export_neutron(filename: str, temperatures: list[float]) -> dict[str, str]:
    source = NEUTRON_DIR / filename
    require(source)
    stem = source.stem.replace("n-", "")
    work = LOG_DIR / f"njoy_{stem}"
    work.mkdir(parents=True, exist_ok=True)
    print(
        f"[{datetime.now().isoformat(timespec='seconds')}] "
        f"processing neutron {filename} at {len(temperatures)} temperatures",
        flush=True,
    )
    data = openmc.data.IncidentNeutron.from_njoy(
        str(source),
        temperatures=temperatures,
        output_dir=str(work),
        njoy_exec=NJOY_EXEC,
        stdout=os.environ.get("OPENMC_CANDU_NJOY_STDOUT") == "1",
        heatr=False,
        gaspr=False,
        purr=os.environ.get("OPENMC_CANDU_PURR") == "1",
        input_filename=str(work / "input.njoy"),
    )
    out = OUT_DIR / f"{data.name}.h5"
    data.export_to_hdf5(out, mode="w")
    if os.environ.get("OPENMC_CANDU_CLEAN_WORK") == "1":
        shutil.rmtree(work)
    print(
        f"[{datetime.now().isoformat(timespec='seconds')}] wrote {out}",
        flush=True,
    )
    return {
        "kind": "neutron",
        "source": str(source),
        "name": data.name,
        "output": str(out),
        "temperatures": " ".join(data.temperatures),
        "status": "PASS",
        "message": "",
    }


def export_thermal(neutron_file: str, thermal_file: str, table_name: str) -> dict[str, str]:
    neutron_source = NEUTRON_DIR / neutron_file
    thermal_source = THERMAL_DIR / thermal_file
    require(neutron_source)
    require(thermal_source)
    work = LOG_DIR / f"njoy_{thermal_source.stem}"
    work.mkdir(parents=True, exist_ok=True)
    print(
        f"[{datetime.now().isoformat(timespec='seconds')}] "
        f"processing thermal {thermal_file}",
        flush=True,
    )
    data = openmc.data.ThermalScattering.from_njoy(
        str(neutron_source),
        str(thermal_source),
        output_dir=str(work),
        njoy_exec=NJOY_EXEC,
        stdout=False,
        table_name=table_name,
    )
    out = OUT_DIR / f"{table_name}.h5"
    data.export_to_hdf5(out, mode="w")
    if os.environ.get("OPENMC_CANDU_CLEAN_WORK") == "1":
        shutil.rmtree(work)
    print(
        f"[{datetime.now().isoformat(timespec='seconds')}] wrote {out}",
        flush=True,
    )
    return {
        "kind": "thermal",
        "source": f"{neutron_source}; {thermal_source}",
        "name": data.name,
        "output": str(out),
        "temperatures": " ".join(data.temperatures),
        "status": "PASS",
        "message": "",
    }


def register_library(records: list[dict[str, str]]) -> None:
    library = openmc.data.DataLibrary()
    if os.environ.get("OPENMC_CANDU_REGISTER_EXISTING") == "1":
        for path in sorted(OUT_DIR.glob("*.h5")):
            library.register_file(path)
    else:
        for record in records:
            if record["status"] == "PASS":
                library.register_file(record["output"])
    library.export_to_xml(XS_XML)


def write_manifest(records: list[dict[str, str]]) -> None:
    fields = ["kind", "source", "name", "output", "temperatures", "status", "message"]
    with MANIFEST.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def write_report(records: list[dict[str, str]], temperatures: list[float]) -> None:
    failures = [r for r in records if r["status"] != "PASS"]
    purr_enabled = os.environ.get("OPENMC_CANDU_PURR") == "1"
    lines = [
        "# OpenMC CANDU CE Data Generation Report",
        "",
        f"- ENDF root: `{ENDF_ROOT}`",
        f"- Output directory: `{OUT_DIR}`",
        f"- cross_sections.xml: `{XS_XML}`",
        f"- NJOY executable: `{NJOY_EXEC}`",
        f"- Neutron processing temperatures: `{', '.join(f'{t:g}' for t in temperatures)} K`",
        "- Temperature mode: `bracketing grid for OpenMC interpolation` "
        "(set `OPENMC_CANDU_FULL_TEMPS=1` for the exact full branch grid).",
        f"- Records processed: `{len(records)}`",
        f"- Failures: `{len(failures)}`",
        "",
        "## Coverage",
        "",
    ]
    for record in records:
        lines.append(
            f"- `{record['kind']}` `{record['name']}`: `{record['status']}` "
            f"-> `{record['output']}`"
        )
    lines += [
        "",
        "## Notes",
        "",
        "- D2O thermal scattering is registered from the available ENDF thermal "
        "files. The available S(a,b) temperatures are limited by ENDF content; "
        "OpenMC will use its normal thermal-scattering temperature treatment.",
        f"- PURR unresolved-resonance probability tables enabled: `{str(purr_enabled).lower()}`.",
        "- The neutron temperature grid brackets the production matrix but is "
        "not the exact full branch grid unless `OPENMC_CANDU_FULL_TEMPS=1` is "
        "used.",
        "- This library is a focused Phase F validation library, not a full "
        "reactor depletion library.",
    ]
    if failures:
        lines += ["", "## Failures", ""]
        for record in failures:
            lines.append(f"- `{record['source']}`: {record['message']}")
    REPORT.write_text("\n".join(lines) + "\n")


def main() -> int:
    require(Path(NJOY_EXEC))
    require(ENDF_ROOT)
    temperatures = branch_temperatures()
    if os.environ.get("OPENMC_CANDU_PRESERVE_OUTPUT") != "1":
        reset_output()
    else:
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        LOG_DIR.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, str]] = []

    for filename in selected_neutron_files():
        try:
            records.append(export_neutron(filename, temperatures))
        except Exception as exc:  # noqa: BLE001 - retain conversion failure evidence
            output = getattr(exc, "output", "")
            records.append({
                "kind": "neutron",
                "source": str(NEUTRON_DIR / filename),
                "name": filename,
                "output": "",
                "temperatures": "",
                "status": "FAIL",
                "message": (
                    f"{type(exc).__name__}: {exc}\n"
                    f"--- NJOY output ---\n{output}\n"
                    f"--- traceback ---\n{traceback.format_exc()}"
                ),
            })

    for neutron_file, thermal_file, table_name in selected_thermal_files():
        try:
            records.append(export_thermal(neutron_file, thermal_file, table_name))
        except Exception as exc:  # noqa: BLE001
            output = getattr(exc, "output", "")
            records.append({
                "kind": "thermal",
                "source": f"{NEUTRON_DIR / neutron_file}; {THERMAL_DIR / thermal_file}",
                "name": table_name,
                "output": "",
                "temperatures": "",
                "status": "FAIL",
                "message": (
                    f"{type(exc).__name__}: {exc}\n"
                    f"--- NJOY output ---\n{output}\n"
                    f"--- traceback ---\n{traceback.format_exc()}"
                ),
            })

    write_manifest(records)
    if any(record["status"] != "PASS" for record in records):
        write_report(records, temperatures)
        return 1

    register_library(records)
    write_report(records, temperatures)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
