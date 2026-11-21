#!/usr/bin/env python3
"""Persist and export DRAGON Phase 3 depletion inventories.

The original Phase 3 run kept the BURNUP object in memory.  This reruns the
same nominal depletion with BURNUP declared as an external XSM object, then
exports the saved DEPL-DATnnnn/ISOTOPESDENS records to CSV for OpenMC Phase F.
"""

from __future__ import annotations

import csv
import os
import pathlib
import shutil
import subprocess
import sys
from datetime import datetime, timezone

import numpy as np

import implement_phase3_dragon_depletion as phase3


def first_existing(*paths: pathlib.Path) -> pathlib.Path:
    for path in paths:
        if path.exists():
            return path
    return paths[0]


DRAGON_ROOT = pathlib.Path(os.environ.get(
    "DRAGON_ROOT",
    str(first_existing(
        pathlib.Path("/Volumes/Alok_Jha_Orico/DRAGON"),
        pathlib.Path("/Users/alokjha/DRAGON"),
    )),
))
XS = pathlib.Path(os.environ.get("DRAGON_XS_ROOT", str(DRAGON_ROOT / "XS_FILES")))
DRAGON = pathlib.Path(os.environ.get("DRAGON_CODE_ROOT", str(DRAGON_ROOT / "5.1" / "Dragon")))
PKG = XS / "phase3_dragon_37el_candu_depletion"
CASE = "CANDU37_NU_PHASE3_DEPL_XSM"
XSM_NAME = "CANDU37_NU_PHASE3_BURNUP.xsm"
BURN_ASCII_NAME = "BURNEXP"
PYGAN = pathlib.Path(os.environ.get(
    "PYGAN_PATH",
    str(DRAGON_ROOT / "5.1" / "PyGan" / "lib" / "Darwin_arm64" / "python"),
))


def write(path: pathlib.Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def build_xsm_deck() -> str:
    deck = phase3.build_deck()
    deck = deck.replace(
        "LINKED_LIST GEOM GEOMT LIBRAW LIBSS SYS FLUX BURNUP ;",
        (
            "LINKED_LIST GEOM GEOMT LIBRAW LIBSS SYS FLUX ;\n"
            f"XSM_FILE BURNUP :: FILE './{XSM_NAME}' ;\n"
            f"SEQ_ASCII {BURN_ASCII_NAME} :: FILE './{BURN_ASCII_NAME}' ;"
        ),
    )
    deck = deck.replace(
        "BURNUP LIBSS LIBRAW SYS FLUX GEOM GEOMT TRACK := DELETE: BURNUP LIBSS LIBRAW SYS FLUX GEOM GEOMT TRACK ;",
        (
            f"{BURN_ASCII_NAME} := BURNUP :: EDIT 10 ;\n"
            "LIBSS LIBRAW SYS FLUX GEOM GEOMT TRACK := DELETE: LIBSS LIBRAW SYS FLUX GEOM GEOMT TRACK ;"
        ),
    )
    return deck


def decode_isotope_names(raw: np.ndarray) -> list[str]:
    if isinstance(raw, str):
        return [raw[i : i + 12].strip() for i in range(0, len(raw), 12) if raw[i : i + 12].strip()]
    arr = np.asarray(raw)
    if arr.shape == () and arr.dtype.kind in {"S", "U"}:
        text = str(arr.item())
        return [text[i : i + 12].strip() for i in range(0, len(text), 12) if text[i : i + 12].strip()]
    if arr.dtype.kind in {"S", "U"}:
        return [str(item).strip().replace("b'", "").strip("'") for item in arr]
    flat = arr.reshape(-1)
    if flat.size % 3 != 0:
        return [str(item) for item in flat]
    names: list[str] = []
    for i in range(0, flat.size, 3):
        chunks = []
        for word in flat[i : i + 3]:
            value = int(word)
            b = value.to_bytes(4, byteorder=sys.byteorder, signed=False)
            text = b.decode("ascii", errors="ignore")
            if not text.strip():
                b = value.to_bytes(4, byteorder="big", signed=False)
                text = b.decode("ascii", errors="ignore")
            chunks.append(text)
        names.append("".join(chunks).strip())
    return names


def load_lcm():
    sys.path.insert(0, str(PYGAN))
    import lcm  # type: ignore[import-not-found]

    return lcm


def find_xsm() -> pathlib.Path:
    candidates = [
        DRAGON / XSM_NAME,
        DRAGON / "Darwin_arm64" / XSM_NAME,
        DRAGON / "data" / XSM_NAME,
        PKG / "burnup_isotopics" / XSM_NAME,
    ]
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError(f"Could not locate {XSM_NAME}; checked {candidates}")


def build_save_script(iso_dir: pathlib.Path) -> str:
    return f"""#!/bin/sh
set -eu
dest='{iso_dir}'
mkdir -p "$dest"
if [ -f './{XSM_NAME}' ]; then
  cp './{XSM_NAME}' "$dest/{XSM_NAME}"
else
  ls -la > "$dest/{CASE}.save_directory_listing.txt"
  echo "Missing ./{XSM_NAME} in DRAGON run directory" >&2
  exit 3
fi
if [ -f './{BURN_ASCII_NAME}' ]; then
  cp './{BURN_ASCII_NAME}' "$dest/{BURN_ASCII_NAME}"
else
  ls -la > "$dest/{CASE}.save_directory_listing.txt"
  echo "Missing ./{BURN_ASCII_NAME} in DRAGON run directory" >&2
  exit 4
fi
"""


def export_inventory(xsm_path: pathlib.Path, out_dir: pathlib.Path) -> dict[str, object]:
    lcm = load_lcm()
    ascii_path = out_dir / BURN_ASCII_NAME
    if not ascii_path.exists():
        raise FileNotFoundError(ascii_path)
    short_name = "burnexp"
    short_ascii = pathlib.Path("/tmp") / f"_{short_name}"
    shutil.copy2(ascii_path, short_ascii)
    cwd = pathlib.Path.cwd()
    try:
        os.chdir("/tmp")
        burn = lcm.new("LCM_INP", short_name, impx=0)
    finally:
        os.chdir(cwd)
    state = np.asarray(burn["STATE-VECTOR"])
    mixes = np.asarray(burn["ISOTOPESMIX"]).astype(int)
    names = decode_isotope_names(burn["ISOTOPESUSED"])
    dragon_times_1e8s = np.asarray(burn["DEPL-TIMES"], dtype=float)
    times_days = dragon_times_1e8s * 1.0e8 / 86400.0
    nbiso = len(names)
    if mixes.size != nbiso:
        raise ValueError(f"ISOTOPESMIX length {mixes.size} != names length {nbiso}")

    out_dir.mkdir(parents=True, exist_ok=True)
    long_csv = out_dir / "phase3_dragon_burnup_isotopic_densities_long.csv"
    fuel_csv = out_dir / "phase3_dragon_fuel_isotopic_densities_openmc.csv"
    summary_csv = out_dir / "phase3_dragon_inventory_summary.csv"

    burnups = phase3.BURNUPS
    with long_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "state_index",
                "burnup_GWd_tU",
                "time_days",
                "mix",
                "isotope",
                "atom_density_barn_cm",
            ]
        )
        for idx, (time_days, burnup) in enumerate(zip(times_days, burnups), start=1):
            dat = burn[f"DEPL-DAT{idx:04d}"]
            dens = np.asarray(dat["ISOTOPESDENS"], dtype=float)
            if dens.size != nbiso:
                raise ValueError(
                    f"DEPL-DAT{idx:04d} density length {dens.size} != {nbiso}"
                )
            for mix, isotope, density in zip(mixes, names, dens):
                writer.writerow(
                    [
                        idx - 1,
                        f"{burnup:.8e}",
                        f"{time_days:.8e}",
                        int(mix),
                        isotope,
                        f"{density:.12e}",
                    ]
                )

    with fuel_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "state_index",
                "burnup_GWd_tU",
                "time_days",
                "fuel_zone",
                "isotope",
                "atom_density_barn_cm",
            ]
        )
        for idx, (time_days, burnup) in enumerate(zip(times_days, burnups), start=1):
            dat = burn[f"DEPL-DAT{idx:04d}"]
            dens = np.asarray(dat["ISOTOPESDENS"], dtype=float)
            for mix, isotope, density in zip(mixes, names, dens):
                if int(mix) in {11, 12, 13, 14} and density > 0.0:
                    writer.writerow(
                        [
                            idx - 1,
                            f"{burnup:.8e}",
                            f"{time_days:.8e}",
                            int(mix),
                            isotope,
                            f"{density:.12e}",
                        ]
                    )

    with summary_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["key", "value"])
        writer.writerow(["xsm_path", str(xsm_path)])
        writer.writerow(["ascii_export_path", str(ascii_path)])
        writer.writerow(["state_vector", " ".join(str(int(x)) for x in state)])
        writer.writerow(["number_isotopes", nbiso])
        writer.writerow(["number_saved_times", len(times_days)])
        writer.writerow(["first_dragon_time_1e8s", f"{dragon_times_1e8s[0]:.12e}"])
        writer.writerow(["last_dragon_time_1e8s", f"{dragon_times_1e8s[-1]:.12e}"])
        writer.writerow(["first_time_days", f"{times_days[0]:.12e}"])
        writer.writerow(["last_time_days", f"{times_days[-1]:.12e}"])
        writer.writerow(["long_csv", str(long_csv)])
        writer.writerow(["fuel_csv", str(fuel_csv)])

    return {
        "nbiso": nbiso,
        "ntimes": len(times_days),
        "long_csv": long_csv,
        "fuel_csv": fuel_csv,
        "summary_csv": summary_csv,
    }


def main() -> int:
    deck_dir = PKG / "decks"
    proc_dir = PKG / "proc"
    log_dir = PKG / "logs"
    iso_dir = PKG / "burnup_isotopics"
    for directory in (deck_dir, proc_dir, log_dir, iso_dir):
        directory.mkdir(parents=True, exist_ok=True)

    phase2_geom = phase3.PHASE2 / "proc" / "GEOM37_NXT.c2m"
    if not phase2_geom.exists():
        raise FileNotFoundError(phase2_geom)

    deck = build_xsm_deck()
    access = phase3.ACCESS.replace(phase3.CASE, CASE)

    write(deck_dir / f"{CASE}.x2m", deck)
    write(deck_dir / f"{CASE}.access", access)
    (deck_dir / f"{CASE}.access").chmod(0o755)
    write(proc_dir / "MATDPLD2O.c2m", phase3.MAT_DEPL)
    shutil.copy2(phase2_geom, proc_dir / "GEOM37_NXT.c2m")

    dragon_deck = DRAGON / "data" / f"{CASE}.x2m"
    dragon_access = DRAGON / "data" / f"{CASE}.access"
    dragon_save = DRAGON / "data" / f"{CASE}.save"
    dragon_proc = DRAGON / "data" / f"{CASE}_proc"
    write(dragon_deck, deck)
    write(dragon_access, access)
    dragon_access.chmod(0o755)
    write(dragon_save, build_save_script(iso_dir))
    dragon_save.chmod(0o755)
    dragon_proc.mkdir(parents=True, exist_ok=True)
    shutil.copy2(phase2_geom, dragon_proc / "GEOM37_NXT.c2m")
    write(dragon_proc / "MATDPLD2O.c2m", phase3.MAT_DEPL)

    for stale in (DRAGON / XSM_NAME, DRAGON / "Darwin_arm64" / XSM_NAME):
        stale.unlink(missing_ok=True)

    env = os.environ.copy()
    env["PYTHONPATH"] = str(PYGAN)
    proc = subprocess.run(
        ["./rdragon", f"{CASE}.x2m"],
        cwd=DRAGON,
        text=True,
        capture_output=True,
        timeout=2400,
        check=False,
        env=env,
    )
    write(log_dir / f"{CASE}.stdout_stderr.log", proc.stdout + proc.stderr)
    result_src = DRAGON / "Darwin_arm64" / f"{CASE}.result"
    result_text = ""
    if result_src.exists():
        result_text = result_src.read_text(encoding="utf-8", errors="ignore")
        shutil.copy2(result_src, log_dir / f"{CASE}.result")

    normal = "normal end of execution for dragon 5" in result_text
    completed = "PHASE3 CANDU37 depletion completed" in result_text
    if proc.returncode != 0 or not normal or not completed:
        report = [
            "# Phase 3 DRAGON Inventory Export",
            "",
            f"Generated UTC: {datetime.now(timezone.utc).isoformat()}",
            f"rdragon exit code: `{proc.returncode}`",
            f"DRAGON normal end: `{normal}`",
            f"Completion marker: `{completed}`",
            "",
            "Inventory export failed before XSM extraction.",
        ]
        write(iso_dir / "phase3_dragon_inventory_export_report.md", "\n".join(report) + "\n")
        print(iso_dir / "phase3_dragon_inventory_export_report.md")
        return 2

    xsm_path = find_xsm()
    copied_xsm = iso_dir / XSM_NAME
    if xsm_path != copied_xsm:
        shutil.copy2(xsm_path, copied_xsm)
        xsm_path = copied_xsm

    info = export_inventory(xsm_path, iso_dir)
    report = [
        "# Phase 3 DRAGON Inventory Export",
        "",
        f"Generated UTC: {datetime.now(timezone.utc).isoformat()}",
        f"DRAGON case: `{CASE}`",
        f"rdragon exit code: `{proc.returncode}`",
        f"DRAGON normal end: `{'PASS' if normal else 'FAIL'}`",
        f"Completion marker: `{'PASS' if completed else 'FAIL'}`",
        f"BURNUP XSM object: `{xsm_path}`",
        f"Saved depletion states: `{info['ntimes']}`",
        f"Isotopes per state: `{info['nbiso']}`",
        "",
        "## Outputs",
        "",
        f"- Long isotope-density CSV: `{info['long_csv']}`",
        f"- Fuel-zone OpenMC inventory CSV: `{info['fuel_csv']}`",
        f"- Summary CSV: `{info['summary_csv']}`",
        "",
        "## Notes",
        "",
        "- Densities are exported directly from DRAGON `DEPL-DATnnnn/ISOTOPESDENS`.",
        "- Units are DRAGON microlib atom density units, treated as atom/barn-cm for OpenMC material reconstruction.",
        "- Fuel-zone CSV is restricted to depleted fuel mixtures 11, 12, 13, and 14.",
    ]
    write(iso_dir / "phase3_dragon_inventory_export_report.md", "\n".join(report) + "\n")
    print(iso_dir / "phase3_dragon_inventory_export_report.md")
    print(f"states={info['ntimes']} isotopes={info['nbiso']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
