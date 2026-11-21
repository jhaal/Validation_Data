#!/usr/bin/env python3
"""Shared DRAGON deck helpers for the CANDU37 library.

step_deck() writes the stepwise nominal depletion (0.5-GWd/tU steps) and, at the burnups in
BRANCH_BU, create-mode temperature branches on the depleted inventory: the branch library LIBBR
receives the depleted densities through EVO: ... NOSA SET, and LIBB2 := USS: LIBBR ... rebuilds
every isotope at the branch temperatures. install() writes a deck, its access script and procedures.
Used by make_matched_inventory_decks.py, make_corrected_matrix.py, make_corrected_2group_and_sweeps.py
and make_mesh_probe.py.
"""
from __future__ import annotations

import pathlib
import shutil
import sys

DRAGON = pathlib.Path("/Users/alokjha/DRAGON/5.1/Dragon")
XS = pathlib.Path("/Users/alokjha/DRAGON/XS_FILES")
PROC = XS / "phase4_dragon_37el_candu_temperature_branches" / "proc"

LIBS = {
    "J2": "draglibJef2p2_v5p1",
    "E81": "draglibendfb8r1_v5p1",
}
POWER = 36.8
NOM = (900.0, 343.15, 343.15)
STEP_BU = [0.5 * i for i in range(41)]
BRANCH_BU = [0.0, 5.0, 10.0, 15.0, 20.0]
BRANCHES = [
    ("fuel", 293.15, 343.15, 343.15),
    ("fuel", 900.0, 343.15, 343.15),
    ("fuel", 1593.15, 343.15, 343.15),
    ("fuel", 1993.15, 343.15, 343.15),
    ("mod", 900.0, 343.15, 293.15),
    ("mod", 900.0, 343.15, 363.15),
    ("cool", 900.0, 293.15, 343.15),
    ("cool", 900.0, 363.15, 343.15),
]


def access(case: str, libkey: str) -> str:
    lib = LIBS[libkey]
    return f"""#!/bin/sh
pos=$1/../../libraries/l_endian
if [ -f "$pos"/{lib}.gz ] && [ ! -f "$pos"/{lib} ]; then
  gunzip -k "$pos"/{lib}.gz
fi
ln -s "$pos"/{lib} DLIB_{libkey}
ln -sf "$1"/data/{case}_proc/*.c2m .
echo "{case} access script terminated"
"""


def days(bu: float) -> float:
    return bu * 1000.0 / POWER


def step_deck(libkey: str, modes: str = "B", only_bu: float | None = None) -> str:
    """modes is kept for caller compatibility; only create-mode branches ("B") are written."""
    assert modes == "B", modes
    L = [
        "* Stepwise depletion and create-mode temperature branches",
        "MODULE NXT: USS: ASM: FLU: EVO: GREP: DELETE: END: ;",
        "PROCEDURE GEOM37_NXT MATDPLD2O ;",
        "LINKED_LIST GEOM GEOMT LIBRAW LIBSS SYS FLUX BURNUP LIBBR LIBB2 SYSB FLUXB ;",
        "SEQ_BINARY TRACK ;",
        "REAL Keff ;",
        f'STRING LibName := "DLIB_{libkey}" ;',
        "GEOM := GEOM37_NXT :: ;",
        "GEOMT TRACK := NXT: GEOM :: ALLG TISO 3 3.0 ;",
        f"LIBRAW := MATDPLD2O :: {NOM[0]:.5f} {NOM[1]:.5f} {NOM[2]:.5f} <<LibName>> ;",
        "LIBSS := USS: LIBRAW GEOMT TRACK :: EDIT 0 GRMIN 0 PIJ ;",
        "SYS := ASM: LIBSS GEOMT TRACK :: PIJ EDIT 0 ;",
        "FLUX := FLU: SYS LIBSS GEOMT TRACK :: EDIT 0 TYPE K THER 5 5E-5 ;",
        "GREP: FLUX :: GETVAL 'K-EFFECTIVE' 1 >>Keff<< ;",
        "ECHO 'VRF_NOMINAL' 0.0 Keff ;",
    ]
    prev = 0.0
    for i, bu in enumerate(STEP_BU):
        if i > 0:
            t0, t1 = days(prev), days(bu)
            src = "LIBSS FLUX GEOMT" if i == 1 else "BURNUP LIBSS FLUX GEOMT"
            L += [
                f"BURNUP LIBSS := EVO: {src} ::",
                "  EDIT 0",
                f"  SAVE {t0:.8f} DAY POWR {POWER:.8f}",
                *(["  NOEX"] if i > 1 else []),
                f"  DEPL {t0:.8f} {t1:.8f} DAY POWR {POWER:.8f}",
                f"  SET {t1:.8f} DAY ;",
                "LIBSS := USS: LIBSS LIBRAW GEOMT TRACK ;",
                "SYS := DELETE: SYS ;",
                "SYS := ASM: LIBSS GEOMT TRACK :: PIJ EDIT 0 ;",
                "FLUX := FLU: FLUX SYS LIBSS GEOMT TRACK :: EDIT 0 TYPE K THER 5 5E-5 ;",
                "GREP: FLUX :: GETVAL 'K-EFFECTIVE' 1 >>Keff<< ;",
                f"ECHO 'VRF_NOMINAL' {bu:.1f} Keff ;",
                f"BURNUP LIBSS := EVO: BURNUP LIBSS FLUX GEOMT :: SAVE {t1:.8f} DAY POWR {POWER:.8f} ;",
            ]
        if bu in BRANCH_BU and (only_bu is None or bu == only_bu):
            t = days(bu)
            for tag, tf, tc, tm in BRANCHES:
                lab = f"{bu:.1f} {tf:.2f} {tc:.2f} {tm:.2f}"
                # Create-mode branch on the depleted inventory.
                L += [
                    f"LIBBR := MATDPLD2O :: {tf:.5f} {tc:.5f} {tm:.5f} <<LibName>> ;",
                ]
                if bu > 0.0:
                    L += [f"BURNUP LIBBR := EVO: BURNUP LIBBR FLUX GEOMT :: EDIT 0 NOSA SET {t:.8f} DAY ;"]
                L += [
                    "LIBB2 := USS: LIBBR GEOMT TRACK :: EDIT 0 GRMIN 0 PIJ ;",
                    "SYSB := ASM: LIBB2 GEOMT TRACK :: PIJ EDIT 0 ;",
                    "FLUXB := FLU: SYSB LIBB2 GEOMT TRACK :: EDIT 0 TYPE K THER 5 5E-5 ;",
                    "GREP: FLUXB :: GETVAL 'K-EFFECTIVE' 1 >>Keff<< ;",
                    f"ECHO 'VRF_BRANCH_B' {lab} Keff ;",
                    "SYSB FLUXB LIBBR LIBB2 := DELETE: SYSB FLUXB LIBBR LIBB2 ;",
                ]
        prev = bu
        if only_bu is not None and bu >= only_bu:
            break
    L += ["ECHO 'VRF_COMPLETED' ;", "END: ;", 'QUIT "LIST" .', ""]
    return "\n".join(L)


def install(case: str, text: str, libkey: str) -> None:
    (DRAGON / "data" / f"{case}.x2m").write_text(text)
    acc = DRAGON / "data" / f"{case}.access"
    acc.write_text(access(case, libkey))
    acc.chmod(0o755)
    pdir = DRAGON / "data" / f"{case}_proc"
    pdir.mkdir(exist_ok=True)
    shutil.copy2(PROC / "GEOM37_NXT.c2m", pdir / "GEOM37_NXT.c2m")
    mat = (PROC / "MATDPLD2O.c2m").read_text()
    if libkey == "E81":
        mat = zr_isotopic(mat)
    (pdir / "MATDPLD2O.c2m").write_text(mat)
    print(case)


# Natural Zr atom fractions (IUPAC); ENDF/B-VIII.1 DRAGLIB has no elemental Zr0 entry.
ZR_NAT = [("Zr90", 0.5145), ("Zr91", 0.1122), ("Zr92", 0.1715), ("Zr94", 0.1738), ("Zr96", 0.0280)]


def zr_isotopic(mat: str) -> str:
    out = []
    for line in mat.splitlines():
        s = line.split()
        if len(s) >= 4 and s[0] == "ZR" and s[2] == "<<IsoZr>>":
            n = float(s[3])
            tail = " ".join(s[4:])
            for j, (iso, f) in enumerate(ZR_NAT):
                out.append(f'    ZR{j} = "{iso}" {n * f:.5E} {tail}')
        else:
            out.append(line)
    return "\n".join(out) + "\n"


def main() -> None:
    """Usage: python make_verification_decks.py step:J2 [step:E81]  (stepwise depletion with branches)."""
    for arg in sys.argv[1:]:
        if arg.startswith("step:"):
            key = arg.split(":", 1)[1]
            install(f"VRF_STEPB_{key}", step_deck(key), key)


if __name__ == "__main__":
    main()
