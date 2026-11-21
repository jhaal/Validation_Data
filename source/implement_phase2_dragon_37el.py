#!/usr/bin/env python3
"""Implement Phase 2 DRAGON 37-element CANDU/PHWR fresh lattice artifacts."""

from __future__ import annotations

import pathlib
import re
import shutil
import subprocess
from datetime import datetime, timezone


DRAGON = pathlib.Path("/Users/alokjha/DRAGON/5.1/Dragon")
XS = pathlib.Path("/Users/alokjha/DRAGON/XS_FILES")
PKG = XS / "phase2_dragon_37el_candu"
CASE = "CANDU37_NU_PHASE2"


GEOM = r"""*DECK GEOM37_NXT
PARAMETER GEOM ::
  EDIT 1
  ::: LINKED_LIST GEOM ;
;
MODULE GEO: ;

* Mixtures:
* 1 coolant D2O, 2 pressure tube, 3 helium gap, 4 calandria tube,
* 5 moderator D2O, 6 fuel sheath, 11-14 natural-U fuel radial zones.
REAL cell_pitch := 28.575 ;
REAL Rci := 0.6122 ;
REAL Rgi := 0.6540 ;
REAL Rci1 Rci2 Rci3 := 0.3061 0.4898 0.5816 ;
REAL Rcool := 5.1689 ;
REAL Rpt := 5.6032 ;
REAL Rgap := 6.4478 ;
REAL Rct := 6.5875 ;
REAL Rout := 14.2875 ;

GEOM := GEO: :: CARCEL 8
  EDIT 5
  X- REFL X+ REFL
  Y- REFL Y+ REFL
  MESHX 0.0 <<cell_pitch>>
  MESHY 0.0 <<cell_pitch>>
  RADIUS 0.0 0.8670 2.4545 3.7295 <<Rcool>> <<Rpt>> <<Rgap>> <<Rct>> <<Rout>>
  MIX 1 1 1 1 2 3 4 5 5
  CLUSTER ROD1 ROD2 ROD3 ROD4
  ::: ROD1 := GEO: TUBE 5
    MIX 11 12 13 14 6
    RADIUS 0.0 <<Rci1>> <<Rci2>> <<Rci3>> <<Rci>> <<Rgi>>
    NPIN 1 RPIN 0.0000 APIN 0.0000 ;
  ::: ROD2 := GEO: TUBE 5
    MIX 11 12 13 14 6
    RADIUS 0.0 <<Rci1>> <<Rci2>> <<Rci3>> <<Rci>> <<Rgi>>
    NPIN 6 RPIN 1.4880 APIN 0.0000 ;
  ::: ROD3 := GEO: TUBE 5
    MIX 11 12 13 14 6
    RADIUS 0.0 <<Rci1>> <<Rci2>> <<Rci3>> <<Rci>> <<Rgi>>
    NPIN 12 RPIN 2.8750 APIN 0.0000 ;
  ::: ROD4 := GEO: TUBE 5
    MIX 11 12 13 14 6
    RADIUS 0.0 <<Rci1>> <<Rci2>> <<Rci3>> <<Rci>> <<Rgi>>
    NPIN 18 RPIN 4.3300 APIN 0.0000 ;
;
QUIT "LIST" .
"""


MAT = r"""*DECK MAT_NU_D2O
PARAMETER LIBRARY ::
  EDIT 1
  ::: LINKED_LIST LIBRARY ;
;
REAL TfuelK TcoolK TmodK ;
STRING LibName ;
:: >>TfuelK<< >>TcoolK<< >>TmodK<< >>LibName<< ;

MODULE LIB: ;

STRING IsoU5 IsoU8 IsoO16 IsoD2 IsoH1 IsoZr IsoNb IsoFe4 IsoFe6 IsoFe7 IsoFe8
       IsoCr0 IsoCr2 IsoCr3 IsoCr4 IsoNi8 IsoNi0 IsoNi1 IsoNi2 IsoNi4 IsoHe4
  :=
       "U235" "U238" "O16" "H2_D2O" "H1_H2O" "Zr0" "Nb93" "Fe54" "Fe56" "Fe57" "Fe58"
       "Cr50" "Cr52" "Cr53" "Cr54" "Ni58" "Ni60" "Ni61" "Ni62" "Ni64" "He4" ;

* Natural uranium approximation, atom fractions in heavy metal based on 0.711% U-235.
REAL NU5 := 1.6100E-04 ;
REAL NU8 := 2.2480E-02 ;
REAL NO16F := 4.5280E-02 ;

LIBRARY := LIB: ::
  EDIT 1
  NMIX 14
  CTRA APOL
  PTSL
  ANIS 2
  MIXS LIB: DRAGON FIL: <<LibName>>
* coolant D2O
  MIX 1 <<TcoolK>>
    D2  = <<IsoD2>>  6.4886E-02
    H1  = <<IsoH1>>  6.4950E-05
    O16 = <<IsoO16>> 3.2398E-02
* pressure tube Zr-Nb approximation
  MIX 2 <<TcoolK>>
    ZR  = <<IsoZr>>  4.1780E-02 2 IRSET PTSL NONE
    NB  = <<IsoNb>>  1.0877E-03
* helium gap
  MIX 3 <<TcoolK>>
    HE4 = <<IsoHe4>> 2.0877E-04
* calandria tube
  MIX 4 <<TmodK>>
    ZR  = <<IsoZr>>  4.2364E-02 2 IRSET PTSL NONE
    FE4 = <<IsoFe4>> 5.7861E-06
    FE6 = <<IsoFe6>> 8.7588E-05
    FE7 = <<IsoFe7>> 1.9873E-06
    FE8 = <<IsoFe8>> 2.5991E-07
    CR0 = <<IsoCr0>> 3.4408E-06
    CR2 = <<IsoCr2>> 6.3805E-05
    CR3 = <<IsoCr3>> 7.0983E-06
    CR4 = <<IsoCr4>> 1.7342E-06
    NI8 = <<IsoNi8>> 2.5562E-05
    NI0 = <<IsoNi0>> 9.5185E-06
    NI1 = <<IsoNi1>> 4.0700E-07
    NI2 = <<IsoNi2>> 1.2765E-06
    NI4 = <<IsoNi4>> 3.1510E-07
* moderator D2O
  MIX 5 <<TmodK>>
    D2  = <<IsoD2>>  6.4886E-02
    H1  = <<IsoH1>>  6.4950E-05
    O16 = <<IsoO16>> 3.2398E-02
* sheath
  MIX 6 <<TcoolK>>
    ZR  = <<IsoZr>>  4.9656E-02 2 IRSET PTSL NONE
    FE4 = <<IsoFe4>> 1.0341E-05
    FE6 = <<IsoFe6>> 1.5654E-04
    FE7 = <<IsoFe7>> 3.5509E-06
    FE8 = <<IsoFe8>> 4.6451E-07
    CR0 = <<IsoCr0>> 3.9532E-06
    CR2 = <<IsoCr2>> 7.3306E-05
    CR3 = <<IsoCr3>> 8.1553E-06
    CR4 = <<IsoCr4>> 1.9925E-06
    NI8 = <<IsoNi8>> 3.7378E-06
    NI0 = <<IsoNi0>> 1.3918E-06
    NI1 = <<IsoNi1>> 5.9514E-08
    NI2 = <<IsoNi2>> 1.8668E-07
    NI4 = <<IsoNi4>> 4.6062E-08
* four radial fuel zones
  MIX 11 <<TfuelK>>
    O16 = <<IsoO16>> <<NO16F>>
    U5  = <<IsoU5>>  <<NU5>> 1 IRSET PTSL NONE
    U8  = <<IsoU8>>  <<NU8>> 1 IRSET PTSL NONE
  MIX 12 COMB 11 1.0
  MIX 13 COMB 11 1.0
  MIX 14 COMB 11 1.0
;
QUIT "LIST" .
"""


DECK = r"""*----
* Phase 2 fresh 37-element natural-uranium CANDU/PHWR lattice probe.
* This deck is a fresh nominal-temperature model only. It is not depleted
* and not yet a validated production safety library.
*----
MODULE NXT: USS: ASM: FLU: EDI: GREP: DELETE: END: ;
PROCEDURE GEOM37_NXT MAT_NU_D2O ;
LINKED_LIST GEOM GEOMT LIB LIBSS SYS FLUX EDIT ;
SEQ_BINARY TRACK ;
REAL Keff ;
REAL TfuelK TcoolK TmodK := 900.0 343.15 343.15 ;
STRING LibName := "DLIB_J2" ;

GEOM := GEOM37_NXT :: ;
GEOMT TRACK := NXT: GEOM :: ALLG TISO 3 3.0 ;
LIB := MAT_NU_D2O :: <<TfuelK>> <<TcoolK>> <<TmodK>> <<LibName>> ;
LIBSS := USS: LIB GEOMT TRACK :: EDIT 1 GRMIN 0 PIJ ;
SYS := ASM: LIBSS GEOMT TRACK :: PIJ EDIT 1 ;
FLUX := FLU: SYS LIBSS GEOMT TRACK :: EDIT 1 TYPE K THER 5 5E-5 ;
GREP: FLUX :: GETVAL 'K-EFFECTIVE' 1 >>Keff<< ;
ECHO 'PHASE2_FRESH_K_EFFECTIVE' Keff ;
EDIT := EDI: FLUX LIBSS GEOMT ::
  EDIT 2
  MERG MIX 2 3 3 3 2 3 0 0 0 0 1 1 1 1
  COND 4.0
  SAVE
;
SYS FLUX EDIT LIBSS LIB GEOM GEOMT TRACK := DELETE: SYS FLUX EDIT LIBSS LIB GEOM GEOMT TRACK ;
ECHO "PHASE2 CANDU37 fresh nominal calculation completed" ;
END: ;
QUIT "LIST" .
"""


ACCESS = r"""#!/bin/sh
if [ $# = 0 ]; then
  echo "usage: CANDU37_NU_PHASE2.access directory" 1>&2
  exit 1
fi
MACH=`uname -s`
Sysx="`echo $MACH | cut -b -6`"
if [ $Sysx = "CYGWIN" ]; then
  MACH=`uname -o`
elif [ $Sysx = "Darwin" ]; then
  MACH=`uname -sm | sed 's/[ ]/_/'`
elif [ $Sysx = "SunOS" ]; then
  MACH=`uname -sm | sed 's/[ ]/_/'`
fi
if [ "$MACH" = "Linux" -o "$MACH" = "OSF1" -o "$MACH" = "Cygwin" -o "$MACH" = "SunOS_i86pc" \
-o "$MACH" = "Darwin_i386" -o "$MACH" = "Darwin_x86_64" -o "$MACH" = "Darwin_arm64" \
-o "$MACH" = "Linux_aarch64" ]; then
  pos=$1/../../libraries/l_endian
else
  pos=$1/../../libraries/b_endian
fi
if [ -f "$pos"/draglibJef2p2_v5p1.gz ]; then
  chmod 755 "$pos"
  gunzip "$pos"/draglibJef2p2_v5p1.gz
fi
if [ -f "$pos"/draglibJef2p2_v5p1 ]; then
  ln -s "$pos"/draglibJef2p2_v5p1 DLIB_J2
fi
ln -s "$1"/data/CANDU37_NU_PHASE2_proc/*.c2m .
echo "CANDU37_NU_PHASE2 access script terminated"
"""


def write(path: pathlib.Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def main() -> int:
    phase2 = PKG
    deck_dir = phase2 / "decks"
    proc_dir = phase2 / "proc"
    logs = phase2 / "logs"
    for d in (deck_dir, proc_dir, logs):
        d.mkdir(parents=True, exist_ok=True)

    write(deck_dir / f"{CASE}.x2m", DECK)
    write(proc_dir / "GEOM37_NXT.c2m", GEOM)
    write(proc_dir / "MAT_NU_D2O.c2m", MAT)
    write(deck_dir / f"{CASE}.access", ACCESS)

    dragon_deck = DRAGON / "data" / f"{CASE}.x2m"
    dragon_proc = DRAGON / "data" / f"{CASE}_proc"
    dragon_access = DRAGON / "data" / f"{CASE}.access"
    write(dragon_deck, DECK)
    write(dragon_proc / "GEOM37_NXT.c2m", GEOM)
    write(dragon_proc / "MAT_NU_D2O.c2m", MAT)
    write(dragon_access, ACCESS)
    dragon_access.chmod(0o755)

    proc = subprocess.run(
        ["./rdragon", f"{CASE}.x2m"],
        cwd=DRAGON,
        text=True,
        capture_output=True,
        timeout=900,
        check=False,
    )
    write(logs / "rdragon.stdout_stderr.log", proc.stdout + proc.stderr)
    result = DRAGON / "Darwin_arm64" / f"{CASE}.result"
    result_text = ""
    if result.exists():
        shutil.copy2(result, logs / f"{CASE}.result")
        result_text = result.read_text(encoding="utf-8", errors="ignore")
    success = (
        proc.returncode == 0
        and "PHASE2_FRESH_K_EFFECTIVE" in result_text
        and "PHASE2 CANDU37 fresh nominal calculation completed" in result_text
        and "READ ERROR" not in result_text
        and "NXT   :" not in result_text
    )
    k_match = re.search(r"PHASE2_FRESH_K_EFFECTIVE\s+([0-9.Ee+-]+)", result_text)
    k_eff = k_match.group(1) if k_match else "not_found"

    report = [
        "# Phase 2 DRAGON 37-Element CANDU Fresh Lattice Status",
        "",
        f"Generated UTC: {datetime.now(timezone.utc).isoformat()}",
        f"DRAGON case: `{CASE}`",
        f"rdragon exit code: `{proc.returncode}`",
        f"physics completion: `{'PASS' if success else 'FAIL'}`",
        f"fresh nominal k-effective: `{k_eff}`",
        "condensation boundary used in Phase 2 probe: `4.0 eV`",
        "",
        "## Artifacts",
        "",
        f"- Deck: `{deck_dir / (CASE + '.x2m')}`",
        f"- Access script: `{deck_dir / (CASE + '.access')}`",
        f"- Geometry procedure: `{proc_dir / 'GEOM37_NXT.c2m'}`",
        f"- Material procedure: `{proc_dir / 'MAT_NU_D2O.c2m'}`",
        f"- Run log: `{logs / 'rdragon.stdout_stderr.log'}`",
        f"- Result listing: `{logs / (CASE + '.result')}`",
        "",
        "## Status",
        "",
    ]
    if success:
        report.append("- Fresh nominal 37-element model run completed normally.")
    else:
        report.append("- Fresh nominal 37-element model did not complete normally. Review the DRAGON listing before Phase 3.")
    report.extend([
        "",
        "## Notes",
        "",
        "- This is a Phase 2 fresh lattice implementation only.",
        "- It is not yet depleted and is not a production safety library.",
        "- Geometry dimensions are benchmark-confirmed separately in `phase2_benchmark_confirmation.md`; structural alloy densities remain baseline values pending plant/material certificate confirmation.",
        "- Natural uranium and D2O coolant/moderator are specified directly in the material procedure.",
        "- The geometry is a 37-pin cluster model with rings 1, 6, 12, and 18.",
        "",
    ])
    write(phase2 / "phase2_fresh_nominal_summary.md", "\n".join(report))
    print(f"phase2={phase2}")
    print(f"exit={proc.returncode}")
    return 0 if success else 2


if __name__ == "__main__":
    raise SystemExit(main())
