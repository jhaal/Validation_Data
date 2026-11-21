#!/usr/bin/env python3
"""Moderator-mesh and tracking convergence study (fresh fuel, nominal temperatures, corrected branch).

The production geometry (GEOM37_NXT.c2m) treats the moderator, 6.5875-14.2875 cm, as one flat-flux
region. Here it is divided into n equal-width rings (n = 1 is production) and tracked at three densities.
Each case is followed by the production whole-cell condensation (EDI: MERG all mixtures, COND 4.0 eV).
Decks VRF_MESHC_<lib>_R<n>_<trk>; run with rdragon, analyse with analyze_mesh_probe.py.
"""
import sys
import pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import make_verification_decks as m
from make_corrected_matrix import PROD_RAD, fine_radius_block

RINGS = (1, 5, 10, 20, 40)
TRACK = {"T3D3": "TISO 3 3.0", "T8D10": "TISO 8 10.0", "T12D20": "TISO 12 20.0"}


def deck(lib: str, trk: str) -> str:
    return f"""MODULE NXT: USS: ASM: FLU: EDI: GREP: END: ;
PROCEDURE GEOM37_NXT MATDPLD2O ;
LINKED_LIST GEOM GEOMT LIBBR LIBB2 SYSB FLUXB EDIT ;
SEQ_BINARY TRACK ;
REAL Keff ;
STRING LibName := "DLIB_{lib}" ;
GEOM := GEOM37_NXT :: ;
GEOMT TRACK := NXT: GEOM :: ALLG {TRACK[trk]} ;
LIBBR := MATDPLD2O :: 900.00000 343.15000 343.15000 <<LibName>> ;
LIBB2 := USS: LIBBR GEOMT TRACK :: EDIT 0 GRMIN 0 PIJ ;
SYSB := ASM: LIBB2 GEOMT TRACK :: PIJ EDIT 0 ;
FLUXB := FLU: SYSB LIBB2 GEOMT TRACK :: EDIT 0 TYPE K THER 5 5E-5 ;
GREP: FLUXB :: GETVAL 'K-EFFECTIVE' 1 >>Keff<< ;
ECHO 'VRF_EDI_BEGIN' 0.0 900.0 343.15 343.15 Keff ;
EDIT := EDI: FLUXB LIBB2 GEOMT ::
  EDIT 2
  MERG MIX 1 1 1 1 1 1 0 0 0 0 1 1 1 1
  COND 4.0
  SAVE
;
ECHO 'VRF_BRANCH_B' 0.0 900.0 343.15 343.15 Keff ;
ECHO 'VRF_COMPLETED' ;
END: ;
QUIT "LIST" .
"""


def main() -> None:
    for lib in ("J2", "E81"):
        for n in RINGS:
            for trk in TRACK:
                name = f"VRF_MESHC_{lib}_R{n:02d}_{trk}"
                m.install(name, deck(lib, trk), lib)
                if n > 1:
                    g = m.DRAGON / "data" / f"{name}_proc" / "GEOM37_NXT.c2m"
                    t = g.read_text()
                    assert t.count(PROD_RAD) == 1 and t.count("CARCEL 8") == 1
                    g.write_text(t.replace(PROD_RAD, fine_radius_block(n)).replace("CARCEL 8", f"CARCEL {7 + n}"))
                print(name)


if __name__ == "__main__":
    main()
