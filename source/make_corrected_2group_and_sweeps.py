#!/usr/bin/env python3
"""Corrected-branch decks for (a) whole-cell two-group constants and (b) full coolant/moderator sweeps.

VRF_2G_<lib>:    corrected (create-mode) branches at fuel 293.15/900/1993.15 K (nominal coolant and
                 moderator) at 0, 5, 10, 15, 20 GWd/tU, each followed by the production whole-cell
                 EDI: condensation (MERG all mixtures to 1, COND 4.0 eV), as in Phase 6.
VRF_SWEEP_<lib>: corrected branches at all 8 coolant and all 8 moderator temperatures (other
                 temperatures nominal) at 0, 5, 10, 15, 20 GWd/tU.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import make_verification_decks as m

BU = [0.0, 5.0, 10.0, 15.0, 20.0]
SWEEP_T = [293.15, 303.15, 313.15, 323.15, 333.15, 343.15, 353.15, 363.15]
EDI = ("EDIT := EDI: FLUXB LIBB2 GEOMT ::\n  EDIT 2\n  MERG MIX 1 1 1 1 1 1 0 0 0 0 1 1 1 1\n"
       "  COND 4.0\n  SAVE\n;\n")


def two_group_deck(lib: str) -> str:
    m.BRANCH_BU = BU
    m.BRANCHES = [("fuel", t, 343.15, 343.15) for t in (293.15, 900.0, 1993.15)]
    d = m.step_deck(lib, "B")
    d = d.replace("MODULE NXT: USS: ASM: FLU: EVO: GREP: DELETE: END: ;",
                  "MODULE NXT: USS: ASM: FLU: EVO: EDI: GREP: DELETE: END: ;")
    d = d.replace("LIBBR LIBB2 SYSB FLUXB ;", "LIBBR LIBB2 SYSB FLUXB EDIT ;")
    out = []
    for line in d.splitlines():
        if line.startswith("ECHO 'VRF_BRANCH_B'"):
            out.append(line.replace("VRF_BRANCH_B", "VRF_EDI_BEGIN"))
            out.append(EDI.rstrip("\n"))
            out.append(line)
        elif line == "SYSB FLUXB LIBBR LIBB2 := DELETE: SYSB FLUXB LIBBR LIBB2 ;":
            out.append("SYSB FLUXB LIBBR LIBB2 EDIT := DELETE: SYSB FLUXB LIBBR LIBB2 EDIT ;")
        else:
            out.append(line)
    return "\n".join(out) + "\n"


def sweep_deck(lib: str) -> str:
    m.BRANCH_BU = BU
    m.BRANCHES = ([("cool", 900.0, t, 343.15) for t in SWEEP_T] +
                  [("mod", 900.0, 343.15, t) for t in SWEEP_T if t != 343.15])
    return m.step_deck(lib, "B")


if __name__ == "__main__":
    for lib in ("J2", "E81"):
        m.install(f"VRF_2G_{lib}", two_group_deck(lib), lib)
        m.install(f"VRF_SWEEP_{lib}", sweep_deck(lib), lib)
