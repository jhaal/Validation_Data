#!/usr/bin/env python3
"""Matched-inventory branch decks: deplete exactly as the production nominal history (JEF-2.2 DRAGLIB,
production geometry, 0.5-GWd/tU steps; the inventories supplied to OpenMC), then evaluate corrected
(create-mode) branches with a chosen library and moderator mesh on that same inventory.

Library switch: the branch library LIBBR is built from the target DRAGLIB (MATDPLE81 for ENDF/B-VIII.1),
and EVO: ... NOSA SET copies the depleted densities into it. Mesh switch: branches are tracked on
GEOM37F (moderator in n equal rings) while depletion uses the production GEOM37_NXT.
Usage: python make_matched_inventory_decks.py test            (10 GWd/tU nominal, 4 combinations)
       python make_matched_inventory_decks.py full <branch_state_grid.csv> <lib> <rings>
"""
import csv
import pathlib
import sys
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import make_verification_decks as m
from make_corrected_matrix import PROD_RAD, fine_radius_block
from make_corrected_2group_and_sweeps import EDI


def deck(bu: float, combos: list[tuple[str, int]], states: list[tuple[float, float, float]], edi: bool) -> str:
    m.BRANCH_BU = []
    base = m.step_deck("J2", "B", only_bu=bu)  # depletion only (no branches since BRANCH_BU is empty)
    head, tail = base.split("ECHO 'VRF_COMPLETED' ;")
    head = head.replace("MODULE NXT: USS: ASM: FLU: EVO: GREP: DELETE: END: ;",
                        "MODULE NXT: USS: ASM: FLU: EVO: EDI: GREP: DELETE: END: ;")
    head = head.replace("PROCEDURE GEOM37_NXT MATDPLD2O ;", "PROCEDURE GEOM37_NXT GEOM37F MATDPLD2O MATDPLE81 ;")
    head = head.replace("LIBBR LIBB2 SYSB FLUXB ;", "LIBBR LIBB2 SYSB FLUXB EDIT GEOMF GEOMTF ;\nSEQ_BINARY TRACKF ;")
    out = [head.rstrip("\n")]
    out += ["GEOMF := GEOM37F :: ;", "GEOMTF TRACKF := NXT: GEOMF :: ALLG TISO 3 3.0 ;"]
    t = m.days(bu)
    for lib, rings in combos:
        proc = "MATDPLD2O" if lib == "J2" else "MATDPLE81"
        g, tr = ("GEOMT", "TRACK") if rings == 1 else ("GEOMTF", "TRACKF")
        for tf, tc, tm in states:
            lab = f"{bu:.1f} {tf:.2f} {tc:.2f} {tm:.2f}"
            out += [f'LIBBR := {proc} :: {tf:.5f} {tc:.5f} {tm:.5f} "DLIB_{lib}" ;']
            if bu > 0.0:
                out += [f"BURNUP LIBBR := EVO: BURNUP LIBBR FLUX GEOMT :: EDIT 0 NOSA SET {t:.8f} DAY ;"]
            out += [f"LIBB2 := USS: LIBBR {g} {tr} :: EDIT 0 GRMIN 0 PIJ ;",
                    f"SYSB := ASM: LIBB2 {g} {tr} :: PIJ EDIT 0 ;",
                    f"FLUXB := FLU: SYSB LIBB2 {g} {tr} :: EDIT 0 TYPE K THER 5 5E-5 ;",
                    "GREP: FLUXB :: GETVAL 'K-EFFECTIVE' 1 >>Keff<< ;"]
            if edi:
                out += [f"ECHO 'VRF_EDI_BEGIN' {lab} Keff ;",
                        EDI.rstrip("\n").replace("EDI: FLUXB LIBB2 GEOMT", f"EDI: FLUXB LIBB2 {g}")]
            out += [f"ECHO 'VRF_BRANCH_B' {lab} Keff ;  ! {lib} R{rings}",
                    "SYSB FLUXB LIBBR LIBB2 := DELETE: SYSB FLUXB LIBBR LIBB2 ;" if not edi else
                    "SYSB FLUXB LIBBR LIBB2 EDIT := DELETE: SYSB FLUXB LIBBR LIBB2 EDIT ;"]
    out += ["ECHO 'VRF_COMPLETED' ;" + tail]
    return "\n".join(out)


def install(name: str, text: str, rings: int) -> None:
    m.install(name, text, "J2")
    pdir = m.DRAGON / "data" / f"{name}_proc"
    mat = (pdir / "MATDPLD2O.c2m").read_text()
    (pdir / "MATDPLE81.c2m").write_text(m.zr_isotopic(mat))
    acc = m.DRAGON / "data" / f"{name}.access"
    a = acc.read_text()
    a = a.replace('ln -sf "$1"', 'if [ -f "$pos"/draglibendfb8r1_v5p1.gz ] && [ ! -f "$pos"/draglibendfb8r1_v5p1 ]; then\n'
                  '  gunzip -k "$pos"/draglibendfb8r1_v5p1.gz\nfi\nln -s "$pos"/draglibendfb8r1_v5p1 DLIB_E81\nln -sf "$1"', 1)
    acc.write_text(a)
    g = (pdir / "GEOM37_NXT.c2m").read_text()
    n = max(rings, 20)
    assert g.count(PROD_RAD) == 1 and g.count("CARCEL 8") == 1
    (pdir / "GEOM37F.c2m").write_text(g.replace(PROD_RAD, fine_radius_block(n)).replace("CARCEL 8", f"CARCEL {7 + n}"))


def main() -> None:
    if sys.argv[1] == "test":
        nom = [(900.0, 343.15, 343.15)]
        install("VRF_MATCH_TEST_BU10p0", deck(10.0, [("J2", 1), ("J2", 20), ("E81", 1), ("E81", 20)], nom, True), 20)
        return
    _, _, status, lib, rings = sys.argv
    rings = int(rings)
    by_bu: dict[float, list[tuple[float, float, float]]] = {}
    for r in csv.DictReader(open(status)):
        by_bu.setdefault(round(float(r["burnup_GWd_tU"]), 1), []).append(
            (float(r["tfuel_K"]), float(r["tcool_K"]), float(r["tmod_K"])))
    for bu, states in sorted(by_bu.items()):
        name = f"VRF_MATCH_{lib}_R{rings:02d}_BU{bu:04.1f}".replace(".", "p")
        install(name, deck(bu, [(lib, rings)], sorted(states), True), rings)


if __name__ == "__main__":
    main()
