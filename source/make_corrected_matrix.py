#!/usr/bin/env python3
"""Full corrected DRAGON branch matrix (all 2132 production states), one deck per burnup and library.

Each deck depletes stepwise (0.5 GWd/tU, as Phase 3) to its burnup and evaluates every production
branch state at that burnup with the corrected create-mode procedure, followed by the production
whole-cell two-group condensation (EDI: MERG all mixtures, COND 4.0 eV).
Usage: python make_corrected_matrix.py <branch_state_grid.csv>
"""
import csv
import sys
import pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import make_verification_decks as m
from make_corrected_2group_and_sweeps import EDI


def deck_for(lib: str, bu: float, states: list[tuple[float, float, float]]) -> str:
    m.BRANCH_BU = [bu]
    m.BRANCHES = [("x", tf, tc, tm) for tf, tc, tm in states]
    d = m.step_deck(lib, "B", only_bu=bu)
    d = d.replace("MODULE NXT: USS: ASM: FLU: EVO: GREP: DELETE: END: ;",
                  "MODULE NXT: USS: ASM: FLU: EVO: EDI: GREP: DELETE: END: ;")
    d = d.replace("LIBBR LIBB2 SYSB FLUXB ;", "LIBBR LIBB2 SYSB FLUXB EDIT ;")
    out = []
    for line in d.splitlines():
        if line.startswith("ECHO 'VRF_BRANCH_B'"):
            out += [line.replace("VRF_BRANCH_B", "VRF_EDI_BEGIN"), EDI.rstrip("\n"), line]
        elif line == "SYSB FLUXB LIBBR LIBB2 := DELETE: SYSB FLUXB LIBBR LIBB2 ;":
            out.append("SYSB FLUXB LIBBR LIBB2 EDIT := DELETE: SYSB FLUXB LIBBR LIBB2 EDIT ;")
        else:
            out.append(line)
    return "\n".join(out) + "\n"


def main() -> None:
    by_bu: dict[float, list[tuple[float, float, float]]] = {}
    for r in csv.DictReader(open(sys.argv[1])):
        by_bu.setdefault(round(float(r["burnup_GWd_tU"]), 1), []).append(
            (float(r["tfuel_K"]), float(r["tcool_K"]), float(r["tmod_K"])))
    assert len(by_bu) == 41 and all(len(v) == 52 for v in by_bu.values())
    for lib in ("J2", "E81"):
        for bu, states in sorted(by_bu.items()):
            m.install(f"VRF_FULL_{lib}_BU{bu:04.1f}".replace(".", "p"), deck_for(lib, bu, sorted(states)), lib)




# --- Refined-moderator-mesh variant (20 equal-width moderator rings) ---------------------------
PROD_RAD = "RADIUS 0.0 0.8670 2.4545 3.7295 <<Rcool>> <<Rpt>> <<Rgap>> <<Rct>> <<Rout>>\n  MIX 1 1 1 1 2 3 4 5 5"


def fine_radius_block(n: int = 20) -> str:
    r = [f"{6.5875 + (14.2875 - 6.5875) * i / n:.4f}" for i in range(1, n)]
    lines = ["    " + " ".join(r[i:i + 6]) for i in range(0, len(r), 6)]
    mix = ["1", "1", "1", "1", "2", "3", "4"] + ["5"] * (n + 1)
    # Lists longer than 30 entries are wrapped to stay within the input line length (n <= 20 unchanged).
    mix_txt = " ".join(mix) if len(mix) <= 30 else "\n    ".join(" ".join(mix[i:i + 20]) for i in range(0, len(mix), 20))
    return ("RADIUS 0.0 0.8670 2.4545 3.7295 <<Rcool>> <<Rpt>> <<Rgap>> <<Rct>>\n" + "\n".join(lines)
            + " <<Rout>>\n  MIX " + mix_txt)


def install_fine_mesh(case: str, n: int = 20) -> None:
    g = m.DRAGON / "data" / f"{case}_proc" / "GEOM37_NXT.c2m"
    t = g.read_text()
    assert t.count(PROD_RAD) == 1 and t.count("CARCEL 8") == 1
    g.write_text(t.replace(PROD_RAD, fine_radius_block(n)).replace("CARCEL 8", f"CARCEL {7 + n}"))


def main_fine() -> None:
    by_bu: dict[float, list[tuple[float, float, float]]] = {}
    for r in csv.DictReader(open(sys.argv[2])):
        by_bu.setdefault(round(float(r["burnup_GWd_tU"]), 1), []).append(
            (float(r["tfuel_K"]), float(r["tcool_K"]), float(r["tmod_K"])))
    for lib in ("J2", "E81"):
        for bu, states in sorted(by_bu.items()):
            name = f"VRF_FULLM_{lib}_BU{bu:04.1f}".replace(".", "p")
            m.install(name, deck_for(lib, bu, sorted(states)), lib)
            install_fine_mesh(name)


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--fine-mesh":
        main_fine()
    else:
        main()
