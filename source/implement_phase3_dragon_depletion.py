#!/usr/bin/env python3
"""Implement Phase 3 DRAGON depletion for the Phase 2 CANDU 37-element model."""

from __future__ import annotations

import csv
import math
import pathlib
import re
import shutil
import subprocess
from datetime import datetime, timezone


XS = pathlib.Path("/Users/alokjha/DRAGON/XS_FILES")
PHASE2 = XS / "phase2_dragon_37el_candu"
PKG = XS / "phase3_dragon_37el_candu_depletion"
DRAGON = pathlib.Path("/Users/alokjha/DRAGON/5.1/Dragon")
CASE = "CANDU37_NU_PHASE3_DEPL"
POWER_MW_PER_TU = 36.8
BURNUPS = [0.5 * i for i in range(41)]


MAT_DEPL = """*DECK MATDPLD2O
PARAMETER LIBRARY ::
  EDIT 1
  ::: LINKED_LIST LIBRARY ;
;
REAL TfuelK TcoolK TmodK ;
STRING LibName ;
:: >>TfuelK<< >>TcoolK<< >>TmodK<< >>LibName<< ;

MODULE LIB: ;

STRING IsoU5 IsoU8
       IsoO16 IsoD2 IsoH1 IsoZr IsoNb IsoFe4 IsoFe6 IsoFe7 IsoFe8
       IsoCr0 IsoCr2 IsoCr3 IsoCr4 IsoNi8 IsoNi0 IsoNi1 IsoNi2 IsoNi4 IsoHe4
  :=
       "U235" "U238"
       "O16" "H2_D2O" "H1_H2O" "Zr0" "Nb93" "Fe54" "Fe56" "Fe57" "Fe58"
       "Cr50" "Cr52" "Cr53" "Cr54" "Ni58" "Ni60" "Ni61" "Ni62" "Ni64" "He4" ;

REAL NU5 := 1.6100E-04 ;
REAL NU8 := 2.2480E-02 ;
REAL NO16F := 4.5280E-02 ;

LIBRARY := LIB: ::
  EDIT 1
  NMIX 14
  CTRA APOL
  PTSL
  ANIS 2
  DEPL LIB: DRAGON FIL: <<LibName>>
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
* four radial fuel zones. Depletion daughters are created by DRAGON chain data.
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


ACCESS = """#!/bin/sh
if [ $# = 0 ]; then
  echo "usage: CANDU37_NU_PHASE3_DEPL.access directory" 1>&2
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
ln -s "$1"/data/CANDU37_NU_PHASE3_DEPL_proc/*.c2m .
echo "CANDU37_NU_PHASE3_DEPL access script terminated"
"""


def burnup_to_days(bu_gwd_tu: float) -> float:
    return bu_gwd_tu * 1000.0 / POWER_MW_PER_TU


def write(path: pathlib.Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def build_deck() -> str:
    lines = [
        "*----",
        "* Phase 3 nominal DRAGON depletion for benchmark-confirmed 37-element CANDU geometry.",
        "* Burnup grid: 0.0 to 20.0 GWd/tU in 0.5 GWd/tU steps.",
        "*----",
        "MODULE NXT: USS: ASM: FLU: EVO: GREP: DELETE: END: ;",
        "PROCEDURE GEOM37_NXT MATDPLD2O ;",
        "LINKED_LIST GEOM GEOMT LIBRAW LIBSS SYS FLUX BURNUP ;",
        "SEQ_BINARY TRACK ;",
        "REAL Keff ;",
        "REAL TfuelK TcoolK TmodK := 900.0 343.15 343.15 ;",
        "STRING LibName := \"DLIB_J2\" ;",
        "",
        "GEOM := GEOM37_NXT :: ;",
        "GEOMT TRACK := NXT: GEOM :: ALLG TISO 3 3.0 ;",
        "LIBRAW := MATDPLD2O :: <<TfuelK>> <<TcoolK>> <<TmodK>> <<LibName>> ;",
        "LIBSS := USS: LIBRAW GEOMT TRACK :: EDIT 1 GRMIN 0 PIJ ;",
        "SYS := ASM: LIBSS GEOMT TRACK :: PIJ EDIT 1 ;",
        "FLUX := FLU: SYS LIBSS GEOMT TRACK :: EDIT 1 TYPE K THER 5 5E-5 ;",
        "GREP: FLUX :: GETVAL 'K-EFFECTIVE' 1 >>Keff<< ;",
        "ECHO 'PHASE3_BURNUP_STATE' 0.000000 0.000000 Keff ;",
    ]
    for i in range(1, len(BURNUPS)):
        bu0 = BURNUPS[i - 1]
        bu1 = BURNUPS[i]
        t0 = burnup_to_days(bu0)
        t1 = burnup_to_days(bu1)
        rhs = "LIBSS FLUX GEOMT" if i == 1 else "BURNUP LIBSS FLUX GEOMT"
        noex = "" if i == 1 else "      NOEX\n"
        lines.extend(
            [
                f"* Burnup step {i}: {bu0:.3f} -> {bu1:.3f} GWd/tU",
                f"BURNUP LIBSS := EVO: {rhs} ::",
                "      EDIT 1",
                f"      SAVE {t0:.8f} DAY POWR {POWER_MW_PER_TU:.8f}",
                f"{noex}      DEPL {t0:.8f} {t1:.8f} DAY POWR {POWER_MW_PER_TU:.8f}",
                f"      SET {t1:.8f} DAY ;",
                "LIBSS := USS: LIBSS LIBRAW GEOMT TRACK ;",
                "SYS := DELETE: SYS ;",
                "SYS := ASM: LIBSS GEOMT TRACK :: PIJ EDIT 0 ;",
                "FLUX := FLU: FLUX SYS LIBSS GEOMT TRACK :: EDIT 0 TYPE K THER 5 5E-5 ;",
                "GREP: FLUX :: GETVAL 'K-EFFECTIVE' 1 >>Keff<< ;",
                f"ECHO 'PHASE3_BURNUP_STATE' {bu1:.6f} {t1:.6f} Keff ;",
                f"BURNUP LIBSS := EVO: BURNUP LIBSS FLUX GEOMT :: SAVE {t1:.8f} DAY POWR {POWER_MW_PER_TU:.8f} ;",
            ]
        )
    lines.extend(
        [
            "BURNUP LIBSS LIBRAW SYS FLUX GEOM GEOMT TRACK := DELETE: BURNUP LIBSS LIBRAW SYS FLUX GEOM GEOMT TRACK ;",
            "ECHO \"PHASE3 CANDU37 depletion completed\" ;",
            "END: ;",
            "QUIT \"LIST\" .",
            "",
        ]
    )
    return "\n".join(lines)


def parse_states(result: str) -> list[dict[str, float]]:
    states = []
    pattern = re.compile(
        r"PHASE3_BURNUP_STATE\s+([0-9.Ee+-]+)\s+([0-9.Ee+-]+)\s+([0-9.Ee+-]+)"
    )
    for bu, days, keff in pattern.findall(result):
        states.append(
            {
                "burnup_GWd_tU": float(bu),
                "time_days": float(days),
                "k_eff": float(keff),
            }
        )
    return states


def write_csv(path: pathlib.Path, rows: list[dict[str, float]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["burnup_GWd_tU", "time_days", "k_eff"])
        writer.writeheader()
        for row in rows:
            writer.writerow({key: f"{value:.8e}" for key, value in row.items()})


def write_manifest(path: pathlib.Path, rows: list[dict[str, float]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["state_index", "burnup_GWd_tU", "time_days", "state_status"])
        have = {round(row["burnup_GWd_tU"], 6): row for row in rows}
        for idx, bu in enumerate(BURNUPS):
            row = have.get(round(bu, 6))
            writer.writerow(
                [
                    idx,
                    f"{bu:.3f}",
                    f"{burnup_to_days(bu):.8f}",
                    "completed" if row else "missing",
                ]
            )


def svg_plot(path: pathlib.Path, rows: list[dict[str, float]]) -> None:
    width, height = 900, 520
    margin = 70
    xs = [row["burnup_GWd_tU"] for row in rows]
    ys = [row["k_eff"] for row in rows]
    if not xs or not ys:
        write(path, "<svg xmlns='http://www.w3.org/2000/svg' width='900' height='520'></svg>\n")
        return
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    if math.isclose(ymin, ymax):
        ymin -= 0.01
        ymax += 0.01

    def px(x: float) -> float:
        return margin + (x - xmin) / (xmax - xmin or 1.0) * (width - 2 * margin)

    def py(y: float) -> float:
        return height - margin - (y - ymin) / (ymax - ymin) * (height - 2 * margin)

    points = " ".join(f"{px(x):.2f},{py(y):.2f}" for x, y in zip(xs, ys))
    circles = "\n".join(
        f"<circle cx='{px(x):.2f}' cy='{py(y):.2f}' r='3.5'><title>BU={x:.3f}, k={y:.6f}</title></circle>"
        for x, y in zip(xs, ys)
    )
    y_ticks = []
    for i in range(6):
        y = ymin + (ymax - ymin) * i / 5.0
        yy = py(y)
        y_ticks.append(
            f"<line x1='{margin-5}' y1='{yy:.2f}' x2='{width-margin}' y2='{yy:.2f}' stroke='#e6e6e6'/>"
            f"<text x='{margin-10}' y='{yy+4:.2f}' text-anchor='end' font-size='12'>{y:.5f}</text>"
        )
    x_ticks = []
    for x in [0, 5, 10, 15, 20]:
        xx = px(x)
        x_ticks.append(
            f"<line x1='{xx:.2f}' y1='{margin}' x2='{xx:.2f}' y2='{height-margin+5}' stroke='#e6e6e6'/>"
            f"<text x='{xx:.2f}' y='{height-margin+24}' text-anchor='middle' font-size='12'>{x:g}</text>"
        )
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="100%" height="100%" fill="white"/>
<text x="{width/2}" y="32" text-anchor="middle" font-family="Arial" font-size="20">Phase 3 DRAGON k-effective vs burnup</text>
<g font-family="Arial" fill="#222">
{''.join(y_ticks)}
{''.join(x_ticks)}
<line x1="{margin}" y1="{margin}" x2="{margin}" y2="{height-margin}" stroke="#222"/>
<line x1="{margin}" y1="{height-margin}" x2="{width-margin}" y2="{height-margin}" stroke="#222"/>
<polyline points="{points}" fill="none" stroke="#005f8f" stroke-width="2.4"/>
<g fill="#c43d2f">{circles}</g>
<text x="{width/2}" y="{height-20}" text-anchor="middle" font-size="14">Burnup (GWd/tU)</text>
<text x="22" y="{height/2}" text-anchor="middle" transform="rotate(-90 22 {height/2})" font-size="14">k-effective</text>
</g>
</svg>
"""
    write(path, svg)


def main() -> int:
    deck_dir = PKG / "decks"
    proc_dir = PKG / "proc"
    log_dir = PKG / "logs"
    plot_dir = PKG / "plots"
    isotopic_dir = PKG / "burnup_isotopics"
    for directory in (deck_dir, proc_dir, log_dir, plot_dir, isotopic_dir):
        directory.mkdir(parents=True, exist_ok=True)
    for stale in (
        proc_dir / "MAT_DEPL_NU_D2O.c2m",
        DRAGON / "data" / f"{CASE}_proc" / "MAT_DEPL_NU_D2O.c2m",
    ):
        stale.unlink(missing_ok=True)

    phase2_geom = PHASE2 / "proc" / "GEOM37_NXT.c2m"
    if not phase2_geom.exists():
        raise FileNotFoundError(phase2_geom)

    deck = build_deck()
    write(deck_dir / "dragon_37el_candu_depletion.x2m", deck)
    write(deck_dir / f"{CASE}.x2m", deck)
    write(deck_dir / f"{CASE}.access", ACCESS)
    (deck_dir / f"{CASE}.access").chmod(0o755)
    shutil.copy2(phase2_geom, proc_dir / "GEOM37_NXT.c2m")
    write(proc_dir / "MATDPLD2O.c2m", MAT_DEPL)

    dragon_deck = DRAGON / "data" / f"{CASE}.x2m"
    dragon_access = DRAGON / "data" / f"{CASE}.access"
    dragon_proc = DRAGON / "data" / f"{CASE}_proc"
    write(dragon_deck, deck)
    write(dragon_access, ACCESS)
    dragon_access.chmod(0o755)
    dragon_proc.mkdir(parents=True, exist_ok=True)
    shutil.copy2(phase2_geom, dragon_proc / "GEOM37_NXT.c2m")
    write(dragon_proc / "MATDPLD2O.c2m", MAT_DEPL)

    proc = subprocess.run(
        ["./rdragon", f"{CASE}.x2m"],
        cwd=DRAGON,
        text=True,
        capture_output=True,
        timeout=1800,
        check=False,
    )
    write(log_dir / f"{CASE}.stdout_stderr.log", proc.stdout + proc.stderr)
    result_src = DRAGON / "Darwin_arm64" / f"{CASE}.result"
    result_text = ""
    if result_src.exists():
        shutil.copy2(result_src, log_dir / "depletion_reference.result")
        shutil.copy2(result_src, log_dir / f"{CASE}.result")
        result_text = result_src.read_text(encoding="utf-8", errors="ignore")

    states = parse_states(result_text)
    write_csv(PKG / "depletion_reference_summary.csv", states)
    write_manifest(PKG / "burnup_state_manifest.csv", states)
    svg_plot(plot_dir / "phase3_k_effective_vs_burnup.svg", states)

    completed = {round(row["burnup_GWd_tU"], 6) for row in states}
    missing = [bu for bu in BURNUPS if round(bu, 6) not in completed]
    normal = "normal end of execution for dragon 5" in result_text
    phase_done = "PHASE3 CANDU37 depletion completed" in result_text
    monotone_nonincreasing = all(
        states[i + 1]["k_eff"] <= states[i]["k_eff"] + 5.0e-3
        for i in range(len(states) - 1)
    )
    listing_observations = []
    for pattern in (
        r"Error on total volume\s*=.*",
        r"RMS Error on total volume\s*=.*",
        r"Global RMS, maximum and average errors \(%\) on region volumes\s*:.*",
        r"Error on total surface area\s*=.*",
        r"RMS Error on total surface area\s*=.*",
        r"Global RMS, maximum and average errors \(%\) on surface areas\s*:.*",
        r"LIBFQD: UNABLE TO SET THE ELEMENTARY LETHARGY WIDTH FOR ISOTOPE.*",
    ):
        listing_observations.extend(re.findall(pattern, result_text))
    error_like_patterns = re.findall(
        r"(?im)\b(error|abort|failed|not enough|invalid)\b", result_text
    )
    success = proc.returncode == 0 and normal and phase_done and not missing and monotone_nonincreasing

    report = [
        "# Phase 3 DRAGON Depletion Status",
        "",
        f"Generated UTC: {datetime.now(timezone.utc).isoformat()}",
        f"DRAGON case: `{CASE}`",
        f"rdragon exit code: `{proc.returncode}`",
        f"depletion run: `{'PASS' if success else 'FAIL'}`",
        "",
        "## Normalization",
        "",
        f"- Specific power: `{POWER_MW_PER_TU:.8f}` MW/tU.",
        "- Time conversion: `days = burnup_GWd_tU * 1000 / specific_power_MW_per_tU`.",
        f"- 0.5 GWd/tU step length: `{burnup_to_days(0.5):.8f}` days.",
        f"- 20.0 GWd/tU end time: `{burnup_to_days(20.0):.8f}` days.",
        "",
        "## Outputs",
        "",
        f"- Depletion deck: `{deck_dir / 'dragon_37el_candu_depletion.x2m'}`",
        f"- Active material procedure: `{proc_dir / 'MATDPLD2O.c2m'}`",
        f"- Burnup state manifest: `{PKG / 'burnup_state_manifest.csv'}`",
        f"- Summary CSV: `{PKG / 'depletion_reference_summary.csv'}`",
        f"- DRAGON result: `{log_dir / 'depletion_reference.result'}`",
        f"- k-effective plot: `{plot_dir / 'phase3_k_effective_vs_burnup.svg'}`",
        "",
        "## Status Checks",
        "",
        f"- DRAGON normal end: `{'PASS' if normal else 'FAIL'}`.",
        f"- Completion marker: `{'PASS' if phase_done else 'FAIL'}`.",
        f"- Burnup states completed: `{len(states)}` of `{len(BURNUPS)}`.",
        f"- Missing burnup states: `{', '.join(f'{x:.1f}' for x in missing) if missing else 'none'}`.",
        f"- Smooth k-effective trend check: `{'PASS' if monotone_nonincreasing else 'FAIL'}`.",
        f"- k-effective range: `{min((row['k_eff'] for row in states), default=float('nan')):.8e}` to `{max((row['k_eff'] for row in states), default=float('nan')):.8e}`.",
        f"- Final burnup k-effective: `{states[-1]['k_eff']:.8e}`." if states else "- Final burnup k-effective: `not available`.",
        f"- Listing observations requiring QA review: `{len(listing_observations)}`.",
        f"- Error-like keyword hits in listing: `{len(error_like_patterns)}`. Review the observations below before using downstream.",
        "",
        "## Listing Observations",
        "",
        *(
            [f"- `{item.strip()}`" for item in listing_observations[:20]]
            if listing_observations
            else ["- `none detected`"]
        ),
        "",
        "## Limitations",
        "",
        "- This is the Phase 3 nominal depletion implementation. It is not yet the final branch library for safety use.",
        "- Isotopic inventory extraction is staged in `burnup_isotopics/`; detailed isotope table export will be added after the depletion object format is confirmed from this run.",
        "- Temperature branch cross-section generation remains Phase 4.",
        "",
    ]
    write(PKG / "phase3_depletion_status.md", "\n".join(report))

    print(PKG / "phase3_depletion_status.md")
    print(f"states={len(states)} missing={len(missing)} exit={proc.returncode}")
    print("verdict", "PASS" if success else "FAIL")
    return 0 if success else 2


if __name__ == "__main__":
    raise SystemExit(main())
