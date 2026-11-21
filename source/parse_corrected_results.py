#!/usr/bin/env python3
"""Parse corrected-branch DRAGON listings: eigenvalues (VRF_BRANCH_B) and whole-cell two-group
constants from the EDI: printout that follows each VRF_EDI_BEGIN marker.

Usage: python parse_corrected_results.py <glob-prefix> <out_states.csv>
e.g.   python parse_corrected_results.py VRF_FULL_J2 corrected_states_J2.csv
"""
import csv
import pathlib
import re
import sys

RES = pathlib.Path("/Users/alokjha/DRAGON/5.1/Dragon/Darwin_arm64")
NUM = r"([-+0-9.Ee]+)"
LINE = re.compile(r"^\s+1\s+" + r"\s+".join([NUM] * 8))


def case_id(bu, tf, tc, tm):
    t = lambda v: f"{v:07.2f}".replace(".", "p")
    return f"BU{bu:04.1f}".replace(".", "p") + f"_TFUEL{t(tf)}_TCOOL{t(tc)}_TMOD{t(tm)}"


def parse(text: str):
    out = []
    cur = None
    group = None
    for line in text.splitlines():
        if line.startswith(">|VRF_EDI_BEGIN"):
            v = [float(x) for x in line[2:].split("|")[0].split()[1:5]]
            cur = {"key": tuple(v), "g": {}}
            continue
        if line.startswith(">|VRF_BRANCH_B"):
            tok = line[2:].split("|")[0].split()
            key = tuple(float(x) for x in tok[1:5])
            rec = {"key": key, "k": float(tok[5]), "g": cur["g"] if cur and cur["key"] == key else {}}
            out.append(rec)
            cur = None
            continue
        if cur is None:
            continue
        g = re.search(r"G R O U P\s+:\s+(\d+)", line)
        if g:
            group = int(g.group(1))
            continue
        mm = LINE.match(line)
        if mm and group in (1, 2):
            f = [float(x) for x in mm.groups()]
            flux = f[1]
            cur["g"][group] = {"Sa": f[3] / flux, "nuSf": f[4] / flux, "S_out": f[7] / flux}
    return out


def main() -> None:
    prefix, out = sys.argv[1], pathlib.Path(sys.argv[2])
    rows = []
    for f in sorted(RES.glob(f"{prefix}*.result")):
        txt = f.read_text(errors="ignore")
        if ">|VRF_COMPLETED" not in txt:
            print("INCOMPLETE:", f.name)
        for r in parse(txt):
            bu, tf, tc, tm = r["key"]
            g = r["g"]
            rows.append({"case_id": case_id(bu, tf, tc, tm), "burnup_GWd_tU": bu, "tfuel_K": tf, "tcool_K": tc,
                         "tmod_K": tm, "k_corrected": r["k"],
                         "Sa1": g.get(1, {}).get("Sa", ""), "Sa2": g.get(2, {}).get("Sa", ""),
                         "S12": g.get(1, {}).get("S_out", ""),
                         "nuSf1": g.get(1, {}).get("nuSf", ""), "nuSf2": g.get(2, {}).get("nuSf", "")})
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"{out}: {len(rows)} states, {sum(1 for r in rows if r['Sa2'] != '')} with two-group constants")


if __name__ == "__main__":
    main()
