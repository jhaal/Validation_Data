#!/usr/bin/env python3
"""Analyse the moderator-mesh and tracking convergence study (make_mesh_probe.py).

For each VRF_MESHC_<lib>_R<n>_<trk> listing: k, whole-cell thermal-to-fast flux ratio and two-group
constants, compared with the OpenMC fresh nominal state (data/openmc_reference_states.csv,
data/openmc_two_group_rate_tallies.csv). Writes data/dragon_moderator_mesh_convergence.csv.
Usage: python analyze_mesh_probe.py <DRAGON listing dir>
"""
import csv
import pathlib
import re
import sys

DATA = pathlib.Path(__file__).resolve().parent.parent / "data"
NUM = r"([-+0-9.Ee]+)"
LINE = re.compile(r"^\s+1\s+" + r"\s+".join([NUM] * 8))
CASE = "BU00p0_TFUEL0900p00_TCOOL0343p15_TMOD0343p15"


def parse(text):
    k, g, grp, on = None, {}, None, False
    for line in text.splitlines():
        if line.startswith(">|VRF_EDI_BEGIN"):
            on = True
        elif line.startswith(">|VRF_BRANCH_B"):
            k = float(line[2:].split("|")[0].split()[5])
            on = False
        elif on:
            mg = re.search(r"G R O U P\s+:\s+(\d+)", line)
            if mg:
                grp = int(mg.group(1))
                continue
            mm = LINE.match(line)
            if mm and grp in (1, 2):
                f = [float(x) for x in mm.groups()]
                g[grp] = {"flux": f[1], "Sa": f[3] / f[1], "nuSf": f[4] / f[1]}
    return k, g


def main() -> None:
    res = pathlib.Path(sys.argv[1])
    o = next(r for r in csv.DictReader((DATA / "openmc_reference_states.csv").open()) if r["case_id"] == CASE)
    t = next(r for r in csv.DictReader((DATA / "openmc_two_group_rate_tallies.csv").open()) if r["case_id"] == CASE)
    ko = float(o["openmc_keff"])
    ref = {"ratio": float(t["thermal_to_fast_flux_ratio_from_tally"]), "Sa1": float(t["Sa1_openmc_rate"]),
           "Sa2": float(t["Sa2_openmc_rate"]), "nuSf1": float(t["nuSf1_openmc_rate"]),
           "nuSf2": float(t["nuSf2_openmc_rate"])}
    rows = []
    for f in sorted(res.glob("VRF_MESHC_*.result")):
        _, _, lib, rings, trk = f.stem.split("_")
        k, g = parse(f.read_text(errors="ignore"))
        ratio = g[2]["flux"] / g[1]["flux"]
        rows.append({"library": lib, "moderator_rings": int(rings[1:]), "tracking": trk, "k_dragon": k,
                     "delta_k_pcm": (ko - k) * 1e5, "flux_ratio_g2_g1": ratio,
                     "openmc_flux_ratio_g2_g1": ref["ratio"],
                     **{f"rel_{m}_pct": 100 * (ref[m] / g[gr][q] - 1) for m, gr, q in
                        (("Sa1", 1, "Sa"), ("Sa2", 2, "Sa"), ("nuSf1", 1, "nuSf"), ("nuSf2", 2, "nuSf"))}})
    rows.sort(key=lambda r: (r["library"], r["moderator_rings"], r["tracking"]))
    with (DATA / "dragon_moderator_mesh_convergence.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"OpenMC k {ko:.5f}, flux ratio {ref['ratio']:.3f}")
    for r in rows:
        print(f"{r['library']:4s} R{r['moderator_rings']:02d} {r['tracking']:7s} k {r['k_dragon']:.5f} dk {r['delta_k_pcm']:+6.0f} "
              f"phi2/phi1 {r['flux_ratio_g2_g1']:.3f} Sa2 {r['rel_Sa2_pct']:+.2f}% nuSf2 {r['rel_nuSf2_pct']:+.2f}% "
              f"Sa1 {r['rel_Sa1_pct']:+.2f}%")


if __name__ == "__main__":
    main()
