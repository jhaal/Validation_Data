#!/usr/bin/env python3
"""Check that every nuclide and S(a,b) table requested by the OpenMC cases has data
covering the requested material temperatures, so the temperature tolerance is never used.

Usage: python audit_openmc_library_temperatures.py <cross_sections.xml> <cases_dir> <out.csv>
Requires h5py. Reads materials.xml in every case directory below <cases_dir>.
"""
import csv
import pathlib
import sys
import xml.etree.ElementTree as ET

import h5py


def main() -> None:
    xs, cases, out = map(pathlib.Path, sys.argv[1:4])
    paths = {m: xs.parent / e.get("path") for e in ET.parse(xs).getroot() for m in e.get("materials").split()}
    requested: dict[str, set[float]] = {}
    n_cases = 0
    for mfile in sorted(cases.glob("*/materials.xml")):
        n_cases += 1
        for m in ET.parse(mfile).getroot().findall("material"):
            t = float(m.get("temperature"))
            for el in m.findall("nuclide") + m.findall("sab"):
                requested.setdefault(el.get("name"), set()).add(t)
    rows = []
    for name, temps in sorted(requested.items()):
        with h5py.File(paths[name], "r") as f:
            avail = sorted(int(k[:-1]) for k in f[name]["kTs"].keys())
        # Interpolation needs data bracketing T; labels are rounded to whole kelvin.
        covered = all(min(avail) - 0.5 <= t <= max(avail) + 0.5 for t in temps)
        rows.append({
            "name": name,
            "requested_min_K": min(temps), "requested_max_K": max(temps), "n_requested": len(temps),
            "library_min_K": min(avail), "library_max_K": max(avail), "n_library": len(avail),
            "bracketed": covered,
        })
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"cases {n_cases}; names {len(rows)}; not bracketed {sum(not r['bracketed'] for r in rows)}")


if __name__ == "__main__":
    main()
