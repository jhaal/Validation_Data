#!/usr/bin/env python3
"""Check all operating inputs before interpreting transport results."""
import csv
import json
import math
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "source"))
from prepare_operating_coefficient_cases import BASE_DENS, BASE_MASS_KGM3, MASS, AVOGADRO, WORK, inventories


def main():
    rows = list(csv.DictReader((ROOT / "data/operating_coefficient_states.csv").open(newline="")))
    assert len(rows) == 63 and len({r["seed"] for r in rows}) == 63
    assert len({r["case_id"] for r in rows}) == 63
    available = list(csv.DictReader((ROOT / "data/operating_library_temperature_audit.csv").open(newline="")))
    assert len(available) == 260 and all(r["bracketed"] == "True" for r in available)
    for row in rows:
        case = WORK / "openmc_cases" / row["case_id"]
        settings = ET.parse(case / "settings.xml").getroot()
        assert settings.findtext("seed") == row["seed"]
        assert settings.findtext("temperature_method") == "interpolation"
        assert settings.find("temperature_range") is None
        mats = ET.parse(case / "materials.xml").getroot().findall("material")
        fluids = 0
        for mat in mats:
            name, temp = mat.get("name"), float(mat.get("temperature"))
            dens = float(mat.find("density").get("value"))
            assert mat.find("density").get("units") == "atom/b-cm"
            fractions = {n.get("name"): float(n.get("ao")) for n in mat.findall("nuclide")}
            assert math.isclose(sum(fractions.values()), 1., abs_tol=1e-12)
            if name.startswith("fuel_zone_"):
                zone = int(name.split("_")[2]); target = inventories()[float(row["burnup_GWd_tU"]), zone]
                assert set(fractions) == set(target)
                for nuclide in target:
                    assert math.isclose(dens*fractions[nuclide], target[nuclide], rel_tol=1e-12)
                assert temp == float(row["tfuel_K"])
            elif name in ("coolant_D2O", "moderator_D2O"):
                fluid = "coolant" if name == "coolant_D2O" else "moderator"
                fluids += 1
                expected_t = float(row["tcool_K"] if fluid == "coolant" else row["tmod_K"])
                assert temp == expected_t
                actual_mass = dens*sum(fractions[n]*MASS[n] for n in fractions)*1e27/AVOGADRO
                assert math.isclose(actual_mass, float(row[fluid+"_density_kg_m3"]), rel_tol=1e-12)
                for n in fractions:
                    assert math.isclose(dens*fractions[n], BASE_DENS[n]*float(row[fluid+"_scale"]), rel_tol=1e-12)
            else:
                assert temp == (343.15 if name.startswith("calandria") else 561.15)
        assert fluids == 2
    print("PASS 63 unique inputs/seeds, exact fuel inventories, fluid mass densities, fixed structural temperatures, and temperature coverage")


if __name__ == "__main__": main()
