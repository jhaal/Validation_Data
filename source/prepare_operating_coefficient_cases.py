#!/usr/bin/env python3
"""Prepare an illustrative operating-state lattice coefficient study.

User-selected reference: Tf=900 K, Tc=561.15 K at 10.1 MPa,
Tm=343.15 K at 0.101325 MPa. Pure-heavy-water IAPWS density is used
as a surrogate while preserving the benchmark's trace-H isotope ratios.
Structural temperatures, densities and geometry stay fixed at the reference.
Fuel inventories stay those of the existing benchmark (not operating depletion).

Run with an interpreter containing iapws. Generates a state manifest and OpenMC
XML files using the original models; no OpenMC Python API is required.
--dragon additionally installs new, uniquely named DRAGON decks/procedures.
"""
import argparse
import copy
import csv
import hashlib
import json
import math
import os
import re
import shutil
import sys
import xml.etree.ElementTree as ET
from functools import lru_cache
from pathlib import Path

from iapws import D2O
import iapws

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "non_upload/operating_coefficients"
BASE = Path(os.environ.get("OPERATING_ORIGINAL_CASES", "/Volumes/Alok_Jha_Orico/DRAGON/XS_FILES/phaseF_openmc_fullgrid_purr_benchmark/full_matrix_openmc_p40k/openmc_cases"))
XS = Path(os.environ.get("OPENMC_CANDU_OUT_DIR", "/Volumes/Alok_Jha_Orico/OpenMC/endfb-viii.1-hdf5-candu-fullgrid-purr"))
TF, TC, TM, PC, PM = 900., 561.15, 343.15, 10.1, .101325
BASE_DENS = {"H2": .064886, "H1": .000064950, "O16": .032398}
MASS = {"H2": 2.014101777844, "H1": 1.007825031898, "O16": 15.99491461926}
AVOGADRO = 6.02214076e23
BASE_MASS_KGM3 = sum(BASE_DENS[n]*MASS[n] for n in BASE_DENS)*1e27/AVOGADRO


def write_csv(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def eos_checks():
    # IAPWS R16-17(2018), Table 7: P at (T, molar density).
    tests = [(300., 55.126, .0529123711), (500., 46.26, 8.35329492),
             (500., .05, .206052588), (800., .25, 1.64466177)]
    out = []
    for t, molar, expected in tests:
        actual = D2O(T=t, rho=molar*D2O.M).P
        assert math.isclose(actual, expected, rel_tol=2e-7, abs_tol=2e-8), (t, actual, expected)
        out.append({"T_K": t, "density_mol_dm3": molar, "reference_P_MPa": expected,
                    "computed_P_MPa": actual})
    return out


def density(t, p):
    state = D2O(T=t, P=p)
    assert state.status == 1 and state.phase == "Liquid", (t, p, state.phase, state.msg)
    return float(state.rho)


@lru_cache(maxsize=1)
def inventories():
    result = {}
    for row in csv.DictReader((ROOT / "data/fuel_isotopic_inventories.csv").open(newline="")):
        den = float(row["atom_density_barn_cm"])
        if den <= 0: continue
        name = row["isotope"].split()[0]
        name = {"U5": "U235", "U8": "U238", "D2": "H2"}.get(name, name)
        if name.endswith("m") and name[-2].isdigit(): name = name[:-1]+"_m1"
        key = (round(float(row["burnup_GWd_tU"]), 8), int(row["fuel_zone"]))
        zone = result.setdefault(key, {})
        zone[name] = zone.get(name, 0.) + den
    return result


def make_states(burnups):
    rc, rm = density(TC, PC), density(TM, PM)
    states = []
    for bu in burnups:
        definitions = [("nominal", "none", "reference", 0., 0.)]
        for variable, steps, modes in (("fuel", (50., 100.), ("fixed",)),
                                       ("coolant", (5., 10.), ("fixed", "pressure")),
                                       ("moderator", (5., 10.), ("fixed", "pressure"))):
            for mode in modes:
                for h in steps:
                    for sign in (-1., 1.):
                        definitions.append((f"{variable}_{mode}_{'p' if sign>0 else 'm'}{h:g}", variable, mode, h, sign))
        for tag, var, mode, h, sign in definitions:
            tf, tc, tm = TF + (sign*h if var == "fuel" else 0), TC + (sign*h if var == "coolant" else 0), TM + (sign*h if var == "moderator" else 0)
            dc = density(tc, PC) if var == "coolant" and mode == "pressure" else rc
            dm = density(tm, PM) if var == "moderator" and mode == "pressure" else rm
            case = f"OP_BU{bu:04.1f}_{tag}".replace(".", "p")
            states.append({"case_id": case, "burnup_GWd_tU": bu, "variable": var, "density_mode": mode,
                           "step_K": h, "sign": sign, "tfuel_K": tf, "tcool_K": tc, "tmod_K": tm,
                           "coolant_reference_pressure_MPa": PC, "moderator_reference_pressure_MPa": PM,
                           "coolant_density_kg_m3": dc, "moderator_density_kg_m3": dm,
                           "coolant_scale": dc / BASE_MASS_KGM3, "moderator_scale": dm / BASE_MASS_KGM3,
                           "seed": 2026112001 + 2*len(states),
                           "structural_temperatures": "fixed_at_reference", "geometry": "fixed"})
    return states


def case_xml(row, particles, batches, inactive):
    bu = row["burnup_GWd_tU"]
    # Reuse the common fresh geometry/IDs; populate every zone explicitly from
    # the deposited inventory table, including external-worker burnups.
    name = "BU00p0_TFUEL0900p00_TCOOL0343p15_TMOD0343p15"
    src = BASE / name
    if not src.is_dir():
        raise FileNotFoundError(src)
    dst = WORK / "openmc_cases" / row["case_id"]
    dst.mkdir(parents=True, exist_ok=True)
    if list(dst.glob("statepoint.*.h5")):
        raise FileExistsError(f"Refusing to overwrite completed case {dst}")
    tree = ET.parse(src / "materials.xml"); mats = tree.getroot()
    mats.find("cross_sections").text = str(XS / "cross_sections.xml")
    for mat in mats.findall("material"):
        name = mat.get("name")
        if name.startswith("fuel_zone_"):
            mat.set("temperature", str(row["tfuel_K"]))
            zone = int(name.split("_")[2])
            composition = inventories()[(bu, zone)]
            total = sum(composition.values())
            mat.find("density").set("value", str(total))
            mat.set("name", f"fuel_zone_{zone}_BU{bu:.3f}")
            for nuc in mat.findall("nuclide"): mat.remove(nuc)
            for nuc, den in composition.items():
                ET.SubElement(mat, "nuclide", {"name": nuc, "ao": str(den/total)})
        elif name in ("coolant_D2O", "moderator_D2O_outer_square", "moderator_D2O"):
            variable = "coolant" if name == "coolant_D2O" else "moderator"
            mat.set("temperature", str(row["tcool_K"] if variable == "coolant" else row["tmod_K"]))
            den = mat.find("density")
            den.set("value", str(sum(BASE_DENS.values())*row[f"{variable}_scale"]))
            assert den.get("units") == "atom/b-cm"
        else:
            mat.set("temperature", str(TM if name.startswith("calandria") else TC))
    tree.write(dst / "materials.xml", encoding="utf-8", xml_declaration=True)
    # This extension evaluates eigenvalue derivatives only. Whole-cell rate
    # tallies are unnecessary and create severe shared-counter contention.
    for file in ("geometry.xml",):
        shutil.copy2(src / file, dst / file)
    if (dst / "tallies.xml").exists():
        (dst / "tallies.xml").rename(dst / "tallies_omitted.xml")
    tree = ET.parse(src / "settings.xml"); settings = tree.getroot()
    for tag, value in (("particles", particles), ("batches", batches), ("inactive", inactive),
                        ("seed", row["seed"]), ("temperature_tolerance", .5)):
        node = settings.find(tag)
        if node is None: node = ET.SubElement(settings, tag)
        node.text = str(value)
    # Load only data bracketing actual material temperatures, not the whole grid.
    node = settings.find("temperature_range")
    if node is not None: settings.remove(node)
    ET.SubElement(settings, "entropy_mesh").text = "1"
    mesh = ET.SubElement(settings, "mesh", {"id": "1"})
    ET.SubElement(mesh, "dimension").text = "8 8 1"
    ET.SubElement(mesh, "lower_left").text = "-14.2875 -14.2875 -10000.0"
    ET.SubElement(mesh, "upper_right").text = "14.2875 14.2875 10000.0"
    tree.write(dst / "settings.xml", encoding="utf-8", xml_declaration=True)


def install_dragon(states, rings=20, prefix="OP_CONV_J2",
                   solver="EXTE 200 1E-6 UNKT 5E-5 THER 100 1E-5"):
    import make_verification_decks as m
    import make_matched_inventory_decks as match
    from make_corrected_matrix import PROD_RAD, fine_radius_block
    from make_corrected_2group_and_sweeps import EDI
    grouped = {}
    for r in states: grouped.setdefault(r["burnup_GWd_tU"], []).append(r)
    for bu, rows in grouped.items():
        m.BRANCH_BU = []
        base = m.step_deck("J2", "B", only_bu=bu)
        head, tail = base.split("ECHO 'VRF_COMPLETED' ;")
        head = head.replace("MODULE NXT:", "MODULE LIB: NXT:")
        head = head.replace("PROCEDURE GEOM37_NXT MATDPLD2O ;", "PROCEDURE GEOM37_NXT GEOM37F MATDPLD2O MATOPER ;")
        head = head.replace("LIBBR LIBB2 SYSB FLUXB ;", "LIBBR LIBB2 SYSB FLUXB GEOMF GEOMTF EDIT ;\nSEQ_BINARY TRACKF ;")
        head = head.replace("REAL Keff ;", "REAL Keff CheckD CheckM ;")
        head = head.replace("MODULE LIB: NXT: USS: ASM: FLU: EVO: GREP:", "MODULE LIB: NXT: USS: ASM: FLU: EVO: EDI: GREP:")
        lines = [head, "GEOMF := GEOM37F :: ;", "GEOMTF TRACKF := NXT: GEOMF :: ALLG TISO 8 10.0 ;"]
        for r in rows:
            tf, tc, tm = r["tfuel_K"], r["tcool_K"], r["tmod_K"]
            lines.append(f'LIBBR := MATOPER :: {tf:.8f} {tc:.8f} {tm:.8f} "DLIB_J2" ;')
            if bu > 0:
                lines.append(f"BURNUP LIBBR := EVO: BURNUP LIBBR FLUX GEOMT :: EDIT 0 NOSA SET {m.days(bu):.8f} DAY ;")
            # EVO SET restores saved nonfuel densities too. Reapply fluid densities
            # AFTER inventory transfer, before creating the shielded library.
            lines.append('LIBBR := LIB: LIBBR :: EDIT 0 MIXS LIB: DRAGON FIL: "DLIB_J2"')
            for mix, temp, scale in ((1, tc, r["coolant_scale"]), (5, tm, r["moderator_scale"])):
                lines.append(f" MIX {mix} {temp:.8f}")
                for alias, isotope, den in (("D2", "H2_D2O", BASE_DENS["H2"]),
                                            ("H1", "H1_H2O", BASE_DENS["H1"]),
                                            ("O16", "O16", BASE_DENS["O16"])):
                    lines.append(f'  {alias} = "{isotope}" {den*scale:.12E}')
            lines.append(" ;")
            lab = f"{bu:.1f} {tf:.2f} {tc:.2f} {tm:.2f}"
            lines += ["GREP: LIBBR :: GETVAL 'ISOTOPESDENS' 1 >>CheckD<< ;",
                      "GREP: LIBBR :: GETVAL 'ISOTOPESDENS' 21 >>CheckM<< ;",
                      f"ECHO 'OP_DENSITY' '{r['case_id']}' CheckD CheckM ;",
                      "LIBB2 := USS: LIBBR GEOMTF TRACKF :: EDIT 0 GRMIN 0 PIJ ;",
                      "SYSB := ASM: LIBB2 GEOMTF TRACKF :: PIJ EDIT 0 ;",
                      "FLUXB := FLU: SYSB LIBB2 GEOMTF TRACKF :: EDIT 1 TYPE K",
                      f"  {solver} ;",
                      "GREP: FLUXB :: GETVAL 'K-EFFECTIVE' 1 >>Keff<< ;",
                      f"ECHO 'OP_RESULT' '{r['case_id']}' Keff ;",
                      f"ECHO 'VRF_EDI_BEGIN' {lab} Keff ;",
                      EDI.rstrip().replace("EDI: FLUXB LIBB2 GEOMT", "EDI: FLUXB LIBB2 GEOMTF"),
                      f"ECHO 'VRF_BRANCH_B' {lab} Keff ;",
                      "SYSB FLUXB LIBBR LIBB2 EDIT := DELETE: SYSB FLUXB LIBBR LIBB2 EDIT ;"]
        lines += ["ECHO 'VRF_COMPLETED' ;" + tail]
        name = f"{prefix}_BU{bu:04.1f}".replace(".", "p")
        m.install(name, "\n".join(lines), "J2")
        proc = m.DRAGON / "data" / f"{name}_proc"
        geometry = (proc / "GEOM37_NXT.c2m").read_text()
        assert geometry.count(PROD_RAD) == 1
        (proc / "GEOM37F.c2m").write_text(geometry.replace(PROD_RAD, fine_radius_block(rings)).replace("CARCEL 8", f"CARCEL {7+rings}"))
        mat = (proc / "MATDPLD2O.c2m").read_text().replace("*DECK MATDPLD2O", "*DECK MATOPER")
        for mix in (2, 3, 6): mat = mat.replace(f"MIX {mix} <<TcoolK>>", f"MIX {mix} {TC:.8f}")
        mat = mat.replace("MIX 4 <<TmodK>>", f"MIX 4 {TM:.8f}")
        (proc / "MATOPER.c2m").write_text(mat)
        target = WORK / "dragon_inputs" / name
        target.mkdir(parents=True, exist_ok=True)
        for file in proc.glob("*.c2m"): shutil.copy2(file, target / file.name)
        shutil.copy2(m.DRAGON / "data" / f"{name}.x2m", target / f"{name}.x2m")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--burnups", nargs="+", type=float, default=[0., 7.5, 20.])
    parser.add_argument("--particles", type=int, default=40000)
    parser.add_argument("--batches", type=int, default=200)
    parser.add_argument("--inactive", type=int, default=50)
    parser.add_argument("--dragon", action="store_true")
    args = parser.parse_args()
    WORK.mkdir(parents=True, exist_ok=True)
    checks = eos_checks()
    rows = make_states(args.burnups)
    for row in rows: case_xml(row, args.particles, args.batches, args.inactive)
    write_csv(ROOT / "data/operating_coefficient_states.csv", rows)
    manifest = {"reference": {"Tf_K": TF, "Tc_K": TC, "Tm_K": TM, "Pc_MPa": PC, "Pm_MPa": PM},
                "illustrative_reference_authorized_by_user": True, "iapws_version": iapws.__version__,
                "EOS": "IAPWS R16-17(2018); pure D2O density surrogate at fixed benchmark isotope ratios",
                "EOS_source": "https://iapws.org/technical-guidance/release/Heavy",
                "EOS_verification": checks,
                "base_composition_mass_density_kg_m3": BASE_MASS_KGM3,
                "structural_temperatures": "Held at Tc0 for sheath, pressure tube and helium; Tm0 for calandria",
                "structural_densities_and_geometry": "Fixed; no thermal expansion",
                "inventories": "Existing shared benchmark histories; not independently validated operating depletion",
                "scattering_correction": "Existing uncorrected resonance-scattering treatment retained",
                "cases": len(rows), "particles": args.particles, "batches": args.batches, "inactive": args.inactive,
                "active_histories_per_case": args.particles*(args.batches-args.inactive),
                "seeds": "Distinct odd seeds for every state"}
    (ROOT / "qa/operating_coefficient_manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    if args.dragon: install_dragon(rows)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__": main()
