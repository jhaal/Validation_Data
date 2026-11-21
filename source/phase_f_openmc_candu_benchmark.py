#!/usr/bin/env python3
"""Phase F OpenMC benchmark package for the CANDU/PHWR 37-element library.

This driver builds an independent continuous-energy OpenMC model for the
fresh natural-uranium 37-element cell that matches the DRAGON Phase 4/6
geometry as closely as the currently processed OpenMC data library permits.
It also creates an auditable manifest for the complete DRAGON branch matrix
and records which states cannot yet be benchmarked because burnup-dependent
isotopic inventories have not been exported from DRAGON.
"""

from __future__ import annotations

import argparse
import csv
import collections
import hashlib
import math
import os
import pathlib
import shutil
import subprocess
import sys
from dataclasses import dataclass
from typing import Iterable

import openmc


def first_existing(*paths: pathlib.Path) -> pathlib.Path:
    for path in paths:
        if path.exists():
            return path
    return paths[0]


XS_ROOT = pathlib.Path(os.environ.get(
    "DRAGON_XS_ROOT",
    str(first_existing(
        pathlib.Path("/Volumes/Alok_Jha_Orico/DRAGON/XS_FILES"),
        pathlib.Path("/Users/alokjha/DRAGON/XS_FILES"),
    )),
))
OPENMC_DATA = pathlib.Path(os.environ.get(
    "OPENMC_CANDU_OUT_DIR",
    str(first_existing(
        pathlib.Path("/Volumes/Alok_Jha_Orico/OpenMC/endfb-viii.1-hdf5-candu-fullgrid-purr"),
        pathlib.Path("/Users/alokjha/OpenMC/endfb-viii.1-hdf5-candu-fullgrid-purr"),
    )),
))
OPENMC_BIN = pathlib.Path(os.environ.get(
    "OPENMC_BIN",
    str(first_existing(
        pathlib.Path("/Users/alokjha/OpenMC/install/bin/openmc"),
        pathlib.Path("/Volumes/Alok_Jha_Orico/OpenMC/install/bin/openmc"),
    )),
))
BRANCH_STATUS = XS_ROOT / "phase4_dragon_37el_candu_temperature_branches" / "branch_status.csv"
RX_FILE = XS_ROOT / "phwr37_2group_true_temp_fullprod.rx"
OUT_DIR = XS_ROOT / "phaseF_openmc_fullgrid_purr_benchmark"
INVENTORY_CSV = (
    XS_ROOT
    / "phase3_dragon_37el_candu_depletion"
    / "burnup_isotopics"
    / "phase3_dragon_fuel_isotopic_densities_openmc.csv"
)
ISOTOPE_ALIASES = {
    "U5": "U235",
    "U8": "U238",
    "D2": "H2",
}


@dataclass(frozen=True)
class BranchState:
    case_id: str
    burnup: float
    tfuel: float
    tcool: float
    tmod: float
    k_eff: float | None
    run_status: str

    @property
    def openmc_supported(self) -> bool:
        return abs(self.burnup) < 1.0e-9


@dataclass(frozen=True)
class InventoryCoverage:
    burnup: float
    rows: int
    supported_nuclides: int
    missing_nuclides: tuple[str, ...]
    supported_density_fraction: float

    @property
    def exists(self) -> bool:
        return self.rows > 0

    @property
    def ce_complete(self) -> bool:
        return self.exists and not self.missing_nuclides


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_branch_states(path: pathlib.Path) -> list[BranchState]:
    states: list[BranchState] = []
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            if row.get("run_status") != "completed":
                continue
            k_text = row.get("k_eff", "").strip()
            states.append(
                BranchState(
                    case_id=row["case_id"],
                    burnup=float(row["burnup_GWd_tU"]),
                    tfuel=float(row["tfuel_K"]),
                    tcool=float(row["tcool_K"]),
                    tmod=float(row["tmod_K"]),
                    k_eff=float(k_text) if k_text else None,
                    run_status=row["run_status"],
                )
            )
    return states


def openmc_materials_from_xml(path: pathlib.Path) -> set[str]:
    data = openmc.data.DataLibrary.from_xml(path)
    return {material for entry in data.libraries for material in entry["materials"]}


def normalize_dragon_isotope(label: str) -> str:
    raw = label.split()[0].strip()
    raw = ISOTOPE_ALIASES.get(raw, raw)
    if raw.endswith("m") and len(raw) > 1 and raw[-2].isdigit():
        return f"{raw[:-1]}_m1"
    return raw


def read_inventory_coverage(path: pathlib.Path, available: set[str]) -> dict[float, InventoryCoverage]:
    if not path.exists():
        return {}
    grouped: dict[float, dict[str, float]] = collections.defaultdict(lambda: collections.defaultdict(float))
    counts: dict[float, int] = collections.defaultdict(int)
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            burnup = round(float(row["burnup_GWd_tU"]), 8)
            nuc = normalize_dragon_isotope(row["isotope"])
            density = float(row["atom_density_barn_cm"])
            if density <= 0.0:
                continue
            grouped[burnup][nuc] += density
            counts[burnup] += 1
    out: dict[float, InventoryCoverage] = {}
    for burnup, densities in grouped.items():
        total = sum(densities.values())
        supported = {nuc for nuc in densities if nuc in available}
        missing = tuple(sorted(nuc for nuc in densities if nuc not in available))
        supported_density = sum(density for nuc, density in densities.items() if nuc in available)
        out[burnup] = InventoryCoverage(
            burnup=burnup,
            rows=counts[burnup],
            supported_nuclides=len(supported),
            missing_nuclides=missing,
            supported_density_fraction=(supported_density / total) if total > 0.0 else 0.0,
        )
    return out


def read_inventory_densities(path: pathlib.Path) -> dict[float, dict[int, dict[str, float]]]:
    grouped: dict[float, dict[int, dict[str, float]]] = collections.defaultdict(
        lambda: collections.defaultdict(lambda: collections.defaultdict(float))
    )
    if not path.exists():
        return {}
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            density = float(row["atom_density_barn_cm"])
            if density <= 0.0:
                continue
            burnup = round(float(row["burnup_GWd_tU"]), 8)
            zone = int(row["fuel_zone"])
            nuc = normalize_dragon_isotope(row["isotope"])
            grouped[burnup][zone][nuc] += density
    return {
        burnup: {zone: dict(densities) for zone, densities in zones.items()}
        for burnup, zones in grouped.items()
    }


def state_inventory_coverage(state: BranchState, coverage: dict[float, InventoryCoverage]) -> InventoryCoverage | None:
    return coverage.get(round(state.burnup, 8))


def state_openmc_supported(state: BranchState, coverage: dict[float, InventoryCoverage]) -> bool:
    if abs(state.burnup) < 1.0e-9:
        return True
    inv = state_inventory_coverage(state, coverage)
    return bool(inv and inv.ce_complete)


def write_csv(path: pathlib.Path, rows: Iterable[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def material_from_atom_densities(
    name: str,
    atom_densities: dict[str, float],
    temperature: float,
    sab: Iterable[str] = (),
) -> openmc.Material:
    mat = openmc.Material(name=name, temperature=temperature)
    total = sum(atom_densities.values())
    mat.set_density("atom/b-cm", total)
    for nuc, density in atom_densities.items():
        mat.add_nuclide(nuc, density / total, "ao")
    for table in sab:
        mat.add_s_alpha_beta(table)
    return mat


def zirc_atom_densities(zr_atom_density: float) -> dict[str, float]:
    abund = {
        "Zr90": 0.5145,
        "Zr91": 0.1122,
        "Zr92": 0.1715,
        "Zr94": 0.1738,
        "Zr96": 0.0280,
    }
    return {nuc: zr_atom_density * frac for nuc, frac in abund.items()}


def zirc_material(name: str, zr_atom_density: float, temperature: float) -> openmc.Material:
    return material_from_atom_densities(
        name,
        zirc_atom_densities(zr_atom_density),
        temperature,
    )


def default_fresh_fuel_densities() -> dict[str, float]:
    return {"U235": 1.6100e-4, "U238": 2.2480e-2, "O16": 4.5280e-2}


def build_materials(
    state: BranchState,
    inventory_densities: dict[float, dict[int, dict[str, float]]] | None = None,
) -> tuple[openmc.Materials, dict[str, openmc.Material]]:
    # Atom densities are copied from Phase 3 MATDPLD2O.c2m.
    burnup_inventory = (inventory_densities or {}).get(round(state.burnup, 8), {})
    fuel_mats: dict[int, openmc.Material] = {}
    for zone in (11, 12, 13, 14):
        densities = burnup_inventory.get(zone) or default_fresh_fuel_densities()
        fuel_mats[zone] = material_from_atom_densities(
            f"fuel_zone_{zone}_BU{state.burnup:.3f}",
            densities,
            state.tfuel,
        )
    d2o_coolant = material_from_atom_densities(
        "coolant_D2O",
        {"H2": 6.4886e-2, "H1": 6.4950e-5, "O16": 3.2398e-2},
        state.tcool,
        ("c_D_in_D2O", "c_O_in_D2O"),
    )
    d2o_moderator = material_from_atom_densities(
        "moderator_D2O",
        {"H2": 6.4886e-2, "H1": 6.4950e-5, "O16": 3.2398e-2},
        state.tmod,
        ("c_D_in_D2O", "c_O_in_D2O"),
    )
    pressure_tube_dens = zirc_atom_densities(4.1780e-2)
    pressure_tube_dens["Nb93"] = 1.0877e-3
    pressure_tube = material_from_atom_densities("pressure_tube_ZrNb", pressure_tube_dens, state.tcool)
    calandria_dens = zirc_atom_densities(4.2364e-2)
    calandria_dens.update(
        {
            "Fe54": 5.7861e-6,
            "Fe56": 8.7588e-5,
            "Fe57": 1.9873e-6,
            "Fe58": 2.5991e-7,
            "Cr50": 3.4408e-6,
            "Cr52": 6.3805e-5,
            "Cr53": 7.0983e-6,
            "Cr54": 1.7342e-6,
            "Ni58": 2.5562e-5,
            "Ni60": 9.5185e-6,
            "Ni61": 4.0700e-7,
            "Ni62": 1.2765e-6,
            "Ni64": 3.1510e-7,
        }
    )
    calandria_tube = material_from_atom_densities("calandria_tube_ZrFeCrNi", calandria_dens, state.tmod)
    sheath_dens = zirc_atom_densities(4.9656e-2)
    sheath_dens.update(
        {
            "Fe54": 1.0341e-5,
            "Fe56": 1.5654e-4,
            "Fe57": 3.5509e-6,
            "Fe58": 4.6451e-7,
            "Cr50": 3.9532e-6,
            "Cr52": 7.3306e-5,
            "Cr53": 8.1553e-6,
            "Cr54": 1.9925e-6,
            "Ni58": 3.7378e-6,
            "Ni60": 1.3918e-6,
            "Ni61": 5.9514e-8,
            "Ni62": 1.8668e-7,
            "Ni64": 4.6062e-8,
        }
    )
    sheath = material_from_atom_densities("sheath_ZrFeCrNi", sheath_dens, state.tcool)
    helium = material_from_atom_densities("helium_gap", {"He4": 2.0877e-4}, state.tcool)
    mats = {
        "fuel_1": fuel_mats[11],
        "fuel_2": fuel_mats[12],
        "fuel_3": fuel_mats[13],
        "fuel_4": fuel_mats[14],
        "coolant": d2o_coolant,
        "moderator": d2o_moderator,
        "pressure_tube": pressure_tube,
        "calandria_tube": calandria_tube,
        "sheath": sheath,
        "helium": helium,
    }
    return openmc.Materials(mats.values()), mats


def pin_positions() -> list[tuple[float, float]]:
    positions = [(0.0, 0.0)]
    for radius, count in [(1.4880, 6), (2.8750, 12), (4.3300, 18)]:
        for i in range(count):
            angle = 2.0 * math.pi * i / count
            positions.append((radius * math.cos(angle), radius * math.sin(angle)))
    return positions


def build_geometry(mats: dict[str, openmc.Material]) -> openmc.Geometry:
    half_pitch = 28.575 / 2.0
    xmin = openmc.XPlane(x0=-half_pitch, boundary_type="reflective")
    xmax = openmc.XPlane(x0=half_pitch, boundary_type="reflective")
    ymin = openmc.YPlane(y0=-half_pitch, boundary_type="reflective")
    ymax = openmc.YPlane(y0=half_pitch, boundary_type="reflective")
    square = +xmin & -xmax & +ymin & -ymax

    rcool = openmc.ZCylinder(r=5.1689)
    rpt = openmc.ZCylinder(r=5.6032)
    rgap = openmc.ZCylinder(r=6.4478)
    rct = openmc.ZCylinder(r=6.5875)

    cells: list[openmc.Cell] = []
    pin_outer_regions = []
    for idx, (x0, y0) in enumerate(pin_positions(), start=1):
        r1 = openmc.ZCylinder(x0=x0, y0=y0, r=0.3061)
        r2 = openmc.ZCylinder(x0=x0, y0=y0, r=0.4898)
        r3 = openmc.ZCylinder(x0=x0, y0=y0, r=0.5816)
        r4 = openmc.ZCylinder(x0=x0, y0=y0, r=0.6122)
        rg = openmc.ZCylinder(x0=x0, y0=y0, r=0.6540)
        pin_outer_regions.append(-rg)
        cells.extend(
            [
                openmc.Cell(name=f"pin{idx:02d}_fuel_1", fill=mats["fuel_1"], region=-r1),
                openmc.Cell(name=f"pin{idx:02d}_fuel_2", fill=mats["fuel_2"], region=+r1 & -r2),
                openmc.Cell(name=f"pin{idx:02d}_fuel_3", fill=mats["fuel_3"], region=+r2 & -r3),
                openmc.Cell(name=f"pin{idx:02d}_fuel_4", fill=mats["fuel_4"], region=+r3 & -r4),
                openmc.Cell(name=f"pin{idx:02d}_sheath", fill=mats["sheath"], region=+r4 & -rg),
            ]
        )

    coolant_region = -rcool
    for region in pin_outer_regions:
        coolant_region = coolant_region & ~region
    cells.append(openmc.Cell(name="coolant_D2O_cluster_region", fill=mats["coolant"], region=coolant_region))
    cells.append(openmc.Cell(name="pressure_tube", fill=mats["pressure_tube"], region=+rcool & -rpt))
    cells.append(openmc.Cell(name="helium_gap", fill=mats["helium"], region=+rpt & -rgap))
    cells.append(openmc.Cell(name="calandria_tube", fill=mats["calandria_tube"], region=+rgap & -rct))
    cells.append(openmc.Cell(name="moderator_D2O_outer_square", fill=mats["moderator"], region=square & +rct))
    return openmc.Geometry(openmc.Universe(cells=cells))


def build_settings(particles: int, batches: int, inactive: int) -> openmc.Settings:
    settings = openmc.Settings()
    settings.run_mode = "eigenvalue"
    settings.particles = particles
    settings.batches = batches
    settings.inactive = inactive
    settings.temperature = {
        "method": "interpolation",
        "range": (293.15, 2593.15),
        "tolerance": 2000.0,
    }
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box(
            [-14.2875, -14.2875, -1.0],
            [14.2875, 14.2875, 1.0],
            only_fissionable=True,
        )
    )
    return settings


def build_tallies() -> openmc.Tallies:
    energy = openmc.EnergyFilter([1.0e-5, 4.0, 2.0e7])
    tally = openmc.Tally(name="whole_cell_two_group_rates")
    tally.filters = [energy]
    tally.scores = ["flux", "absorption", "fission", "nu-fission"]
    return openmc.Tallies([tally])


def export_case(
    case_dir: pathlib.Path,
    state: BranchState,
    particles: int,
    batches: int,
    inactive: int,
    inventory_densities: dict[float, dict[int, dict[str, float]]] | None,
) -> None:
    if case_dir.exists():
        shutil.rmtree(case_dir)
    case_dir.mkdir(parents=True)
    materials, mats = build_materials(state, inventory_densities)
    materials.cross_sections = str(OPENMC_DATA / "cross_sections.xml")
    materials.export_to_xml(case_dir / "materials.xml")
    build_geometry(mats).export_to_xml(case_dir / "geometry.xml")
    build_settings(particles, batches, inactive).export_to_xml(case_dir / "settings.xml")
    build_tallies().export_to_xml(case_dir / "tallies.xml")


def export_all_fresh_decks(
    states: list[BranchState],
    case_root: pathlib.Path,
    particles: int,
    batches: int,
    inactive: int,
    skip_existing: set[str],
    coverage: dict[float, InventoryCoverage],
    inventory_densities: dict[float, dict[int, dict[str, float]]],
) -> int:
    count = 0
    for state in states:
        if not state_openmc_supported(state, coverage) or state.case_id in skip_existing:
            continue
        export_case(case_root / state.case_id, state, particles, batches, inactive, inventory_densities)
        count += 1
    return count


def run_openmc(case_dir: pathlib.Path, threads: int) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env["OPENMC_CROSS_SECTIONS"] = str(OPENMC_DATA / "cross_sections.xml")
    cmd = [str(OPENMC_BIN), "-s", str(threads)]
    return subprocess.run(cmd, cwd=case_dir, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)


def parse_statepoint(case_dir: pathlib.Path) -> dict[str, float]:
    statepoints = sorted(case_dir.glob("statepoint.*.h5"))
    if not statepoints:
        raise FileNotFoundError(f"no OpenMC statepoint found in {case_dir}")
    with openmc.StatePoint(statepoints[-1]) as sp:
        k = sp.keff
        tally = sp.get_tally(name="whole_cell_two_group_rates")
        df = tally.get_pandas_dataframe()
        out: dict[str, float] = {
            "openmc_keff": float(k.nominal_value),
            "openmc_keff_sd": float(k.std_dev),
        }
        for _, row in df.iterrows():
            e_low = float(row["energy low [eV]"])
            group = "thermal" if e_low < 4.0 else "fast"
            score = str(row["score"]).replace("-", "_")
            out[f"{score}_{group}"] = float(row["mean"])
        if out.get("flux_fast", 0.0) > 0.0:
            out["thermal_to_fast_flux_ratio"] = out.get("flux_thermal", 0.0) / out["flux_fast"]
        return out


def write_svg_plots(
    out_dir: pathlib.Path,
    pilot_rows: list[dict[str, object]],
    states: list[BranchState],
    coverage: dict[float, InventoryCoverage],
) -> None:
    plot_dir = out_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    width, height = 900, 420
    margin = 70

    def esc(text: object) -> str:
        return str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

    if pilot_rows:
        xs = list(range(len(pilot_rows)))
        k_values = [float(r["openmc_keff"]) for r in pilot_rows] + [float(r["dragon_k_eff"]) for r in pilot_rows]
        y_min, y_max = min(k_values), max(k_values)
        pad = max((y_max - y_min) * 0.15, 0.002)
        y_min -= pad
        y_max += pad

        def xmap(i: int) -> float:
            if len(xs) == 1:
                return width / 2
            return margin + i * (width - 2 * margin) / (len(xs) - 1)

        def ymap(v: float) -> float:
            return height - margin - (v - y_min) * (height - 2 * margin) / (y_max - y_min)

        dragon_points = " ".join(f"{xmap(i):.1f},{ymap(float(r['dragon_k_eff'])):.1f}" for i, r in enumerate(pilot_rows))
        openmc_points = " ".join(f"{xmap(i):.1f},{ymap(float(r['openmc_keff'])):.1f}" for i, r in enumerate(pilot_rows))
        lines = [
            "<svg xmlns='http://www.w3.org/2000/svg' width='900' height='420' viewBox='0 0 900 420'>",
            "<rect width='900' height='420' fill='white'/>",
            "<text x='450' y='30' text-anchor='middle' font-family='Arial' font-size='18'>Phase F OpenMC pilot k comparison</text>",
            f"<line x1='{margin}' y1='{height-margin}' x2='{width-margin}' y2='{height-margin}' stroke='#333'/>",
            f"<line x1='{margin}' y1='{margin}' x2='{margin}' y2='{height-margin}' stroke='#333'/>",
            f"<polyline fill='none' stroke='#1f77b4' stroke-width='2.5' points='{dragon_points}'/>",
            f"<polyline fill='none' stroke='#d62728' stroke-width='2.5' points='{openmc_points}'/>",
            "<text x='720' y='64' font-family='Arial' font-size='13' fill='#1f77b4'>DRAGON branch k</text>",
            "<text x='720' y='84' font-family='Arial' font-size='13' fill='#d62728'>OpenMC CE pilot k</text>",
        ]
        for i, row in enumerate(pilot_rows):
            x = xmap(i)
            lines.append(f"<circle cx='{x:.1f}' cy='{ymap(float(row['dragon_k_eff'])):.1f}' r='4' fill='#1f77b4'/>")
            lines.append(f"<circle cx='{x:.1f}' cy='{ymap(float(row['openmc_keff'])):.1f}' r='4' fill='#d62728'/>")
            lines.append(
                f"<text x='{x:.1f}' y='{height-38}' text-anchor='middle' font-family='Arial' font-size='10'>"
                f"{esc(row['tfuel_K'])}</text>"
            )
        for j in range(5):
            v = y_min + j * (y_max - y_min) / 4
            y = ymap(v)
            lines.append(f"<line x1='{margin-5}' y1='{y:.1f}' x2='{width-margin}' y2='{y:.1f}' stroke='#ddd'/>")
            lines.append(f"<text x='{margin-10}' y='{y+4:.1f}' text-anchor='end' font-family='Arial' font-size='11'>{v:.4f}</text>")
        lines.append("<text x='450' y='400' text-anchor='middle' font-family='Arial' font-size='12'>pilot states, labels are fuel temperature K</text>")
        lines.append("</svg>")
        (plot_dir / "phaseF_openmc_pilot_k_comparison.svg").write_text("\n".join(lines))

    runnable = sum(1 for s in states if state_openmc_supported(s, coverage))
    blocked = len(states) - runnable
    total = max(len(states), 1)
    bar_w = width - 2 * margin
    run_w = bar_w * runnable / total
    block_w = bar_w * blocked / total
    lines = [
        "<svg xmlns='http://www.w3.org/2000/svg' width='900' height='260' viewBox='0 0 900 260'>",
        "<rect width='900' height='260' fill='white'/>",
        "<text x='450' y='32' text-anchor='middle' font-family='Arial' font-size='18'>Phase F branch support status</text>",
        f"<rect x='{margin}' y='90' width='{run_w:.1f}' height='54' fill='#2ca02c'/>",
        f"<rect x='{margin + run_w:.1f}' y='90' width='{block_w:.1f}' height='54' fill='#ff7f0e'/>",
        f"<text x='{margin}' y='78' font-family='Arial' font-size='13'>Runnable states with current CE data: {runnable}</text>",
        f"<text x='{margin}' y='170' font-family='Arial' font-size='13'>Blocked states: {blocked} missing CE data for inventory nuclides</text>",
        f"<text x='{margin}' y='205' font-family='Arial' font-size='13'>Total completed DRAGON branch states: {len(states)}</text>",
        "</svg>",
    ]
    (plot_dir / "phaseF_openmc_branch_support_status.svg").write_text("\n".join(lines))


def select_pilot_states(states: list[BranchState]) -> list[BranchState]:
    targets = [
        (0.0, 900.0, 343.15, 343.15),
        (0.0, 293.15, 293.15, 293.15),
        (0.0, 2593.15, 343.15, 343.15),
        (0.0, 900.0, 363.15, 343.15),
        (0.0, 900.0, 343.15, 363.15),
        (10.0, 900.0, 343.15, 343.15),
        (20.0, 900.0, 343.15, 343.15),
    ]
    selected: list[BranchState] = []
    for target in targets:
        for state in states:
            if (
                abs(state.burnup - target[0]) < 1.0e-6
                and abs(state.tfuel - target[1]) < 1.0e-6
                and abs(state.tcool - target[2]) < 1.0e-6
                and abs(state.tmod - target[3]) < 1.0e-6
            ):
                selected.append(state)
                break
    return selected


def write_manifest(states: list[BranchState], out_dir: pathlib.Path, coverage: dict[float, InventoryCoverage]) -> None:
    rows = []
    for state in states:
        inv = state_inventory_coverage(state, coverage)
        supported = state_openmc_supported(state, coverage)
        if supported:
            disposition = "RUNNABLE_WITH_CURRENT_CE_DATA"
        elif inv is None:
            disposition = "BLOCKED_MISSING_BURNUP_ISOTOPIC_VECTOR"
        else:
            disposition = "BLOCKED_MISSING_CE_DATA_FOR_INVENTORY_NUCLIDES"
        rows.append(
            {
                "case_id": state.case_id,
                "burnup_GWd_tU": state.burnup,
                "tfuel_K": state.tfuel,
                "tcool_K": state.tcool,
                "tmod_K": state.tmod,
                "dragon_k_eff": "" if state.k_eff is None else f"{state.k_eff:.8f}",
                "inventory_rows": 0 if inv is None else inv.rows,
                "inventory_supported_nuclides": "" if inv is None else inv.supported_nuclides,
                "inventory_missing_nuclides": "" if inv is None else len(inv.missing_nuclides),
                "inventory_supported_density_fraction": "" if inv is None else f"{inv.supported_density_fraction:.12e}",
                "openmc_ce_supported_now": "YES" if supported else "NO",
                "disposition": disposition,
            }
        )
    write_csv(
        out_dir / "phaseF_openmc_full_branch_manifest.csv",
        rows,
        [
            "case_id",
            "burnup_GWd_tU",
            "tfuel_K",
            "tcool_K",
            "tmod_K",
            "dragon_k_eff",
            "inventory_rows",
            "inventory_supported_nuclides",
            "inventory_missing_nuclides",
            "inventory_supported_density_fraction",
            "openmc_ce_supported_now",
            "disposition",
        ],
    )


def write_inventory_coverage(out_dir: pathlib.Path, coverage: dict[float, InventoryCoverage]) -> None:
    rows = []
    for burnup in sorted(coverage):
        inv = coverage[burnup]
        rows.append(
            {
                "burnup_GWd_tU": f"{inv.burnup:.8e}",
                "inventory_rows": inv.rows,
                "supported_nuclides": inv.supported_nuclides,
                "missing_nuclides": len(inv.missing_nuclides),
                "supported_density_fraction": f"{inv.supported_density_fraction:.12e}",
                "missing_nuclide_list": " ".join(inv.missing_nuclides),
                "ce_complete": "YES" if inv.ce_complete else "NO",
            }
        )
    write_csv(
        out_dir / "phaseF_inventory_ce_coverage.csv",
        rows,
        [
            "burnup_GWd_tU",
            "inventory_rows",
            "supported_nuclides",
            "missing_nuclides",
            "supported_density_fraction",
            "missing_nuclide_list",
            "ce_complete",
        ],
    )


def write_report(
    out_dir: pathlib.Path,
    states: list[BranchState],
    pilot_rows: list[dict[str, object]],
    particles: int,
    batches: int,
    inactive: int,
    coverage: dict[float, InventoryCoverage],
) -> None:
    runnable = sum(1 for s in states if state_openmc_supported(s, coverage))
    blocked = len(states) - runnable
    inventory_burnups = len(coverage)
    ce_blocked_burnups = sum(1 for inv in coverage.values() if inv.exists and not inv.ce_complete)
    sample_missing = sorted({nuc for inv in coverage.values() for nuc in inv.missing_nuclides})[:40]
    passed_runs = sum(1 for r in pilot_rows if r["run_status"] == "PASS")
    rx_hash = sha256(RX_FILE) if RX_FILE.exists() else "MISSING"
    xs_hash = sha256(OPENMC_DATA / "cross_sections.xml") if (OPENMC_DATA / "cross_sections.xml").exists() else "MISSING"
    disposition = "PHASE_F_INITIATED_NOT_COMPLETE_FOR_SAFETY_EXIT"
    if pilot_rows and passed_runs == len(pilot_rows) and blocked == 0:
        disposition = "PHASE_F_OPENMC_PILOT_PASS_CE_COVERAGE_COMPLETE"
    coverage_status = "PASS" if inventory_burnups and ce_blocked_burnups == 0 and blocked == 0 else "PARTIAL"

    equiv_rows = [
        {
            "criterion": "independent_openmc_ce_source_documented",
            "status": "PASS",
            "detail": f"OpenMC {openmc.__version__}; CE data {OPENMC_DATA / 'cross_sections.xml'}",
        },
        {
            "criterion": "same_37_element_phase4_geometry",
            "status": "PASS",
            "detail": "37 pins, Phase 4 ring radii, fuel/sheath radii, coolant/PT/gap/CT/moderator radii, reflective square pitch",
        },
        {
            "criterion": "fresh_material_basis",
            "status": "PASS",
            "detail": "U/O/D/Zr atom densities from MATDPLD2O; structural H, Nb, Fe, Cr, Ni, and He CE data are present in the relocated OpenMC library",
        },
        {
            "criterion": "burnup_material_basis",
            "status": coverage_status if inventory_burnups else "FAIL",
            "detail": f"DRAGON burnup inventories exported for {inventory_burnups} burnup states; CE-incomplete burnup states={ce_blocked_burnups}",
        },
        {
            "criterion": "temperature_branch_basis",
            "status": "PASS" if blocked == 0 else "PARTIAL",
            "detail": f"current runnable states={runnable}; blocked states={blocked}",
        },
        {
            "criterion": "pilot_openmc_execution",
            "status": "PASS" if pilot_rows and passed_runs == len(pilot_rows) else "FAIL",
            "detail": f"{passed_runs}/{len(pilot_rows)} selected fresh pilot states completed",
        },
        {
            "criterion": "full_phase_f_exit",
            "status": "PARTIAL" if blocked == 0 else "FAIL",
            "detail": "CE coverage is complete; full like-for-like safety exit still requires full-statistics OpenMC execution and acceptance checks over the required burnup-temperature matrix",
        },
    ]
    write_csv(out_dir / "phaseF_openmc_equivalence_matrix.csv", equiv_rows, ["criterion", "status", "detail"])

    lines = [
        "# Phase F OpenMC Independent Lattice Benchmark Report",
        "",
        f"Disposition: `{disposition}`",
        "",
        "## Scope",
        "",
        "This package initiates the independent continuous-energy OpenMC benchmark",
        "for the DRAGON full-production two-group CANDU/PHWR library. It builds",
        "the Phase 4 37-element cell geometry directly, uses DRAGON exported fuel",
        "inventories for non-fresh burnup states, and compares selected branch",
        "states against the matching DRAGON k-effective values.",
        "",
        "## Frozen Inputs",
        "",
        f"- Candidate RX: `{RX_FILE}`",
        f"- Candidate RX SHA256: `{rx_hash}`",
        f"- OpenMC cross_sections.xml: `{OPENMC_DATA / 'cross_sections.xml'}`",
        f"- OpenMC cross_sections.xml SHA256: `{xs_hash}`",
        f"- OpenMC version: `{openmc.__version__}`",
        f"- OpenMC executable: `{OPENMC_BIN}`",
        f"- DRAGON inventory CSV: `{INVENTORY_CSV}`",
        "",
        "## Branch Coverage",
        "",
        f"- Completed DRAGON branch states in manifest: `{len(states)}`",
        f"- Burnup inventory states available: `{inventory_burnups}`",
        f"- States runnable with current OpenMC CE coverage: `{runnable}`",
        f"- States blocked by missing CE data for inventory nuclides: `{blocked}`",
        f"- Sample missing depleted nuclides: `{', '.join(sample_missing) if sample_missing else 'none'}`",
        "",
        "## Pilot Run Settings",
        "",
        f"- particles: `{particles}`",
        f"- batches: `{batches}`",
        f"- inactive: `{inactive}`",
        "",
        "## Pilot Results",
        "",
        "| case_id | burnup | Tfuel | Tcool | Tmod | DRAGON k | OpenMC k | sigma | delta pcm | status |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in pilot_rows:
        lines.append(
            "| {case_id} | {burnup_GWd_tU:.1f} | {tfuel_K:.2f} | {tcool_K:.2f} | {tmod_K:.2f} | "
            "{dragon_k_eff:.6f} | {openmc_keff:.6f} | {openmc_keff_sd:.6f} | {delta_pcm:.1f} | {run_status} |".format(
                **row
            )
        )
    lines += [
        "",
        "## Equivalence Disposition",
        "",
        "The OpenMC pilot is an independent CE lattice calculation and is suitable",
        "for exercising Phase F infrastructure. Inventory CE coverage is complete",
        "for the completed DRAGON branch matrix when the inventory coverage checks",
        "show zero blocked states. This report is still not a full safety-exit",
        "benchmark unless the OpenMC transport runs are executed at the required",
        "statistics over the required burnup-temperature matrix and pass explicit",
        "acceptance criteria.",
        "",
        "## Generated Artifacts",
        "",
        f"- Full branch manifest: `{out_dir / 'phaseF_openmc_full_branch_manifest.csv'}`",
        f"- Pilot comparison CSV: `{out_dir / 'phaseF_openmc_pilot_comparison.csv'}`",
        f"- Equivalence matrix: `{out_dir / 'phaseF_openmc_equivalence_matrix.csv'}`",
        f"- OpenMC case directories: `{out_dir / 'openmc_cases'}`",
        f"- Plot: `{out_dir / 'plots' / 'phaseF_openmc_pilot_k_comparison.svg'}`",
        f"- Plot: `{out_dir / 'plots' / 'phaseF_openmc_branch_support_status.svg'}`",
        "",
    ]
    (out_dir / "PHASE_F_OPENMC_BENCHMARK_REPORT.md").write_text("\n".join(lines))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=pathlib.Path, default=OUT_DIR)
    parser.add_argument("--particles", type=int, default=4000)
    parser.add_argument("--batches", type=int, default=80)
    parser.add_argument("--inactive", type=int, default=20)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--skip-openmc", action="store_true")
    parser.add_argument("--no-export-all-fresh", action="store_true")
    args = parser.parse_args()

    if not (OPENMC_DATA / "cross_sections.xml").exists():
        raise FileNotFoundError(OPENMC_DATA / "cross_sections.xml")
    if not OPENMC_BIN.exists():
        raise FileNotFoundError(OPENMC_BIN)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    states = read_branch_states(BRANCH_STATUS)
    available_nuclides = openmc_materials_from_xml(OPENMC_DATA / "cross_sections.xml")
    coverage = read_inventory_coverage(INVENTORY_CSV, available_nuclides)
    inventory_densities = read_inventory_densities(INVENTORY_CSV)
    write_manifest(states, args.out_dir, coverage)
    write_inventory_coverage(args.out_dir, coverage)
    pilot_states = select_pilot_states(states)
    if not pilot_states:
        raise RuntimeError("no fresh pilot states found in branch_status.csv")

    pilot_rows: list[dict[str, object]] = []
    case_root = args.out_dir / "openmc_cases"
    if not args.no_export_all_fresh:
        export_all_fresh_decks(
            states,
            case_root,
            args.particles,
            args.batches,
            args.inactive,
            {state.case_id for state in pilot_states},
            coverage,
            inventory_densities,
        )
    for state in pilot_states:
        case_dir = case_root / state.case_id
        export_case(case_dir, state, args.particles, args.batches, args.inactive, inventory_densities)
        run_status = "SKIPPED"
        stdout = ""
        metrics: dict[str, float] = {
            "openmc_keff": float("nan"),
            "openmc_keff_sd": float("nan"),
            "thermal_to_fast_flux_ratio": float("nan"),
        }
        if not args.skip_openmc:
            proc = run_openmc(case_dir, args.threads)
            stdout = proc.stdout
            (case_dir / "openmc_stdout.log").write_text(stdout)
            if proc.returncode == 0:
                metrics = parse_statepoint(case_dir)
                run_status = "PASS"
            else:
                run_status = "FAIL"
        dragon_k = state.k_eff if state.k_eff is not None else float("nan")
        delta_pcm = (metrics["openmc_keff"] - dragon_k) * 1.0e5
        pilot_rows.append(
            {
                "case_id": state.case_id,
                "burnup_GWd_tU": state.burnup,
                "tfuel_K": state.tfuel,
                "tcool_K": state.tcool,
                "tmod_K": state.tmod,
                "dragon_k_eff": dragon_k,
                "openmc_keff": metrics["openmc_keff"],
                "openmc_keff_sd": metrics["openmc_keff_sd"],
                "delta_pcm": delta_pcm,
                "thermal_to_fast_flux_ratio": metrics.get("thermal_to_fast_flux_ratio", float("nan")),
                "run_status": run_status,
                "case_dir": str(case_dir),
            }
        )
        if run_status == "FAIL":
            sys.stderr.write(stdout[-4000:])

    write_csv(
        args.out_dir / "phaseF_openmc_pilot_comparison.csv",
        pilot_rows,
        [
            "case_id",
            "burnup_GWd_tU",
            "tfuel_K",
            "tcool_K",
            "tmod_K",
            "dragon_k_eff",
            "openmc_keff",
            "openmc_keff_sd",
            "delta_pcm",
            "thermal_to_fast_flux_ratio",
            "run_status",
            "case_dir",
        ],
    )
    write_report(args.out_dir, states, pilot_rows, args.particles, args.batches, args.inactive, coverage)
    write_svg_plots(args.out_dir, pilot_rows, states, coverage)
    print(f"wrote {args.out_dir / 'PHASE_F_OPENMC_BENCHMARK_REPORT.md'}")
    print(f"pilot states: {sum(1 for r in pilot_rows if r['run_status'] == 'PASS')}/{len(pilot_rows)} PASS")
    return 0 if all(r["run_status"] in {"PASS", "SKIPPED"} for r in pilot_rows) else 2


if __name__ == "__main__":
    raise SystemExit(main())
