# Validation_Data

Data and inputs for the benchmark of a DRAGON5 two-group branch library for a 37-element natural-uranium PHWR lattice against OpenMC: 2132 lattice states, 41 burnups (0–20 GWd/tU) × 52 temperature states (fuel 293–2593 K; coolant and moderator 293–363 K at fixed density). Both codes use the same stepwise depletion inventories.

## Files

| Path | Content |
|---|---|
| `data/openmc_reference_states.csv` | OpenMC k-eff and standard deviation for every state |
| `data/corrected_matrix_dragon_states.csv` | DRAGON k-eff and two-group constants (cm⁻¹) for every state and variant (below) |
| `data/openmc_two_group_rate_tallies.csv` | OpenMC whole-cell two-group tallies (boundary 4 eV) and rate coefficients |
| `data/fuel_isotopic_inventories.csv` | Atom densities (b⁻¹cm⁻¹) in the four fuel zones at the 41 burnups, used by both codes |
| `data/dragon_moderator_mesh_convergence.csv` | Moderator-mesh (1–40 rings) and tracking study, fresh nominal state |
| `dragon_inputs/benchmark_library/` | DRAGON 5.1 decks of the benchmarked library (one per burnup) with procedures and DRAGLIB access script |
| `openmc_inputs/` | OpenMC 0.15 inputs for the fresh nominal state; other states change compositions (`fuel_isotopic_inventories.csv`) and temperatures |
| `source/` | Core result-analysis scripts and supporting model/run/plotting tools (dependencies below) |

DRAGON variants (`variant` column): `corrected_refined_mesh_J2_matched_inventory` is the benchmarked library (JEF-2.2, 20-ring moderator, OpenMC inventories). `corrected_J2` and `corrected_E81` use a single-region moderator (JEF-2.2, ENDF/B-VIII.1). `corrected_refined_mesh_E81` uses ENDF/B-VIII.1 with a 20-ring moderator. `corrected_refined_mesh_J2_own_depletion` uses JEF-2.2 depleted on the 20-ring mesh. The ENDF/B-VIII.1 variants deplete with their own data.

## Reproducing the results

```bash
python source/analyze_benchmark_library.py   # Section 3 results -> data/benchmark_library_*.csv
python source/analyze_corrected_matrix.py    # mesh and evaluated-data sensitivity -> data/corrected_matrix_*.csv
```

Δk = 10⁵ (k_OpenMC − k_DRAGON). Transport calculations used DRAGON 5.1.0 (DRAGLIBs `draglibJef2p2_v5p1`, `draglibendfb8r1_v5p1`) and OpenMC 0.15.3/0.15.4-dev169 with ENDF/B-VIII.1.

License: MIT.

## Current paper revision (21 November 2026)

This bundle includes the corrected tally moments/statistics, fitted temperature-response windows, and the 63-state operating-reference extension. The Zenodo deposit is maintained separately from this Git repository. The manuscript and submission documents remain in the separate local submission directory.

`data/` contains the accepted paper tables and the derived tables needed by its checks; partial-run exports are excluded. `source/` contains model generators, analysis and plotting scripts. `qa/` contains numerical checks and diagnostics. `non_upload/operating_coefficients/` contains only selected accepted input files and DRAGON deck/procedure snapshots; this historical directory name is retained so input-hash records and scripts resolve correctly. It contains no statepoint binaries or revision backups.

Run these checks from the bundle root with Python 3.10 or newer:

```bash
python3 qa/check_tally_statistics.py
python3 qa/check_temperature_response_windows.py
python3 qa/check_operating_results.py
```

The first two checks and the operating-results check use the Python standard library. Core benchmark analysis also uses the standard library. Plotting needs numpy and matplotlib; transport-output analysis/audits need numpy and h5py; operating input generation needs iapws. Transport runs additionally require separately installed DRAGON 5.1 and its DRAGLIBs, OpenMC and evaluated cross-section libraries. Some original run scripts retain machine-specific paths; configure those before transport reruns. Complete statepoints, processed-data binaries, and the original transport workspace are not included, so the bundle does not by itself reproduce all transport calculations.

`qa/paper_bundle_manifest.json` records SHA-256 hashes of the files copied from the current local submission revision. The accepted input hashes are retained separately in `qa/operating_input_hashes.json`.

The retained `dragon_inputs/moderator_mesh_study/` snapshots support the moderator-mesh sensitivity study; their hashes and legacy archive provenance are in `qa/retained_mesh_input_manifest.json`. The unrelated legacy reactivity-device/NJOY inputs remain only in the recovery archive outside this project.
