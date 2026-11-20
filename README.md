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
| `data/njoy_dragon_comparison.csv` | NJOY21 and DRAGON two-group constants of the reactivity-device states |
| `dragon_inputs/benchmark_library/` | DRAGON 5.1 decks of the benchmarked library (one per burnup) with procedures and DRAGLIB access script |
| `openmc_inputs/` | OpenMC 0.15 inputs for the fresh nominal state; other states change compositions (`fuel_isotopic_inventories.csv`) and temperatures |
| `source/` | Scripts that compute the published results from `data/` (Python 3 standard library) |

DRAGON variants (`variant` column): `corrected_refined_mesh_J2_matched_inventory` is the benchmarked library (JEF-2.2, 20-ring moderator, OpenMC inventories). `corrected_J2` and `corrected_E81` use a single-region moderator (JEF-2.2, ENDF/B-VIII.1). `corrected_refined_mesh_E81` uses ENDF/B-VIII.1 with a 20-ring moderator. `corrected_refined_mesh_J2_own_depletion` uses JEF-2.2 depleted on the 20-ring mesh. The ENDF/B-VIII.1 variants deplete with their own data.

## Reproducing the results

```bash
python source/analyze_benchmark_library.py   # Section 3 results -> data/benchmark_library_*.csv
python source/analyze_corrected_matrix.py    # mesh and evaluated-data sensitivity -> data/corrected_matrix_*.csv
```

Δk = 10⁵ (k_OpenMC − k_DRAGON). Transport calculations used DRAGON 5.1.0 (DRAGLIBs `draglibJef2p2_v5p1`, `draglibendfb8r1_v5p1`), OpenMC 0.15.3/0.15.4-dev169 with ENDF/B-VIII.1, and NJOY21 1.2.74.

License: MIT.
