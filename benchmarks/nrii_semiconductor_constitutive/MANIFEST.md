# NRII semiconductor constitutive benchmark manifest

This directory supports `NRII_SEMICONDUCTOR_CONSTITUTIVE_CAMPAIGN.md`.

## Stored evidence

- `nrii_semiconductor_constitutive_summary.csv`: model-level statistics aggregated from the 1,800-point validated main sweep.
- `nrii_semiconductor_scale_summary.csv`: model-by-grid-size statistics from the same 1,800-point main sweep.
- `nrii_semiconductor_focus_impact_ionization_s18.csv`: complete 72-point impact-ionization boundary sweep, including points beyond the NRII accuracy boundary.
- `nrii_semiconductor_focus_combined_pin_electrothermal_traps_s18.csv`: complete 72-point dense electrothermal/trap/impact boundary sweep.
- `nrii_semiconductor_boundary_summary.csv`: valid/invalid counts and boundary statistics for the two focused sweeps.
- `nrii_semiconductor_constitutive_decision.json`: machine-readable campaign summary and current interpretation.

The associated root report is `NRII_SEMICONDUCTOR_CONSTITUTIVE_CAMPAIGN.md`.

## Campaign scope

The complete campaign contains 1,944 measured points: 1,800 validated main-sweep points and 144 boundary-stress points. The main sweep spans five cubic 3D grid sizes, twelve constitutive families, five constitutive-strength levels, and six implicit-step levels.

The main acceptance gates are:

- relative NRII residual <= 1e-8;
- relative NRII/Newton difference <= 1e-7.

All 1,800 main-sweep points pass both gates. Boundary-stress points that fail accuracy remain in the stored boundary tables and must not be counted as valid speedups.

## Interpretation boundary

These are controlled normalized constitutive/transport benchmarks designed to expose solver-mechanism behavior. They are not production TCAD validation. The current evidence does not include a self-consistent elliptic Poisson solve, Scharfetter-Gummel flux discretization, calibrated silicon/GaN material parameters, production AMG/ILU preconditioning, or device-level IV/breakdown/thermal validation.

The principal structural result is that low coefficient-tail ratio and locality are not sufficient to predict NRII speed. Per-order constitutive coefficient-algebra cost is an independent performance coordinate. Large explicit trap-state stacks and dense electrothermal/trap/impact combinations can erase the advantage even while q_eff remains small.
