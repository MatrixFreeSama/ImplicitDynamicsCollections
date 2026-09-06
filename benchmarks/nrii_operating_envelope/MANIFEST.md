# NRII operating-envelope benchmark manifest

This directory supports `NRII_OPERATING_ENVELOPE.md`.

## Files stored directly

- `nrii_final_comprehensive_summary.csv`: benchmark-group summary.
- `nrii_operating_envelope_reproducibility_bundle.zip`: complete reproducibility snapshot.

## Contents of the reproducibility bundle

The archive contains:

- the consolidated 572-point phase table and machine-readable decision-rule snapshot;
- all parameter-sweep CSV files used by the report;
- the six C benchmark sources used for the baseline, sparsity, curvature, nonnormality, stiffness/family, and multidimensional/preconditioner studies;
- the report source and paper-oriented SVG figures.

The benchmark programs are research microbenchmarks intended to expose solver-regime behavior under controlled conditions. They do not substitute for production AMG/ILU libraries or application-specific finite-element codes.
