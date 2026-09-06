# NRII golden-region evidence

This directory supports `NRII_GOLDEN_REGION.md`.

The files in this directory are compact evidence summaries rather than a replacement for the complete benchmark sources and reproducibility bundles already retained by the project.

- `nrii_golden_region_evidence.csv` consolidates measured speed ranges and qualifications from the new 3D localization campaign, PCG fairness controls, the localized-hotspot study, and the existing operating-envelope suite.
- `nrii_3d_globalization_summary.csv` summarizes the new 3D step-sweep and concrete-model measurements by scenario.
- `nrii_3d_localization_classification.csv` records which tested 3D structures can remove solver-induced globality and which remain intrinsically nonlocal.

The speed ratio convention is `reference solver runtime / NRII runtime`; values above one favor NRII. Controlled pathological results and microbenchmark timing ranges must retain their stated qualifications.
