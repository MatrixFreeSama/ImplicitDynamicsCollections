# NRII semiconductor Poisson and physical Si PIN follow-up manifest

This directory supports two root-level reports:

- `NRII_SEMICONDUCTOR_POISSON_FOLLOWUP.md`
- `NRII_SI_PIN_PHYSICAL_DEVICE_FOLLOWUP.md`

## Stored readable evidence

- `nrii_si_pin_device_breakdown_summary.csv`
- `nrii_si_pin_transient_summary.csv`
- `nrii_semiconductor_poisson_decision.json`
- `nrii_si_pin_physical_device_decision.json`

The reports summarize the complete normalized Poisson campaign (390 points) and the physical-unit Si PIN device campaign, including breakdown refinement, coefficient-radius measurements, contact-mode condensation, and NRII-to-Newton relay tests.

## Evidence rules

1. Physical Poisson globality is retained; NRII is not described as global-free.
2. The physical Si PIN breakdown voltage is an Okuto-Crowell-style ionization-integral crossing estimate, not a full avalanche-multiplication IV solution.
3. Full-domain NRII timing points that fail the complete-system residual gate are not valid speedups.
4. Contact-condensed NRII is classified as an approximate seed until the complete uncondensed residual passes.
5. Relay timing includes both the NRII seed and the full Newton correction.
6. The 1D Newton baseline uses structured block-tridiagonal direct linear algebra and is intentionally strong.

A complete local reproducibility bundle was generated with the experiment and is retained separately from this lightweight repository summary.
