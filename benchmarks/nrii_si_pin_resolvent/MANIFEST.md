# Si PIN resolvent fast-mode elimination manifest

This directory supports `NRII_SI_PIN_RESOLVENT_FAST_MODE_ELIMINATION.md`.

## Stored files

- `nrii_si_pin_transport_resolvent.c`: exact contact + SG transport resolvent implementation used for the final result.
- `nrii_si_pin_resolvent_summary.csv`: measured operating-radius summary against bare and contact-only NRII.
- `nrii_si_pin_tc_resolvent_long_combined.csv`: complete long-step transport/contact resolvent sweep for N=61,121,241,481 through 1e-9 s.
- `nrii_si_pin_resolvent_theta_summary.csv`: physical charge-relaxation matching-factor summary.
- `nrii_si_pin_resolvent_decision.json`: machine-readable interpretation and limits.

The earlier bare physical-device NRII data are stored with the physical Si PIN follow-up already on `main`.

## Validation rules

Every claimed resolvent result is checked against the complete original backward-Euler residual and against the direct coupled block-Newton solution. The transformed equation is algebraically identical to the original implicit relation when the fixed resolvent operator is invertible.

## Important boundary

The exact 1D transport resolvent uses tridiagonal linear solves at every coefficient order. This resolves the tested convergence-radius problem but does not establish a 1D timing win, and a direct 3D analogue could reintroduce an expensive global solve. The next target is a locality-preserving approximate resolvent.
