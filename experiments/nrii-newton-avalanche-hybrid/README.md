# NRII -> Newton hybrid avalanche root-basin cloud study

## Purpose

This branch tests a hybrid interpretation of NRII: NRII is not required to finish as the final root solver. Instead, its coefficient closure + Chebyshev audit + Padé continuation is used to construct the branch-informed state that Newton needs as an initial seed. Newton then performs the residual/Jacobian correction.

This is deliberately different from standalone NRII. The predictor is allowed to have a residual far above the final tolerance, as long as it passes coefficient/rational safety and positivity gates.

## Model

Dimensionless 1D PN drift-diffusion-Poisson model with smooth reduced impact ionization, identical to the prior avalanche root-basin experiment.

- nominal global grid: 2^34 cells
- local materialized window: 1537 cells
- dt = 0.005
- final residual tolerance = 1e-10
- Newton: analytic 3x3 block-tridiagonal Jacobian, block Thomas solve, positivity cap, Armijo backtracking

## Locked hybrid policy

Declared NRII predictor orders:

16, 24, 32, 48, 64, 80, 96, 112, 128, 160, 192.

The handoff order is the first declared order satisfying all of the following without consulting the original-equation residual:

1. coefficient bank remains finite;
2. Chebyshev tail <= sqrt(tolerance);
3. Padé denominator gate is safe;
4. predicted electron and hole concentrations remain positive.

Only after the seed has been selected is its original-equation residual measured for reporting. Therefore the hybrid handoff is not residual-tuned.

## Main result

Pure old-state Newton passes through reverse drive 77 and fails from drive 78 onward in this sweep. The NRII-seeded Newton hybrid passes through drive 87 and fails at drive 88.

- pure Newton last passing drive: 77
- hybrid last passing drive: 87

The important observation is that the NRII seed does not need a small residual. At drive 78, pure Newton starts from residual 2.35607 and fails after 40 updates, while the locked NRII seed has an even larger residual 3.19066; nevertheless Newton converges from that seed in 3 full updates to 5.976e-12.

The same effect strengthens deeper into breakdown:

- drive 80: seed residual 3.31661, hybrid converges in 3 updates to 1.945e-14; pure Newton fails.
- drive 84: seed residual 10.2702, hybrid converges in 4 updates to 3.434e-14; pure Newton fails.
- drive 87: seed residual 78.8785, hybrid converges in 5 updates to 6.938e-11; pure Newton fails.
- drive 88: the first coefficient/rational-safe seed no longer places Newton in a convergent basin, and the hybrid fails.

## Interpretation

This experiment separates "being close in residual norm" from "being in a useful Newton basin". NRII can produce a state whose residual is numerically worse than the old-state initial guess but whose branch information places Newton in a basin where full Newton steps work.

That supports a hybrid role for NRII: produce structured branch-aware information for Newton, rather than require the raw NRII output itself to satisfy the final root tolerance.

This does not prove a universal basin advantage. The extension is finite in this model (78 -> 87), and the hybrid fails by drive 88.

## Timing note

The cloud predictor is a dense CPU reference implementation and is not the optimized GPU NRII path. Predictor-generation wall time is therefore included only for reproducibility, not as a performance claim about the production NRII implementation. The robust-basin result is the primary result of this study.
