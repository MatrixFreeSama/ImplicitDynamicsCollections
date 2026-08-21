# General nonlinear dynamics Newton reciprocal-alternation tournament, bare C17

## Purpose

This benchmark keeps the same 128-DOF manufactured nonlinear structural-dynamics implicit problem and the same solver pairings as the previous Newton-hybrid tournament, but changes the composition rule from one-way/terminal hybridization to strict reciprocal alternation.

For paired methods the schedule is exactly 1:1:

```text
A -> B -> A -> B -> ...
```

The experiment measures the largest contiguous operational convergence envelope on the same 121-point logarithmic stress grid from `0.05` to `1e12`.

## Fixed problem

The discrete nonlinear equation is

```text
R(q) = q + f_int(q) - f_ext = 0
```

with the same 128-DOF cubic-hardening chain and manufactured root used by the directional tournament. Success requires both relative residual <= `1e-12` and relative error to the manufactured root <= `1e-8`.

## Solver schedules

- `Newton`: Newton only.
- `Newton<->LineSearch`: raw Newton alternates with Armijo line-search Newton.
- `Newton<->Aitken`: raw Newton alternates with an Aitken-relaxed Newton-map stage.
- `Newton<->Anderson`: raw Newton alternates with depth-3 Anderson combination of Newton-map history.
- `Broyden<->Newton`: inverse-Broyden alternates with exact Newton; each Newton stage refreshes the inverse Jacobian used by the next Broyden stage.
- `PseudoTransient<->Newton`: shifted-Jacobian pseudo-transient stage alternates with exact Newton.
- `NRII<->Newton`: re-anchored order-12 NRII construction alternates with exact Newton. The NRII stage uses `q(h)=S+h[(f_ext-S)-f_int(q(h))]`, which recovers the original target at `h=1`; residual is an accept/reject audit only and never generates an NRII correction.

## Primary result: equal total-stage budget

At 24 total solver stages, every method receives the same total number of stage opportunities. A paired method therefore receives at most 12 stages from each member, while standalone Newton receives 24 Newton stages.

| Method | contiguous limit | relative to Newton |
|---|---:|---:|
| Newton | 38.0942859 | 1x |
| Newton<->LineSearch | 38.0942859 | 1x |
| Newton<->Aitken | 136.480128 | 3.58269x |
| Newton<->Anderson | 2.29931649 | 0.0603586x |
| Broyden<->Newton | 38.0942859 | 1x |
| PseudoTransient<->Newton | 29.51331 | 0.774744x |
| NRII<->Newton | 2.96784127 | 0.0779078x |

The largest contiguous envelope under this equal-stage budget is `Newton<->Aitken` at `136.480128`, about `3.58269x` the Newton baseline. `Broyden<->Newton` reaches the same contiguous limit as Newton. `LineSearch`, pseudo-transient, Anderson, and NRII reciprocal alternation do not improve the Newton envelope under this schedule and budget.

## Budget sensitivity

The same 1:1 schedules were also swept with 48 total stages. This is not cost-equivalent to the 24-stage experiment; it is included to show how strongly an "operational convergence limit" depends on the allowed stage budget.

| Method | 24 stages | 48 stages | 48/24 |
|---|---:|---:|---:|
| Newton | 38.0942859 | 4862.46236 | 127.643x |
| Newton<->LineSearch | 38.0942859 | 4862.46236 | 127.643x |
| Newton<->Aitken | 136.480128 | 103982.26 | 761.886x |
| Newton<->Anderson | 2.29931649 | 2.29931649 | 1x |
| Broyden<->Newton | 38.0942859 | 4862.46236 | 127.643x |
| PseudoTransient<->Newton | 29.51331 | 3767.16234 | 127.643x |
| NRII<->Newton | 2.96784127 | 38.0942859 | 12.8357x |

At 48 stages, `Newton<->Aitken` remains the widest alternating envelope at `103982.26`, about `21.3847x` the 48-stage Newton baseline.

This budget sweep also demonstrates that this benchmark is measuring an **operational envelope under a declared iteration/stage budget**, not a theorem about the mathematical basin of Newton. In this monotone cubic-hardening family, increasing the allowed Newton stages itself moves the frontier substantially.

## Directional versus reciprocal composition

The previous 24-iteration directional tournament is retained as a reference. Reciprocal alternation is not automatically better:

- `Broyden<->Newton` expands the contiguous limit from the directional Broyden->Newton result `2.29931649` to `38.0942859`.
- `Newton<->LineSearch` collapses the huge directional line-search envelope because every other stage is an unguarded raw Newton step.
- `PseudoTransient<->Newton` likewise loses the large directional pseudo-transient globalization envelope when raw Newton is forced every other stage.
- `NRII<->Newton` reaches only `2.96784127` under 24 total stages, below both standalone Newton and the previous one-way NRII->Newton result. The re-anchored NRII stages do not compensate for halving the number of Newton stages in this benchmark.

Therefore the data support a narrower but useful conclusion:

> Reciprocal solver alternation changes the convergence envelope in a strongly pair-dependent way. Alternation can enlarge, preserve, or shrink the operational convergence region; it is not itself a universal globalization mechanism.

## Reproducibility

Primary build:

```bash
gcc -O3 -march=native -std=c17 -Wall -Wextra -pedantic newton_reciprocal_alternation_dynamics_tournament.c -lm -o tournament
./tournament
```

For the 48-stage sweep:

```bash
gcc -O3 -march=native -std=c17 -DMAX_STAGES=48 -DNO_TIMING newton_reciprocal_alternation_dynamics_tournament.c -lm -o tournament48
./tournament48
```

A 24-stage ASan + UBSan sweep completed with zero sanitizer stderr bytes. The warning-clean optimized build produced zero compiler-warning bytes.

## Files

- `newton_reciprocal_alternation_dynamics_tournament.c`: complete C17 source.
- `CONVERGENCE_ENVELOPE_24_STAGES.csv`: primary equal-stage envelope.
- `CONVERGENCE_ENVELOPE_48_STAGES.csv`: budget-sensitivity envelope.
- `DIRECTIONAL_VS_ALTERNATING_24.csv`: direct composition comparison against the previous tournament.
- `BUDGET_SENSITIVITY.csv`: 24-vs-48 stage comparison.
- `SPEED_SWEEP_24_STAGES.csv`, `SPEED_SWEEP_48_STAGES.csv`: all grid points.
- `SELECTED_TRACES_24_STAGES.csv`: selected stress checkpoints.
- `TIMING_24_STAGES.csv`: secondary timing data.
- `RESULTS.json`: combined machine-readable summary.

## Limits

This is a manufactured reduced structural-dynamics benchmark. The rankings are not universal solver rankings. In particular, the line-search and pseudo-transient results show that forcing strict 1:1 alternation can deliberately remove the persistence that makes a globalization method effective. This is evidence about composition contracts and scheduling, not evidence that those established methods are intrinsically weak.
