# General nonlinear dynamics Newton-hybrid tournament, bare C17

## Purpose

This branch compares Newton with several globalization, acceleration, continuation, quasi-Newton, and independent-construction hybrids on the **same** 128-DOF nonlinear structural-dynamics implicit step.

This is a manufactured-solution convergence benchmark, not a reproduction of one published numerical example and not a universal ranking of nonlinear solvers.

The benchmark is ISO C17 and uses no BLAS, LAPACK, PETSc, NumPy, or SciPy.

## Literature basis

The solver families are established in the literature:

- Crisfield, *Nonlinear structural dynamics via Newton and quasi-Newton methods*, Nuclear Engineering and Design 58(3), 339-348, 1980, DOI `10.1016/0029-5493(80)90147-8`.
- Gee, Kelley, and Lehoucq, *Pseudo-transient continuation for nonlinear transient elasticity*, International Journal for Numerical Methods in Engineering 78(10), 1209-1219, 2009, DOI `10.1002/nme.2527`.
- Dallas and Pollock, *Newton-Anderson at Singular Points*, International Journal of Numerical Analysis and Modeling 20(5), 667-692, 2023, DOI `10.4208/ijnam2023-1029`.

## Dynamics problem

A unit-mass chain has cubic hardening foundation forces and nonlinear nearest-neighbor springs:

```text
f_i(q) = k q_i + alpha q_i^3 + edge contributions
edge(d) = kc d + beta d^3
```

with

```text
NDOF  = 128
k     = 0.10
alpha = 1.00
kc    = 0.20
beta  = 0.02
dt    = 1.00
q0    = 0
v0    = 0
```

Backward Euler reduces the step to

```text
R(q) = q + dt^2 f_int(q) - dt^2 f_ext = 0.
```

For every stress level a deterministic target `q* = stress * shape` is manufactured and `f_ext` is generated from the same discrete equation. Therefore the declared root is known independently.

Success requires both

```text
relative residual <= 1e-12
relative forward error to q* <= 1e-8.
```

## Methods

- Newton: full-step exact-Jacobian Newton.
- Newton + LineSearch: Newton with Armijo-type residual backtracking.
- Newton + Aitken: Aitken-relaxed Newton; accelerated candidate is accepted only if its residual improves over raw Newton.
- Newton + Anderson: depth-3 Anderson acceleration on the Newton fixed-point map; candidate is residual-safeguarded.
- Broyden -> Newton: inverse-Broyden phase followed by exact Newton terminal convergence.
- PseudoTransient -> Newton: `(J + mu I) dx = -R`, with `mu` reduced toward exact Newton.
- NRII -> Newton: order-12 NRII constructs an implicit-Euler coefficient candidate. Residual is used only as a handoff audit. It never drives NRII correction.

## Operational convergence envelope

`stress` is a dimensionless nonlinear step-severity / predictor-distance multiplier. It is not a physical constant.

The scan uses 121 logarithmically spaced values from `0.05` to `1e12` with at most 24 nonlinear iterations. Line search permits up to 80 backtracks per Newton iteration. Anderson depth is 3 and NRII order is 12.

The principal metric is the **largest contiguous successful stress starting from the easy end**. This avoids treating isolated reacquisition at a harder point as a continuous basin.

| Method | Contiguous success limit | First failed stress | Farthest isolated success |
| --- | ---: | ---: | ---: |
| Newton | 38.0943 | 49.1702 | 38.0943 |
| Newton + LineSearch | 7.74744e11 | 1.0e12 | 7.74744e11 |
| Newton + Aitken | 488.966 | 631.133 | 488.966 |
| Newton + Anderson | 2.96784 | 3.83074 | 1.0e12 |
| Broyden -> Newton | 2.29932 | 2.96784 | 1.0e12 |
| PseudoTransient -> Newton | 801114 | 1.03404e6 | 801114 |
| NRII -> Newton | 38.0943 | 49.1702 | 38.0943 |

Relative to plain Newton's contiguous envelope in this benchmark:

```text
Newton + Aitken           ~12.84 x
PseudoTransient -> Newton ~2.103e4 x
Newton + LineSearch       ~2.034e10 x
NRII -> Newton             1.00 x
```

These are operational benchmark boundaries under the locked budgets, not mathematical convergence limits.

## Speed observations

Repeated same-process timings are secondary to operation counts because microsecond timings are sensitive to container scheduling.

At `stress = 1` the recorded optimized run gives approximately:

```text
Newton                   7 iterations, 36.7 us
Newton + LineSearch      5 iterations, 28.9 us
PseudoTransient -> Newton 8 iterations, 42.4 us
NRII -> Newton            7 iterations, 119.7 us
```

At `stress = 5`:

```text
Newton                   14 iterations, 71.2 us
Newton + LineSearch       7 iterations, 42.3 us
PseudoTransient -> Newton 9 iterations, 46.8 us
Newton + Anderson         failed within 24 iterations
Broyden -> Newton         failed within 24 iterations
NRII -> Newton            14 iterations, 152.4 us
```

At low stress, the audited NRII candidate is accepted and can reduce terminal Newton iterations, for example from 3 to 2 at `stress = 0.1`. In this small CPU case the coefficient-construction overhead makes total NRII->Newton wall time slower, so this is an iteration-count reduction, not a speed victory.

The NRII handoff gate accepts candidates only through approximately `stress = 0.4972`; outside that region it falls back to the original Newton predictor.

## Interpretation

The benchmark does not produce one universal winner.

Plain Newton is cheapest in the easy regime. Line-search Newton has by far the largest contiguous convergence envelope in this monotone hardening family. Pseudo-transient continuation is the second-strongest globalization under the declared budget. Aitken extends the Newton envelope moderately. Anderson and Broyden show non-monotone success islands, demonstrating that acceleration is not automatically equivalent to globalization.

The NRII/Newton relay is intentionally retained as a negative control for universal superiority: order-12 NRII reduces terminal Newton iterations only in the low-stress region and does not enlarge the contiguous envelope in this benchmark.

## Build

```bash
gcc -O3 -march=native -std=c17 -Wall -Wextra -pedantic \
  newton_hybrid_dynamics_tournament.c -lm -o tournament
./tournament
```

The source wrapper includes `source_part1.inc`, `source_part2.inc`, and `source_part3.inc`; splitting is only for repository transport and does not change the C translation unit.

The optimized build was warning-clean. A separate AddressSanitizer + UndefinedBehaviorSanitizer run completed with exit code 0 and zero sanitizer stderr bytes.

Recorded optimized full-tournament runtime in the current cloud container was about `0.69 s`.
