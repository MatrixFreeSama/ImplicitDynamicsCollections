# Bare-C reduced thermo-electromagneto-hydrodynamic solver-relay benchmark

## Status

This is a **literature-grounded reduced multiphysics stress benchmark** implemented in ISO C17 with no external numerical library. It is **not** a bitwise reproduction of the geometry, material table, mesh, industrial furnace, or numerical code of the cited paper.

The benchmark is motivated by the induction-furnace model of Bermúdez, Gómez, Muñiz, Salgado, and Vázquez, *Numerical simulation of a thermo-electromagneto-hydrodynamic problem in an induction heating furnace*, Applied Numerical Mathematics 59(9), 2082–2104 (2009), DOI `10.1016/j.apnum.2008.12.005`. That paper studies an axisymmetric coupled problem containing electromagnetics, heat transfer, hydrodynamics, and phase change. Its reported numerical strategy combines BEM/FEM for the electromagnetic model, the method of characteristics for flow, and a fixed-point-like algorithm for the overall coupling.

The present code preserves the **coupled pathology class** while deliberately using a different reduced discretization and a much broader solver set in order to test explicit solver handoff.

## Purpose

The benchmark asks a different question from a monolithic multiphysics code:

> Can a coupled simulation continue when the solver currently assigned to one subproblem loses its declared numerical contract, without forcing every other subproblem to change solver at the same time?

Five stress cases are used. Each activates a different failure surface while the physical coupling remains active.

## Reduced coupled model

The computational domain is a `48 x 64` structured cross-section with 3072 cells.

### 1. Electrical coil state

A reduced RL coil state is advanced by implicit midpoint,

```text
I1 = I0 + dt/L [ V - R(Tbar) (I0 + I1)/2 ].
```

NRII constructs the implicit candidate with a fixed coefficient cascade of order 64. The implicit residual is evaluated only after construction and does not create a correction direction.

### 2. Time-harmonic eddy-current surrogate

A complex magnetic-potential surrogate `A = Ar + i Ai` is represented as two coupled real fields. The matrix-free finite-difference operator contains spatial diffusion, temperature-dependent conductivity, and the skew real/imaginary eddy-current coupling.

The source amplitude depends on the current coil state and the coupled temperature field.

From the solved field, the benchmark forms reduced Joule heating and Lorentz-force densities,

```text
Q_J ~ sigma(T) |A|^2
F_L ~ sigma(T) |A|^2.
```

These are intentionally reduced constitutive expressions, not a claim to reproduce the exact industrial material law of the paper.

### 3. Nonlinear convective-flow surrogate

The reduced momentum field `u` satisfies a nonlinear diffusion/convection problem,

```text
u - nu0 Laplacian(u) + gamma u d_x u = F_L + F_buoyancy(T).
```

It is used to preserve the key numerical ingredients needed for solver stress: diffusion, convection, electromagnetic forcing, and thermal feedback.

### 4. Thermal diffusion, convection, Joule heating, and phase change

Temperature uses an enthalpy form,

```text
H(T1) - H(T0)
  - dt [ k Laplacian(T1) - u d_x T1 + Q_J ] = 0.
```

The liquid fraction is piecewise linear between the solidus and liquidus. Therefore the ordinary derivative contract is undefined exactly at the phase corner, while a generalized semismooth derivative remains available.

### 5. Global multiphysics coupling

Temperature changes conductivity and electrical loading. Electromagnetics produces Joule heating and Lorentz forcing. Flow advects heat. Phase state changes thermal response. The entire map is iterated until the relative coupling tolerance is met.

## Solver set

The implementation contains ten distinct solver roles:

| Role | Primary | Handoff target |
|---|---|---|
| RL coil implicit step | NRII | residual audit only |
| eddy-current linear system | restarted GMRES | BiCGSTAB |
| nonlinear flow | Picard | Newton-Krylov |
| thermal nonlinear system | ordinary Newton | semismooth active-set Newton |
| global coupling | fixed point | Aitken |
| accelerated global coupling | Aitken | Anderson acceleration |

GMRES, BiCGSTAB, Newton-Krylov, the nonlinear residuals, the small Anderson least-squares system, and all matrix-free operators are implemented directly in C. No BLAS, LAPACK, PETSc, Trilinos, SciPy, or NumPy is used by the benchmark executable.

## Recorded cloud run

Build:

```bash
gcc -O3 -march=native -std=c17 nrii_induction_multiphysics_relay.c -lm -o bench
./bench
```

Recorded environment:

```text
gcc 14.2.0
x86_64
AMD EPYC 9V74
5 visible CPU cores
```

The optimized run completed the five cases in approximately `1.29 s` in the recorded container run.

AddressSanitizer + UndefinedBehaviorSanitizer completed with an empty sanitizer error stream.

## Results

| Case | Terminal solver path | Outer iterations | Handoffs | Result |
|---|---|---:|---:|---|
| weak coupling | GMRES + Picard + ordinary Newton + fixed point | 11 | 0 | completed |
| skin-effect contrast | GMRES -> BiCGSTAB once | 11 | 1 | completed |
| strong convection | Picard -> Newton-Krylov | 11 | 10 | completed |
| phase boundary | ordinary Newton -> semismooth active set | 17 | 17 | completed |
| strong global coupling | fixed point -> Aitken -> Anderson | 6 Anderson iterations after escalation | 2 | completed |

Total recorded handoffs across the five cases: `30`.

### Weak-coupling control

The control case completes with no handoff. This demonstrates that the scheduler does not switch solvers merely because relay capability exists.

### Electromagnetic contrast case

Restarted GMRES exhausts its declared operational iteration budget once. The same matrix-free electromagnetic operator and right-hand side are handed to BiCGSTAB, which converges. No thermal or flow solver is replaced.

### Strong-convection case

Picard fails its nonlinear convergence contract in 10 of the 11 coupling evaluations. Newton-Krylov accepts all 10 handoffs. The electromagnetic and thermal solvers remain unchanged.

### Phase-boundary case

The stress case deliberately places an interior temperature exactly on the solidus. Ordinary Newton therefore encounters an undefined classical derivative contract in every coupling evaluation. The semismooth active-set formulation accepts all 17 handoffs and the coupled solve reaches its target tolerance.

This is a solver-domain transition, not a claim that Newton is generally incapable of phase-change problems.

### Strong-global-coupling case

The plain outer fixed-point method exhausts its declared operational budget. Aitken relaxation also exhausts the budget assigned to it in this stress case. Anderson acceleration then reaches the coupling tolerance.

The fixed-point and Aitken events in this case are **operational-budget failures**, not proofs of mathematical divergence. Their purpose is to exercise an explicit global-coupling escalation policy under a common coupled model.

## Numerical audits from the recorded run

The final reported NRII implicit residuals are in the `4e-17` to `1.4e-16` range.

Final electromagnetic relative residuals are approximately `4.6e-8` to `5.0e-8`, consistent with the declared `5e-8` Krylov tolerance.

Final nonlinear-flow infinity residuals are approximately `4e-9` to `1e-7` depending on case.

Final thermal infinity residuals are approximately `3e-13` to `1.8e-11`.

These values audit the numerical contracts used in this reduced benchmark. They do not establish discretization accuracy against the full industrial furnace experiment.

## Interpretation

Supported conclusion:

> In one coupled thermo-electromagneto-hydrodynamic reduced benchmark, different solver families can lose their local numerical contracts in different physical regimes while the global computation continues by transferring the same declared subproblem data to a compatible fallback solver.

The central result is therefore not that one universal solver dominates all fields. It is nearly the opposite: the benchmark intentionally preserves solver heterogeneity.

```text
NRII
GMRES -> BiCGSTAB
Picard -> Newton-Krylov
ordinary Newton -> semismooth active set
fixed point -> Aitken -> Anderson
```

The failure of one stage does not automatically invalidate the state of every other stage.

## Limitations

This benchmark does **not** support the following claims:

- that it reproduces the industrial furnace results of Bermúdez et al. numerically;
- that its reduced momentum equation is a complete Navier-Stokes solver;
- that its electromagnetic finite-difference surrogate is equivalent to the paper's BEM/FEM formulation;
- that the chosen solver escalation order is universally optimal;
- that NRII solves phase-change nonsmoothness;
- that Anderson acceleration always succeeds when fixed point or Aitken fails;
- that a successful coupled algebraic solve proves spatial or temporal convergence of the underlying industrial PDE model.

The next fidelity upgrade would be a full axisymmetric eddy-current + incompressible-flow + enthalpy discretization with published furnace geometry/material data and mesh-convergence studies.

## Files

- `nrii_induction_multiphysics_relay.c` — C17 translation unit that includes the three source parts
- `nrii_multiphysics_part1.inc`, `nrii_multiphysics_part2.inc`, `nrii_multiphysics_part3.inc` — complete benchmark implementation
- `CASE_RESULTS.csv` — terminal result of each stress case
- `SOLVER_COUNTS.csv` — solver invocation/failure counts
- `RELAY_EVENTS.csv` — explicit handoff events and interpretations
- `STRESS_CASES.csv` — locked stress parameters
- `MODEL_PARAMETERS.csv` — grid, solver budgets, tolerances, and phase parameters
- `RESULTS.json` — aggregate summary
- `RUN_STDOUT.txt` — recorded optimized run output
- `RUN_TIME.txt` — recorded optimized wall time
- `ENVIRONMENT.txt` — compiler and CPU environment
- `BUILD_AND_RUN.txt` — reproduction commands
- `SANITIZER_STDERR.txt` — sanitizer error stream, empty in the recorded run
