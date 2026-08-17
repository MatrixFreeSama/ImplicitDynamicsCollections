# Bare-C constrained-dynamics singular-KKT continuation benchmark

## Status

Literature-grounded benchmark adaptation. This is **not** a bitwise reproduction of the numerical setup or source code of any cited paper.

The implementation is ISO C17 and uses no external numerical library. Dense diagnostic operations, Gaussian elimination, symmetric Jacobi eigendecomposition, numerical-rank estimation, and the pseudoinverse specialist control are implemented directly in C.

Build used in the recorded cloud run:

```bash
gcc -O3 -march=native -std=c17 nrii_constrained_dynamics_singular_kkt_continuation.c -lm -o nrii_constrained_dynamics_singular_kkt_continuation
./nrii_constrained_dynamics_singular_kkt_continuation
```

A separate AddressSanitizer + UndefinedBehaviorSanitizer run completed without a reported sanitizer error.

## Literature basis

The benchmark targets a documented failure class in constrained multibody dynamics.

- S. K. Ider and F. M. L. Amirouche, “Numerical stability of the constraints near singular positions in the dynamics of multibody systems,” *Computers & Structures* 33(1), 129–137, 1989, DOI `10.1016/0045-7949(89)90135-1`. The paper states that instantaneous singular configurations can make constraint equations linearly dependent and produce singularities in conventional simulation algorithms, and demonstrates methods on a three-link planar manipulator.
- I. M. Khan and K. S. Anderson, “Performance investigation and constraint stabilization approach for the orthogonal complement-based divide-and-conquer algorithm,” *Mechanism and Machine Theory* 67, 111–121, 2013, DOI `10.1016/j.mechmachtheory.2013.04.009`. The paper states that geometric singularities can make system matrices rank deficient or ill-conditioned and can rapidly lead to simulation failure.
- N. C. Parida and S. Raha, “Regularized numerical integration of multibody dynamics with the generalized alpha method,” *Applied Mathematics and Computation* 215(3), 1224–1243, 2009, DOI `10.1016/j.amc.2009.06.063`. The paper studies ill-conditioned discretization Jacobians in singular/high-index multibody DAEs and regularization that increases the smallest singular values.

The present benchmark reproduces the **pathology class**, not the exact published mechanism parameters.

## Model

A planar 3R manipulator is the singularity-producing core. It is embedded in 24 generalized coordinates and mixed by a fixed dense orthogonal basis. Eight auxiliary independent constraints are added to the two end-effector Cartesian constraints.

- generalized coordinates: 24
- constraints: 10
- dense KKT dimension: 34 x 34
- auxiliary generalized masses: `1e-2` to `1e2`
- step size: `0.0025 s`
- singular step: 200
- singular time: `0.5 s`
- target step: 400
- target time: `1.0 s`
- NRII order: 16

At the collinear configuration `q2 = q3 = 0`, the 2 x 3 end-effector Jacobian loses one row-space dimension.

The acceleration-level matrix-dependent formulation is

```text
[ M  G^T ] [a     ] = [M a*]
[ G   0  ] [lambda]   [G a*]
```

The right-hand side is deliberately consistent with a smooth physical acceleration `a*`. Thus the singular KKT can lose multiplier uniqueness without making the physical acceleration undefined.

## Matrix-dependent baseline contract

The baseline represents a conventional Newton/KKT correction route whose linearized correction matrix must be numerically full rank. Numerical rank uses

```text
tol = DBL_EPSILON * N * sigma_max
```

where the singular values of the symmetric KKT matrix are the absolute values of its eigenvalues.

When rank drops below 34, this baseline terminates rather than treating an arbitrary floating-point LU output as a unique Newton correction.

## Matrix-Free NRII path

The smooth physical state is `y = [u,v]` with uncoupled oscillator operator action

```text
A[u,v] = [v, -omega^2 u].
```

For implicit midpoint,

```text
y1 = y0 + dt A (y0 + y1)/2.
```

NRII constructs

```text
C1 = dt A y0
Ck = (dt/2) A C{k-1},  k >= 2
```

and sums through fixed order 16. No global KKT matrix is assembled or solved on the NRII time-advance path. The implicit-equation residual is audited after construction and never creates a correction direction inside NRII.

An independent component-wise closed-form implicit-midpoint solve is used as the state oracle.

## Recorded result

The full-rank KKT baseline commits through step 199 and terminates when attempting step 200:

```text
last committed time = 0.4975 s
attempted singular time = 0.5 s
KKT numerical rank = 33 / 34
sigma_min ~= 6.64e-17
cond_2 ~= 1.51e18
```

The Matrix-Free NRII path commits all 400 steps and reaches `t = 1.0 s`.

```text
max NRII implicit-midpoint residual (inf norm) ~= 2.56e-16
max absolute state error vs closed-form midpoint oracle ~= 1.94e-14
```

At step 200, the task Jacobian numerical rank is 1 and its smallest singular value is approximately `9.37e-16` on the actual computed trajectory.

A raw Gaussian elimination solve can still return a vector at the rank-deficient point with a tiny backward error. This is not interpreted as a unique solution certificate. The pseudoinverse specialist control remains included and recovers the physical acceleration accurately, showing that the result is a **formulation-bypass result**, not a claim that all established singular-system methods fail.

## Singularity sweep

`SINGULARITY_SWEEP.csv` varies

```text
q2 = 1e-1, 1e-2, ..., 1e-12, 0
q3 = 0.6 q2
```

The KKT condition number grows from approximately `6.63e3` to the `1e18` range and numerical rank drops from 34 to 33 before the exact collinear endpoint.

## Interpretation

Supported claim:

> In this constrained-dynamics benchmark, a full-rank Newton/KKT continuation contract terminates at the singular configuration, while a Matrix-Free NRII construction that does not depend on that KKT solve continues the same independently verified smooth discrete state to the declared target time.

Not supported:

- that all conventional multibody methods fail at singular configurations;
- that pseudoinverse, null-space, regularized, or projection formulations cannot continue;
- that every physically singular dynamics problem has a uniquely continuable state;
- that NRII universally resolves rank loss.

## Versioned files

- `nrii_constrained_dynamics_singular_kkt_continuation.c`: complete C17 reproducer
- `RESULTS.json`: aggregate result and singular snapshot
- `SINGULARITY_SWEEP.csv`: conditioning/rank sweep
- `MODEL_PARAMETERS.csv`: locked parameters

Running the reproducer also generates `CONTINUATION_TRACE.csv` containing all 400 time steps. The recorded trace is included in the downloadable benchmark bundle accompanying this run rather than duplicated in the Git branch.
