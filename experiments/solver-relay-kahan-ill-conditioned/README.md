# Solver Relay: Kahan-128 Ill-Conditioned Linear System

## Target

This branch tests solver handoff on a canonical ill-conditioned linear-algebra problem rather than on an implicit-dynamics model.

Matrix family:

\[
K_n(\theta)=\operatorname{diag}(1,s,\ldots,s^{n-1})U,\qquad s=\sin\theta,\ c=\cos\theta,
\]

where `U` is upper triangular with diagonal `1` and every strict upper-triangular entry `-c`.

Locked case:

- `n = 128`
- `theta = 1.2`
- `pert = 0`
- deterministic hidden-oracle solution pattern in the reproducibility harness
- mathematical input is the generator, so arbitrary-precision reconstruction is allowed

A dense orthogonal mask is intentionally not used in the binary64 challenge. For singular values far below machine epsilon, explicitly forming a dense rotated matrix can inject rounding perturbations that alter the sub-epsilon conditioning being tested. The canonical Kahan structure is retained.

## Conditioning

The Kahan family is a standard counterexample for rank-revealing QR with column pivoting. In this locked `128 x 128` case, the experiment estimated

- `sigma_max ≈ 1.0728719264e1`
- `sigma_min ≈ 2.1537625220e-21`
- `cond_2 ≈ 4.9813845095e21`

using high-precision power/inverse-power iteration.

Double CPQR kept every column in its original order and reported all 128 diagonal entries of `R` above its usual rank threshold. Its smallest `|R_ii|` was about `1.3126e-4`, roughly `6.09e16` times larger than the high-precision estimate of `sigma_min`.

## Residual versus forward error

A small backward residual is not sufficient for this problem.

The double-LU solution had

- backward error: `5.93e-17`
- forward error: `9.68e-1`

so a residual-only acceptance condition would accept a solution whose vector is almost one full relative unit away from the oracle.

Double QR/CPQR improved the forward error only to about `6.76e-1`. Double SVD recognized one numerically lost direction (`rank = 127`) and reduced forward error to about `2.69e-2`, but still did not recover the hidden solution accurately.

## Mixed-precision handoff result

The useful transfer was not another approximate solution. It was higher-precision residual information while retaining the original double-precision factorization.

1. Double LU builds and retains its factorization.
2. A higher-precision stage recomputes the residual of the same mathematical problem.
3. That residual becomes the correction right-hand side.
4. The already-built double LU factorization solves one correction equation.

This is a standard mixed-precision iterative-refinement mechanism. The branch does not claim this algorithm as a new numerical method. Its role here is to test the Solver Relay architecture: one numerical stage can provide the specific information required by another stage without replacing the entire solver.

One correction produced:

| residual precision | forward error after correction |
|---:|---:|
| 18 dps | `5.66e-5` |
| 20 dps | `6.84e-6` |
| 22 dps | `2.41e-8` |
| 24 dps | `1.89e-10` |
| 26 dps | `3.40e-12` |
| **30 dps** | **`4.01e-16`** |

Thus the same double factorization that initially produced a `~0.97` forward error became essentially binary64-accurate after receiving a sufficiently accurate residual evaluation.

By contrast, ordinary double-precision iterative refinement did not reliably converge: after five corrections its forward error was still about `4.84e-1`.

## Solver-handoff interpretation

This branch deliberately does not force NRII into a linear-algebra role where it is unnecessary. The purpose is broader: the relay architecture remains meaningful when the participating numerical stages are unrelated algorithms.

Transferable information exposed by this benchmark includes:

- current approximate state `x`;
- reusable LU factors and pivots;
- rank diagnostics from CPQR and SVD;
- numerically unresolved singular directions;
- required residual-evaluation precision;
- higher-precision residual data;
- independent forward/backward error verification.

The experiment therefore complements the NRII/Newton branches: NRII is one possible producer of solver handoff data, not a mandatory stage in every relay graph.

## Evaluation criteria

A solver should not receive credit for residual alone. Report all of:

1. normalized backward error;
2. oracle forward error;
3. reported numerical rank at binary64 precision;
4. condition-number decade;
5. whether the solver explicitly warns that residual-only acceptance is unsafe;
6. total cost and maximum precision used.

For this locked case, the declared acceptance conditions are:

- forward error `<= 1e-12`;
- backward error `<= 1e-14`;
- report the effective binary64 rank loss (`127` under the SVD tolerance used by the harness) or an equivalent conditioning warning;
- no hidden use of the oracle solution in the solve path.

## Reproduction

Requires Python 3, NumPy, SciPy and mpmath.

```bash
python kahan_ill_conditioned_relay.py
```

Generated files:

- `BASELINES.csv`
- `DOUBLE_REFINEMENT.csv`
- `PRECISION_RELAY.csv`
- `MULTIPRECISION_DIRECT.csv`
- `MATRIX_DIAGNOSTICS.csv`
- `SUMMARY.json`
