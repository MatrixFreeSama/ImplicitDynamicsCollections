# Newton information-deficiency redirection tournament, bare C17

## Purpose

This experiment keeps the same 128-DOF manufactured nonlinear structural-dynamics implicit problem and the same Newton partner families used by the preceding directional-hybrid and fixed 1:1 reciprocal-alternation tournaments. The composition rule is changed again.

Newton is now the primary solver. A separate information gate observes Newton progress and redirects to the paired specialist only when the currently observed Newton contraction is operationally insufficient under the declared remaining dispatch budget, or when Newton progress is nondecreasing. After one specialist dispatch, control returns to Newton.

The scheduler therefore has the form

```text
Newton -> Newton -> ... -> information gate -> specialist -> Newton -> ...
```

rather than unconditional 1:1 alternation.

This is an operational scheduling experiment, not a theorem about the mathematical basin of Newton or any partner solver.

## Fixed nonlinear dynamics problem

The discrete implicit equation is

```text
R(q) = q + f_int(q) - f_ext = 0
```

for the same 128-DOF cubic-hardening nearest-neighbor chain used by the previous tournaments. A manufactured root `q* = stress * shape` is used to verify the selected branch independently.

Success requires both

```text
relative residual <= 1e-12
relative error to manufactured root <= 1e-8
```

on a 121-point logarithmic stress grid from `0.05` to `1e12`.

## Information-deficiency gate

A redirect is considered only after at least two Newton stages since the previous redirect. Let

```text
rho_k = ||R_{k+1}|| / ||R_k||
```

for the latest Newton step. If `rho_k >= 1`, the step provides nondecreasing progress and the gate requests a redirect. If `0 < rho_k < 1`, the scheduler estimates the number of additional Newton steps required at the observed geometric contraction rate:

```text
n_need = log(tol / rel_residual) / log(rho_k)
```

A redirect is requested when `n_need` exceeds the number of dispatch stages remaining in the declared budget.

This criterion is intentionally budget-aware. It should be read as **operational Newton information deficiency under the declared budget**, not as proof that Newton is mathematically unable to converge given unlimited iterations.

## Solver pairings

The comparison families are unchanged:

- `Newton`: standalone baseline.
- `Newton->LineSearch->Newton`: Armijo line-search Newton specialist.
- `Newton->Aitken->Newton`: Aitken-relaxed Newton-map specialist.
- `Newton->Anderson->Newton`: depth-3 Anderson Newton-map specialist.
- `Newton->Broyden->Newton`: inverse-Broyden specialist with secant-updated inverse information accumulated along Newton states.
- `Newton->PseudoTransient->Newton`: shifted-Jacobian pseudo-transient specialist.
- `Newton->NRII->Newton`: re-anchored order-12 NRII specialist. NRII residual evaluation is accept/reject audit only and never drives an NRII correction.

A specialist dispatch counts against the same outer dispatch budget as a Newton stage. Internal residual, Jacobian, and linear-solve counts are also recorded because native specialist work is not cost-equivalent across methods.

## Primary 24-dispatch result

| Method | contiguous operational envelope | relative to Newton | redirects | successful redirects |
|---|---:|---:|---:|---:|
| Newton | 38.0942859 | 1x | 0 | 0 |
| Newton->LineSearch->Newton | 38.0942859 | 1x | 808 | 808 |
| **Newton->Aitken->Newton** | **81.9191782** | **2.15043x** | 793 | 793 |
| Newton->Anderson->Newton | 2.29931649 | 0.06036x | 815 | 813 |
| Newton->Broyden->Newton | 13.7243622 | 0.36027x | 821 | 102 |
| Newton->PseudoTransient->Newton | 38.0942859 | 1x | 808 | 808 |
| Newton->NRII->Newton | 8.23774486 | 0.21625x | 831 | 0 |

Under this gate and 24-dispatch budget, triggered Aitken redirection is the only tested pairing that enlarges the contiguous envelope beyond standalone Newton, reaching about `2.15x` the Newton frontier.

The gate is dormant in the easy regime. At stresses `0.1`, `0.5`, and `1.0`, all paired methods complete with zero redirects. At `stress=2`, Aitken uses one redirect and reduces the total dispatch count from Newton's 10 stages to 9 in the recorded timing run.

## What redirect success means

`successful_redirects` counts specialist dispatches that reduce the residual from the state handed to the specialist. It does **not** mean that the specialist beats the raw Newton candidate, improves wall time, or enlarges the convergence envelope.

This distinction is visible in the results:

- line search and pseudo-transient specialists reduce residual on every requested redirect but do not enlarge the 24-stage Newton envelope under this one-dispatch-return policy;
- Aitken both accepts every redirect and enlarges the envelope;
- Broyden accepts only about 12.4% of requested redirects and shrinks the envelope under repeated retry;
- NRII accepts none of the 831 requested redirects in the primary 24-stage sweep and therefore only consumes dispatch budget in this particular cubic-hardening family.

The NRII negative result is retained intentionally.

## Comparison with the two previous composition rules

The same solver families have now been tested under three distinct scheduling rules.

| Family | directional / terminal hybrid | fixed 1:1 alternation | information-triggered redirect |
|---|---:|---:|---:|
| Newton | 38.0943 | 38.0943 | 38.0943 |
| LineSearch | 7.74744e11 | 38.0943 | 38.0943 |
| Aitken | 488.966 | 136.480 | 81.9192 |
| Anderson | 2.96784 | 2.29932 | 2.29932 |
| Broyden | 2.29932 | 38.0943 | 13.7244 |
| PseudoTransient | 801114 | 29.5133 | 38.0943 |
| NRII | 38.0943 | 2.96784 | 8.23774 |

The data therefore reject a simple ordering such as `triggered redirect > alternation > directional hybrid`. Scheduling changes the operational convergence frontier in a strongly solver-dependent way.

In particular, line search and pseudo-transient continuation need persistence to act as globalization mechanisms. A single specialist dispatch followed immediately by raw Newton does not reproduce the large envelope of their directional forms. Aitken benefits from triggered use but reaches a smaller frontier than fixed 1:1 alternation in this benchmark.

## 48-dispatch sensitivity

The same triggered policy was rerun with 48 total dispatch stages.

| Method | 24 stages | 48 stages |
|---|---:|---:|
| Newton | 38.0943 | 4862.46 |
| LineSearch redirect | 38.0943 | 4862.46 |
| **Aitken redirect** | **81.9192** | **37462.1** |
| Anderson redirect | 2.29932 | 136.480 |
| Broyden redirect | 13.7244 | 176.162 |
| PseudoTransient redirect | 38.0943 | 4862.46 |
| NRII redirect | 8.23774 | 227.381 |

At 48 stages Aitken redirection reaches about `7.70x` the 48-stage Newton frontier. The large movement of every frontier with budget again confirms that these are **operational convergence envelopes under a declared budget**, not absolute mathematical basins.

## Failed-specialist policy sensitivity

The primary scheduler permits later retries after a specialist fails to improve the handed-off state. A secondary control disables all later redirects after the first failed specialist dispatch.

This matters for poorly matched specialists:

- Broyden improves from `13.7244` to `17.7147` when repeated failed redirects are disabled.
- NRII improves from `8.23774` to `29.5133`, despite still producing zero accepted specialist candidates.

That control shows a second scheduling principle: repeatedly redirecting to a specialist that is demonstrably unproductive can itself destroy the remaining Newton budget.

## Main architectural finding

The benchmark supports a narrower scheduling statement:

> Solver redirection is useful only when the trigger signal and the specialist's native contract are compatible. Detecting that the primary solver is operationally information-poor is not enough; the scheduler must also choose a specialist capable of converting the available handoff data into a better state.

The results also separate three questions that should not be collapsed:

1. Did Newton trigger a redirect?
2. Did the specialist reduce the handed-off residual?
3. Did the overall solver pair enlarge the contiguous operational convergence envelope?

A positive answer to one does not imply a positive answer to the others.

## Reproducibility

Primary build:

```bash
gcc -O3 -march=native -std=c17 -Wall -Wextra -pedantic \
  newton_information_deficiency_redirect_tournament.c -lm -o tournament
./tournament
```

48-stage sweep:

```bash
gcc -O3 -march=native -std=c17 -DMAX_STAGES=48 -DNO_TIMING \
  newton_information_deficiency_redirect_tournament.c -lm -o tournament48
./tournament48
```

Failed-specialist circuit-breaker sensitivity:

```bash
gcc -O3 -march=native -std=c17 -DDISABLE_AFTER_SPECIALIST_FAILURE=1 -DNO_TIMING \
  newton_information_deficiency_redirect_tournament.c -lm -o tournament_faildisable
./tournament_faildisable
```

The warning-clean optimized build emits zero compiler-warning bytes. A 24-stage AddressSanitizer + UndefinedBehaviorSanitizer sweep completed with exit code 0 and zero sanitizer stderr bytes.

## Files

The GitHub branch keeps the complete compilable C17 source as a small wrapper plus three source include parts, together with the compact reproducibility summaries:

- `CONVERGENCE_ENVELOPE_24_STAGES.csv`: primary convergence frontiers.
- `CONVERGENCE_ENVELOPE_48_STAGES.csv`: 48-stage budget sensitivity.
- `COMPOSITION_COMPARISON_24.csv`: directional vs fixed-alternation vs triggered-redirection comparison.
- `GATE_DIAGNOSTICS.csv`: trigger reasons and redirect acceptance statistics.
- `FAILED_SPECIALIST_POLICY_SENSITIVITY.csv`: repeated-retry vs fail-disable control.
- `BUDGET_SENSITIVITY.csv`: 24-vs-48 triggered-redirection comparison.
- `METHODS.csv`, `MODEL_PARAMETERS_24_STAGES.csv`: machine-readable contracts and parameters.

The downloadable experiment archive additionally contains the complete monolithic C17 source, full 121-point 24/48-stage sweeps, repeated timing data, sanitizer records, and the combined JSON result. All omitted large tables are regenerated directly by the committed source.

## Limits

This is a manufactured reduced nonlinear structural-dynamics benchmark. The gate is deliberately simple and uses only the most recent Newton contraction rate plus the declared remaining dispatch budget. It is not an optimal scheduler, not a learned policy, and not a universal information metric. Rankings should not be generalized to unrelated mechanics, multiphysics, nonsmooth events, singular systems, or different cost models without separate tests.
