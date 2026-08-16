# NRII <-> Newton reciprocal avalanche experiment

## Purpose

This experiment asks whether NRII and Newton can be used as reciprocal information producers rather than forcing either method to solve the entire nonlinear step alone.

The core cycle is:

1. NRII coefficient closure constructs a branch-informed partial continuation state.
2. Newton performs a short correction burst toward the final implicit root.
3. If the root is not yet accepted, Newton's corrected state is fed back as the next NRII anchor.
4. The cycle repeats.

This is an **external hybrid solver**. It is not a redefinition of standalone NRII. The original NRII coefficient recursion remains non-residual-iterative. Residual information is allowed only in the outer NRII-Newton controller when it compares Newton-corrected candidates.

## Same-target re-anchoring

A naive replacement of NRII C0 would change the implicit problem. This experiment therefore uses a re-anchored homotopy that preserves the original final equation.

For the carrier state U and Newton-produced anchor S:

    U(h) = S + h[(U_old - S) + dt F(U(h))]

At h=0 the construction starts from S. At h=1:

    U = U_old + dt F(U)

which is exactly the original implicit target.

Consequently:

    C1 = U_old - S + dt F_0
    C{k+1} = dt F_k, k >= 1

For Poisson, the corresponding anchor residual is removed continuously:

    R_phi(U(h)) = (1-h) R_phi(S)

so the h=1 endpoint again satisfies the original Poisson equation.

## Partial NRII step

The experiment does not force the Padé continuation all the way to h=1 on every rung. A declared fraction theta is used by scaling the coefficient sequence for U(theta t), t in [0,1], and then applying the same Chebyshev/Padé construction on that scaled interval.

This matters in the avalanche model. At drive 88, the full h=1 one-shot predictor had an original-equation residual of about 8.73e5 and Newton failed. A partial continuation can land in a useful Newton basin instead.

## Model

The model is the same dimensionless 1D PN drift-diffusion-Poisson avalanche stress test used in the preceding hybrid branch.

- nominal global grid: 2^34 cells
- locally materialized window: 1537 cells
- dt = 0.005
- final residual tolerance = 1e-10
- impact ionization: smooth reduced Chynoweth-type field law
- Newton: analytic 3x3 block-tridiagonal Jacobian, block Thomas solve, positivity cap, Armijo backtracking
- NRII orders: 16, 24, 32, 48, 64, 80, 96, 112, 128, 160, 192

## Experiment A: locked alternating policy

A calibration sweep at drive 94 was used only to choose a fixed policy. The locked policy was then:

- theta = 0.7
- two Newton updates per rung
- at most 16 rungs

The locked sweep shows a qualitatively different behavior from the one-shot hybrid. Pure old-state Newton fails at drive 78. The previous full-h one-shot NRII->Newton hybrid passes through drive 87 and fails at drive 88. The fixed alternating policy passes every point that was explicitly included in the locked sweep through drive 107. It fails at drive 108, passes drive 109, and fails at drive 110, showing that Padé/continuation safety is not monotone in the drive parameter.

Representative cases:

- drive 88: pure Newton fails; one-shot hybrid fails; fixed alternating passes in 9 rungs / 17 Newton burst updates, residual 3.13e-13.
- drive 94: pure Newton fails; one-shot hybrid fails; fixed alternating passes in 6 rungs / 11 Newton burst updates, residual 1.07e-13.
- drive 104: pure Newton fails; one-shot hybrid has no declared safe full-step seed; fixed alternating passes in 4 rungs / 7 Newton burst updates, residual 1.90e-13.
- drive 107: fixed alternating passes in 4 rungs / 7 Newton burst updates, residual 1.02e-12.

## Experiment B: reciprocal feedback controller

The stronger exploratory controller declares theta candidates:

    0.7, 0.6, 0.5, 0.4, 0.3

For each rung, every structurally safe NRII partial predictor is given a two-update Newton burst. If none reaches the final tolerance, the Newton-corrected candidate with the smallest post-burst residual becomes the next NRII anchor. This is the explicit "Newton feeds NRII back" step.

Again, this residual selection belongs to the outer hybrid and is not part of NRII coefficient generation.

Selected stress results:

- drive 108: pass, 10 rungs, 46 candidate trials, 92 Newton burst updates in all trials, residual 1.99e-13.
- drive 120: pass, 6 rungs, 24 trials, residual 4.38e-13.
- drive 140: pass, 8 rungs, 29 trials, residual 4.21e-13.
- drive 160: pass, 6 rungs, 20 trials, residual 7.98e-13.
- drive 180: pass, 5 rungs, 15 trials, residual 8.28e-13.
- drive 190: pass, 6 rungs, 21 trials, residual 4.04e-12.
- drive 195: fail under the declared order/theta set; first accepted feedback rung remains at about 5.54e11 residual.
- drive 198: no structurally safe declared candidate was found.
- drive 200: fail under the declared policy.

Thus the reciprocal construction pushes the tested old-state root-basin stress far beyond both pure old-state Newton and the prior one-shot hybrid, but it still has a finite domain.

## Important control

Standard Newton parameter continuation is a strong alternative. In a separate control, Newton was given the previous converged drive solution as the next initial guess. That continuation solved the tested ladder through drive 800 and failed at drive 1000.

Therefore this experiment does **not** establish that NRII is uniquely capable of branch continuation. Its narrower claim is that, starting each stress case from the original old state rather than a previous-drive root oracle, NRII and Newton can exchange intermediate information repeatedly and recover roots that neither old-state Newton nor the full-step one-shot hybrid reached under their locked policies.

## Timing note

The cloud implementation is a dense CPU reference and repeatedly rebuilds coefficient banks, Chebyshev transforms, and Padé fits. Reciprocal feedback is therefore much more expensive than the local Newton corrector in wall time. The present experiment is a basin/robustness study, not a performance claim. GPU conclusions require the optimized matrix-free production path.

## NRII historical context

The repository main README records the historical design intent: NRII was originally conceived as a direct, general-purpose constructive route for transcendental equations, not specifically as an implicit-dynamics time integrator. The immediate emitted construction was not required to be a converged standalone final solution; it could instead supply structured information to a downstream solver. This branch deliberately explores that original "information producer" role.
