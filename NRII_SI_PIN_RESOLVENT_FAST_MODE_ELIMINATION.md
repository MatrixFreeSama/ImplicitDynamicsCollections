# Resolvent Fast-Mode Elimination for Physical-Unit Si PIN NRII

**Date:** 2026-09-07  
**Implementation:** bare C, FP64, single-threaded CPU  
**Device:** the same physical-unit 24 um Si PIN benchmark used in `NRII_SI_PIN_PHYSICAL_DEVICE_FOLLOWUP.md`

## 1. Motivation

The previous physical-device study identified a full-domain NRII radius controlled by the few-femtosecond charge-relaxation mode of the heavily doped contacts:

q_eff approximately dt / 4.94e-15 s.

This follow-up tests whether that restriction belongs to the implicit device equation itself or to the bare power-series representation.

The result is that the dominant restrictions are representation-level fast linear spectra. They can be removed by an exactly equivalent resolvent transformation.

## 2. Exact transformation

For the backward-Euler implicit relation

Y - X = h G(Y),

introduce any fixed linear operator A and write

(I - h A)(Y - X) = h [G(Y) - A(Y - X)].

This equation is algebraically identical to the original implicit relation whenever I-hA is invertible.

An auxiliary continuation parameter lambda is then introduced:

Y(lambda) - X =
lambda h (I - h A)^(-1)
[G(Y(lambda)) - A(Y(lambda)-X)].

NRII coefficients are generated in lambda and evaluated at lambda=1.

For a purely linear mode G(Y)=A Y, choosing the matching A removes that mode from the coefficient depth exactly: the solution is obtained through the resolvent rather than a geometric h-series.

## 3. Contact-only resolvent

The first resolvent contains only a diagonal contact charge-relaxation model derived from physical conductivity:

tau_q = epsilon / sigma,

with the tested default using the physical matching factor theta=1.

This transformation preserves the original implicit equation. It does not freeze contact densities and does not replace the complete residual.

Measured N=121 result:

- bare NRII maximum validated dt in the original sweep: **3.206e-15 s**
- contact-only resolvent maximum validated dt in the refined sweep: **6.876e-13 s**
- measured radius increase: **215x**

The contact-only resolvent therefore removes the previously identified femtosecond charge-relaxation limitation, but a second grid/transport-related boundary appears at larger dt.

## 4. Transport + contact resolvent

The second operator A adds the fixed-potential linearized Scharfetter-Gummel carrier-transport Jacobian to the physical contact relaxation term.

Each coefficient order now requires:

- one tridiagonal electron resolvent solve;
- one tridiagonal hole resolvent solve;
- one Poisson coefficient closure.

The transformation remains algebraically equivalent to the original backward-Euler device equation. The tridiagonal transport resolvent changes the representation, not the target root.

For N=121:

- dt = 1e-9 s remains valid;
- q_eff = **0.144**
- retained order = **11**
- full coupled residual = **9.210e-10**
- relative difference from direct full Newton = **4.620e-13**

Relative to the original measured full-domain NRII limit, this is already a lower-bound usable-step expansion of at least **3.12e+05x**. The upper boundary was not reached in the 1e-9 s sweep.

The same 1e-9 s test remains valid for N=61, 241 and 481. At that step:

| Nodes | q_eff | order | full residual | relative difference from Newton |
|---:|---:|---:|---:|---:|
| 61 | 0.143 | 11 | 7.024e-10 | 4.386e-13 |
| 121 | 0.144 | 11 | 9.210e-10 | 4.620e-13 |
| 241 | 0.133 | 10 | 3.883e-09 | 1.420e-12 |
| 481 | 0.089 | 9 | 1.184e-09 | 3.284e-13 |

The coefficient depth is therefore bounded over a step range spanning six orders of magnitude beyond the original femtosecond scale.

## 5. Physical matching factor

The contact relaxation operator was tested with theta = 0, 0.5, 1.0 and 1.5.

The theta=1 case, which uses the physical epsilon/sigma charge-relaxation scale without an empirical multiplier, gives the strongest tested contraction. Under-resolving or over-resolving that fast mode increases q and the required coefficient order.

This provides evidence that the preconditioner is capturing a physical fast spectrum rather than acting as an arbitrary numerical damping factor.

## 6. Interpretation

The physical PIN study now separates three notions that were previously mixed together.

### 6.1 Bare NRII radius

The bare h-series directly sees every fast eigenmode. Heavily doped contact charge relaxation therefore forces a few-femtosecond radius.

### 6.2 Resolvent-transformed NRII radius

A stiff linear mode can be moved into (I-hA)^(-1). It no longer determines coefficient depth when A captures that mode.

### 6.3 Remaining nonlinear radius

After the tested contact and SG transport spectra are absorbed, the measured q_eff stays approximately 0.09-0.14 over the complete tested range through 1e-9 s. No new nonlinear branch boundary was reached in this campaign.

The statement "physical Si PIN forces NRII to femtosecond steps" is therefore false for the resolvent representation. It is true only for bare full-domain power-series NRII.

## 7. Runtime result remains negative in one dimension

The new representation solves the convergence-radius problem, but it does not beat the structured one-dimensional Newton baseline.

Direct Newton uses an O(N) block-tridiagonal factorization. The transport-resolvent NRII also performs two tridiagonal linear solves per coefficient order, plus coefficient algebra and Poisson closure.

Consequently the measured speed ratios remain below one in this one-dimensional benchmark.

This distinction is important:

**representation success does not imply performance success.**

The present result establishes that fast-mode elimination can preserve exactness and enlarge the NRII operating radius. It does not establish a 1D speed advantage.

## 8. Consequence for three-dimensional design

The transport resolvent introduces a new tradeoff. In one dimension, tridiagonal resolvent solves are cheap. In three dimensions, an exact global transport resolvent can itself become a global solve and erase the locality advantage that motivated NRII.

The next design problem is therefore not convergence theory. It is the construction of a resolvent that is both spectrally effective and locality-preserving, for example:

- element/block-local fast-mode resolvents;
- line or directional factorizations;
- local Schwarz resolvents;
- Chebyshev approximations to the fast inverse;
- sparse approximate inverse / short polynomial resolvents;
- relay of only the stiff subdomains to a specialist solver.

The exact 1D resolvent result supplies a reference target for evaluating those local approximations.

## 9. Main conclusion

The experiment establishes a hierarchy:

bare NRII
-> physical contact resolvent
-> contact + SG transport resolvent.

The first transition removes the contact dielectric-relaxation radius. The second removes the remaining linearized transport stiffness over the tested range.

The core equation is not changed. The improvement comes from changing the representation of the implicit branch.

A practical future NRII preflight should therefore estimate not only q and Gamma, but also identify which part of the predicted q is generated by removable linear fast spectra. Those spectra should be sent to a resolvent or relay layer before coefficient construction begins.
