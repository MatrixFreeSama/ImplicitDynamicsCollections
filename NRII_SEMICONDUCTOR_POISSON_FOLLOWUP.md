# NRII Semiconductor Poisson Follow-Up

**Date:** 2026-09-06  
**Implementation:** bare C, FP64, single-threaded CPU  
**Scope:** controlled normalized 3D drift-diffusion/Poisson benchmark, not production TCAD validation.

## Objective

The previous constitutive campaign excluded the electrostatic Poisson equation. This follow-up restores a genuine global elliptic field. Poisson is never treated as local and is never removed.

The NRII organization is:

carrier coefficient construction -> scalar Poisson solve -> next coefficient order.

The comparison asks whether this is cheaper than solving the complete potential/electron/hole correction globally.

## Controlled discretization

The carrier operator uses a 3D finite-neighborhood Scharfetter-Gummel-style Bernoulli flux. The Bernoulli function is represented by a smooth sixth-order local polynomial over the small normalized potential differences used in this campaign. Electrostatics uses a seven-point Dirichlet Poisson operator with Jacobi-preconditioned CG.

Tested model families:

- SG-style carrier transport;
- radiative recombination;
- localized heterojunction transport/recombination;
- impact ionization;
- combined heterojunction + impact ionization.

Every NRII result is checked against the complete coupled residual.

## Solver organizations

The same implicit relation is solved four ways:

1. monolithic Newton-GMRES over potential, electrons and holes;
2. reduced Newton-GMRES after Poisson elimination, where each reduced Jacobian-vector product invokes a scalar Poisson solve;
3. segregated Poisson/carrier fixed-point iteration;
4. NRII + Poisson, with one scalar Poisson closure per retained coefficient order.

## Broad sweep: 6^3 through 16^3

Measured points: **360**

Maximum residuals:

- monolithic Newton: 9.642e-10
- reduced Newton: 9.642e-10
- segregated solve: 1.420e-11
- NRII + Poisson: 1.426e-12

At 16^3:

| Model | vs monolithic Newton | vs reduced Newton | vs segregated |
|---|---:|---:|---:|
| sg_transport | 1.84x | 1.06x | 1.03x |
| radiative | 2.56x | 1.33x | 1.17x |
| heterojunction | 3.10x | 1.56x | 1.16x |
| impact_ionization | 2.37x | 1.22x | 1.11x |
| combined_hetero_impact | 3.71x | 1.50x | 1.09x |

The self-consistent Poisson field narrows the acceleration substantially relative to purely local hotspot experiments, but the relative NRII position generally improves as the 3D grid grows.

## Larger-grid focus: 18^3 through 24^3

Measured points: **30**

At 24^3:

### Heterojunction

- vs monolithic Newton-GMRES: **3.72x**
- vs reduced Newton-GMRES: **1.39x**
- vs segregated Poisson/carrier iteration: **1.09x**

### Combined heterojunction + impact ionization

- vs monolithic Newton-GMRES: **4.07x**
- vs reduced Newton-GMRES: **1.33x**
- vs segregated Poisson/carrier iteration: **1.26x**

The median q_eff at 24^3 is approximately **0.034**, with median retained NRII order **6**.

## Poisson work is compressed, not removed

At 24^3 the median cumulative Poisson-PCG iteration counts per implicit step are:

- NRII + Poisson: **597**
- reduced Newton: **995**
- segregated iteration: **984**

The physical Poisson inverse remains global. The difference is organizational:

- reduced Newton embeds Poisson solves inside Krylov directions;
- segregated iteration invokes Poisson once per outer sweep;
- NRII invokes Poisson once per retained coefficient order.

The measured mechanism is therefore:

**global Poisson retained + global carrier/block correction removed + cumulative Poisson inversion depth reduced.**

## Revised semiconductor golden region

A favorable Poisson-coupled NRII case requires:

1. low or moderate q_eff, keeping coefficient order small;
2. controlled transient amplification;
3. finite-neighborhood carrier and constitutive algebra;
4. low or moderate coefficient-algebra density;
5. a scalar Poisson solve materially cheaper than the original coupled block correction;
6. enough device size or coupling complexity for global block overhead to matter.

Poisson creates a hard performance floor. This explains why the 10-20x gains of purely local active-sparse 3D tests do not survive unchanged.

## Negative result

On small grids and the simplest SG transport model, NRII + Poisson can lose to conventional organizations even with very small q_eff.

Therefore:

**low q is not sufficient once a genuine global Poisson solve is present.**

The cost and iteration depth of Poisson closure must be included in the preflight cost model.

## Evidence boundary

This remains a controlled normalized solver-structure benchmark. It is not production TCAD.

The current carrier flux is SG-style rather than a production finite-volume semiconductor implementation. The benchmark does not yet include calibrated silicon or GaN parameters, contact models, Fermi-Dirac statistics, production AMG/ILU/geometric-multigrid preconditioning, commercial continuation strategies, or device-level IV and breakdown-voltage validation.

The next production-oriented step is a calibrated PN/PIN device with exact finite-volume Scharfetter-Gummel flux, self-consistent Poisson, physical contacts, and production-strength Newton/Gummel/preconditioned baselines.
