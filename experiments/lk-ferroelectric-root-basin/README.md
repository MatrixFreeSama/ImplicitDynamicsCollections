# NRII vs Newton: Landau–Khalatnikov Root-Basin Cloud Study

## Purpose
Isolate a Newton root-basin / initial-Jacobian failure without simultaneously making the physical continuation branch singular.

## Model
Dimensionless homogeneous Landau–Khalatnikov switching law

`dP/dt = P - P^3 + E`

with `P0=-0.3`, `E=0.2`. The implicit one-step family is

`P(h)=P0+h*dt*(P(h)-P(h)^3+E)`, `h in [0,1]`.

The deliberately diagnostic timestep is

`dt_crit = 1/(1-3 P0^2) = 1.36986301369863`.

At this value the Newton Jacobian evaluated at the old-state initial guess is zero, but an independent 100000-substep continuation audit finds the connected solution branch stays nonsingular with minimum absolute branch Jacobian `0.612661208`.

## Critical point
- Newton: `singular_jacobian_at_iterate`, iteration `0`, initial residual `0.1`.
- Full NRII chain: declared p=74, Pade M=27, residual `4.35435704e-09`, candidate P `-0.551756714948`.
- Independent connected-branch reference P: `-0.551756710007`.
- NRII absolute branch error: `4.94116115e-09`.
- Chebyshev tail ratio: `2.56667459e-19`.
- Pade minimum sampled denominator magnitude: `0.0058590052`.

## Above the critical timestep
For factors above 1, direct Newton often converges successfully but to a different real root from the h-connected branch. In this sweep there are 4 such tested points.

This distinguishes root selection from simple nonlinear-solver failure: the final implicit algebraic equation has multiple valid roots, while NRII's coefficient construction plus continuation follows the branch connected to the old state.

## Scope
This is a constitutive root-basin diagnostic, not a calibrated material/device performance model. It uses the homogeneous Landau–Khalatnikov form and dimensionless parameters to isolate solver behavior. It does not claim all ferroelectric switching states have this geometry.
