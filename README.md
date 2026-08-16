# ImplicitDynamicsCollections

## Experimental branch: NRII Landau-Khalatnikov root-basin study

Branch: `experiment/nrii-lk-root-basin-cloud`

This branch isolates a root-basin / root-selection difference between a Newton solve started from the old state and the NRII coefficient-construction path on a dimensionless homogeneous Landau-Khalatnikov ferroelectric switching model.

### Critical diagnostic point

Model:

`dP/dt = P - P^3 + E`

with `P0=-0.3`, `E=0.2`, and

`dt_crit = 1/(1-3 P0^2) = 1.36986301369863`.

At `dt_crit` the Newton Jacobian evaluated at the old-state initial guess is approximately zero, so the tested Newton method fails immediately with `singular_jacobian_at_iterate`. An independent 100000-substep continuation audit shows that the connected physical branch itself remains nonsingular, with minimum absolute branch Jacobian about `0.612661208`.

The full NRII chain used in this experiment is:

`NRII coefficient closure -> Chebyshev audit -> Pade continuation -> original implicit residual gate`.

At the critical point the declared NRII configuration `p=74`, Padé `M=27` reaches residual `4.35435704e-09`, with candidate `P=-0.551756714948`, matching the connected-branch reference `P=-0.551756710007` to about `4.94e-09` absolute error.

Above the critical timestep, the tested Newton method can again converge but may converge to a different real root than the branch connected to the old state. The experiment therefore distinguishes simple nonlinear-solver convergence from branch ancestry / root selection.

### Scope

This is a constitutive solver diagnostic, not a calibrated ferroelectric material model. The parameters are dimensionless and chosen to isolate solver behavior. It does not claim that all Landau-Khalatnikov switching states have the same geometry.

Experiment files are under `experiments/lk-ferroelectric-root-basin/`.
