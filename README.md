# ImplicitDynamicsCollections

## NRII historical origin

NRII (Non-Residual-Iterative Implicit) was originally conceived as a direct, general-purpose constructive route for transcendental equations, not as an implicit-dynamics time integrator.

In that original design, the primary objective was to construct useful equation information directly. Convergence of the immediately emitted expression as a standalone final solution was not a required design target. The construction could instead be consumed by a downstream solver.

The implicit-dynamics work in this repository is a later application of the same coefficient-closure idea. Accordingly, experiments may study NRII in more than one role: as a standalone implicit construction, as a branch/trajectory information generator, or as a producer of structured input for another solver such as Newton.

## Experimental status

Experiment branches are exploratory and should not be read as universal performance or convergence claims. Each branch records its own model, numerical contract, validation gates, and comparison conditions.
