# ImplicitDynamicsCollections

## NRII historical origin

NRII (Non-Residual-Iterative Implicit) was originally conceived as a direct, general-purpose constructive route for transcendental equations, not as an implicit-dynamics time integrator.

In that original design, the primary objective was to construct useful equation information directly. Convergence of the immediately emitted expression as a standalone final solution was not a required design target. The construction could instead be consumed by a downstream solver.

The implicit-dynamics work in this repository is a later application of the same coefficient-closure idea. Accordingly, experiments may study NRII in more than one role: as a standalone implicit construction, as a branch or trajectory information generator, or as a producer of structured input for another solver such as Newton.

## Universal Solver Relay Interface (USRI)

The **Universal Solver Relay Interface (USRI)** is an interface-level architecture for composing heterogeneous numerical solvers without requiring their internal update laws, convergence criteria, numerical representations, precision models, or implementation strategies to be unified.

A solver participating in USRI does not need to produce the final solution. It may instead produce structured solver-handoff data that another solver can consume.

A relay is written abstractly as

\[
S_i \xrightarrow{\mathcal{H}_{i\rightarrow j}} S_j,
\]

where \(S_i\) and \(S_j\) are independently defined solver cores and \(\mathcal{H}_{i\rightarrow j}\) is the structured handoff data produced by \(S_i\) and consumed by \(S_j\).

The defining principle is:

> **Solver compatibility is determined by transferable information at solver boundaries, not by similarity of internal update laws.**

USRI is intended to support heterogeneous solver composition without solver homogenization. A residual-driven solver, a non-residual constructive solver, an interval method, a continuation method, a factorization-based method, a matrix-free method, or a mixed-precision module may remain internally unchanged while participating through a declared handoff contract.

### Solver-handoff data

Solver-handoff data is a structured information package passed between solver cores. It is not restricted to a candidate solution vector. Depending on the producer and consumer, it may contain one or more of the following:

- **State data**: candidate state, corrected state, anchor, or initial point.
- **Domain data**: interval, trust region, candidate box, admissible neighborhood, or convergence region.
- **Branch data**: branch identity, ancestry, mode identity, or other information required to preserve solution selection.
- **Direction data**: tangent, continuation direction, descent direction, local trajectory direction, or other directional information.
- **Local-geometry data**: Jacobian information, approximate inverse, local metric, curvature information, or other local geometric structure.
- **Subspace data**: Krylov space, coarse space, low-rank space, null-space estimate, or other reusable search subspace.
- **Factorization data**: LU, QR, RRQR, preconditioner, or other reusable decomposition or operator structure.
- **Residual data**: residual or defect evaluated by a solver or precision regime better suited to producing it.
- **Precision data**: precision requirements, error budgets, conditioning estimates, or information about which numerical components require increased precision.
- **Verification data**: interval inclusion information, uniqueness information, admissibility conditions, forward/backward error information, or other numerical or rigorous verification results.

The handoff schema may be extended when a solver pair requires information not covered by these categories. USRI defines the interface principle rather than a closed list of solver-specific payloads.

### Relay modes

USRI permits several composition modes.

**One-way relay**

\[
S_A \rightarrow S_B.
\]

One solver constructs information needed by another solver that performs the next stage or terminal convergence.

**Reciprocal relay**

\[
S_A \rightarrow S_B \rightarrow S_A.
\]

The output of the second solver becomes a new anchor or information source for the first.

**Alternating relay**

\[
S_A \rightarrow S_B \rightarrow S_A \rightarrow S_B \rightarrow \cdots.
\]

The solver pair repeatedly exchanges information while each solver retains its native update rule.

**Multi-solver relay graph**

\[
S_1 \rightarrow S_2 \rightarrow \cdots \rightarrow S_n,
\]

or, more generally, a directed graph in which a scheduler selects the next solver according to the currently available handoff data and the information still required by the problem state.

### Interface-level hybridization

USRI treats hybridization as an interface problem rather than a requirement to derive a single blended update formula.

Two solvers may be mathematically unrelated, internally incompatible, or based on different numerical formulations and still cooperate if one can produce information the other can consume. Their native solver cores remain independently optimizable and independently testable.

This repository refers to this as **interface-level solver hybridization** or **non-invasive solver hybridization**.

The intended advantage is not that handoff has zero cost. Relay scheduling, representation conversion, verification, precision escalation, and handoff-data construction may all introduce overhead. The intended advantage is that these costs occur at the interface without requiring the internal algorithms to be merged into a single homogenized method.

## Universal Implicit Solver Relay Interface (UISRI)

The **Universal Implicit Solver Relay Interface (UISRI)** is the implicit-solver specialization of USRI.

UISRI permits implicit solvers with unrelated or mutually incompatible internal formulations to cooperate through structured handoff data while preserving the mathematical contract of each solver core.

Typical UISRI handoffs include:

- a non-residual constructive solver producing a state, branch, direction, or candidate domain for a residual-driven corrector;
- a local corrector returning a converged state that becomes the next construction anchor;
- an interval or Krawczyk-style method validating or narrowing a candidate domain before handing it to a local solver;
- a pseudo-arclength or other continuation method carrying a branch across a fold and returning a new state and tangent;
- a Newton-Krylov or other large-scale corrector consuming a state, preconditioner, subspace, or local structure generated upstream;
- a higher-order local method such as Halley or Householder consuming a sufficiently accurate state and local high-order information for terminal convergence.

An alternating UISRI relay may therefore take forms such as

\[
\text{NRII} \rightarrow \text{Newton} \rightarrow \text{NRII} \rightarrow \text{Newton} \rightarrow \cdots,
\]

while a multi-solver relay may take forms such as

\[
\text{NRII} \rightarrow \text{Krawczyk} \rightarrow \text{Newton} \leftrightarrow \text{Pseudo-arclength}.
\]

These expressions describe information handoff, not a redefinition of the internal algorithms. In particular, using Newton, Krawczyk, continuation, or another consumer does not modify the NRII coefficient-construction core unless a separate experiment explicitly defines a different method.

## Relationship between NRII, UISRI, and USRI

NRII is not defined as the universal interface itself. It is one possible producer of solver-handoff data and, in some experiments, one possible solver core inside UISRI.

The relationship is therefore

\[
\text{USRI} \supset \text{UISRI},
\]

with NRII participating in UISRI when its constructed coefficient, trajectory, branch, state, domain, or related information is handed to another implicit solver.

The broader USRI abstraction is intentionally independent of NRII. For example, a factorization-based solver may hand a low-precision solution and reusable factorization to a high-precision residual module, which then returns a more accurate residual for correction by the original factorization. This remains a USRI relay even though NRII is not involved.

## Scope of the word "Universal"

"Universal" refers to the **scope of the interface architecture**, not to a theorem that every pair of numerical solvers is mathematically compatible.

A solver can participate when a meaningful producer-consumer relation can be defined between its outputs and another solver's inputs. Some solver pairs may have no useful handoff relation, and some problems may terminate a relay because of singularity, discontinuity, branch termination, nonexistence of a valid domain, loss of numerical meaning, or other mathematical obstructions.

USRI and UISRI therefore make no claim that relay composition guarantees convergence, uniqueness, correctness, or improved performance for arbitrary problems.

## Experimental evidence in this repository

The current experimental branches explore different parts of this interface concept rather than constituting a universal proof.

- `experiment/nrii-newton-avalanche-hybrid` studies NRII as a branch-aware seed producer for Newton in a PN drift-diffusion-Poisson impact-ionization model.
- `experiment/nrii-newton-alternating-continuation` studies reciprocal and alternating NRII/Newton handoff in the same nonlinear impact-ionization model.
- `experiment/nrii-solver-relay-transcendental-fold` studies state, domain, branch, direction, Krawczyk, Newton, and pseudo-arclength handoff across a transcendental fold.
- `experiment/solver-relay-kahan-ill-conditioned` studies a non-NRII relay in which reusable factorization and higher-precision residual information are transferred between numerical modules on a severely ill-conditioned Kahan system.
- `experiment/nrii-semismooth-published-impacting-bar-relay` studies smooth NRII bulk evolution and semismooth unilateral-contact resolution on a published-parameter elastic-bar impact adaptation.
- `experiment/nrii-zero-hysteresis-buck-sliding-relay` studies smooth fixed-switch NRII branches and projected semismooth resolution of a zero-hysteresis electrical sliding singularity.
- `experiment/nrii-constrained-dynamics-singular-kkt-bare-c` studies a bare-C17 dense constrained multibody adaptation in which a full-rank Newton/KKT continuation contract loses numerical rank at a three-link manipulator singular configuration while a Matrix-Free NRII implicit-midpoint path continues the same independently verified smooth discrete state; a pseudoinverse specialist is retained as a singularity-aware control.
- `experiment/nrii-oscillating-field-moving-area-flux-bare-c` studies a bare-C17 combined electromagnetic/moving-surface stress benchmark in which a time-periodic magnetic-field vector, a rotating surface normal, and a periodically changing area generate analytically predictable carrier and sideband flux components; fixed-order Matrix-Free NRII reproduces the implicit-midpoint trajectory and Faraday-voltage spectrum without a global matrix solve.
- `experiment/nrii-induction-multiphysics-solver-relay-bare-c` studies a bare-C17 reduced thermo-electromagneto-hydrodynamic induction-furnace benchmark in which NRII, GMRES, BiCGSTAB, Picard, Newton-Krylov, ordinary Newton, semismooth active-set Newton, fixed-point coupling, Aitken relaxation, and Anderson acceleration retain separate solver contracts; electromagnetic contrast, strong convection, phase-boundary nondifferentiability, and global-coupling stress trigger different explicit solver handoffs while all five cases complete.
- `experiment/newton-hybrid-dynamics-tournament-bare-c` studies a bare-C17 128-DOF nonlinear structural-dynamics implicit-step family in which Newton is compared under the same manufactured branch oracle with line-search, Aitken, Anderson, Broyden-to-Newton, pseudo-transient-to-Newton, and audit-gated NRII-to-Newton hybrids; the branch records convergence speed and contiguous operational convergence envelopes rather than assuming that any hybrid must dominate.

These experiments support the architectural motivation for structured solver handoff and solver relay graphs, but they should not be read as evidence that USRI or UISRI dominates every standalone solver or every established continuation, interval, mixed-precision, projection, null-space, pseudoinverse, regularized, electromagnetic, structural-dynamics, multiphysics, or high-precision method.

## Experimental status

Experiment branches are exploratory and should not be read as universal performance or convergence claims. Each branch records its own model, numerical contract, validation conditions, comparison conditions, and known limitations.
