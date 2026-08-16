# ImplicitDynamicsCollections

## NRII historical origin

NRII (Non-Residual-Iterative Implicit) was originally conceived as a direct, general-purpose constructive route for transcendental equations, not as an implicit-dynamics time integrator.

In that original design, the primary objective was to construct useful equation information directly. Convergence of the immediately emitted expression as a standalone final solution was not a required design target. The construction could instead be consumed by a downstream solver.

The implicit-dynamics work in this repository is a later application of the same coefficient-closure idea. Accordingly, experiments may study NRII in more than one role: as a standalone implicit construction, as a branch/trajectory information generator, or as a producer of structured input for another solver such as Newton.

## Universal Solver Relay Interface (USRI)

The **Universal Solver Relay Interface (USRI)** is an interface-level architecture for composing heterogeneous numerical solvers without requiring their internal update laws, convergence philosophies, numerical representations, precision models, or implementation strategies to be unified.

A solver participating in USRI does not need to produce the final solution. It may instead produce structured information that another solver can consume more effectively. Such transferable information is called a **Solver Baton**.

A relay is written abstractly as

\[
S_i \xrightarrow{\mathcal{B}_{i\rightarrow j}} S_j,
\]

where \(S_i\) and \(S_j\) are independently defined solver cores and \(\mathcal{B}_{i\rightarrow j}\) is the Solver Baton produced by \(S_i\) and consumed by \(S_j\).

The defining principle is:

> **Solver compatibility is determined by transferable information at solver boundaries, not by similarity of internal update laws.**

USRI is therefore intended to support heterogeneous solver composition without solver homogenization. A residual-driven solver, a non-residual constructive solver, an interval method, a continuation method, a factorization-based method, a matrix-free method, or a mixed-precision module may remain internally unchanged while participating through a common relay contract.

### Solver Baton

A **Solver Baton** is a structured information package passed between solver cores. It is not restricted to a candidate solution vector. Depending on the producer and consumer, a Baton may contain one or more of the following:

- **State Baton**: a candidate state, corrected state, anchor, or initial point.
- **Domain Baton**: an interval, trust region, candidate box, admissible neighborhood, or convergence region.
- **Branch Baton**: branch identity, ancestry, mode identity, or other information needed to preserve solution selection.
- **Direction Baton**: tangent, continuation direction, descent direction, local trajectory direction, or other directional information.
- **Geometry Baton**: Jacobian information, approximate inverse, local metric, curvature information, or other local geometric structure.
- **Subspace Baton**: Krylov space, coarse space, low-rank space, null-space estimate, or other reusable search subspace.
- **Factorization Baton**: an LU, QR, RRQR, preconditioner, or other reusable decomposition or operator structure.
- **Residual Baton**: a residual or defect evaluated by a solver or precision regime better suited to producing it.
- **Precision Baton**: precision requirements, error budgets, conditioning estimates, or information about which numerical components require increased precision.
- **Certificate Baton**: interval inclusion information, uniqueness information, admissibility gates, or other numerical or rigorous certificates.

A Baton may be extended with additional fields when a solver pair requires information not covered by the categories above. USRI defines the handoff concept rather than a closed list of solver-specific payloads.

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

The solver pair repeatedly hands information back and forth while each solver retains its native update rule.

**Multi-solver relay graph**

\[
S_1 \rightarrow S_2 \rightarrow \cdots \rightarrow S_n,
\]

or, more generally, a directed graph in which a scheduler selects the next solver according to the currently available Baton and the information still missing from the problem state.

### Interface-level hybridization

USRI treats hybridization as an interface problem rather than a requirement to derive a single blended update formula.

Two solvers may therefore be mathematically unrelated, internally incompatible, or based on opposing numerical philosophies and still cooperate if one can produce information the other can consume. Their native solver cores remain independently optimizable and independently testable.

This is referred to in this repository as **interface-level solver hybridization** or **non-invasive solver hybridization**.

The intended advantage is not that handoff has zero cost. Relay scheduling, representation conversion, certification, precision escalation, and Baton construction may all introduce overhead. The intended advantage is that these costs occur at the interface without requiring the internal strengths of participating solvers to be diluted into a single homogenized method.

## Universal Implicit Solver Relay Interface (UISRI)

The **Universal Implicit Solver Relay Interface (UISRI)** is the implicit-solver specialization of USRI.

UISRI permits implicit solvers with unrelated or even mutually incompatible internal formulations to cooperate through Solver Batons while preserving the native mathematical authority of each solver core.

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

NRII is not defined as the universal interface itself. It is one possible **Baton Producer** and, in some experiments, one possible solver core inside UISRI.

The relationship is therefore

\[
\text{USRI} \supset \text{UISRI},
\]

with NRII participating in UISRI when its constructed coefficient, trajectory, branch, state, domain, or related information is handed to another implicit solver.

The broader USRI abstraction is intentionally independent of NRII. For example, a factorization-based solver may hand a low-precision solution and reusable factorization to a high-precision residual module, which then returns a Residual Baton for correction by the original factorization. This remains a USRI relay even though NRII is not involved.

## Scope of the word "Universal"

"Universal" refers to the **scope of the interface architecture**, not to a theorem that every pair of numerical solvers is mathematically compatible.

A solver can participate when a meaningful producer-consumer relation can be defined between its outputs and another solver's inputs. Some solver pairs may have no useful Baton, and some problems may terminate a relay because of singularity, discontinuity, branch termination, nonexistence of a valid domain, loss of numerical meaning, or other mathematical obstructions.

USRI and UISRI therefore make no claim that relay composition guarantees convergence, uniqueness, correctness, or improved performance for arbitrary problems.

## Experimental evidence in this repository

The current experimental branches explore different parts of this interface concept rather than constituting a universal proof.

- `experiment/nrii-newton-avalanche-hybrid` studies NRII as a branch-aware seed producer for Newton.
- `experiment/nrii-newton-alternating-continuation` studies reciprocal and alternating NRII/Newton handoff in a nonlinear avalanche model.
- `experiment/nrii-solver-relay-transcendental-fold` studies state, domain, branch, direction, Krawczyk, Newton, and pseudo-arclength handoff across a transcendental fold.
- `experiment/solver-relay-kahan-black` studies a non-NRII relay in which reusable factorization and higher-precision residual information are handed between numerical modules on a severely ill-conditioned Kahan system.

These experiments support the architectural motivation for Solver Batons and solver relay graphs, but they should not be read as evidence that USRI or UISRI dominates every standalone solver or every established continuation, interval, mixed-precision, or high-precision method.

## Experimental status

Experiment branches are exploratory and should not be read as universal performance or convergence claims. Each branch records its own model, numerical contract, validation gates, comparison conditions, and known limitations.
