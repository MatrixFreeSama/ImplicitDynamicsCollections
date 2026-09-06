# NRII Golden-Region Identification and Measured Fast-Regime Characteristics

**Benchmark snapshot:** 2026-09-06  
**Scope:** smooth, non-destructive implicit systems  
**Arithmetic in the reported C microbenchmarks:** FP64  

## Abstract

This note defines a practical procedure for identifying the operating region in which NRII is both numerically suitable and computationally advantageous. The region is not defined by a single time-step threshold and is not equivalent to the condition \(q<1\). It is the intersection of four properties: rapid coefficient contraction, bounded transient amplification, locality of the underlying operator, and a sufficiently expensive competing global correction solve.

The evidence combines the existing 572-point operating-envelope suite with a new 158-point three-dimensional localization campaign and the previously measured localized-hotspot study. The new campaign includes reaction-diffusion, a local fourth-order phase-field proxy, smooth vector coupling, two-field local multiphysics, nonlinear heat conduction with temperature-dependent conductivity, Allen-Cahn dynamics with a Newton-PCG baseline, and a deliberately nonlocal global-mean negative control.

The measured pattern is consistent across these studies. NRII is fastest when the original physics is local, the coefficient order remains small, and the conventional implicit method converts that local physics into a costly globally coupled correction solve. Spatially localized activity can increase the advantage further because the set of nodes that require high-order coefficient propagation may shrink with coefficient order. Conversely, low \(q\) alone is not sufficient: a cheap direct solve or an unusually favorable PCG system can erase the runtime advantage, and true physical nonlocality cannot be removed by coefficient propagation alone.

## 1. Definition of the golden region

For the implicit relation

\[
Y=X+hG(Y),
\]

NRII constructs

\[
Y(h)=\sum_{k=0}^{p}C_k h^k,
\qquad
C_{k+1}=[h^k]G\left(\sum_j C_jh^j\right).
\]

Define the coefficient term

\[
T_k=C_kh^k.
\]

The measured tail ratio is

\[
q_{\mathrm{eff}}
\approx
\operatorname{median}
\left(
\frac{\|T_{k+1}\|}{\|T_k\|}
\right)
\]

over the last several usable coefficient ratios. The transient-amplification diagnostic is

\[
\Gamma_{\mathrm{transient}}
=
\frac{\max_k\|T_k\|}{\|T_0\|+\varepsilon}.
\]

The practical NRII golden region is defined here as the subset of the problem-state space for which

\[
q_{\mathrm{eff}}<q_*,
\qquad
\Gamma_{\mathrm{transient}}\epsilon_{\mathrm{machine}}\ll\epsilon_{\mathrm{target}},
\qquad
T_{\mathrm{NRII}}<T_{\mathrm{reference}},
\]

and the coefficient recurrence does not introduce a global operation equivalent in cost to the global solve that it is intended to replace.

This definition deliberately separates three questions:

1. **Convergence suitability:** does the coefficient tail contract rapidly enough?
2. **Numerical safety:** does transient coefficient growth remain compatible with the requested precision?
3. **Computational economy:** is coefficient propagation cheaper than the best appropriate competing solve for the same problem?

A low tail ratio answers only the first question.

## 2. Preflight estimation before a full NRII solve

A complete solve is not required to estimate whether the current state lies near the favorable region.

### 2.1 Local spectral estimate

With

\[
J_0=J_G(X),
\]

a first preflight estimate is

\[
q_{\mathrm{pred}}^{(1)}\approx |h|\rho(J_0).
\]

For a linear system this reduces to the usual Neumann-series condition. In a nonlinear problem, a more conservative estimate uses a local Jacobian norm or Lipschitz bound over the predicted state neighborhood,

\[
q_{\mathrm{safe}}=|h|L,
\qquad
L\ge \sup_{Y\in\mathcal B}\|J_G(Y)\|.
\]

A practical neighborhood radius can be initialized from

\[
r_0=\eta\|hG(X)\|,
\qquad \eta>1,
\]

with the Jacobian sampled at the initial state and one or more predicted states inside this neighborhood.

For matrix-free implementations, neither estimate requires assembling a Jacobian. Short power or Arnoldi probes using only \(v\mapsto Jv\) are sufficient.

### 2.2 Transient-amplification probe

The spectral radius does not detect nonnormal transient growth. For a linearized operator, propagate several probe vectors,

\[
v_j^{(k+1)}=hJ_0v_j^{(k)},
\]

and record

\[
\widehat\Gamma
=
\max_{j,k}
\frac{\|v_j^{(k)}\|}{\|v_j^{(0)}\|}.
\]

The probe set should include the current-state direction when meaningful, the current residual or defect direction when available, and several independent generic directions. A small number of vectors and a depth of roughly six to twelve operator applications is sufficient for a low-cost warning test; it is not a rigorous worst-case bound.

For nonlinear problems, the same probe can be repeated at the initial state, a half-step predictor, and a full explicit predictor. A useful companion diagnostic is the Jacobian drift

\[
\eta_J=
\frac{\|J_G(X+\Delta X)-J_G(X)\|}
{\|J_G(X)\|+\varepsilon},
\qquad
\Delta X=hG(X).
\]

Large \(\eta_J\) indicates that a single-state spectral estimate is unlikely to predict the complete coefficient sequence reliably.

### 2.3 Coefficient pilot

The preferred runtime gate is a short coefficient pilot. Generate four to eight NRII orders without committing to a full solve, then measure

\[
q_{\mathrm{pilot}}
=\operatorname{median}
\left(
\frac{\|T_{k+1}\|}{\|T_k\|}
\right),
\qquad
\Gamma_{\mathrm{pilot}}
=\frac{\max_k\|T_k\|}{\|T_0\|+\varepsilon}.
\]

If NRII is retained, these coefficients are reused by the full solve, so the pilot is not discarded work. The pilot also provides an order estimate. For an approximately geometric tail,

\[
p_{\mathrm{req}}
\approx
\frac{\ln(\epsilon_{\mathrm{target}}/A)}{\ln q_{\mathrm{pilot}}}.
\]

The pilot therefore provides a direct estimate of both convergence depth and cost.

## 3. Cost-side test

NRII is not selected from \(q\) and \(\Gamma\) alone. The competing solver must also be modeled.

A simplified Newton-Krylov cost model is

\[
T_N
\approx
n_N
\left[
 n_K(C_{Jv}+C_{\mathrm{red}}+C_{\mathrm{prec}})
 +C_R
\right],
\]

where \(n_N\) is the number of nonlinear corrections, \(n_K\) is the Krylov depth, \(C_{\mathrm{red}}\) accounts for global reductions and orthogonalization, and \(C_{\mathrm{prec}}\) accounts for preconditioning.

A corresponding NRII model is

\[
T_R\approx p_{\mathrm{req}}C_{\mathrm{coeff}}.
\]

For spatially adaptive coefficient propagation,

\[
T_R^{\mathrm{active}}
\propto
N\sum_k f_k,
\]

where \(f_k\) is the fraction of nodes that remain active at coefficient order \(k\).

The principal performance coordinate is therefore not \(q\) alone but the ratio

\[
\mathcal A
\sim
\frac{n_N n_K C_{\mathrm{global}}}
{p_{\mathrm{req}}C_{\mathrm{local}}},
\]

with a further gain possible when \(f_k\) decreases strongly with \(k\).

## 4. Structural criterion: solver-induced globality

The three-dimensional localization study distinguishes solver-induced globality from physical nonlocality.

A local operator has the form

\[
G_i(U)=G_i(U_{\mathcal N(i)}),
\]

where \(\mathcal N(i)\) is a finite local neighborhood. A Newton or Newton-Krylov discretization can remain matrix-free and still require a globally coupled correction solve,

\[
J(U_k)\Delta U=-F(U_k).
\]

When the NRII coefficient recurrence preserves the original finite-neighborhood structure,

\[
C_{k+1,i}=\mathcal G_i(C_{0:k,\mathcal N(i)}),
\]

the global correction solve is replaced by local coefficient propagation. This is stronger than avoiding matrix assembly; it removes solver-induced globality from the coefficient stage.

The distinction has a hard boundary. If the physical operator itself contains a true nonlocal term, such as

\[
\bar u=\frac1N\sum_i u_i,
\]

then each NRII order must also evaluate that global quantity. NRII cannot create locality that the original operator does not possess.

## 5. New 158-point 3D localization campaign

The new campaign contains four groups of measurements: a size probe, a continuous implicit-step sweep, a symmetric-problem PCG control, and two concrete smooth 3D models. A consolidated evidence table and scenario-level localization summary are stored under `benchmarks/nrii_golden_region/`.

### 5.1 Summary by scenario

| Scenario | Points | Min speed ratio | Median | Max | \(q_{\mathrm{eff}}\) range | Main observation |
|---|---:|---:|---:|---:|---:|---|
| Fourth-order phase-field proxy | 14 | 1.49 | 1.90 | 2.20 | 0.00465-0.01994 | Local recurrence; no NRII global reductions |
| Reaction-diffusion, 7-point | 14 | 0.88 | 2.06 | 2.77 | 0.00569-0.02468 | Local recurrence; speed advantage is not uniform at every timing point |
| Two-field local multiphysics | 14 | 1.98 | 2.60 | 3.24 | 0.00403-0.01702 | Local block coupling remains local under coefficient propagation |
| Smooth vector solid proxy | 14 | 1.82 | 2.57 | 3.25 | 0.00529-0.02277 | Three-field local coupling remains favorable in the tested region |
| Nonlinear heat conduction \(k(T)\) | 12 | 1.64 | 2.34 | 3.09 | 0.000120-0.001404 | Concrete local face-flux model; Newton-GMRES baseline |
| Allen-Cahn | 12 | 0.38 | 1.15 | 1.60 | 0.00730-0.02246 | Concrete counterexample: locality is preserved but Newton-PCG can be cheaper |
| Global-mean negative control | 14 | 1.24 | 1.43 | 2.13 | 0.01881-0.08077 | True nonlocal term forces a global reduction at every NRII coefficient order |

The reported speed ratio is

\[
S_T=T_{\mathrm{Newton-family}}/T_{\mathrm{NRII}}.
\]

Values greater than one favor NRII. The measurements are controlled single-node C microbenchmarks and should not be interpreted as application-wide production solver rankings.

### 5.2 PCG fairness control

Reaction-diffusion, the fourth-order local operator, and the smooth vector-coupling proxy were repeated with Newton-PCG where the linearized structure permits it. This removes the artificial advantage that could arise from forcing a general GMRES solver onto a symmetric problem.

Measured ranges in the PCG control include approximately:

- reaction-diffusion: 1.30-2.50x;
- fourth-order local operator: 1.01-1.77x;
- smooth vector coupling: 1.25-2.58x.

The comparison shows that NRII localization can survive a more favorable Krylov baseline, but the margin contracts when the conventional linear system is inexpensive.

### 5.3 Nonlinear heat conduction

The concrete heat model uses a temperature-dependent conductivity and local face fluxes. The measured Newton-GMRES/NRII runtime ratio spans 1.64-3.09x over 24^3 through 48^3 grids and the tested implicit-step multipliers. NRII uses only three or four coefficient orders in these runs, with \(q_{\mathrm{eff}}\) between approximately \(1.2\times10^{-4}\) and \(1.4\times10^{-3}\).

This is a representative favorable structure: the original operator is local, the coefficient depth is very small, and the Newton formulation still requires a globally coupled correction.

### 5.4 Allen-Cahn counterexample

The Allen-Cahn runs use Newton-PCG. NRII preserves strict local coefficient propagation and remains accurate, but the measured speed ratio spans 0.38-1.60x. The conventional linearized system is sufficiently favorable to PCG that eliminating the global correction structure does not guarantee lower wall time.

This is the clearest measured counterexample to the rule "localizable implies faster."

### 5.5 Global-mean negative control

The global-mean model adds a true nonlocal coupling to an otherwise local operator. In the step sweep, Newton requires approximately ten global reductions in the tested configuration, while NRII requires one global mean evaluation per retained coefficient order, increasing from seven to eleven reductions as the required order rises.

This negative control establishes the structural limit:

\[
\boxed{\text{NRII can remove solver-induced globality; it does not remove intrinsic physical nonlocality.}}
\]

## 6. Previously measured high-acceleration cases

The existing operating-envelope suite and the localized-hotspot study identify two regimes in which the runtime separation becomes much larger than in ordinary smooth local PDE tests.

### 6.1 Localized 3D hotspot with active-sparse coefficient propagation

For the smooth localized-hotspot benchmark, the active fraction at coefficient order \(k\) is well approximated by

\[
f_k\approx f_0r^k,
\qquad r\approx0.79
\]

for the larger tested grids, with high coefficient of determination in the measured fits. An online active-sparse implementation propagates only the current active nodes and a local safety halo.

Against Newton-GMRES on the same controlled problem, the measured active-sparse NRII speed ratio spans approximately

\[
13.7\times\text{ to }22.1\times.
\]

At \(64^3\), the measured full-voxel-evaluation to active-evaluation ratio is approximately 11.55. Full NRII on the same hotspot is only about 1.2-2.4x faster than Newton-GMRES, showing that the large gain comes from the combination of low coefficient depth and shrinking spatial support rather than from coefficient propagation alone.

### 6.2 Controlled nonnormal Krylov-stagnation case

The existing nonnormal benchmark reaches a maximum measured runtime ratio of approximately 114x in a deliberately pathological regime where restarted GMRES stagnates while the NRII coefficient sequence remains numerically usable. This result demonstrates separation between Krylov difficulty and coefficient-propagation difficulty. It is not a representative engineering speedup and should not be quoted without the nonnormal-pathology qualification.

## 7. Empirical characteristics of the fastest NRII regimes

The fastest measured cases share the following properties.

### 7.1 Low coefficient depth

The effective tail ratio is small enough that the requested tolerance is reached in a small number of coefficient orders. This keeps

\[
p_{\mathrm{NRII}}\approx O(1)
\]

over the tested size range.

### 7.2 Local original physics

The governing residual or update at a node depends on a finite spatial neighborhood. NRII therefore retains local operator semantics instead of replacing the local physics with a global correction equation.

### 7.3 Expensive conventional global correction

The reference implicit method requires repeated global Krylov work, preconditioning, orthogonalization, reductions, or a difficult block solve. The strongest relative speedups occur when this cost is large compared with a coefficient sweep.

### 7.4 Controlled transient amplification

Small \(q_{\mathrm{eff}}\) is useful only while \(\Gamma_{\mathrm{transient}}\) remains compatible with the arithmetic precision. Nonnormality can create a low final tail ratio and still consume the floating-point error budget through intermediate growth.

### 7.5 Spatial activity contraction

When the set of nodes that require high-order terms decreases with coefficient order, active-sparse NRII can avoid most of the full-grid coefficient work. This mechanism produced the largest non-pathological speed ratios measured so far.

### 7.6 Absence of an unusually cheap specialized solve

A low \(q\) does not imply a performance win when the reference problem has a very cheap direct solve, a highly favorable PCG system, or another specialized global solver with a small constant factor. The scalar, tridiagonal, analytic-family, and Allen-Cahn controls all demonstrate this boundary.

These observations can be summarized as

\[
\boxed{\text{local operator} + \text{low }q + \text{controlled }\Gamma + \text{expensive global correction}}
\]

with a further acceleration term when the active spatial fraction decreases with coefficient order.

## 8. Runtime selection procedure

A practical solver selector can use the following sequence.

1. Estimate \(|h|\rho(J)\) or a conservative local Jacobian norm using matrix-free probes.
2. Probe short-term amplification in several directions to detect nonnormal transient growth.
3. Estimate Jacobian drift over a predicted state displacement.
4. Generate four to eight NRII coefficient orders.
5. Measure \(q_{\mathrm{pilot}}\), \(\Gamma_{\mathrm{pilot}}\), and the predicted remaining order.
6. Determine whether the recurrence is local or contains unavoidable global reductions or nonlocal inverses.
7. Estimate the competing direct/Krylov/preconditioned solve cost using the current problem structure.
8. Continue NRII only when the coefficient representation is safe and its predicted cost is lower; otherwise select Newton, another specialist, or a relay.

The selector should use continuous cost estimates rather than permanent fixed thresholds. Values such as \(q<0.05\) are useful empirical markers in the present experiments, but they are not universal constants.

## 9. Current boundary of the claim

The current evidence supports a structural statement rather than a universal complexity claim:

\[
\boxed{\text{NRII can convert some globally solved implicit workflows into bounded-order local propagation.}}
\]

This does not imply that all implicit problems become local, that all Newton methods require an assembled global matrix, or that all localizable problems run faster under NRII. Production AMG/ILU, geometric multigrid, unstructured finite elements, production multiphysics block preconditioners, distributed-memory scaling, and GPU implementations remain necessary comparison targets.

The current golden-region hypothesis is therefore falsifiable: performance should deteriorate when coefficient depth rises, transient amplification exceeds the precision budget, the original operator becomes intrinsically nonlocal, or the competing global solver becomes sufficiently cheap.
