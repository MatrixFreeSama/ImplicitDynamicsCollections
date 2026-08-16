# NRII versus Newton–Krylov on a Kármán Bridge-Wake Problem: An Experimental Cloud-CPU Reference

## Abstract

This experimental note compares a forward coefficient-construction implicit method, NRII (Non-Residual-Iterative Implicit), with a conventional Newton–Krylov baseline on a two-dimensional bridge-deck wake problem. The comparison is deliberately restricted to a cloud CPU reference implementation. It does not report CUDA performance and does not use measurements from the author's local machine.

The benchmark uses a common discrete state \(\mathbf U=(u,v,p)\), a common matrix-free spatial operator, and a common implicit residual. NRII constructs a finite sequence of time-series coefficients and then applies Chebyshev/rational representation machinery to form a candidate state without a residual-driven correction loop. The baseline uses inexact Newton iterations, analytic matrix-free Jacobian-vector products, restarted GMRES, a Jacobi preconditioner, and Armijo line search.

Two experiment styles are recorded. A 12-step median benchmark measures end-to-end wall time under a fixed residual gate. A one-step floor-seeking benchmark instead records residual as a function of computational work and wall time. The latter avoids treating a single user-chosen tolerance as the sole definition of performance and makes visible a structural cost of one-shot coefficient methods: once the numerical floor has been reached, additional expansion order produces little or no useful accuracy.

On the cloud CPU reference, NRII reaches the observed \(\sim 4.4\times10^{-9}\) residual plateau substantially earlier than Newton–Krylov at 256² and 512² in the one-step experiment. In the 12-step benchmark, however, the current general Chebyshev/rational implementation is sufficiently expensive that Newton is faster at several grid sizes. These results are therefore evidence about the present reference implementation and experimental metric, not a claim of universal superiority.

## 1. Experimental status and provenance

This document belongs to an experimental Git branch and is not part of the repository's `main` branch.

All CSV data committed with this note were produced in the assistant's cloud environment. The corresponding native C++ source is included verbatim except that one absolute `/mnt/data/` include path was replaced by a relative include for repository portability.

The measured environment reported by the cloud run was:

- CPU: AMD EPYC 9V74
- available logical CPU allocation used by the experiment: 5 vCPU
- OpenMP threads: 5
- compiler mode: `g++ -O3 -march=native -fopenmp -std=c++17`
- GPU: none used for these measurements

Accordingly, no result in this branch should be presented as RTX, CUDA, Blackwell, or local-machine performance.

## 2. Benchmark problem

The state vector is

\[
\mathbf U=(u,v,p).
\]

The bridge-wake scene uses a fixed bridge-deck obstacle and a two-dimensional velocity-pressure field intended to exhibit alternating wake vortices. The common implicit step is tested by the discrete residual

\[
\mathbf R(\mathbf U_{n+1})
=
\mathbf U_{n+1}
-
\mathbf U_n
-
\Delta t\,\mathbf F(\mathbf U_{n+1}).
\]

The reference parameters used by the 12-step cloud benchmark are:

\[
\Delta t=5\times10^{-5},
\qquad
\nu=2.5\times10^{-4},
\qquad
\beta^2=1,
\qquad
\sigma_p=0.05.
\]

The fixed residual gate in that benchmark is

\[
\|\mathbf R\|_{\mathrm{RMS}}\le 5\times10^{-6}.
\]

The one-step floor-seeking experiment does not use this tolerance to adapt NRII order. Each NRII order is launched as an independently declared run; residual is measured afterward.

## 3. Methods

### 3.1 NRII

The NRII core is represented conceptually by

\[
\mathbf C_0=\mathbf U_n,
\qquad
\mathbf C_k=[h^{k-1}]\mathbf G(\mathbf U(h)).
\]

The defining property relevant to this experiment is the absence of a residual-driven Newton-style correction loop. For a declared order \(p\), the coefficient sequence is generated forward and used to construct a candidate state.

In the reference code, the declared order also controls the available Chebyshev and rational representation information. The rational denominator degree is component-adaptive within the information supplied by the declared order.

### 3.2 Newton–Krylov baseline

The baseline uses:

- inexact Newton iterations,
- analytic matrix-free Jacobian-vector products,
- restarted GMRES,
- Jacobi preconditioning,
- Armijo line search.

The baseline solves the same discrete implicit residual used to audit NRII candidates.

## 4. Metrics

Two metrics are emphasized.

### 4.1 Residual at equal wall time

Rather than defining success only by a single arbitrary tolerance, the experiment records

\[
R(t)=\|\mathbf R(\mathbf U(t))\|_{\mathrm{RMS}}.
\]

This makes excess accuracy visible as useful progress until a numerical plateau is reached.

### 4.2 Time to observed numerical floor

For each grid, the observed floor is taken from the residual plateau actually reached by the reference implementation. A method is considered to have reached the floor when further work produces negligible improvement on that plateau.

This is an empirical measurement. It is not identified with raw IEEE machine epsilon.

## 5. One-step floor-seeking results

### 5.1 Grid 256²

NRII samples:

| Order \(p\) | Wall time (ms) | RMS residual |
|---:|---:|---:|
| 2 | 1.7206 | \(1.3740\times10^{-6}\) |
| 3 | 4.0880 | \(2.1683\times10^{-8}\) |
| 5 | 4.7846 | \(4.4717\times10^{-9}\) |
| 7 | 6.2615 | \(4.4675\times10^{-9}\) |
| 12 | 58.3931 | \(4.4675\times10^{-9}\) |

Newton samples:

| Newton updates | Wall time (ms) | RMS residual |
|---:|---:|---:|
| 0 | 0 | \(1.1638\times10^{-4}\) |
| 1 | 8.7051 | \(1.2928\times10^{-6}\) |
| 2 | 22.6171 | \(4.4714\times10^{-9}\) |
| 3 | 33.6968 | \(4.4674\times10^{-9}\) |

The first NRII point essentially on the observed floor occurs at approximately 4.78 ms with \(p=5\). Newton reaches a comparable floor after approximately 22.62 ms. The ratio is approximately

\[
22.62/4.78\approx 4.73.
\]

### 5.2 Grid 512²

NRII samples:

| Order \(p\) | Wall time (ms) | RMS residual |
|---:|---:|---:|
| 2 | 6.8418 | \(4.3040\times10^{-6}\) |
| 3 | 13.2888 | \(1.4315\times10^{-7}\) |
| 5 | 19.4380 | \(5.7562\times10^{-9}\) |
| 7 | 37.9238 | \(4.3586\times10^{-9}\) |
| 10 | 97.8765 | \(4.3572\times10^{-9}\) |
| 14 | 202.4071 | \(4.3572\times10^{-9}\) |

Newton samples:

| Newton updates | Wall time (ms) | RMS residual |
|---:|---:|---:|
| 0 | 0 | \(1.7603\times10^{-4}\) |
| 1 | 49.8579 | \(4.1614\times10^{-6}\) |
| 2 | 90.8150 | \(5.6054\times10^{-9}\) |
| 3 | 148.0502 | \(4.3559\times10^{-9}\) |

NRII first reaches the observed floor neighborhood at approximately 37.92 ms with \(p=7\). Newton reaches a comparable plateau at approximately 148.05 ms. The ratio is approximately

\[
148.05/37.92\approx 3.90.
\]

At roughly 50 ms of wall time, the NRII reference is already near \(4.36\times10^{-9}\), while the first Newton update is at \(4.16\times10^{-6}\), a residual difference of roughly three orders of magnitude at that stage.

## 6. Twelve-step median benchmark

A separate experiment runs 12 implicit steps and compares median wall time under the common \(5\times10^{-6}\) residual gate.

| Grid | DOF | NRII order | Padé \(M_{u/v/p}\) | NRII median (ms) | Newton median (ms) | Newton / NRII | NRII residual | Newton residual |
|---:|---:|---:|:---:|---:|---:|---:|---:|---:|
| 128² | 49,152 | 4 | 2/2/2 | 13.76 | 8.57 | 0.623× | \(5.50\times10^{-9}\) | \(4.29\times10^{-7}\) |
| 256² | 196,608 | 5 | 2/2/2 | 57.55 | 91.27 | 1.586× | \(4.65\times10^{-9}\) | \(1.24\times10^{-6}\) |
| 512² | 786,432 | 6 | 3/3/3 | 489.36 | 335.61 | 0.686× | \(5.44\times10^{-9}\) | \(3.60\times10^{-6}\) |
| 1024² | 3,145,728 | 7 | 3/3/3 | 2611.41 | 1518.08 | 0.581× | \(5.42\times10^{-9}\) | \(3.84\times10^{-8}\) |

This benchmark does not show monotonic NRII superiority. In the present CPU implementation, Newton is faster at 128², 512², and 1024². NRII is faster at 256².

The reference code makes a likely implementation cost visible: the general Chebyshev transform scales approximately as \(O(p^2N)\), and the rational identification performs repeated whole-field reductions. These costs grow sharply when the representation order increases.

The attempted 2048² run did not complete within the 120 s cloud execution window and is not assigned a fabricated result.

## 7. Accuracy over-resolution as a structural cost

The one-step data expose a characteristic cost of one-shot coefficient construction.

At 512², for example,

\[
p=7
\]

already gives approximately

\[
4.36\times10^{-9},
\]

and increasing the order to \(p=14\) leaves the residual effectively unchanged while wall time grows substantially.

If NRII is required to remain a one-shot forward construction rather than a residual-driven outer iteration, the exact minimal sufficient order is generally not known by repeatedly solving and correcting within the same attempt. A conservative order can therefore produce accuracy beyond what the application needs.

This branch calls that phenomenon an **accuracy over-resolution cost** or **one-shot accuracy tax**. It is an efficiency limitation, not a CFL stability restriction.

## 8. Limitations

The following limitations are essential:

1. This is a CPU reference, not a GPU benchmark.
2. The bridge-wake model is a benchmark scene, not a validated full-scale bridge engineering CFD model.
3. The observed residual floor is implementation- and discretization-dependent.
4. The present Chebyshev/rational implementation contains substantial overhead and should not be treated as an optimized endpoint.
5. The data are too limited to establish asymptotic superiority of either solver family.
6. No claim is made that NRII dominates Newton–Krylov for arbitrary nonlinear implicit systems.
7. The cloud source is a research reference and should not be conflated with the separate local Windows/CUDA executable.

## 9. Reproducibility files

Raw experimental data:

- `experiments/karman-vortex-cloud-cpu/equal_time_results.csv`
- `experiments/karman-vortex-cloud-cpu/median_results.csv`

Source:

- `experiments/karman-vortex-cloud-cpu/nrii_cpu_cloud_bench.cpp`
- `experiments/karman-vortex-cloud-cpu/nrii_karman_equal_time.cpp`

The equal-time driver includes the benchmark source through a relative include so that the repository version is portable.

## 10. Current interpretation

The cloud experiment supports three narrow observations.

First, the CPU reference implementation can produce NRII candidates whose measured implicit residual reaches the same numerical plateau as Newton–Krylov on the tested 256² and 512² one-step cases.

Second, the residual-versus-time view is more informative for this one-shot method than a single fixed tolerance, because it distinguishes useful additional accuracy from work performed after the numerical floor has already been reached.

Third, the present generalized Chebyshev/rational implementation has enough overhead to reverse the performance ranking in the longer multi-step benchmark at several grid sizes. Removing that implementation overhead is therefore a separate engineering problem from evaluating the NRII coefficient construction itself.

No stronger conclusion is claimed in this experimental branch.
