# ImplicitDynamicsCollections

## Experimental branch: NRII vs Newton–Krylov on a Kármán bridge-wake benchmark

> **Branch status:** experimental only. This branch is intentionally not merged into `main`.
>
> **Provenance:** the benchmark data and C++ sources in this branch were generated and executed in the assistant's cloud CPU environment. They are **not** measurements from the author's local RTX system.

This branch records an early controlled comparison between a forward coefficient-construction implicit method, **NRII (Non-Residual-Iterative Implicit)**, and a conventional **inexact Newton + matrix-free analytic Jv + restarted GMRES + Jacobi preconditioner + Armijo line search** baseline.

The benchmark is a two-dimensional bridge-deck wake problem designed to produce a Kármán-style alternating vortex street. The state is

\[
\mathbf U=(u,v,p),
\]

and both methods are evaluated against the same discrete implicit residual

\[
\mathbf R(\mathbf U_{n+1})
=
\mathbf U_{n+1}-\mathbf U_n-\Delta t\,\mathbf F(\mathbf U_{n+1}).
\]

The purpose of this branch is not to claim a production CFD solver or a definitive GPU result. It is a reproducible **cloud CPU reference experiment** for studying:

1. residual versus wall time,
2. time to the observed numerical residual floor,
3. the one-shot order/accuracy cost of NRII,
4. iteration growth in Newton–Krylov,
5. implementation overhead introduced by Chebyshev and rational continuation.

### Headline equal-time observations

For the one-step floor-seeking experiment:

| Grid | NRII first reaches the observed residual floor | Newton first reaches the observed residual floor | Approx. NRII time advantage |
|---:|---:|---:|---:|
| 256² | 4.78 ms, \(p=5\) | 22.62 ms, 2 Newton / 3 GMRES | 4.73× |
| 512² | 37.92 ms, \(p=7\) | 148.05 ms, 3 Newton / 6 GMRES | 3.90× |

The observed floor is around \(4.4\times10^{-9}\) in this CPU reference. This is an empirical plateau of this implementation and discretization, **not** a universal machine-precision theorem.

The longer 12-step median benchmark shows a different picture: NRII wins at 256², while Newton is faster at 128², 512², and 1024² in the current CPU implementation. This exposes substantial overhead in the present general Chebyshev/rational machinery and is one reason the branch is explicitly experimental.

### Repository layout on this branch

```text
README.md
paper/
  NRII_Karman_Cloud_CPU_Experiment.md
experiments/
  karman-vortex-cloud-cpu/
    README.md
    equal_time_results.csv
    median_results.csv
    nrii_cpu_cloud_bench.cpp
    nrii_karman_equal_time.cpp
```

### Important scope boundary

NRII itself is treated here as the coefficient-generation core. Chebyshev and Padé-family rational continuation are auxiliary representation tools, not definitions of NRII. The Newton–Krylov implementation is a comparison baseline.

No local Windows/CUDA executable, local GPU measurement, or unpublished local-machine result is included in this branch.

See [`paper/NRII_Karman_Cloud_CPU_Experiment.md`](paper/NRII_Karman_Cloud_CPU_Experiment.md) for the experimental note and [`experiments/karman-vortex-cloud-cpu/`](experiments/karman-vortex-cloud-cpu/) for raw data and source.
