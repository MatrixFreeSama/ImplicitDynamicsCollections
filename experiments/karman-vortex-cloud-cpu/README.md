# Kármán bridge-wake cloud CPU experiment

This directory contains only **cloud-generated CPU reference material**. It intentionally does not contain the author's local CUDA executable or local GPU measurements.

## Files

- `nrii_cpu_cloud_bench.cpp`  
  Twelve-step native CPU reference comparing adaptive-order NRII + Chebyshev/Padé continuation against inexact Newton + analytic matrix-free Jv + restarted GMRES + Jacobi + Armijo.

- `nrii_karman_equal_time.cpp`  
  One-step residual-versus-wall-time driver. Each NRII order is an independent declared run; residual is measured after the run rather than used to increase the order inside the same attempt.

- `median_results.csv`  
  Median timing and final residual data from the 12-step experiment.

- `equal_time_results.csv`  
  Per-order NRII and per-Newton-update traces for 256² and 512².

## Cloud build

The recorded cloud configuration used:

```bash
g++ -O3 -march=native -fopenmp -std=c++17 nrii_cpu_cloud_bench.cpp -o nrii_cpu_cloud_bench
g++ -O3 -march=native -fopenmp -std=c++17 nrii_karman_equal_time.cpp -o nrii_karman_equal_time
```

The run used five OpenMP threads on an AMD EPYC 9V74 cloud CPU allocation.

## Reproducibility warning

Timing is hardware- and runtime-dependent. These CSV files are the measurements from the recorded cloud run, not canonical performance constants. A rerun on a different machine should be expected to produce different wall times while ideally preserving the qualitative residual trajectories if the numerical implementation is unchanged.
