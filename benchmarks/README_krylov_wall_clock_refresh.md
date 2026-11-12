# Krylov wall-clock refresh comparison

The numerical benchmark `benchmark_krylov_quadrature_refresh.f90` uses
midpoint quadrature inside nonlinear residual and analytic Jacobian-action
callbacks. `wall_timed` is an additional opt-in policy; the other modes
are unchanged. The included CSV records 84 GNU Fortran 14.2 `-O2` runs:
42 using `OMP_NUM_THREADS=1` and 42 using `OMP_NUM_THREADS=4`.

To compile with the SciFort sources and their dependencies, use Fortran
2023-compatible GNU Fortran 14.2 or another compiler supporting OpenMP
and `system_clock` int64 counters. Compile the benchmark with `-O2
-fopenmp` and run it twice with `OMP_NUM_THREADS=1` and `4` (and
`OMP_DYNAMIC=FALSE`). The script requires the SciFort root modules
and their normal dependencies. Benchmark builds used a compile-only
MINPACK interface stub; the benchmark does not invoke MINPACK solvers.

For the 128-variable problems with 2,048 quadrature points per variable,
mean elapsed wall seconds for the two nonlinearity rates in this run:

| Policy | 1 OpenMP thread | 4 OpenMP threads |
| --- | ---: | ---: |
| Initial preconditioner only | 0.1295 | 0.0591 |
| CPU-timed adaptive | 0.1002 | 0.0439 |
| Wall-timed adaptive | 0.0655 | 0.0340 |
| Fixed every four iterations | 0.0686 | 0.0266 |
| Measured products (1x transpose) | 0.0666 | 0.0258 |

Every run converged to the known solution and passed the benchmark's
scaled residual and accuracy checks. These are single-run, machine-
dependent timings, not evidence of universal superiority. Process CPU
time sums work across threads, while elapsed time accounts for parallel
execution and scheduling delays. Test repeatability across operating
systems, OpenMP runtimes, and machine loads before relying on a
wall-timed policy for production performance.
