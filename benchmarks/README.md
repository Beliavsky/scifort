## Cost-aware FFT convolution padding (2026-10-26)

`bench_signal_fft_smooth_padding.f90` measures four one-dimensional and four
two-dimensional convolution workloads. Compile the same benchmark against the
old and new `src/signal_fftconvolve.f90` implementations, both with `-O3`,
along with `kinds.f90`, `signal_dpss.f90`, `signal_windows_extra.f90`, and
`signal_spectral.f90`. All workloads produce deterministic checksums.

In eight interleaved local GNU Fortran 14.2 runs, the median times for the
two-dimensional case with inputs of shape `(101,251)` and `(29,59)` were
0.03359 s (power-of-two padding) and 0.02486 s (cost-aware padding), a 1.35x
speedup. The selected padded shape shrank from `(256,512)` to `(135,315)`.
Most other cases saw substantially smaller differences. The recursive
mixed-radix FFT has higher per-element cost than radix-2; using the smallest
possible smooth length everywhere is **not** a reliable speed optimization.
Timings depend on hardware, compiler, and problem shape and should not be
used as pass/fail thresholds.

## Analytic spline-gradient benchmark (2026-10-24)

`bench_ndimage_spline_gradients.f90` compares the cached spline plan's
analytic value-and-gradient evaluator with 2-sided coordinate differences,
using a 64x60 real image, 900 points, fifth-order splines and four repeats.
A local GNU Fortran 14.2.0 `-O2` run measured 9.589 ms analytic versus
20.844 ms for two-dimensional finite differences (2.17x speedup); maximum
absolute derivative difference was 4.33e-10. Results are workload-specific.

# Elementwise typed-argument performance benchmark

`benchmark_elementwise_typed_args.f90` compares all four elementwise packed
solvers using the existing copied typed callback and the indexed, zero-copy
typed callback. Each objective reads a scalar field from a derived-type value
containing a large allocatable array. The unused array is intentional: it
measures the cost of copying the whole payload for a small objective.

Compile from the project root without needing any external FPM dependencies:

```sh
mkdir -p _bench_build
cd _bench_build
gfortran -O2 -std=f2018 -ffree-line-length-none -fimplicit-none -Werror=implicit-interface \
    -c ../src/kinds.f90 ../src/optimize/elementwise_batch.f90 \
    ../src/differentiate.f90 ../src/optimize/root_scalar.f90 \
    ../src/optimize/minimize_scalar.f90
gfortran -O2 -std=f2018 -ffree-line-length-none -fimplicit-none \
    -Werror=implicit-interface ../benchmarks/benchmark_elementwise_typed_args.f90 \
    ./*.o -o benchmark_elementwise_typed_args
./benchmark_elementwise_typed_args 512 4096 10
```

The optional integer arguments are the number of independent problem lanes,
the number of `real(dp)` values allocated inside each context object, and the
number of timing repetitions. Defaults are `256 512 3`.

The program builds and boxes its argument list **before** timing, warms up each
method, alternates the order of copied and indexed measurements, and reports
mean process CPU seconds per solver invocation. It checks that every lane
succeeds and that both callback styles produce equal result checksums and
function-evaluation counts. The four methods are `find_root`, `find_minimum`,
`bracket_root`, and `bracket_minimum`.

For meaningful comparisons, use an optimized build and increase repetitions
when indexed times are under a millisecond. Results vary with memory bandwidth,
compiler, CPU frequency and payload size. Do not use speed ratios as CI
pass/fail thresholds. The main benefit is avoidance of repeated derived-type
payload copying; differences caused by the reusable integer-offset work buffer
are substantially smaller and may be within run-to-run timing noise.

## LGMRES pruning convergence/cost comparison

`test/test_optimize_root_krylov_pruning_comparison.f90` is a deterministic
regression and comparison driver for seven analytic augmentation policies: unpruned
LGMRES, `outer_prune_tol=1e-7` (solution directions),
`outer_image_prune_tol=1e-7` (current Jacobian images), and
`outer_image_prune_auto=.true.` (automatically selected image-space cutoff),
retained-image reuse, candidate-image reuse, and both reuse options together.

From the SciFort root, with FPM and the declared dependencies installed:

```sh
fpm test --target test_optimize_root_krylov_pruning_comparison
```

The program prints CSV-compatible records before its final `PASS` line. It
uses 16 variables and a dense matrix constructed as `Q1*D*Q2**T`, where
`Q1` and `Q2` are deterministic Householder reflections and `D` has
geometrically spaced singular values. The nominal *linear-factor* condition
numbers are `1e3`, `1e4`, and `1e7`. The nonlinear residual is

```
F(x) = A * [e + 0.3*e**2 + 0.03*q*(dot_product(q,e))**2], e = x - x_star.
```

Its Jacobian changes with `x`; the analytic `J*v` is passed as an explicit
callback. Each policy uses `inner_maxiter=12`, `outer_k=4`, selective Arnoldi
reorthogonalization, an Armijo line search, `fatol=1e-9`, and `maxiter=80`.
The driver also repeats five applicable policies for the `1e3` case using two-point numerical derivatives
with `fatol=1e-8`, checking that fixed and automatic image pruning's extra function
evaluations are included in `nfev`.

Columns are `condition`, `product`, `policy`, `success`, `nit`, `nfev`,
`njvp`, `max_residual`, and `max_x_error`. `njvp` counts calls to the
*analytic* Jacobian-vector callback; it is zero for finite-difference runs,
whose extra samples instead contribute to `nfev`.

`lgmres_pruning_gfortran14_o2.csv` captures the GNU Fortran 14.2 `-O2`
results from the checked focused build (with `-std=f2018`, `-fcheck=all`,
`-ffpe-trap=invalid,zero,overflow`). Results may differ with compiler,
optimization, and floating-point ordering; this file is reference output,
**not** a golden file for assertions of exact counts or rankings. This stress
case is not a statistical claim of superiority for any pruning policy.
On the GNU Fortran 14.2 `-O2` reference cases, automatic and fixed
image pruning happened to produce identical nonlinear iterations, final
residuals, and work counts; this is not a universal equivalence. The
separate automatic-pruning unit tests verify that the projection-dependent
cutoff differs from fixed thresholds on near-dependent image vectors.
The full dependency-backed FPM test suite was not run in the snapshot
that generated these results.

## LGMRES Jacobian-image reuse scaling benchmark

`benchmark_lgmres_image_cache.f90` compares four analytic-Jacobian-vector
configurations for **32, 64, and 128 variables**: automatic image pruning
without caching, retained-image caching (`outer_image_reuse`), candidate-step
caching (`outer_candidate_image_reuse`), and both. Each mode is repeated 30
times and reports mean `cpu_time`, iterations, residual evaluations,
Jacobian-vector callback calls, and achieved maximum residual. The benchmark
uses the same deterministic dense Householder/singular-value construction
as the 16-variable regression, nominal matrix condition number `1e3`,
`inner_maxiter=12`, `outer_k=4`, and Armijo line search. Do not use raw elapsed
times or performance rankings as regression pass/fail requirements.

With FPM and the declared external dependencies installed, one simple way to
compile and run this *optional* benchmark is to temporarily copy the program
into the project's automatically discovered executable directory:

```sh
cp benchmarks/benchmark_lgmres_image_cache.f90 app/
fpm run --target benchmark_lgmres_image_cache --flag '-O2'
rm app/benchmark_lgmres_image_cache.f90
```

The GNU Fortran 14.2 `-O2` **checked-build** reference results are in
`lgmres_image_cache_scaling_gfortran14_o2.csv`. They were measured using
`-fcheck=all`, `-ffpe-trap=invalid,zero,overflow` and repeated runs; they are
not optimized production benchmarks. Compiler, matrix multiplication library,
CPU and runtime checking flags can change timings and even nonlinear iteration
counts. The measured savings in callback counts on the reference build ranged
from 23% to 34% when both caches were enabled, but some cases also followed
different nonlinear paths. Only the separate focused regressions enforce
numerical correctness and callback accounting.


## LGMRES image-cache memory benchmark

`benchmark_lgmres_cache_memory.f90` solves a large, diagonally scaled
nonlinear system using an analytic, matrix-free Jacobian-vector callback.
It performs no dense matrix allocation, so the process's resident memory
reflects Krylov basis/augmentation workspaces rather than an `n*n` Jacobian.
The test enables automatic image pruning, retained-image reuse, and
candidate-image reuse simultaneously. An optional command-line argument
selects the number of variables (default `200000`). Both this snapshot and
the immediately preceding candidate-image-reuse snapshot have the same
public API and can compile the *same* benchmark source.

After compiling the benchmark against each version (for example by copying
it into `app/` temporarily and running `fpm build`), measure peak RSS on
Linux with:

```sh
/usr/bin/time -f 'peak_kb=%M elapsed_sec=%e' ./benchmark_lgmres_cache_memory 250000
```

The `%M` output is peak resident memory in KiB for the whole process, not an
allocation profiler. GNU Fortran 14.2, `-O2 -fcheck=all`, and floating-point
traps produced these single-run reference measurements:

| Variables | Previous peak RSS | Compact peak RSS | Reduction |
| ---: | ---: | ---: | ---: |
| 250,000 | 184,992 KiB | 145,812 KiB | 21.2% |
| 500,000 | 337,112 KiB | 294,192 KiB | 12.7% |

Both implementations converged in six nonlinear iterations with six residual
evaluations and 46 analytic J*v callbacks at each size, with the same reported
maximum residual. Timings and peak memory depend on compiler, runtime flags,
allocator, OS, and system load. These figures demonstrate reduced memory for
this benchmark, not a guaranteed improvement for every solver workload.
The raw results are in `lgmres_cache_workspace_gfortran14_o2.csv`. The
`test_optimize_root_krylov_cache_compaction.f90` regression checks numerical
and callback behavior, including preconditioning and one-column history.

## GMRES/LGMRES Arnoldi workspace comparison (2026-10-11)

The same `benchmark_lgmres_cache_memory.f90` benchmark can compare the
original compact-image-cache snapshot with the per-solve Arnoldi-workspace
version. Build both snapshots with identical compiler and flags, then use
`/usr/bin/time -f 'peak_kb=%M elapsed_sec=%e'` on the executable with
`250000` and `500000` as its command-line argument. The benchmark needs
no dense Jacobian, and checks the nonlinear solution and callback counts.
The source uses analytic J*v and both image-cache options. This is a
whole-process peak-RSS measurement, not a byte-accurate allocation profile.

`lgmres_arnoldi_workspace_gfortran14_o2.csv` records matched single-run
GNU Fortran 14.2, checked `-O2` results. Peak RSS changed from 145,820 to
141,876 KiB (250k) and from 294,168 to 282,528 KiB (500k), while the
iterations and Jacobian/residual callback counts remained unchanged.
Measurements from different machines, compilers, allocator states, or run
orders may disagree; do not assert exact performance measurements in tests.

## Lazy Arnoldi and Jacobian-image workspaces (2026-10-11)

`benchmark_lgmres_cache_memory.f90` now optionally accepts an `outer_k`
capacity as its second argument (default 8). A third argument, `precond`,
enables a deterministic diagonal approximate-inverse preconditioner. For
example, after building the benchmark with the same compiler settings against
both the previous `lgmres-arnoldi-workspace-reuse` and the updated source:

```sh
/usr/bin/time -f 'peak_rss_kib=%M elapsed_sec=%e' ./benchmark_lgmres_cache_memory 250000 64
/usr/bin/time -f 'peak_rss_kib=%M elapsed_sec=%e' ./benchmark_lgmres_cache_memory 150000 64 precond
```

With an oversized `outer_k`, the previous version reserved all possible
Arnoldi basis columns and zero-filled `n*outer_k` Jacobian image elements.
The updated version allocates an initial two extra Arnoldi columns, grows
the basis in powers of two as history expands, and allocates image space for
only the current retained directions plus a candidate. It also lets a
preconditioner write directly into the persistent Krylov direction. When the
retention history is full, both versions can require the same memory.

`lgmres_lazy_workspace_gfortran14_o2.csv` contains matched single-run
measurements from GNU Fortran 14.2 with `-O2 -fcheck=all` and floating-point
traps. For 250,000 variables and `outer_k=64`, peak RSS changed from 366,532
to 251,408 KiB (31.4% less). At 500,000 variables and `outer_k=64`,
the corresponding values were 731,840 and 512,956 KiB (29.9% less).
With `outer_k=8`, the change was within normal
run-to-run variation. All compared solves took six nonlinear iterations,
with identical residual and J*v evaluation counts for the corresponding
preconditioner setting. RSS is whole-process memory, not an allocation
profile; wall times and exact memory usage are not portable regression goals.


### Reusing augmentation-pruning workspaces

`benchmark_lgmres_pruning_workspace.f90` compares repeated image-space
orthogonalization with local temporary allocation to reusing an explicit
`lgmres_pruning_workspace`. It depends only on `src/kinds.f90` and
`src/optimize/krylov_augmentation.f90` and is not an FPM test target:

```sh
gfortran -O2 -std=f2018 -fcheck=all src/kinds.f90 \
  src/optimize/krylov_augmentation.f90 \
  benchmarks/benchmark_lgmres_pruning_workspace.f90 -o prune_bench
./prune_bench 64 20000
./prune_bench 512 4000
```

The two optional arguments are the vector length and the number of repeated
span tests; each run uses eight retained directions. Exact results depend
on allocator and CPU. Some short-vector runs gain from reuse, but for longer
vectors the orthogonalization work dominates and an individual run can be
slower. This microbenchmark measures CPU time, not peak process memory.

`lgmres_pruning_workspace_gfortran14_o2.csv` instead records single-run
process-memory comparisons for the full nonlinear solver using
`benchmark_lgmres_cache_memory.f90` and `outer_k=64`. Both snapshots were
built with the same GNU Fortran flags. The new implementation uses retired
Arnoldi basis columns as orthogonalization scratch, so retaining the pruning
workspace does not introduce another large dense basis during Arnoldi.


### Incremental LGMRES image-span screening (2026-10-11)

`benchmark_lgmres_incremental_pruning.f90` measures the difference between
rebuilding the Jacobian-image span for each retained candidate and reusing an
incrementally extended two-pass orthonormal basis during one nonlinear-iterate
filtering pass. Both paths use reusable work arrays, so allocation overhead
is not the primary source of the measured gain. The source depends only on
`src/kinds.f90` and `src/optimize/krylov_augmentation.f90`:

```sh
gfortran -O2 -std=f2018 -fcheck=all \
  -ffpe-trap=invalid,zero,overflow src/kinds.f90 \
  src/optimize/krylov_augmentation.f90 \
  benchmarks/benchmark_lgmres_incremental_pruning.f90 -o prune_stream_bench
./prune_stream_bench
```

With GNU Fortran 14.2 on one Linux machine, eight passes of 96 proposals
(length 192, 84 retained) required 0.55291 seconds with reconstruction and
0.01994 seconds with the incremental basis. These are single-run CPU-time
figures, not a whole-solver speedup guarantee or a required test threshold.

## Incremental LGMRES solution-space pruning

`benchmark_lgmres_solution_span.f90` measures the benefit of preserving an
orthogonal solution-direction span while LGMRES augmentation history is
append-only, versus rebuilding that span for each candidate. It uses 96
proposed vectors of length 192, a capacity of 112, and 12 repetitions.
It verifies identical retained directions and counts before reporting CPU
seconds; no full root solve is timed.

Build independently of FPM, from the project root:

```sh
mkdir -p _bench_solution_span
cd _bench_solution_span
gfortran -std=f2018 -O2 -ffree-line-length-none -fimplicit-none \
  -Werror=implicit-interface -fcheck=all -ffpe-trap=invalid,zero,overflow \
  -c ../src/kinds.f90 ../src/optimize/krylov_augmentation.f90
gfortran -std=f2018 -O2 -ffree-line-length-none -fimplicit-none \
  -Werror=implicit-interface -fcheck=all -ffpe-trap=invalid,zero,overflow \
  ../benchmarks/benchmark_lgmres_solution_span.f90 \
  kinds.o krylov_augmentation.o -o benchmark_lgmres_solution_span
./benchmark_lgmres_solution_span
```

`lgmres_solution_span_gfortran14_o2.csv` records a single matched run.
The speedup applies to candidate-screening work, not end-to-end solve time.
The integrated solver uses incremental screening only while the retained
history grows and caps its additional basis storage at 8 MiB; full histories
continue with the previous FIFO reconstruction.

## LGMRES full-history FIFO QR benchmark

`benchmark_lgmres_fifo_pruning.f90` isolates solution-space pruning after a
48-vector retained history fills. The new opt-in full-history QR downdate
(`outer_fifo_prune_fast=.true.` in the public root APIs) is compared with
complete span reconstruction over 200 proposed directions of length 192.
The benchmark checks that final normalized retained directions match exactly.

Build directly, without the full FPM dependency graph, from the project root:

```sh
mkdir -p _bench_fifo
cd _bench_fifo
gfortran -O2 -std=f2018 -ffree-line-length-none -fimplicit-none \
    -fcheck=all -ffpe-trap=invalid,zero,overflow -c \
    ../src/kinds.f90 ../src/optimize/krylov_augmentation.f90
gfortran -O2 -std=f2018 -ffree-line-length-none -fimplicit-none \
    -fcheck=all -ffpe-trap=invalid,zero,overflow \
    ../benchmarks/benchmark_lgmres_fifo_pruning.f90 ./*.o -o fifo_bench
./fifo_bench
```

The reference GNU Fortran 14.2 `-O2` checked-build timing is recorded in
`lgmres_fifo_gfortran14_o2.csv`. CPU time is machine dependent; this is
not an end-to-end nonlinear-solver benchmark.

## LGMRES candidate image-basis reuse benchmark

`benchmark_lgmres_image_basis_reuse.f90` times the append-only current-
Jacobian filtering pass plus one new candidate, comparing the original
basis-reconstructing insertion with the optional reuse of the image basis
built earlier during that pass. It runs both fixed and automatic image-
dependence thresholds with 192 components, a 48-direction history, and 250
repetitions. It intentionally excludes Jacobian-vector evaluations and
Arnoldi work. Compile from the project root without FPM dependencies:

```sh
mkdir -p _bench_image_basis
cd _bench_image_basis
gfortran -std=f2018 -O2 -fcheck=all -ffpe-trap=invalid,zero,overflow \
  -c ../src/kinds.f90 ../src/optimize/krylov_augmentation.f90
gfortran -std=f2018 -O2 -fcheck=all -ffpe-trap=invalid,zero,overflow \
  ../benchmarks/benchmark_lgmres_image_basis_reuse.f90 \
  kinds.o krylov_augmentation.o -o benchmark_image_basis
./benchmark_image_basis
```

`lgmres_image_basis_gfortran14_o2.csv` records one matched run. Numerical
regressions do not enforce speedup; the extra cached basis can increase
memory, up to the documented 8 MiB cap, and overall nonlinear-root speedup
may be much smaller than the isolated filtering-stage speedup.


## End-to-end LGMRES in-place image-basis benchmark

`benchmark_lgmres_inplace_image_basis.f90` runs complete matrix-free
nonlinear root solves and verifies successful convergence, a small true
residual, and stable iteration counts. It exercises automatic Jacobian-image
pruning with an analytic `jacobian_product` callback and deliberately does
not enable `outer_image_reuse`, allowing the solver to store the pruning
basis in the already allocated image workspace. Optional integer arguments
are system dimension, number of repeated solves, and the
`outer_image_basis_reuse` flag (`0` or `1`).

Build using the project dependencies with FPM if available, or compile the
root/nonlinear modules as in the focused tests. Example execution:

```sh
/usr/bin/time -f 'peak_kb=%M' ./benchmark_lgmres_inplace_image_basis 100000 4 1
```

To compare implementations, compile the *same* benchmark with the same
compiler flags against the previous release and the revised code, using
separate build directories. Run both variants with `1` for the basis-reuse
flag, and run flag `0` as a no-reuse control. The benchmark prints nonlinear
iteration and evaluation counts, a solution checksum, maximum residual,
and process CPU seconds. End-to-end timings include all residual and
Jacobian-vector work; do not extrapolate standalone pruning speedups.
`lgmres_inplace_image_basis_gfortran14_o2.csv` records two checked GNU
Fortran 14.2 Linux process runs per version and size. Peak-RSS savings were
about 3.6%-3.8% for the supplied test sizes; timing differences were small
and inconsistent. No performance threshold is part of the regression.


## Packed retained-image cache benchmark

`benchmark_lgmres_cache_memory.f90` also measures packed retained-image
pruning buffers. Compile the same source file with GNU Fortran 14.2 `-O2`,
`-fcheck=all`, and floating-point traps, against the preceding and current
SciFort source snapshots in separate directories (plus their dependencies).
Run `/usr/bin/time -f 'peak_kib=%M' ./benchmark_memory 250000 8` for
example. The two repeated Linux RSS measurements at 100,000, 250,000, and
500,000 unknowns are recorded in `lgmres_packed_images_gfortran14_o2.csv`.
The optimized build uses no additional solver flags: the saved arrays are
internal only, with no change to the API or Jacobian-product counts.


### Memory-limited raw candidate-image caching

`benchmark_lgmres_candidate_cache_budget.f90` accepts positional arguments
`n outer_k cache_limit_mib`, where `cache_limit_mib=0` is unlimited
and a positive integer bounds the separate raw Arnoldi-product cache.
It uses an analytic Jacobian-vector callback, automatic Jacobian-image
pruning and candidate-image reuse. Compile against SciFort and measure peak
process memory, for example with Linux `/usr/bin/time -f 'peak_kb=%M'`:

```sh
/usr/bin/time -f 'peak_kb=%M elapsed_s=%e' ./benchmark_lgmres_candidate_cache_budget 250000 4 0
/usr/bin/time -f 'peak_kb=%M elapsed_s=%e' ./benchmark_lgmres_candidate_cache_budget 250000 4 1
```

The 1 MiB case can invoke an additional analytic Jv callback per nonlinear
step; `nfev` counts residual evaluations and is not increased by those
analytic callback calls. A raw-cache budget does not limit total RSS.
The checked GNU Fortran 14.2 reference results are recorded in
`lgmres_candidate_cache_budget_gfortran14_o2.csv`.


## Memory-bounded LGMRES restart benchmark

Build `benchmark_krylov_restarted_lgmres.f90` with the SciFort optimization
modules and dependencies. Run the two modes independently to measure peak
resident process memory (Linux example):

```sh
/usr/bin/time -f 'peak_kib=%M' ./benchmark_lgmres_restart full 100000
/usr/bin/time -f 'peak_kib=%M' ./benchmark_lgmres_restart restart 100000
```

The benchmark uses an exactly specified diagonal, matrix-free linear residual
and analytic Jacobian-vector products; full LGMRES has `inner_maxiter=64`,
while restarted LGMRES uses `inner_restart=8` with the same total budget and
retains `outer_k=4` augmentation vectors. Both modes verify convergence,
reported residual evaluation counts, and the true final residual.

`krylov_lgmres_restart_gfortran14_o2.csv` holds two measured runs per mode
for 20,000, 100,000, and 250,000 variables. Peak resident memory includes
all process allocations, not only Arnoldi storage. Work counts and CPU times
can differ between modes, and timings are not portable speedup guarantees.


## Large, ill-conditioned nonlinear LGMRES restart benchmark (2026-10-11)

`benchmark_krylov_restart_scaling.f90` compares **full-width LGMRES**,
**fixed restarting**, and **adaptive restarting** on an explicitly
constructed matrix-free nonlinear system with a known root and a changing
Jacobian. For `n` even, the linear part is `D*Q`, where `Q` is an orthogonal
rotation between the first and second halves of the vector and `D` has
geometrically spaced positive diagonal entries from 1 to `1/condition`.
At the exact root, the Jacobian has 2-norm condition number equal to
`condition`; away from the root, the `0.04*D*(x-target)**2` contribution
changes the Jacobian. No dense `n*n` matrix is allocated.

The three policies use `inner_maxiter=48`, `outer_k=2` and, respectively,
no restart, `inner_restart=6`, and adaptive restart widths 6 through 24.
The nonlinear solver uses `maxiter=12`, `fatol=1e-8`, an analytic
Jacobian-vector callback, and optional approximate inverse preconditioning
`M*v = D**(-0.5)*v`.

With FPM and the declared external dependencies installed, copy this
standalone optional benchmark temporarily into `app/` and compile/run it:

```sh
cp benchmarks/benchmark_krylov_restart_scaling.f90 app/
fpm build --profile release
fpm run --target benchmark_krylov_restart_scaling --profile release -- 65536 1e6 adaptive yes raw
rm app/benchmark_krylov_restart_scaling.f90
```

Or, once the benchmark executable is built, compare separate **processes**
so peak resident memory can be measured independently:

```sh
/usr/bin/time -f 'peak_kb=%M' ./benchmark_krylov_restart_scaling 65536 1e6 full yes raw
/usr/bin/time -f 'peak_kb=%M' ./benchmark_krylov_restart_scaling 65536 1e6 fixed yes raw
/usr/bin/time -f 'peak_kb=%M' ./benchmark_krylov_restart_scaling 65536 1e6 adaptive yes raw
```

On Windows, run the corresponding compiled `.exe` from CMD and use a
process-memory monitor if desired; the benchmark's own CSV includes
CPU seconds and callback counts without requiring Linux `time`.
Arguments are `N CONDITION full|fixed|adaptive [yes|no] [raw|scaled]`.
The default preconditioner is `yes`; the default norm is `raw`. The optional
`scaled` norm calls the already implemented SciFort `tol_norm=` API,
using `maxval(abs(residual/D))` instead of the ordinary max norm for
stopping. It is **not** a new SciPy-compatible solver option.

`krylov_restart_scaling_gfortran14_o2.csv` contains 78 checked
GNU Fortran 14.2 `-O2` Linux process runs: sizes 64 through 65,536,
condition numbers `1e4`, `1e6`, and `1e8`, all three policies, and selected
scaled-norm cases. It reports `n`, condition number, preconditioner, norm,
restart policy, convergence flag, nonlinear iterations, residual/Jacobian
evaluations, preconditioner calls, **unscaled** max residual, **scaled** max
residual, maximum solution error against the known root, CPU time, and peak
process RSS in KiB.

For `n=65536`, `condition=1e6`, raw norm and preconditioning, the GNU
Fortran reference run gave 59,832 KiB and 164 Jv calls for full width,
17,912 KiB and 301 calls for fixed restart, and 24,004 KiB and 276 calls
for adaptive restart. All three met the requested absolute residual
tolerance. Adaptive restart was about 60% smaller in peak RSS than full
width, but it used more Jv calls. Exact timings and memory figures are
environment dependent and must **not** be used as pass/fail criteria.

**Accuracy caveat:** for `condition=1e8`, raw residual convergence may not
imply a small solution error. At `n=65536`, all three policies returned
`success=.true.` but maximum solution errors ranged from approximately
`5.3e-3` to `9.5e-3`. On the 64-variable problem, a scaled stopping norm
reduced the full-width solution error from `2.64e-5` to `2.79e-10`, at
the cost of additional work. On harder cases, the scaled norm can prevent
a premature success flag, but convergence may then require more nonlinear
iterations or a different preconditioner. No policy is universally best.

The regression `test/test_optimize_root_krylov_restart_scaling.f90` checks
the analytic derivative against centered differences, scaled and unscaled
convergence at sizes 64 and 256, accuracy against the known root, and
callback accounting. The focused offline checks do not replace a full
dependency-backed `fpm test`.

### Diagonal preconditioner for restarted Krylov root solving

`benchmark_krylov_diagonal_preconditioner.f90` compares an unpreconditioned
restarted LGMRES solve with (1) `inner_m_diag=` and (2) the existing
`inner_m=` callback on the same diagonally scaled nonlinear residual. The
Jacobian is diagonal with a known condition number of `1e8`, and a custom
residual norm measures scaled residual error rather than reporting a
premature success when the smallest rows have tiny absolute residuals.

Compile this benchmark with the SciFort modules (as with other standalone
Krylov benchmarks), then run `benchmark_krylov_diagonal_preconditioner.exe
256` or `... 1024` to print CSV with convergence, iterations, residual and
Jacobian callback counts, approximate CPU time, scaled residual, and maximum
solution error. Reference output from GNU Fortran 14.2 `-O2` runs is in
`krylov_diagonal_preconditioner_gfortran14_o2.csv`. These diagonal examples
show the benefit of a *good* preconditioner; they do not predict performance
for dense or strongly coupled Jacobians. Test-only MINPACK facade objects
used in offline builds are not distributed as part of the project.

## Automatic Krylov diagonal preconditioning

`benchmark_krylov_auto_diagonal.f90` evaluates coupled 32- and 128-variable
nonlinear systems with a diagonal scaling range of six orders of magnitude.
It compares no preconditioning, a supplied static inverse diagonal,
automatic initial-only estimation, and automatic updates at every nonlinear
step. It prints CSV with success flags, iteration counts, function and
analytic Jacobian-vector evaluation counts, elapsed process CPU time, and
maximum solution error. The stopping norm rescales residuals by the known
per-coordinate linear coefficients, making the reported solution accuracy
more meaningful than an unscaled absolute residual alone.

Run with the FPM dependencies installed, or link the benchmark against the
SciFort root library. The raw GNU Fortran 14.2 `-O2` checked-build output
is recorded in `krylov_auto_diagonal_gfortran14_o2.csv`. CPU times are
single-run measurements and are **not** regression acceptance thresholds.
Automatic estimation costs n analytic Jv calls per refresh on these
problems; with finite differences it would cost at least n or 2n
additional residual evaluations per refresh, depending on the stencil.

See `test/test_optimize_root_krylov_auto_diagonal.f90` for preflight
validation, failure handling, zero/tiny diagonal behavior, update schedules,
MINRES positive preconditioning and `nfev` accounting.

## Deterministic diagonal probes for large Krylov systems

`benchmark_krylov_diagonal_probes.f90` compares exact automatic diagonal
estimation (`inner_m_diag_probes` omitted) against 8 and 16 deterministic
Walsh probes on 128-, 512-, and 2,048-variable nonlinear tridiagonal systems
with five orders of magnitude of diagonal scaling. It prints a CSV with
nonlinear iteration and residual counts, analytic Jv calls, solution error,
and CPU time. Recorded checked GNU Fortran 14.2 `-O2` output is in
`krylov_diagonal_probes_gfortran14_o2.csv`.

The short-range couplings are chosen so 8 probes cancel them exactly. This
benchmark does **not** establish estimator accuracy on arbitrary sparse or
dense Jacobians: see the companion regression for an explicit long-range
aliasing counterexample. The timing values are single-run diagnostics.

## Seeded and verified diagonal-probe comparison (October 2026)

`benchmark_krylov_random_diagonal.f90` compares exact coordinate probing,
8 Walsh probes, 8 seeded Rademacher probes, and both methods with four
coordinate verification samples on 128-, 512-, and 2,048-variable
irregularly coupled Jacobians. It records relative inverse-preconditioner
errors, analytic Jv evaluations, residual evaluations, and CPU time.
The original absolute diagonal has four orders of magnitude of scaling;
long-range couplings deliberately alias an 8-probe Walsh pattern.
Results are in `krylov_random_diagonal_gfortran14_o2.csv`.

These timings are single-run diagnostics. Exact fallback can make a
verified method more expensive than using the exact mode immediately.
Randomized probes sometimes lower mean error while increasing maximum
error; neither convergence nor accuracy is guaranteed by a small probe
count or a passing finite sample of diagonal checks.

## Low-rank Woodbury Krylov comparison (October 2026)

`benchmark_krylov_woodbury.f90` compares no preconditioning, exact automatic
diagonal, checked rank one, checked rank two, and the new opt-in checked
**rank-two sketch** on 32- and 128-variable structured nonlinear problems
with scaling ranges of `1e2` and `1e6`. The two system variants use
`J = D + U*V^T` with one or two nonzero low-rank columns, respectively.
An analytic Jacobian-vector callback and a row-scaled stopping norm make
iteration and evaluation counts reproducible, and prevent small raw
residuals in ill-scaled rows from concealing solution error. The CSV
reports solver success, nonlinear iteration count, residual calls, analytic
Jv calls (including all setup/verification), CPU time, scaled residual,
and error against the known root.

Compile with the SciFort library and run `benchmark_krylov_woodbury` to
print all 40 comparisons. The GNU Fortran 14.2 `-O2` reference output is
`krylov_woodbury_gfortran14_o2.csv`. On the structured model, the sketched
rank-two setup requires seven Jv products, instead of `n+4`. In the
128-variable comparisons, this saves 125 setup Jv calls per refresh.
Benchmarks **do not establish** that rank-two is always faster than plain
diagonal preconditioning: the model here is purposely low rank, CPU times
are single-run noise-prone values, and an unmodeled Jacobian may trigger
the slower exact fallback. A full dependency-backed `fpm test` was not run
in the offline validation environment.

## Actual quadrature workload: CPU-timed versus work-count refresh

`benchmark_krylov_quadrature_refresh.f90` uses 128 independent local
exponential-response integrals coupled by two nonlinear global moments.
Each local integral and its derivative are evaluated with midpoint quadrature
(order 96, 640, or 2048), producing a diagonal-plus-rank-two Jacobian with
an exact *discrete* forward Jacobian-vector callback and its transpose.
OpenMP `parallel do` directives optionally evaluate the quadrature across
variables. The model itself is synthetic, but the callbacks perform real
numerical integration rather than delays or busy waiting.

The six policies are initial-only, every-iteration, interval four, measured
work with transpose weights 1 and 16, and CPU-timed work-aware refresh.
The benchmark covers two nonlinearity rates and prints 36 CSV rows per run:
`success`, iterations, `nfev`, forward and transpose Jv counts, base callback
counts, process CPU time, elapsed wall time, scaled residual, and error
against a known root. It checks the exact discrete Jacobian using centered
finite differences and the adjoint identity, and aborts on nonconvergence.
OpenMP is optional: without `-fopenmp`, directives are ignored and the
benchmark runs serially. With an OpenMP-enabled SciFort build, run the
executable with `OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE` and then with
`OMP_NUM_THREADS=4 OMP_DYNAMIC=FALSE` to compare execution. GNU Fortran 14.2
`-O2` reference data for both configurations are combined in
`krylov_quadrature_refresh_gfortran14_o2.csv` (72 rows). A checked build with
`-fcheck=all -ffpe-trap=invalid,zero,overflow` also converged all 72 cases;
all 67 existing focused Krylov tests passed in both builds using a temporary
MINPACK stub (not included in the distributed source archive). The full
FPM dependency-backed suite has not been run.

**Interpretation:** `cpu_time` measures process CPU usage, which can sum
worker-thread execution, whereas `system_clock` measures elapsed time.
Consequently the CPU-timed heuristic can make different refresh decisions
for different OpenMP thread counts even when the mathematical problem is
identical. Reference timings are single runs, susceptible to OpenMP startup,
CPU contention, and clock resolution. Work counts, validation, and
within-case correctness should be weighted more heavily than small
single-run wall-time differences. These quadrature problems are not a
substitute for validation on independent applications.

## Labeled image callback lookup

`bench_ndimage_labeled_comprehension.f90` compares grouped C-order callbacks
with a direct, independently implemented full-image scan for each requested
label (80,000 pixels, 1,200 requests). From the repository root, compile
`src/kinds.f90`, `src/ndimage_measurements.f90`, and the benchmark with a
Fortran 2018 compiler and an optimization level such as `-O2`. It prints both
elapsed CPU times, their ratio, and a reduction-agreement check. Times will
vary with compiler, architecture and input distribution.

## Batched multidimensional neighborhood callbacks (2026-10-23)

`bench_ndimage_vectorized_multi.f90` compares two separate scalar
`vectorized_filter_nd` passes with a single two-channel
`vectorized_filter_multi_nd` pass on a 160 by 128 real image with 7 by 7
windows. It checks numerical equality before printing CPU times and their
ratio. This exercises reuse of neighborhood extraction across statistics,
not a general benchmark of arbitrary callback complexity. Compile against
the core SciFort `ndimage` modules with a chosen GNU Fortran optimization
profile; the measured speedup varies by machine and compiler.

## Repeated spline interpolation using coefficient plans (2026-10-23)

`bench_ndimage_spline_plan.f90` compares 32 order-5 `grid-wrap` coordinate
queries on the same 80 x 88 image using a cached plan versus separate
`map_coordinates_nd(...,prefilter=.true.)` calls. It verifies that the outputs
are equal before printing CPU times and their ratio. Compile with
`src/kinds.f90`, `src/signal_spline.f90` and `src/signal_spline_nd.f90`.
On one GNU Fortran 14.2.0 `-O2` run, it measured 0.002469 versus
0.012985 CPU seconds (about 5.3x faster). Actual ratios vary with output point count, source
shape, order, compiler settings and number of repeated queries.
