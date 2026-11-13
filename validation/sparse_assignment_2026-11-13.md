# Sparse weighted matching compatibility validation

Date: 2026-11-13. Compiler: GNU Fortran (GCC) 17.0.0 20260510
(experimental), Windows. Reference: SciPy 1.15.3.

## Scope

The existing solver already provided weighted rectangular sparse matching as
`minimum_weight_full_bipartite_matching`. This change adds a forwarding entry
point named `min_weight_full_bipartite_matching`, with required status and
optional `maximize` and `total_weight`. No numerical algorithm changed.
The original entry point and argument ordering remain supported.

The graph module has no library callers under `src/`; its consumers are tests
and reference drivers. All six graph test programs and the independent dense
assignment regression were selected. The full repository suite was not run,
following the impact-based policy in `AGENTS.md`.

## Results

- Seven focused programs passed at both `-O0` and `-O3`, with bounds/runtime
  checks and invalid, zero-division, and overflow traps enabled.
- `test_sparse_csgraph_assignment`: 22,389 assertions and 1,378 exhaustive
  cases per build. Every sparsity pattern through 3x3, including empty sides,
  was tested for minimization and maximization. Optima and feasibility agree
  with exhaustive enumeration and native dense assignment. Original/new APIs
  return identical results. Additional checks cover default/keyword options,
  parallel edges, explicit zero costs, and invalid nonfinite input.
- The existing weighted regression also covers a 2,000-by-2,000 sparse
  diagonal graph and weighted graph constructors.
- Randomized reference checks passed in checked and optimized builds. Per
  build: 500 weighted cases (458 feasible, 42 infeasible); 161 had exhaustive
  small-problem oracles. The same driver also passed 454 Laplacian cases and
  rejected 46 invalid negative-degree normalization cases.
- The updated guide example compiled and ran with debug runtime checks.
- Diff whitespace and new-source ASCII/line-length checks passed. The new
  local gap-inventory link points to the existing graph guide.

## Commands

The native checks compiled these sources separately into each build directory:
`src/kinds.f90`, `src/sparse_csgraph.f90`, and `src/optimize/assignment.f90`.
Each test was linked against the resulting three objects.

Debug flags:

```text
-O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface
```

Optimized flags:

```text
-O3 -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface
```

Each compile/link also used `-Jbuild/matching-debug -Ibuild/matching-debug`
or the corresponding `build/matching-release` directory. For example, the
debug compile/link commands were:

```cmd
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/matching-debug -Ibuild/matching-debug -c src/kinds.f90 -o build/matching-debug/kinds.o
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/matching-debug -Ibuild/matching-debug -c src/sparse_csgraph.f90 -o build/matching-debug/sparse_csgraph.o
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/matching-debug -Ibuild/matching-debug -c src/optimize/assignment.f90 -o build/matching-debug/assignment.o
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/matching-debug -Ibuild/matching-debug test/test_sparse_csgraph_assignment.f90 build/matching-debug/kinds.o build/matching-debug/sparse_csgraph.o build/matching-debug/assignment.o -o build/matching-debug/test_sparse_csgraph_assignment.exe
build\matching-debug\test_sparse_csgraph_assignment.exe
```

The same link/run sequence was used for `test_sparse_csgraph`,
`test_sparse_csgraph_advanced`, `test_sparse_csgraph_extensions`,
`test_sparse_csgraph_weighted`, `test_sparse_csgraph_yen`, and
`test_optimize_assignment`. Both modes passed.

Reference command:

```cmd
python tools/csgraph_reference/compare_csgraph_weighted_scipy.py --cases 500 --mode both
```

That script uses `-std=f2018 -ffree-line-length-none`, `-O0 -g -fcheck=all
-fbacktrace` for checked builds and `-O3` for optimized builds, with isolated
temporary module/object directories. Its duplicate-edge oracle first selects
the cheapest/dearest parallel edge; all costs are then shifted by 100 for
SciPy to retain explicit zero-cost edges. Compare feasibility and total cost,
not tied pairings.

## Compatibility limits

Native CSR input, one-based indices, and integer error statuses follow
SciFort conventions. Explicit zero edges, parallel-edge selection, and the
total-cost output are documented SciFort extensions. The original binary64
solver's extreme-range arithmetic/total-overflow limits remain; this change
does not claim full SciPy input-format, algorithm, or performance parity.
C ABI work remains deferred.
