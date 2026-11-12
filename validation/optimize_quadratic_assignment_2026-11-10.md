# Quadratic assignment validation, 2026-11-10

New SciFort module: `src/optimize/quadratic_assignment.f90`.
Compiler: GNU Fortran 14.2.0. Comparison library: SciPy 1.17.0.

Commands run from the extracted repository:

```sh
python validation/compare_optimize_quadratic_assignment_scipy.py
python validation/compare_optimize_quadratic_assignment_scipy.py --optimize
```

Both modes passed 5,156 independent checks across 288 solves:
144 FAQ and 144 2-opt cases. All permutations were valid, retained all fixed
matches, and had exactly the independently recomputed objectives. All 2-opt
results were locally optimal with respect to allowed pairwise swaps.

| Comparison | Checked build | `-O3` build |
|---|---:|---:|
| FAQ cases equal to SciPy objective | 142/144 | 141/144 |
| FAQ cases with different objective | 2/144 | 3/144 |
| 2-opt cases better than SciPy stochastic run | 41 | 41 |
| 2-opt cases equal to SciPy stochastic run | 80 | 80 |
| 2-opt cases worse than SciPy stochastic run | 23 | 23 |
| Independent assertions | 5,156 | 5,156 |

Quality comparisons are not exact-solver claims and cannot establish global
optimality. SciPy 2-opt requires sorted `partial_match` pairs; the independent
comparison sorts those pairs on its side while the Fortran API accepts them
in any order. SciFort 2-opt uses a deterministic initial permutation.

Native `test/test_optimize_quadratic_assignment.f90`: **469 assertions** in
both checked and optimized builds. The SciPy documentation 4-node example
produced FAQ objective 178 and 2-opt objective 176 in both builds.

Prior module smoke/regression tests rerun with GNU Fortran runtime checking:

- `test_optimize_assignment`: PASS
- `optimize_direct_reference`: 109 assertions passed
- `optimize_shgo_reference`: 61 assertions passed
- `optimize_shgo_multiconstraints_reference`: 40 assertions passed

FPM is not installed in this environment, so the full `fpm test` target,
including the external SLSQP dependency from prior SHGO work, was not run.
