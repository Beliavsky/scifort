# Hadamard, Leslie and DFT matrix validation — 2026-11-09

## Environment

- GNU Fortran 14.2.0, checked: `-std=f2018 -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`
- GNU Fortran 14.2.0, optimized: `-std=f2018 -O3 -march=native -fcheck=bounds`
- SciPy 1.17.0 and NumPy 2.3.5 as independent numerical references
- FPM not installed in the test environment; **full `fpm test` suite not run**

## Results

| Check | Checked | Optimized |
|---|---:|---:|
| `test_linalg_hadamard_leslie_dft` | PASS | PASS |
| Independent SciPy comparisons | 1473 / 1473 | 1473 / 1473 |
| Matrix entries compared | 1,308,441 | 1,308,441 |
| Expected invalid-input cases | 89 | 89 |
| Maximum scaled numerical difference | 1.234e-13 | 1.234e-13 |
| Existing `test_linalg_pascal` | PASS | PASS |
| Existing `test_linalg_special_matrices` | PASS | PASS |
| Existing `test_linalg_companion` | PASS | PASS |
| New example | PASS | PASS |

The comparison script uses SciPy's constructors as direct references, testing real and complex Leslie matrices, real/complex mixed types, and batched/broadcast output. It also checks integer Leslie matrices and all three Hadamard output types. DFT cases cover all three normalization modes and sizes from `0` to `128`, including non-power-of-two orders. Random coefficient magnitudes range from `1e-100` to `1e100`, and error cases require an explicit invalid-input status.

Native tests independently check exact Hadamard orthogonality up to order `256`, Leslie matrix structure, correct broadcasting of real and complex batch dimensions, DFT unitarity through order `65`, scaling consistency, invalid input statuses, and output clearing. No additional dependencies or third-party code are required.

## Reproduce

From the repository root, see [the API guide](linalg_hadamard_leslie_dft.md#running-the-example-and-independent-tests) for the complete checked-build commands. For the optimized build, compile the same sources with `-O3 -march=native -fcheck=bounds` and run the same native program and SciPy comparison harness. Run `fpm test` separately wherever FPM and the existing project dependencies are installed.
