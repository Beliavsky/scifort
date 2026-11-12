# Pascal matrix validation — 2026-11-09

## Environment

- GNU Fortran 14.2.0, checked: `-std=f2008 -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow`
- GNU Fortran 14.2.0, optimized: `-std=f2008 -O3 -march=native -fcheck=bounds`
- SciPy 1.17.0 as an independently computed comparison oracle
- FPM unavailable here; **the full `fpm test` suite has not been run**

## Results

| Check | Checked | Optimized |
|---|---:|---:|
| `test_linalg_pascal` | PASS | PASS |
| Independent SciPy comparisons | 474/474 | 474/474 |
| Floating matrix entries compared | 192,222 | 192,222 |
| Exact integer / decimal entries compared | 200,223 | 200,223 |
| Expected integer-overflow cases | 13 | 13 |
| Maximum scaled floating-point difference | 1.036e-15 | 1.036e-15 |
| Existing `test_linalg_special_matrices` | PASS | PASS |
| Existing `test_linalg_companion` | PASS | PASS |

The reference script checks real64, exact int64, and exact decimal outputs
against SciPy's `pascal` and `invpascal` (`exact=True`) for all three matrix
kinds, including order zero, small matrices, randomized orders, orders up to
120 for real64, and orders up to 80 for decimal strings. At `int64` overflow
boundaries, the script expects an explicit overflow status, rather than
silently accepting wrapped values. Native tests additionally verify exact
small-order inverses via integer matrix multiplication, shape errors, invalid
kinds, zero dimensions, output clearing, and undersized decimal buffers.

### Reproduce the standalone checks

From the repository root (GNU Fortran and Python with SciPy installed):

```sh
mkdir -p build_pascal && cd build_pascal
gfortran -std=f2008 -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow \
  -c ../src/kinds.f90 ../src/linalg.f90 ../src/linalg_pascal.f90
gfortran -std=f2008 -O0 -g -fcheck=all -o native_pascal \
  kinds.o linalg.o linalg_pascal.o ../test/test_linalg_pascal.f90
gfortran -std=f2008 -O0 -g -fcheck=all -o driver_pascal \
  kinds.o linalg.o linalg_pascal.o ../reference/linalg_pascal_driver.f90
./native_pascal
python ../tools/compare_linalg_pascal.py ./driver_pascal
```

`fpm test` remains the recommended full integration check where the Fortran
Package Manager and project dependencies are available.
