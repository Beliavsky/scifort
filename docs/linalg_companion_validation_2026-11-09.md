# Companion and Fiedler constructors — validation, 2026-11-09

## Environment

- GNU Fortran 14.2.0 (both checked and optimized `-O3 -march=native` builds)
- SciPy 1.17.0 (independent numerical reference), NumPy 2.3.5
- Native test: `test/test_linalg_companion.f90`
- Reference driver: `reference/linalg_companion_driver.f90`
- Python comparison: `tools/compare_linalg_companion.py`

## Results

| Validation | Checked (`-O0 -g -fcheck=all`) | Optimized (`-O3 -march=native`) |
|---|---:|---:|
| Native `test_linalg_companion` | PASS | PASS |
| Direct SciPy matrix comparisons | 1,710 / 1,710 | 1,710 / 1,710 |
| Valid matrices | 1,566 | 1,566 |
| Invalid leading coefficients/degrees | 144 | 144 |
| Polynomial-root/eigenvalue checks | 58 | 58 |
| Maximum relative matrix Frobenius difference | 1.73e-16 | 2.11e-16 |
| Previously existing dense LU test | PASS | PASS |
| Previously existing indefinite LDL test | PASS | PASS |

The reference cases cover degrees 0–15 plus 17, 18, 23, 24, 31 and 32,
real/complex/integer coefficients, extremely small or large coefficient
scales, zero leading coefficients and both parities. Root tests compare
matrix eigenvalues to `numpy.roots` with a minimum-cost permutation match.
The native tests also cover SciPy's documented matrices, one-dimensional
and batched inputs, malformed dimensions, empty arrays, atomic batched
failure and complex absolute-value distances.

The extra SciPy comparisons compare **matrix entries**, not just resulting
eigenvalues. The matrix errors are measured relative to
`max(1, norm(reference, 'fro'))`. A good companion matrix can still be
ill-conditioned for root finding; these checks do not imply accurate roots
for arbitrary polynomials or floating-point dynamic ranges.

The complete `fpm test` suite was **not** run here because FPM is not
installed in this environment. Only the new test and two existing dense
linear-algebra regression programs were compiled and run manually.

### Reproduction without FPM

From the repository root, after creating a writable build directory:

```bash
mkdir -p /tmp/scifort-companion
cd /tmp/scifort-companion
gfortran -O0 -g -fcheck=all -c /path/to/scifort/src/kinds.f90
gfortran -O0 -g -fcheck=all -c /path/to/scifort/src/linalg.f90
gfortran -O0 -g -fcheck=all -c /path/to/scifort/src/linalg_companion.f90
gfortran -O0 -g -fcheck=all -c /path/to/scifort/src/linalg_fiedler.f90
gfortran -O0 -g -fcheck=all /path/to/scifort/test/test_linalg_companion.f90 \
    kinds.o linalg.o linalg_companion.o linalg_fiedler.o -o test_companion
./test_companion
gfortran -O0 -g -fcheck=all /path/to/scifort/reference/linalg_companion_driver.f90 \
    kinds.o linalg.o linalg_companion.o linalg_fiedler.o -o companion_driver
python /path/to/scifort/tools/compare_linalg_companion.py ./companion_driver
```
