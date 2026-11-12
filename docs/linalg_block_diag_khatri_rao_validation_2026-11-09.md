# Block-diagonal and Khatri–Rao matrix validation (2026-11-09)

## Environment

- GNU Fortran 14.2.0 (`-std=f2018`)
- Python 3.13, NumPy 2.3.5, SciPy 1.17.0
- Checked build: `-O0 -g -fcheck=all -fbacktrace`
- Optimized build: `-O3 -march=native`
- No extra dependencies for the new Fortran modules

## Results

| Test | Checked | Optimized |
|---|---:|---:|
| Randomized SciPy/NumPy comparisons | 3,200 / 3,200 | 3,200 / 3,200 |
| Expected invalid-input/shape failures (included above) | 532 | 532 |
| Largest scaled numerical difference | 1.716e-16 | 1.412e-16 |
| New native `test_linalg_block_diag_khatri_rao` | PASS | PASS |
| Existing Hadamard/Leslie/DFT native regression | PASS | PASS |
| Existing Pascal native regression | PASS | PASS |
| Existing companion/Fiedler native regression | PASS | PASS |
| Existing Helmert/Hilbert native regression | PASS | PASS |

**Total: 6,400 SciPy comparisons**, including 1,064 expected invalid
configurations across both builds. The reference harness tests randomized
block counts, zero-dimensional cases, scalar/vector/matrix blocks,
integer32/integer64/real64/complex64 combinations, mixed output types,
singleton batch broadcasting, incompatible batches, malformed caller-allocated
outputs, and all 16 `khatri_rao` kind pairs across 2D/3D inputs. It uses
SciPy solely as an external numerical reference, not as runtime Fortran code.

The independent tests use seed `921091` and 1,600 block-diagonal plus
1,600 Khatri–Rao examples per build, with normalized absolute comparison
`max(abs(expected - actual))/max(1,max(abs(expected)))`, where the maximum
is taken over the nonempty output. A threshold of `1e-12` was required.

The native tests additionally verify exact preservation of integer64 values
above `2**53`, int64-to-int32 out-of-range rejection, invalid-shape output
zeroing, zero-sized batches, empty block semantics, and examples of mixed
real/complex products. Integer comparisons use exact equality where practical.

## Reproduce (without FPM)

From the repository root, compile the four module dependencies into a
working directory `build-check` and build the reference driver:

```sh
mkdir -p build-check
FLAGS='-std=f2018 -O0 -g -fcheck=all -fbacktrace'
gfortran $FLAGS -J build-check -I build-check -c src/kinds.f90 -o build-check/kinds.o
gfortran $FLAGS -J build-check -I build-check -c src/linalg.f90 -o build-check/linalg.o
gfortran $FLAGS -J build-check -I build-check -c src/linalg_block_diag.f90 -o build-check/block.o
gfortran $FLAGS -J build-check -I build-check -c src/linalg_khatri_rao.f90 -o build-check/khatri.o
gfortran $FLAGS -J build-check -I build-check reference/linalg_block_diag_khatri_rao_driver.f90 \
    build-check/kinds.o build-check/linalg.o build-check/block.o \
    build-check/khatri.o -o build-check/driver
gfortran $FLAGS -J build-check -I build-check test/test_linalg_block_diag_khatri_rao.f90 \
    build-check/kinds.o build-check/linalg.o build-check/block.o \
    build-check/khatri.o -o build-check/test
./build-check/test
python tools/compare_linalg_block_diag_khatri_rao.py --exe ./build-check/driver \
    --blocks 1600 --khatri 1600 --seed 921091
```

Repeat in a separate directory with `FLAGS='-std=f2018 -O3 -march=native'`
for the optimized validation. The program in
[`example/linalg_block_diag_khatri_rao_example.f90`](../example/linalg_block_diag_khatri_rao_example.f90)
can be linked against these same four module objects.

The whole-project `fpm test` suite was not run, because `fpm` was not
installed in the validation environment. Existing unrelated modules were
not changed.
