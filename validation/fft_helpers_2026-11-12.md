# FFT helper and 2D validation, 2026-11-12

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Transform/shift reference: SciPy 1.15.3. NumPy RNG seed: 20261112.

All 8,227 native checks passed in checked and optimized builds. All 3,579
external/helper cases passed in both builds. Maximum transform scaled error:
`2.931e-15`; shifts and integer helper results matched exactly.

## Native checks

Create `build\fft-helpers-native` before running these CMD commands:

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Werror=implicit-interface -Jbuild\fft-helpers-native -Ibuild\fft-helpers-native -Ibuild\gfortran_2806A042DAA8BBF7 src\fft.f90 test\test_fft_helpers.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\fft-helpers-native\test_fft_helpers.exe
build\fft-helpers-native\test_fft_helpers.exe
```

The optimized command substitutes `-O3` for `-O0 -g` and
`fft-helpers-optimized` for `fft-helpers-native`, retaining the other flags.
FPM module/archive hashes are local build artifacts.

Native coverage includes brute-force factorization searches for targets
1-2048, both real flag values, invalid/maximum-integer targets, shifts with
odd/even extents, selected/reordered/empty axes, empty arrays, rank 15,
strided input/output, input preservation, NaN/infinity preservation,
wrong output shapes, invalid axes, all 2D transform pairs and real/complex
overloads, rank-three batches, and invalid 2D axis/length/rank arguments.
The intentional strided call produces array-temporary warnings.

## Reference checks

```cmd
python validation\check_fft_helpers.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
```

The script compiles the changed source with
`validation/fft_nd_reference_driver.f90`, outside unattended FPM tests.
Reference flags are `-O0 -g -fcheck=all
-ffpe-trap=invalid,zero,overflow -fbacktrace -Werror=implicit-interface`
and `-O3 -fcheck=bounds -Werror=implicit-interface`.

SciPy comparisons cover all six 2D FFT names, real/complex overloads,
rank-four batches, singleton dimensions, default/explicit/reordered axes,
cropping/padding, all normalizations, odd/even half spectra, and amplitudes
`1e-100` and `1e100`. Shift permutations must agree exactly. Transform
errors are divided by the larger of maximum input magnitude, maximum
reference magnitude, and minimum normal binary64 value; acceptance is
`3e-12`, with no unit absolute floor.

Fast lengths intentionally follow SciFort's backend rather than SciPy's
radix choices. The independent Python oracle uses a priority queue to
generate all 5,911 positive 7-smooth integers representable by the tested
32-bit default integer. It checks small targets, 500 seeded random large
targets, prime-power boundaries, and both sides of the final representable
smooth length. Results and status codes must match exactly.

## Full suite and repository checks

Full-suite command:

```cmd
cmd /d /c "fpm test > build\fft-helpers-full-test.log 2>&1"
```

The full `fpm test` run passed with exit code 0, including the 8,227 new
helper checks, 5,743 existing standalone FFT checks, and 139 multidimensional
FFT checks. Array-temporary warnings were emitted by intentional strided
tests and existing interpolation/spatial code; there were no test failures.

`python tools\check_test_inputs.py` passed. Changed Fortran/Python files
are ASCII with no lines exceeding 132 columns. Affected local documentation
links, whitespace, provenance, the public API diff, and the changed-file list
were checked. No dependency, imported third-party implementation, C ABI,
or generated build artifact is added.
