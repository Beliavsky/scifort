# Complex Riemann zeta validation

Date: 2026-11-13. Compiler: GNU Fortran (GCC) 17.0.0 20260510
(experimental), Windows. Python references: mpmath 1.3.0 (90 decimal digits)
and SciPy 1.15.3.

## Results and scope

The new complex overload, existing real overload, and special-function umbrella
were compiled at `-O0` and `-O3`. Both builds used `-fcheck=all`,
`-ffpe-trap=invalid,zero,overflow`, and `-Werror=implicit-interface`.
The debug build also used `-g -fbacktrace`.

All seven focused programs passed in both configurations:

- `test_special_complex_zeta`: 159 assertions, including reference values,
  continuation, real-axis agreement, trivial and nontrivial zeros, conjugate
  symmetry, arrays, a functional equation, the reciprocal pole down to an
  imaginary offset of `1e-300`, NaNs/infinities, and range/overflow guards.
- `test_special_zeta_spence`: existing real Riemann/Hurwitz and Spence behavior.
- `test_special_gamma_family`: the complex log-gamma dependency.
- `test_special_integer_sequences`: an existing real-zeta consumer.
- `test_special_gamma_combinatorics`: existing Hurwitz-zeta consumers.
- `test_betanbinom_yulesimon_zipf_zipfian`: probability distributions using
  real Hurwitz zeta, including existing C ABI checks.
- `test_score`: likelihood scores using real Hurwitz zeta and its derivative.

The direct importing modules `special/integer_sequences`,
`special/gamma_combinatorics`, `stats/zipf`, and `score` also compiled separately
against the new generic interface with debug checks. The real specific retains
its existing symbol and implementation. The example in `docs/complex_zeta.md`
compiled and ran, returning `(0.79802198514627576,-0.11374430805293852)`.

No existing numerical kernel was changed; only the complex overload and its
new kernel were added. The selected regressions cover the direct imports and
new dependency. The full repository suite was not run, consistent with the
impact-based testing policy.

## Independent references

`check_complex_zeta.py` passed 248 points in each build against 90-digit
mpmath values. It also independently checked all 25 adapted Bernoulli
coefficients against `bernoulli(2*k)/factorial(2*k)`, including agreement after
rounding to binary64. The set covers positive and negative real parts, the
critical strip, the continuation boundary, near-zero and near-pole arguments,
five nontrivial zeros, deterministic random inputs, small imaginary parts
down to `1e-300`, and imaginary arguments up to both signs of 100,000.

Errors use `abs(result-reference)/(1+abs(reference))`:

| Region | Acceptance bound | Maximum observed, both builds |
| --- | --- | --- |
| `abs(aimag(s)) <= 200` | `5e-13` | `2.44414e-13` |
| Higher imaginary cases | `3e-9` | `2.33864e-10` at `s=-10+100000i` |

These sampled bounds are not a uniform guarantee throughout the supported
rectangle. Relative accuracy near zeros and phase accuracy at large imaginary
arguments have the limitations described in the guide.

SciPy comparisons were also performed. Five SciPy 1.15.3 values exceeded
the applicable mpmath bound: `1e-8i`, `1+1e-8i`, `1.0001+0.0001i`,
`1.00000001+1e-8i`, and `1.000000000001+1e-12i` (decimal strings denote their
binary64 inputs). The largest normalized native/SciPy difference was
`1.47451e-5`. The native/multiprecision bounds remained unchanged for these
cases. For all other finite SciPy values, the native/SciPy difference was
required to be at most twice the applicable bound.

## Commands

The focused builds reused the existing FPM archive for unchanged modules;
new kernels and `src/special.f90` were compiled explicitly ahead of it.
The debug object commands were:

```cmd
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/zeta-debug -Ibuild/zeta-debug -Ibuild/gfortran_2806A042DAA8BBF7 -c src/special/complex_zeta.f90 -o build/zeta-debug/complex_zeta.o
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/zeta-debug -Ibuild/zeta-debug -Ibuild/gfortran_2806A042DAA8BBF7 -c src/special/zeta.f90 -o build/zeta-debug/zeta.o
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/zeta-debug -Ibuild/zeta-debug -Ibuild/gfortran_2806A042DAA8BBF7 -c src/special.f90 -o build/zeta-debug/special.o
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/zeta-debug -Ibuild/zeta-debug -Ibuild/gfortran_2806A042DAA8BBF7 test/test_special_complex_zeta.f90 build/zeta-debug/special.o build/zeta-debug/zeta.o build/zeta-debug/complex_zeta.o build/gfortran_3C159EB802C22FB9/scifort/libscifort.a -o build/zeta-debug/test_special_complex_zeta.exe
build\zeta-debug\test_special_complex_zeta.exe
```

The same link/run command was used for the other six tests. The Zipf-family
test additionally compiled `test/test_checks.f90` and
`test/betanbinom_yulesimon_zipf_zipfian_reference.f90` before the test program.
An initial manual link omitted those two test helpers and failed; linking
them resolved the harness error, and the test passed in both modes.

The optimized commands replaced `-O0 -g -fbacktrace` with `-O3` and used
`build/zeta-release`. Runtime checks and traps remained enabled.
The direct importer compile checks used the debug object command above,
substituting each importer source and a separate `dependent_*.o` output.

The reference driver was linked with the same flags, `zeta.o`,
`complex_zeta.o`, and the existing archive. Reference checks were:

```cmd
python validation/check_complex_zeta.py build/zeta-debug/reference.exe
python validation/check_complex_zeta.py build/zeta-release/reference.exe
```

Whitespace, local documentation links, new Fortran ASCII and line lengths,
and retained license text were checked. No build artifacts or downloaded raw
upstream source are part of this change.
