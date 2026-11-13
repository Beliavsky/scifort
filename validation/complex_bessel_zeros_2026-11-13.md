# Complex Bessel Y zero validation

Date: 2026-11-13. Compiler: GNU Fortran (GCC) 17.0.0 20260510
(experimental), Windows. References: SciPy 1.15.3 and mpmath 1.3.0,
with 80 decimal digits for independent root solves.

## Scope and results

The change adds complex-output overloads of `y0_zeros`, `y1_zeros`, and
`y1p_zeros`. The real specifics retain their implementations and symbols.
Only `src/special/bessel_zeros.f90` changes executable library behavior;
its sole importing library module is the `scifort_special` umbrella, which
was recompiled to expose the overloads. Existing complex Bessel and Hankel
kernels are reused without modification.

Four focused programs passed at both `-O0` and `-O3`, with bounds/runtime
checks and invalid, zero-division, and overflow traps enabled:

- `test_special_complex_bessel_zeros`: 203 assertions covering the first
  16 roots of each family, multiprecision first roots, finite outputs,
  residuals, derivative/value identities, ordering, the exceptional
  positive-real-part `Y_1'` root, prefix agreement, omitted status, and
  invalid/oversized counts.
- `test_special_bessel_zeros`: existing real roots and companions, including
  tables beyond 140 roots and the merged J/J-prime API.
- `test_special_bessel_complex`: the complex J/Y/I/K dependency.
- `test_special_hankel_complex`: the underlying complex Hankel dependency.

The reference checker generated 160 native roots per family in each build:
480 root/companion pairs. It compared 140 pairs per family with SciPy, then
independently refined 21 roots per family using mpmath (63 pairs total).
The mpmath selection includes the first 16, spread-out tail entries through
160, and every significant SciPy disagreement. Multiprecision initial seeds
come from mathematical phases and coarse exceptional-root guesses, not from
either library's computed roots. Conjugate-root residuals were also checked.

Both builds produced the same reported maxima:

| Check against mpmath | Acceptance bound | Maximum observed |
| --- | --- | --- |
| `abs(root-reference)/(1+abs(reference))` | `5e-13` | `3.24679e-15` |
| Relative companion error | `3e-11` | `1.42599e-14` |
| `abs(equation(root))*sqrt(abs(root))` | `2e-12*(1+abs(root))` | `1.48916e-13` |

Separately, the checked debug driver returned status zero and 1000 finite
table entries for each family. This maximum-count convergence stress check
took about 172 seconds for all three families combined on this machine. It
was not repeated in the optimized build. The 1000-root stress output was
not independently compared root by root with mpmath; the detailed reference
accuracy measurements above cover the selected roots through 160.

The new guide example compiled and ran with debug checks, returning the
first `Y_0` root and actual derivative:

```text
(-2.4030166319937760, 0.53988231304266410)
(-0.10074768928855733, 0.88196771009944841)
```

The full repository suite was not run. The selected programs cover the
changed module, its public umbrella interface, and its numerical dependencies,
consistent with the impact-based test policy.

## Observed SciPy discrepancies

- The first complex `y1p_zeros` root is approximately
  `+0.576785128606845+0.903984792166411i`. It has positive real part despite
  SciPy's documented description. SciFort follows the numerical root sequence.
- SciPy's complex `y0_zeros` companion equals `+Y_1`, whereas the documented
  actual derivative is `-Y_1`. SciFort returns the derivative, consistent with
  its existing real overload. The checker negates this SciPy companion.
- On this SciPy version, requesting 141 or 160 roots produced only 140 finite
  roots for all three families. The reference checker explicitly limits the
  SciPy comparison to 140 and uses mpmath for selected native tail roots.
- Among those 420 SciPy pairs, 16 exceeded either `5e-12` normalized root
  difference or `2e-10` relative companion difference. Each such entry was
  independently checked with mpmath, which agreed with SciFort within the
  stricter native bounds above. No native tolerance was relaxed to match SciPy.
  The maximum normalized native/SciPy root difference was `3.67877e-6`; the
  maximum relative companion difference was `2.35792e-6`.
- For example, the 140th `Y_0` root from SciFort is approximately
  `-439.037858051447+0.549305788115904i`, agreeing with mpmath to rounding
  accuracy. SciPy returned approximately
  `-439.036267449817+0.549305782962784i`, an absolute difference of `0.00159`.

These are observations about SciPy 1.15.3, not claims about every version.

## Commands

Both configurations used:

```text
-fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface
```

Debug additionally used `-O0 -g -fbacktrace`; optimized used `-O3`.
Debug module search/output flags were:

```text
-Jbuild/bessel-zeros-debug -Ibuild/bessel-zeros-debug -Ibuild/zeta-debug -Ibuild/gfortran_2806A042DAA8BBF7
```

The `zeta-debug` directory supplies the current generic zeta interface from
the preceding change; the existing FPM archive supplies unchanged modules.
Debug commands (CMD-compatible spelling) included:

```cmd
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/bessel-zeros-debug -Ibuild/bessel-zeros-debug -Ibuild/zeta-debug -Ibuild/gfortran_2806A042DAA8BBF7 -c src/special/bessel_zeros.f90 -o build/bessel-zeros-debug/bessel_zeros.o
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/bessel-zeros-debug -Ibuild/bessel-zeros-debug -Ibuild/zeta-debug -Ibuild/gfortran_2806A042DAA8BBF7 -c src/special.f90 -o build/bessel-zeros-debug/special.o
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/bessel-zeros-debug -Ibuild/bessel-zeros-debug -Ibuild/zeta-debug -Ibuild/gfortran_2806A042DAA8BBF7 test/test_special_complex_bessel_zeros.f90 build/bessel-zeros-debug/special.o build/bessel-zeros-debug/bessel_zeros.o build/zeta-debug/zeta.o build/zeta-debug/complex_zeta.o build/gfortran_3C159EB802C22FB9/scifort/libscifort.a -o build/bessel-zeros-debug/test_special_complex_bessel_zeros.exe
build\bessel-zeros-debug\test_special_complex_bessel_zeros.exe
```

The same link/run command was used for each of the other three test names.
Optimized commands used `-O3` and `build/bessel-zeros-release` for the changed
objects, executable outputs, and first module search directory; unchanged
dependency objects retained the paths above. Each reference driver was linked
with its corresponding `bessel_zeros.o` and the existing FPM archive.

```cmd
python validation/check_complex_bessel_zeros.py build/bessel-zeros-debug/reference.exe
python validation/check_complex_bessel_zeros.py build/bessel-zeros-release/reference.exe
```

The maximum-count stress check invoked the debug reference driver once per
family with inputs `0 1000`, `1 1000`, and `2 1000`.

Diff whitespace, new-source ASCII/line lengths, documentation links, and
the guide example were checked. No build artifacts or new third-party
implementation code are included.
