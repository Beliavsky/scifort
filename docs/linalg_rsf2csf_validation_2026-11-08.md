# Real-to-complex Schur conversion validation (2026-11-08)

Implementation: `src/linalg_schur.f90`, procedure `rsf2csf`.

The input is a **precomputed** real Schur form `T` and its real orthogonal
Schur vectors `Z`. The converter triangularizes each real 2-by-2 diagonal
block through a local complex unitary similarity and updates the complex
Schur vectors accordingly. It does not recalculate the Schur decomposition.

The independent reference is `scipy.linalg.rsf2csf` from SciPy 1.17.0.
The test matrix generator includes arbitrary real matrices, real symmetric
and skew-symmetric inputs, upper/lower triangular matrices, repeated spectra,
block matrices with conjugate eigenvalues, near-real conjugate pairs,
rectangular-free zero-order cases, and scales from 1e-260 to 1e130. Tests use
SciPy's real Schur factors as direct inputs, and separately route a subset
of real matrices through SciFort's real `schur` routine first.

For each case we verify:

1. The normalized reconstruction `Zc*Tc*Zc^H = A`.
2. Unitarity of `Zc`.
3. Exact zero entries strictly below the diagonal of `Tc`.
4. The unordered eigenvalue spectrum, using Hungarian assignment because
   ordering and phases of Schur factors are not unique.

**Important reference scaling caveat:** SciPy's `rsf2csf` on matrices
with entries around `1e-260` can lose the imaginary parts of 2-by-2 block
eigenvalues because the tiny-matrix eigensolver in the reference may
underflow or rescale poorly. To compare eigenvalues fairly, the Python
harness normalizes the *real Schur factor* before calling the SciPy
converter. The Fortran implementation internally scales its inputs and
preserves the imaginary parts in native `1e-260` tests.

Reproduce the tests:

```bash
python tools/linalg_reference/run_rsf2csf_checks.py --mode both \
    --cases 480 --max-order 64 --seed 19083
python tools/linalg_reference/run_rsf2csf_checks.py --mode both \
    --cases 240 --max-order 64 --seed 20261108
```

The harness builds the native regression test and a separate Fortran driver
with GNU Fortran (`-fcheck=all` and IEEE floating-point traps in checked
mode; `-O3` in optimized mode), and compares with SciPy. Native tests
also cover isolated real and conjugate-pair blocks, large and tiny scales,
empty matrices, malformed shapes, nonfinite inputs, invalid real Schur
structure, and real-Schur-to-complex-Schur composition.

This algorithm is independent of the existing Schur QR order-40 bound:
`rsf2csf` can convert valid precomputed real Schur factors of larger orders.

Limitations: inputs must already be *real Schur factors* and a corresponding
orthogonal basis, not an arbitrary dense matrix. The converter does not
calculate real Schur factors, batch multiple matrices, or produce a rigorous
backward-error bound. Extremely unbalanced 2-by-2 blocks, malformed block
partitions, or floating-point overflows may produce nonconvergence status.

The complete FPM suite is not run in this environment; independent native
and SciPy comparisons are used instead.
