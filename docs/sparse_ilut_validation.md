# SciFort sparse ILUT validation — 2026-11-05

This extension adds `spilu` / `spilu_solve` to the real-binary64
`scifort_sparse_linalg` module, using canonical one-based CSR or CSC matrices.
No external sparse factorization implementation is linked or copied.

## Verification performed

The standalone sparse test runner completed **all 19 native programs** under
both `-O0 -g -fcheck=all` and `-O3` using GNU Fortran 14.2, including
`test_sparse_spilu` and the previous 18 sparse regression programs.

The final version of `compare_spilu_scipy.py` passed **200 independent
randomized problems per mode**, a total of 400 cases, using SciPy 1.17 and
NumPy references. Of these, 80 total cases used zero dropping and generous
fill and were checked for near-exact no-pivot LU factorization and solve.
The other 320 used nonzero dropping/limited fill; all 400 passed independent
Python dictionary-ILUT factor checks, factor application checks, and original
matrix GMRES solves. CSR and CSC inputs, duplicate COO coordinates, stored
zeros, multiple right-hand sides and diagonal shifts are represented.

Additional random seeds completed 120 cases in checked mode and 240 in
optimized mode during development. All previously existing native sparse
regression programs also passed in both modes. Checked-build legacy
comparisons passed: 75 sparse direct-solver cases and 75 iterative-solver
matrices (350 SciPy solver comparisons).

Test the updated code with:

```bash
python tools/sparse_reference/run_native_sparse.py --mode both
python tools/sparse_linalg_reference/compare_spilu_scipy.py --mode both --cases 200
```

## Meaning of the reference comparisons

`scipy.sparse.linalg.spilu` uses SuperLU ILUTP, with different drop rules,
permutations, optional pivoting, and storage accounting. **Its factors
are not required to match SciFort's natural-order ILUT factors.** The new
reference harness checks the independently reconstructed `L*U`, true
linear-system residuals, exact solutions when no dropping is enabled,
preconditioned GMRES, and a separate Python ILUT with the *specified*
SciFort dropping policy. It also calls SciPy `spilu` on the same numerical
inputs to verify independent factor applicability.

For very aggressive dropping (`drop_tol=0.12`, `fill_factor=1`), the
approximate solve can have a relative residual around 0.33 on these tests;
this is expected for a deliberately coarse preconditioner, **not** a direct
solve accuracy guarantee. Preconditioned GMRES still achieved its specified
true residual tolerance for every randomized test. For difficult sparse
matrices, `spilu` can fail on a zero or tiny unpivoted diagonal even if
SuperLU can pivot. Use a diagonal shift where appropriate, or the pivoted
`splu` for direct solutions.

## Known unrelated limitations

The full project-wide FPM test suite was not run in this environment because
FPM was unavailable. The pre-existing executable-stack linker warning for
internal LSMR callbacks remains; this update did not change those callbacks.
