# Independent SciPy reference comparisons for `envelope` and `deconvolve`

Run from the repository root after installing GNU Fortran, NumPy, and SciPy:

```bash
python reference/run_checks.py
```

The script builds **only** SciFort's four dependency-free signal modules into
`build/signal_reference/`, compiles the four Fortran data drivers, then compares
1,900 deterministic seeded cases with `scipy.signal.deconvolve` and
`scipy.signal.envelope`. The data drivers are intentionally outside `test/`
because they read their input from stdin and are not FPM standalone tests.

Fortran-native regressions are in `test/test_signal_envelope_deconvolve.f90`
and `test/test_signal_envelope_deconvolve_public.f90`. The public test needs
the complete SciFort signal layer, whereas the numerical comparison only
needs `kinds`, `signal_spectral`, `signal_envelope`, and `signal_deconvolve`.

These verification scripts depend on SciPy, but the Fortran runtime does not.
