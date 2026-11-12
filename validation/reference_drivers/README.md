# Input-driven reference programs

These Fortran programs produce data for independent numerical comparisons.
They are outside `test/` so FPM does not execute them as standalone tests.
Their input/output protocols, program names, and license notices are unchanged.

- `driver_linalg_*.f90` consume cases from standard input supplied by the
  corresponding `test/compare_linalg_*_scipy.py` scripts.
- `reference_signal_short_time_fft*.f90` produce reference values for the
  corresponding `test/compare_signal_short_time_fft*.py` scripts. Most consume
  case files; the tight-window producer generates its own cases.
- `reference_signal_spectral_nd.f90` and `reference_signal_timefreq_nd.f90`
  consume input/output paths and are maintained by the corresponding source
  generators in `scripts/`.
- `reference_signal_whittaker_reml.f90` consumes the case file produced by
  `test/compare_signal_whittaker_reml.py`.
- `signal_*_reference_driver.f90` consume input/output paths supplied by the
  corresponding FFT-padding, Hilbert, and tiled-convolution comparison scripts.
- `test_signal_missing_extra_reference.f90` consumes standard input from
  `test/compare_signal_missing_extra_scipy.py`.

Compile these sources against SciFort using the same compiler and flags as
the library. Comparison scripts that build their own drivers use this directory;
scripts accepting an executable still require a separately compiled driver.
For example, from cmd at the repository root:

```cmd
python test\compare_signal_fft_padding_scipy.py --checked
python test\compare_linalg_estimate_rank_scipy.py C:\path\to\rank_driver.exe
```

Keep self-checking, unattended regression programs in `test/` and run them with
`fpm test`. The console-input layout check is `python tools/check_test_inputs.py`.
