"""Optional cross-language signal validation (requires gfortran, numpy, scipy).

Run from any directory: python tools/signal_reference/verify.py
The three drivers are intentionally outside the FPM test/ tree.
"""
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.signal import convolve, correlate, correlation_lags, detrend, find_peaks

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = Path(__file__).resolve().parent
RNG = np.random.default_rng(20261014)


def compile_drivers(output):
    sources = ["kinds", "signal_convolve", "signal_convolve2d", "signal_detrend", "signal_peaks", "signal_savgol", "signal_lfilter", "signal_sos", "signal_filtfilt", "signal_filtfilt_nd", "signal_dpss", "signal_windows_extra", "signal_spectral", "signal_fir_advanced", "signal_spectral_complex", "signal_fftconvolve", "signal_analog_design", "signal_design", "signal_elliptic", "signal_analog_sos", "signal_iir_extra", "signal_fourier_resample", "signal_resample", "signal_analog", "signal_zpk_conversion", "signal_state_space", "signal_cont2discrete", "signal_discrete_simulation", "signal_continuous_simulation", "signal_lti_frequency", "signal_partial_fraction", "signal_sequences", "signal_hilbert", "signal_waveforms", "signal_periodic", "signal_irregular", "signal_remez", "signal_czt", "signal"]
    objects = []
    for stem in sources:
        obj = output / (stem + ".o")
        subprocess.run(
            ["gfortran", "-O2", "-std=f2018", "-c", str(ROOT / "src" / (stem + ".f90")),
             "-o", str(obj)], cwd=output, check=True
        )
        objects.append(str(obj))
    exes = {}
    for stem in ["peak", "detrend", "convolve"]:
        exe = output / (stem + (".exe" if __import__("os").name == "nt" else ""))
        subprocess.run(
            ["gfortran", "-O2", "-std=f2018", "-o", str(exe),
             str(SCRIPTS / (stem + "_driver.f90")), *objects], cwd=output, check=True
        )
        exes[stem] = exe
    return exes


def run_driver(exe, inp):
    res = subprocess.run([str(exe)], input=inp, text=True, capture_output=True)
    if res.returncode:
        raise RuntimeError(res.stderr)
    return iter(res.stdout.splitlines())


def check_peaks(exe):
    cases = []
    for k in range(700):
        n = int(RNG.integers(2, 170))
        if k % 4 == 0:
            x = RNG.uniform(-5, 5, n)
        elif k % 4 == 1:
            x = RNG.integers(-4, 5, n).astype(float)
        elif k % 4 == 2:
            x = np.sin(np.arange(n) * RNG.uniform(.15, 2.0)) + RNG.normal(0, .2, n)
        else:
            x = np.round(RNG.uniform(-4, 4, n) / 2) * 2
        # Stable distance filtering in ties is not part of SciPy's contract.
        x = x + RNG.normal(0, 1.e-6, n)
        d = int(RNG.integers(1, 14))
        wlen = [0, 2, 3, 4, 7, 11, 31, 101][k % 8]
        ht = [-1e20, 0.0, 1.0][k % 3]
        threshold = [-1e20, 0., .5][(k // 3) % 3]
        pm = [0., .5, 2.][(k // 5) % 3]
        wm = [-1., .8, 1.5, 2.5][(k // 7) % 4]
        rh = [0., .2, .5, 1., 1.8][(k // 11) % 5]
        cases.append((x, d, wlen, ht, threshold, pm, wm, rh))
    inp = "".join(
        f"{len(x)} {d} {wlen} {ht} {thr} {pm} {wm} {rh}\n"
        + " ".join(f"{v:.17g}" for v in x) + "\n"
        for x, d, wlen, ht, thr, pm, wm, rh in cases
    )
    output = run_driver(exe, inp)
    keys = ["peaks", "prominences", "left_bases", "right_bases", "widths",
            "width_heights", "left_ips", "right_ips", "left_edges", "right_edges",
            "left_thresholds", "right_thresholds"]
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for idx, (x, d, wlen, ht, thr, pm, wm, rh) in enumerate(cases):
            peaks, props = find_peaks(
                x, distance=d, height=(ht, None), threshold=(thr, None),
                prominence=(pm, None), width=(wm, None), wlen=wlen or None,
                rel_height=rh, plateau_size=(1, None))
            expected = dict(props, peaks=peaks)
            for key in ["peaks", "left_bases", "right_bases", "left_edges", "right_edges",
                        "left_ips", "right_ips"]:
                expected[key] = expected[key] + 1
            status, count = map(int, next(output).split())
            assert (status, count) == (0, len(peaks)), (idx, "peak count")
            for key in keys:
                got = np.fromstring(next(output), sep=" ")
                assert np.allclose(got, expected[key], atol=1e-10, rtol=1e-10), (idx, key)
    return len(cases)


def check_detrend(exe):
    cases = []
    for k in range(450):
        n = int(RNG.integers(1, 300))
        bp = np.unique(RNG.integers(0, n + 1, int(RNG.integers(0, 8))))
        kind = k % 2
        x = RNG.uniform(-2, 2, n) + .23 * np.arange(n) + .00001 * np.arange(n) ** 2
        if k % 5 == 0:
            x = 7 + 2 * np.arange(n).astype(float)
        cases.append((x, bp, kind))
    inp = "".join(
        f"{len(x)} {len(bp)} {kind}\n"
        + " ".join(f"{v:.17g}" for v in x) + "\n"
        + (" ".join(str(int(v + 1)) for v in bp) + "\n" if len(bp) else "")
        for x, bp, kind in cases
    )
    output = run_driver(exe, inp)
    for idx, (x, bp, kind) in enumerate(cases):
        expected = detrend(x, type="constant" if kind else "linear", bp=bp)
        status, count = map(int, next(output).split())
        got = np.fromstring(next(output), sep=" ")
        assert (status, count) == (0, len(x)), idx
        assert np.allclose(got, expected, atol=3e-11, rtol=1e-11), idx
    return len(cases)


def check_convolve(exe):
    cases = []
    for k in range(900):
        n = int(RNG.integers(1, 65))
        m = int(RNG.integers(1, 65))
        mode = ["full", "same", "valid"][k % 3]
        a = RNG.normal(size=n)
        b = RNG.normal(size=m)
        if k % 2:
            a = a + 1j * RNG.normal(size=n)
            b = b + 1j * RNG.normal(size=m)
        cases.append((a, b, mode))
    lines = []
    for a, b, mode in cases:
        lines.append(f"{len(a)} {len(b)} {1 + ['full', 'same', 'valid'].index(mode)}")
        for values in (a.real, a.imag, b.real, b.imag):
            lines.append(" ".join(f"{v:.17g}" for v in values))
    output = run_driver(exe, "\n".join(lines) + "\n")
    for idx, (a, b, mode) in enumerate(cases):
        status, count = map(int, next(output).split())
        expected_conv = convolve(a, b, mode, method="direct")
        expected_corr = correlate(a, b, mode, method="direct")
        expected_lags = correlation_lags(len(a), len(b), mode)
        assert (status, count) == (0, len(expected_conv)), idx
        got = np.array([[float(v) for v in next(output).split()] for _ in range(count)])
        assert np.allclose(got[:, 0] + 1j * got[:, 1], expected_conv, atol=1e-12,
                           rtol=2e-14), (idx, "convolve")
        assert np.allclose(got[:, 2] + 1j * got[:, 3], expected_corr, atol=1e-12,
                           rtol=2e-14), (idx, "correlate")
        assert np.array_equal(got[:, 4].astype(int), expected_lags), (idx, "lags")
    return len(cases)


if __name__ == "__main__":
    with tempfile.TemporaryDirectory(prefix="scifort_signal_") as tmp:
        executables = compile_drivers(Path(tmp))
        print("peak cases:", check_peaks(executables["peak"]))
        print("detrend cases:", check_detrend(executables["detrend"]))
        print("convolution/correlation cases:", check_convolve(executables["convolve"]))
    print("all scipy reference comparisons passed")
