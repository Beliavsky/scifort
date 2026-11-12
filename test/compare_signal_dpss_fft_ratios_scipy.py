"""Compare FFT/direct DPSS concentration ratios against independent SciPy 1.17.

Both Fortran and SciPy compute the full, untrimmed taper's concentration;
no SciPy code is embedded in SciFort. Invocations:
    python test/compare_signal_dpss_fft_ratios_scipy.py checked
    python test/compare_signal_dpss_fft_ratios_scipy.py optimized
"""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from scipy.signal.windows import dpss

ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / "reference" / "signal_dpss_ratio_reference_driver.f90"


def main() -> None:
    mode = sys.argv[1] if len(sys.argv) > 1 else "checked"
    if mode not in ("checked", "optimized"):
        raise SystemExit("mode should be checked or optimized")
    cases: list[tuple[int, float, int, int, int, int, int]] = []
    for n in (2, 3, 4, 7, 16, 32, 65, 127, 255, 256, 257, 509,
              1024, 2048, 4096, 8192, 16384):
        for nw in (0.25, min(1.75, n * 0.3), min(3.5, n * 0.35), n * 0.41):
            if not 0 < nw < n / 2:
                continue
            k = min(n if n > 2 else 1, 6)
            # Force the FFT for short windows and exercise default dispatch.
            for ratio_method in (0, 2):
                cases.append((n, nw, k, int(n % 2 == 0), 0, 0, ratio_method))
            if n <= 2048:
                cases.append((n, nw, k, int(n % 2 == 1), 2, 0, 1))
    rng = np.random.default_rng(20261027)
    for _ in range(64):
        n = int(rng.integers(3, 1400))
        nw = float(rng.uniform(0.13, n * 0.45))
        k = int(rng.integers(1, min(8, n) + 1))
        cases.append((n, nw, k, int(rng.integers(2)),
                      int(rng.integers(3)), 0, 2))

    with tempfile.TemporaryDirectory(prefix="scifort-dpss-ratio-") as directory:
        work = Path(directory)
        binary = work / ("ratio.exe" if os.name == "nt" else "ratio")
        flags = ["-std=f2018", "-Werror=implicit-interface"]
        flags += (["-O0", "-g", "-fcheck=all", "-fbacktrace"]
                  if mode == "checked" else ["-O3"])
        subprocess.run([os.environ.get("FC", "gfortran"), *flags,
                        str(ROOT / "src" / "kinds.f90"),
                        str(ROOT / "src" / "signal_dpss.f90"),
                        str(DRIVER), "-o", str(binary)], cwd=work, check=True)
        request = "".join(" ".join(map(str, case)) + "\n" for case in cases)
        process = subprocess.run([str(binary)], input=request, text=True,
                                 capture_output=True, check=True)
    output = iter(process.stdout.splitlines())
    worst = 0.0
    number = 0
    for n, nw, k, isym, inorm, eigensolver, ratio_method in cases:
        header = next(output).split()
        assert header == ["CASE", str(n), str(k), "0"], (n, nw, header)
        result = np.fromstring(next(output), sep=" ")
        ref, concentration = dpss(n, nw, k, sym=bool(isym),
                                  norm=(2, "approximate", "subsample")[inorm],
                                  return_ratios=True)
        concentration = np.asarray(concentration)
        assert result.shape == concentration.shape, (n, nw, k)
        difference = float(np.max(np.abs(result - concentration)))
        worst = max(worst, difference)
        if difference > 5e-10:
            raise AssertionError((n, nw, k, isym, inorm, ratio_method, difference))
        number += len(result)
    assert next(output, None) is None
    print(f"PASS {mode}: {len(cases)} DPSS SciPy comparisons, {number} ratios, "
          f"maximum absolute error={worst:.3g}")


if __name__ == "__main__":
    main()
