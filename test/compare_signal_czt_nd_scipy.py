"""Independent SciPy 1.17 numerical comparisons for assumed-rank CZT/zoom FFT.

Run from any directory:
    python test/compare_signal_czt_nd_scipy.py checked
    python test/compare_signal_czt_nd_scipy.py optimized

The reference Fortran driver uses rank-three lanes; the native regression
separately tests ranks 1-15, strided arrays, and transactional failures.
"""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from scipy.signal import czt, zoom_fft

ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / "reference" / "signal_czt_nd_reference_driver.f90"


def cases() -> list[tuple[int, int, int, int, int]]:
    result = []
    for kind in range(2):
        for n in (1, 2, 3, 7, 8, 17, 31, 64, 97):
            for m in sorted({1, 2, 5, n, n+4}):
                for ctype in (0, 1):
                    for variant in range(3):
                        if kind == 0 and variant == 2 and n > 31:
                            continue
                        # SciPy 1.17 raises ZeroDivisionError for endpoint=True, m=1.
                        if kind == 1 and variant == 1 and m == 1:
                            continue
                        result.append((kind, n, m, ctype, variant))
    return result


def expected(kind: int, n: int, m: int, ctype: int, variant: int) -> np.ndarray:
    i = np.arange(1, 3)[:, None, None]
    j = np.arange(1, n+1)[None, :, None]
    k = np.arange(1, 4)[None, None, :]
    real = np.cos(0.21*(i+3*j+7*k)) + 0.11*(i-j+2*k)
    x = real + 1j*0.35*np.sin(0.11*(2*i+j+5*k)) if ctype else real
    if kind == 0:
        if variant == 0:
            return czt(x, m=m, axis=1)
        if variant == 1:
            w = np.exp(-2j*np.pi/(m+4))
            a = np.exp(0.17j)
            return czt(x, m=m, w=w, a=a, axis=1)
        return czt(x, m=m, w=0.985-0.009j, a=1.016+0.035j, axis=1)
    if variant == 0:
        return zoom_fft(x, [-1.2, 3.1], m=m, fs=12.0, axis=1)
    if variant == 1:
        return zoom_fft(x, [-0.8, 1.75], m=m, fs=7.0, endpoint=True, axis=1)
    return zoom_fft(x, 2.8, m=m, fs=14.0, axis=1)


def main() -> None:
    mode = sys.argv[1] if len(sys.argv) > 1 else "checked"
    if mode not in ("checked", "optimized"):
        raise SystemExit("mode must be checked or optimized")
    options = cases()
    with tempfile.TemporaryDirectory(prefix="scifort-czt-nd-") as tmp:
        work = Path(tmp)
        binary = work / ("reference.exe" if os.name == "nt" else "reference")
        flags = ["-std=f2018", "-Werror=implicit-interface"]
        flags += (["-O0", "-g", "-fcheck=all", "-fbacktrace"] if mode == "checked" else ["-O3"])
        cmd = [os.getenv("FC", "gfortran"), *flags,
               str(ROOT / "src" / "kinds.f90"),
               str(ROOT / "src" / "signal_windows_extra.f90"),
               str(ROOT / "src" / "signal_dpss.f90"),
               str(ROOT / "src" / "signal_spectral.f90"),
               str(ROOT / "src" / "signal_czt.f90"),
               str(ROOT / "src" / "signal_czt_nd.f90"),
               str(DRIVER), "-o", str(binary)]
        subprocess.run(cmd, cwd=work, check=True, capture_output=True)
        request = "".join(" ".join(map(str, c)) + "\n" for c in options)
        process = subprocess.run([str(binary)], cwd=work, input=request,
                                 text=True, capture_output=True, check=True)
    lines = iter(process.stdout.splitlines())
    worst = 0.0
    total = 0
    for kind, n, m, ctype, variant in options:
        header = next(lines).split()
        assert header == ["CASE", str(kind), str(n), str(m), str(ctype), str(variant), "0"], header
        values = np.empty(6*m, dtype=np.complex128)
        for idx in range(6*m):
            a,b = map(float, next(lines).split())
            values[idx] = a + 1j*b
        actual = values.reshape((2,m,3), order="F")
        reference = expected(kind,n,m,ctype,variant)
        error = np.max(np.abs(actual-reference))
        scaled = error/max(1.0, float(np.max(np.abs(reference))))
        if scaled > 5e-10:
            raise AssertionError((kind,n,m,ctype,variant,error,scaled))
        worst = max(worst, scaled)
        total += values.size
    assert next(lines, None) is None
    print(f"PASS {mode}: {len(options)} SciPy 1.17 CZT/zoom FFT comparisons, "
          f"{total} coefficients; max magnitude-scaled error={worst:.3g}")


if __name__ == "__main__":
    main()
