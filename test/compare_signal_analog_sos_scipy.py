"""Reproducible randomized SciFort/SciPy analog IIR ZPK and SOS verification.

Run from any directory with `python test/compare_signal_analog_sos_scipy.py`.
Requires a `gfortran` on PATH and SciPy/NumPy, but not FPM.
"""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy import signal
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
FAMILIES = ("butter", "cheby1", "cheby2", "ellip")
SHAPES = ("lowpass", "highpass", "bandpass", "bandstop")
RNG = np.random.default_rng(20261016)


def main() -> None:
    cases = []
    for family in range(1, 5):
        for shape in range(1, 5):
            for order in range(1, 9):
                for _ in range(10):
                    lo = float(RNG.uniform(0.6, 8.0))
                    hi = float(lo * RNG.uniform(1.2, 3.5))
                    rp = float(RNG.uniform(0.2, 2.5))
                    rs = float(RNG.uniform(30.0, 72.0))
                    cases.append((order, family, shape, lo, hi, rp, rs))

    with tempfile.TemporaryDirectory() as tmp:
        directory = Path(tmp)
        flags = ["gfortran", "-O2", "-ffree-line-length-none", "-Werror=implicit-interface"]
        objects = []
        for name in (
            "kinds", "signal_dpss", "signal_windows_extra", "signal_spectral", "signal_analog_design", "signal_design",
            "signal_elliptic", "signal_analog_sos",
        ):
            subprocess.run(
                flags + ["-c", str(SRC / f"{name}.f90")], cwd=directory, check=True,
                capture_output=True, text=True,
            )
            objects.append(str(directory / f"{name}.o"))
        executable = directory / "analog_driver"
        subprocess.run(
            flags + [str(ROOT / "scripts" / "analog_sos_reference_driver.f90"),
                     *objects, "-o", str(executable)],
            cwd=directory, check=True, capture_output=True, text=True,
        )
        inputs = "\n".join(" ".join(str(x) for x in case) for case in cases) + "\n"
        output = subprocess.run(
            [str(executable)], input=inputs, capture_output=True,
            text=True, check=True,
        ).stdout.splitlines()

    if len(output) != 2 * len(cases):
        raise AssertionError(f"expected {2 * len(cases)} output lines, found {len(output)}")
    max_response_difference = 0.0
    max_root_difference = 0.0
    max_gain_difference = 0.0
    for j, (order, family, shape, lo, hi, rp, rs) in enumerate(cases):
        if output[2*j] == "FAILED":
            raise AssertionError(f"Fortran design rejected case {j}: {cases[j]}")
        nz_text, np_text, gain_text = output[2*j].split()
        nz, npoles, gain = int(nz_text), int(np_text), float(gain_text)
        raw = np.fromstring(output[2*j + 1], sep=" ")
        expected = 2*(nz+npoles+10)
        if raw.size != expected:
            raise AssertionError(f"case {j}: {raw.size} numeric values, expected {expected}")
        roots = raw.reshape((-1,2))
        z = roots[:nz, 0] + 1j*roots[:nz, 1]
        p = roots[nz:nz+npoles, 0] + 1j*roots[nz:nz+npoles, 1]
        response = roots[nz+npoles:,0] + 1j*roots[nz+npoles:,1]
        wn = lo if shape < 3 else [lo, hi]
        ref_z, ref_p, ref_k = signal.iirfilter(
            order, wn, btype=SHAPES[shape-1], ftype=FAMILIES[family-1],
            rp=rp, rs=rs, analog=True, output="zpk",
        )
        if nz != len(ref_z) or npoles != len(ref_p):
            raise AssertionError(f"case {j}: root count mismatch")
        for actual, reference in ((z, ref_z), (p, ref_p)):
            if actual.size:
                errors = abs(actual[:,None]-reference[None,:])
                row, col = linear_sum_assignment(errors)
                difference = max(errors[row,col] / np.maximum(1.0, abs(reference[col])))
                max_root_difference = max(max_root_difference, difference)
                if difference > 2e-7:
                    raise AssertionError(f"case {j}: ZPK root mismatch {difference:g}")
        gain_difference = abs(gain-ref_k) / max(1.0, abs(ref_k))
        max_gain_difference = max(max_gain_difference, gain_difference)
        if gain_difference > 3e-7:
            raise AssertionError(f"case {j}: gain mismatch {gain_difference:g}")
        w = np.array([0,lo*0.1,lo*0.3,lo*0.8,lo,lo*1.2,
                      hi*0.8,hi,hi*1.2,hi*10])
        ref_response = signal.freqs_zpk(ref_z, ref_p, ref_k, worN=w)[1]
        difference = max(abs(response-ref_response) / np.maximum(1.0,abs(ref_response)))
        max_response_difference = max(max_response_difference, difference)
        if difference > 3e-7:
            raise AssertionError(f"case {j}: response mismatch {difference:g}")

    print(f"PASS {len(cases)} designs, {len(cases)*10} response points")
    print(f"max scaled response difference: {max_response_difference:.3e}")
    print(f"max scaled root difference: {max_root_difference:.3e}")
    print(f"max scaled gain difference: {max_gain_difference:.3e}")


if __name__ == "__main__":
    main()
