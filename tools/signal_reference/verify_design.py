#!/usr/bin/env python3
"""Optional independent SciPy 1.17 comparisons for SciFort signal design.

Requires NumPy, SciPy, and gfortran. Run from any current directory:
    python tools/signal_reference/verify_design.py
"""
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[2]
SEED = 20261015


def pairs(values: np.ndarray) -> str:
    return " ".join(f"{x.real:.17e} {x.imag:.17e}" for x in np.ravel(values))


def real_values(values: np.ndarray) -> str:
    return " ".join(f"{x:.17e}" for x in np.ravel(values))


def main() -> None:
    rng = np.random.default_rng(SEED)
    input_lines: list[str] = []
    cases: list[tuple] = []
    count = 0
    # Valid and invalid windowed-sinc designs. Include the non-finite Hann(2)
    # SciPy case separately as an invalid-normalization status check.
    for _ in range(480):
        mode = str(rng.choice(["lowpass", "highpass", "bandpass", "bandstop"]))
        numtaps = int(rng.integers(2, 64))
        fs = float(rng.choice([2.0, 24.0, 200.0]))
        window = str(rng.choice(["boxcar", "hann", "hamming", "blackman"]))
        if mode in ("highpass", "bandstop") and numtaps % 2 == 0:
            numtaps += 1
        if numtaps == 2 and window == "hann":
            numtaps = 4
        scale = bool(rng.integers(0, 2))
        cutoff = (
            np.array([float(rng.uniform(.1, .85))])
            if mode in ("lowpass", "highpass")
            else np.sort(rng.uniform(.07, .9, size=2))
        ) * (fs / 2)
        expected = signal.firwin(numtaps, cutoff, pass_zero=mode,
                                 window=window, fs=fs, scale=scale)
        cases.append(("FIR", expected))
        input_lines.extend([
            f"FIR {numtaps} {len(cutoff)} 0 {fs}",
            f"{mode} {window} {int(scale)}",
            real_values(cutoff),
        ])

    for _ in range(240):
        mode = str(rng.choice(["lowpass", "highpass", "bandpass", "bandstop"]))
        order = int(rng.integers(1, 9))
        fs = float(rng.choice([2, 16, 100]))
        cutoff = (
            np.array([float(rng.uniform(.12, .8))])
            if mode in ("lowpass", "highpass")
            else np.sort(rng.uniform(.15, .78, size=2))
        ) * fs / 2
        b, a = signal.butter(order, cutoff, btype=mode, fs=fs, output="ba")
        cases.append(("BUT", b, a, mode, order, cutoff, fs))
        input_lines.extend([
            f"BUT {order} {len(cutoff)} 0 {fs}",mode,real_values(cutoff),
        ])

    for _ in range(500):
        complex_input = bool(rng.integers(0, 2))
        nb = int(rng.integers(1, 9))
        na = int(rng.integers(1, 7))
        b = rng.normal(size=nb)
        a = np.r_[1, rng.normal(size=na-1)*.12]
        if complex_input:
            b = b + 1j*rng.normal(size=nb)*.25
            a = a + 1j*np.r_[0, rng.normal(size=na-1)*.08]
        fs = float(rng.choice([2, 30, 100]))
        freq = rng.uniform(-.25*fs, .85*fs, size=int(rng.integers(2,15)))
        expected = signal.freqz(b, a, worN=freq, fs=fs)[1]
        cases.append(("RESPONSE", expected))
        input_lines.extend([
            f"{'CFREQ' if complex_input else 'FREQ'} {nb} {na} {len(freq)} {fs}",
            real_values(freq), pairs(b), pairs(a),
        ])

    for _ in range(400):
        complex_input = bool(rng.integers(0, 2))
        ns = int(rng.integers(1, 7))
        sos = np.zeros((ns,6),dtype=np.complex128)
        sos[:,0] = rng.uniform(.3,1,size=ns)
        sos[:,1:3] = rng.normal(size=(ns,2))*.2
        sos[:,3] = 1
        sos[:,4:6] = rng.normal(size=(ns,2))*.2
        if complex_input:
            sos[:,0:3] += 1j*rng.normal(size=(ns,3))*.1
            sos[:,4:6] += 1j*rng.normal(size=(ns,2))*.1
        fs = float(rng.choice([2, 30, 100]))
        freq = rng.uniform(-.25*fs, .85*fs, size=int(rng.integers(2,15)))
        expected = signal.sosfreqz(sos,worN=freq,fs=fs)[1]
        cases.append(("RESPONSE", expected))
        input_lines.extend([
            f"{'CSOS' if complex_input else 'SOS'} {ns} 0 {len(freq)} {fs}",
            real_values(freq), pairs(sos),
        ])

    with tempfile.TemporaryDirectory(prefix="scifort-design-") as tmp:
        tmp_path = Path(tmp)
        exe = tmp_path / "design_driver"
        compile_cmd = [
            "gfortran", "-std=f2018", "-O2", "-J", str(tmp_path),
            "-I", str(tmp_path),
            str(ROOT / "src/kinds.f90"),
            str(ROOT / "src/signal_windows_extra.f90"),
            str(ROOT / "src/signal_dpss.f90"),
            str(ROOT / "src/signal_spectral.f90"),
            str(ROOT / "src/signal_analog_design.f90"),
            str(ROOT / "src/signal_design.f90"),
            str(ROOT / "tools/signal_reference/design_driver.f90"),
            "-o", str(exe),
        ]
        subprocess.run(compile_cmd, cwd=tmp_path, check=True)
        result = subprocess.run([str(exe)],input="\n".join(input_lines)+"\n",
                                text=True,capture_output=True,check=True)

    output_lines = iter(result.stdout.splitlines())
    worst_error = 0.0
    for case in cases:
        status, length = map(int, next(output_lines).split())
        if status != 0:
            raise AssertionError(f"case {count}: unexpected Fortran status {status}; case={case}")
        if case[0] == "FIR":
            got = np.fromstring(next(output_lines),sep=" ")
            reference = case[1]
            if len(got)!=length or len(got)!=len(reference):
                raise AssertionError(f"case {count}: FIR length mismatch")
            err = np.max(np.abs(got-reference))
            tolerance = 2e-11
        elif case[0] == "BUT":
            reference_b,reference_a,mode,order,cutoff,fs = case[1:]
            b = np.fromstring(next(output_lines),sep=" ")
            a = np.fromstring(next(output_lines),sep=" ")
            sections = (order+1)//2 if len(cutoff)==1 else order
            sos = np.array([np.fromstring(next(output_lines),sep=" ")
                            for _ in range(sections)])
            if len(b)!=length or len(a)!=len(reference_a):
                raise AssertionError(f"case {count}: Butterworth coefficient shape")
            w = np.linspace(.01,fs/2-.01,128)
            reference_sos = signal.butter(order, cutoff, btype=mode,
                                          fs=fs, output="sos")
            href = signal.sosfreqz(reference_sos,worN=w,fs=fs)[1]
            hsos = signal.sosfreqz(sos,worN=w,fs=fs)[1]
            coeff_error = max(np.max(np.abs(b-reference_b)) / max(1,np.max(np.abs(reference_b))),
                              np.max(np.abs(a-reference_a)) / max(1,np.max(np.abs(reference_a))))
            response_error = np.max(np.abs(hsos-href))
            err = max(coeff_error,response_error)
            tolerance = 2e-8
        else:
            vals = np.array([tuple(map(float,next(output_lines).split()))
                             for _ in range(length)])
            got = vals[:,0]+1j*vals[:,1]
            reference = case[1]
            err = np.max(np.abs(got-reference))
            tolerance = 1e-10 * max(1.0,float(np.max(np.abs(reference))))
        if err > tolerance or not np.isfinite(err):
            raise AssertionError(f"case {count}: error {err:.3g} > {tolerance:.3g}")
        count += 1
        worst_error = max(worst_error,float(err))
    print(f"signal design: {count} SciPy comparisons passed; worst absolute error {worst_error:.3g}")


if __name__ == "__main__":
    main()
