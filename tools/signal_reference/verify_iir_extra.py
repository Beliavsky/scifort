#!/usr/bin/env python3
"""Compare new SciFort IIR algorithms with SciPy 1.17.

Requires NumPy, SciPy and GNU Fortran. Run from any directory:
    python tools/signal_reference/verify_iir_extra.py

The comparison uses polynomial coefficients and whole-filter responses rather
than requiring identical orders for the biquad sections.
"""
from __future__ import annotations

from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[2]
RNG = np.random.default_rng(20261015)


def build(directory: Path, driver: str) -> Path:
    exe = directory / (Path(driver).stem + ".exe")
    cmd = [
        "gfortran", "-std=f2018", "-O2", "-J", str(directory),
        "-I", str(directory),
        str(ROOT / "src/kinds.f90"),
        str(ROOT / "src/signal_spectral.f90"),
        str(ROOT / "src/signal_analog_design.f90"),
            str(ROOT / "src/signal_design.f90"),
        str(ROOT / "src/signal_elliptic.f90"),
        str(ROOT / "src/signal_iir_extra.f90"),
        str(ROOT / "tools/signal_reference" / driver),
        "-o", str(exe),
    ]
    subprocess.run(cmd, cwd=directory, check=True, capture_output=True)
    return exe


def run(exe: Path, lines: list[str]) -> list[str]:
    result = subprocess.run([str(exe)], input="\n".join(lines)+"\n",
                            text=True, capture_output=True, check=True)
    return result.stdout.splitlines()


def test_filters(exe: Path) -> tuple[int, float]:
    cases: list[tuple[int, float, float, int, np.ndarray]] = []
    lines = []
    for _ in range(450):
        order = int(RNG.integers(1, 9))
        rs = float(RNG.uniform(8, 95))
        fs = float(RNG.choice([2, 100, 200]))
        mode = int(RNG.integers(1, 5))
        edges = np.sort(RNG.uniform(.14, .84, size=2))*fs/2
        cases.append((order, rs, fs, mode, edges))
        lines.append(f"{order} {rs:.17g} {fs:.17g} {mode} {edges[0]:.17g} {edges[1]:.17g}")
    for _ in range(300):
        fs = float(RNG.choice([2, 100, 1000]))
        w0 = float(RNG.uniform(.05, .9)*fs/2)
        q = float(RNG.uniform(.9, 100))
        mode = int(RNG.integers(5, 7))
        cases.append((1, q, fs, mode, np.array([w0, 0.0])))
        lines.append(f"1 {q:.17g} {fs:.17g} {mode} {w0:.17g} 0")

    rows = run(exe, lines)
    i = 0
    worst = 0.0
    for order, rs, fs, mode, edges in cases:
        status = int(rows[i].split()[0])
        i += 1
        if status != 0:
            raise AssertionError(f"unexpected filter status {status}")
        size = int(rows[i]); i += 1
        b = np.fromstring(rows[i], sep=" "); i += 1
        a = np.fromstring(rows[i], sep=" "); i += 1
        ns = int(rows[i]); i += 1
        sos = np.vstack([np.fromstring(rows[i+k], sep=" ") for k in range(ns)])
        i += ns
        if len(a) != size or len(b) != size:
            raise AssertionError("filter coefficient shape")
        if mode <= 4:
            kind = ["lowpass", "highpass", "bandpass", "bandstop"][mode-1]
            cutoff = edges[:1] if mode <= 2 else edges
            br, ar = signal.cheby2(order, rs, cutoff, btype=kind, fs=fs)
            sos_ref = signal.cheby2(order, rs, cutoff, btype=kind, fs=fs, output="sos")
            w = np.linspace(.001*fs, .499*fs, 257)
            _, actual = signal.sosfreqz(sos, worN=w, fs=fs)
            _, expected = signal.sosfreqz(sos_ref, worN=w, fs=fs)
            worst = max(worst, float(np.max(abs(actual-expected))))
            if np.max(abs(actual-expected)) > 2e-6:
                raise AssertionError(f"Chebyshev-II response, order={order}, mode={kind}, cutoff={list(cutoff)}, fs={fs}, rs={rs}, err={np.max(abs(actual-expected)):.3g}, num={sos[:2]}")
            if max(np.max(abs(b-br)), np.max(abs(a-ar))) > 2e-5:
                raise AssertionError(f"Chebyshev-II coefficients, order={order}, mode={kind}")
        else:
            if mode == 5:
                br, ar = signal.iirnotch(edges[0], rs, fs=fs)
            else:
                br, ar = signal.iirpeak(edges[0], rs, fs=fs)
            error = max(np.max(abs(b-br)), np.max(abs(a-ar)))
            worst = max(worst, float(error))
            if error > 2e-12:
                raise AssertionError(f"notch/peak coefficients mode={mode}")
    return len(cases), worst


def test_orders(exe: Path) -> tuple[int, float]:
    cases = []
    lines = []
    while len(cases) < 600:
        family = int(RNG.integers(1,4))
        fs = float(RNG.choice([2, 100, 200]))
        wp = float(RNG.uniform(.13, .77)*fs/2)
        ws = float(RNG.uniform(.13, .77)*fs/2)
        if abs(wp-ws) < .09*fs/2:
            if wp < .55*fs/2:
                ws = min(.96*fs/2, wp+.11*fs/2)
            else:
                ws = max(.04*fs/2, wp-.11*fs/2)
        gp = float(RNG.uniform(.3,3))
        gs = float(RNG.uniform(22,60))
        family_name = ["butter", "cheby1", "cheby2"][family-1]
        n, wn = getattr(signal, ["buttord", "cheb1ord", "cheb2ord"][family-1])(
            wp, ws, gp, gs, fs=fs)
        if n > 20:
            continue
        cases.append((family, wp, ws, gp, gs, fs, int(n), float(wn), family_name))
        lines.append(f"{family} {wp:.17g} {ws:.17g} {gp:.17g} {gs:.17g} {fs:.17g} 1")
    rows = run(exe, lines)
    idx = 0
    largest_frequency_error = 0.0
    for family, wp, ws, gp, gs, fs, n, wn, family_name in cases:
        status, n_actual, wn_actual = rows[idx].split(); idx += 1
        if int(status) != 0 or int(n_actual) != n:
            raise AssertionError(f"{family_name} filter order")
        frequency_error = abs(float(wn_actual)-wn)
        largest_frequency_error = max(largest_frequency_error,frequency_error)
        if frequency_error > 1e-9:
            raise AssertionError(f"{family_name} critical frequency")
        status, ns = map(int, rows[idx].split()); idx += 1
        b = np.fromstring(rows[idx], sep=" "); idx += 1
        a = np.fromstring(rows[idx], sep=" "); idx += 1
        b_ref, a_ref = signal.iirdesign(wp,ws,gp,gs,ftype=family_name,fs=fs)
        if status != 0 or max(np.max(abs(b-b_ref)),np.max(abs(a-a_ref))) > 1e-6:
            raise AssertionError(f"iirdesign {family_name} coefficients")
    return len(cases), largest_frequency_error


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="scifort-iir-") as tmp:
        directory = Path(tmp)
        filters = build(directory, "iir_extra_filter_driver.f90")
        orders = build(directory, "iir_extra_order_driver.f90")
        count1, error1 = test_filters(filters)
        count2, error2 = test_orders(orders)
    print(f"SciPy IIR extra: {count1+count2} randomized cases passed")
    print(f"Largest SOS-response/coefficient difference: {error1:.3e}")
    print(f"Largest order critical-frequency difference: {error2:.3e}")


if __name__ == "__main__":
    main()
