#!/usr/bin/env python3
"""Compile two independent Fortran drivers and verify them against SciPy.

Run at the project root: python tools/signal_reference/verify_cwt_place_poles_scipy.py
Requires GNU Fortran, NumPy and SciPy 1.17 (or compatible versions).
"""
from pathlib import Path
from tempfile import TemporaryDirectory
import subprocess

import numpy as np
from scipy.signal import find_peaks_cwt, place_poles

ROOT = Path(__file__).resolve().parents[2]


def build_driver(work, module, driver, name):
    exe = work / name
    sources = [str(ROOT / "src/kinds.f90")]
    if module == "signal_place_poles":
        sources.append(str(ROOT / "src/linalg.f90"))
        sources.append(str(ROOT / "src/signal_place_poles_mimo.f90"))
    sources.append(str(ROOT / f"src/{module}.f90"))
    cmd = ["gfortran", "-std=f2018", "-O1", "-fcheck=all", "-J", str(work),
           "-I", str(work), *sources,
           str(ROOT / f"tools/signal_reference/{driver}.f90"), "-o", str(exe)]
    subprocess.run(cmd, cwd=work, check=True, capture_output=True)
    return exe


def wavelet_cases():
    rng = np.random.default_rng(4289)
    cases = []
    for i in range(300):
        n = int(rng.integers(20, 115))
        nw = int(rng.integers(2, 18))
        t = i % 4
        if t == 0:
            widths = np.arange(1, nw+1, dtype=float)
        elif t == 1:
            widths = np.linspace(.8, 9, nw)
        elif t == 2:
            widths = rng.uniform(.7, 13., nw).round(1)
        else:
            widths = rng.choice([1.5, 2.5, 3.5, 5.5, 8.5], nw).astype(float)
        x = rng.normal(size=n) * rng.uniform(.1, 2.)
        for _ in range(int(rng.integers(0, 5))):
            z = rng.uniform(0, n)
            x += rng.uniform(.1, 5)*np.exp(-.5*((np.arange(n)-z)/rng.uniform(.9, 9))**2)
        ml = int(rng.integers(1, nw+1))
        win = int(rng.integers(3, 16))
        noise = rng.uniform(0, 95)
        gap = float(rng.integers(1, 5))
        snr = float(rng.uniform(.3, 2.5))
        custom = i % 2 == 0
        dist = rng.uniform(.1, 5., size=nw) if custom else widths/4
        cases.append((x, widths, gap, snr, noise, win, ml, custom, dist))
    return cases


def check_wavelets(exe):
    lines = []
    expected = []
    for x, w, gap, snr, noise, win, ml, custom, dist in wavelet_cases():
        lines.append(f"{len(x)} {len(w)} {ml} {win} {int(custom)}")
        lines.append(" ".join(f"{v:.17g}" for v in x))
        lines.append(" ".join(f"{v:.17g}" for v in w))
        lines.append(f"{gap:.17g} {snr:.17g} {noise:.17g}")
        if custom:
            lines.append(" ".join(f"{v:.17g}" for v in dist))
        p = find_peaks_cwt(x, w, max_distances=dist if custom else None,
                           gap_thresh=gap, min_length=ml, min_snr=snr,
                           noise_perc=noise, window_size=win)
        expected.append(np.asarray(p, dtype=int)+1)
    output = subprocess.run([str(exe)], input="\n".join(lines)+"\n",
                            text=True, check=True, capture_output=True).stdout.splitlines()
    if len(output) != len(expected):
        raise AssertionError("wrong number of CWT responses")
    for i, (line, ref) in enumerate(zip(output, expected)):
        v = np.fromstring(line, sep=" ", dtype=int)
        if v.size < 2 or v[0] != 0 or v[1] != len(ref) or not np.array_equal(v[2:], ref):
            raise AssertionError(f"CWT case {i}: Fortran {v}, SciPy {ref}")
    print(f"CWT peaks: {len(expected)}/{len(expected)} exact SciPy matches")


def placement_cases():
    rng = np.random.default_rng(8402)
    records = []
    for i in range(280):
        n = int(rng.integers(1, 7))
        m = 1 if i % 2 == 0 else n + int(rng.integers(0, 3))
        p = -(rng.random(n)*2 + .25)
        if n >= 2 and i % 3 == 0:
            p = p.astype(complex)
            p[:2] = -(rng.random()+.5) + 1j * np.array([1, -1])*rng.uniform(.4, 2.)
            p[1] = p[0].conjugate()
        if m == 1:
            a = np.zeros((n, n))
            for j in range(n-1):
                a[j, j+1] = 1
            a[-1, :] = rng.uniform(-1, 1, size=n)
            b = np.zeros((n, 1))
            b[-1, 0] = 1
            q, _ = np.linalg.qr(rng.normal(size=(n, n)))
            a = q @ a @ q.T
            b = q @ b
        else:
            a = rng.uniform(-1.3, 1.3, size=(n, n))
            b = rng.uniform(-1, 1, size=(n, m))
        # Use SciPy as an independent verification of system feasibility.
        place_poles(a, b, p)
        records.append((a, b, p))
    return records


def check_poles(exe):
    lines = []
    cases = placement_cases()
    for a, b, p in cases:
        lines.append(f"{a.shape[0]} {b.shape[1]}")
        lines.append(" ".join(f"{v:.17g}" for v in a.flatten(order="F")))
        lines.append(" ".join(f"{v:.17g}" for v in b.flatten(order="F")))
        lines.append(" ".join(f"({v.real:.17g},{v.imag:.17g})" for v in p))
    output = subprocess.run([str(exe)], input="\n".join(lines)+"\n",
                            text=True, check=True, capture_output=True).stdout.splitlines()
    if len(output) != len(cases):
        raise AssertionError("wrong number of pole-placement responses")
    max_error = 0.
    for i, (line, (a, b, p)) in enumerate(zip(output, cases)):
        v = np.fromstring(line, sep=" ")
        n, m = b.shape
        if len(v) != 1+n*m or int(v[0]) != 0:
            raise AssertionError(f"pole placement case {i}: status={v}")
        k = v[1:].reshape((m, n), order="F")
        computed = np.linalg.eigvals(a-b@k)
        error = np.max(np.abs(np.sort_complex(computed)-np.sort_complex(p)))
        max_error = max(error, max_error)
        if error > 1e-6:
            raise AssertionError(f"pole placement case {i} has eigenvalue error {error}")
    print(f"Pole placement: {len(cases)}/{len(cases)} valid systems, "
          f"maximum eigenvalue error {max_error:.3g}")


def main():
    with TemporaryDirectory(prefix="scifort_cwt_poles_") as dirname:
        work = Path(dirname)
        cwt = build_driver(work, "signal_cwt_peaks", "cwt_peaks_driver", "cwt_check")
        poles = build_driver(work, "signal_place_poles", "place_poles_driver", "poles_check")
        check_wavelets(cwt)
        check_poles(poles)


if __name__ == "__main__":
    main()
