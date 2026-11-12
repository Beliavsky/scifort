#!/usr/bin/env python3
"""Compare SciFort analytic signals with SciPy 1.17+ across random sizes.

Requires NumPy, SciPy, and gfortran. Run from any current directory:
    python tools/signal_reference/verify_hilbert.py
"""
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[2]
RNG = np.random.default_rng(20261016)


def numbers(x: np.ndarray) -> str:
    return " ".join(f"{float(v):.17e}" for v in np.ravel(x, order="F"))


def main() -> None:
    inputs: list[str] = []
    expected: list[np.ndarray] = []
    for _ in range(460):
        nx = int(RNG.integers(0, 42))
        x = RNG.normal(size=nx)
        if nx == 0:
            n = int(RNG.integers(1, 18))
        elif RNG.random() < 0.33:
            n = -1
        else:
            n = int(RNG.integers(1, 47))
        reference = signal.hilbert(x, None if n == -1 else n)
        inputs.append(f"1 {nx} 0 {n} 0")
        if nx:
            inputs.append(numbers(x))
        expected.append(reference)

    for _ in range(320):
        nr = int(RNG.integers(0, 13))
        nc = int(RNG.integers(0, 14))
        a = RNG.normal(size=(nr, nc))
        if nr == 0 or nc == 0:
            mode = 2
            nf, ng = int(RNG.integers(1, 15)), int(RNG.integers(1, 15))
            reference = signal.hilbert2(a, (nf, ng))
        elif RNG.random() < .3:
            mode = 2
            nf, ng = -1, -1
            reference = signal.hilbert2(a)
        elif RNG.random() < .5:
            mode = 3
            nf, ng = int(RNG.integers(1, 18)), 0
            reference = signal.hilbert2(a, nf)
        else:
            mode = 2
            nf, ng = int(RNG.integers(1, 18)), int(RNG.integers(1, 18))
            reference = signal.hilbert2(a, (nf, ng))
        inputs.append(f"{mode} {nr} {nc} {nf} {ng}")
        if a.size:
            inputs.append(numbers(a))
        expected.append(reference)

    with tempfile.TemporaryDirectory(prefix="scifort-hilbert-") as tmp:
        cwd = Path(tmp)
        exe = cwd / "hilbert_driver"
        cmd = [
            "gfortran", "-std=f2018", "-O2", "-fcheck=all",
            "-Werror=implicit-interface", "-J", str(cwd), "-I", str(cwd),
            str(ROOT / "src/kinds.f90"),
            str(ROOT / "src/signal_spectral.f90"),
            str(ROOT / "src/signal_hilbert.f90"),
            str(ROOT / "tools/signal_reference/hilbert_driver.f90"),
            "-o", str(exe),
        ]
        subprocess.run(cmd, cwd=cwd, check=True, capture_output=True, text=True)
        result = subprocess.run(
            [str(exe)], input="\n".join(inputs)+"\n", capture_output=True,
            check=True, text=True, cwd=cwd,
        )

    lines = iter(result.stdout.splitlines())
    worst = 0.0
    for case, ref in enumerate(expected):
        dims = tuple(map(int, next(lines).split()))
        if len(dims) != ref.ndim+1 or dims[0] != 0 or dims[1:] != ref.shape:
            raise AssertionError(f"case {case}: wrong status/dimensions {dims} vs {ref.shape}")
        got = np.array([complex(*map(float, next(lines).split()))
                        for _ in range(ref.size)]).reshape(ref.shape,order="F")
        error = float(np.max(np.abs(got-ref))) if ref.size else 0.0
        worst = max(worst,error)
        if not np.allclose(got,ref,atol=4e-12,rtol=3e-12):
            raise AssertionError(f"case {case} differs from SciPy by {error:.3e}")

    print(f"Passed {len(expected)} randomized Hilbert comparisons with SciPy ")
    print(f"Maximum absolute error: {worst:.3e}")


if __name__ == "__main__":
    main()
