#!/usr/bin/env python3
"""SciFort periodic signal and maximum-length sequence regression vs SciPy.

Builds a test-only Fortran driver; it adds no SciPy dependency to SciFort.
The random generator uses a fixed seed for reproducibility.
"""
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from scipy import signal

ROOT = Path(__file__).resolve().parents[2]
RNG = np.random.default_rng(20261016)


def generate_cases() -> list[tuple]:
    cases = []
    # Boundary locations, small cycles, and negative phases are deliberate.
    for j in range(800):
        op = j % 4 + 1
        n = int(RNG.integers(0, 64))
        t = RNG.uniform(-100.0, 100.0, size=n)
        if n > 0 and j % 7 == 0:
            t[: min(n, 6)] = np.pi * np.array([0, 0.5, 1, 1.5, 2, -2])[: min(n, 6)]
        choices = np.array([0.0,0.01,0.25,0.5,0.75,0.99,1.0])
        if op in (1, 3):
            width = float(RNG.choice(choices)) if j % 3 else float(RNG.random())
        else:
            width = RNG.random(n)
            if n > 0:
                width[::5] = RNG.choice(choices, size=len(width[::5]))
        func = signal.sawtooth if op <= 2 else signal.square
        expected = np.asarray(func(t,width),dtype=np.float64)
        cases.append((op, t, width, expected))

    for j in range(500):
        nbits = int(RNG.integers(2,33))
        length = int(RNG.integers(0, 180))
        initial = RNG.integers(0, 2, size=nbits,dtype=np.int8)
        if j % 13 == 0:
            initial = np.ones(nbits,dtype=np.int8)
        if not np.any(initial):
            initial[0] = 1
        taps = None
        if j % 2:
            count = int(RNG.integers(1,min(5,nbits)))
            taps = RNG.integers(1, nbits, size=count)
        seq,state = signal.max_len_seq(nbits,state=initial,length=length,taps=taps)
        cases.append((5, nbits, length, initial, taps, seq, state))
    return cases


def encode(cases: list[tuple]) -> str:
    buf = []
    for case in cases:
        op = case[0]
        if op <= 4:
            _, t, width, _ = case
            ctl = width if np.isscalar(width) else 0.0
            buf.append(f"{op} {len(t)} {ctl:.17g}")
            if len(t):
                buf.append(" ".join(f"{v:.17g}" for v in t))
                if not np.isscalar(width):
                    buf.append(" ".join(f"{v:.17g}" for v in width))
        else:
            _,nbits,length,initial,taps,*_ = case
            buf.append(f"5 {nbits} {length}")
            buf.append(" ".join(map(str, initial)))
            buf.append(str(0 if taps is None else len(taps)))
            if taps is not None:
                buf.append(" ".join(map(str,taps)))
    return "\n".join(buf) + "\n"


def check_output(cases: list[tuple], output: str) -> None:
    lines = iter(output.splitlines())
    max_error = 0.0
    for i, case in enumerate(cases):
        try:
            header = list(map(int,next(lines).split()))
            if case[0] <= 4:
                assert header == [0, len(case[1])],(i,header)
                raw = next(lines).split()
                got = np.array([float(x) for x in raw])
                ref = case[3]
                assert_allclose(got,ref,rtol=0,atol=2.0e-13,
                                err_msg=f"waveform case {i}, op={case[0]}")
                if len(ref):
                    max_error = max(max_error,float(np.max(np.abs(got-ref))))
            else:
                assert header == [0, len(case[5]),case[1]],(i,header)
                gotseq = np.fromstring(next(lines),sep=' ',dtype=int)
                gotstate = np.fromstring(next(lines),sep=' ',dtype=int)
                assert_array_equal(gotseq,case[5],err_msg=f"MLS bits case {i}")
                assert_array_equal(gotstate,case[6],err_msg=f"MLS state case {i}")
        except StopIteration as exc:
            raise AssertionError(f"missing driver output at case {i}") from exc
    assert next(lines,None) is None,"unexpected extra output"
    print(f"PASS {len(cases)} SciPy comparisons; max periodic-wave error: {max_error:.3g}")


def main() -> None:
    cases = generate_cases()
    with tempfile.TemporaryDirectory(prefix='scifort-signal-periodic-') as temp:
        temp = Path(temp)
        cmd = [
            'gfortran','-std=f2018','-O0','-g','-fcheck=all',
            '-Werror=implicit-interface',
            '-J',str(temp),'-I',str(temp),'-o',str(temp/'driver'),
            str(ROOT/'src/kinds.f90'),
            str(ROOT/'src/signal_periodic.f90'),
            str(ROOT/'src/signal_sequences.f90'),
            str(ROOT/'tools/signal_reference/periodic_sequences_driver.f90'),
        ]
        subprocess.run(cmd,check=True,capture_output=True,text=True)
        completed = subprocess.run([str(temp/'driver')],input=encode(cases),
                                   capture_output=True,text=True,check=True)
        check_output(cases,completed.stdout)


if __name__ == '__main__':
    main()
