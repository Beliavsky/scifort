#!/usr/bin/env python3
"""Build and compare SciFort's O(log n) MLS skip-ahead with SciPy 1.17.

Run from any directory with: python tools/run_signal_sequence_skip.py
Uses SciPy only for independent regression data, not as a Fortran dependency.
Large offsets are evaluated by exponentiating the *SciPy* one-clock state
transition (in Python), rather than generating billions of discarded bits.
"""
from __future__ import annotations

import subprocess
import tempfile
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[1]
RNG = np.random.default_rng(20261028)


@dataclass(frozen=True)
class Case:
    nbits: int
    offset: int
    length: int
    state: tuple[int, ...]
    taps: tuple[int, ...] = ()  # empty => SciPy default feedback polynomial


def encode(cases: list[Case]) -> str:
    lines = [str(len(cases))]
    for c in cases:
        lines.append(f"{c.nbits} {c.length} {c.offset} {len(c.taps)}")
        lines.append(' '.join(map(str, c.state)))
        if c.taps:
            lines.append(' '.join(map(str, c.taps)))
    return '\n'.join(lines) + '\n'


def make_cases() -> list[Case]:
    cases: list[Case] = []
    lengths = [0, 1, 2, 5, 17, 43, 61]
    offsets = [0, 1, 2, 7, 23, 89, 257, 987, -1, -2, -25, -127]
    for nbits in range(2, 33):
        initial = tuple(int(x) for x in RNG.integers(0, 2, size=nbits))
        if not any(initial):
            initial = (1,) + initial[1:]
        for j, offset in enumerate(offsets):
            cases.append(Case(nbits, offset, lengths[j % len(lengths)], initial))
    for j in range(270):
        nbits = int(RNG.integers(2,33))
        state = RNG.integers(0, 2, size=nbits)
        if not np.any(state):
            state[0] = 1
        taps: tuple[int, ...] = ()
        if j % 3:
            num_taps = int(RNG.integers(1,min(nbits,5)))
            taps = tuple(sorted(map(int,RNG.choice(np.arange(1,nbits),size=num_taps,replace=False))))
        offset = int(RNG.integers(-250,250))
        cases.append(Case(nbits,offset,int(RNG.integers(0,80)),tuple(map(int,state)),taps))
    # Exercise large integers (including the minimum int64). Use both custom
    # and primitive taps, and compare distinct initial states.
    big = [2**20+1, 2**40+13, 2**62+99, 2**63-1,
           -(2**20+7), -(2**40+11), -(2**62+17), -(2**63)]
    for nbits in [2, 4, 7, 11, 17, 25, 32]:
        for j, offset in enumerate(big):
            taps = () if j % 2 == 0 else (nbits-1,)
            cases.append(Case(nbits,offset,24,tuple([1]*nbits),taps))
    return cases


def to_bits(state: tuple[int, ...]) -> int:
    return sum((int(bit) != 0) << j for j, bit in enumerate(state))


def from_bits(state: int, n: int) -> np.ndarray:
    return np.array([(state >> j) & 1 for j in range(n)], dtype=np.int8)


def apply_columns(columns: tuple[int, ...], state: int) -> int:
    result = 0
    while state:
        bit = state & -state
        result ^= columns[bit.bit_length() - 1]
        state ^= bit
    return result


def inverse_columns(columns: tuple[int,...]) -> tuple[int,...]:
    """Independent bit-packed Gauss-Jordan inverse of the SciPy transition."""
    n = len(columns)
    rows = []
    for row in range(n):
        arow = sum(((col >> row) & 1) << j for j, col in enumerate(columns))
        rows.append(arow | (1 << (n + row)))
    for col in range(n):
        pivot = next((k for k in range(col,n) if (rows[k] >> col) & 1), None)
        if pivot is None:
            raise AssertionError('noninvertible LFSR transition')
        rows[col],rows[pivot] = rows[pivot],rows[col]
        for row in range(n):
            if row != col and ((rows[row] >> col) & 1):
                rows[row] ^= rows[col]
    return tuple(sum(((rows[row] >> (n+col)) & 1) << row for row in range(n))
                 for col in range(n))


@lru_cache(maxsize=128)
def scipy_columns(nbits: int, taps: tuple[int, ...]) -> tuple[int,...]:
    kwargs = {} if not taps else {'taps':np.array(taps,dtype=np.int64)}
    images = []
    for col in range(nbits):
        seed = from_bits(1<<col, nbits)
        _, next_state = signal.max_len_seq(nbits, state=seed, length=1, **kwargs)
        images.append(to_bits(tuple(map(int,next_state))))
    return tuple(images)


def scipy_jump_state(c: Case) -> np.ndarray:
    base = scipy_columns(c.nbits,c.taps)
    columns = base if c.offset >= 0 else inverse_columns(base)
    remaining = abs(c.offset)
    bits = to_bits(c.state)
    while remaining:
        if remaining & 1:
            bits = apply_columns(columns,bits)
        remaining >>= 1
        if remaining:
            columns = tuple(apply_columns(columns,x) for x in columns)
    return from_bits(bits,c.nbits)


def expected(c: Case) -> tuple[np.ndarray, np.ndarray]:
    kwargs = {} if not c.taps else {'taps':np.asarray(c.taps,dtype=np.int64)}
    state = scipy_jump_state(c)
    return signal.max_len_seq(c.nbits,state=state,length=c.length,**kwargs)


def run_driver(executable: Path, cases: list[Case]) -> None:
    completed = subprocess.run([str(executable)],input=encode(cases),text=True,
                               capture_output=True,check=True)
    lines = iter(completed.stdout.splitlines())
    for i,c in enumerate(cases):
        header = list(map(int,next(lines).split()))
        assert header == [0,c.length,c.nbits],(i,c,header)
        got_seq = np.fromstring(next(lines),sep=' ',dtype=np.int8)
        got_state = np.fromstring(next(lines),sep=' ',dtype=np.int8)
        reference_seq,reference_state = expected(c)
        if not np.array_equal(got_seq,reference_seq):
            raise AssertionError(f'case {i}: incorrect sequence, {c}')
        if not np.array_equal(got_state,reference_state):
            raise AssertionError(f'case {i}: incorrect state, {c}')
    assert next(lines,None) is None, 'extra output'
    print(f'  PASS {len(cases)} SciPy skip/rewind comparisons')


def build_and_run(temp: Path, optimized: bool, cases: list[Case]) -> None:
    mode = 'optimized' if optimized else 'checked'
    folder = temp / mode
    folder.mkdir()
    flags = ['-std=f2018','-O3','-march=native'] if optimized else [
        '-std=f2018','-O0','-g','-fcheck=all','-fbacktrace']
    common = ['gfortran',*flags,'-Wall','-Wextra','-Werror=implicit-interface',
              '-J',str(folder),'-I',str(folder)]
    source = [str(ROOT/'src/kinds.f90'),str(ROOT/'src/signal_sequences.f90')]
    tests = [ROOT/'test/test_signal_sequences.f90',
             ROOT/'test/test_signal_sequences_skip.f90']
    for idx,test in enumerate(tests):
        exe = folder/f'native_{idx}'
        subprocess.run([*common,*source,str(test),'-o',str(exe)],
                       cwd=folder,check=True,capture_output=True,text=True)
        result = subprocess.run([str(exe)],cwd=folder,
                                check=True,capture_output=True,text=True)
        print(' ',result.stdout.strip())
    exe = folder/'reference_driver'
    subprocess.run([*common,*source,
                    str(ROOT/'tools/signal_reference/sequence_skip_driver.f90'),
                    '-o',str(exe)],cwd=folder,check=True,capture_output=True,text=True)
    run_driver(exe,cases)
    print(' ',mode,'PASS')


def main() -> None:
    cases = make_cases()
    with tempfile.TemporaryDirectory(prefix='scifort-mls-skip-') as d:
        root = Path(d)
        for optimized in [False,True]:
            build_and_run(root,optimized,cases)
    print('PASS all targeted maximum-length-sequence regressions')


if __name__ == '__main__':
    main()
