#!/usr/bin/env python3
"""Compare SciFort Hadamard, Leslie, and DFT matrices against SciPy.

Usage: python tools/compare_linalg_hadamard_leslie_dft.py PATH_TO_DRIVER
Compile the Fortran driver using kinds.f90, linalg.f90, and the three new
linalg_{hadamard,leslie,dft}.f90 modules before running this script.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import subprocess

import numpy as np
from scipy.linalg import dft, hadamard, leslie


@dataclass
class Case:
    mode: str
    dtype: str
    n: int
    scale: str = "none"
    f: np.ndarray | None = None
    s: np.ndarray | None = None
    valid: bool = True


def generate_cases() -> list[Case]:
    rng = np.random.default_rng(20261109)
    cases: list[Case] = []
    hadamard_orders = [0, 1, 2, 3, 4, 6, 7, 8, 16, 32, 64, 128]
    for dtype in ("I", "R", "C"):
        for n in hadamard_orders:
            cases.append(Case("H", dtype, n, valid=n > 0 and n & (n - 1) == 0))
        for _ in range(75):
            n = int(rng.choice(hadamard_orders))
            cases.append(Case("H", dtype, n, valid=n > 0 and n & (n - 1) == 0))

    for scale in ("none", "sqrtn", "n"):
        for n in (0, 1, 2, 3, 4, 5, 7, 8, 13, 16, 31, 32, 47, 64, 97, 128):
            cases.append(Case("D", "C", n, scale=scale))
        for _ in range(70):
            cases.append(Case("D", "C", int(rng.integers(0, 85)), scale=scale))
    for n in (0, 1, 4, 16):
        cases.append(Case("D", "C", n, scale="bad", valid=False))

    def random_array(shape: tuple[int, ...], kind: str) -> np.ndarray:
        if kind == "I":
            return rng.integers(-100, 101, size=shape, dtype=np.int64)
        # Include small and large coefficients, and both signs.
        magnitude = float(rng.choice([1e-100, 1e-9, 1.0, 1e6, 1e100]))
        result = magnitude * rng.standard_normal(shape)
        if kind == "C":
            result = result + 1j * magnitude * rng.standard_normal(shape)
        return result

    for _ in range(950):
        dtype = str(rng.choice(["R", "C", "X", "Y", "I"]))
        n = int(rng.integers(2, 17))
        mode = "L" if dtype == "I" or rng.integers(2) == 0 else "B"
        if mode == "L":
            fshape, sshape = (n,), (n - 1,)
        else:
            batch = int(rng.integers(1, 6))
            bf = batch if rng.integers(2) else 1
            bs = batch if rng.integers(2) else 1
            if bf == 1 and bs == 1:
                bf = batch
            fshape, sshape = (bf, n), (bs, n - 1)
        fk = "C" if dtype in ("C", "Y") else dtype
        sk = "C" if dtype in ("C", "X") else dtype
        cases.append(Case(mode, dtype, n,
                          f=random_array(fshape, fk),
                          s=random_array(sshape, sk)))
    return cases


def encode_row(row: np.ndarray) -> str:
    def encode(value: object) -> str:
        if isinstance(value, (complex, np.complexfloating)):
            c = complex(value)
            return f"({c.real:.17g},{c.imag:.17g})"
        if isinstance(value, (int, np.integer)):
            return str(int(value))
        return f"{float(value):.17g}"
    return " ".join(encode(v) for v in row)


def build_payload(cases: list[Case]) -> str:
    rows = [str(len(cases))]
    for case in cases:
        bf = case.f.shape[0] if case.mode == "B" else 0
        bs = case.s.shape[0] if case.mode == "B" else 0
        rows.append(f"{case.mode} {case.dtype} {case.n} {bf} {bs} {case.scale}")
        if case.mode == "L":
            rows.append(encode_row(case.f))
            rows.append(encode_row(case.s))
        elif case.mode == "B":
            rows.extend(encode_row(row) for row in case.f)
            rows.extend(encode_row(row) for row in case.s)
    return "\n".join(rows) + "\n"


def evaluate(executable: str) -> None:
    cases = generate_cases()
    process = subprocess.run([executable], input=build_payload(cases), text=True,
                             capture_output=True, check=True)
    lines = iter(process.stdout.splitlines())
    max_scaled_error = 0.0
    entry_count = 0
    expected_failures = 0
    for case_id, case in enumerate(cases, start=1):
        status = int(next(lines))
        if not case.valid:
            assert status == 1, (case_id, case, "expected invalid shape", status)
            expected_failures += 1
            continue
        assert status == 0, (case_id, case, "unexpected failure", status)
        if case.mode == "H":
            reference = hadamard(case.n, dtype=complex if case.dtype == "C" else int)
        elif case.mode == "D":
            reference = dft(case.n, scale=None if case.scale == "none" else case.scale)
        else:
            reference = leslie(case.f, case.s)
        original_shape = reference.shape
        nentries = reference.size
        actual = np.empty(nentries, dtype=np.complex128)
        for index in range(nentries):
            fields = next(lines).split()
            assert len(fields) == 2, (case_id, index, fields)
            actual[index] = float(fields[0]) + 1j * float(fields[1])
        actual = actual.reshape(original_shape)
        if case.mode == "H":
            assert np.array_equal(actual, reference), (case_id, case)
        else:
            difference = np.max(np.abs(actual - reference)) if nentries else 0.0
            denominator = max(1.0, float(np.max(np.abs(reference)))) if nentries else 1.0
            error = float(difference) / denominator
            max_scaled_error = max(max_scaled_error, error)
            assert error < 2.5e-12, (case_id, case.mode, case.dtype, case.n, error)
        entry_count += nentries
    assert not list(lines), "Unexpected extra output from Fortran driver"
    print(f"Independent SciPy comparisons: {len(cases)}/{len(cases)} PASS")
    print(f"Matrix entries compared: {entry_count}")
    print(f"Expected invalid input detections: {expected_failures}")
    print(f"Maximum scaled difference: {max_scaled_error:.3e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", help="compiled Fortran reference driver")
    evaluate(parser.parse_args().executable)
