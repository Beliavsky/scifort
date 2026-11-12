#!/usr/bin/env python3
"""Check SciFort's rectangular CSR constructors and operations against SciPy.

Builds an independent Fortran reference driver and sends it randomized COO
matrices, including duplicate coordinates and explicitly stored zeros.
No SciPy sources are imported into the Fortran implementation.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[2]
FIELDS = ("A", "B", "C", "T", "P", "S", "H", "V")


def sample(rng: np.random.Generator, nr: int, nc: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    m = int(rng.integers(0, max(nr * nc * 2, 1))) if nr * nc else 0
    rows = rng.integers(0, nr, size=m, dtype=np.int32)
    cols = rng.integers(0, nc, size=m, dtype=np.int32)
    vals = rng.integers(-6, 7, size=m).astype(float) / 4.0
    # Some cases deliberately have repeated coordinates and explicit zeros.
    if m:
        vals[rng.random(m) < 0.15] = 0.0
    return rows, cols, vals


def scipy_matrix(shape: tuple[int, int], sample_data: tuple[np.ndarray, np.ndarray, np.ndarray]):
    rows, cols, vals = sample_data
    result = sparse.coo_matrix((vals, (rows, cols)), shape=shape).tocsr()
    result.sum_duplicates()
    result.sort_indices()
    return result


def format_cases(rng: np.random.Generator, count: int):
    input_lines = [str(count)]
    cases = []
    for case_id in range(count):
        # Non-square, empty, and repeated-coordinate graphs all occur.
        m, k, n = [int(x) for x in rng.integers(0, 15, size=3)]
        if case_id % 13 == 0:
            m, k, n = (0, 0, 0)
        elif case_id % 13 == 1:
            m, k, n = (13, 0, 8)
        elif case_id % 13 == 2:
            m, k, n = (0, 5, 4)
        triples = [sample(rng, m, k), sample(rng, k, n), sample(rng, m, k)]
        a, b, c = [scipy_matrix(shape, data) for shape, data in
                   zip(((m, k), (k, n), (m, k)), triples)]
        x = (np.arange(1, k + 1) % 7 - 3) / 3.0
        expected = {
            "A": a.toarray(), "B": b.toarray(), "C": c.toarray(),
            "T": a.T.toarray(), "P": (a @ b).toarray(),
            "S": (2.0 * a - 3.0 * c).toarray(),
            "H": a.multiply(c).toarray(), "V": (a @ x)[:, None],
        }
        cases.append(expected)
        input_lines.append(f"{m} {k} {n}")
        for rows, cols, vals in triples:
            input_lines.append(str(len(vals)))
            input_lines.extend(f"{i+1} {j+1} {v:.17g}" for i, j, v in zip(rows, cols, vals))
    return "\n".join(input_lines) + "\n", cases


def parse_and_check(output: str, cases: list[dict[str, np.ndarray]]) -> None:
    lines = iter(output.splitlines())
    for case_id, expected in enumerate(cases, start=1):
        heading = next(lines)
        assert heading == f"CASE {case_id}", (case_id, heading)
        for label in FIELDS:
            hdr = next(lines).split()
            assert hdr[0] == label, (case_id, label, hdr)
            n_rows, n_cols = map(int, hdr[1:])
            arr = np.empty((n_rows, n_cols))
            for i in range(n_rows):
                row = next(lines).strip()
                numbers = [float(item) for item in row.split()]
                assert len(numbers) == n_cols, (case_id, label, numbers)
                arr[i, :] = numbers
            ref = expected[label]
            assert ref.shape == arr.shape, (case_id, label, ref.shape, arr.shape)
            try:
                np.testing.assert_allclose(arr, ref, rtol=2e-12, atol=2e-12)
            except AssertionError as e:
                raise AssertionError(f"case={case_id} operation={label}\n{e}") from e
    assert next(lines, None) is None, 'Unexpected output beyond last case'


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=('checked', 'optimized', 'both'), default='both')
    parser.add_argument('--cases', type=int, default=400)
    parser.add_argument('--seed', type=int, default=20261103)
    args = parser.parse_args()
    stdin, cases = format_cases(np.random.default_rng(args.seed), args.cases)
    for mode in (('checked', 'optimized') if args.mode == 'both' else (args.mode,)):
        with tempfile.TemporaryDirectory(prefix='scifort-sparse-') as temp:
            build = Path(temp)
            flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
            base = ['gfortran', '-std=f2018', '-ffree-line-length-none', *flags,
                    '-I', str(build), '-J', str(build)]
            objs = []
            for name in ('kinds', 'sparse'):
                obj = build / (name + '.o')
                subprocess.run([*base, '-c', str(ROOT / 'src' / (name + '.f90')),
                                '-o', str(obj)], cwd=build, check=True)
                objs.append(obj)
            exe = build / 'sparse_reference'
            subprocess.run([*base, str(ROOT / 'tools/sparse_reference/sparse_reference_driver.f90'),
                            *(str(o) for o in objs), '-o', str(exe)], cwd=build, check=True)
            result = subprocess.run([str(exe)], input=stdin, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, text=True, cwd=build, check=True)
            parse_and_check(result.stdout, cases)
            print(f'{mode}: {args.cases} randomized sparse matrix cases PASS '
                  f'({len(FIELDS)*args.cases} matrix comparisons)')


if __name__ == '__main__':
    main()
