#!/usr/bin/env python3
"""Randomized SciPy checks for SciFort's sparse diagonal and random constructors.

Matches numeric diagonals and the distributional/structural properties of random
matrices. Exact random samples intentionally differ because SciFort uses MRG32k3a,
not NumPy's bit generator.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path
import random
import subprocess
import tempfile

import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[2]


def compile_driver(mode: str, build: Path) -> Path:
    flags = ['-std=f2018', '-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-std=f2018', '-O3']
    cc = ['gfortran', *flags, '-J', str(build), '-I', str(build)]
    objects: list[Path] = []
    for name, source in [('kinds', 'src/kinds.f90'), ('random', 'src/random.f90'),
                         ('sparse', 'src/sparse.f90'), ('constructors', 'src/sparse_constructors.f90')]:
        obj = build / f'{name}.o'
        subprocess.run([*cc, '-c', str(ROOT / source), '-o', str(obj)], cwd=build, check=True)
        objects.append(obj)
    binary = build / 'sparse_constructors_reference_driver'
    subprocess.run([*cc, str(ROOT / 'tools/sparse_reference/sparse_constructors_reference_driver.f90'),
                    *(str(obj) for obj in objects), '-o', str(binary)], cwd=build, check=True)
    return binary


def run_cases(binary: Path, test_cases: list[dict]) -> list[tuple[str, np.ndarray | None, tuple[int, int]]]:
    lines: list[str] = []
    for case in test_cases:
        if case['type'] == 'D':
            lines.append('D')
            rows, cols = case['shape'] if case['shape'] is not None else (-1, -1)
            lines.append(f"{rows} {cols} {len(case['offsets'])}")
            for offset, values in zip(case['offsets'], case['diags']):
                lines.append(f'{offset} {len(values)}')
                if len(values):
                    lines.append(' '.join(f'{v:.17e}' for v in values))
        else:
            lines.append('R')
            lines.append(f"{case['shape'][0]} {case['shape'][1]} {case['density']:.17e} {case['seed']} {case['dist']}")
    output = subprocess.run([str(binary)], input='\n'.join(lines) + '\n', text=True,
                            check=True, capture_output=True).stdout.splitlines()
    pos = 0
    results = []
    for _case in test_cases:
        fields = output[pos].split()
        pos += 1
        if fields[0] == 'ERR':
            results.append(('ERR', None, (0, 0)))
            continue
        assert fields[0] == 'OK', fields
        m, n, nnz = map(int, fields[1:])
        coords = np.empty((nnz, 3))
        for i in range(nnz):
            coords[i] = np.fromstring(output[pos], sep=' ')
            pos += 1
        results.append(('OK', coords, (m, n)))
    assert pos == len(output), f'unparsed Fortran output: {len(output)-pos} lines'
    return results


def generate_cases(rng: random.Random, count: int) -> list[dict]:
    tests: list[dict] = []
    for iteration in range(count):
        shape = (rng.randrange(0, 16), rng.randrange(0, 16))
        use_infer = rng.random() < 0.25
        offsets: list[int] = []
        if use_infer:
            offset0 = rng.randrange(-5, 6)
            width0 = rng.randrange(1, 9)
            shape = (width0 + abs(offset0), width0 + abs(offset0))
            offsets.append(offset0)
        term_count = rng.randint(1, 5)
        all_offsets = list(range(-max(0, shape[0]), max(0, shape[1]) + 1))
        if not all_offsets:
            all_offsets = [0]
        rng.shuffle(all_offsets)
        if use_infer:
            offsets.extend(x for x in all_offsets if x != offsets[0])
        else:
            offsets.extend(all_offsets)
        offsets = offsets[:term_count]
        diags = []
        for index, offset in enumerate(offsets):
            first = max(0, -offset)
            last = min(shape[0], shape[1] - offset)
            needed = max(0, last-first)
            if use_infer and index == 0:
                length = width0
            elif needed == 0:
                length = rng.choice([0,1,2,3])
            elif not use_infer and rng.random() < 0.20:
                length = 1  # SciPy scalar broadcasting with explicit shape.
            else:
                length = needed + rng.choice([0,0,0,1,2])
            values = [float(rng.randrange(-40, 41)) / 10 for _ in range(length)]
            diags.append(values)
        tests.append({'type': 'D', 'shape': None if use_infer else shape,
                      'offsets': offsets, 'diags': diags})

        if iteration % 2 == 0:
            rr, cc = rng.choice([0,1,2,3,4,5,10,15,30,60]), rng.choice([0,1,3,4,7,20,70])
            density = rng.choice([0,0.01,0.1,0.2,0.5,0.8,1.0, rng.random()])
            tests.append({'type': 'R', 'shape': (rr,cc), 'density': density,
                          'seed': rng.randrange(1,100000000), 'dist': rng.choice(['U','N'])})
    return tests


def compare(test_cases: list[dict], results: list[tuple[str,np.ndarray | None,tuple[int,int]]]) -> tuple[int,int]:
    diag_count = rand_count = 0
    for i, (case, (status, coords, shape)) in enumerate(zip(test_cases, results)):
        if case['type'] == 'D':
            if case['shape'] is None:
                try:
                    ref = sparse.diags_array(case['diags'], offsets=case['offsets'],
                                             format='csr').toarray()
                except ValueError:
                    assert status == 'ERR', f'case {i} expected error, got {shape}'
                    continue
            else:
                try:
                    ref = sparse.diags_array(case['diags'], offsets=case['offsets'],
                                             shape=case['shape'], format='csr').toarray()
                except ValueError:
                    assert status == 'ERR', f'case {i} expected error, got {shape}'
                    continue
            assert status == 'OK', f'case {i} expected valid diags, got ERR, case={case}'
            assert tuple(ref.shape) == shape, f'case {i}: {ref.shape} vs {shape}'
            ours = sparse.coo_matrix((coords[:,2],(coords[:,0].astype(int)-1,
                                                  coords[:,1].astype(int)-1)),
                                     shape=shape).toarray()
            np.testing.assert_allclose(ours, ref, rtol=0, atol=1e-13,
                                       err_msg=f'diags case {i} {case}')
            diag_count += 1
        else:
            assert status == 'OK', f'case {i} invalid random'
            m, n = shape
            nnz = len(coords)
            expected = round(case['density'] * m * n)
            assert nnz == expected, f'case {i}: observed {nnz} != expected {expected}'
            assert len(set(map(tuple, coords[:,:2]))) == nnz, f'case {i} duplicate coords'
            assert np.all((coords[:,0] >= 1) & (coords[:,0] <= m)) if nnz else True
            assert np.all((coords[:,1] >= 1) & (coords[:,1] <= n)) if nnz else True
            assert np.isfinite(coords[:,2]).all()
            if case['dist'] == 'U':
                assert (coords[:,2] > 0).all() and (coords[:,2] < 1).all()
            rand_count += 1
    return diag_count,rand_count


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=['checked','optimized','both'], default='both')
    parser.add_argument('--cases', type=int, default=180)
    parser.add_argument('--seed', type=int, default=20261105)
    args = parser.parse_args()
    cases = generate_cases(random.Random(args.seed), args.cases)
    for mode in (['checked', 'optimized'] if args.mode == 'both' else [args.mode]):
        with tempfile.TemporaryDirectory(prefix='scifort-sparse-constructors-') as d:
            binary = compile_driver(mode, Path(d))
            results = run_cases(binary, cases)
            counts = compare(cases,results)
            print(f'{mode}: {counts[0]} SciPy diagonal cases and {counts[1]} random CSR cases passed',flush=True)


if __name__ == '__main__':
    main()
