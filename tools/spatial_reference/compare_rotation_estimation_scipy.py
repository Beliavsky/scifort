#!/usr/bin/env python3
"""Independent SciPy reference for vector alignment and chordal rotation means.

The test covers weighted Wahba/Kabsch alignment, sensitivity, root residuals,
shortest single-vector alignment, one infinite primary weight, and weighted
quaternion chordal means. The Fortran executable is compiled from local code.

    python tools/spatial_reference/compare_rotation_estimation_scipy.py --mode both
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile
import warnings

import numpy as np
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parents[2]


def compile_program(temp: Path, mode: str) -> Path:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    command = ['gfortran', '-std=f2018', '-ffree-line-length-none', *flags]
    objects = []
    for name in ('kinds', 'spatial_transform'):
        obj = temp / (name + '.o')
        subprocess.run([*command, '-I', str(temp), '-J', str(temp), '-c',
                        str(ROOT / 'src' / (name + '.f90')), '-o', str(obj)],
                       cwd=temp, check=True, capture_output=True)
        objects.append(str(obj))
    exe = temp / 'rotation_estimation_driver'
    subprocess.run([*command, '-I', str(temp),
                    str(ROOT / 'tools/spatial_reference/rotation_estimation_driver.f90'),
                    *objects, '-o', str(exe)], cwd=temp, check=True, capture_output=True)
    return exe


def make_cases(seed: int, number: int):
    rng = np.random.default_rng(seed)
    cases = []
    for z in range(number):
        n = [1, 2, 3, 5, 10, 20][z % 6]
        m = [1, 2, 3, 10][z % 4]
        actual = Rotation.random(random_state=rng)
        b = rng.normal(size=(n, 3)) * 10**rng.uniform(-3, 3)
        a = actual.apply(b)
        if z % 5 == 0:
            a += rng.normal(size=(n, 3)) * 0.1
        if z % 7 == 0 and n > 1:
            a[1, :] *= 1.2
        if z % 29 == 0 and n == 1:
            a[0, :] = [-1.0, 1e-8, 0.0]
            b[0, :] = [1.0, 0.0, 0.0]
        w = np.exp(rng.uniform(-3, 3, size=n))
        if z % 13 == 0:
            w[0] = np.inf
        q = np.atleast_2d(Rotation.random(m, random_state=rng).as_quat())
        wm = np.exp(rng.uniform(-3, 3, size=m))
        if z % 11 == 0 and m > 1:
            wm[0] = 0.0
        if z % 23 == 0:
            q[0] *= -1.0
        sensitivity = int(n > 1 and np.all(np.isfinite(w)) and z % 3 == 0)
        cases.append((a, b, w, q, wm, sensitivity))
    return cases


def run_comparison(exe: Path, seed: int, count: int):
    cases = make_cases(seed, count)
    lines = [str(len(cases))]
    for a, b, w, q, wm, sensitivity in cases:
        lines.append(f'{len(a)} {len(q)} {sensitivity}')
        for ai, bi, wi in zip(a, b, w):
            lines.append(' '.join(f'{float(v):.17e}' for v in (*ai, *bi, wi)))
        for quat, weight in zip(q, wm):
            lines.append(' '.join(f'{float(v):.17e}' for v in (*quat, weight)))
    output = subprocess.run([str(exe)], input='\n'.join(lines) + '\n',
                            text=True, capture_output=True, check=True)
    rows = output.stdout.splitlines()
    if len(rows) != 2 * len(cases):
        raise AssertionError((len(rows), len(cases)))
    largest = dict(align_angle=0.0, mean_angle=0.0,
                   residual=0.0, sensitivity=0.0)
    checked_sensitivity = 0
    infinite_cases = 0
    for k, (a, b, w, q, wm, want_sens) in enumerate(cases):
        codes = [int(x) for x in rows[2*k].split()]
        if codes != [0, 0]:
            raise AssertionError((k, 'Fortran status', codes))
        result = np.array([float(x) for x in rows[2*k+1].split()])
        if len(result) != 18:
            raise AssertionError((k, 'Fortran result size', len(result)))
        got_align = Rotation.from_quat(result[:4])
        got_rssd = result[4]
        got_mean = Rotation.from_quat(result[5:9])
        got_sens = result[9:].reshape(3, 3)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            if want_sens:
                target_align, target_rssd, target_sens = Rotation.align_vectors(
                    a, b, weights=w, return_sensitivity=True)
                checked_sensitivity += 1
            else:
                target_align, target_rssd = Rotation.align_vectors(a, b, weights=w)
            target_mean = Rotation.from_quat(q).mean(weights=wm)
        align_error = float(np.linalg.norm((got_align.inv() * target_align).as_rotvec()))
        mean_error = float(np.linalg.norm((got_mean.inv() * target_mean).as_rotvec()))
        if align_error > 2e-8 or mean_error > 2e-8:
            raise AssertionError((k, 'rotation', align_error, mean_error))
        # SciPy computes rssd from two large squared norms minus a dot product,
        # which can introduce an O(sqrt(machine_epsilon))*vector-scale error
        # for perfectly aligned observations. Also validate against residuals
        # computed independently with the rotation returned by SciPy.
        finite = np.isfinite(w)
        target_direct = np.sqrt(np.sum(w[finite] * np.sum(
            (a[finite] - target_align.apply(b[finite]))**2, axis=1)))
        residual_error = abs(got_rssd - target_direct)
        scale = max(1.0, float(np.max(abs(a))), float(np.max(abs(b))), target_direct)
        if residual_error > 3e-10 * scale:
            raise AssertionError((k, 'direct residual', got_rssd, target_direct))
        if abs(got_rssd - target_rssd) > 4e-7 * scale:
            raise AssertionError((k, 'SciPy rssd', got_rssd, target_rssd))
        if want_sens:
            sensitivity_error = float(np.linalg.norm(got_sens - target_sens) /
                                      (1.0 + np.linalg.norm(target_sens)))
            if sensitivity_error > 2e-8:
                raise AssertionError((k, 'sensitivity', sensitivity_error))
            largest['sensitivity'] = max(largest['sensitivity'], sensitivity_error)
        if np.any(~finite):
            infinite_cases += 1
        largest['align_angle'] = max(largest['align_angle'], align_error)
        largest['mean_angle'] = max(largest['mean_angle'], mean_error)
        largest['residual'] = max(largest['residual'], residual_error / scale)
    return checked_sensitivity, infinite_cases, largest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=('both','checked','optimized'), default='both')
    parser.add_argument('--seed', type=int, default=84211)
    parser.add_argument('--cases', type=int, default=360)
    args = parser.parse_args()
    modes = ('checked','optimized') if args.mode == 'both' else (args.mode,)
    for mode in modes:
        with tempfile.TemporaryDirectory(prefix='scifort-rotation-estimation-') as td:
            exe = compile_program(Path(td), mode)
            sensitivity, infinite, err = run_comparison(exe, args.seed, args.cases)
            print(f'{mode}: {args.cases} SciPy vector-fit and mean cases PASS; '
                  f'{sensitivity} sensitivity matrices, '
                  f'{infinite} infinite-primary cases; '
                  + ', '.join(f'{key}={value:.3g}' for key,value in err.items()),
                  flush=True)


if __name__ == '__main__':
    main()
