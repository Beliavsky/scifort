#!/usr/bin/env python3
"""Compare high-dimensional supplied-simplex interpolation against SciPy/Qhull.

The geometry is constructed by scipy.spatial.Delaunay; Fortran only consumes
its simplex vertex indices, so this does NOT test native >=4D triangulation.
Checks both scalar/multichannel real and complex values, 2D--7D, shape and fill.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay

ROOT = Path(__file__).resolve().parents[2]


def run(command: list[str], *, cwd: Path, stdin: str | None = None) -> str:
    proc = subprocess.run(command, input=stdin, text=True, capture_output=True, cwd=cwd)
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(command)}\n{proc.stdout[-1200:]}\n{proc.stderr[-1800:]}")
    return proc.stdout


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed', type=int, default=315024)
    parser.add_argument('--cases', type=int, default=75)
    parser.add_argument('--only', choices=('checked', 'optimized', 'both'), default='both')
    args = parser.parse_args()
    modes = (('checked', ['-O0', '-g', '-fcheck=all', '-fbacktrace']), ('optimized', ['-O3']))
    for name, flags in modes:
        if args.only not in (name, 'both'):
            continue
        rng = np.random.default_rng(args.seed)
        worst = 0.0
        nvalues = 0
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            compiler = ['gfortran', '-std=f2018', '-ffree-line-length-none', *flags]
            modules = ('kinds', 'interpolate_kdtree', 'interpolate_scattered')
            for module in modules:
                run([*compiler, '-I' + str(tmp), '-J' + str(tmp), '-c',
                     str(ROOT / 'src' / f'{module}.f90'),
                     '-o', str(tmp / f'{module}.o')], cwd=tmp)
            executable = tmp / 'simplex_nd_driver'
            run([*compiler, '-I' + str(tmp),
                 str(ROOT / 'tools/interpolate_reference/simplex_nd_driver.f90'),
                 *(str(tmp / f'{m}.o') for m in reversed(modules)),
                 '-o', str(executable)], cwd=tmp)
            for case in range(args.cases):
                dim = (2, 3, 4, 5, 6, 7)[case % 6]
                n = dim + 5 + (case % 9)
                channels = 1 + case % 3
                pts = rng.standard_normal((n, dim))
                if case % 13 == 0:
                    pts += 120.0
                if case % 7 == 0:
                    pts[:, -1] *= 0.1
                tri = Delaunay(pts)
                if tri.simplices.size == 0:
                    continue
                yr = rng.standard_normal((n, channels))
                yi = rng.standard_normal((n, channels))
                # Points inside random mesh simplices plus exact data sites.
                ni = 9
                ids = rng.integers(len(tri.simplices), size=ni)
                weights = rng.dirichlet(np.ones(dim+1), size=ni)
                q = np.einsum('ij,ijk->ik', weights, pts[tri.simplices[ids]])
                q = np.vstack([q, pts[case % n], pts.mean(axis=0),
                               pts.max(axis=0) + 2.0])
                real_ref = np.asarray(LinearNDInterpolator(pts, yr, fill_value=-91.0)(q))
                complex_ref = np.asarray(LinearNDInterpolator(pts, yr + 1j*yi,
                                                             fill_value=-91.0 + 17.0j)(q))
                lines = [f'{n} {dim} {channels} {len(tri.simplices)} {len(q)}']
                lines.extend(' '.join(f'{val:.17g}' for val in row) for row in pts)
                lines.extend(' '.join(f'{val:.17g}' for val in [*a, *b]) for a, b in zip(yr, yi))
                lines.extend(' '.join(str(v+1) for v in row) for row in tri.simplices)
                lines.extend(' '.join(f'{val:.17g}' for val in row) for row in q)
                out = run([str(executable)], cwd=tmp, stdin='\n'.join(lines)+'\n').splitlines()
                if [int(s) for s in out[:4]] != [0, 0, 0, 0]:
                    raise AssertionError(f'{name} case {case} invalid status: {out[:4]}')
                observed = np.array([[float(s) for s in line.split()] for line in out[4:]])
                if observed.shape != (len(q), 3 * channels):
                    raise AssertionError(f'{name} case {case} bad shape {observed.shape}')
                real_got = observed[:, :channels]
                complex_got = observed[:, channels:2*channels] + 1j*observed[:, 2*channels:]
                err = max(float(np.max(np.abs(real_ref - real_got))),
                          float(np.max(np.abs(complex_ref - complex_got))))
                worst = max(worst, err)
                if err > 2e-9:
                    raise AssertionError(f'{name} case {case}, dimension {dim}, error {err}')
                nvalues += 3 * real_got.size
        print(f'{name}: {args.cases} cases (dimensions 2-7), {nvalues} values, '
              f'max error {worst:.3g}: PASS', flush=True)


if __name__ == '__main__':
    main()
