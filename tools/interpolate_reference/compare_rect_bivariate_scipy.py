#!/usr/bin/env python3
"""Compare SciFort rectangular-grid tensor splines to SciPy FITPACK.

Builds locally with gfortran; SciPy is a numeric oracle only, never linked.
"""
from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy.interpolate import RectBivariateSpline

ROOT = Path(__file__).resolve().parents[2]
FMT = lambda v: ' '.join(f'{float(a):.17g}' for a in np.ravel(v))


def compile_driver(work: Path, flags: list[str]) -> Path:
    exe = work / 'rect_bivariate_driver'
    sources = ['kinds.f90', 'interpolate_bspline.f90',
               'interpolate_fitpack_bivariate.f90']
    cmd = ['gfortran', '-std=f2018', '-ffree-line-length-none', *flags,
           '-J' + str(work), '-I' + str(work),
           *(str(ROOT / 'src' / src) for src in sources),
           str(ROOT / 'tools/interpolate_reference/rect_bivariate_driver.f90'),
           '-o', str(exe)]
    subprocess.run(cmd, check=True, text=True, capture_output=True)
    return exe


def compare_case(exe: Path, rng: np.random.Generator, case: int) -> np.ndarray:
    kx = 1 + case % 5
    ky = 1 + (case // 5 + 2 * case) % 5
    nx = kx + 1 + (case * 11) % 21
    ny = ky + 1 + (case * 7) % 18
    x = np.cumsum(rng.uniform(0.1, 1.1, nx))
    y = np.cumsum(rng.uniform(0.15, 1.3, ny))
    x = (x - x[0]) / (x[-1] - x[0]) * 2.7 - 1.0
    y = (y - y[0]) / (y[-1] - y[0]) * 3.1 - 1.4
    # Smooth and noisy grid data; no separability assumed in comparisons.
    z = (np.sin(1.7 * x[:, None] - 0.8 * y[None, :]) +
         0.12 * x[:, None]**2 * y[None, :] +
         0.05 * rng.standard_normal((nx, ny)))
    nq = 25
    qx = rng.uniform(x[0] - 0.2, x[-1] + 0.2, nq)
    qy = rng.uniform(y[0] - 0.3, y[-1] + 0.3, nq)
    bounds = np.empty((8, 4))
    for j in range(8):
        bounds[j] = [*rng.uniform(x[0] - 0.5, x[-1] + 0.5, 2),
                     *rng.uniform(y[0] - 0.5, y[-1] + 0.5, 2)]
    bounds[0] = [x[0] - 2, x[-1] + 2, y[0] - 2, y[-1] + 2]
    bounds[1] = [x[-1]+1, x[-1]+2, y[0], y[-1]]
    lines = [f'{nx} {ny} {kx} {ky} {nq} {len(bounds)}',
             FMT(x), FMT(y), FMT(z.ravel(order='F')), FMT(qx), FMT(qy)]
    lines.extend(FMT(r) for r in bounds)
    result = subprocess.run([str(exe)], input='\n'.join(lines)+'\n',
                            capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(f'Fortran exit {result.returncode}: {result.stderr[-1200:]}')
    output = result.stdout.splitlines()
    if int(output[0]):
        raise AssertionError(f'case {case}: invalid Fortran status {output[0]}')
    nt_x,nt_y,nc = map(int,output[1].split())
    tx = np.fromstring(output[2],sep=' ')
    ty = np.fromstring(output[3],sep=' ')
    coeff = np.fromstring(output[4],sep=' ')
    assert (len(tx),len(ty),len(coeff)) == (nt_x,nt_y,nc)
    vals = np.array([np.fromstring(line,sep=' ') for line in output[5:5+nq]])
    integral = np.array([float(line) for line in output[5+nq:5+nq+len(bounds)]])
    spline = RectBivariateSpline(x,y,z,kx=kx,ky=ky,s=0)
    tx_ref,ty_ref = spline.get_knots()
    cr = spline.get_coeffs()
    refval = np.zeros_like(vals)
    for j,(dx,dy) in enumerate([(0,0),(1,0),(0,1),(2,0),(0,2),(1,1)]):
        if dx < kx and dy < ky:
            # SciPy/FITPACK's derivative evaluator only accepts nu < degree.
            refval[:,j] = spline.ev(np.clip(qx,x[0],x[-1]),
                                    np.clip(qy,y[0],y[-1]),dx=dx,dy=dy)
    refint = np.array([spline.integral(*b) for b in bounds])
    # Knot vectors and coefficient arrays should match FITPACK for exact
    # rectangular interpolation; use scale-aware errors for all outcomes.
    diff = np.array([
        np.max(np.abs(tx-tx_ref)),
        np.max(np.abs(ty-ty_ref)),
        np.max(np.abs(coeff-cr)/np.maximum(1.,np.abs(cr))),
        np.max(np.abs(vals-refval)/np.maximum(1.,np.abs(refval))),
        np.max(np.abs(integral-refint)/np.maximum(1.,np.abs(refint)))
    ])
    # Dense not-a-knot collocation may differ from FITPACK in the final ulps.
    limits = [3e-13,3e-13,2e-7,2e-6,2e-7]
    if np.any(diff > limits):
        raise AssertionError(f'case {case}: nx,ny={nx},{ny} kx,ky={kx},{ky} errors {diff}')
    # The interpolation property is separately tested, independent of SciPy.
    return diff


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--cases',type=int,default=90)
    parser.add_argument('--seed',type=int,default=20261101)
    parser.add_argument('--mode',choices=['both','checked','optimized'],default='both')
    args=parser.parse_args()
    modes=[]
    if args.mode!='optimized':
        modes.append(('checked',['-O0','-g','-fcheck=all','-fbacktrace']))
    if args.mode!='checked':
        modes.append(('optimized',['-O3']))
    for name,flags in modes:
        with tempfile.TemporaryDirectory() as dirname:
            exe=compile_driver(Path(dirname),flags)
            rng=np.random.default_rng(args.seed)
            errors=np.zeros(5)
            for i in range(args.cases):
                errors=np.maximum(errors,compare_case(exe,rng,i))
            print(f'{name}: {args.cases}/{args.cases} scipy rectangular B-spline cases PASS; '
                  f'max errors knot_x/knot_y/coeff/value/integral: {errors}',flush=True)


if __name__ == '__main__':
    main()
