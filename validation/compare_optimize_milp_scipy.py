#!/usr/bin/env python3
"""Compare SciFort's independent branch-and-bound MILP to SciPy HiGHS.

Compiles only scifort_kinds, linprog and milp; does not require FPM.
Use --optimize for -O3 instead of a fully checked build.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp as scipy_milp


def cases(seed: int = 20261110, count: int = 240):
    rng = np.random.default_rng(seed)
    items = []
    for _ in range(count):
        n = int(rng.integers(1, 6))
        m = int(rng.integers(0, 6))
        types = rng.integers(0, 4, size=n, dtype=np.int32)
        lower = np.empty(n)
        upper = np.empty(n)
        feasible_seed = np.empty(n)
        for j, typ in enumerate(types):
            if typ >= 2:
                lower[j] = float(rng.integers(1, 5)) + (0.2 if typ == 3 else 0.0)
                upper[j] = lower[j] + float(rng.integers(2, 5))
                if typ == 3:
                    feasible_seed[j] = 0.0 if rng.random() < 0.4 else float(np.ceil(lower[j]))
                else:
                    feasible_seed[j] = 0.0 if rng.random() < 0.4 else lower[j] + 0.5
            elif typ == 1:
                lower[j] = float(rng.integers(-2, 2))
                upper[j] = lower[j] + float(rng.integers(1, 5))
                feasible_seed[j] = float(rng.integers(int(lower[j]), int(upper[j])+1))
            else:
                lower[j] = float(rng.integers(-2, 2))
                upper[j] = lower[j] + float(rng.integers(1, 5))
                feasible_seed[j] = lower[j] + (upper[j]-lower[j]) * rng.random()
        c = rng.integers(-5, 6, size=n).astype(float) + rng.normal(0.0, 0.1, size=n)
        a = rng.integers(-3, 4, size=(m, n)).astype(float)
        v = a @ feasible_seed
        bl = np.full(m, -np.inf)
        bu = np.full(m, np.inf)
        for row in range(m):
            kind = int(rng.integers(0, 5))
            if kind == 0:
                bl[row] = v[row]
                bu[row] = v[row]
            elif kind == 1:
                bl[row] = v[row] - float(rng.integers(0, 4))
            elif kind == 2:
                bu[row] = v[row] + float(rng.integers(0, 4))
            elif kind == 3:
                bl[row] = v[row] - float(rng.integers(0, 4))
                bu[row] = v[row] + float(rng.integers(0, 4))
        items.append((c, types, lower, upper, a, bl, bu))
    return items


def encode(items):
    lines = [str(len(items))]
    for c, types, lower, upper, a, bl, bu in items:
        n = len(c)
        lines.append(f"{n} {len(bl)}")
        vectors = [c, types, lower, upper]
        vectors.extend(row for row in a)
        if len(bl):
            vectors.extend([bl, bu])
        for v in vectors:
            lines.append(" ".join(str(int(x)) if isinstance(x, np.integer) else
                                  f"{x:.17g}" if np.isfinite(x) else
                                  "1e100" if x > 0 else "-1e100" for x in v))
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--optimize', action='store_true', help='compile with -O3')
    parser.add_argument('--count', type=int, default=240, help='random problems')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    items = cases(count=args.count)
    with tempfile.TemporaryDirectory(prefix='scifort-milp-') as tmp:
        tmp = Path(tmp)
        sources = [root / 'src/kinds.f90', root / 'src/optimize/linprog.f90',
                   root / 'src/optimize/milp.f90', root / 'validation/optimize_milp_driver.f90']
        flags = (['-O3'] if args.optimize else
                 ['-O0', '-g', '-fcheck=all', '-ffpe-trap=invalid,zero,overflow'])
        command = ['gfortran', '-std=f2018', '-Wall', '-Wextra', '-Werror',
                   '-ffree-line-length-none', *flags, *map(str, sources), '-o', str(tmp / 'driver')]
        subprocess.run(command, cwd=tmp, check=True, capture_output=True, text=True)
        process = subprocess.run([str(tmp / 'driver')], input=encode(items),
                                 cwd=tmp, check=True, capture_output=True, text=True)
    lines = process.stdout.strip().splitlines()
    assert len(lines) == 2*len(items), f'Expected 2 lines per case; got {len(lines)}'
    checks = 0
    mismatch = []
    for index, item in enumerate(items):
        c, types, lower, upper, a, bl, bu = item
        row = lines[2*index].split()
        xline = lines[2*index+1]
        assert int(row[0]) == index+1
        status, fun = int(row[1]), float(row[2])
        x = np.fromstring(xline, sep=' ') if xline != 'NONE' else None
        ref = scipy_milp(c, integrality=types, bounds=Bounds(lower, upper),
                         constraints=[LinearConstraint(a, bl, bu)],
                         options={'mip_rel_gap': 0.0, 'node_limit': 100000})
        if x is not None:
            assert len(x)==len(c), (index, x)
            assert abs(float(c@x)-fun) < 2e-6 * max(1.0,abs(fun)), (index,x,fun)
            for j, t in enumerate(types):
                if t >= 2 and abs(x[j]) <= 1e-6:
                    continue
                assert lower[j] - 2e-6 <= x[j] <= upper[j]+2e-6, (index,j,'bounds',x)
                if t in (1,3):
                    assert abs(x[j]-np.rint(x[j])) < 2e-6, (index,j,'integrality',x)
            if len(bl):
                ax = a @ x
                assert np.all(ax >= bl-2e-6) and np.all(ax <= bu+2e-6), (index,'constraints',x,ax,bl,bu)
            checks += 1
        if ref.status == 0:
            if status != 0 or x is None or abs(float(ref.fun)-fun) > 3e-5*max(1.0,abs(float(ref.fun))):
                mismatch.append((index, status, fun, ref.status, ref.fun, str(ref.message), item))
            else:
                checks += 1
        elif ref.status == 2:
            if status != 2:
                mismatch.append((index, status, fun, ref.status, ref.fun, str(ref.message), item))
            else:
                checks += 1
        else:
            # HiGHS failure/inconclusive comparisons are not counted as successes.
            mismatch.append((index, status, fun, ref.status, ref.fun, str(ref.message), item))
    if mismatch:
        for index,status,fun,refstatus,reffun,message,item in mismatch[:8]:
            print(f'Case {index}: Fortran status={status}, fun={fun}, SciPy status={refstatus}, fun={reffun}, message={message}')
            print('c, integrality, lower, upper, a, bl, bu:', item)
        raise AssertionError(f'{len(mismatch)} of {len(items)} SciPy reference cases disagreed')
    print(f'MILP SciPy comparisons passed: {len(items)} problems, {checks} checks ' +
          ('(-O3)' if args.optimize else '(runtime-checked)'))


if __name__ == '__main__':
    main()
