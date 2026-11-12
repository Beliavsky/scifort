#!/usr/bin/env python3
"""Independent SciPy validation of seeded random rotations and 2D broadcast grids.

Fortran and NumPy deliberately use different PRNGs; we test the documented
Fortran stream against an independent Park-Miller/Shoemake implementation and
compare distributions statistically, NOT seed-for-seed with SciPy's RNG.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.spatial.transform import Rotation
from scipy.stats import kstest, ks_2samp

ROOT = Path(__file__).resolve().parents[2]
MOD = 2147483647
A = 16807


def build(tmp: Path, mode: str) -> Path:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    fc = ['gfortran', '-std=f2018', '-ffree-line-length-none', *flags]
    objects = []
    for name in ('kinds', 'spatial_transform', 'spatial_rotation_stack', 'spatial_rotation_grid'):
        obj = tmp / (name + '.o')
        subprocess.run([*fc, '-I', str(tmp), '-J', str(tmp), '-c',
                        str(ROOT / 'src' / (name + '.f90')), '-o', str(obj)],
                       cwd=tmp, check=True, capture_output=True)
        objects.append(str(obj))
    exe = tmp / 'rotation_grid_driver'
    subprocess.run([*fc, '-I', str(tmp),
                    str(ROOT / 'tools/spatial_reference/rotation_grid_driver.f90'),
                    *objects, '-o', str(exe)], cwd=tmp, check=True, capture_output=True)
    return exe


def expected_random(seed: int, n: int) -> np.ndarray:
    state = seed % (MOD - 1) + 1

    def uniform():
        nonlocal state
        state = A * state % MOD
        first = state - 1
        state = A * state % MOD
        second = state - 1
        return min((first + second / (MOD - 1)) / (MOD - 1), 1 - np.finfo(float).eps)

    quats = np.empty((n, 4))
    for i in range(n):
        u1, u2, u3 = uniform(), uniform(), uniform()
        a = 2 * np.pi * u2
        b = 2 * np.pi * u3
        q = np.array([np.sqrt(1-u1)*np.sin(a), np.sqrt(1-u1)*np.cos(a),
                      np.sqrt(u1)*np.sin(b), np.sqrt(u1)*np.cos(b)])
        quats[i] = q / np.linalg.norm(q)
    return quats


def qerr(actual: np.ndarray, expected: np.ndarray) -> float:
    return float(np.max(np.minimum(np.max(np.abs(actual-expected), axis=-1),
                                   np.max(np.abs(actual+expected), axis=-1))))


def random_cases(seed: int, ncases: int):
    rng = np.random.default_rng(seed)
    cases = []
    for j in range(ncases):
        n, m = rng.integers(1, 6, size=2)
        k = n if j % 3 else 1
        l = m if j % 4 else 1
        vn = n if j % 5 else 1
        vm = m if j % 2 else 1
        a = Rotation.random(shape=(n,m), rng=rng)
        b = Rotation.random(shape=(k,l), rng=rng)
        qa = a.as_quat()
        qb = b.as_quat()
        # Scalar last, sign arbitrary, and scalar multiplication allowed.
        if j % 3 == 0:
            qa[::2] *= -1
        if j % 5 == 0:
            qb[:,::2] *= 1e70
        if j % 11 == 0:
            qa *= 1e-80
        v = rng.normal(size=(vn,vm,3))
        w0 = rng.uniform(0.1,3.0,size=n)
        w1 = rng.uniform(0.1,3.0,size=m)
        cases.append((qa,qb,v,w0,w1))
    return cases


def run_comparison(exe: Path, seed: int, ncases: int, sample_count: int) -> tuple[float, int]:
    cases = random_cases(seed,ncases)
    lines = [f'{sample_count} {seed}', str(ncases)]
    for qa,qb,v,w0,w1 in cases:
        n,m = qa.shape[:2]
        k,l = qb.shape[:2]
        vn,vm = v.shape[:2]
        lines.append(f'{n} {m} {k} {l} {vn} {vm}')
        for arr in (qa,qb,v):
            for j in range(arr.shape[1]):
                for i in range(arr.shape[0]):
                    lines.append(' '.join(f'{x:.17e}' for x in arr[i,j,:]))
        lines.append(' '.join(f'{x:.17e}' for x in w0))
        lines.append(' '.join(f'{x:.17e}' for x in w1))
    proc = subprocess.run([str(exe)], input='\n'.join(lines)+'\n',
                          capture_output=True, text=True, check=True)
    output = iter(proc.stdout.splitlines())
    assert int(next(output)) == 0, 'Random generator failed'
    qr = np.array([[float(v) for v in next(output).split()] for _ in range(sample_count)])
    expected = expected_random(seed,sample_count)
    worst = qerr(qr,expected)
    assert worst < 2e-14, ('reproducible stream',worst)
    rotated = Rotation.from_quat(qr).apply(np.array([0.,0.,1.]))
    for dim in range(3):
        ks = kstest((rotated[:,dim]+1)/2,'uniform')
        assert ks.pvalue > 1e-6, (dim,ks)
    ref_rotated = Rotation.random(sample_count, rng=np.random.default_rng(seed)).apply([0,0,1])
    for dim in range(3):
        ks = ks_2samp(rotated[:,dim], ref_rotated[:,dim])
        assert ks.pvalue > 1e-6, (dim,ks)
    count = 4 * sample_count

    def get_array(shape: tuple[int, ...], ncol: int) -> np.ndarray:
        nonlocal count
        arr = np.empty(shape)
        if len(shape) == 3:
            n,m,c = shape
            assert c == ncol
            for j in range(m):
                for i in range(n):
                    row = np.fromstring(next(output),sep=' ')
                    assert row.size == c, (i,j,row,c)
                    arr[i,j,:] = row
        else:
            n,c = shape
            assert c == ncol
            for i in range(n):
                row = np.fromstring(next(output),sep=' ')
                assert row.size == c
                arr[i,:] = row
        count += arr.size
        return arr

    for case_id,(qa,qb,v,w0,w1) in enumerate(cases):
        assert int(next(output)) == 0, case_id
        n,m = qa.shape[:2]
        k,l = qb.shape[:2]
        vn,vm = v.shape[:2]
        q = get_array((n,m,4),4)
        rv = get_array((n,m,3),3)
        flat_m = get_array((n,m,9),9)
        mats = np.empty((n,m,3,3))
        for j in range(m):
            for i in range(n):
                mats[i,j,:,:] = flat_m[i,j,:].reshape(3,3,order='F')
        mrps = get_array((n,m,3),3)
        eulers = get_array((n,m,3),3)
        back_mrp = get_array((n,m,4),4)
        back_matrix = get_array((n,m,4),4)
        back_euler = get_array((n,m,4),4)
        mags = np.array([np.fromstring(next(output),sep=' ') for _ in range(m)]).T
        assert mags.shape == (n,m)
        count += n*m
        applied = get_array((max(n,vn),max(m,vm),3),3)
        inverse_applied = get_array((max(n,vn),max(m,vm),3),3)
        composed_q = get_array((max(n,k),max(m,l),4),4)
        mean0_q = get_array((m,4),4)
        mean1_q = get_array((n,4),4)
        a = Rotation.from_quat(qa/np.max(np.abs(qa),axis=-1,keepdims=True))
        b = Rotation.from_quat(qb/np.max(np.abs(qb),axis=-1,keepdims=True))
        errors = (
            qerr(q,a.as_quat(canonical=True)),
            float(np.max(np.abs(rv-a.as_rotvec()))),
            float(np.max(np.abs(mats-a.as_matrix()))),
            float(np.max(np.abs(mrps-a.as_mrp()))),
            float(np.max(np.abs(eulers-a.as_euler('ZYX')))),
            qerr(back_mrp,a.as_quat()),
            qerr(back_matrix,a.as_quat()),
            qerr(back_euler,a.as_quat()),
            float(np.max(np.abs(mags-a.magnitude()))),
            float(np.max(np.abs(applied-a.apply(v)))),
            float(np.max(np.abs(inverse_applied-a.apply(v,inverse=True)))),
            qerr(composed_q,(a*b).as_quat()),
            qerr(mean0_q,a.mean(weights=w0[:,None],axis=0).as_quat()),
            qerr(mean1_q,a.mean(weights=w1,axis=1).as_quat()),
        )
        error = max(errors)
        worst = max(error,worst)
        assert error < 2e-10, (case_id, errors)
    assert next(output,None) is None, 'unexpected driver output'
    return worst,count


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--mode',choices=('checked','optimized','both'),default='both')
    parser.add_argument('--count',type=int,default=140)
    parser.add_argument('--random-count',type=int,default=12000)
    args=parser.parse_args()
    modes=['checked','optimized'] if args.mode=='both' else [args.mode]
    for mode in modes:
        with tempfile.TemporaryDirectory(prefix='scifort-rotation-grid-') as td:
            exe=build(Path(td),mode)
            worst=0.0
            checked=0
            for seed in (20261102,20261103):
                err,num=run_comparison(exe,seed,args.count,args.random_count)
                worst=max(worst,err)
                checked+=num
                print(f'{mode}: seed {seed}: {args.count} grids, {args.random_count} random rotations: '
                      f'max error {err:.3e}',flush=True)
            print(f'{mode}: {2*args.count} grids, {2*args.random_count} random rotations, '
                  f'{checked} values: PASS (max error {worst:.3e})',flush=True)

if __name__=='__main__':
    main()
