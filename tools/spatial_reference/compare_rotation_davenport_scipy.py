#!/usr/bin/env python3
"""Independent SciPy comparisons for scalar Davenport-angle conversions.

Checks constructor quaternions, extracted angles, singularity convention,
large and small input axes, and rotation round trips.
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


def compile_driver(tmp: Path, mode: str) -> Path:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    fc = ['gfortran', '-std=f2018', '-ffree-line-length-none', *flags]
    objs = []
    for name in ('kinds', 'spatial_transform'):
        obj = tmp / f'{name}.o'
        subprocess.run([*fc, '-I', str(tmp), '-J', str(tmp), '-c',
                        str(ROOT / 'src' / f'{name}.f90'), '-o', str(obj)],
                       cwd=tmp, check=True, capture_output=True)
        objs.append(str(obj))
    exe = tmp / 'rotation_davenport_driver'
    subprocess.run([*fc, '-I', str(tmp),
                    str(ROOT / 'tools/spatial_reference/rotation_davenport_driver.f90'),
                    *objs, '-o', str(exe)], cwd=tmp, check=True, capture_output=True)
    return exe


def qerror(a, b):
    return float(min(np.max(np.abs(a - b)), np.max(np.abs(a + b))))


def cases_for(seed: int, ncases: int):
    rng = np.random.default_rng(seed)
    cases = []
    for i in range(ncases):
        n = (i % 13 == 0) + 2 if i % 13 == 0 else (1 if i % 17 == 0 else (2 if i % 19 == 0 else 3))
        order = 'e' if i % 2 == 0 else 'i'
        rot = Rotation.random(random_state=rng)
        a, b, _ = rot.as_matrix()
        lam = rng.uniform(-np.pi, np.pi)
        if i % 29 == 0:
            lam = 0.0
        elif i % 29 == 1:
            lam = np.pi/2
        elif i % 29 == 2:
            lam = -np.pi/2
        elif i % 29 == 3:
            lam = np.pi
        elif i % 29 == 4:
            lam = -np.pi
        c = np.cos(lam)*a + np.sin(lam)*np.cross(b, a)
        axes = np.stack((a,b,c))
        if i % 11 == 0:
            scales = 10.0**rng.uniform(-120, 120, size=3)
            axes = axes*scales[:, None]
        elif i % 11 == 1:
            axes = axes*np.array([-2.0,3.0,-4.0])[:, None]
        angles = rng.uniform(-3*np.pi,3*np.pi,size=3)
        if i % 23 == 0:
            angles[:] = 0.0
        elif i % 23 == 1:
            angles[1] = lam
        elif i % 23 == 2:
            angles[1] = lam + np.pi
        elif i % 23 == 3:
            angles[1] = lam + 1e-11
        elif i % 23 == 4:
            angles[1] = lam + np.pi - 1e-11
        # Axis arrays with subnormal/huge magnitudes are fine when normalized.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            ref = Rotation.from_davenport(axes[:n],order,angles[:n]); q = ref.as_quat()
        if i % 3 == 0:
            q = -q
        cases.append((n,order,axes,angles,q))
    return cases


def compare(exe: Path, seed: int, ncases: int):
    cases = cases_for(seed,ncases)
    rows = [str(len(cases))]
    for n, order, axes, angles, q in cases:
        rows.append(f'{n} {order} ' + ' '.join(f'{float(x):.17e}' for x in (*axes.flat,*angles,*q)))
    run = subprocess.run([str(exe)], input='\n'.join(rows)+'\n',
                         text=True, check=True, capture_output=True)
    output = run.stdout.strip().splitlines()
    if len(output) != 2*len(cases):
        raise AssertionError((len(output),len(cases), run.stderr[-300:]))
    errors = {'construct':0.0,'roundtrip':0.0,'angle_regular':0.0}
    locked = 0
    ambiguous = 0
    for i, (n,order,axes,angles,q) in enumerate(cases):
        code = list(map(int,output[2*i].split()))
        nums = np.array([float(s) for s in output[2*i+1].split()])
        if nums.shape != (11,):
            raise AssertionError((i,nums))
        qout, aout, qrecon = nums[:4],nums[4:7],nums[7:]
        expected = Rotation.from_davenport(axes[:n],order,angles[:n]).as_quat()
        errors['construct'] = max(errors['construct'],qerror(qout,expected))
        if code[0] != 0 or qerror(qout,expected) > 5e-13:
            raise AssertionError((i,'from',code, qout,expected))
        if n != 3:
            continue
        if code[1:4] != [0,0,0]:
            raise AssertionError((i,'invalid extraction',code))
        errors['roundtrip'] = max(errors['roundtrip'], qerror(qrecon,q))
        if qerror(qrecon,q) > 5e-9:
            raise AssertionError((i,'roundtrip',axes,angles,aout, qrecon,q))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            scipy_angles = Rotation.from_quat(q).as_davenport(axes,order)
        diff = (aout-scipy_angles+np.pi)%(2*np.pi)-np.pi
        normalized = axes/np.max(abs(axes),axis=1)[:,None]
        normalized /= np.linalg.norm(normalized,axis=1)[:,None]
        degenerate_axis = abs(abs(normalized[0]@normalized[2])-1.0) < 2e-13
        if not degenerate_axis:
            errors['angle_regular'] = max(errors['angle_regular'],float(np.max(abs(diff))))
        if max(abs(diff)) > 2e-7:
            if not degenerate_axis:
                raise AssertionError((i,'as',order,axes,angles,aout,scipy_angles,diff,code))
            ambiguous += 1
        if code[4]:
            locked += 1
            if abs(aout[2]) > 1e-10:
                raise AssertionError((i,'nonzero singular last angle',aout))
    return errors, locked, ambiguous


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--mode', choices=('checked','optimized','both'),default='both')
    parser.add_argument('--seed',type=int,default=20261103)
    parser.add_argument('--cases',type=int,default=1200)
    args=parser.parse_args()
    modes=('checked','optimized') if args.mode=='both' else (args.mode,)
    for mode in modes:
        with tempfile.TemporaryDirectory(prefix='scifort-davenport-') as tmp:
            exe=compile_driver(Path(tmp),mode)
            errors, locks, ambiguous = compare(exe,args.seed,args.cases)
            print(f'{mode}: {args.cases} Davenport cases PASS, {locks} gimbal locks, '
                  f'{ambiguous} equivalent angle branch choices, errors={errors}',flush=True)

if __name__=='__main__':
    main()
