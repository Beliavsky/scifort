#!/usr/bin/env python3
"""Compare SciFort's Euler and modified Rodrigues conversions to SciPy.

Runs all 24 three-angle conventions, random orientations, exact/near gimbal
locks, high-angle rotations and MRPs in both shadow and principal sets.

    python tools/spatial_reference/compare_rotation_euler_mrp_scipy.py --mode both
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
SEQUENCES = (
    'xyz', 'xzy', 'yxz', 'yzx', 'zxy', 'zyx',
    'xyx', 'xzx', 'yxy', 'yzy', 'zxz', 'zyz',
    'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX',
    'XYX', 'XZX', 'YXY', 'YZY', 'ZXZ', 'ZYZ',
)


def compile_driver(tmp: Path, mode: str) -> Path:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    fc = ['gfortran', '-std=f2018', '-ffree-line-length-none', *flags]
    objs = []
    for name in ('kinds', 'spatial_transform'):
        obj = tmp / (name + '.o')
        subprocess.run([*fc, '-I', str(tmp), '-J', str(tmp), '-c',
                        str(ROOT / 'src' / (name + '.f90')), '-o', str(obj)],
                       cwd=tmp, check=True, capture_output=True)
        objs.append(str(obj))
    exe = tmp / 'rotation_euler_mrp_driver'
    subprocess.run([*fc, '-I', str(tmp),
                    str(ROOT / 'tools/spatial_reference/rotation_euler_mrp_driver.f90'),
                    *objs, '-o', str(exe)], cwd=tmp, check=True, capture_output=True)
    return exe


def qerror(a: np.ndarray, b: np.ndarray) -> float:
    return float(min(np.max(abs(a - b)), np.max(abs(a + b))))


def make_cases(seed: int, per_seq: int):
    rng = np.random.default_rng(seed)
    cases = []
    for seq in SEQUENCES:
        repeated = seq[0] == seq[2]
        for i in range(per_seq):
            angles = rng.uniform(-np.pi, np.pi, size=3)
            if i % 13 in (0,1,2,3):
                angles[1] = (0.0 if i % 13 == 0 else np.pi) if repeated else \
                    (np.pi/2 if i % 13 <= 1 else -np.pi/2)
            elif i % 13 == 4:
                angles[1] = (np.pi if repeated else np.pi/2) + 1e-12
            elif i % 13 == 5:
                angles[1] = (0.0 if repeated else -np.pi/2) - 1e-12
            elif i % 13 == 6:
                angles[1] = (0.0 if repeated else np.pi/2) + 2e-6
            elif i % 13 == 7:
                angles *= 0.0
            elif i % 13 == 8:
                angles[0] += 40*np.pi
            ref = Rotation.from_euler(seq, angles)
            quat = ref.as_quat()
            if i % 4 == 0:
                quat = -quat
            if i % 7 == 0:
                quat *= 1e150
            elif i % 7 == 1:
                quat *= 1e-150
            sigma = rng.normal(size=3)*10.0**rng.uniform(-5.0,5.0)
            if i % 17 == 0:
                sigma = np.zeros(3)
            if i % 17 == 1:
                sigma = rng.normal(size=3)*1e100
            if i % 17 == 2:
                sigma = rng.normal(size=3)*1e-100
            if i % 17 == 3:
                sigma = np.array([0.0,0.0,1.0])
            cases.append((seq, angles, quat, sigma))
    return cases


def compare(exe: Path, seed: int, per_seq: int):
    cases = make_cases(seed, per_seq)
    rows = [str(len(cases))]
    for seq, a, q, sigma in cases:
        rows.append(seq + ' ' + ' '.join(f'{float(x):.17e}' for x in (*a,*q,*sigma)))
    output = subprocess.run([str(exe)], input='\n'.join(rows) + '\n',
                            text=True, capture_output=True, check=True)
    lines = output.stdout.splitlines()
    assert len(lines) == 2*len(cases), (len(lines),len(cases))
    errors = dict(euler=0.0, mrp=0.0, euler_angles=0.0, mrp_params=0.0)
    locks = 0
    for idx, (seq, angles, q, sigma) in enumerate(cases):
        codes = [int(x) for x in lines[2*idx].split()]
        assert codes[:5] == [0]*5, (idx,codes)
        nums = np.array([float(x) for x in lines[2*idx+1].split()])
        assert len(nums) == 22, (idx,len(nums))
        q_e = nums[0:4]
        a_out = nums[4:7]
        q_back = nums[7:11]
        q_m = nums[11:15]
        mrp_out = nums[15:18]
        q_out = nums[18:22]
        scipy_e = Rotation.from_euler(seq,angles)
        scipy_base = Rotation.from_quat(q / np.max(abs(q)))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',UserWarning)
            scipy_angles = scipy_base.as_euler(seq)
        scipy_m = Rotation.from_mrp(sigma)
        errors['euler'] = max(errors['euler'],qerror(q_e,scipy_e.as_quat()),
                              qerror(q_back,scipy_base.as_quat()),
                              qerror(q_out,scipy_base.as_quat()))
        e = qerror(q_m,scipy_m.as_quat())
        errors['mrp'] = max(errors['mrp'],e)
        d = np.mod(a_out - scipy_angles+np.pi,2*np.pi)-np.pi
        # At singularities different gauges can give different angles but
        # equivalent rotations. SciFort follows SciPy's third-angle-zero gauge.
        errors['euler_angles'] = max(errors['euler_angles'],float(max(abs(d))))
        errors['mrp_params'] = max(errors['mrp_params'],
                                   float(max(abs(mrp_out-scipy_base.as_mrp()))))
        if codes[5]:
            locks += 1
            assert abs(a_out[-1]) < 1e-14, (seq,idx,a_out)
        if errors['euler'] > 3e-8:
            raise AssertionError((idx,seq,angles,'Euler quat',errors))
        if e > 2e-12:
            raise AssertionError((idx,seq,sigma,'MRP quat',e))
        if max(abs(d)) > 3e-8:
            raise AssertionError((idx,seq,angles,'Euler angle',d,a_out,scipy_angles,codes))
        if max(abs(mrp_out-scipy_base.as_mrp())) > 1e-12:
            raise AssertionError((idx,seq,'MRP roundtrip',mrp_out,scipy_base.as_mrp()))
    return len(cases),locks,errors


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=('both','checked','optimized'), default='both')
    parser.add_argument('--seed', type=int, default=45617)
    parser.add_argument('--per-sequence', type=int, default=80)
    args = parser.parse_args()
    modes = ('checked','optimized') if args.mode == 'both' else (args.mode,)
    for mode in modes:
        with tempfile.TemporaryDirectory(prefix='scifort-euler-mrp-') as td:
            exe = compile_driver(Path(td),mode)
            count, locks, err = compare(exe,args.seed,args.per_sequence)
            print(f'{mode}: {count} Euler/MRP scenarios across 24 sequences, '
                  f'{locks} gimbal locks PASS; '
                  + ', '.join(f'{k}={v:.3g}' for k,v in err.items()),flush=True)


if __name__ == '__main__':
    main()
