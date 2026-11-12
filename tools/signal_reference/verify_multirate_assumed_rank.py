#!/usr/bin/env python3
"""Randomized SciPy comparisons for rank-4..15 Fortran assumed-rank multirate APIs.

Usage: python verify_multirate_assumed_rank.py path/to/multirate_assumed_rank_driver
"""
import sys
import subprocess
import numpy as np
from scipy import signal


def check(exe):
    rng = np.random.default_rng(10252026)
    taps_r = np.array([0.2, 0.6, 0.2])
    taps_c = taps_r + 1j * np.array([0.07, -0.02, 0.11])
    cases = 0
    values = 0
    max_abs = 0.0
    # Test first, middle and final axes across a selection of high ranks.
    for rank in (4, 5, 6, 8, 14, 15):
        for a in sorted(set((0, rank // 2, rank - 1))):
            for dtype in ('r', 'c'):
                for op in ('resample', 'poly', 'upfirdn', 'decimate'):
                    for taps_type in ('r', 'c') if op in ('poly', 'upfirdn') else ('r',):
                        dims = [1] * rank
                        dims[0] = 2
                        dims[-1] = 2
                        dims[a] = 9
                        if rank >= 5 and a != rank // 2:
                            dims[rank // 2] = 2
                        mode_list = {
                            'resample': ('hann', 'time', 'freq') if dtype == 'c' else ('hann', 'none'),
                            'poly': ('constant', 'antireflect', 'median', 'line'),
                            'upfirdn': ('constant', 'antisymmetric', 'smooth', 'wrap'),
                            'decimate': ('fir', 'iir'),
                        }[op]
                        for mode in mode_list:
                            shape = tuple(dims)
                            orig = rng.normal(size=shape)
                            if dtype == 'c':
                                orig = orig + 1j * rng.normal(size=shape)
                            p1,p2 = {
                                'resample': (11, 1), 'poly': (3, 2),
                                'upfirdn': (2, 3), 'decimate': (2, 1),
                            }[op]
                            h = taps_r if taps_type == 'r' else taps_c
                            if op == 'resample':
                                kw = {'window':'hann'} if mode == 'hann' else {}
                                if dtype == 'c' and mode == 'freq':
                                    kw = {'domain': 'freq'}
                                target = signal.resample(orig,p1,axis=a,**kw)
                            elif op == 'poly':
                                target = signal.resample_poly(orig,p1,p2,window=h,padtype=mode,axis=a)
                            elif op == 'upfirdn':
                                target = signal.upfirdn(h,orig,up=p1,down=p2,mode=mode,axis=a)
                            else:
                                target = signal.decimate(orig,p1,n=8,ftype=mode,zero_phase=False,axis=a)
                            dims_out = list(target.shape)
                            axis_fortran = (a + 1) if (cases % 2 == 0) else (a-rank)
                            real_flat = np.real(orig).ravel(order='F')
                            imag_flat = np.imag(orig).ravel(order='F')
                            inp = '\n'.join((
                                f'{op} {dtype} {taps_type} {mode} {axis_fortran} {p1} {p2} {rank}',
                                ' '.join(map(str,dims)),
                                ' '.join(map(str,dims_out)),
                                ' '.join(f'{v:.17g}' for v in real_flat),
                                ' '.join(f'{v:.17g}' for v in imag_flat),
                                '',
                            ))
                            proc = subprocess.run([exe],input=inp,text=True,capture_output=True,timeout=20)
                            if proc.returncode:
                                raise AssertionError(f'{rank=} {op=} {dtype=} {taps_type=} {mode=} {a=}\n{proc.stderr}')
                            out = proc.stdout.splitlines()
                            if out[0] != 'STATUS 0':
                                raise AssertionError(f'bad status {out[0]!r}: {rank=} {op=} {mode=} {a=}')
                            got = np.array([complex(*[float(x) for x in line.split()]) for line in out[1:]], dtype=np.complex128)
                            if got.size != target.size:
                                raise AssertionError(f'wrong size {got.size} vs {target.size}: {rank=} {op=} {a=}')
                            exp = target.ravel(order='F')
                            err = float(np.max(np.abs(got-exp)))
                            if not np.isfinite(err) or err > 1e-10:
                                raise AssertionError(f'error {err}: {rank=} {op=} {dtype=} {taps_type=} {mode=} {a=}')
                            max_abs = max(max_abs,err)
                            cases += 1
                            values += len(got)
    print(f'PASS: {cases} SciPy cases, {values} values; max_abs_error={max_abs:.3e}')


if __name__ == '__main__':
    check(sys.argv[1])
