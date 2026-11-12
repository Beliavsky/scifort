#!/usr/bin/env python3
"""Compare SciFort multidimensional resampling driver output with SciPy 1.17+.

Usage: python check_multirate_nd.py /path/to/multirate_nd_driver
"""
import re
import subprocess
import sys

import numpy as np
from scipy import signal


def arrays():
    i2, j2 = np.indices((12, 4)) + 1
    r2 = np.sin(0.17*i2*j2) + 0.03*(i2+j2)
    c2 = r2 + 1j*np.cos(0.2*(i2+2*j2))
    i3, j3, k3 = np.indices((5, 12, 3)) + 1
    r3 = np.sin((i3+j3+k3)*0.12)
    c3 = r3 + 1j*np.cos((i3+2*j3+3*k3)*0.13)
    return r2, c2, r3, c3


def reference(case):
    r2, c2, r3, c3 = arrays()
    taps = np.array([0.2, 0.6, 0.2])
    ctaps = taps + 1j*np.array([0.1,-0.2,0.3])
    if case == 'poly_r2_edge':
        return signal.resample_poly(r2, 2, 3, axis=1, window=taps, padtype='edge')
    if case == 'poly_c2_const':
        return signal.resample_poly(c2, 2, 3, axis=1, window=taps, padtype='constant', cval=0.75-0.25j)
    if case == 'dec_r2_iir':
        return signal.decimate(r2, 2, n=4, ftype='iir', axis=0, zero_phase=False)
    m = re.fullmatch(r'(fourier|poly|up|dec)_(r|c|rr|rc|cr|cc)([23])_a([123])', case)
    if not m:
        raise ValueError('unknown case '+case)
    operation, kind, rank, axis = m.groups()
    rank, axis = int(rank), int(axis)-1
    x = ((r2 if 'r' == kind[-1] else c2) if rank == 2
         else (r3 if 'r' == kind[-1] else c3))
    if operation == 'fourier':
        return signal.resample(x,9 if rank==2 and kind=='r' else
                               7 if rank==2 else 7 if kind=='r' else 6,
                               axis=axis, window='hann' if (rank==2 and kind=='c') or
                               (rank==3 and kind=='r') else None)
    if operation == 'poly':
        if rank==2:
            if kind=='r':
                return signal.resample_poly(x,2,3,axis=axis,window=taps)
            return signal.resample_poly(x,3,2,axis=axis,window=taps,padtype='mean')
        if kind=='r':
            return signal.resample_poly(x,2,3,axis=axis)
        return signal.resample_poly(x,3,2,axis=axis,window=taps)
    if operation == 'up':
        h = taps if kind[0]=='r' else ctaps
        mode = {'rr':'edge','rc':'symmetric','cr':'reflect','cc':'constant'}[kind]
        return signal.upfirdn(h,x,up=2,down=3,axis=axis,mode=mode)
    if operation == 'dec':
        return signal.decimate(x,2,n=12,ftype='fir',axis=axis,zero_phase=(kind=='r'))
    raise AssertionError(operation)


def main(path):
    text = subprocess.run([path], capture_output=True, text=True, check=True).stdout.splitlines()
    count, values, largest_abs, largest_scaled = 0, 0, 0.0, 0.0
    i = 0
    while i < len(text):
        parts = text[i].split()
        if len(parts) < 5 or parts[0] != 'CASE':
            raise RuntimeError(f'malformed header at line {i+1}: {text[i]}')
        name, rank = parts[1], int(parts[2])
        shape=tuple(map(int,parts[3:]))
        if rank != len(shape):
            raise RuntimeError(f'rank mismatch {name} {shape}')
        n = int(np.prod(shape))
        lines=text[i+1:i+1+n]
        if len(lines)<n: raise RuntimeError('short result for '+name)
        got=np.array([complex(*map(float,line.split())) for line in lines]).reshape(shape,order='F')
        want=reference(name)
        if got.shape != want.shape:
            raise AssertionError(f'{name}: shape {got.shape} != {want.shape}')
        err=np.max(np.abs(got-want))
        scale=max(1,np.max(np.abs(want)))
        scaled=err/scale
        largest_abs=max(largest_abs,err)
        largest_scaled=max(largest_scaled,scaled)
        if err>3e-8*scale:
            raise AssertionError(f'{name}: max absolute error {err:.5g}, scale {scale:.5g}')
        count += 1
        values += n
        i += n+1
    print(f'PASS {count} SciPy 1.17+ cases, {values} complex values; '
          f'max abs {largest_abs:.4g}, max scaled {largest_scaled:.4g}')


if __name__ == '__main__':
    main(sys.argv[1])
