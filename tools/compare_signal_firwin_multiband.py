#!/usr/bin/env python3
"""Independent SciPy reference comparisons for scalar/logical/multiband firwin.

Usage: python tools/compare_signal_firwin_multiband.py /path/to/reference_exe
Requires NumPy and SciPy 1.17 or newer.  The reference output is produced
by the independent Fortran FIR code, never derived from Python coefficients.
"""
from __future__ import annotations

import random
import subprocess
import sys
import numpy as np
from scipy import signal


def cases():
    rng = random.Random(271828)
    base = {1:[.29], 2:[.13,.57], 3:[.09,.28,.65],
            4:[.07,.22,.46,.79], 5:[.05,.19,.35,.53,.81]}
    # Deliberately include the full pass_zero alphabet for each edge parity.
    for ncut in range(1,6):
        for pass_id in range(7):
            for k in range(16):
                n = [11,12,17,18,27,28,33,34][k%8]
                win_id = (k+pass_id) % 4
                scale = (k // 2) % 2
                width = (k // 4) % 2
                scalar = int(ncut == 1 and ((k+pass_id) % 2 == 0))
                fs = [2.0, 8.0, 150.0][(k+pass_id)//2%3]
                bounds = np.array(base[ncut], dtype=float)
                # Small independent perturbations keep the knots ordered.
                bounds += np.array([rng.uniform(-.004,.004) for _ in bounds])
                bounds *= fs/2
                yield n,ncut,pass_id,win_id,scale,width,scalar,fs,bounds
    # Short and long tap lengths and other aliasing corner cases.
    for ncut in range(1,6):
        for pass_id in range(7):
            for n in (1,2,3,4,49,50,61):
                fs=2.0
                yield n,ncut,pass_id,n%4,1,n%2, int(ncut==1 and n%2==0),fs,                    np.array(base[ncut])*fs/2


def compare(executable):
    samples = list(cases())
    payload = []
    for n,ncut,p,w,scale,width,scalar,fs,bounds in samples:
        limits = list(bounds) + [0.0]*(5-ncut)
        payload.append(' '.join(map(str,[n,ncut,p,w,scale,width,scalar,fs,*limits])))
    result = subprocess.run([executable],input='\n'.join(payload)+'\n',
                            text=True, capture_output=True, check=True)
    lines = result.stdout.splitlines()
    assert len(lines) == len(samples), (len(lines),len(samples),result.stderr)
    failed = 0
    passed = 0
    coefficients = 0
    max_error = 0.0
    names={2:'lowpass',3:'highpass',4:'bandpass',5:'bandstop'}
    taper_names=['hamming','hann','blackman',('kaiser',3.2)]
    for idx,(line,case) in enumerate(zip(lines,samples)):
        n,ncut,p,w,scale,width,scalar,fs,bounds = case
        args = line.split()
        assert len(args) >= 2, (idx,line)
        status, count = map(int,args[:2])
        output = np.array(list(map(float,args[2:])))
        assert len(output)==count, (idx,count,len(output))
        kwargs={'fs':fs, 'scale':bool(scale),
                'window':taper_names[w]}
        if width:kwargs['width']=0.04*fs
        if p in (0,1):kwargs['pass_zero']=bool(p)
        elif p in names:kwargs['pass_zero']=names[p]
        cutoff = float(bounds[0]) if scalar else bounds
        try:
            truth = signal.firwin(n,cutoff,**kwargs)
        except ValueError:
            failed += 1
            if status==0 or count!=0:
                raise AssertionError(f'case {idx}: expected SciPy rejection, got status={status}, len={count}, {case}')
            continue
        if status!=0 or count!=len(truth):
            raise AssertionError(f'case {idx}: valid design rejected status={status}, len={count}, {case}')
        error = float(np.max(np.abs(output-truth)))
        if not np.allclose(output,truth,rtol=5e-12,atol=1e-12):
            raise AssertionError(f'case {idx}: error {error:g}, {case}, Fortran={output}, SciPy={truth}')
        max_error = max(max_error,error)
        passed += 1
        coefficients += len(truth)
    print(f'PASS: {passed} valid FIR designs, {failed} invalid designs, '
          f'{coefficients} coefficients vs SciPy {signal.__version__ if hasattr(signal,"__version__") else __import__("scipy").__version__}; '
          f'max error {max_error:.3e}')


if __name__=='__main__':
    compare(sys.argv[1])
