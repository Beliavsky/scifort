#!/usr/bin/env python3
"""Compare Fortran 1D upfirdn/resample_poly/FIR decimate with SciPy.

Compile resample_driver.f90 with kinds, signal_spectral, signal_design,
signal_lfilter, signal_sos, signal_filtfilt, and signal_resample
then run this script with its path as the first command-line argument.
"""
from __future__ import annotations
import os
import subprocess
import sys
from pathlib import Path
import numpy as np
from scipy import signal

rng=np.random.default_rng(4819)
mode_list=['constant','edge','reflect','symmetric','wrap','line','smooth']
cases=[]
for typ in (1,2,3,4,5,6,7):
    for j in range(160):
        n=int(rng.integers(2,60))
        m=int(rng.integers(1,13)) if typ in (1,2,3) else (int(rng.integers(1,22)) if j%3==0 and typ in (4,5) else 0)
        up=int(rng.integers(1,8))
        down=int(rng.integers(2,9)) if typ in (6,7) else int(rng.integers(1,8))
        mode=mode_list[j%len(mode_list)] if typ in (1,2,3,4,5) else 'constant'
        window=['kaiser','hann','hamming','blackman','boxcar'][j%5]
        order=int(rng.integers(12,65))
        centered=j%2
        x=rng.standard_normal(n)
        xi=rng.standard_normal(n)
        h=rng.standard_normal(m) if m>0 else np.empty(0)
        z=x+1j*xi
        if typ==1:
            reference=signal.upfirdn(h,x,up,down,mode=mode)
        elif typ==2:
            reference=signal.upfirdn(h,z,up,down,mode=mode)
        elif typ==3:
            reference=signal.upfirdn(h.astype(complex),z,up,down,mode=mode)
        elif typ==4:
            if m>0:
                reference=signal.resample_poly(x,up,down,window=h,padtype=mode)
            else:
                reference=signal.resample_poly(x,up,down,window=('kaiser',5.0) if window=='kaiser' else window,padtype=mode)
        elif typ==5:
            if m>0:
                reference=signal.resample_poly(z,up,down,window=h,padtype=mode)
            else:
                reference=signal.resample_poly(z,up,down,window=('kaiser',5.0) if window=='kaiser' else window,padtype=mode)
        elif typ==6:
            reference=signal.decimate(x,down,n=order,ftype='fir',zero_phase=bool(centered))
        else:
            reference=signal.decimate(z,down,n=order,ftype='fir',zero_phase=bool(centered))
        cases.append((typ,n,m,up,down,mode,window,order,centered,x,xi,h,reference))

lines=[]
for typ,n,m,up,down,mode,window,order,centered,x,xi,h,reference in cases:
    lines.append(f'{typ} {n} {m} {up} {down} {mode} {window} {order} {centered}')
    lines.append(' '.join(f'{v:.17g}' for v in x))
    lines.append(' '.join(f'{v:.17g}' for v in xi))
    if m: lines.append(' '.join(f'{v:.17g}' for v in h))
proc=subprocess.run([sys.argv[1]],input='\n'.join(lines)+'\n',text=True,capture_output=True,check=False)
if proc.returncode:
    print(proc.stderr[-3000:]); raise SystemExit(proc.returncode)
rows=iter(proc.stdout.splitlines())
counts={t:[0,0] for t in range(1,8)}
worst=(0.0,None)
fail=[]
for j,c in enumerate(cases):
    typ,n,m,up,down,mode,window,order,centered,x,xi,h,ref=c
    hdr=next(rows).split()
    stat,nv=map(int,hdr)
    r=[]
    for _ in range(nv):
        item=next(rows).split()
        if len(item)==1:r.append(float(item[0]))
        else:r.append(float(item[0])+1j*float(item[1]))
    r=np.asarray(r)
    counts[typ][0]+=1
    if stat != 0 or nv != len(ref):
        fail.append((j,typ,'shape/status',stat,nv,len(ref),mode))
        continue
    diff=np.max(np.abs(r-ref)) if nv else 0.
    if diff>worst[0]:worst=(diff,(j,typ,mode,window,up,down,n,m))
    if not np.allclose(r,ref,atol=2e-10,rtol=2e-10):
        fail.append((j,typ,'diff',diff,mode,window,up,down,n,m))
    else: counts[typ][1]+=1
print('counts:',counts,'worst:',worst,'failures:',len(fail))
for item in fail[:20]:print(item)
assert not fail
