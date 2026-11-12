"""Seeded numerical parity test for SciFort versus scipy.signal.deconvolve."""
import subprocess
import os
from pathlib import Path
import numpy as np
from scipy.signal import deconvolve
def executable(name):
    p = Path(__file__).resolve().parent.parent / 'build' / 'signal_reference' / name
    return p.with_suffix('.exe') if os.name == 'nt' else p

rng=np.random.default_rng(20261019)
cases=[]
for i in range(600):
    n=int(rng.integers(1,21)); m=int(rng.integers(1, n+5)); typ=int(rng.integers(0,4))
    x=rng.normal(size=n)+1j*rng.normal(size=n)*(typ>1)
    d=rng.normal(size=m)+1j*rng.normal(size=m)*(typ%2==1)
    d[0]=d[0]+(2.0+0.0j if d[0].real>=0 else -2.0-0j)
    cases.append((n,m,typ,x,d))
lines=[str(len(cases))]
for n,m,typ,x,d in cases:
    lines.append(f'{n} {m} {typ}')
    lines.extend(f'{v.real:.17g} {v.imag:.17g}' for v in x)
    lines.extend(f'{v.real:.17g} {v.imag:.17g}' for v in d)
p=subprocess.run([str(executable('ref_deconvolve'))],input='\n'.join(lines)+'\n',text=True,capture_output=True,check=True)
outs=iter(p.stdout.splitlines()); failures=[]; mx=0.0
for i,(n,m,typ,x,d) in enumerate(cases):
    code,nq,nr=map(int,next(outs).split())
    q=np.array([complex(*map(float,next(outs).split())) for _ in range(nq)])
    r=np.array([complex(*map(float,next(outs).split())) for _ in range(nr)])
    qe,re=deconvolve(x if typ>1 else x.real,d if typ%2 else d.real)
    qe=np.asarray(qe); re=np.asarray(re)
    err=max(np.max(abs(q-qe),initial=0),np.max(abs(r-re),initial=0))
    mx=max(mx,err)
    # Numerical cancellation can increase absolute error for ill-conditioned polynomials.
    thresh=2e-8 * max(1.0,np.max(abs(qe),initial=0),np.max(abs(re),initial=0))
    if code!=0 or err>thresh: failures.append((i,n,m,typ,err,thresh))
print(f'deconvolve: {len(cases)-len(failures)}/{len(cases)} pass, max abs {mx:.3g}')
print('failed',failures[:10])
if failures: raise SystemExit(1)
