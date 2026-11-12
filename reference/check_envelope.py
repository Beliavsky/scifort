"""Seeded independent comparisons against scipy.signal.envelope."""
import numpy as np
from scipy.signal import envelope
from pathlib import Path
import subprocess
import os
def executable(name):
    p = Path(__file__).resolve().parent.parent / 'build' / 'signal_reference' / name
    return p.with_suffix('.exe') if os.name == 'nt' else p

rng=np.random.default_rng(314159)
cases=[]
for index in range(700):
    n = int(rng.choice([1,2,3,4,5,6,7,8,9,11,12,15,16,17,29,32,33,64,65]))
    nout=int(rng.choice([1,2,3,max(1,n//2),n,n+1,2*n+1]))
    low=int(rng.integers(-(n//2), (n+1)//2))
    high=int(rng.integers(low+1,(n+1)//2+1))
    cx=int(rng.integers(0,2)); square=int(rng.integers(0,2)); res=int(rng.integers(0,3))
    z=rng.normal(size=n)+1j*rng.normal(size=n)*(cx!=0)
    cases.append((n,nout,low,high,cx,square,res,z))
lines=[str(len(cases))]
for n,nout,low,high,cx,square,res,z in cases:
    lines.append(f'{n} {nout} {low} {high} {cx} {square} {res}')
    lines.extend(f'{w.real:.17g} {w.imag:.17g}' for w in z)
p=subprocess.run([str(executable('ref_envelope'))],input='\n'.join(lines)+'\n',text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True)
outs=iter(p.stdout.splitlines()); max_err=0; failures=[]
for idx,(n,nout,low,high,cx,square,res,z) in enumerate(cases):
    status, n_y, n_res=map(int,next(outs).split())
    got=np.array([list(map(float,next(outs).split())) for i in range(n_y)])
    a=envelope(z if cx else z.real, bp_in=(low,high),n_out=nout,squared=bool(square),residual=[None,'all','lowpass'][res]); y=a if res==0 else a[0]; r=np.zeros(nout) if res==0 else a[1]
    err=max(np.max(abs(got[:,0]-y)) if n_y else np.inf, np.max(abs(got[:,1]+1j*got[:,2]-r)) if n_res else (0 if res==0 else np.inf))
    max_err=max(err,max_err)
    if err>1e-9 or status!=0: failures.append((idx,n,nout,low,high,cx,square,res,err,np.max(abs(got[:,0]-y))));
print('cases',len(cases),'max error',max_err,'failures',len(failures))
print('examples',failures[:18])
if failures: raise SystemExit(1)
