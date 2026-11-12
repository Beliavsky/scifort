import numpy as np
from scipy.signal import envelope
from pathlib import Path
import subprocess
import os
def executable(name):
    p = Path(__file__).resolve().parent.parent / 'build' / 'signal_reference' / name
    return p.with_suffix('.exe') if os.name == 'nt' else p

rng=np.random.default_rng(9001)
cases=[]
for i in range(380):
    nr=int(rng.choice([1,2,3,4,8,9,16])); nt=int(rng.choice([1,2,3,7,8,9,16])); ax=int(rng.integers(1,3)); n=[nr,nt][ax-1]
    nout=int(rng.choice([1,2,max(1,n//2),n,n+3,2*n]))
    lo=int(rng.integers(-n//2 if n%2==0 else -(n//2),(n+1)//2))
    # -n//2 in Python means floor of negative number; use -(n//2).
    lo=int(rng.integers(-(n//2),(n+1)//2))
    hi=int(rng.integers(lo+1,(n+1)//2+1))
    cx=int(rng.integers(0,2)); square=int(rng.integers(0,2)); res=int(rng.integers(1,3))
    z=rng.normal(size=(nr,nt))+1j*rng.normal(size=(nr,nt))*(cx!=0)
    cases.append((nr,nt,ax,nout,lo,hi,cx,square,res,z))
lines=[str(len(cases))]
for nr,nt,ax,nout,lo,hi,cx,square,res,z in cases:
    lines.append(f'{nr} {nt} {ax} {nout} {lo} {hi} {cx} {square} {res}')
    lines.extend(f'{v.real:.17g} {v.imag:.17g}' for v in z.ravel(order='F'))
p=subprocess.run([str(executable('ref_envelope_2d'))],input='\n'.join(lines)+'\n',text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True)
outs=iter(p.stdout.splitlines()); max_err=0; failures=[]
for idx,(_,_,ax,nout,lo,hi,cx,square,res,z) in enumerate(cases):
    status,r,c,rr,rc=map(int,next(outs).split())
    got=np.array([list(map(float,next(outs).split())) for j in range(r*c)]).reshape((r,c,3),order='F')
    reference=envelope(z if cx else z.real,bp_in=(lo,hi),n_out=nout,axis=ax-1,squared=bool(square),residual=['none','all','lowpass'][res])
    y=reference[0]; expected_res=reference[1]
    err=max(np.max(abs(got[:,:,0]-y)),np.max(abs(got[:,:,1]+1j*got[:,:,2]-expected_res)))
    max_err=max(max_err,err)
    if err>1e-9 or status !=0: failures.append((idx,err,(r,c),y.shape))
print('2D cases',len(cases),'max err',max_err,'fail',failures[:12]);
if failures: raise SystemExit(1)
