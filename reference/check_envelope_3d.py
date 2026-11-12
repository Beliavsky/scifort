"""Seeded rank-3 randomized comparison with SciPy envelope."""
from pathlib import Path
import numpy as np
from scipy.signal import envelope
import subprocess
import os
def executable(name):
    p = Path(__file__).resolve().parent.parent / 'build' / 'signal_reference' / name
    return p.with_suffix('.exe') if os.name == 'nt' else p

rng=np.random.default_rng(1230992)
cases=[]
for t in range(220):
    shape=tuple(int(rng.integers(1,7)) for j in range(3))
    ax=int(rng.integers(1,4)); n=shape[ax-1]
    nout=int(rng.choice([1,2,n,2*n,n+1]))
    lo=int(rng.integers(-(n//2),(n+1)//2))
    hi=int(rng.integers(lo+1,(n+1)//2+1))
    cx=int(rng.integers(0,2)); square=int(rng.integers(0,2)); res=int(rng.integers(0,3))
    z=rng.normal(size=shape)+1j*rng.normal(size=shape)*(cx!=0)
    cases.append((shape,ax,nout,lo,hi,cx,square,res,z))
lines=[str(len(cases))]
for shape,ax,nout,lo,hi,cx,square,res,z in cases:
    lines.append(' '.join(map(str,(*shape,ax,nout,lo,hi,cx,square,res))))
    lines.extend(f'{v.real:.17g} {v.imag:.17g}' for v in z.ravel(order='F'))
p=subprocess.run([str(executable('ref_envelope_3d'))],input='\n'.join(lines)+'\n',text=True,capture_output=True,check=True)
outs=iter(p.stdout.splitlines()); failed=[]; max_err=0
for index,(shape,ax,nout,lo,hi,cx,square,res,z) in enumerate(cases):
    data=list(map(int,next(outs).split()))
    status=data[0]; out_shape=tuple(data[1:4]); res_shape=tuple(data[4:7])
    got=np.array([list(map(float,next(outs).split())) for j in range(np.prod(out_shape))]).reshape((*out_shape,3),order='F')
    expected=envelope(z if cx else z.real,bp_in=(lo,hi),n_out=nout, squared=bool(square),residual=[None,'all','lowpass'][res],axis=ax-1)
    y=expected if res==0 else expected[0]
    residue=np.zeros(out_shape) if res==0 else expected[1]
    err=max(np.max(abs(got[:,:,:,0]-y)),np.max(abs(got[:,:,:,1]+1j*got[:,:,:,2]-residue)))
    max_err=max(max_err,err)
    if err>1e-9 or status!=0 or (res==0 and res_shape!=(0,0,0)):
        failed.append((index,status,out_shape,err,res_shape))
print('3D cases',len(cases),'max abs',max_err,'failed',failed[:10])
if failed: raise SystemExit(1)
