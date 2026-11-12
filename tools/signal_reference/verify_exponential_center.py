#!/usr/bin/env python3
"""Compare arbitrary-center periodic exponential/Poisson taper samples against SciPy."""
from __future__ import annotations
import sys, subprocess
import numpy as np
from scipy.signal import windows

def main(executable):
    rng=np.random.default_rng(20261031)
    cases=[]
    for n in [0,1,2,3,4,7,15,32,33,100]:
        for _ in range(20):
            c=rng.uniform(-n-2,2*n+2)
            tau=np.exp(rng.uniform(-3.0,4.0))
            cases.append((n,c,tau,'exponential'))
            cases.append((n,c,tau,'poisson'))
    data=''.join(f'{n} {c:.17g} {t:.17g} {name}\n' for n,c,t,name in cases)
    out=subprocess.check_output([executable],input=data.encode()).decode().splitlines()
    k=0
    worst=0.0
    for n,c,t,name in cases:
        status,got_n=map(int,out[k].split()); k+=1
        if status!=0 or got_n!=n: raise AssertionError('wrong shape or status')
        actual=np.array([float(value) for value in out[k:k+n]])
        k+=n
        ref=windows.exponential(n,center=c,tau=t,sym=False)
        err=float(np.max(np.abs(actual-ref),initial=0))
        worst=max(worst,err)
        if err>3e-14: raise AssertionError(f'mismatch {name} n={n} center={c}: {err}')
    if k!=len(out): raise AssertionError('extra native output')
    print(f'PASS: {len(cases)} SciPy arbitrary-center exponential window cases; max abs {worst:.3g}')
if __name__=='__main__': main(sys.argv[1])
