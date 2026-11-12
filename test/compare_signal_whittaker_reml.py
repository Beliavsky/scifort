#!/usr/bin/env python3
"""Compare Whittaker-Henderson REML penalty and output to an independent
profiled-likelihood implementation using SciPy linear algebra and optimization.

SciPy's public whittaker_henderson was added in 1.18, so this reference also
runs with SciPy 1.17 by using the documented profiled REML criterion directly.
Usage: python test/compare_signal_whittaker_reml.py EXECUTABLE
"""
import itertools
import pathlib
import subprocess
import sys
import tempfile

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize_scalar


def reference(y, w, p):
    n = len(y)
    ok = w > 0
    yzero = np.where(ok,y,0)
    d = np.diff(np.eye(n),n=p,axis=0)
    q = d.T@d
    positive = int(n-p)

    def cost(t, return_x=False):
        lam = np.exp(t)
        a = np.diag(w)+lam*q
        try:
            c, low = cho_factor(a,lower=True,check_finite=False)
            x = cho_solve((c,low),w*yzero,check_finite=False)
        except np.linalg.LinAlgError:
            return 1e300
        r2 = np.sum(w[ok]*(yzero[ok]-x[ok])**2)+lam*np.sum((d@x)**2)
        if not np.isfinite(r2) or r2 <= 0:
            return 1e300
        if return_x:
            return x
        return positive*(np.log(r2/positive)-t)+2*np.log(np.diag(c)).sum()

    grid = np.linspace(-24,24,33)
    i = np.argmin([cost(t) for t in grid])
    left,right = grid[max(0,i-1)],grid[min(len(grid)-1,i+1)]
    sol = minimize_scalar(cost,bounds=(left,right),method='bounded',options={'xatol':2e-9})
    return np.exp(sol.x), cost(sol.x,return_x=True)


def run(exe):
    cases = list(itertools.product((18,35,55),(1,2,3),(0,1,2)))
    with tempfile.TemporaryDirectory() as tmp:
        d=pathlib.Path(tmp)
        (d/'whittaker_reml_cases.txt').write_text(f'{len(cases)}\n'+
            '\n'.join(' '.join(map(str,x)) for x in cases)+'\n')
        run = subprocess.run([str(pathlib.Path(exe).resolve())],cwd=d,
                             capture_output=True,text=True)
        if run.returncode:
            raise RuntimeError(f'Fortran solver failed: {run.stdout}\n{run.stderr}')
        err_x = 0.
        err_loglam = 0.
        count=0
        with (d/'whittaker_reml_results.txt').open() as fh:
            for i,(n,p,mode) in enumerate(cases,1):
                h=list(map(int,fh.readline().split()))
                assert h==[i,n,p,mode],(h,i,n,p,mode)
                lam=float(fh.readline())
                x=np.array([float(fh.readline()) for _ in range(n)])
                t=np.arange(n)
                y=0.03*t+np.sin(.12*t)+.23*np.cos(.51*t)+.2*np.sin(2.1*t+.15*i)
                w=np.ones(n)
                for k in range(1,n+1):
                    if mode==1 and (k+i)%7==0:w[k-1]=0.
                    if mode==2 and (k+i)%7==0:w[k-1]=.3
                    if mode==2 and (k+i)%11==0:w[k-1]=0.
                y[w==0]=np.nan
                targetlam,targetx=reference(y,w,p)
                logerr=abs(np.log(lam/targetlam))
                xerr=np.max(abs(x-targetx))
                err_loglam=max(logerr,err_loglam)
                err_x=max(xerr,err_x)
                count+=n
                if logerr>2e-2 or xerr>3e-6:
                    raise AssertionError(f'case {i}, n={n}, p={p}, mode={mode}: '
                                         f'lambda {lam} vs {targetlam} (logerr {logerr}); '
                                         f'error {xerr}')
            assert not fh.readline()
    print(f'Whittaker REML: {len(cases)} cases, {count} fitted values, '
          f'max abs x error={err_x:.3e}, max |log(lambda ratio)|={err_loglam:.3e}')


if __name__=='__main__':
    if len(sys.argv)!=2:
        raise SystemExit('Usage: compare_signal_whittaker_reml.py EXECUTABLE')
    run(sys.argv[1])
