#!/usr/bin/env python3
"""Independent SciPy comparison for continuous band-stop order objectives.

Compare against scipy.signal.band_stop_obj in the well-conditioned regime.
For very large gstop/small gpass, SciPy's direct 10**x-1 and K(1-m)
expressions can lose precision; compare to a stable formulation using
scipy.special.ellipkm1 (a separately implemented elliptic-integral kernel).
"""
from __future__ import annotations
import math
import sys
import subprocess
import numpy as np
from scipy import signal, special


def stable_objective(wp, ind, p, s, gp, gs, mode):
    edges=p.copy()
    edges[ind]=wp
    nat=np.min(np.abs(s*(edges[0]-edges[1])/(s*s-edges[0]*edges[1])))
    logratio=np.log(np.expm1(gs*np.log(10)/10))-np.log(np.expm1(gp*np.log(10)/10))
    if mode=='butter':
        return logratio/(2*np.log(nat))
    if mode=='cheby':
        return np.arccosh(np.exp(0.5*logratio))/np.arccosh(nat)
    m0=1/(nat*nat)
    if gp==gs:
        return 0.0
    m1=np.exp(-logratio)
    return (special.ellipk(m0)/special.ellipkm1(m0))*\
           (special.ellipkm1(m1)/special.ellipk(m1))


def main(executable: str) -> None:
    rng=np.random.default_rng(20261031)
    cases=[]
    for _ in range(180):
        # Valid analog band-stop design with both stop edges inside pass edges.
        a = rng.uniform(0.1,1.5)
        b = rng.uniform(2.8,7.5)
        s0 = rng.uniform(a*1.1,(a+b)/2)
        s1 = rng.uniform(s0+(b-s0)*0.12,b*0.9)
        if not a<s0<s1<b:
            continue
        passb=np.array([a,b])
        stopb=np.array([s0,s1])
        gp = 10**rng.uniform(-2,0.7)
        gs = gp + 10**rng.uniform(0,2)
        for i in [0,1]:
            for mode in ['butter','cheby','ellip']:
                wp=passb[i] * rng.uniform(0.91,1.09)
                modified=passb.copy()
                modified[i]=wp
                nat=np.min(abs(stopb*(modified[0]-modified[1]) / (stopb**2 - np.prod(modified))))
                if nat<=1.0:
                    continue
                cases.append((wp,i,passb,stopb,gp,gs,mode))
    # Difficult but valid cases: tiny pass ripple and high stop attenuation.
    for gp,gs in [(0.0000001,40),(0.0001,150),(0.001,0.00101),(0.05,80),(1,40),(3,25),(3,3)]:
        for mode in ['butter','cheby','ellip']:
            cases.append((0.5,0,np.array([0.5,4.0]),np.array([1.5,2.5]),gp,gs,mode))
    # SciPy permits negative Butterworth objective values in trial regions.
    cases.append((2.,1,np.array([1.,3.]),np.array([0.5,4.]),3.,30.,'butter'))
    assert len(cases)>400
    text=''.join(f'{wp:.17g} {idx} {p[0]:.17g} {p[1]:.17g} {s[0]:.17g} '
                 f'{s[1]:.17g} {gp:.17g} {gs:.17g} {mode}\n'
                 for wp,idx,p,s,gp,gs,mode in cases)
    out = subprocess.check_output([executable], input=text.encode()).decode().splitlines()
    if len(out)!=len(cases):
        raise AssertionError(f'output count {len(out)} != {len(cases)}')
    worst_abs=0.0
    worst_rel=0.0
    worst_scipy=0.0
    scipy_count=0
    for i,(line,case) in enumerate(zip(out,cases)):
        actual_s,status_s=line.split()
        actual=float(actual_s)
        if int(status_s)!=0:
            raise AssertionError(f'fortran rejected case {i}: {case}')
        wp,ind,p,s,gp,gs,mode=case
        expected=float(stable_objective(wp,ind,p,s,gp,gs,mode))
        abserr=abs(actual-expected)
        relerr=abserr/max(1.0,abs(expected))
        worst_abs=max(worst_abs,abserr)
        worst_rel=max(worst_rel,relerr)
        if not math.isfinite(actual) or not math.isfinite(expected) or relerr>5.e-11:
            raise AssertionError(f'stable mismatch case {i}: actual {actual}, expected {expected}: {case}')
        if gp>=0.01 and gs<=60.0:
            scipy_expected=float(signal.band_stop_obj(wp,ind,p,s,gp,gs,mode))
            scipy_rel=abs(actual-scipy_expected)/max(1.0,abs(scipy_expected))
            worst_scipy=max(worst_scipy,scipy_rel)
            if scipy_rel>5e-9:
                raise AssertionError(f'SciPy mismatch case {i}: actual {actual}, SciPy {scipy_expected}: {case}')
            scipy_count+=1
    print(f'PASS: {len(cases)} band-stop reference cases; max stable absolute error {worst_abs:.3g}; '
          f'max stable relative error {worst_rel:.3g}')
    print(f'PASS: {scipy_count} direct scipy.signal.band_stop_obj cases; '
          f'max relative error {worst_scipy:.3g}')

if __name__=='__main__':
    main(sys.argv[1])
