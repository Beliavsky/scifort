#!/usr/bin/env python3
"""Independent reference test for scifort_optimize_linprog.

Usage: python3 test/compare_linprog_scipy.py /path/to/linprog_reference_driver

Only dense continuous LPs are tested. No HiGHS executable or library is called
from the Fortran implementation; the installed SciPy is an independent oracle.
"""
import sys
import subprocess
import numpy as np
from scipy.optimize import linprog as scipy_linprog


def gen_cases():
    rng = np.random.default_rng(20261110)
    cases = []
    for k in range(850):
        n = int(rng.integers(1, 9))
        mu = int(rng.integers(0, 13))
        me = int(rng.integers(0, min(n + 4, 9)))
        c = rng.normal(size=n)
        lo = np.zeros(n)
        hi = np.zeros(n)
        x0 = np.zeros(n)
        kind = k % 5
        for j in range(n):
            z = int(rng.integers(0, 6)) if kind == 0 else int(rng.integers(0, 4))
            if z == 0:
                lo[j], hi[j] = -np.inf, np.inf
                x0[j] = rng.normal()
            elif z == 1:
                lo[j], hi[j] = rng.uniform(-5, 0), np.inf
                x0[j] = lo[j] + rng.uniform(0, 5)
            elif z == 2:
                lo[j], hi[j] = -np.inf, rng.uniform(0, 5)
                x0[j] = hi[j] - rng.uniform(0, 5)
            elif z == 3:
                lo[j], hi[j] = rng.uniform(-5, 0), rng.uniform(0, 5)
                x0[j] = rng.uniform(lo[j],hi[j])
            elif z == 4:
                lo[j] = rng.uniform(-3, 3)
                hi[j] = lo[j]
                x0[j] = lo[j]
            else:
                lo[j], hi[j] = 0, np.inf
                x0[j] = rng.uniform(0, 3)
        if kind == 3:
            lo[:] = -rng.uniform(0, 3, n)
            hi[:] = rng.uniform(0, 3, n)
            x0 = rng.uniform(lo, hi)
        if kind == 4:
            lo[:] = 0
            hi[:] = np.inf
            x0 = rng.uniform(0, 4, n)
        aub = rng.normal(size=(mu,n))
        bub = aub@x0 + rng.uniform(0, 3, mu)
        aeq = rng.normal(size=(me,n))
        beq = aeq@x0
        if me > 1 and k % 9 == 0:
            aeq[-1] = aeq[0]
            beq[-1] = beq[0]
        if me > 1 and k % 19 == 0:
            aeq[-1] = 2*aeq[0]
            beq[-1] = 2*beq[0]
        if k % 11 == 0 and mu > 1:
            aub[-1] = aub[0]
            bub[-1] = bub[0]
        if k % 15 == 0:
            c[:] = 0.0
        cases.append((c, aub, bub, aeq, beq, lo, hi, f'feasible_{k}'))
    # Known infeasible and unbounded LPs, plus degenerate test cases.
    for k in range(100):
        n = int(rng.integers(1, 8))
        c = rng.normal(size=n)
        lo = np.zeros(n)
        hi = np.full(n,np.inf)
        aub = np.zeros((2,n))
        aub[0,0], aub[1,0] = 1, -1
        bub = np.array([0.0, -1.0-float(k%3)])
        cases.append((c,aub,bub,np.zeros((0,n)),np.zeros(0),lo,hi,f'infeasible_{k}'))
    for k in range(100):
        n = int(rng.integers(1, 8))
        c = rng.uniform(-3,3,n)
        c[0] = -1-float(k%3)
        lo = np.zeros(n)
        hi = np.full(n,np.inf)
        cases.append((c,np.zeros((0,n)),np.zeros(0),np.zeros((0,n)),np.zeros(0),lo,hi,f'unbounded_{k}'))
    for k in range(80):
        n = int(rng.integers(1, 8))
        c = rng.normal(size=n)
        lo = -rng.uniform(0,5,n)
        hi = rng.uniform(0,5,n)
        x = np.zeros(n)
        mu, me = 12, n+4
        aub = rng.choice([-1.0,0.0,1.0], size=(mu,n))
        bub = aub@x
        aeq = rng.choice([-1.0,0.0,1.0],size=(me,n))
        beq = aeq@x
        cases.append((c,aub,bub,aeq,beq,lo,hi,f'degenerate_{k}'))
    for k in range(80):
        n = int(rng.integers(1, 8))
        c = rng.normal(size=n)
        lo = -rng.uniform(0,5,n)
        hi = rng.uniform(0,5,n)
        x0 = rng.uniform(lo,hi)
        aeq = rng.normal(size=(1,n))
        beq = aeq@x0
        aeq = np.vstack([aeq,aeq,2*aeq])
        beq = np.concatenate([beq,beq,2*beq])
        aub = rng.normal(size=(4,n))
        bub = aub@x0+rng.uniform(0,3,4)
        cases.append((c,aub,bub,aeq,beq,lo,hi,f'redundant_{k}'))
    # Conditioning checks: equivalent LPs with independently scaled rows and
    # a scaled objective. Avoid pathologically tiny feasibility tolerances in
    # the independent HiGHS oracle by restricting row factors to 1e-3..1e3.
    for k,source in enumerate(cases[:240]):
        c, aub, bub, aeq, beq, lo, hi, _ = source
        su = 10.0**rng.uniform(-3,3,len(bub))
        se = 10.0**rng.uniform(-3,3,len(beq))
        cs = 10.0**rng.uniform(-3,3)
        cases.append((c*cs, aub*su[:,None], bub*su, aeq*se[:,None],
                      beq*se, lo, hi, f'scaled_{k}'))
    return cases


def to_lines(case):
    c, aub, bub, aeq, beq, lo, hi, _ = case
    def line(a):
        return ' '.join('Infinity' if z==np.inf else '-Infinity' if z==-np.inf else f'{z:.17g}' for z in a) + '\n'
    s = f'{len(c)} {len(bub)} {len(beq)}\n' + line(c) + line(lo) + line(hi)
    for row in aub: s += line(row)
    if len(bub): s += line(bub)
    for row in aeq: s += line(row)
    if len(beq): s += line(beq)
    return s


def main(path):
    cases = gen_cases()
    text = ''.join(to_lines(case) for case in cases) + '-1 0 0\n'
    process = subprocess.run([path], input=text, text=True, capture_output=True, check=True)
    output = process.stdout.splitlines()
    assert len(output) == 2*len(cases), (len(output),len(cases),process.stderr[:1000])
    max_obj = 0.0
    max_feas = 0.0
    count = 0
    bad = []
    counts = {}
    for k,case in enumerate(cases):
        c,aub,bub,aeq,beq,lo,hi,tag = case
        status, val, nit = output[2*k].split()
        status, val, nit = int(status), float(val), int(nit)
        xraw = output[2*k+1].split()
        x = None if 'no-solution' in xraw else np.array([float(z) for z in xraw])
        try:
            ref = scipy_linprog(c, A_ub=aub, b_ub=bub, A_eq=aeq, b_eq=beq,
                              bounds=list(zip(lo,hi)), method='highs',
                              options={'presolve':False, 'primal_feasibility_tolerance':1e-8,
                                       'dual_feasibility_tolerance':1e-8})
        except Exception as exc:
            bad.append((tag,'SciPy exception',str(exc)))
            continue
        counts[ref.status] = counts.get(ref.status,0)+1
        if status!=ref.status:
            bad.append((tag,'status',status,ref.status,val,ref.fun,nit))
            continue
        if status == 0:
            if x is None or len(x)!=len(c):
                bad.append((tag,'missing x'))
                continue
            scale = max(1,abs(ref.fun))
            gap = abs(val-ref.fun)/scale
            feas = 0.
            if len(bub): feas = max(feas,np.max(np.maximum(aub@x-bub,0)))
            if len(beq): feas = max(feas,np.max(np.abs(aeq@x-beq)))
            feas = max(feas,np.max(np.maximum(lo-x,0)),np.max(np.maximum(x-hi,0)))
            bound_scale = max(1,np.linalg.norm(bub,ord=np.inf) if len(bub) else 0,
                              np.linalg.norm(beq,ord=np.inf) if len(beq) else 0)
            feas/=bound_scale
            max_obj = max(max_obj,gap)
            max_feas = max(max_feas,feas)
            if gap>1e-6 or feas>1e-7:
                bad.append((tag,'numerics',gap,feas,val,ref.fun,nit))
        count += 1
    print('SciFort/SciPy cases:',count,'/',len(cases),'SciPy statuses:',counts)
    print('max relative objective gap:',max_obj,'max normalized infeasibility:',max_feas)
    print('failures:',len(bad))
    if bad:
        for b in bad[:30]: print('  ',b)
        raise SystemExit(1)

if __name__ == '__main__':
    main(sys.argv[1])
