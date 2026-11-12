#!/usr/bin/env python3
"""Independent SciPy 1.17 comparisons for real/complex PPoly/BPoly.extend.

Runs checked and optimized GNU Fortran builds; no nonstandard Fortran dependencies.
SciPy 1.17 has a known issue for prepending PPoly segments with a different
polynomial order; in that case an explicitly constructed reference is used.
"""
from __future__ import annotations
import argparse
import subprocess
import tempfile
from pathlib import Path
import numpy as np
from scipy.interpolate import PPoly, BPoly

ROOT=Path(__file__).resolve().parents[2]
SOURCES=['src/kinds.f90','src/interpolate_piecewise.f90',
         'src/interpolate_bspline.f90','src/interpolate_regular_grid.f90','src/interpolate_barycentric.f90', 'src/interpolate_global.f90', 'src/interpolate_kdtree.f90', 'src/interpolate_rbf.f90', 'src/interpolate_legacy.f90', 'src/interpolate_scattered.f90', 'src/interpolate_clough_tocher.f90', 'src/interpolate_fitpack_bivariate.f90', 'src/interpolate_fitpack_compat.f90', 'src/interpolate.f90']

def fmt(a):
    return ' '.join(f'{float(z):.17g}' for z in np.ravel(a,order='F'))

def compile_driver(wd, flag):
    command=['gfortran','-std=f2018','-ffree-line-length-none',*flag,'-J'+str(wd),'-I'+str(wd),
             *(str(ROOT/s) for s in SOURCES),
             str(ROOT/'tools/interpolate_reference/piecewise_extend_driver.f90'),'-o',str(wd/'piecewise_extend')]
    cp=subprocess.run(command, cwd=wd, capture_output=True,text=True)
    if cp.returncode: raise RuntimeError(cp.stderr[-7000:])
    return wd/'piecewise_extend'

def make_cases(seed, count):
    rng=np.random.default_rng(seed)
    out=[]
    for i in range(count):
        kind=1+i%4
        deg0=i%11
        deg1=(i*7+3)%12
        if i%13 == 0: deg1=deg0
        seg0=1+i%5
        seg1=1+(i//5)%5
        nch=1+(i//7)%4
        oldx=np.r_[0,np.cumsum(rng.uniform(.4,1.25,seg0))]+rng.uniform(-3,3)
        append=(i%4<2)
        if append:
            newx=oldx[-1]+np.cumsum(rng.uniform(.15,1.3,seg1))
        else:
            newx=oldx[0]-np.cumsum(rng.uniform(.15,1.3,seg1))[::-1]
        if i%8 == 5:
            oldx=oldx[::-1].copy()
            # Append/prepend selection is in index order, not physical direction.
            if append:
                newx=oldx[-1]-np.cumsum(rng.uniform(.15,1.3,seg1))
            else:
                newx=oldx[0]+np.cumsum(rng.uniform(.15,1.3,seg1))[::-1]
        c0=rng.normal(size=(deg0+1,seg0,nch))*0.7
        c1=rng.normal(size=(deg1+1,seg1,nch))*0.7
        if kind>=3:
            c0=c0+1j*rng.normal(size=c0.shape)*0.7
            c1=c1+1j*rng.normal(size=c1.shape)*0.7
        fullx=np.r_[oldx,newx] if append else np.r_[newx,oldx]
        q=rng.uniform(min(fullx)-.3,max(fullx)+.3,14)
        q=np.r_[q,fullx[0],fullx[-1],oldx[0],oldx[-1]]
        out.append((kind,oldx,c0,newx,c1,q))
    return out

def serial(cases):
    out=[str(len(cases))]
    for kind,x,c,xn,cn,q in cases:
        out += [f'{kind} {len(c)-1} {len(x)-1} {len(cn)-1} {len(xn)} {c.shape[-1]} {len(q)}',
                fmt(x),fmt(c.real),fmt(c.imag),fmt(xn),fmt(cn.real),fmt(cn.imag),fmt(q)]
    return '\n'.join(out)+'\n'

def expected_polynomial(kind,x,c,xn,cn):
    cls=PPoly if kind%2 else BPoly
    original=cls(c,x)
    append = (xn[0] >= x[-1]) if x[-1]>x[0] else (xn[0]<=x[-1])
    if kind%2:
        k=max(c.shape[0],cn.shape[0])
        old=np.pad(c,((k-c.shape[0],0),(0,0),(0,0)))
        new=np.pad(cn,((k-cn.shape[0],0),(0,0),(0,0)))
    else:
        k=max(c.shape[0],cn.shape[0])
        old=BPoly._raise_degree(c,k-c.shape[0])
        new=BPoly._raise_degree(cn,k-cn.shape[0])
    ec=np.concatenate([old,new],axis=1) if append else np.concatenate([new,old],axis=1)
    ex=np.r_[x,xn] if append else np.r_[xn,x]
    constructed=cls(ec,ex)
    try:
        original.extend(cn,xn)
    except (ValueError,TypeError):
        # A SciPy 1.17 prepend mismatch may raise due incorrect leading index.
        if kind%2==1 and not append and c.shape[0]!=cn.shape[0]:
            return constructed,True
        raise
    if not (np.array_equal(original.x,ex) and np.allclose(original.c,ec)):
        # SciPy 1.17 prepend with mismatched orders can write the wrong data.
        if kind%2==1 and not append and c.shape[0]!=cn.shape[0]:
            return constructed,True
        raise AssertionError('independent assembled coefficients disagree with SciPy.extend')
    return original,False

def compare(binary,cases):
    cp=subprocess.run([str(binary)],input=serial(cases),text=True,capture_output=True)
    if cp.returncode: raise RuntimeError('Fortran failed: '+cp.stderr[-3000:])
    lines=iter(cp.stdout.splitlines())
    counter=0; maxerr=0.; scipy_workarounds=0
    for index,(kind,x,c,xn,cn,q) in enumerate(cases,1):
        expected,workaround=expected_polynomial(kind,x,c,xn,cn)
        scipy_workarounds+=workaround
        head=next(lines).split()
        assert head[:2]==['CASE',str(index)],(index,head)
        k=int(head[2]); m=int(head[3]); nch=c.shape[-1]
        actual_knots=np.fromstring(next(lines),sep=' ')
        assert np.array_equal(actual_knots,expected.x),(index,actual_knots,expected.x)
        real_coeff=np.fromstring(next(lines),sep=' ').reshape((k,m,nch),order='F')
        actual_coeff=real_coeff
        if kind>=3:
            imag_coeff=np.fromstring(next(lines),sep=' ').reshape((k,m,nch),order='F')
            actual_coeff=real_coeff+1j*imag_coeff
        assert actual_coeff.shape == expected.c.shape,(index,actual_coeff.shape,expected.c.shape)
        np.testing.assert_allclose(actual_coeff,expected.c,rtol=2e-12,atol=2e-12,err_msg=f'coefficients case={index}')
        for deriv in range(3):
            for val in q:
                got=np.fromstring(next(lines),sep=' ')
                if kind>=3:
                    imag=np.fromstring(next(lines),sep=' ')
                    got=got+1j*imag
                should=(np.zeros(nch, dtype=expected.c.dtype) if (kind%2==0 and deriv>=len(expected.c))
                        else np.asarray(expected(val,nu=deriv)).reshape(-1))
                np.testing.assert_allclose(got,should,rtol=2e-9,atol=2e-9,
                                           err_msg=f'values case={index}, deriv={deriv}, x={val}')
                maxerr=max(maxerr,float(np.max(np.abs(got-should))))
                counter+=len(got)
    assert next(lines,None) is None
    return counter,maxerr,scipy_workarounds

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--seed',type=int,default=20261031)
    parser.add_argument('--cases',type=int,default=480)
    args=parser.parse_args()
    cases=make_cases(args.seed,args.cases)
    for name,flags in [('checked',['-O0','-g','-fcheck=all','-fbacktrace']),('optimized',['-O3'])]:
        with tempfile.TemporaryDirectory() as folder:
            wd=Path(folder)
            exe=compile_driver(wd,flags)
            n,err,workarounds=compare(exe,cases)
        print(f'{name}: PASS {len(cases)} cases / {n} values; max_abs_error={err:.3g}; '+
              f'prepending-mismatched-order SciPy workarounds={workarounds}')

if __name__=='__main__': main()
