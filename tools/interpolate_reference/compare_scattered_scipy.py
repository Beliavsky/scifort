#!/usr/bin/env python3
"""Compare native SciFort scattered interpolation against SciPy 1.17.

NearestNDInterpolator is exercised in 2–6 dimensions, while LinearNDInterpolator
is tested in 2D with real/complex, multichannel, unit-cube rescaling,
nonuniform coordinate scales, convex-hull fill, and sample-site queries.
No FPM, Qhull, external Fortran libraries, or network access required.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator, griddata

ROOT=Path(__file__).resolve().parents[2]


def run(cmd,cwd):
    r=subprocess.run(cmd,cwd=cwd,text=True,capture_output=True)
    if r.returncode:
        raise RuntimeError(f"Command failed: {' '.join(map(str,cmd))}\n{r.stdout[-500:]}\n{r.stderr[-2000:]}")
    return r.stdout


def build(tmp,flags):
    base=['gfortran','-std=f2018','-ffree-line-length-none',*flags]
    for name in ('kinds','interpolate_kdtree','interpolate_scattered'):
        run([*base,f'-J{tmp}',f'-I{tmp}','-c',str(ROOT/'src'/f'{name}.f90'),'-o',str(tmp/f'{name}.o')],tmp)
    exe=tmp/'scattered_driver'
    run([*base,f'-I{tmp}',str(ROOT/'tools/interpolate_reference/scattered_driver.f90'),
        str(tmp/'interpolate_scattered.o'),str(tmp/'interpolate_kdtree.o'),str(tmp/'kinds.o'),'-o',str(exe)],tmp)
    return exe


def run_case(exe,mode,x,y,q,fill,rescale,complex_data):
    n,d=x.shape
    c=y.shape[1]
    m=len(q)
    lines=[f'{mode} {n} {d} {c} {m} {int(rescale)} {int(complex_data)}',
           f'{fill.real:.17g} {fill.imag:.17g}']
    lines+=[' '.join(f'{v:.17g}' for v in row) for row in x]
    lines += [f'{v.real:.17g} {v.imag:.17g}' for row in y for v in row]
    lines+=[' '.join(f'{v:.17g}' for v in row) for row in q]
    r=subprocess.run([str(exe)],input='\n'.join(lines)+'\n',text=True,capture_output=True)
    if r.returncode or not r.stdout:
        raise RuntimeError(f'driver failed: {r.stdout[-500:]} {r.stderr[-500:]}')
    lines=r.stdout.splitlines()
    if lines[0].strip()!='0':
        raise AssertionError(f'native returned status {lines[0]}')
    arr=np.array([list(map(float,t.split())) for t in lines[1:]])
    if arr.shape!=(m*c,2):
        raise AssertionError(f'invalid shape {arr.shape} vs {(m*c,2)}')
    return (arr[:,0]+1j*arr[:,1]).reshape(m,c)


def case(rng,k):
    mode= k%2 # alternate nearest / linear
    dim = (2+k//2%5) if mode==0 else 2
    n = (8+k%26) if mode==0 else (8+k%28)
    channels=1+k//7%5
    complex_data=bool(k//3%2)
    rescale=bool(k//5%2)
    scale = np.exp(rng.uniform(-2, 2, size=(1,dim)))
    offset = rng.uniform(-3,3,size=(1,dim))
    x= rng.uniform(-1,1,size=(n,dim))*scale+offset
    y=rng.normal(size=(n,channels))
    if complex_data: y=y+1j*rng.normal(size=(n,channels))
    y=np.asarray(y,dtype=complex)
    q=rng.uniform(-1.15,1.15,size=(12,dim))*scale+offset
    q[:3]=x[:3]
    if mode==1:
        # Include points known to lie within the data's convex hull, as
        # well as queries outside it; Qhull and Bowyer-Watson may choose
        # different triangles but agree on barycentric interpolation.
        weights=rng.dirichlet(np.ones(n),size=8)
        q=np.vstack([q,weights@x,x])
    fill=complex(-2.75,0.65) if complex_data else complex(-2.75)
    return mode,x,y,q,fill,rescale,complex_data


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--seed',type=int,default=20261031)
    parser.add_argument('--cases',type=int,default=240)
    parser.add_argument('--mode',choices=('checked','optimized','both'),default='both')
    a=parser.parse_args()
    modes=[]
    if a.mode!='optimized': modes.append(('checked',['-O0','-g','-fcheck=all','-fbacktrace']))
    if a.mode!='checked': modes.append(('optimized',['-O3']))
    for title,flags in modes:
        rng=np.random.default_rng(a.seed)
        with tempfile.TemporaryDirectory() as td:
            exe=build(Path(td),flags)
            total=0
            maxerror=0.
            for k in range(a.cases):
                mode,x,y,q,fill,rescale,complex_data=case(rng,k)
                use_griddata=(k%7==0)
                if use_griddata:
                    ref=griddata(x,y,q,method=('nearest' if mode==0 else 'linear'),
                        rescale=rescale,fill_value=fill)
                elif mode==0:
                    ref=NearestNDInterpolator(x,y,rescale=rescale)(q)
                else:
                    ref=LinearNDInterpolator(x,y,rescale=rescale,fill_value=fill)(q)
                if ref.ndim==1: ref=ref[:,None]
                try:
                    got=run_case(exe,mode+2 if use_griddata else mode,x,y,q,fill,rescale,complex_data)
                except Exception as exc:
                    raise RuntimeError(f'failed {title} seed {a.seed} case {k}') from exc
                missing_ref=np.isnan(ref.real)
                missing_got=np.isnan(got.real)
                if not np.array_equal(missing_ref,missing_got):
                    raise AssertionError(f'{title} case {k}: mismatched NaN locations\nref={ref}\ngot={got}')
                err=np.max(np.abs((got-ref)[~missing_ref]),initial=0.)
                if err>5e-10:
                    raise AssertionError(f'{title} case {k}: {err=:.3g}\nref={ref}\ngot={got}')
                maxerror=max(maxerror,float(err))
                total+=ref.size
            print(f'{title}: PASS {a.cases} randomized cases, {total} scalar/complex entries, max abs={maxerror:.3g}',flush=True)

if __name__=='__main__': main()
