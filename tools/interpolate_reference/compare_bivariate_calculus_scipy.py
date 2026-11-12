#!/usr/bin/env python3
"""Cross-validate bivariate TCK differentiation and Boehm knot refinement.

SciPy FITPACK provides the derivative spline reference and the 1-D BSpline
class provides an independently compiled tensor-axis knot-insertion oracle.
No external compiled library is linked into the Fortran implementation.
"""
from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy.interpolate import BSpline, RectBivariateSpline, bisplev

ROOT = Path(__file__).resolve().parents[2]


def as_text(v):
    return ' '.join(f'{float(n):.17g}' for n in np.ravel(v))


def compile_driver(work: Path, flags: list[str]):
    exe = work / 'bivariate_calculus_driver'
    command = ['gfortran', '-std=f2018', '-ffree-line-length-none', *flags,
               '-J' + str(work), '-I' + str(work),
               *(str(ROOT / 'src' / name) for name in (
                   'kinds.f90','interpolate_bspline.f90',
                   'interpolate_fitpack_bivariate.f90')),
               str(ROOT / 'tools/interpolate_reference/bivariate_calculus_driver.f90'),
               '-o',str(exe)]
    subprocess.run(command,check=True,text=True,capture_output=True)
    return exe


def read_tck(lines, offset):
    nx,ny,nc,kx,ky = (int(n) for n in lines[offset].split())
    tx = np.fromstring(lines[offset+1],sep=' ')
    ty = np.fromstring(lines[offset+2],sep=' ')
    c = np.fromstring(lines[offset+3],sep=' ')
    assert (len(tx),len(ty),len(c)) == (nx,ny,nc)
    return (tx,ty,c,kx,ky),offset+4


def scipy_insert(tck,axis,knot,m):
    tx,ty,c,kx,ky=tck
    nx=len(tx)-kx-1
    ny=len(ty)-ky-1
    mat=c.reshape((nx,ny),order='C')
    if m == 0:
        return tck
    if axis==1:
        bs=BSpline(tx,mat,kx,axis=0).insert_knot(knot,m=m)
        return (bs.t,ty,bs.c.ravel(order='C'),kx,ky)
    bs=BSpline(ty,mat.T,ky,axis=0).insert_knot(knot,m=m)
    return (tx,bs.t,bs.c.T.ravel(order='C'),kx,ky)


def scipy_derivative(tck, dx, dy):
    """Derive tensor coefficients through SciPy's independently tested BSpline.

    SciPy 1.17's FITPACK pardtc signature mistakenly tests dx < ky in
    some asymmetric-degree cases; the 1D BSpline derivative avoids that bug.
    """
    tx,ty,c,kx,ky=tck
    mat=c.reshape((len(tx)-kx-1,len(ty)-ky-1),order='C')
    bsx=BSpline(tx,mat,kx,axis=0).derivative(dx)
    cx=bsx.c[:len(bsx.t)-(kx-dx)-1,:]
    bsy=BSpline(ty,cx.T,ky,axis=0).derivative(dy)
    cy=bsy.c[:len(bsy.t)-(ky-dy)-1,:]
    return (bsx.t,bsy.t,cy.T.ravel(order='C'),kx-dx,ky-dy)


def scipy_tensor_integral(tck, bounds):
    """Independent exact integral by successive SciPy BSpline integrations."""
    tx,ty,c,kx,ky=tck
    nx=len(tx)-kx-1
    ny=len(ty)-ky-1
    xa,xb,ya,yb=bounds
    lo_x=max(min(xa,xb),tx[kx])
    hi_x=min(max(xa,xb),tx[nx])
    lo_y=max(min(ya,yb),ty[ky])
    hi_y=min(max(ya,yb),ty[ny])
    if lo_x>=hi_x or lo_y>=hi_y:
        return 0.0
    cmat=c.reshape((nx,ny),order='C')
    cy=BSpline(ty,cmat.T,ky,axis=0).integrate(lo_y,hi_y)
    result=BSpline(tx,cy,kx).integrate(lo_x,hi_x)
    return float(result)*(-1 if xa>xb else 1)*(-1 if ya>yb else 1)


def relative_error(a,b):
    a=np.asarray(a)
    b=np.asarray(b)
    assert a.shape==b.shape, (a.shape,b.shape)
    return float(np.max(np.abs(a-b)/np.maximum(1.0,np.abs(b)),initial=0))


def compare_tck(actual, expected, stage):
    tx,ty,c,kx,ky=actual
    ex,ey,ec,ekx,eky=expected
    assert (kx,ky)==(ekx,eky), f'{stage}: degrees {(kx,ky)} != {(ekx,eky)}'
    return np.array([relative_error(tx,ex),relative_error(ty,ey),relative_error(c,ec)])


def compare_case(exe, rng, case):
    kx = 1+case%5
    ky = 1+(case*3//5 + 2*case)%5
    nx = kx+1+(case*13)%18
    ny = ky+1+(case*7)%17
    x = np.cumsum(rng.uniform(.25,1.2,nx))
    y = np.cumsum(rng.uniform(.21,1.3,ny))
    x = (x-x[0])/(x[-1]-x[0])*3.0-1.4
    y = (y-y[0])/(y[-1]-y[0])*2.3-.8
    z = (np.sin(1.2*x[:,None]+1.4*y[None,:])+
         .15*x[:,None]*y[None,:]**2 +
         .07*rng.standard_normal((nx,ny)))
    dx = int(rng.integers(kx))
    dy = int(rng.integers(ky))
    axis=1+case%2
    original=RectBivariateSpline(x,y,z,kx=kx,ky=ky,s=0)
    tx,ty=original.get_knots()
    original_tck=(tx,ty,original.get_coeffs(),kx,ky)
    knots=tx if axis==1 else ty
    degree=kx if axis==1 else ky
    base_x=x if axis==1 else y
    # Include insertion at a pre-existing interior knot when possible.
    available = knots[degree+1:len(knots)-degree-1]
    if case%3==0 and len(available):
        knot=float(available[int(rng.integers(len(available)))])
    else:
        knot=float(base_x[0]+(base_x[-1]-base_x[0])*rng.uniform(.1,.9))
    multiplicity=int(np.count_nonzero(knots==knot))
    # Keep the refined derivative representable as a B-spline. Inserting a
    # discontinuity of multiplicity degree+1 makes differentiation singular.
    order_axis=dx if axis==1 else dy
    max_reps=degree+1-order_axis-multiplicity
    m=int(case%min(max_reps+1,4))
    # Ensure both actual insertions and m=0 cases are represented.
    if case%5==0 and max_reps>0: m=max_reps if case%2 else 0
    nq=20
    qx=rng.uniform(x[0]-.3,x[-1]+.3,nq)
    qy=rng.uniform(y[0]-.3,y[-1]+.3,nq)
    bounds=[*rng.uniform(x[0]-.5,x[-1]+.5,2),
            *rng.uniform(y[0]-.4,y[-1]+.4,2)]
    lines=[f'{nx} {ny} {kx} {ky} {nq} {dx} {dy} {axis} {m} {knot:.17g}',
           as_text(x),as_text(y),as_text(z.ravel(order='F')),
           as_text(qx),as_text(qy),as_text(bounds)]
    proc=subprocess.run([str(exe)],input='\n'.join(lines)+'\n',capture_output=True,text=True)
    if proc.returncode:
        raise RuntimeError(f'case {case}: {proc.stderr[-1500:]}')
    output=proc.stdout.splitlines()
    statuses=tuple(int(n) for n in output[0].split())
    if statuses!=(0,0,0,0):
        raise AssertionError(f'case {case}: statuses={statuses} input={kx,ky,dx,dy,axis,m,knot}')
    tcks=[]
    offset=1
    for _ in range(4):
        tck,offset=read_tck(output,offset)
        tcks.append(tck)
    values=np.array([np.fromstring(line,sep=' ') for line in output[offset:offset+nq]])
    integral=np.fromstring(output[offset+nq],sep=' ')
    assert values.shape==(nq,6), values.shape
    assert integral.shape==(4,),integral.shape
    ref_der=scipy_derivative(original_tck,dx,dy)
    # Also cross-check against SciPy's native bivariate derivative when its
    # FITPACK wrapper accepts this asymmetric-degree combination.
    try:
        native_derived=original.partial_derivative(dx,dy)
    except Exception as exc:
        if 'pardtc:nuy=' not in str(exc):
            raise
        # SciPy 1.17.0 mistakenly checks nux against ky in some cases.
        native_derived=None
    if native_derived is not None:
        native_tck=(*native_derived.get_knots(),native_derived.get_coeffs(),*native_derived.degrees)
        assert np.all(compare_tck(native_tck,ref_der,'SciPy derivative consistency')<2e-8)
    inserted=scipy_insert(original_tck,axis,knot,m)
    errors=np.zeros(8)
    for i,(ours,ref) in enumerate(zip(tcks[:3],(original_tck,ref_der,inserted))):
        errs=compare_tck(ours,ref,stage=['fit','derivative','insertion'][i])
        errors[:3]=np.maximum(errors[:3],errs)
    # Check the derivative of the refined surface also via SciPy's TCK
    # evaluator, without needing to use SciPy's private FITPACK classes.
    cx=np.clip(qx,x[0],x[-1]); cy=np.clip(qy,y[0],y[-1])
    ref_undiff=original.ev(cx,cy)
    ref_dxdy=original.ev(cx,cy,dx=dx,dy=dy)
    ref_insert=np.array([bisplev(a,b,inserted) for a,b in zip(cx,cy)])
    ref_insert_dxdy=np.array([bisplev(a,b,inserted,dx=dx,dy=dy) for a,b in zip(cx,cy)])
    targets=[ref_undiff,ref_dxdy,ref_insert,ref_insert_dxdy,ref_dxdy,ref_insert_dxdy]
    for j,target in enumerate(targets):
        errors[3]=max(errors[3],relative_error(values[:,j],target))
    # The derivative surface generated after insertion should also preserve
    # the original FITPACK derivative coefficients at the numerical level.
    errors[4]=relative_error(values[:,3],values[:,1])
    # SciPy FITPACK verifies exact integral of the original and unchanged
    # integral after knot insertion. The derivative satisfies the mixed-
    # partial fundamental theorem of calculus (when dx=dy=1).
    xlo,xhi,ylo,yhi=bounds
    ref_integral=original.integral(xlo,xhi,ylo,yhi)
    errors[5]=relative_error(integral[[0,2]],np.full(2,ref_integral))
    ref_der_integral=scipy_tensor_integral(ref_der,bounds)
    errors[6]=relative_error(integral[[1,3]],np.full(2,ref_der_integral))
    if dx==1 and dy==1:
        x1,x2=sorted((max(xlo,x[0]),min(xhi,x[-1])))
        y1,y2=sorted((max(ylo,y[0]),min(yhi,y[-1])))
        # Only compare the four corners when the clipping intervals overlap.
        if x1<x2 and y1<y2:
            corners=(original.ev(x2,y2)-original.ev(x1,y2)-
                     original.ev(x2,y1)+original.ev(x1,y1))
            # FITPACK's integral is non-oriented in some SciPy versions:
            # compare magnitudes, and only when given bounds were ordered.
            if xlo<xhi and ylo<yhi:
                errors[7]=relative_error(integral[1],corners)
    tolerance=np.array([5e-13,5e-13,2e-7,3e-6,3e-6,2e-7,3e-6,2e-6])
    if np.any(errors>tolerance):
        raise AssertionError(f'case={case} nx={nx} ny={ny} degrees={kx,ky} '
                             f'nu={dx,dy} insert={axis,knot,m} errors={errors}')
    return errors


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--cases',type=int,default=100)
    ap.add_argument('--seed',type=int,default=20261101)
    ap.add_argument('--mode',choices=('both','checked','optimized'),default='both')
    args=ap.parse_args()
    modes=[]
    if args.mode!='optimized': modes.append(('checked',['-O0','-g','-fcheck=all','-fbacktrace']))
    if args.mode!='checked': modes.append(('optimized',['-O3']))
    for mode,flags in modes:
        with tempfile.TemporaryDirectory() as d:
            exe=compile_driver(Path(d),flags)
            rng=np.random.default_rng(args.seed)
            totals=np.zeros(8)
            for n in range(args.cases):
                totals=np.maximum(totals,compare_case(exe,rng,n))
            print(f'{mode}: {args.cases}/{args.cases} SciPy mixed-derivative / '
                  f'knot-refinement cases PASS; max errors={totals}',flush=True)


if __name__=='__main__': main()
