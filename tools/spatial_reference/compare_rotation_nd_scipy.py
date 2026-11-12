#!/usr/bin/env python3
"""SciPy reference tests for rank-four rotation/rigid-transform broadcasting.

Inputs use Fortran-major batch ordering; outputs are flattened in the same
order. Random rotations and translations are generated independently of SciFort.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.spatial.transform import Rotation, RigidTransform

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ['kinds','spatial_transform','spatial_rigid_transform','spatial_transform_nd']

def build(mode: str, tmp: Path) -> Path:
    flags = ['-O0','-g','-fcheck=all','-fbacktrace'] if mode == 'checked' else ['-O3']
    cmd = ['gfortran','-std=f2018','-ffree-line-length-none',*flags]
    objs=[]
    for name in SOURCES:
        obj=tmp/(name+'.o')
        subprocess.run([*cmd,'-I',str(tmp),'-J',str(tmp),'-c',str(ROOT/'src'/f'{name}.f90'),'-o',str(obj)],
                       cwd=tmp,check=True,capture_output=True,text=True)
        objs.append(str(obj))
    exe=tmp/'rotation_nd_driver'
    subprocess.run([*cmd,'-I',str(tmp),str(ROOT/'tools/spatial_reference/rotation_nd_driver.f90'),*objs,'-o',str(exe)],
                   cwd=tmp,check=True,capture_output=True,text=True)
    return exe

def flat_line(arr: np.ndarray) -> str:
    return ' '.join(f'{v:.17e}' for v in arr.ravel(order='F'))

def trial_cases(seed: int, ncases: int):
    rng=np.random.default_rng(seed)
    cases=[]
    chunks=[str(ncases)]
    for _ in range(ncases):
        # Every axis may be shared, absent (length 1), or expanded.
        # Require compatible nontrivial cases, including cross broadcasting.
        shape_a=[];shape_b=[]
        for k in range(4):
            common=int(rng.integers(1,5))
            kind=int(rng.integers(0,4))
            if kind==0:
                shape_a.append(1);shape_b.append(common)
            elif kind==1:
                shape_a.append(common);shape_b.append(1)
            else:
                shape_a.append(common);shape_b.append(common)
        sa=tuple(shape_a);sb=tuple(shape_b)
        qa=Rotation.random(int(np.prod(sa)),random_state=rng).as_quat().reshape((*sa,4))
        qb=Rotation.random(int(np.prod(sb)),random_state=rng).as_quat().reshape((*sb,4))
        ta=rng.normal(size=(*sa,3))
        tb=rng.normal(size=(*sb,3))
        points=rng.normal(size=(*sb,3))
        for v in (sa,sb): chunks.append(' '.join(map(str,v)))
        for arr in (qa,qb,ta,tb,points): chunks.append(flat_line(arr))
        ra=Rotation.from_quat(qa); rb=Rotation.from_quat(qb)
        rc=ra*rb
        aa=RigidTransform.from_components(ta,ra)
        bb=RigidTransform.from_components(tb,rb)
        cc=aa*bb
        expected_shape=np.broadcast_shapes(sa,sb)
        expected = (rc.as_matrix(), cc.as_matrix(), rc.apply(points), cc.apply(points))
        cases.append((expected_shape,expected))
    return '\n'.join(chunks)+'\n', cases

def rows_from_expected(exp):
    shape,(rm,tm,rp,tp)=exp
    count=int(np.prod(shape))
    # Each reference output has batch axes in Fortran-major order.
    # Matrix components are arranged in column-major order within each row.
    matrices=np.reshape(rm,(-1,3,3),order='F')
    transforms=np.reshape(tm,(-1,4,4),order='F')
    rpts=np.reshape(rp,(-1,3),order='F')
    tpts=np.reshape(tp,(-1,3),order='F')
    out=np.empty((count,31))
    for i in range(count):
        out[i,:9]=matrices[i].ravel(order='F')
        out[i,9:25]=transforms[i].ravel(order='F')
        out[i,25:28]=rpts[i]
        out[i,28:31]=tpts[i]
    return out

def compare(mode: str, seed: int, ncases: int):
    with tempfile.TemporaryDirectory(prefix='scifort-rotation-nd-') as directory:
        exe=build(mode,Path(directory))
        input_, cases = trial_cases(seed,ncases)
        proc=subprocess.run([str(exe)],input=input_,text=True,capture_output=True,check=True)
        lines=iter(proc.stdout.splitlines())
        count=0; worst=0.0
        for expected_shape, expected in cases:
            header=next(lines).split()
            if header[0]!='CASE' or tuple(map(int,header[1:]))!=expected_shape:
                raise AssertionError(f'shape mismatch: {header} vs {expected_shape}')
            ref=rows_from_expected((expected_shape,expected))
            actual=np.array([[float(x) for x in next(lines).split()] for i in range(len(ref))])
            if actual.shape!=ref.shape:
                raise AssertionError((actual.shape,ref.shape))
            diff=np.max(np.abs(actual-ref))
            scale=max(1.0,np.max(np.abs(ref)))
            worst=max(worst,diff/scale)
            if diff/scale>5e-11:
                i,j=np.unravel_index(np.argmax(np.abs(actual-ref)),actual.shape)
                raise AssertionError(f'case={count} error={diff/scale} i={i} j={j} actual={actual[i,j]} reference={ref[i,j]}')
            count+=1
        if next(lines,None) is not None:
            raise AssertionError('unexpected trailing driver output')
        print(f'{mode} seed={seed}: {count} cases {sum(np.prod(c[0])*31 for c in cases)} values; max scale-relative error {worst:.3e}')

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--mode',choices=['checked','optimized','both'],default='both')
    p.add_argument('--cases',type=int,default=120)
    p.add_argument('--seed',type=int,default=20261103)
    a=p.parse_args()
    for mode in (['checked','optimized'] if a.mode=='both' else [a.mode]):
        compare(mode,a.seed,a.cases)
