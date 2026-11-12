#!/usr/bin/env python3
"""Compare SciFort SE(3) transforms and stacks with scipy.spatial.transform.

Two reproducible seeds; scalar conversions, screw powers, composition,
point application, and one-dimensional stack means/broadcast operations.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.spatial.transform import RigidTransform as Tf, Rotation as R

ROOT = Path(__file__).resolve().parents[2]


def build(tmp: Path, mode: str) -> Path:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    fc = ['gfortran', '-std=f2018', '-ffree-line-length-none', *flags]
    objs = []
    for name in ('kinds', 'spatial_transform', 'spatial_rigid_transform'):
        obj = tmp / (name + '.o')
        subprocess.run([*fc, '-I', str(tmp), '-J', str(tmp), '-c',
                        str(ROOT/'src'/f'{name}.f90'), '-o', str(obj)],
                       cwd=tmp,check=True,capture_output=True)
        objs.append(str(obj))
    exe = tmp/'rigid_transform_driver'
    subprocess.run([*fc, '-I', str(tmp), str(ROOT/'tools/spatial_reference/rigid_transform_driver.f90'),
                    *objs, '-o', str(exe)],cwd=tmp,check=True,capture_output=True)
    return exe


def generate(seed: int,count: int) -> list[tuple]:
    rng = np.random.default_rng(seed)
    result=[]
    for i in range(count):
        q1=R.random(rng=rng).as_quat()
        q2=R.random(rng=rng).as_quat()
        if i%7==0:
            q1=R.from_rotvec(rng.normal(size=3)*1e-10).as_quat()
        if i%11==0:
            axis=rng.normal(size=3); axis/=np.linalg.norm(axis)
            q1=R.from_rotvec((np.pi-1e-9)*axis).as_quat()
        if i%13==0:
            q1=R.from_rotvec([0,0,0]).as_quat()
        if i%17==0:
            q1=-q1
        if i%37==0:
            q1=R.from_rotvec(rng.normal(size=3)*1e-7).as_quat()
        t1=rng.normal(size=3)
        t2=rng.normal(size=3)
        point=rng.normal(size=3)
        exponent=float(rng.uniform(-2.5,3.0))
        if i%6==0:
            exponent=0.5
        elif i%6==1:
            exponent=-1.0
        elif i%6==2:
            exponent=0.0
        elif i%6==3:
            exponent=2.0
        if i%19==0:
            t1*=1e-7
            t2*=1e5
        if i%41==0:
            t1*=1e5
            t2*=1e-7
        result.append((q1,t1,q2,t2,point,exponent))
    return result


def read_vector(records, size: int, label: str) -> np.ndarray:
    row=next(records)
    arr=np.fromstring(row,sep=' ')
    assert arr.size==size,(label,arr.size,size,row)
    assert np.all(np.isfinite(arr)),label
    return arr


def run(exe: Path, seed: int, count: int) -> tuple[float,int]:
    cases=generate(seed,count)
    lines=[str(len(cases))]
    for case in cases:
        q1,t1,q2,t2,pt,power=case
        lines.append(' '.join(f'{float(x):.17e}' for x in np.r_[q1,t1,q2,t2,pt,power]))
    proc=subprocess.run([str(exe)],input='\n'.join(lines)+'\n',text=True,capture_output=True,check=True)
    records=iter(proc.stdout.splitlines())
    worst=0.0
    values=0
    for i,(q1,t1,q2,t2,pt,power) in enumerate(cases):
        a=Tf.from_components(t1,R.from_quat(q1))
        b=Tf.from_components(t2,R.from_quat(q2))
        stack=Tf.from_components(np.stack([t1,t2]),R.from_quat([q1,q2]))
        expected=[
            ('matrix',a.as_matrix().reshape(-1,order='F')),
            ('log',a.as_exp_coords()),
            ('dual',a.as_dual_quat()),
            ('compose',(a*b).as_matrix().reshape(-1,order='F')),
            ('inverse',a.inv().as_matrix().reshape(-1,order='F')),
            ('power',Tf.from_exp_coords(a.as_exp_coords()*power).as_matrix().reshape(-1,order='F')),
            ('apply',np.r_[a.apply(pt),a.apply(pt,inverse=True)]),
            ('exp_roundtrip',a.as_matrix().reshape(-1,order='F')),
            ('dual_roundtrip',a.as_matrix().reshape(-1,order='F')),
            ('matrix_roundtrip',a.as_matrix().reshape(-1,order='F')),
            ('mean',stack.mean(weights=[1,2]).as_matrix().reshape(-1,order='F')),
            ('stack_apply',stack.apply(pt).ravel(order='F')),
            ('stack_compose',np.broadcast_to(np.eye(4),(2,4,4)).ravel(order='F')),
            ('stack_exp',stack.as_matrix().ravel(order='F'))
        ]
        for label,ref in expected:
            got=read_vector(records,ref.size,label)
            if label=='dual':
                err=min(np.max(np.abs(got-ref)),np.max(np.abs(got+ref)))
            elif label=='log' and np.linalg.norm(ref[:3])>=np.pi-1e-8:
                # At exactly pi the logarithm can switch equivalent axis sign.
                err=max(0.0,float(np.max(np.abs(got[3:]-ref[3:]))))
            else:
                err=float(np.max(np.abs(got-ref)))
            norm=max(1.0,float(np.max(np.abs(ref))))
            scaled=err/norm
            worst=max(worst,scaled)
            tol = 3e-9 if label=='log' else 6e-10
            assert scaled<tol,(seed,i,label,scaled,got,ref)
            values+=ref.size
    assert next(records,None) is None,'unexpected extra output'
    return worst,values


def run_twist_stress(tmp: Path, mode: str, seed: int=20261106, count: int=600) -> tuple[float,int]:
    """Compare full SE(3) exponential coordinates, including nonprincipal angles."""
    fc = ['gfortran','-std=f2018','-ffree-line-length-none',
          *(['-O0','-g','-fcheck=all'] if mode=='checked' else ['-O3'])]
    objs=[]
    for name in ('kinds','spatial_transform','spatial_rigid_transform'):
        obj=tmp/f'{name}.o'
        subprocess.run([*fc,'-I',str(tmp),'-J',str(tmp),'-c',
                        str(ROOT/'src'/f'{name}.f90'),'-o',str(obj)],
                       check=True,capture_output=True,cwd=tmp)
        objs.append(str(obj))
    exe=tmp/'rigid_twist_driver'
    subprocess.run([*fc,'-I',str(tmp),str(ROOT/'tools/spatial_reference/rigid_twist_driver.f90'),
                    *objs,'-o',str(exe)],check=True,capture_output=True,cwd=tmp)
    rng=np.random.default_rng(seed)
    data=rng.normal(size=(count,6))
    angles=[0.0,1e-12,1e-8,1e-5,1e-3,0.1,np.pi-1e-10,
            np.pi,np.pi+1e-10,2*np.pi-1e-6,2*np.pi,2*np.pi+1e-6,
            3*np.pi,5*np.pi,10.0]
    for i in range(count):
        if i%2==0:
            theta=angles[(i//2)%len(angles)]
            axis=data[i,:3]; axis/=np.linalg.norm(axis)
            data[i,:3]=theta*axis
        if i%11==0:
            data[i,3:]*=1000.0
    lines=[str(count)]+[' '.join(f'{x:.17e}' for x in row) for row in data]
    p=subprocess.run([str(exe)],input='\n'.join(lines)+'\n',text=True,
                     capture_output=True,check=True)
    received=np.array([[float(x) for x in line.split()] for line in p.stdout.splitlines()])
    assert received.shape==(count,16),received.shape
    actual=received.reshape((count,4,4),order='F')
    expected=Tf.from_exp_coords(data).as_matrix()
    norm=np.maximum(1.0,np.max(np.abs(expected),axis=(1,2)))
    err=np.max(np.abs(actual-expected),axis=(1,2))/norm
    worst=float(err.max())
    assert np.all(err<3e-10),(worst,np.argmax(err),data[np.argmax(err)])
    return worst, count*16

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--mode',choices=('checked','optimized','both'),default='both')
    parser.add_argument('--count',type=int,default=260)
    args=parser.parse_args()
    for mode in (['checked','optimized'] if args.mode=='both' else [args.mode]):
        with tempfile.TemporaryDirectory(prefix='scifort-rigid-') as td:
            exe=build(Path(td),mode)
            total=0
            worst=0.0
            for seed in (20261102,20261105):
                disc,nvalues=run(exe,seed,args.count)
                total+=nvalues
                worst=max(worst,disc)
                print(f'{mode}: seed {seed}: {args.count} SE(3) cases, {nvalues} values, '
                      f'max scale-relative error {disc:.3e}',flush=True)
            with tempfile.TemporaryDirectory(prefix='scifort-rigid-twist-') as twd:
                stress_error, nstress=run_twist_stress(Path(twd),mode)
            total += nstress
            worst = max(worst,stress_error)
            print(f'{mode}: {2*args.count} SE(3) cases + 600 twist cases, '
                  f'{total} values PASS; max error {worst:.3e}',flush=True)
if __name__=='__main__':
    main()
