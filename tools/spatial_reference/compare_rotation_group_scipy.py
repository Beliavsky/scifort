#!/usr/bin/env python3
"""Validate SO(3) symmetry groups and reduction against SciPy 1.17+.

Group element order may differ from SciPy; compare rotation sets using
q/-q equivalence. Reduction indices are one-based and refer to SciFort's
own ordering, so check their reconstruction rather than raw index equality.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parents[2]


def build(tmp: Path, mode: str) -> Path:
    flags = ['-O0','-g','-fcheck=all','-fbacktrace'] if mode == 'checked' else ['-O3']
    fc = ['gfortran','-std=f2018','-ffree-line-length-none',*flags]
    objects = []
    for name in ('kinds','spatial_transform','spatial_rotation_stack','spatial_rotation_group'):
        obj = tmp/(name+'.o')
        subprocess.run([*fc,'-I',str(tmp),'-J',str(tmp),'-c',str(ROOT/'src'/f'{name}.f90'),'-o',str(obj)],
                       cwd=tmp,check=True,capture_output=True)
        objects.append(str(obj))
    exe = tmp/'rotation_group_driver'
    subprocess.run([*fc,'-I',str(tmp),str(ROOT/'tools/spatial_reference/rotation_group_driver.f90'),
                    *objects,'-o',str(exe)],cwd=tmp,check=True,capture_output=True)
    return exe


def qerror(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.minimum(np.max(np.abs(a-b),axis=-1), np.max(np.abs(a+b),axis=-1))


def cases(seed: int, count: int) -> list[tuple]:
    rng=np.random.default_rng(seed)
    allgroups=[f'{p}{n}' for p in 'CD' for n in (1,2,3,4,5,7,11,17)] + ['T','O','I']
    axis=['X','Y','Z','x','y','z']
    generated=[]
    for i in range(count):
        gname=allgroups[i%len(allgroups)]
        ga=axis[(i//len(allgroups))%len(axis)]
        lname=allgroups[int(rng.integers(len(allgroups)))] if i%7 else 'NONE'
        rname=allgroups[int(rng.integers(len(allgroups)))] if i%9 else 'NONE'
        if i%21==0:
            lname, rname='NONE','NONE'
        # Large two-sided finite groups periodically exercise O(60^2).
        if i%37==0:
            lname,rname='I','O'
        if i%41==0:
            lname,rname='O','I'
        la=axis[int(rng.integers(len(axis)))]
        ra=axis[int(rng.integers(len(axis)))]
        n=int(rng.integers(1,9))
        rot=Rotation.random(n,rng=rng)
        q=rot.as_quat()
        if i%5==0:
            q[0]=Rotation.from_rotvec([0.,0.,1e-13]).as_quat()
        if i%11==0:
            q[-1]=Rotation.from_rotvec([np.pi-1e-9,0.,0.]).as_quat()
        generated.append((gname,ga,lname,la,rname,ra,q))
    return generated


def run(exe: Path, seed: int, count: int) -> tuple[float,int,int]:
    data=cases(seed,count)
    lines=[str(len(data))]
    for gname,ga,lname,la,rname,ra,q in data:
        lines.append(f'{gname} {ga} {lname} {la} {rname} {ra} {len(q)}')
        lines.extend(' '.join(f'{x:.17e}' for x in row) for row in q)
    process=subprocess.run([str(exe)],input='\n'.join(lines)+'\n',text=True,capture_output=True,check=True)
    records=iter(process.stdout.splitlines())
    worst=0.0
    compared=0
    for i,(gname,ga,lname,la,rname,ra,q) in enumerate(data):
        n,k=map(int,next(records).split())
        assert k==len(q),(i,k,len(q))
        actual=np.array([np.fromstring(next(records),sep=' ') for _ in range(n)])
        expected=Rotation.create_group(gname,axis=ga).as_quat()
        assert actual.shape==expected.shape,(i,actual.shape,expected.shape)
        # Directed Hausdorff distance under quaternion sign equivalence.
        dot=np.abs(actual@expected.T)
        mismatch=max(float(1-dot.max(axis=0).min()),float(1-dot.max(axis=1).min()))
        worst=max(worst,mismatch)
        assert mismatch<2e-13,(i,gname,ga,mismatch)
        compared+=4*n
        left_n=int(next(records))
        left_q=np.array([np.fromstring(next(records),sep=' ') for _ in range(left_n)])
        right_n=int(next(records))
        right_q=np.array([np.fromstring(next(records),sep=' ') for _ in range(right_n)])
        out=np.array([np.fromstring(next(records),sep=' ') for _ in range(k)])
        assert out.shape==(k,6),(i,out.shape)
        il=out[:,0].astype(int); ir=out[:,1].astype(int)
        result=Rotation.from_quat(out[:,2:6]); p=Rotation.from_quat(q)
        left=None if lname=='NONE' else Rotation.create_group(lname,axis=la)
        right=None if rname=='NONE' else Rotation.create_group(rname,axis=ra)
        ref=p.reduce(left,right)
        # Multiple symmetry pairs can yield different optimum rotations with
        # exactly the same minimal angle. Compare the invariant magnitude,
        # and independently reconstruct the returned representatives.
        err=float(np.max(np.abs(result.magnitude()-ref.magnitude())))
        worst=max(worst,err)
        assert err<2e-12,(i,lname,rname,err)
        if left is None:
            assert left_n==0
        else:
            assert left_n==len(left)
            ld=np.abs(left_q@left.as_quat().T)
            assert max(1-ld.max(0).min(),1-ld.max(1).min())<2e-13
        if right is None:
            assert right_n==0
        else:
            assert right_n==len(right)
            rd=np.abs(right_q@right.as_quat().T)
            assert max(1-rd.max(0).min(),1-rd.max(1).min())<2e-13
        rebuilt=p
        if left is not None:
            rebuilt=Rotation.from_quat(left_q[il-1])*rebuilt
        if right is not None:
            rebuilt=rebuilt*Rotation.from_quat(right_q[ir-1])
        reconstruction_err=float(np.max(qerror(result.as_quat(),rebuilt.as_quat())))
        assert reconstruction_err<2e-12,(i,lname,rname,'reconstruction',reconstruction_err)
        # Verify 1-based indices of each factor, with absent groups indexed 0.
        assert np.all(il==0) if left is None else np.all((il>=1)&(il<=len(left)))
        assert np.all(ir==0) if right is None else np.all((ir>=1)&(ir<=len(right)))
        compared+=4*(k+left_n+right_n)
    assert next(records,None) is None,'Unexpected extra output'
    return worst,compared,len(data)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--mode',choices=('checked','optimized','both'),default='both')
    parser.add_argument('--count',type=int,default=180)
    args=parser.parse_args()
    for mode in (['checked','optimized'] if args.mode=='both' else [args.mode]):
        with tempfile.TemporaryDirectory(prefix='scifort-group-') as td:
            exe=build(Path(td),mode)
            total=0
            worst=0.0
            for seed in (20261102,20261103):
                discrepancy,nvalues,ncases=run(exe,seed,args.count)
                total+=nvalues
                worst=max(worst,discrepancy)
                print(f'{mode}: seed {seed}: {ncases} group/reduction cases, {nvalues} values: '
                      f'max discrepancy {discrepancy:.3e}',flush=True)
            print(f'{mode}: {2*args.count} comparisons, {total} values PASS; '
                  f'max discrepancy {worst:.3e}',flush=True)

if __name__=='__main__':
    main()
