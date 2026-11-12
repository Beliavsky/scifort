#!/usr/bin/env python3
"""Cross-check SciFort's one-based cluster leaders with SciPy 1.17."""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.cluster import hierarchy as sch

ROOT = Path(__file__).resolve().parents[2]


def compile_native(folder, mode):
    flags = ['-O0','-g','-fcheck=all','-fbacktrace'] if mode=='checked' else ['-O3']
    base=['gfortran','-std=f2018','-ffree-line-length-none',*flags,
          '-J',str(folder),'-I',str(folder)]
    objects=[]
    for module in ('kinds','spatial_distance','cluster_hierarchy'):
        obj=folder/(module+'.o')
        subprocess.run([*base,'-c',str(ROOT/'src'/(module+'.f90')),
                        '-o',str(obj)],check=True,cwd=folder,capture_output=True)
        objects.append(str(obj))
    exe=folder/'leaders_driver'
    subprocess.run([*base,str(ROOT/'tools/cluster_reference'/'cluster_hierarchy_leaders_driver.f90'),
                    *objects,'-o',str(exe)],check=True,cwd=folder,capture_output=True)
    return exe


def run_case(exe,z,labels,expected_valid):
    n=len(labels)
    z1=z.copy(); z1[:,:2]+=1
    data=[str(n), *[' '.join(f'{v:.17e}' for v in row) for row in z1],
          ' '.join(str(int(v)) for v in labels)]
    proc=subprocess.run([str(exe)],input='\n'.join(data)+'\n',text=True,
                        capture_output=True,check=True)
    lines=proc.stdout.splitlines()
    status,nfound=map(int,lines[0].split())
    if not expected_valid:
        if status != 1 or nfound != 0:
            raise AssertionError('invalid assignment not detected')
        return
    if status != 0:
        raise AssertionError(f'valid assignment rejected: {proc.stdout}')
    pairs=[tuple(map(int,ln.split())) for ln in lines[1:]]
    if len(pairs)!=nfound:
        raise AssertionError('wrong leader count')
    scipy_nodes,scipy_labels=sch.leaders(z,labels.astype(np.int32))
    expected={(int(node)+1,int(label)) for node,label in zip(scipy_nodes,scipy_labels)}
    if set(pairs) != expected:
        raise AssertionError(f'leader pairs differ: {set(pairs)} vs {expected}')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--mode',choices=('checked','optimized','both'),default='both')
    parser.add_argument('--cases',type=int,default=180)
    args=parser.parse_args()
    for mode in (('checked','optimized') if args.mode=='both' else (args.mode,)):
        rng=np.random.default_rng(20261103)
        with tempfile.TemporaryDirectory(prefix='scifort-leaders-') as tmp:
            exe=compile_native(Path(tmp),mode)
            for i in range(args.cases):
                n=int(rng.integers(2,35))
                x=rng.normal(size=(n,int(rng.integers(1,7))))
                method=('single','complete','average','ward','centroid','median')[i%6]
                z=sch.linkage(x,method=method)
                target=int(rng.integers(1,n+1))
                # Build a row-order forest after n-target merges. SciPy's
                # cut_tree can create non-subtree partitions if the linkage
                # has inversions (notably centroid/median methods).
                active = {leaf: [leaf] for leaf in range(n)}
                for row_no in range(n-target):
                    a,b = map(int,z[row_no,:2])
                    active[n+row_no] = active.pop(a)+active.pop(b)
                T = np.zeros(n,dtype=np.int32)
                for label,member in enumerate(active.values()):
                    T[member] = label
                ids,inv=np.unique(T,return_inverse=True)
                relabel=rng.choice(np.arange(-1000,1000,7),size=len(ids),replace=False)
                T=relabel[inv].astype(np.int32)
                run_case(exe,z,T,True)
                if n>=4 and len(ids)>=2:
                    # Construct a deliberately impossible flat partition:
                    # merge two disjoint singleton labels while preserving
                    # the other labels (if the resulting set is nonsubtree).
                    invalid=np.arange(n,dtype=np.int32)
                    invalid[-1]=invalid[0]
                    try:
                        sch.leaders(z,invalid)
                    except ValueError:
                        run_case(exe,z,invalid,False)
        print(f'{mode}: {args.cases} randomized leader mappings PASS',flush=True)

if __name__ == '__main__': main()
