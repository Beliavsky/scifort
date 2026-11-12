#!/usr/bin/env python3
"""SciPy comparisons for Fortran discrete LTI series and parallel MIMO systems."""
from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy import signal

from verify_dlti_objects import SOURCES, build, read_complex, recurrence, write_complex

ROOT = Path(__file__).resolve().parents[2]


def create_cases():
    rng=np.random.default_rng(24680301)
    cases=[]
    for mode in (0,1):
        for idx in range(55):
            n1=1+idx%3
            n2=1+(idx//3)%3
            ni=1+idx%3
            no=1+(idx//7)%3
            mid=1+(idx//5)%3
            nu1 = mid if mode==0 else ni
            ny2 = mid if mode==0 else no
            nu2=ni
            def make(n,nu,ny):
                a=np.triu(rng.normal(size=(n,n))*0.07)
                np.fill_diagonal(a,rng.uniform(-.6,.6,n))
                b=rng.normal(size=(n,nu))*.5
                c=rng.normal(size=(ny,n))*.5
                d=rng.normal(size=(ny,nu))*.2
                if idx%3==1:
                    a=a+0.02j*rng.normal(size=(n,n))
                    b=b+0.2j*rng.normal(size=(n,nu))
                    c=c+0.2j*rng.normal(size=(ny,n))
                    d=d+0.2j*rng.normal(size=(ny,nu))
                return a,b,c,d
            one=make(n1,nu1,no)
            two=make(n2,nu2,ny2)
            dt=(.1,.25,.5,1.0)[idx%4]
            sys1=signal.StateSpace(*one,dt=dt)
            sys2=signal.StateSpace(*two,dt=dt)
            composed=sys1*sys2 if mode==0 else sys1+sys2
            nt=15+idx%7
            u=rng.normal(size=(nt,nu2))
            if idx%3==1:
                u=u+0.2j*rng.normal(size=(nt,nu2))
            w=np.linspace(.03,2.8,14)
            cases.append(dict(mode=mode,one=one,two=two,dt=dt,nt=nt,nu1=nu1,
                              nu2=nu2,ny1=no,ny2=ny2,n1=n1,n2=n2,u=u,w=w,
                              sys=composed))
    return cases


def verify(executable,cases):
    data=[str(len(cases))]
    for c in cases:
        data.append(' '.join(map(str,[c['mode'],c['n1'],c['n2'],c['nu1'],
                                      c['ny1'],c['nu2'],c['ny2'],c['nt'],
                                      len(c['w']),c['dt']])))
        data.append(' '.join(f'{w:.17e}' for w in c['w']))
        for values in (*c['one'],*c['two'],c['u']):
            write_complex(data,values)
    completed=subprocess.run([str(executable)],input='\n'.join(data)+'\n',
                             capture_output=True,text=True,check=True)
    lines=completed.stdout.splitlines()
    pos=0
    num_values=0
    max_scaled=0.0
    max_abs=0.0
    for index,c in enumerate(cases,1):
        expected_status=[index,0,0,0,0,0]
        status=list(map(int,lines[pos].split()))
        pos+=1
        if status!=expected_status:
            raise AssertionError(f'{index}: statuses {status}')
        n,nu,ny=c['nt'],c['nu2'],c['ny1']
        w=c['w']
        sys=c['sys']
        aa,bb,cc,dd=(np.asarray(v,dtype=complex) for v in (sys.A,sys.B,sys.C,sys.D))
        eye=np.eye(len(aa),dtype=complex)
        fr_ref=np.stack([cc@np.linalg.solve(np.exp(1j*wi)*eye-aa,bb)+dd for wi in w],axis=2)
        y_ref=recurrence(sys,c['u'])
        step_ref=np.empty((n,ny,nu),dtype=complex)
        impulse_ref=np.empty_like(step_ref)
        for input_index in range(nu):
            ex=np.zeros((n,nu),dtype=complex)
            ex[:,input_index]=1
            step_ref[:,:,input_index]=recurrence(sys,ex)
            ex[:]=0
            ex[0,input_index]=1
            impulse_ref[:,:,input_index]=recurrence(sys,ex)
        fr=np.empty((ny,nu,len(w)),dtype=complex)
        y=np.empty((n,ny),dtype=complex)
        sr=np.empty((n,ny,nu),dtype=complex)
        ir=np.empty_like(sr)
        for k in range(len(w)):
            for j in range(nu):
                for i in range(ny):
                    fr[i,j,k],pos=read_complex(lines,pos)
        for k in range(n):
            for i in range(ny):
                y[k,i],pos=read_complex(lines,pos)
        for k in range(n):
            for j in range(nu):
                for i in range(ny):
                    sr[k,i,j],pos=read_complex(lines,pos)
                    ir[k,i,j],pos=read_complex(lines,pos)
        for label,actual,expected in (('frequency',fr,fr_ref),('output',y,y_ref),
                                      ('step',sr,step_ref),('impulse',ir,impulse_ref)):
            delta=np.max(np.abs(actual-expected))
            scaled=np.max(np.abs(actual-expected)/(1+np.abs(expected)))
            max_scaled=max(max_scaled,float(scaled))
            max_abs=max(max_abs,float(delta))
            num_values+=actual.size
            if scaled>2e-10:
                raise AssertionError(f'{index} {label}: max_scaled={scaled}, max_abs={delta}')
    if pos!=len(lines):
        raise AssertionError(f'{len(lines)-pos} trailing lines')
    return len(cases),num_values,max_abs,max_scaled


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--optimized',action='store_true')
    args=parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='scifort-dlti-composition-') as tmp:
        work=Path(tmp)
        build(work,args.optimized)
        exe=work/'dlti_composition_driver'
        flags=['-std=f2018','-O3' if args.optimized else '-O0']
        if not args.optimized:
            flags.extend(['-fcheck=all','-g','-fbacktrace'])
        objects=[str(work/f'{stem}.o') for stem in SOURCES]
        subprocess.run(['gfortran',*flags,'-I',str(work),
                        str(ROOT/'tools/signal_reference/dlti_composition_driver.f90'),
                        *objects,'-o',str(exe)],check=True,cwd=work)
        count,values,diff,scaled=verify(exe,create_cases())
        print(f'PASS: {count} SciPy MIMO compositions, {values} values, '
              f'max_abs={diff:.4g}, max_scaled={scaled:.4g} '
              f'({"-O3" if args.optimized else "-O0"})')


if __name__=='__main__':
    main()
