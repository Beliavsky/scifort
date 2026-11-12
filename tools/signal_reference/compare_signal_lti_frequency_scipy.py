#!/usr/bin/env python3
"""Seeded SciPy 1.17 comparisons for SciFort LTI frequency-response APIs.

Fortran test driver is in lti_frequency_driver.f90.  This script can use a
prebuilt executable with --driver; otherwise it compiles the minimal signal
module stack using the standalone test runner's MODULES list.
"""
from __future__ import annotations
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path
import warnings
import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[2]


def make_cases():
    rng = np.random.default_rng(20261017)
    cases = []
    for mode in range(1,7):
        for _ in range(120):
            order = int(rng.integers(1,6))
            degree = int(rng.integers(0,order+1))
            digital = mode in (2,4,6)
            dt = float(rng.uniform(0.1,1.5))
            if digital:
                roots = rng.uniform(-0.75,0.75,size=order)
                zeros = rng.uniform(-0.9,0.9,size=degree)
                w = np.sort(rng.uniform(0,2*np.pi,size=18))
            else:
                roots = -rng.uniform(0.4,3.0,size=order)
                zeros = -rng.uniform(0.3,4.0,size=degree)
                w = np.exp(rng.uniform(-2.4,2.0,size=18))
                w.sort()
            if mode in (3,4):
                b = zeros.astype(complex)
                a = roots.astype(complex)
                dt = float(rng.uniform(0.3,2.0))  # ZPK gain (probe convention)
            else:
                b = np.poly(zeros)*rng.uniform(0.2,2.0)
                a = np.poly(roots)
            cases.append((mode,np.atleast_1d(np.array(b,dtype=complex)),np.atleast_1d(np.array(a,dtype=complex)),w,dt))
    return cases


def reference(mode, b, a, w, dt):
    if mode in (3,4):
        system=(b,a,dt) if mode==3 else (b,a,dt,1.0)
    else:
        system=(b,a) if mode in (1,5) else (b,a,dt)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        if mode in (1,3): return signal.freqresp(system,w=w)
        if mode in (2,4): return signal.dfreqresp(system,w=w)
        if mode==5: return signal.bode(system,w=w)
        return signal.dbode(system,w=w)


def probe(driver, cases):
    buf=[]
    for mode,b,a,w,dt in cases:
        buf.append(f'{mode} {len(b)} {len(a)} {len(w)} {dt:.17g}')
        buf.extend(f'{v.real:.17g} {v.imag:.17g}' for v in b)
        buf.extend(f'{v.real:.17g} {v.imag:.17g}' for v in a)
        buf.extend(f'{v:.17g}' for v in w)
    process=subprocess.run([str(driver)], input='\n'.join(buf)+'\n',text=True,
        capture_output=True,check=True,timeout=60)
    lines=iter(process.stdout.splitlines())
    results=[]
    for mode,b,a,w,dt in cases:
        status,count=map(int,next(lines).split())
        values=np.array([[float(t.replace('D','E')) for t in next(lines).split()]
                         for _ in range(count)],dtype=float)
        results.append((status,values))
    return results


def check(driver):
    cases=make_cases()
    results=probe(driver,cases)
    max_error=0.0
    errors=[]
    for index,((mode,b,a,w,dt),(status,data)) in enumerate(zip(cases,results)):
        ref=reference(mode,b,a,w,dt)
        if status!=0 or data.shape!=(len(w),3):
            errors.append((index,'status',status,data.shape));continue
        expected_w=np.real(ref[0])
        expected2=np.asarray(ref[1])
        if mode<5:
            actual2=data[:,1]+1j*data[:,2]
            err=np.max(np.abs(actual2-expected2)/(1+np.abs(expected2)))
        else:
            expected3=np.asarray(ref[2])
            err=max(np.max(np.abs(data[:,1]-expected2)),
                    np.max(np.abs(data[:,2]-expected3)))*1.0
        err=max(err,np.max(np.abs(data[:,0]-expected_w)))
        max_error=max(max_error,float(err))
        if err > (3e-8 if mode>=5 else 1e-8):
            errors.append((index,mode,err))
    if errors:
        for x in errors[:10]:print('FAIL',x)
        raise SystemExit(f'{len(errors)}/{len(cases)} mismatches')
    print(f'PASS: {len(cases)} seeded SciPy LTI comparisons; worst error {max_error:.3g}')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--driver',type=Path)
    args=parser.parse_args()
    if args.driver is not None:
        check(args.driver)
        return
    sys.path.insert(0,str(ROOT/'tools'/'signal_reference'))
    from run_signal_suite import MODULES
    with tempfile.TemporaryDirectory(prefix='scifort-lti-frequency-') as temp:
        work=Path(temp)
        objs=[]
        for name in MODULES:
            obj=work/(name+'.o')
            subprocess.run(['gfortran','-std=f2018','-O2','-I',str(work),'-J',str(work),
                '-c',str(ROOT/'src'/(name+'.f90')),'-o',str(obj)],check=True,cwd=work)
            objs.append(str(obj))
        exe=work/'lti_frequency_driver'
        subprocess.run(['gfortran','-std=f2018','-O2','-I',str(work),'-J',str(work),
            str(ROOT/'tools/signal_reference/lti_frequency_driver.f90'),*objs,
            '-o',str(exe)],check=True,cwd=work)
        check(exe)


if __name__=='__main__': main()
