#!/usr/bin/env python3
"""Compile the public scifort_signal dependency closure, then test mixed generics.

No FPM/third-party dependency is required. Explicit source dependency ordering
lets this operate on the signal-only part of the project.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
MODULE = re.compile(r'^\s*module\s+(?!procedure\b|function\b|subroutine\b)(\w+)',re.I|re.M)
USE = re.compile(r'^\s*use\s*(?:,\s*(?:intrinsic|non_intrinsic)\s*)?(?:::)?\s*(\w+)',re.I|re.M)

def closure() -> list[Path]:
    mapping={}
    reqs={}
    for p in (ROOT/'src').rglob('*.f90'):
        code=p.read_text(errors='replace')
        modules={n.lower() for n in MODULE.findall(code)}
        for n in modules:
            if n in mapping:
                raise ValueError(f'Duplicate Fortran module {n}: {p} and {mapping[n]}')
            mapping[n]=p
        reqs[p]={n.lower() for n in USE.findall(code)}
    output=[]
    visited=set()
    stack=set()
    def walk(p:Path):
        if p in visited: return
        if p in stack: raise ValueError('Cyclic dependency: '+str(p))
        stack.add(p)
        for name in sorted(reqs[p]):
            dep=mapping.get(name)
            if dep is not None and dep != p: walk(dep)
        stack.remove(p)
        visited.add(p)
        output.append(p)
    walk(ROOT/'src'/'signal.f90')
    return output


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler',default='gfortran')
    parser.add_argument('--build-dir',type=Path,default=ROOT/'build'/'signal_mixed_public')
    parser.add_argument('--count',type=int,default=999)
    opts=parser.parse_args()
    order=closure()
    opts.build_dir.mkdir(parents=True,exist_ok=True)
    flags=['-std=f2018','-ffree-line-length-none','-fimplicit-none','-Werror=implicit-interface',
           '-O0','-J.','-I.']
    print(f'Public module dependency closure has {len(order)} source files',flush=True)
    for i,p in enumerate(order,1):
        obj=p.stem+'.o'
        path=opts.build_dir/obj
        if path.exists(): continue
        if i>opts.count: break
        subprocess.run([opts.compiler,*flags,'-c',str(p),'-o',obj],cwd=opts.build_dir,check=True)
        if i%15==0 or i==len(order): print(f'Compiled {i}/{len(order)}: {p.name}',flush=True)
    else:
        objects=[str(opts.build_dir/(p.stem+'.o')) for p in order]
        # Full legacy convolution and FFT suites import public scifort_signal.
        for t in ('test_signal_convolve','test_signal_fftconvolve'):
            subprocess.run([opts.compiler,*flags,str(ROOT/'test'/f'{t}.f90'),
                            *objects,'-o',t],cwd=opts.build_dir,check=True)
            subprocess.run([str(opts.build_dir/t)],cwd=opts.build_dir,check=True)
        driver=opts.build_dir/'mixed_public_demo.f90'
        driver.write_text('''program demo
 use scifort_signal, only: dp, convolve, correlate, fftconvolve, oaconvolve, choose_conv_method
 implicit none
 real(dp) :: x(3) = [1.0_dp, 2.0_dp, 3.0_dp]
 complex(dp) :: h(2) = [(1.0_dp,2.0_dp),(3.0_dp,-1.0_dp)]
 complex(dp), allocatable :: y(:)
 character(6) :: algo
 y = convolve(x,h,method='direct')
 if (size(y)/=4) error stop 'public convolve'
 y = correlate(h,x)
 if (size(y)/=4) error stop 'public correlate'
 y = fftconvolve(x,h)
 if (size(y)/=4) error stop 'public fftconvolve'
 y = oaconvolve(h,x)
 if (size(y)/=4) error stop 'public oaconvolve'
 algo = choose_conv_method(x,h)
 if (algo/='direct'.and.algo/='fft') error stop 'public selector'
 print *, 'public mixed convolution generic passed'
end program
''')
        subprocess.run([opts.compiler,*flags,str(driver),*objects,'-o','mixed_public_demo'],
                       cwd=opts.build_dir,check=True)
        subprocess.run([str(opts.build_dir/'mixed_public_demo')],cwd=opts.build_dir,check=True)
        print('PASS: public module and legacy one-dimensional convolution tests',flush=True)

if __name__=='__main__': main()
