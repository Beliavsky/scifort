#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compare standalone transforms with SciPy in checked and optimized builds."""
import argparse
from pathlib import Path
import subprocess
import mpmath as mp
import numpy as np
import scipy
from scipy import fft

ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--module-dir',type=Path,required=True)
    parser.add_argument('--library',type=Path,required=True)
    parser.add_argument('--compiler',default='gfortran')
    args=parser.parse_args()
    rng=np.random.default_rng(20261112)
    norms=['backward','forward','ortho']
    names=['fft','ifft','rfft','irfft','hfft','ihfft','dct','idct','dst','idst',
           'fftfreq','rfftfreq','fftshift','ifftshift']
    cases=[]
    inputs=[]
    def add(op,x,n,kind=2,mode=0,orth=-1,d=1.,real_input=False):
        x=np.asarray(x)
        if real_input: x=x.real
        name=names[op-1]
        if op<=6:
            expected=getattr(fft,name)(x,n=n,norm=norms[mode])
        elif op<=10:
            options=dict(type=kind,n=n,norm=norms[mode],
                         orthogonalize=None if orth<0 else bool(orth))
            if np.iscomplexobj(x):
                # SciPy 1.15.3's complex path drops the orthogonalize argument.
                # Apply the public real-transform definition to each component.
                expected=(getattr(fft,name)(x.real,**options)
                          +1j*getattr(fft,name)(x.imag,**options))
            else:
                expected=getattr(fft,name)(x,**options)
        elif op<=12:
            if abs(d)>np.finfo(float).max/n:
                # SciPy/NumPy forms n*d, which overflows for these spacings.
                # Check the defining bin/(n*d) values independently instead.
                bins=np.arange(n//2+1) if op==12 else np.array(
                    [j if j<=(n-1)//2 else j-n for j in range(n)])
                with mp.workdps(100):
                    expected=np.array([float(mp.mpf(int(j))/n/mp.mpf(float(d))) for j in bins])
                name+=' (mpmath)'
            else:
                expected=getattr(fft,name)(n,d=d)
        else:
            expected=getattr(fft,name)(x)
        label=f'{name}: n={n},type={kind},norm={norms[mode]},orth={orth},real={real_input}'
        cases.append((label,expected,max(np.max(np.abs(x),initial=0.),
                                       np.max(np.abs(expected),initial=0.),np.finfo(float).tiny)))
        inputs.append(f'{op} {int(real_input)} {len(x)} {n} {kind} {mode} {orth} {d:.17g}')
        inputs.extend(f'{complex(z).real:.17g} {complex(z).imag:.17g}' for z in x)
    lengths=[1,2,3,4,5,7,8,9,15,16,17,31,64,97,128,257,1009]
    for n in lengths:
        for mode in range(3):
            for nx in sorted(set([0,max(1,n-2),n,n+3])):
                x=rng.normal(size=nx)+1j*rng.normal(size=nx)
                for op in range(1,7):
                    for real_input in ([True] if op in [3,6] else [False,True]):
                        add(op,x,n,mode=mode,real_input=real_input)
    for n in [1,2,3,4,5,8,17,31,64,127,257]:
        for kind in range(1,5):
            for mode in range(3):
                for orth in [-1,0,1]:
                    for op in range(7,11):
                        if n==1 and kind==1 and op in [7,8]: continue
                        nx=n+(1 if (kind+mode)%2 else -1)
                        x=rng.normal(size=nx)+1j*rng.normal(size=nx)
                        for real_input in [False,True]:
                            add(op,x,n,kind,mode,orth,real_input=real_input)
    for n in [1,2,3,8,9,31]:
        for d in [1.,.1,-.25,1e-300,1e300,np.finfo(float).max]:
            for op in [11,12]: add(op,[],n,d=d)
    for n in [0,1,2,3,8,9,31]:
        x=rng.normal(size=n)+1j*rng.normal(size=n)
        for op in [13,14]:
            for real_input in [False,True]: add(op,x,max(1,n),real_input=real_input)
    for amplitude in [1e-100,1e100]:
        x=amplitude*(rng.normal(size=97)+1j*rng.normal(size=97))
        for mode in range(3):
            for op in range(1,11):
                add(op,x,97,mode=mode,real_input=op in [3,6])
    stdin='\n'.join(inputs)+'\n'
    print(f'SciPy {scipy.__version__}; {len(cases)} reference cases',flush=True)
    for mode,flags in [('checked',['-O0','-g','-fcheck=all',
                                  '-ffpe-trap=invalid,zero,overflow','-fbacktrace']),
                       ('optimized',['-O3','-fcheck=bounds'])]:
        directory=ROOT/'build/fft_reference'/mode
        directory.mkdir(parents=True,exist_ok=True)
        exe=directory/'driver.exe'
        subprocess.run([args.compiler,*flags,'-J'+str(directory),'-I'+str(directory),
                        '-I'+str(args.module_dir.resolve()),str(ROOT/'src/fft.f90'),
                        str(ROOT/'validation/fft_reference_driver.f90'),
                        str(args.library.resolve()),'-o',str(exe)],cwd=ROOT,check=True)
        run=subprocess.run([str(exe)],input=stdin,text=True,capture_output=True,check=True)
        lines=iter(run.stdout.splitlines())
        worst=(0.,'')
        for name,expected,scale in cases:
            status,count=map(int,next(lines).split())
            assert status==0 and count==len(expected),(mode,name,status,count,len(expected))
            actual=np.array([complex(*map(float,next(lines).split())) for _ in range(count)])
            assert np.all(np.isfinite(actual)),(mode,name,'nonfinite output')
            error=np.max(np.abs(actual-expected),initial=0.)/scale
            assert error<=2e-12,(mode,name,error,actual,expected)
            worst=max(worst,(error,name))
        assert next(lines,None) is None
        print(f'{mode}: {len(cases)} passed; max scaled error={worst[0]:.3g} ({worst[1]})',flush=True)

if __name__=='__main__':
    main()
