#!/usr/bin/env python3
"""Compare SciFort Gustafsson zero-phase filtering with SciPy 1.17.

Compile driver_signal_filtfilt_gust.f90 with signal_filtfilt_gust.f90 first.
Usage: python compare_signal_filtfilt_gust_scipy.py /path/to/driver
"""
import subprocess
import sys
import warnings

import numpy as np
from scipy import signal, linalg

rng = np.random.default_rng(20261026)
cases = []
for order in range(0, 7):
    for n in (1, 2, 3, 5, 9, 16, 31, 80, 120):
        for irlen in (-1, 1, 2, 4, 8, 13, 24):
            if order == 0:
                a = np.array([1.2])
                b = np.array([0.4])
            else:
                # Real, stable poles; smaller filter gains for higher orders.
                radius = rng.uniform(0.05, 0.7, order)
                a = np.poly(radius)
                b = rng.normal(0, 0.12, size=order+1)
            x = rng.normal(size=n)
            cases.append((0,n,b,a,x,irlen))
            if (n+order)%3 == 0:
                cases.append((1,n,b,a,x+1j*rng.normal(size=n),irlen))

for order in (2, 3, 4):
    b,a = signal.butter(order,0.2)
    for n in (7,15,50,140):
        for irlen in (-1,3,9,35):
            x = rng.normal(size=n)
            cases.append((0,n,b,a,x,irlen))

# Unequal numerator/denominator lengths and more sharply decaying filters.
for n in (2, 3, 8, 20, 80, 240):
    for nb,na in ((1,4),(6,1),(2,5),(5,2),(4,6)):
        for irlen in (-1, 2, 8, 16):
            b = rng.normal(0,0.06,size=nb)
            poles = rng.uniform(0.02,0.65,size=na-1)
            a = np.atleast_1d(np.poly(poles))*rng.uniform(0.4,2.3)
            x = rng.normal(size=n)
            cases.append((0,n,b,a,x,irlen))


def reference_complex_gust(b, a, x, irlen):
    """Independent complex-valued Gustafsson linear-system reference.

    SciPy 1.17 builds its observability matrix as real even for complex filter
    coefficients. This reference uses the mathematically correct complex states.
    """
    n=len(x)
    order=max(len(b),len(a))-1
    if order == 0:
        return (b[0]/a[0])**2*x
    m=n if irlen < 0 or n <= 2*irlen else irlen
    obs=np.empty((m,order), dtype=complex)
    for k in range(order):
        initial=np.zeros(order, dtype=complex)
        initial[k]=1.0
        obs[:,k]=signal.lfilter(b,a,np.zeros(m, dtype=complex),zi=initial)[0]
    s=signal.lfilter(b,a,obs[::-1],axis=0)
    sr=s[::-1]
    obsr=obs[::-1]
    fb=signal.lfilter(b,a,signal.lfilter(b,a,x)[::-1])[::-1]
    bf=signal.lfilter(b,a,signal.lfilter(b,a,x[::-1])[::-1])
    if m == n:
        mat=np.column_stack((sr-obs,obsr-s))
        adj=np.column_stack((sr,obsr))
        delta=bf-fb
    else:
        mat=np.zeros((2*m,2*order),dtype=complex)
        adj=np.zeros_like(mat)
        mat[:m,:order]=sr-obs
        mat[m:,order:]=obsr-s
        adj[:m,:order]=sr
        adj[m:,order:]=obsr
        delta=np.r_[bf[:m]-fb[:m],bf[-m:]-fb[-m:]]
    ic=linalg.lstsq(mat,delta)[0]
    shift=adj@ic
    result=fb.copy()
    if m==n:
        result+=shift
    else:
        result[:m]+=shift[:m]
        result[-m:]+=shift[m:]
    return result

for order in range(1,6):
    for n in (2,3,8,20,80):
        for irlen in (-1,2,7,25):
            roots=rng.uniform(.1,.55,size=order)*np.exp(1j*rng.uniform(-1,1,size=order))
            a=np.poly(roots)
            b=(rng.normal(size=order+1)+1j*rng.normal(size=order+1))*.06
            x=rng.normal(size=n)+1j*rng.normal(size=n)
            cases.append((2,n,b,a,x,irlen))

# Bounded-state estimation on long recordings: only the two short
# edge intervals enter the least-squares solve, not the entire recording.
for mode, length in ((0, 8192), (1, 4096)):
    b, a = signal.butter(4, 0.17)
    x = rng.normal(size=length)
    if mode == 1:
        x = x + 1j*rng.normal(size=length)
    cases.append((mode, length, b, a, x, 48))

def write_values(v, complex_pair=False):
    if complex_pair:
        return ''.join(f'{item.real:.17e} {item.imag:.17e}\n' for item in v)
    return ' '.join(f'{item:.17e}' for item in v)+'\n'

input_text = [str(len(cases))+'\n']
expected = []
for mode,n,b,a,x,irlen in cases:
    input_text.append(f'{mode} {n} {len(b)} {len(a)} {irlen}\n')
    input_text.append(write_values(b,mode==2))
    input_text.append(write_values(a,mode==2))
    input_text.append(write_values(x,mode>=1))
    kw={'method':'gust'}
    if irlen != -1:
        kw['irlen'] = irlen
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        expected.append(reference_complex_gust(b,a,x,irlen) if mode == 2 else
                        signal.filtfilt(b,a,x,**kw))
completed = subprocess.run([sys.argv[1]], input=''.join(input_text),text=True,
                           capture_output=True,check=True,timeout=75)
lines=iter(completed.stdout.splitlines())
worst=(0.0,None)
by_mode={0:[0,0.,0.],1:[0,0.,0.],2:[0,0.,0.]}
failed=[]
count=0
for j,(case,ref) in enumerate(zip(cases,expected)):
    mode,n,b,a,x,irlen=case
    status, length = map(int,next(lines).split())
    samples=np.array([list(map(float,next(lines).split())) for i in range(length)])
    out=(samples[:,0]+1j*samples[:,1]) if mode>=1 else samples[:,0]
    error=np.max(np.abs(out-ref))
    count += n
    by_mode[mode][0]+=1
    by_mode[mode][1]=max(by_mode[mode][1],float(error))
    by_mode[mode][2]=max(by_mode[mode][2],float(error/(1+np.max(np.abs(ref)))))
    if error>worst[0]: worst=(float(error),case)
    tol=2e-7+1e-7*np.max(np.abs(ref))
    if status!=0 or not np.all(np.isfinite(out)) or error>tol:
        failed.append((j,mode,n,len(a),irlen,float(error),float(tol),status))
print(f'cases={len(cases)} samples={count} worst_error={worst[0]:.3g} failures={len(failed)}')
print('case groups (count, absolute error, scaled error):',by_mode)
print('worst case mode,n,num,den,irlen:',
      (worst[1][0],worst[1][1],len(worst[1][2]),len(worst[1][3]),worst[1][-1]))
if failed:
    print('first failures:', *failed[:30],sep='\n')
    raise SystemExit(1)
