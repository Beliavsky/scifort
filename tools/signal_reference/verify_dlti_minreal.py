#!/usr/bin/env python3
"""Cross-check SciFort's complex MIMO minimal-realization algorithm with SciPy.

Usage: python verify_dlti_minreal.py /path/to/dlti_minreal_driver
"""
import sys
import subprocess
import numpy as np
from scipy import linalg, signal


def read_matrix(lines,pos,shape):
    n = int(np.prod(shape))
    data = [complex(*map(float,lines[k].split())) for k in range(pos,pos+n)]
    return np.asarray(data).reshape(shape,order='F'),pos+n


def make_case(rng,i):
    r = i % 5  # required minimal order, including static-feedthrough systems
    nextra = (i * 3 + 1) % 5
    n = r+nextra
    nu = 1+(i*7)%3
    ny = 1+(i*11)%3
    dt = (0.05,0.2,1.,3.)[i%4]
    real_case = i % 3 == 0
    def noise(shape):
        z = rng.normal(size=shape)
        return z if real_case else (z+1j*rng.normal(size=shape))/np.sqrt(2)
    # Four disjoint state classes: minimal, unobservable, unreachable, both.
    # Distinct poles avoid accidental degeneracy in Kalman ranks.
    poles = (0.04+0.05*np.arange(n)).astype(complex)
    if not real_case:
        poles += .025j*(-1)**np.arange(n)
    poles[r:] += np.where(np.arange(n-r)%2==0,.12,.22) if nextra else 0.
    a = np.diag(poles)
    b = np.zeros((n,nu),dtype=complex)
    c = np.zeros((ny,n),dtype=complex)
    if r:
        b[:r,:] = noise((r,nu))
        c[:,:r] = noise((ny,r))
    for j in range(nextra):
        idx = r+j
        if j % 3 == 0:  # reachable but unobservable
            b[idx,:] = noise((nu,))
        elif j % 3 == 1:  # observable but unreachable
            c[:,idx] = noise((ny,))
    q,_ = linalg.qr(noise((n,n))) if n else (np.zeros((0,0)),None)
    a = q @ a @ q.conj().T
    b = q @ b
    c = c @ q.conj().T
    d = .13*noise((ny,nu))
    return (n,nu,ny,dt,r,real_case),(a,b,c,d)


def main():
    if len(sys.argv)!=2:
        raise SystemExit(__doc__)
    rng=np.random.default_rng(20261025)
    cases=[make_case(rng,i) for i in range(200)]
    parts=[str(len(cases))]
    for (n,nu,ny,dt,r,_),arrays in cases:
        parts.append(f'{n} {nu} {ny} {dt:.17g}')
        for x in arrays:
            for val in x.ravel(order='F'):
                parts.append(f'({val.real:.17g},{val.imag:.17g})')
    result=subprocess.run([sys.argv[1]],input='\n'.join(parts)+'\n',text=True,
                          capture_output=True,check=True)
    lines=result.stdout.splitlines()
    pos=0
    values=0
    max_error=0.
    for k,((n,nu,ny,dt,r,isreal),arrays) in enumerate(cases):
        case,status,nred=map(int,lines[pos].split())
        pos+=1
        if (case,status,nred)!=(k+1,0,r):
            raise AssertionError(('status/degree',k,case,status,nred,r))
        a,b,c,d=arrays
        a1,pos=read_matrix(lines,pos,(r,r))
        b1,pos=read_matrix(lines,pos,(r,nu))
        c1,pos=read_matrix(lines,pos,(ny,r))
        d1,pos=read_matrix(lines,pos,(ny,nu))
        orig=(a,b,c,d)
        red=(a1,b1,c1,d1)
        max_error=max(max_error,np.max(abs(d-d1)))
        # Original and reduced transfer functions compared against SciPy's
        # independent ss2tf polynomial conversion for each input channel.
        for idx in range(nu):
            if n:
                num,den=signal.ss2tf(a,b,c,d,input=idx)
            else:
                num=np.array([[d[row,idx]] for row in range(ny)])
                den=np.array([1.])
            for w in (0.,.14,.4,.9,1.8,2.7):
                z=np.exp(1j*w)
                href=c1 @ linalg.solve(z*np.eye(r)-a1,b1,assume_a="gen") + d1 if r else d1
                honly=c @ linalg.solve(z*np.eye(n)-a,b,assume_a="gen") + d if n else d
                expected=np.array([np.polyval(num[row],z)/np.polyval(den,z)
                                   for row in range(ny)])
                err=np.max(abs(href[:,idx]-expected))
                if err>1e-8:
                    raise AssertionError((k,'scipy ss2tf',err))
                err2=np.max(abs(href-honly))
                if err2>1e-8:
                    raise AssertionError((k,'frequency matrix',err2))
                max_error=max(max_error,err,err2)
                values += ny + ny*nu
        # Compare transient simulations with SciPy's own dlsim for real systems.
        if isreal:
            inp=rng.normal(size=(28,nu))
            tout,y0,_=signal.dlsim((a.real,b.real,c.real,d.real,dt),inp)
            t2,y1,_=signal.dlsim((a1.real,b1.real,c1.real,d1.real,dt),inp)
            err=np.max(abs(y0-y1))
            if err>1e-8:
                raise AssertionError((k,'SciPy dlsim',err))
            max_error=max(max_error,err)
            values+=y0.size
        values+=r*r+r*nu+ny*r+ny*nu
    assert pos==len(lines),(pos,len(lines))
    print(f'PASS: {len(cases)} MIMO minimal-realization cases, {values} values, max error {max_error:.3e}')

if __name__=='__main__':
    main()
