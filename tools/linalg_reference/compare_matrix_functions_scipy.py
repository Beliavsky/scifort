#!/usr/bin/env python3
"""Compare standalone SciFort dense sqrtm/signm with SciPy and residuals.

Requires a precompiled driver at --driver. No SciFort Fortran sources are altered.
"""
import argparse
import subprocess
import numpy as np
from scipy.linalg import sqrtm, signm


def cases(seed, count):
    rng = np.random.default_rng(seed)
    for k in range(count):
        n = int(rng.integers(1, 9))
        typ = k % 10
        z = (rng.normal(size=(n,n)) + 1j*rng.normal(size=(n,n)))/np.sqrt(n)
        u = rng.normal(size=(n,n))/np.sqrt(n)
        if typ == 0:
            a = z + 1.5*np.eye(n)
            mode = 1
        elif typ == 1:
            a = z - 1.5*np.eye(n)
            mode = 1
        elif typ == 2:
            a = z@z.conj().T+0.2*np.eye(n)
            mode = 1
        elif typ == 3:
            a = u@u.T+0.5*np.eye(n)
            mode = 3
        elif typ == 4:
            # real symmetric indefinite, including negative real eigenvalues
            q, _ = np.linalg.qr(u)
            a = q@np.diag(rng.choice([-3., -1., 0.5, 2.],n))@q.T
            mode = 3
        elif typ == 5:
            a = np.diag(rng.choice([-4., -1., 1., 3.],n)).astype(float)
            a += 0.02*np.triu(u,1)
            mode = 3
        elif typ == 6:
            a = z + 2*np.eye(n)
            mode = 2
        elif typ == 7:
            a = z - 2*np.eye(n)
            mode = 2
        elif typ == 8:
            a = u
            a += np.diag(rng.choice([-2., -1., 1., 2.],n))
            mode = 4
        else:
            a = np.diag(rng.choice([-2.,-1.,1.,2.],n)).astype(float)
            a += 0.2*np.triu(u,1)
            mode = 4
        yield mode, np.asarray(a,dtype=np.complex128)


def compare(driver, count, seed, verbose):
    examples=list(cases(seed,count))
    lines=[str(count)]
    for mode,a in examples:
        n=a.shape[0];lines.append(f'{mode} {n}')
        for val in a.T.flat:
            lines.append(f'{val.real:.17g} {val.imag:.17g}')
    proc=subprocess.run([driver],input='\n'.join(lines)+'\n',text=True,capture_output=True,check=True)
    out=iter(proc.stdout.splitlines())
    diffs=[]; residuals=[]; failed=[]; count_pass=0; status_counts={}
    for k,(mode,a) in enumerate(examples):
        n=a.shape[0]
        status,res=next(out).split()
        status=int(status); res=float(res)
        status_counts[status]=status_counts.get(status,0)+1
        vals=np.array([complex(*map(float,next(out).split())) for j in range(n*n)])
        y=vals.reshape((n,n),order='F')
        if status:
            failed.append((k,mode,n,status))
            continue
        if mode == 3 and np.allclose(a, a.T, atol=1e-12, rtol=1e-12):
            # Negative real eigenvalues lie exactly on sqrt's branch cut.
            # SciPy's Schur implementation can return +/-i depending on
            # insignificant imaginary roundoff.  Choose the +i principal
            # branch analytically using the symmetric eigendecomposition.
            w, v = np.linalg.eigh(a.real)
            ref = (v * np.sqrt(w.astype(complex))) @ v.T
        else:
            ref=sqrtm(a) if mode in (1,3) else signm(a)
        err=np.linalg.norm(y-ref,ord='fro')/max(1,np.linalg.norm(ref,ord='fro'))
        residual=(np.linalg.norm(y@y-a,ord='fro')/max(1,np.linalg.norm(a,ord='fro')) if mode in (1,3) else
            np.linalg.norm(y@y-np.eye(n),ord='fro')/max(1,n))
        diffs.append(err);residuals.append(residual)
        if err < 2e-7 and residual < 2e-8:
            count_pass+=1
        else:
            failed.append((k,mode,n,err,residual, np.linalg.eigvals(a)))
    if verbose:
        print('failures:',failed[:15])
    print(f'cases={count} passed={count_pass} status_counts={status_counts} max_error={max(diffs) if diffs else 0:.3e} max_residual={max(residuals) if residuals else 0:.3e}')
    if failed:
        raise SystemExit('FAIL, first failures: '+repr(failed[:5]))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--driver',required=True)
    p.add_argument('--cases',type=int,default=300)
    p.add_argument('--seed',type=int,default=20261106)
    p.add_argument('--verbose',action='store_true')
    args=p.parse_args()
    compare(args.driver,args.cases,args.seed,args.verbose)
