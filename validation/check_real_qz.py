# SPDX-License-Identifier: MIT
# Copyright (c) 2026 SciFort contributors
"""Validate real QZ/ordered QZ through reconstruction and homogeneous SciPy spectra."""
from pathlib import Path
import subprocess
import sys
import numpy as np
from scipy.linalg import eigvals, qz, ordqz
from scipy.optimize import linear_sum_assignment


def distance(alpha,beta,ref):
    scale=np.maximum(abs(alpha),abs(beta))
    assert np.all(scale>0)
    aa=alpha/scale;bb=beta/scale
    rs=np.maximum(abs(ref[0]),abs(ref[1]));ra=ref[0]/rs;rb=ref[1]/rs
    cost=abs(aa[:,None]*rb[None,:]-bb[:,None]*ra[None,:])
    i,j=linear_sum_assignment(cost)
    return np.max(cost[i,j],initial=0.)


def main():
    rng=np.random.default_rng(20261124)
    cases=[]
    for n in (0,1,2,3,5,8,16,40):
        for trial in range(4):
            a=rng.normal(size=(n,n));b=rng.normal(size=(n,n))
            if n and trial==1:b[-1]=0.
            if n and trial==2:b=np.diag(np.r_[np.ones(n-1),0.])
            if trial==3:a=a*1e150;b=b*1e-150
            for mode in range(7):cases.append((a,b,mode))
    # Constructed real, conjugate, repeated, zero and infinite eigenvalues.
    for n in (6,12):
        a=np.diag(np.linspace(-3,3,n));b=np.eye(n)
        a[:2,:2]=[[1.,-2.],[2.,1.]]
        a[-1,-1]=2.;b[-1,-1]=0.
        left=np.linalg.qr(rng.normal(size=(n,n)))[0]
        right=np.linalg.qr(rng.normal(size=(n,n)))[0]
        for mode in range(7):cases.append((left@a@right.T,left@b@right.T,mode))
    rows=[]
    for a,b,mode in cases:
        rows.append(f'{len(a)} {mode}')
        rows.extend(' '.join(f'{v:.17g}' for v in row) for row in np.vstack((a,b)))
    proc=subprocess.run([str(Path(sys.argv[1]).resolve())],input='\n'.join(rows)+'\n',
                        text=True,capture_output=True,check=True)
    upstream_lines=None
    if len(sys.argv)>2:
        upstream=subprocess.run([str(Path(sys.argv[2]).resolve())],input='\n'.join(rows)+'\n',
                                text=True,capture_output=True,check=True)
        actual=np.fromstring(proc.stdout,sep=' ')
        expected=np.fromstring(upstream.stdout,sep=' ')
        assert actual.shape==expected.shape
        error=np.max(abs(actual-expected)/(1+abs(expected)),initial=0.)
        print(f'Unchanged upstream raw output difference {error:.6g} (nonunique factors)')
        upstream_lines=iter(upstream.stdout.splitlines())
    lines=iter(proc.stdout.splitlines());spectrum_max=residual_max=orthogonal_max=0.
    sorts=['none','lhp','rhp','iuc','ouc']
    pairs=ordered=reference_rejections=0
    upstream_spectrum_max=0.
    for index,(a,b,mode) in enumerate(cases):
        code,sdim=map(int,next(lines).split());n=len(a)
        assert code==0,(index,n,mode,code,proc.stderr[-1000:])
        mats=[np.array([np.fromstring(next(lines),sep=' ') for _ in range(n)]).reshape(n,n)
              for _ in range(4)]
        aa,bb,q,z=mats
        if mode!=6:
            alpha=np.fromstring(next(lines),sep=' ')+1j*np.fromstring(next(lines),sep=' ')
            beta=np.fromstring(next(lines),sep=' ')
        if upstream_lines is not None:
            refcode,refsdim=map(int,next(upstream_lines).split())
            assert (refcode,refsdim)==(code,sdim)
            original=[np.array([np.fromstring(next(upstream_lines),sep=' ') for _ in range(n)])
                      .reshape(n,n) for _ in range(4)]
            if mode!=6:
                upalpha=np.fromstring(next(upstream_lines),sep=' ')+1j*np.fromstring(next(upstream_lines),sep=' ')
                upbeta=np.fromstring(next(upstream_lines),sep=' ')
        if not n:continue
        sa=max(np.max(abs(a)),1e-300);sb=max(np.max(abs(b)),1e-300)
        for x,y,scale in ((a,aa,sa),(b,bb,sb)):
            error=np.max(abs(q@(y/scale)@z.T-x/scale))
            residual_max=max(residual_max,error)
            assert error<3e-11,(index,'reconstruction',error)
        for basis in (q,z):
            error=np.max(abs(basis.T@basis-np.eye(n)));orthogonal_max=max(orthogonal_max,error)
            assert error<3e-12
        assert np.max(abs(np.tril(aa/sa,-2)),initial=0)==0
        assert np.max(abs(np.tril(bb/sb,-1)),initial=0)==0
        assert np.all(np.diag(bb)>=0)
        for i in range(n-2):assert aa[i+1,i]==0 or aa[i+2,i+1]==0
        for i in range(n-1):
            if aa[i+1,i]!=0:
                pairs+=1
                assert bb[i,i]>0 and bb[i+1,i+1]>0
                assert abs(bb[i,i+1]/sb)<3e-12
        if mode==6:
            alpha,beta=eigvals(aa/sa,bb/sb,homogeneous_eigvals=True)
            ref=eigvals(a/sa,b/sb,homogeneous_eigvals=True)
        else:
            # Compare physical homogeneous outputs in independently scaled units.
            alpha=alpha/sa;beta=beta/sb
            ref=eigvals(a/sa,b/sb,homogeneous_eigvals=True)
        if upstream_lines is not None:
            if mode==6:
                upref=eigvals(original[0]/sa,original[1]/sb,homogeneous_eigvals=True)
            else:upref=np.vstack((upalpha/sa,upbeta/sb))
            up_error=distance(alpha,beta,upref)
            upstream_spectrum_max=max(upstream_spectrum_max,up_error)
            assert up_error<3e-8,(index,'unchanged upstream spectrum',up_error)
        error=distance(alpha,beta,ref);spectrum_max=max(spectrum_max,error)
        assert error<3e-8,(index,'homogeneous spectrum',error)
        if mode:
            if mode==6:sort='lhp'
            elif mode==5:
                sort=lambda ar,be: np.imag(ar)>0
            else:sort=sorts[mode]
            # Select the independent normalized SciPy eigenspectrum in physical units.
            ra,rb=ref[0]*sa,ref[1]*sb
            if callable(sort):chosen=np.imag(ra)!=0
            elif sort=='lhp':chosen=np.real(ra)*rb<0
            elif sort=='rhp':chosen=np.real(ra)*rb>0
            elif sort=='iuc':chosen=abs(ra)<abs(rb)
            else:chosen=abs(ra)>abs(rb)
            assert sdim==np.count_nonzero(chosen),(index,'selection count',sdim,chosen)
            try:
                ordqz(a,b,sort=sort,output='real')
            except ValueError:
                # Some extreme physical scales make SciPy's block reordering fail.
                # Reconstruction, the selected spectrum and counts are still checked.
                reference_rejections+=1
            if sdim:
                first=eigvals(aa[:sdim,:sdim]/sa,bb[:sdim,:sdim]/sb,homogeneous_eigvals=True)
                selected_ref=np.vstack((ra[chosen]/sa,rb[chosen]/sb))
                err=distance(first[0],first[1],selected_ref)
                assert err<3e-8,(index,'selected spectrum',err)
            ordered+=1
    assert next(lines,None) is None
    if upstream_lines is not None:
        assert next(upstream_lines,None) is None
        print(f'Max unchanged upstream homogeneous spectrum difference {upstream_spectrum_max:.6g}')
    print(f'PASS: {len(cases)} real QZ/ordered/factor-reordering cases; {ordered} ordered spectra; '
          f'{pairs} conjugate blocks checked; {reference_rejections} SciPy reorder rejections')
    print(f'Max scaled reconstruction {residual_max:.6g}; orthogonality {orthogonal_max:.6g}; '
          f'homogeneous spectral error {spectrum_max:.6g}')


if __name__=='__main__':main()
