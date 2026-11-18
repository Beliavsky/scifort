#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Multiprecision defining-integral checks for censored noncentral t and Landau."""
import argparse
from functools import lru_cache
import json
import math
from pathlib import Path
import random
import subprocess
import mpmath as mp
from scipy import stats, optimize

ROOT = Path(__file__).resolve().parents[1]


def complement(x):
    return mp.log1p(-mp.exp(x)) if x < -mp.log(2) else mp.log(-mp.expm1(x))


@lru_cache(None)
def nct(x, df, nc, loc=0., scale=1.):
    x, df, nc, loc, scale = map(mp.mpf, (x, df, nc, loc, scale))
    z = (x-loc)/scale
    if z < 0:
        ld, lp, lq = nct(float(-x),float(df),float(-nc),float(-loc),float(scale))
        return ld, lq, lp
    a = df/2
    if 0 < z < mp.mpf('1e-40'):
        # Local defining Taylor limit; error is negligible at binary64 subnormal z.
        ld = mp.loggamma(a+mp.mpf('.5'))-mp.loggamma(a)-mp.log(df*mp.pi)/2-nc*nc/2
        cdf = mp.erfc(nc/mp.sqrt(2))/2+z*mp.exp(ld)
        return ld-mp.log(scale),mp.log(cdf),mp.log1p(-cdf)
    d = df+z*z
    c = nc*z/mp.sqrt(d)
    # Integrating the normal numerator first gives a positive normal moment.
    moment = mp.quad(lambda y: y**df*mp.exp(-y*y/2+c*y),
                     [0, 1, 5, 10, 20, 40, 80], maxdegree=5)
    ld = (mp.log(2)+a*mp.log(a)-mp.loggamma(a)-mp.log(2*mp.pi)/2
          -nc*nc/2-(df+1)*mp.log(d)/2+mp.log(moment)-mp.log(scale))
    if z == 0:
        return ld, mp.log(mp.erfc(nc/mp.sqrt(2))/2), mp.log(mp.erfc(-nc/mp.sqrt(2))/2)
    baseline = -df*max(0,mp.log(z))-min(nc,0)**2/2
    def integrand(w):
        if w == 0:
            return mp.mpf(0)
        h = a*(w/z)**2
        g = mp.gammainc(a,0,h,regularized=True)
        return mp.exp(-(w-nc)**2/2+mp.log(g)-baseline)/mp.sqrt(2*mp.pi)
    lq = baseline+mp.log(mp.quad(integrand,[0,.1,1,5,10,20,40,80],maxdegree=5))
    if lq > -mp.log(2):
        def lower(w):
            h = a*(w/z)**2
            return mp.exp(-(w-nc)**2/2)*mp.gammainc(a,h,mp.inf,regularized=True)
        lp = mp.log(mp.erfc(nc/mp.sqrt(2))/2+
                    mp.quad(lower,[0,.1,1,5,10,20,40,80],maxdegree=5)/mp.sqrt(2*mp.pi))
        lq = complement(lp)
    else:
        lp = complement(lq)
    return ld, lp, lq


@lru_cache(None)
def landau(x, loc=0., scale=1.):
    x, loc, scale = map(mp.mpf,(x,loc,scale))
    z = (x-loc)/scale
    if z <= -1:
        sigma = mp.exp(-mp.pi*z/2-1)*2/mp.pi
        if sigma > mp.mpf('1e30'):
            ld = -sigma+mp.log(sigma)/2+mp.log(mp.pi/2)-mp.log(2*mp.pi)/2
            lp = -sigma-mp.log(sigma)/2-mp.log(2*mp.pi)/2
        else:
            root = mp.sqrt(sigma)
            def core(v):
                u = v/root
                return mp.exp(sigma*((1+1j*u)*mp.log(1+1j*u)-1j*u))
            # Independent adaptive complex contour integration, no tabulated rule.
            sd = mp.quad(lambda v: mp.re(core(v)),[0,1,4,10,20,40,80],maxdegree=5)
            sp = mp.quad(lambda v: mp.re(core(v)/(1+1j*v/root)),
                         [0,1,4,10,20,40,80],maxdegree=5)
            ld = -sigma+mp.log(sigma)/2-mp.log(2)+mp.log(sd)
            lp = -sigma-mp.log(sigma)/2-mp.log(mp.pi)+mp.log(sp)
        return ld-mp.log(scale),lp,complement(lp)
    if z >= 1:
        lam = mp.log(mp.pi/2)+mp.pi*z/2
        def base(u):
            if not u: return mp.mpf(1)
            t = u/lam
            return mp.exp(-u-t*mp.log(t))*mp.sin(mp.pi*t)/(mp.pi*t)
        # Scale out the small powers before adaptive quadrature.
        pdf = mp.pi*mp.quad(lambda u:u*base(u),[0,1,5,20,60,200],maxdegree=5)/(2*lam*lam)
        sf = mp.quad(base,[0,1,5,20,60,200],maxdegree=5)/lam
        lq = mp.log(sf)
        return mp.log(pdf)-mp.log(scale),complement(lq),lq
    def phase(t): return t*(z+2/mp.pi*mp.log(t)) if t else mp.mpf(0)
    pdf = mp.quad(lambda t: mp.exp(-t)*mp.cos(phase(t)),[0,1,5,20,60,200],maxdegree=6)/mp.pi
    cdf = mp.mpf('.5')+mp.quad(lambda t:mp.exp(-t)*mp.sin(phase(t))/t if t else 0,
                                [0,1,5,20,60,200],maxdegree=6)/mp.pi
    return mp.log(pdf)-mp.log(scale),mp.log(cdf),mp.log1p(-cdf)


def interval(oracle, a, b, args, loc, scale):
    a, b = mp.mpf(a),mp.mpf(b)
    left = oracle(float(a),*args,loc,scale)
    right = oracle(float(b),*args,loc,scale)
    mid = (a+b)/2
    # Only use a local density when its relative variation is negligible.
    width = b-a
    z = (mid-mp.mpf(loc))/mp.mpf(scale)
    variation = width/mp.mpf(scale)
    if oracle is landau and z < -1:
        variation *= mp.exp(-mp.pi*z/2-1)
    else:
        variation /= max(1,abs(z))
    if variation < mp.mpf('1e-10'):
        return (left[0]+right[0])/2+mp.log(width)
    if right[1] < -mp.log(2):
        return right[1]+mp.log(-mp.expm1(left[1]-right[1]))
    return left[2]+mp.log(-mp.expm1(right[2]-left[2]))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reuse',action='store_true')
    args=parser.parse_args()
    cache=ROOT/'build'/'specialist-reference'
    cache.mkdir(exist_ok=True)
    cases,inputs=[],[]
    def add(family, shapes, x, part, loc=0.,scale=1.,end=None,label='MP'):
        oracle=nct if family==105 else landau
        expected=(oracle(x,*shapes,loc,scale)[part] if part<3 else
                  interval(oracle,x,end,shapes,loc,scale))
        cases.append((family,label,float(expected)))
        counts=[int(part==j) for j in range(4)]
        first=shapes[0] if shapes else 0.
        inputs.append(f'{family} {first:.17g} {loc:.17g} {scale:.17g} '+' '.join(map(str,counts)))
        if shapes:inputs.append(f'{shapes[1]:.17g}')
        inputs.append(f'{x:.17g}'+(' '+f'{end:.17g}' if end is not None else ''))
    if args.reuse:
        cases=json.loads((cache/'cases.json').read_text())
        inputs=(cache/'input.txt').read_text().splitlines()
    else:
        with mp.workdps(100):
            for shapes in ((.5,-10.),(.5,10.),(1.,-2.),(2.,0.),(10.,2.),(100.,10.),(100.,-10.)):
                for x in (-100.,-1.,0.,1.,100.):
                    for part in range(3):add(105,shapes,x,part)
                for x in (-1e300,1e300):
                    for part in (0,1,2):add(105,shapes,x,part,label='extreme')
                add(105,shapes,.2,3,end=math.nextafter(.2,math.inf),label='adjacent')
                print('NCT fixtures ready:',shapes,flush=True)
            for x in (-440.,-100.,-20.,-10.,-3.,-2.7,-1.,0.,1.,10.,1e5,1e100,1e300):
                for part in range(3):add(106,(),x,part,label='Landau')
                if abs(x)<100: add(106,(),x,3,end=math.nextafter(x,math.inf),label='adjacent')
                print('Landau fixtures ready:',x,flush=True)
            for family,shapes in ((105,(.5,-10.)),(106,())):
                for scale in (1e-200,1.,1e200):
                    for sign in (-1,1):
                        loc=sign*1e12*scale
                        x=loc+.7*scale
                        for part in range(3):add(family,shapes,x,part,loc,scale,label='affine')
                        add(family,shapes,x,3,loc,scale,math.nextafter(x,math.inf),'affine adjacent')
                print('Affine fixtures ready:',family,flush=True)
                add(family,shapes,math.nextafter(0.,1.),3,end=math.nextafter(0.,1.)*2,label='subnormal')
                add(family,shapes,-.4,3,end=.8,label='ordinary interval')
            # Direct high-precision subtraction retains finite logs with an overflowed z.
            for family,shapes in ((105,(10.,2.)),(106,())):
                add(family,shapes,1e300,0,scale=1e-200,label='overflowed coordinate')
                add(family,shapes,1e300,2,scale=1e-200,label='overflowed coordinate')
            for family,shapes in ((105,(.5,-10.)),(105,(100.,10.)),(106,())):
                for sign in (-1,1):
                    if family==106 and sign<0:continue
                    x=sign*1e300
                    for delta in (.01,.0001,1e-8,0.):
                        end=x+abs(x)*delta if delta else math.nextafter(x,math.inf)
                        add(family,shapes,x,3,scale=1e-200,end=end,label='overflowed interval')
            for family,shapes in ((105,(.5,-10.)),(106,())):
                add(family,shapes,1.795e300,3,scale=1e-8,end=1.8e300,
                    label='partially overflowed interval')
        (cache/'cases.json').write_text(json.dumps(cases),encoding='ascii')
        (cache/'input.txt').write_text('\n'.join(inputs)+'\n',encoding='ascii')
        print('MP fixtures ready:',len(cases),flush=True)
    for name,opt in (('checked','-O0'),('optimized','-O3')):
        out=cache/name;out.mkdir(exist_ok=True)
        cmd=['gfortran',opt,'-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow',
             '-J'+str(out),'-I'+str(out),'-Ibuild/gfortran_2806A042DAA8BBF7',
             'src/censored.f90','validation/censored_reference_driver.f90',
             'build/gfortran_3C159EB802C22FB9/scifort/libscifort.a','-o',str(out/'driver.exe')]
        print(subprocess.list2cmdline(cmd),flush=True)
        subprocess.run(cmd,cwd=ROOT,check=True)
        result=subprocess.run([str(out/'driver.exe')],input='\n'.join(inputs)+'\n0\n',text=True,capture_output=True)
        if result.returncode:raise RuntimeError((len(result.stdout.split()),result.stderr))
        actual=list(map(float,result.stdout.split()));assert len(actual)==len(cases)
        worst=0.;failures=[]
        for value,case in zip(actual,cases):
            family,label,expected=case
            error=abs(value-expected)/max(1.,abs(expected))
            if 0<abs(expected)<1e-10:error=abs(value/expected-1)
            if value==expected:error=0.
            if not math.isfinite(error) or error>3e-11:failures.append((case,value,error))
            worst=max(worst,error)
        assert not failures,(name,failures)
        print(name,len(cases),'MP cases passed, maximum error',worst,flush=True)
        scipy_check(out/'driver.exe')
    sample=stats.CensoredData(uncensored=[.2,.4,.6],left=[.1],right=[.8],interval=[[.3,.35],[.5,.55]])
    for dist,shape in ((stats.nct,(4.,1.)),(stats.landau,())):
        opt=optimize.minimize_scalar(lambda s:dist._penalized_nnlf((*shape,0.,s),sample),
                                    bounds=(.3,3.),method='bounded',options={'xatol':1e-13})
        print(dist.name,'scale',opt.x,'NNLF',opt.fun,flush=True)


def scipy_check(driver):
    rng=random.Random(2026111809);inputs=[];references=[]
    for family,dist in ((105,stats.nct),(106,stats.landau)):
        for _ in range(30):
            shapes=(rng.uniform(.5,100),rng.uniform(-2,2)) if family==105 else ()
            x=rng.uniform(-1,3)
            for part,method in enumerate(('logpdf','logcdf','logsf')):
                expected=float(getattr(dist,method)(x,*shapes))
                references.append(expected)
                inputs.append(f'{family} {shapes[0] if shapes else 0} 0 1 '+
                              ' '.join(str(int(part==j)) for j in range(4)))
                if shapes:inputs.append(str(shapes[1]))
                inputs.append(str(x))
    result=subprocess.run([str(driver)],input='\n'.join(inputs)+'\n0\n',capture_output=True,text=True,check=True)
    values=list(map(float,result.stdout.split()));assert len(values)==len(references)
    assert all(abs(v-r)<=3e-10*max(1,abs(r)) for v,r in zip(values,references))
    print('SciPy ordinary comparisons passed:',len(values),flush=True)


if __name__=='__main__':main()
