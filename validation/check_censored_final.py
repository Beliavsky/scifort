#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent references for censored stable and studentized-range methods."""
import argparse
import json
import math
from pathlib import Path
import subprocess
import mpmath as mp
from scipy import stats, integrate, special
ROOT=Path(__file__).resolve().parents[1]

def stable(x,a,b):
    x,a,b=map(mp.mpf,(x,a,b))
    phase=lambda t:x*t-b*mp.tan(mp.pi*a/2)*t**a
    pdf=mp.quad(lambda t:mp.exp(-t**a)*mp.cos(phase(t)),[0,1,3,6,12,30,60])/mp.pi
    cdf=mp.mpf('.5')+mp.quad(lambda t:mp.exp(-t**a)*mp.sin(phase(t))/t if t else 0,[0,1,3,6,12,30,60])/mp.pi
    return mp.log(pdf),mp.log(cdf),mp.log1p(-cdf)

def two_means(q,df):
    q,df=map(mp.mpf,(q,df));t=q/mp.sqrt(2)
    if mp.isinf(df):
        cdf=mp.erf(q/2);sf=mp.erfc(q/2)
        ld=-q*q/4-mp.log(mp.pi)/2
        if cdf<sf:return ld,mp.log(cdf),mp.log1p(-cdf)
        return ld,mp.log1p(-sf),mp.log(sf)
    density=mp.sqrt(2)*mp.gamma((df+1)/2)/(mp.sqrt(df*mp.pi)*mp.gamma(df/2))*(1+t*t/df)**(-(df+1)/2)
    sf=mp.betainc(df/2,mp.mpf('.5'),0,df/(df+t*t),regularized=True)
    cdf=mp.betainc(mp.mpf('.5'),df/2,0,t*t/(df+t*t),regularized=True)
    if cdf<sf:return mp.log(density),mp.log(cdf),mp.log1p(-cdf)
    return mp.log(density),mp.log1p(-sf),mp.log(sf)

def adaptive_range(q,k,df,pdf):
    # Positive adaptive integrals, independent of the fixed log-coordinate rule.
    phi=lambda z:math.exp(-z*z/2)/math.sqrt(2*math.pi)
    def conditional(r):
        def inner(z):
            d=(special.ndtr(-z)-special.ndtr(-z-r) if z>=0 else
               special.ndtr(z+r)-special.ndtr(z))
            if pdf:return k*(k-1)*phi(z)*phi(z+r)*d**(k-2)
            return k*phi(z)*d**(k-1)
        return integrate.quad(inner,-12,12,epsabs=1e-30,epsrel=3e-13)[0]
    if math.isinf(df):return math.log(conditional(q))
    c=math.log(2)+df/2*math.log(df/2)-special.gammaln(df/2)
    def outer(s):
        if s==0:return 0.
        return math.exp(c+(df-1)*math.log(s)-df*s*s/2)*conditional(q*s)*(s if pdf else 1)
    return math.log(integrate.quad(outer,0,20,epsabs=1e-30,epsrel=3e-13)[0])

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reuse',action='store_true');args=parser.parse_args()
    cache=ROOT/'build'/'final-censored-reference';cache.mkdir(exist_ok=True)
    cases=[];inputs=[]
    def add(family,a,b,x,part,expected,loc=0.,scale=1.,end=None,label='reference',tol=3e-9):
        cases.append([label,float(expected),tol])
        inputs.append(f'{family} {a:.17g} {loc:.17g} {scale:.17g} '+' '.join(str(int(part==j)) for j in range(4)))
        inputs.append(f'{b:.17g}')
        inputs.append(f'{x:.17g}'+(' '+f'{end:.17g}' if end is not None else ''))
    if args.reuse:
        cases=json.loads((cache/'cases.json').read_text());inputs=(cache/'input.txt').read_text().splitlines()
    else:
        with mp.workdps(70):
            for a,b in ((1.25,-.5),(1.5,0.),(1.5,.5),(1.75,-.5),(1.75,.5)):
                for x in (-6.,-.1,0.,1.,6.):
                    refs=stable(x,a,b)
                    for part in range(3):add(108,a,b,x,part,refs[part],label='MP Fourier',tol=3e-10)
                    if x==1.:
                        end=math.nextafter(x,math.inf)
                        add(108,a,b,x,3,refs[0]+mp.log(mp.mpf(end)-x),end=end,label='MP adjacent')
                left=stable(-.2,a,b)[1];right=stable(.8,a,b)[1]
                add(108,a,b,-.2,3,right+mp.log(-mp.expm1(left-right)),end=.8,label='MP ordinary interval')
                print('Stable MP fixtures ready',a,b,flush=True)
            for df in (1.,2.,10.,30.,31.,50.,75.,100.,math.inf):
                for q in (1e-100,.01,1.,5.,30.,1e100):
                    if math.isinf(df) and q>40:continue
                    refs=two_means(q,df)
                    for part in range(3):add(109,2.,df,q,part,refs[part],label='MP folded t',tol=3e-10)
                q=1.;end=math.nextafter(q,math.inf)
                add(109,2.,df,q,3,two_means(q,df)[0]+mp.log(mp.mpf(end)-q),end=end,label='MP folded adjacent')
            for k,df in ((3.,2.),(5.,10.),(10.,30.),(2.5,5.),(2.5,31.),(5.,50.),(10.,100.),(3.,math.inf),(10.,math.inf)):
                for q in (.1,1.,4.,10.):
                    for part,name in enumerate(('logpdf','logcdf','logsf')):
                        ref=float(getattr(stats.studentized_range,name)(q,k,df))
                        if (q==.1 and part<2) or (q==10. and part==0 and 30<df<=100):
                            ref=adaptive_range(q,k,df,part==0)
                            add(109,k,df,q,part,ref,label='adaptive range',tol=3e-10)
                        elif math.isfinite(ref) and (part!=2 or ref>-18) and not (-1e-8<ref<0):
                            add(109,k,df,q,part,ref,label='SciPy central',tol=3e-8)
                mass=stats.studentized_range.cdf(2.,k,df)-stats.studentized_range.cdf(1.,k,df)
                add(109,k,df,1.,3,math.log(mass),end=2.,label='SciPy interval',tol=3e-8)
                m=mp.mpf(k)-1
                c=mp.log(k)/2-m*mp.log(2*mp.pi)/2
                if math.isfinite(df):c+=m/2*mp.log(2/mp.mpf(df))+mp.loggamma((df+m)/2)-mp.loggamma(df/2)
                q=1e-100
                add(109,k,df,q,1,c+m*mp.log(q),label='MP lower power',tol=3e-10)
                add(109,k,df,q,0,c+mp.log(m)+(m-1)*mp.log(q),label='MP lower density',tol=3e-10)
            # Affine references use exactly represented binary64 inputs.
            for family,a,b,oracle in ((108,1.5,.5,lambda x:stable(x,1.5,.5)),(109,2.,10.,lambda x:two_means(x,10.))):
                for scale in (1e-200,1e200):
                    loc=.3*scale;x=loc+scale
                    z=(mp.mpf(x)-loc)/scale
                    refs=oracle(z)
                    for part in range(3):add(family,a,b,x,part,refs[part]-(mp.log(scale) if part==0 else 0),loc,scale,label='MP affine')
            with mp.workdps(120):
                x=1e300;scale=1e-200;refs=two_means(mp.mpf(x)/scale,10.)
                for part in range(3):add(109,2.,10.,x,part,refs[part]-(mp.log(scale) if part==0 else 0),scale=scale,label='MP overflowed coordinate')
                end=math.nextafter(x,math.inf)
                add(109,2.,10.,x,3,refs[0]-mp.log(scale)+mp.log(mp.mpf(end)-x),scale=scale,end=end,label='MP overflowed adjacent')
        (cache/'cases.json').write_text(json.dumps(cases));(cache/'input.txt').write_text('\n'.join(inputs)+'\n')
    for name,opt in (('checked','-O0'),('optimized','-O3')):
        out=cache/name;out.mkdir(exist_ok=True)
        flags=[opt,'-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow','-fbacktrace',
               '-J'+str(out),'-I'+str(out),'-I'+str(ROOT/'build/gfortran_2806A042DAA8BBF7')]
        sources=['src/empirical_distribution.f90','src/censored.f90','src/fit.f90','src/stats.f90']
        cmd=['gfortran',*flags,'-c',*[str(ROOT/path) for path in sources]]
        print(subprocess.list2cmdline(cmd),flush=True);subprocess.run(cmd,cwd=out,check=True)
        objects=[str(out/(Path(path).stem+'.o')) for path in sources]
        archive=str(ROOT/'build/gfortran_3C159EB802C22FB9/scifort/libscifort.a')
        for source,exe in [('validation/censored_reference_driver.f90','driver.exe'),
                           ('test/test_censored_final.f90','test.exe'),
                           ('test/test_censored_studentized_range_large_k.f90','large_k.exe')]:
            cmd=['gfortran',*flags,str(ROOT/source),*objects,archive,'-o',str(out/exe)]
            print(subprocess.list2cmdline(cmd),flush=True);subprocess.run(cmd,cwd=ROOT,check=True)
        subprocess.run([str(out/'test.exe')],check=True)
        subprocess.run([str(out/'large_k.exe')],check=True)
        result=subprocess.run([str(out/'driver.exe')],input='\n'.join(inputs)+'\n0\n',text=True,capture_output=True)
        if result.returncode:raise RuntimeError(result.stderr)
        actual=list(map(float,result.stdout.split()));assert len(actual)==len(cases),(len(actual),len(cases))
        worst=0.;bad=[]
        for i,(a,(label,b,tol)) in enumerate(zip(actual,cases)):
            error=0. if a==b else abs(a-b)/max(1.,abs(b))
            if label.startswith('MP') and 0<abs(b)<1e-10:error=abs(a/b-1)
            if not math.isfinite(error) or error>tol:bad.append((i,label,a,b,error))
            worst=max(worst,error)
        assert not bad,(name,bad)
        print(name,len(cases),'independent references passed; maximum error',worst,flush=True)
if __name__=='__main__':main()
