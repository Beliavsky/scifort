# SPDX-License-Identifier: MIT
"""Independent references for the extended stable censored skewness domain.

Run check_censored_final.py --reuse first to compile the numerical drivers.
Use --compile to rebuild the drivers and run the native regression directly.
Reference formulas are shared with check_censored_stable_tails.py; see its
provenance and stable_censored_skew_2026-11-19.md.
"""
import argparse
from pathlib import Path
import json
import math
import subprocess
import mpmath as mp
from scipy import stats, optimize
from check_censored_stable_tails import reference, CACHE


def build_references():
    cases, inputs = [], []
    def add(a,b,x,part,ref,loc=0.,scale=1.,end=None,label='MP skew'):
        cases.append((label,float(ref),3e-10))
        inputs.append(f'108 {a:.17g} {loc:.17g} {scale:.17g} '+
                      ' '.join(str(int(part == j)) for j in range(4)))
        inputs.append(f'{b:.17g}')
        inputs.append(f'{x:.17g}'+(' '+f'{end:.17g}' if end is not None else ''))
    with mp.workdps(90):
        for a,b in ((1.25,-.75),(1.25,.75),(1.5,-.75),(1.5,.75),
                    (1.75,-.75),(1.75,.75),(1.37,-.625),(1.68,.625)):
            for x in (0.,1.,6.,math.nextafter(6.,math.inf),8.,16.,64.,1e100,1e300):
                refs = reference(x,a,b)
                for part in range(3):
                    add(a,b,x,part,refs[part])
                    add(a,-b,-x,part,refs[(0,2,1)[part]],label='MP skew reflection')
                if x:
                    end = math.nextafter(x,math.inf)
                    add(a,b,x,3,refs[0]+mp.log(mp.mpf(end)-x),end=end,
                        label='MP skew adjacent')
                if x >= 64:
                    end = 1.01*x
                    rs = reference(end,a,b)[2]
                    add(a,b,x,3,refs[2]+mp.log(-mp.expm1(rs-refs[2])),end=end,
                        label='MP skew interval')
            left = reference(-.2,a,b)[1]
            right = reference(.8,a,b)[1]
            add(a,b,-.2,3,right+mp.log(-mp.expm1(left-right)),end=.8,
                label='MP skew central interval')
            left = reference(5.99,a,b)[2]
            right = reference(6.01,a,b)[2]
            add(a,b,5.99,3,left+mp.log(-mp.expm1(right-left)),end=6.01,
                label='MP skew join interval')
            print('Stable skew fixtures ready',a,b,flush=True)
        for b in (-.75,.75):
            for x,scale in ((1e300,1e-200),(-1e300,1e-200),(1e200,1e200)):
                loc = .3*scale
                refs = reference((mp.mpf(x)-loc)/scale,1.25,b)
                for part in range(3):
                    add(1.25,b,x,part,refs[part]-(mp.log(scale) if part == 0 else 0),
                        loc,scale,label='MP skew affine')
                end = math.nextafter(x,math.inf)
                add(1.25,b,x,3,refs[0]-mp.log(scale)+mp.log(mp.mpf(end)-x),
                    loc,scale,end,label='MP skew affine adjacent')
    (CACHE/'stable-skew-cases.json').write_text(json.dumps(cases))
    (CACHE/'stable-skew-input.txt').write_text('\n'.join(inputs)+'\n0\n')
    return cases, inputs



def compile_drivers():
    """Build the affected modules and native test with runtime checks in both profiles."""
    root = Path(__file__).resolve().parents[1]
    sources = ['src/empirical_distribution.f90','src/censored.f90',
               'src/fit.f90','src/stats.f90']
    for profile,opt in (('checked','-O0'),('optimized','-O3')):
        out = CACHE/profile
        out.mkdir(parents=True,exist_ok=True)
        flags = [opt,'-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow','-fbacktrace',
                 '-J'+str(out),'-I'+str(out),'-I'+str(root/'build/gfortran_2806A042DAA8BBF7')]
        cmd = ['gfortran',*flags,'-c',*[str(root/path) for path in sources]]
        print(subprocess.list2cmdline(cmd),flush=True)
        subprocess.run(cmd,cwd=out,check=True)
        objects = [str(out/(Path(path).stem+'.o')) for path in sources]
        archive = str(root/'build/gfortran_3C159EB802C22FB9/scifort/libscifort.a')
        for source,exe in (('validation/censored_reference_driver.f90','driver.exe'),
                           ('test/test_censored_final.f90','test.exe')):
            cmd = ['gfortran',*flags,str(root/source),*objects,archive,'-o',str(out/exe)]
            print(subprocess.list2cmdline(cmd),flush=True)
            subprocess.run(cmd,cwd=root,check=True)
        subprocess.run([str(out/'test.exe')],check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reuse',action='store_true')
    parser.add_argument('--compile',action='store_true',help='compile and run both native profiles')
    parser.add_argument('--profile',choices=('both','checked','optimized'),default='both')
    args = parser.parse_args()
    if args.compile:
        compile_drivers()
    if args.reuse:
        cases = json.loads((CACHE/'stable-skew-cases.json').read_text())
        inputs = (CACHE/'stable-skew-input.txt').read_text().splitlines()[:-1]
    else:
        cases, inputs = build_references()
    profiles = ('checked','optimized') if args.profile == 'both' else (args.profile,)
    for profile in profiles:
        result = subprocess.run([str(CACHE/profile/'driver.exe')],
                                input='\n'.join(inputs)+'\n0\n',text=True,
                                capture_output=True,check=True)
        actual = list(map(float,result.stdout.split()))
        assert len(actual) == len(cases)
        worst = 0.
        for value,(label,expected,tol) in zip(actual,cases):
            error = 0. if value == expected else abs(value-expected)/max(1.,abs(expected))
            if 0 < abs(expected) < 1e-10:
                error = abs(value/expected-1)
            assert math.isfinite(error) and error <= tol,(profile,label,value,expected,error)
            worst = max(worst,error)
        print(profile,len(cases),'skew references passed; maximum error',worst,flush=True)
    def objective(scale):
        rv = stats.levy_stable(1.5,.75,scale=scale)
        return -(sum(rv.logpdf([.5,1.,2.,8.]))+rv.logcdf(-10.)+rv.logsf(20.)+
                 math.log(rv.sf(30.)-rv.sf(32.)))
    fit = optimize.minimize_scalar(objective,bounds=(.5,20.),method='bounded',
                                   options={'xatol':1e-10})
    assert fit.success
    print('Independent SciPy skew fit: scale',fit.x,'NNLF',fit.fun,flush=True)


if __name__ == '__main__':
    main()
