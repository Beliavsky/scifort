# SPDX-License-Identifier: MIT
"""Independent Fourier and asymptotic references for stable censored tails.

Run after python validation/check_censored_final.py (or --reuse for its fixtures).
Uses its checked/optimized drivers. No new runtime dependency is introduced.
"""
from pathlib import Path
import json
import argparse
import math
import subprocess
import mpmath as mp
from scipy import stats, optimize

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT/'build/final-censored-reference'


def reference(x, alpha, beta):
    """Real Fourier quadrature near the join; Saenko (2023), (3),(5) remotely."""
    x, a, b = map(mp.mpf, (x, alpha, beta))
    if abs(x) < 64:
        phase = lambda t: x*t-b*mp.tan(mp.pi*a/2)*t**a
        grid = list(range(31))
        pdf = mp.quad(lambda t: mp.exp(-t**a)*mp.cos(phase(t)), grid)/mp.pi
        cdf = mp.mpf('.5')+mp.quad(
            lambda t: mp.exp(-t**a)*mp.sin(phase(t))/t if t else 0, grid)/mp.pi
        return mp.log(pdf), mp.log(cdf), mp.log1p(-cdf)
    if x < 0:
        ld, lc, ls = reference(-x, a, -b)
        return ld, ls, lc
    c = 1+1j*b*mp.tan(mp.pi*a/2)
    phase = mp.pi*a/2+mp.arg(c)
    r = abs(c)/x**a
    # This is an asymptotic series, not an infinite sum. Validate the independent
    # 32- and 48-term truncations before using a fixture.
    def sums(count):
        sf = mp.fsum((-1)**(n+1)*mp.gamma(a*n)/mp.factorial(n)
                     *mp.sin(n*phase)*r**n for n in range(1, count)) / mp.pi
        pdf = mp.fsum((-1)**(n+1)*mp.gamma(a*n+1)/mp.factorial(n)
                      *mp.sin(n*phase)*r**n for n in range(1, count)) / (mp.pi*x)
        return pdf, sf
    pdf, sf = sums(48)
    other = sums(32)
    assert max(abs(pdf/other[0]-1),abs(sf/other[1]-1)) < mp.mpf('1e-35')
    return mp.log(pdf), mp.log1p(-sf), mp.log(sf)


def build_references():
    cases = []
    inputs = []
    def add(a,b,x,part,ref,loc=0.,scale=1.,end=None,label='MP tail'):
        cases.append((label,float(ref),3e-10))
        inputs.append(f'108 {a:.17g} {loc:.17g} {scale:.17g} '+
                      ' '.join(str(int(part == j)) for j in range(4)))
        inputs.append(f'{b:.17g}')
        inputs.append(f'{x:.17g}'+(' '+f'{end:.17g}' if end is not None else ''))
    with mp.workdps(90):
        for a,b in ((1.25,-.5),(1.25,.5),(1.5,0.),(1.5,.5),(1.75,-.5),(1.75,.5),
                    (1.37,-.23),(1.68,.31)):
            for x in (math.nextafter(6.,math.inf),8.,16.,64.,100.,1e10,1e100,1e300):
                refs = reference(x,a,b)
                for part in range(3):
                    add(a,b,x,part,refs[part])
                    reflected = reference(-x,a,-b) if x >= 64 else (refs[0],refs[2],refs[1])
                    add(a,-b,-x,part,reflected[part],label='MP reflected tail')
                end = math.nextafter(x,math.inf)
                add(a,b,x,3,refs[0]+mp.log(mp.mpf(end)-x),end=end,label='MP tail adjacent')
                if x >= 64:
                    end = 1.01*x
                    rs = reference(end,a,b)[2]
                    add(a,b,x,3,refs[2]+mp.log(-mp.expm1(rs-refs[2])),
                        end=end,label='MP tail interval')
            print('Stable tail fixtures ready',a,b,flush=True)
        for x,scale in ((1e300,1e-200),(-1e300,1e-200),(1e200,1e200)):
            loc = .3*scale
            refs = reference((mp.mpf(x)-loc)/scale,1.5,.5)
            for part in range(3):
                add(1.5,.5,x,part,refs[part]-(mp.log(scale) if part == 0 else 0),loc,scale,
                    label='MP tail affine')
            end = math.nextafter(x,math.inf)
            add(1.5,.5,x,3,refs[0]-mp.log(scale)+mp.log(mp.mpf(end)-x),loc,scale,end,
                label='MP tail affine adjacent')
    (CACHE/'stable-tail-cases.json').write_text(json.dumps(cases))
    (CACHE/'stable-tail-input.txt').write_text('\n'.join(inputs)+'\n0\n')
    return cases, inputs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reuse',action='store_true',help='reuse existing tail fixtures')
    args = parser.parse_args()
    if args.reuse:
        cases = json.loads((CACHE/'stable-tail-cases.json').read_text())
        inputs = (CACHE/'stable-tail-input.txt').read_text().splitlines()[:-1]
    else:
        cases, inputs = build_references()
    for profile in ('checked','optimized'):
        result = subprocess.run([str(CACHE/profile/'driver.exe')],
                                input='\n'.join(inputs)+'\n0\n',text=True,capture_output=True,check=True)
        actual = list(map(float,result.stdout.split()))
        assert len(actual) == len(cases)
        worst = 0.
        for value,(label,expected,tol) in zip(actual,cases):
            error = 0. if value == expected else abs(value-expected)/max(1.,abs(expected))
            if 0 < abs(expected) < 1e-10:
                error = abs(value/expected-1)
            assert math.isfinite(error) and error <= tol,(profile,label,value,expected,error)
            worst = max(worst,error)
        print(profile,len(cases),'tail references passed; maximum error',worst,flush=True)
    # Moderate-tail SciPy values suffice for an independent free-scale fit.
    def objective(scale):
        rv = stats.levy_stable(1.5,.5,scale=scale)
        return -(sum(rv.logpdf([.5,1.,2.,8.]))+rv.logcdf(-10.)+rv.logsf(20.)+
                 math.log(rv.sf(30.)-rv.sf(32.)))
    fit = optimize.minimize_scalar(objective,bounds=(.5,20.),method='bounded',
                                   options={'xatol':1e-10})
    print('Independent SciPy tail fit: scale',fit.x,'NNLF',fit.fun,flush=True)


if __name__ == '__main__':
    main()
