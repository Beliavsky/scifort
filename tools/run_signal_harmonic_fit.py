#!/usr/bin/env python3
"""Check SciFort's harmonic coefficients and predictions against NumPy QR.

Run: python tools/run_signal_harmonic_fit.py
No installed Python package is needed beyond NumPy and SciPy for reference.
"""
from __future__ import annotations

from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.signal import lombscargle

ROOT = Path(__file__).resolve().parents[1]
RNG = np.random.default_rng(20261028)


def make_cases():
    cases = []
    for j in range(256):
        nterms = 1 + j % 4
        n = int(RNG.integers(38, 104))
        m = int(RNG.integers(2, 9))
        x = np.sort(RNG.uniform(-4.0, 14.0, size=n))
        freqs = RNG.uniform(0.31, 4.0, size=m)
        xnew = RNG.uniform(-7.0, 18.0, size=7)
        floating = j % 2 == 0
        center = j % 3 == 0
        weighted = j % 5 != 0
        weights = RNG.uniform(0.02, 3.0, size=n)
        if j % 7 == 0:
            weights[:n//7] = 0.0
        if j % 11 == 0:
            weights *= 1e6
        y = RNG.normal(0.0, 0.3, n) + RNG.uniform(-1.5, 2.5)
        for h in range(1, nterms+1):
            y += RNG.normal() * np.cos(h*freqs[0]*x)
            y += RNG.normal() * np.sin(h*freqs[0]*x)
        cases.append((x,y,freqs,weights,xnew,nterms,center,floating,weighted))
    return cases


def expected(case):
    x,y,freqs,w0,xnew,nterms,center,floating,weighted = case
    w = w0.copy() if weighted else np.ones(len(x))
    w /= np.sum(w)
    shift = float(np.mean(y)) if center else 0.0
    centered = y-shift
    baseline = float(w@centered) if floating else 0.0
    residual = centered-baseline
    variance = float(w@(residual*residual))
    out = []
    for omega in freqs:
        terms = [a for h in range(1,nterms+1) for a in
                 (np.cos(h*omega*x),np.sin(h*omega*x))]
        basis = np.column_stack(terms)
        mean_basis = w@basis if floating else np.zeros(2*nterms)
        design = basis-mean_basis
        coefficients = np.linalg.lstsq(np.sqrt(w)[:,None]*design,
                                        np.sqrt(w)*residual, rcond=None)[0]
        offset = shift+baseline-float(mean_basis@coefficients)
        model = offset + basis@coefficients
        future_terms = [a for h in range(1,nterms+1) for a in
                        (np.cos(h*omega*xnew),np.sin(h*omega*xnew))]
        future = offset + np.column_stack(future_terms)@coefficients
        improvement = float(w@((design@coefficients)**2))
        improvement = min(variance, max(0.0, improvement))
        power = 0.5*len(x)*improvement
        normalized_power = improvement/variance
        rms = float(np.sqrt(w@((y-model)**2)))
        out.append((offset,coefficients,power,normalized_power,rms,future))
    if nterms == 1:
        kw = dict(precenter=center, floating_mean=floating, normalize='power')
        if weighted:
            kw['weights'] = w0
        scipy_p = lombscargle(x,y,freqs,**kw)
        np.testing.assert_allclose([r[2] for r in out],scipy_p,atol=1e-9,rtol=1e-9)
    return out


def input_text(cases):
    lines = [str(len(cases))]
    for x,y,f,w,xnew,nterms,center,floating,weighted in cases:
        lines.append(f'{len(x)} {len(f)} {len(xnew)} {nterms} '
                     f'{int(center)} {int(floating)} {int(weighted)}')
        for data in [x,y,f,w,xnew]:
            lines.append(' '.join(f'{v:.17e}' for v in data))
    return '\n'.join(lines)+'\n'


def verify(exe, cases, expectations):
    proc = subprocess.run([str(exe)], input=input_text(cases), text=True,
                          capture_output=True, check=True)
    it = iter(proc.stdout.splitlines())
    worst = 0.0
    nvalues = 0
    for j,(case,expected_fits) in enumerate(zip(cases,expectations)):
        x,y,f,w,xnew,nterms,center,floating,weighted = case
        status,harmonics,nfreq = map(int,next(it).split())
        assert status == 0 and harmonics == nterms and nfreq == len(f), (j,status)
        for k,(offset,coefficients,power,npower,rms,prediction) in enumerate(expected_fits):
            assert int(next(it)) == 1
            got_offset = float(next(it))
            got_coefficients = np.array([float(next(it)) for _ in range(2*nterms)])
            got_power = float(next(it))
            got_norm = float(next(it))
            got_rms = float(next(it))
            assert int(next(it)) == 0
            got_pred = np.array([float(next(it)) for _ in xnew])
            pairs = [(np.asarray(got_offset),np.asarray(offset),'offset'),
                     (got_coefficients,coefficients,'harmonics'),
                     (np.asarray(got_power),np.asarray(power),'power'),
                     (np.asarray(got_norm),np.asarray(npower),'normalized'),
                     (np.asarray(got_rms),np.asarray(rms),'RMS'),
                     (got_pred,prediction,'prediction')]
            for got,ref,name in pairs:
                np.testing.assert_allclose(got,ref,atol=5e-10,rtol=5e-10,
                                           err_msg=f'case={j}, bin={k}, {name}')
                worst = max(worst,float(np.max(np.abs(got-ref))))
                nvalues += got.size
    assert next(it,None) is None
    return nvalues,worst


def run_build(directory,optimized,cases,expectations):
    directory.mkdir()
    flags = (['-O3','-march=native'] if optimized else
             ['-O0','-g','-fcheck=all','-fbacktrace'])
    common = ['gfortran','-std=f2018',*flags,'-Wall','-Wextra',
              '-Werror=implicit-interface','-J',str(directory),'-I',str(directory)]
    sources = [ROOT/'src/kinds.f90',ROOT/'src/signal_irregular.f90']
    for name,tail in [('native',ROOT/'test/test_signal_harmonic_fit.f90'),
                      ('reference',ROOT/'tools/reference_signal_harmonic_fit.f90'),
                      ('existing',ROOT/'test/test_signal_multiharmonic.f90')]:
        exe = directory/name
        subprocess.run([*common,*(str(p) for p in sources),str(tail),'-o',str(exe)],
                       cwd=directory,check=True,text=True,capture_output=True)
        if name == 'reference':
            nvalues,worst = verify(exe,cases,expectations)
            print(f'  {len(cases)} independent cases / {nvalues} values; '
                  f'max absolute error {worst:.3e}')
        else:
            output = subprocess.run([str(exe)],check=True,capture_output=True,text=True)
            print(' ',output.stdout.strip())


def main():
    cases = make_cases()
    expected_values = [expected(case) for case in cases]
    with tempfile.TemporaryDirectory(prefix='scifort_harmonic_fit_') as tmp:
        for optimized in (False,True):
            print('optimized -O3:' if optimized else 'checked -fcheck=all:')
            run_build(Path(tmp)/('opt' if optimized else 'checked'),
                      optimized,cases,expected_values)
    print('All independent NumPy and single-harmonic SciPy comparisons passed.')


if __name__ == '__main__':
    main()
