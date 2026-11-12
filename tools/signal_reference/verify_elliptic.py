#!/usr/bin/env python3
"""Compare SciFort elliptic IIR design/order selection with SciPy 1.17.

Compile and run a temporary Fortran driver; no SciPy dependency is introduced
for users of the Fortran library.

    python tools/signal_reference/verify_elliptic.py
"""
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path
import warnings

import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[2]
MODULES = [
    'kinds', 'signal_dpss', 'signal_windows_extra', 'signal_spectral', 'signal_analog_design', 'signal_design',
    'signal_elliptic', 'signal_iir_extra',
]
SEED = 20261016
RNG = np.random.default_rng(SEED)
BANDS = ('lowpass', 'highpass', 'bandpass', 'bandstop')


def edges_for_order(mode: int) -> tuple[np.ndarray, np.ndarray]:
    """Construct a valid analog/digital band specification in [0, 1]."""
    if mode == 1:
        wp = float(RNG.uniform(.05, .7))
        ws = float(RNG.uniform(wp+.04, .95))
        return np.array([wp]), np.array([ws])
    if mode == 2:
        ws = float(RNG.uniform(.05, .7))
        wp = float(RNG.uniform(ws+.04, .95))
        return np.array([wp]), np.array([ws])
    if mode == 3:
        wp0, wp1 = float(RNG.uniform(.2, .4)), float(RNG.uniform(.55, .75))
        ws0, ws1 = float(RNG.uniform(.05, wp0-.04)), float(RNG.uniform(wp1+.04, .95))
        return np.array([wp0, wp1]), np.array([ws0, ws1])
    wp0, wp1 = float(RNG.uniform(.05, .3)), float(RNG.uniform(.7, .95))
    ws0, ws1 = float(RNG.uniform(wp0+.04, .45)), float(RNG.uniform(.55, wp1-.04))
    return np.array([wp0, wp1]), np.array([ws0, ws1])


def cases():
    """Generate filters, order-selection problems, and minimum-order designs."""
    results = []
    for j in range(1200):
        mode = j % 4 + 1
        analog = bool((j // 4) % 2)
        order = int(RNG.integers(1, 9))
        rp = float(RNG.uniform(.1, 5.0))
        rs = float(RNG.uniform(rp+6, 90))
        fs = 2.0 if analog else float(RNG.choice([2., 100., 1000.]))
        if mode <= 2:
            cutoff = np.array([float(RNG.uniform(.04, .92))])
        else:
            start = float(RNG.uniform(.04, .60))
            cutoff = np.array([start, float(RNG.uniform(start+.05, .95))])
        if analog:
            cutoff = cutoff * float(10**RNG.uniform(-1, 2))
        else:
            cutoff = cutoff * fs/2
        results.append(('F',order,mode,analog,rp,rs,fs,cutoff,np.array([]),np.array([])))

    for j in range(800):
        mode = j % 4 + 1
        analog = bool((j // 4) % 2)
        rp = float(RNG.uniform(.2, 4))
        rs = float(RNG.uniform(rp+14, 90))
        fs = 2.0 if analog else float(RNG.choice([2., 100., 1000.]))
        wp, ws = edges_for_order(mode)
        factor = float(10**RNG.uniform(-1, 2)) if analog else fs/2
        wp, ws = factor*wp, factor*ws
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            ref_n, _ = signal.ellipord(wp,ws,rp,rs,analog=analog,**({} if analog else {'fs':fs}))
        if ref_n <= 20:
            results.append(('O',0,mode,analog,rp,rs,fs,np.array([]),wp,ws))
            if ref_n <= 6 and j % 3 == 0:
                results.append(('D',0,mode,analog,rp,rs,fs,np.array([]),wp,ws))
    return results


def main():
    jobs = cases()
    rows = []
    for op, n, mode, analog, rp, rs, fs, cut, wp, ws in jobs:
        c = np.zeros(2); p = np.zeros(2); s = np.zeros(2)
        c[:len(cut)] = cut
        p[:len(wp)] = wp
        s[:len(ws)] = ws
        vals = [op,n,mode,int(analog),rp,rs,fs,*c,*p,*s]
        rows.append(' '.join(str(x) for x in vals))

    errs = {'F':0.,'O':0.,'D':0.}
    passed = {'F':0,'O':0,'D':0}
    with tempfile.TemporaryDirectory(prefix='scifort-elliptic-') as temp:
        work = Path(temp)
        flags = ['-std=f2018','-Wimplicit-interface','-Werror=implicit-interface',
                 '-O0','-fcheck=all']
        objects = []
        for name in MODULES:
            obj = work/f'{name}.o'
            subprocess.run(['gfortran',*flags,'-I',str(work),'-J',str(work),
                            '-c',str(ROOT/'src'/f'{name}.f90'),'-o',str(obj)],
                           check=True,cwd=work,capture_output=True)
            objects.append(str(obj))
        driver = work/'elliptic_driver'
        subprocess.run(['gfortran',*flags,'-I',str(work),
                        str(ROOT/'tools/signal_reference/elliptic_driver.f90'),
                        *objects,'-o',str(driver)],check=True,cwd=work,capture_output=True)
        proc = subprocess.run([str(driver)],input='\n'.join(rows)+'\n',text=True,
                              check=True,capture_output=True,cwd=work)
        outputs = iter(proc.stdout.splitlines())
        for i, (op, n, mode, analog, rp, rs, fs, cut, wp, ws) in enumerate(jobs):
            status, m, count = map(int,next(outputs).split())
            if status != 0:
                raise AssertionError(f'case {i} {op} returned status {status}')
            values = np.fromstring(next(outputs),sep=' ')
            expected_length = count if op == 'O' else m
            if len(values) != expected_length:
                raise AssertionError(f'case {i} {op} produced unexpected array length')
            kwargs = {'analog':analog}
            if not analog:
                kwargs['fs'] = fs
            if op == 'O':
                ref_n, ref_cut = signal.ellipord(wp,ws,rp,rs,**kwargs)
                ref_cut = np.atleast_1d(ref_cut)
                err = max(abs(values-ref_cut)/(1+abs(ref_cut)))
                # Difference arises from bounded minimization tolerance in
                # bandstop order selection; SciFort uses a tighter bound.
                tol = 5e-5 if mode == 4 else 2e-10
                if m != ref_n or err > tol:
                    raise AssertionError(f'case {i} order: got {m} {values}, expected {ref_n} {ref_cut}')
            else:
                a = np.fromstring(next(outputs),sep=' ')
                if len(a) != count:
                    raise AssertionError(f'case {i}: denominator length mismatch')
                args = (n,rp,rs,cut) if op == 'F' else (wp,ws,rp,rs)
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    if op == 'F':
                        b_ref,a_ref = signal.ellip(*args,btype=BANDS[mode-1],**kwargs)
                    else:
                        b_ref,a_ref = signal.iirdesign(*args,ftype='ellip',**kwargs)
                err = max(float(np.max(abs(values-b_ref)/(1+abs(b_ref)))),
                          float(np.max(abs(a-a_ref)/(1+abs(a_ref)))))
                tol = 2e-8 if op == 'F' else (3e-3 if mode == 4 else 1e-8)
                if err > tol:
                    raise AssertionError(f'case {i} design: error {err:.3e} > {tol:.3e}')
            errs[op] = max(errs[op],float(err))
            passed[op] += 1
        if next(outputs,None) is not None:
            raise AssertionError('unexpected extra output lines from Fortran')
    print(f'PASS: {sum(passed.values())} randomized SciPy elliptic comparisons')
    for op, name in (('F','filter designs'),('O','order selections'),('D','iirdesign calls')):
        print(f'  {passed[op]} {name}; maximum scaled discrepancy {errs[op]:.3e}')


if __name__ == '__main__':
    main()
