#!/usr/bin/env python3
"""Cross-check seven independently implemented signal algorithms against SciPy.

Run without fpm: python tools/signal_reference/compare_signal_missing_algorithms_scipy.py
Requirements: gfortran, NumPy and SciPy. The script builds its own small test
executable in a temporary folder and leaves the working tree unchanged.
"""
from __future__ import annotations

from pathlib import Path
import argparse
import subprocess
import tempfile

import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[2]
SOURCES = (
    'kinds', 'signal_windows_extra', 'signal_dpss', 'signal_spectral',
    'signal_analog_design', 'signal_design', 'signal_analog', 'signal_fir2d',
    'signal_stft_dual', 'signal_czt', 'signal_unique_roots',
    'signal_spatial', 'signal_spline', 'signal_spline2d', 'signal_symiir',
)


def parse_complex(parts):
    data = np.asarray(list(map(float, parts)), dtype=float)
    if len(data) % 2:
        raise AssertionError('Expected paired real and imaginary values')
    return data[::2] + 1j*data[1::2]


def check(actual, expected, description, tracker, tol=2e-10):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    if actual.shape != expected.shape:
        raise AssertionError(f'{description}: shape {actual.shape} != {expected.shape}')
    diff = float(np.max(abs(actual-expected), initial=0))
    rel = float(np.max(abs(actual-expected)/np.maximum(1,abs(expected)), initial=0))
    if not np.isfinite(rel) or rel > tol:
        raise AssertionError(f'{description}: max absolute error {diff:.3g}, relative {rel:.3g}')
    tracker[0] += 1
    tracker[1] += actual.size
    tracker[2] = max(tracker[2],diff)


def validate(line, tracker):
    parts = line.split()
    code = parts[0]
    if code == 'D':
        n,hop,complex_flag,desired_flag,scaled,status,count = map(int,parts[1:8])
        assert status == 0 and n == count
        alpha = complex(float(parts[8]),float(parts[9]))
        result = parse_complex(parts[10:])
        ix = np.arange(1,n+1,dtype=float)
        win = 0.1+np.sin(.17*ix)**2 + 1j*.2*np.cos(.19*ix)
        desired = .7+.25*np.cos(.27*ix) + 1j*.4*np.sin(.23*ix)
        if not complex_flag:
            win = win.real
            desired = desired.real
        if not desired_flag:
            desired = None
        ref, factor = signal.closest_STFT_dual_window(win,hop,desired_dual=desired,scaled=bool(scaled))
        check(result,ref,f'D(n={n},hop={hop},complex={complex_flag},scaled={scaled})',tracker)
        check([alpha],[factor],f'D alpha(n={n},hop={hop},complex={complex_flag})',tracker)
    elif code == 'F':
        r,c,ir,wi,mode,status,nr,nc = map(int,parts[1:9])
        assert status == 0
        arr = np.asarray(list(map(float,parts[9:])),dtype=float).reshape((nr,nc),order='F')
        names = ('hamming','hann','blackman',('kaiser',5.))
        window = names[wi-1]
        cutoff = .1 if mode == 0 else [.17,.32]
        if ir == 2:
            expected = signal.firwin_2d((r,c),window,fc=cutoff,circular=True)
        else:
            if mode == 0:
                expected = signal.firwin_2d((r,c),(window,window),fc=cutoff)
            else:
                # SciPy 1.17 firwin_2d accidentally ignores pass_zero in
                # separable mode. Verify the documented behavior against
                # the product of two correctly configured SciPy firwin calls.
                tap = signal.firwin(r,cutoff,window=window,pass_zero='bandpass')
                expected = np.outer(tap,signal.firwin(c,cutoff,window=window,
                                                    pass_zero='bandpass'))
        check(arr,expected,f'FIR2D {(r,c)} {window},radial={ir==2},mode={mode}',tracker)
    elif code == 'U':
        is_real,case,method,status,count = map(int,parts[1:6])
        assert status == 0
        values = np.asarray(list(map(float,parts[6:])),dtype=float)
        assert len(values) == 3*count
        actual = values[::3]+1j*values[1::3]
        counts = values[2::3].astype(int)
        ix = np.arange(1,9)
        roots = (ix//4)+.007*(ix%4)+1j*.007*((ix+1)%4)
        if is_real == 1:
            roots = roots.real
        tol = .002*(1+5*(case-1))
        how = ('min','max','avg')[method-1]
        ref, multiplicities = signal.unique_roots(roots,tol=tol,rtype=how)
        check(actual,ref,f'unique_roots {is_real,case,how}',tracker)
        if not np.array_equal(counts,multiplicities):
            raise AssertionError('unique_roots: multiplicities differ')
    elif code == 'W':
        case,status,n = map(int,parts[1:4])
        assert status == 0
        data = np.asarray(list(map(float,parts[4:])))
        assert len(data) == n
        args = {
            1: ([1,0],[1,8,25],13,'ba'),
            2: ([1,2],[1,2,3],17,'ba'),
            3: ([0],[-1+2j,-1-2j],11,'zp'),
            4: ([-2+1j],[-2+2j,-2-2j],12,'zp'),
            5: ([0],[-.1+0j],7,'zp'),
        }[case]
        ref = signal.findfreqs(args[0],args[1],args[2],kind=args[3])
        check(data,ref,f'findfreqs case {case}',tracker)
    elif code == 'S':
        family,case,complex_flag,rank,status,nr,nc = map(int,parts[1:8])
        assert status == 0
        actual = parse_complex(parts[8:]).reshape((nr,nc),order='F')
        ix = np.arange(1,121,dtype=float)
        x = np.sin(.1*ix)+.1*np.cos(.023*ix) + 1j*(.2*np.cos(.07*ix)-.08*np.sin(.03*ix))
        if rank == 2:
            x = np.vstack([x*(1+.25*k+1j*.03*k) for k in (1,2)])
        if not complex_flag:
            # SciPy 1.17 symiirorder1 mishandles strided real views of complex data.
            x = x.real.copy()
        if family == 1:
            p1,p2 = ((.64,.2),(.25,.5),(1.69,-.3))[case-1]
            expected = signal.symiirorder1(x,p1,p2,precision=1e-12)
        else:
            p1,p2 = ((.2,.4),(.5,.8),(.7,1.2))[case-1]
            if complex_flag:
                # SciPy 1.17 does not accept complex inputs to symiirorder2.
                # Use linearity to form the complex reference.
                expected = (signal.symiirorder2(x.real.copy(),p1,p2,precision=1e-12) +
                            1j*signal.symiirorder2(x.imag.copy(),p1,p2,precision=1e-12))
            else:
                expected = signal.symiirorder2(x,p1,p2,precision=1e-12)
        expected = np.asarray(expected).reshape((nr,nc))
        check(actual,expected,f'symiir family={family} case={case} complex={complex_flag} rank={rank}',
              tracker, tol=3e-10)
    elif code == 'C':
        n,non_default,status,count = map(int,parts[1:5])
        assert status == 0 and count == n
        actual = parse_complex(parts[5:])
        expected = (signal.czt_points(n,w=np.exp(-.19j),a=1.1+.4j) if non_default else
                    signal.czt_points(n))
        check(actual,expected,f'czt_points({n},non_default={non_default})',tracker)
    else:
        raise AssertionError(f'Unexpected Fortran record: {line[:50]}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--optimized',action='store_true',help='compile kernels with -O3')
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='scifort-signal-missing-') as folder:
        work = Path(folder)
        flags = ['-std=f2018','-O3' if args.optimized else '-O0',
                 '-Werror=implicit-interface','-ffree-line-length-none','-I',str(work),'-J',str(work)]
        if not args.optimized:
            flags.append('-fcheck=all')
        objects = []
        for name in SOURCES:
            path = work/(name+'.o')
            subprocess.run(['gfortran',*flags,'-c',str(ROOT/'src'/(name+'.f90')),
                            '-o',str(path)],cwd=work,check=True)
            objects.append(str(path))
        exe = work/'missing_algorithms_driver'
        subprocess.run(['gfortran',*flags,str(ROOT/'tools/signal_reference'/
                        'missing_algorithms_driver.f90'),*objects,'-o',str(exe)],
                       cwd=work,check=True)
        records = subprocess.check_output([str(exe)],text=True).splitlines()
        tracker = [0,0,0.]
        for row in records:
            validate(row,tracker)
    print(f'PASS: {tracker[0]} SciPy signal algorithm comparisons, '
          f'{tracker[1]} values (max absolute error {tracker[2]:.3e})')


if __name__ == '__main__':
    main()
