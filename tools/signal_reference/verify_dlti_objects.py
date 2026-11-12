#!/usr/bin/env python3
"""Validate Fortran discrete-LTI objects against SciPy frequency and state-space references.

Uses SciPy for the state-space realizations, discrete frequency responses,
and real-input dlsim; the complex-valued time-domain reference evaluates the
same state equation explicitly because SciPy 1.17 casts complex dlsim values.
"""
from __future__ import annotations

import argparse
import subprocess
import tempfile
import warnings
from pathlib import Path

import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[2]
SOURCES = (
    'kinds', 'signal_analog_design', 'signal_windows_extra', 'signal_dpss',
    'signal_spectral', 'signal_design', 'signal_analog',
    'signal_zpk_conversion', 'signal_state_space',
    'signal_discrete_simulation', 'signal_lti_frequency', 'signal_dlti',
)


def build(directory: Path, optimized: bool) -> Path:
    flags = ['-std=f2018', '-O3' if optimized else '-O0', '-Wimplicit-interface',
             '-Werror=implicit-interface']
    if not optimized:
        flags.extend(['-fcheck=all', '-g', '-fbacktrace'])
    objects = []
    for name in SOURCES:
        src = ROOT / 'src' / f'{name}.f90'
        obj = directory / f'{name}.o'
        subprocess.run(['gfortran', *flags, '-I', str(directory), '-J', str(directory),
                        '-c', str(src), '-o', str(obj)], check=True, cwd=directory)
        objects.append(str(obj))
    exe = directory / 'dlti_object_driver'
    subprocess.run(['gfortran', *flags, '-I', str(directory),
                    str(ROOT / 'tools/signal_reference/dlti_object_driver.f90'),
                    *objects, '-o', str(exe)], check=True, cwd=directory)
    return exe


def write_complex(lines: list[str], values: np.ndarray) -> None:
    for value in np.asarray(values, dtype=complex).ravel(order='F'):
        lines.append(f'({value.real:.17e},{value.imag:.17e})')


def recurrence(ss, u):
    a, b, c, d = (np.asarray(v, dtype=complex) for v in (ss.A, ss.B, ss.C, ss.D))
    n, ni = b.shape
    nt = len(u)
    y = np.empty((nt, c.shape[0]), dtype=complex)
    x = np.zeros(n, dtype=complex)
    for t in range(nt):
        y[t] = c @ x + d @ u[t]
        x = a @ x + b @ u[t]
    return y


def reference(c):
    w = c['w']
    sys = c['sys']
    n, ni = c['u'].shape
    no = c['no']
    if no == 1 and ni == 1:
        if c['mode'] == 2:
            # SciPy 1.17 freqz_zpk rejects complex gains in its array-api
            # casting path. Evaluate SciPy's zpk2tf polynomials instead.
            b_poly, a_poly = signal.zpk2tf(c['z'],c['p'],c['gain'])
            freq = np.polyval(b_poly,np.exp(1j*w))/np.polyval(a_poly,np.exp(1j*w))
        else:
            freq = sys.freqresp(w=w)[1]
        freq = freq.reshape(1, 1, -1)
    else:
        a, b, cc, d = (np.asarray(v, dtype=complex) for v in (sys.A, sys.B, sys.C, sys.D))
        eye = np.eye(a.shape[0], dtype=complex)
        freq = np.stack([cc @ np.linalg.solve(np.exp(1j*wi)*eye-a,b)+d
                         for wi in w], axis=2)
    magnitude = 20*np.log10(np.abs(freq))
    phase = np.unwrap(np.angle(freq), axis=2)*180/np.pi
    ss = sys.to_ss()
    y = recurrence(ss,c['u'])
    steps = np.empty((n,no,ni), dtype=complex)
    impulses = np.empty((n,no,ni), dtype=complex)
    for j in range(ni):
        ustep = np.zeros((n,ni),dtype=complex)
        ustep[:,j] = 1
        steps[:,:,j] = recurrence(ss,ustep)
        ustep[:] = 0
        ustep[0,j] = 1
        impulses[:,:,j] = recurrence(ss,ustep)
    return freq, magnitude, phase, y, steps, impulses


def create_cases():
    rng = np.random.default_rng(20261025)
    cases = []
    for mode in (1,2,3):
        for it in range(56):
            complex_case = it % 3 == 0
            dt = [0.05,0.2,0.5,1.0,2.0][it%5]
            nt = 12 + it % 7
            w = np.linspace(0.03,2.6,17) + 0.002*(it%3)
            # Complex 1-D systems and MIMO of sizes 1..3, 1..3.
            if mode < 3:
                order = 1 + it%4
                poles = rng.uniform(-.65,.65,order)
                zeros = rng.uniform(-.7,.7,it%(order+1))
                gain = float(rng.uniform(.1,2.0))
                if complex_case:
                    poles = poles + 1j*rng.uniform(-.15,.15,order)
                    zeros = zeros + 1j*rng.uniform(-.3,.3,len(zeros))
                    gain = gain + 0.15j
                b, a = signal.zpk2tf(zeros,poles,gain)
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    if mode == 1:
                        sys = signal.TransferFunction(b,a,dt=dt)
                    else:
                        if complex_case:
                            # Work around SciPy ZPK class dtype handling.
                            sys = signal.TransferFunction(b,a,dt=dt)
                        else:
                            sys = signal.ZerosPolesGain(zeros,poles,gain,dt=dt)
                ni = no = 1
                case = dict(mode=mode,dt=dt,nt=nt,w=w,ns=order,ni=ni,no=no,
                            b=np.asarray(b),a=np.asarray(a),z=np.asarray(zeros),
                            p=np.asarray(poles),gain=gain,sys=sys)
            else:
                ns = 1 + it%4
                ni = 1 + it%3
                no = 1 + (it//3)%3
                a = np.triu(rng.normal(size=(ns,ns))*0.08)
                np.fill_diagonal(a,rng.uniform(-.65,.65,ns))
                b = rng.normal(size=(ns,ni))*0.8
                cc = rng.normal(size=(no,ns))*0.8
                d = rng.normal(size=(no,ni))*0.15
                if complex_case:
                    a = a + 0.04j*rng.normal(size=(ns,ns))
                    b = b + 0.1j*rng.normal(size=(ns,ni))
                    cc = cc + 0.1j*rng.normal(size=(no,ns))
                    d = d + 0.1j*rng.normal(size=(no,ni))
                sys = signal.StateSpace(a,b,cc,d,dt=dt)
                case = dict(mode=mode,dt=dt,nt=nt,w=w,ns=ns,ni=ni,no=no,
                            aa=a,bb=b,cc=cc,dd=d,sys=sys)
            u = rng.normal(size=(nt,ni))
            if complex_case:
                u = u+0.12j*rng.normal(size=(nt,ni))
            case['u'] = u
            cases.append(case)
    return cases


def read_complex(lines: list[str], pos: int):
    re, im = [float(x) for x in lines[pos].split()]
    return re+1j*im, pos+1


def verify(executable: Path, cases) -> tuple[int,float,float]:
    contents = [str(len(cases))]
    for c in cases:
        nb = len(c.get('b',[])) if c['mode'] == 1 else 0
        na = len(c.get('a',[])) if c['mode'] == 1 else 0
        nz = len(c.get('z',[])) if c['mode'] == 2 else 0
        npoles = len(c.get('p',[])) if c['mode'] == 2 else 0
        contents.append(f"{c['mode']} {c['nt']} {len(c['w'])} {c['ns']} "
                        f"{c['ni']} {c['no']} {nb} {na} {nz} {npoles} {c['dt']}")
        contents.append(' '.join(f'{v:.17e}' for v in c['w']))
        if c['mode'] == 1:
            write_complex(contents,c['b'])
            write_complex(contents,c['a'])
        elif c['mode'] == 2:
            write_complex(contents,c['z'])
            write_complex(contents,c['p'])
            write_complex(contents,np.asarray([c['gain']]))
        else:
            for name in ('aa','bb','cc','dd'):
                write_complex(contents,c[name])
        write_complex(contents,c['u'])
    completed = subprocess.run([str(executable)], input='\n'.join(contents)+'\n',
                               capture_output=True, text=True, check=True)
    lines = completed.stdout.splitlines()
    pos = 0
    comparisons = 0
    max_diff = 0.0
    max_scaled = 0.0
    for idx,c in enumerate(cases,1):
        header = list(map(int,lines[pos].split()))
        pos += 1
        if header[:6] != [idx,0,0,0,0,0]:
            raise AssertionError(f'case {idx} status {header}')
        fr_status = list(map(int,lines[pos].split()))
        pos += 1
        if fr_status != [0,0]:
            raise AssertionError(f'case {idx} frequency status {fr_status}')
        fr_ref, mag_ref, phase_ref, y_ref, step_ref, impulse_ref = reference(c)
        n,ni,no = c['nt'],c['ni'],c['no']
        fr_actual = np.empty((no,ni,len(c['w'])),dtype=complex)
        mag_actual = np.empty_like(fr_actual.real)
        phase_actual = np.empty_like(fr_actual.real)
        for k in range(len(c['w'])):
            for j in range(ni):
                for i in range(no):
                    fr_actual[i,j,k],pos=read_complex(lines,pos)
                    mag_actual[i,j,k],phase_actual[i,j,k] = map(float,lines[pos].split())
                    pos+=1
        y_actual = np.empty((n,no),dtype=complex)
        step_actual = np.empty((n,no,ni),dtype=complex)
        impulse_actual = np.empty((n,no,ni),dtype=complex)
        for k in range(n):
            for i in range(no):
                y_actual[k,i],pos=read_complex(lines,pos)
        for k in range(n):
            for j in range(ni):
                for i in range(no):
                    step_actual[k,i,j],pos=read_complex(lines,pos)
                    impulse_actual[k,i,j],pos=read_complex(lines,pos)
        for label,actual,expected in (("frequency",fr_actual,fr_ref),
                                      ("bode magnitude",mag_actual,mag_ref),
                                      ("bode phase",phase_actual,phase_ref),
                                      ("output",y_actual,y_ref),
                                      ("step",step_actual,step_ref),
                                      ("impulse",impulse_actual,impulse_ref)):
            if not np.all(np.isfinite(expected)):
                raise AssertionError(f'{idx} {label}: nonfinite reference')
            if not np.all(np.isfinite(actual)):
                raise AssertionError(f'{idx} {label}: nonfinite actual')
            difference=np.max(np.abs(actual-expected))
            scaled=np.max(np.abs(actual-expected)/(1+np.abs(expected)))
            max_diff=max(max_diff,float(difference))
            max_scaled=max(max_scaled,float(scaled))
            comparisons += actual.size
            if scaled > (3e-8 if label=='bode phase' else 1e-8):
                raise AssertionError(f'{idx} mode {c["mode"]} {label}: diff={difference} scaled={scaled}')
    if pos!=len(lines):
        raise AssertionError(f'Output trailing lines: {len(lines)-pos}')
    return len(cases),max_diff,max_scaled,comparisons


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--optimized',action='store_true')
    args=parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='scifort-dlti-') as tmp:
        binary=build(Path(tmp),args.optimized)
        count,diff,scaled,values=verify(binary,create_cases())
        print(f'PASS: {count} SciPy LTI systems, {values} values; max_abs={diff:.4g}; '
              f'max_scaled={scaled:.4g} ({"-O3" if args.optimized else "-O0"})')


if __name__=='__main__':
    main()
