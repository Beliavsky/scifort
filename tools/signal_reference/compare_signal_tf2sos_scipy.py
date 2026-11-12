#!/usr/bin/env python3
"""Compare SciFort tf2sos with SciPy 1.17 for digital and analog filters.

Requires Python, NumPy, SciPy and gfortran; no FPM installation needed.
Sections may have different order from SciPy, so compare frequency response.
Analog and minimal-pairing rows are polynomials in positive powers.
"""
from __future__ import annotations

import subprocess
import tempfile
import warnings
from pathlib import Path

import numpy as np
from scipy import signal
from run_signal_suite import ROOT, MODULES


PROBE = r'''program tf2sos_probe
    use scifort_signal, only : dp, tf2sos, filter_sos_result
    implicit none
    real(dp), allocatable :: b(:), a(:)
    type(filter_sos_result) :: ss
    integer :: ncases, i, j, nb, na, analog_flag, pairing_flag
    read (*,*) ncases
    do i = 1, ncases
        read (*,*) nb, na, analog_flag, pairing_flag
        allocate(b(nb),a(na))
        do j = 1, nb
            read (*,*) b(j)
        end do
        do j = 1, na
            read (*,*) a(j)
        end do
        if (analog_flag == 1) then
            ss = tf2sos(b,a,analog=.true.)
        else
            select case (pairing_flag)
            case (0)
                ss = tf2sos(b,a)
            case (1)
                ss = tf2sos(b,a,pairing='keep_odd')
            case (2)
                ss = tf2sos(b,a,pairing='minimal')
            end select
        end if
        if (ss%status == 0) then
            write (*,'(2(i0,1x))') ss%status,size(ss%sos,1)
            do j = 1, size(ss%sos,1)
                write (*,'(12(es25.16e3,1x))') &
                    real(ss%sos(j,:),dp),aimag(ss%sos(j,:))
            end do
        else
            write (*,'(2(i0,1x))') ss%status,0
        end if
        deallocate(a,b)
    end do
end program tf2sos_probe
'''


def make_roots(rng, order, kind, mode):
    values = []
    for k in range(order // 2):
        if kind == 'analog':
            real = rng.uniform(-2.0, -.15)
            imag = rng.uniform(.1, 2.0)
        else:
            real = rng.uniform(-.8, .8)
            imag = rng.uniform(.1, .65)
        if mode % 4 == 0:
            if kind == 'analog':
                values.extend([real, rng.uniform(-2.0, -.15)])
            else:
                values.extend([real, rng.uniform(-.8, .8)])
        else:
            values.extend([real + 1j * imag, real - 1j * imag])
    if order % 2:
        if kind == 'analog':
            values.append(rng.uniform(-2., -.15))
        else:
            values.append(rng.uniform(-.8, .8))
    if len(values) and mode % 19 == 0 and kind == 'digital':
        # Add real roots at the origin without breaking conjugacy.
        values = np.asarray(values, dtype=complex)
        if order % 2:
            values[-1] = 0.
        elif mode % 4 == 0:
            values[:2] = 0.
        return values
    return np.asarray(values,complex)


def cases():
    rng = np.random.default_rng(179013)
    result = []
    for kind in ('digital', 'analog'):
        for order in range(11):
            for i in range(48):
                if kind == 'analog':
                    degree_b = (i * 7 + order) % (order + 1)
                else:
                    degree_b = (i * 11 + order) % (order + 4)
                zeros = make_roots(rng, degree_b,kind,i+3)
                poles = make_roots(rng, order,kind,i)
                b = np.atleast_1d((np.poly(zeros) * rng.uniform(.08, 3.) * (-1 if i % 9 == 0 else 1)).real)
                a = np.atleast_1d((np.poly(poles) * rng.uniform(.3, 2.5)).real)
                if i % 31 == 0:
                    b = np.r_[0.,b]
                pairing = 0 if kind == 'analog' else (i % 3)
                result.append((b,a,kind,pairing))
    return result


def response(rows, frequencies, positive_power):
    q = (1j * frequencies if positive_power == 'analog'
         else np.exp(1j*frequencies))
    if positive_power is False:
        q = 1 / q
    h = np.ones_like(q, dtype=complex)
    for row in rows:
        if positive_power is False:
            num = row[0]+row[1]*q+row[2]*q*q
            den = row[3]+row[4]*q+row[5]*q*q
        else:
            num = np.polyval(row[:3], q)
            den = np.polyval(row[3:], q)
        h *= num / den
    return h


def run():
    examples = cases()
    lines = [str(len(examples))]
    for b, a, kind, pairing in examples:
        lines.append(f'{len(b)} {len(a)} {int(kind == "analog")} {pairing}')
        lines.extend(format(v,'.17e') for v in b)
        lines.extend(format(v,'.17e') for v in a)
    with tempfile.TemporaryDirectory(prefix='scifort-tf2sos-') as temp:
        temp = Path(temp)
        flags = ['-std=f2018','-O0','-fcheck=all','-Wimplicit-interface',
                 '-Werror=implicit-interface']
        objects = []
        for name in MODULES:
            obj = temp / f'{name}.o'
            subprocess.run(['gfortran',*flags,'-J',str(temp),'-I',str(temp),
                '-c',str(ROOT/'src'/f'{name}.f90'),'-o',str(obj)],check=True,cwd=temp)
            objects.append(str(obj))
        source = temp / 'probe.f90'
        source.write_text(PROBE)
        executable = temp / 'probe'
        subprocess.run(['gfortran',*flags,'-I',str(temp),'-J',str(temp),
            str(source),*objects,'-o',str(executable)],check=True,cwd=temp)
        out = subprocess.run([str(executable)],input='\n'.join(lines)+'\n',
            capture_output=True,check=True,text=True,cwd=temp)
    output = out.stdout.splitlines()
    pos = 0
    max_scaled = 0.
    max_polynomial_scaled = 0.
    count = 0
    rejected = 0
    scipy_failed = 0
    for idx,(b,a,kind,pairing) in enumerate(examples):
        status, ns = map(int,output[pos].split()); pos += 1
        mode = {0:None,1:'keep_odd',2:'minimal'}[pairing]
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', signal.BadCoefficients)
                ref = signal.tf2sos(b,a,pairing=mode,analog=(kind=='analog'))
        except ValueError as exc:
            if 'must have len(p)>=len(z)' in str(exc):
                if status == 0:
                    raise AssertionError(f'case {idx}: accepted an improper minimal design')
                rejected += 1
                continue
            ref = None
            scipy_failed += 1
        except IndexError:
            ref = None
            scipy_failed += 1
        if status:
            raise AssertionError(f'case {idx}: status {status} for kind={kind}, pairing={pairing}, b={b}, a={a}')
        rows = []
        for section in range(ns):
            values = np.array(list(map(float,output[pos].split())))
            if values.size != 12:
                raise AssertionError(f'case {idx} row incomplete')
            rows.append(values[:6]+1j*values[6:]); pos += 1
        ours = np.asarray(rows)
        freq = (np.logspace(-2,1.5,31) if kind=='analog'
                else np.linspace(.025,np.pi-.025,31))
        positive_power = ('analog' if kind=='analog' else pairing==2)
        actual = response(ours,freq,positive_power)
        q = 1j*freq if kind=='analog' else np.exp(1j*freq)
        if kind=='digital' and pairing != 2:
            q = 1/q
            bb = np.trim_zeros(b,'f')
            aa = np.trim_zeros(a,'f')
            if len(bb) == 0:
                bb = np.array([0.])
            polynomial = sum(v*q**j for j,v in enumerate(bb)) / \
                sum(v*q**j for j,v in enumerate(aa))
        else:
            polynomial = np.polyval(b,q)/np.polyval(a,q)
        direct_scaled = abs(actual-polynomial) / \
            (1.+np.maximum(abs(actual),abs(polynomial)))
        max_polynomial_scaled = max(max_polynomial_scaled,np.max(direct_scaled))
        if np.max(direct_scaled) > 5.e-6:
            raise AssertionError(f'case {idx}: polynomial mismatch {np.max(direct_scaled)}')
        expected = response(ref,freq,positive_power) if ref is not None else polynomial
        scaled = abs(actual-expected)/(1.+np.maximum(abs(actual),abs(expected)))
        discrepancy = np.max(scaled)
        max_scaled = max(max_scaled,discrepancy)
        if not np.isfinite(discrepancy) or discrepancy > 5.e-6:
            raise AssertionError(f'case {idx}: kind {kind} pairing {pairing}, err {discrepancy}; b={b}, a={a}, ours={ours}, ref={ref}')
        count += 1
    if pos != len(output):
        raise AssertionError('unexpected output')
    print(f'PASS: {count} randomized tf2sos analog/digital cases, {rejected} rejected improper requests; '
          f'{scipy_failed} SciPy reference failures checked against polynomial response; '
          f'max scaled response error {max_scaled:.3e}; max vs original TF {max_polynomial_scaled:.3e}')

if __name__=='__main__':
    run()
