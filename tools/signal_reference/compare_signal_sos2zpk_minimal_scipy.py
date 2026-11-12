#!/usr/bin/env python3
"""Check SciFort SOS/ZPK conversions against SciPy 1.17 and rational responses.

Compile the Fortran modules with gfortran (no FPM required). In minimal mode
SOS coefficients represent polynomials in positive powers of z and may have
a0=0, so do NOT pass those sections to digital sosfreqz/sosfilt.
"""
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy import signal
from scipy.optimize import linear_sum_assignment

from run_signal_suite import ROOT, MODULES


def roots(rng, n, mode):
    parts = []
    if mode % 3 == 0:
        return np.asarray(rng.uniform(-0.75, 0.75, n), complex)
    for _ in range(n // 2):
        re = rng.uniform(-0.7, 0.7)
        im = rng.uniform(0.08, 0.6)
        if mode % 3 == 1:
            im = 0.0
        parts.extend([re + 1j * im, re - 1j * im])
    if n % 2:
        parts.append(rng.uniform(-0.8, 0.8))
    return np.asarray(parts, complex)


def make_cases():
    rng = np.random.default_rng(173418)
    cases = []
    for npoles in range(11):
        for i in range(72):
            nz = (i * 5 + npoles) % (npoles + 1)
            z = roots(rng, nz, i)
            p = roots(rng, npoles, i + 1)
            # Include true roots at the origin so their allocation cannot
            # simply be assumed to be artificial padding.
            if i % 13 == 0 and nz:
                z[-1] = 0.0
            if i % 17 == 0 and npoles:
                p[-1] = 0.0
            # Keep a conjugate pair intact if replacing one complex root.
            if i % 13 == 0 and nz > 1 and z[-2] == np.conj(z[-1]):
                z[-2] = 0.0
            if i % 17 == 0 and npoles > 1 and p[-2] == np.conj(p[-1]):
                p[-2] = 0.0
            k = (0.1 + rng.random() * 1.8) * (-1 if i % 8 == 0 else 1)
            cases.append((z, p, complex(k)))
    sos_cases = []
    for i, (z, p, k) in enumerate(cases):
        if i % 5 == 0:
            pairing = 'minimal'
        elif i % 5 == 1:
            pairing = 'keep_odd'
        else:
            pairing = 'nearest'
        try:
            sos = np.asarray(signal.zpk2sos(z, p, k.real, pairing=pairing), complex)
        except (ValueError, IndexError):
            continue
        # Scale numerator and denominator of a section separately. This
        # probes the section-wise leading-coefficient normalization.
        if i % 11 == 0:
            sos[0, :3] *= 0.8 + 0.2j
            sos[0, 3:] *= 1.3 - 0.1j
        sos_cases.append(sos)
    return cases, sos_cases


def vec(values):
    return [f'{v.real:.17e} {v.imag:.17e}' for v in np.ravel(values)]


def probe_source():
    return r'''program signal_minimal_probe
    use scifort_signal, only : dp, zpk2sos, sos2zpk, &
        filter_sos_result, filter_zpk_result
    implicit none
    integer :: ncases, nsos, i, j, nz, np, ns
    real(dp) :: re, im
    complex(dp), allocatable :: z(:), p(:), sections(:,:)
    complex(dp) :: k
    type(filter_sos_result) :: ss
    type(filter_zpk_result) :: zp
    read (*,*) ncases, nsos
    do i = 1, ncases
        read (*,*) nz, np
        allocate(z(nz),p(np))
        read (*,*) re, im
        k = cmplx(re,im,dp)
        do j = 1, nz
            read (*,*) re, im
            z(j) = cmplx(re,im,dp)
        end do
        do j = 1, np
            read (*,*) re, im
            p(j) = cmplx(re,im,dp)
        end do
        ss = zpk2sos(z,p,k,pairing='minimal')
        write (*,'(2(i0,1x))') ss%status, size(ss%sos,1)
        if (ss%status == 0) then
            do j = 1, size(ss%sos,1)
                call print_complex(ss%sos(j,:))
            end do
        end if
        deallocate(z,p)
    end do
    do i = 1, nsos
        read (*,*) ns
        allocate(sections(ns,6))
        do j = 1, ns
            call read_complex(sections(j,:))
        end do
        zp = sos2zpk(sections)
        write (*,'(2(i0,1x))') zp%status, size(zp%zeros)
        if (zp%status == 0) then
            write (*,'(2(es25.16e3,1x))') real(zp%gain,dp), aimag(zp%gain)
            call print_complex(zp%zeros)
            call print_complex(zp%poles)
        end if
        deallocate(sections)
    end do
contains
    subroutine read_complex(v)
        complex(dp), intent(out) :: v(:) !! coefficients to read
        integer :: t
        do t = 1, size(v)
            read (*,*) re, im
            v(t) = cmplx(re,im,dp)
        end do
    end subroutine
    subroutine print_complex(v)
        complex(dp), intent(in) :: v(:) !! coefficients to output
        integer :: t
        do t = 1, size(v)
            write (*,'(2(es25.16e3,1x))') real(v(t),dp), aimag(v(t))
        end do
    end subroutine
end program signal_minimal_probe
'''


def compare_roots(actual, expected, label):
    if len(actual) != len(expected):
        raise AssertionError(f'{label}: {len(actual)} != {len(expected)} roots')
    if not len(actual):
        return
    costs = abs(actual[:, None] - expected[None, :])
    rows, cols = linear_sum_assignment(costs)
    err = max(costs[rows, cols])
    if err > 3e-7:
        raise AssertionError(f'{label}: roots differ by {err}')


def main():
    cases, sos_cases = make_cases()
    lines = [f'{len(cases)} {len(sos_cases)}']
    for z, p, k in cases:
        lines.append(f'{len(z)} {len(p)}')
        lines.extend(vec([k]))
        lines.extend(vec(z))
        lines.extend(vec(p))
    for sos in sos_cases:
        lines.append(str(len(sos)))
        for row in sos:
            lines.extend(vec(row))
    with tempfile.TemporaryDirectory(prefix='scifort-sos-conversion-') as temp:
        work = Path(temp)
        objects = []
        flags = ['-std=f2018', '-O0', '-fcheck=all', '-Wimplicit-interface',
                 '-Werror=implicit-interface']
        for module in MODULES:
            obj = work / f'{module}.o'
            subprocess.run(['gfortran', *flags, '-J', str(work), '-I', str(work),
                '-c', str(ROOT / 'src' / f'{module}.f90'), '-o', str(obj)],
                cwd=work, check=True)
            objects.append(str(obj))
        source = work / 'probe.f90'
        source.write_text(probe_source())
        exe = work / 'probe'
        subprocess.run(['gfortran', *flags, '-I', str(work), '-J', str(work),
            str(source), *objects, '-o', str(exe)], cwd=work, check=True)
        result = subprocess.run([str(exe)], input='\n'.join(lines) + '\n',
            text=True, capture_output=True, cwd=work, check=True)
    records = result.stdout.splitlines()
    pos = 0

    def fetch(n):
        nonlocal pos
        data = np.array([complex(*map(float, records[pos + j].split()))
                         for j in range(n)])
        pos += n
        return data

    # The minimal sections must represent exactly the given positive-power
    # ZPK transfer function, with no changes in degree.
    for i, (z, p, k) in enumerate(cases):
        status, ns = map(int, records[pos].split())
        pos += 1
        if status:
            raise AssertionError(f'minimal {i}: status {status}, roots {z}, {p}')
        sos = fetch(ns * 6).reshape(ns, 6)
        if ns != max(1, (len(p) + 1) // 2):
            raise AssertionError(f'minimal {i}: unexpected sections {ns}')
        frequencies = [0.03, 0.15, 0.4, 0.9, 1.8, 2.8]
        for w in frequencies:
            q = np.exp(1j * w)
            h = 1.0 + 0j
            for row in sos:
                h *= np.polyval(row[:3], q) / np.polyval(row[3:], q)
            href = k * np.prod(q - z) / np.prod(q - p)
            if not np.allclose(h, href, rtol=3e-10, atol=3e-10):
                raise AssertionError(f'minimal {i}: response {h} != {href} at {w}')
        # Independent SciPy factorization confirms the no-spurious-root
        # property, after removing its documented origin padding.
        scipy_zz, scipy_pp, scipy_gain = signal.sos2zpk(sos)
        expect_n = 2 * ns
        compare_roots(scipy_zz, np.r_[z, np.zeros(expect_n - len(z))], f'minimal zeros {i}')
        compare_roots(scipy_pp, np.r_[p, np.zeros(expect_n - len(p))], f'minimal poles {i}')
        if not np.allclose(scipy_gain, k, rtol=2e-10, atol=2e-10):
            raise AssertionError(f'minimal {i}: gain {scipy_gain} != {k}')

    for i, sos in enumerate(sos_cases):
        status, n = map(int, records[pos].split())
        pos += 1
        if status:
            raise AssertionError(f'sos2zpk {i}: status {status}')
        gain = fetch(1)[0]
        zz = fetch(n)
        pp = fetch(n)
        ez, ep, ek = signal.sos2zpk(sos)
        compare_roots(zz, ez, f'zeros {i}')
        compare_roots(pp, ep, f'poles {i}')
        if not np.allclose(gain, ek, rtol=3e-10, atol=3e-10):
            raise AssertionError(f'sos2zpk {i}: gain {gain} != {ek}')
    if pos != len(records):
        raise AssertionError(f'{len(records)-pos} unused output lines')
    print(f'PASS: {len(cases)} minimal ZPK designs + {len(sos_cases)} SOS inverses vs SciPy 1.17')


if __name__ == '__main__':
    main()
