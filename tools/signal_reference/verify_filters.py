#!/usr/bin/env python3
"""Optional randomized SciPy 1D signal regression (no SciPy runtime dependency).

Compiles a test-only Fortran driver in a temporary directory; compares generated
outputs with the public scipy.signal API. Requires Python, NumPy, SciPy, gfortran.
"""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[2]
DRIVER = r"""
program compare_signal_filters
    use scifort_signal
    implicit none
    character(len=1) :: tag
    character(len=20) :: mode
    real(dp), allocatable :: x(:), a(:), b(:), y(:), zi(:), zf(:)
    complex(dp), allocatable :: xc(:), ac(:), bc(:), yc(:), zic(:), zfc(:)
    real(dp) :: delta, cval
    integer :: ios, w, p, d, pos, n, na, nb, nz, status
    do
        read(*, *, iostat=ios) tag
        if (ios /= 0) exit
        select case (tag)
        case ('G')
            read(*, *) w, p, d, delta, pos, mode
            if (pos == 0) then
                y = savgol_coeffs(w, p, deriv=d, delta=delta, use=trim(mode), status=status)
            else
                y = savgol_coeffs(w, p, deriv=d, delta=delta, pos=pos, &
                    use=trim(mode), status=status)
            end if
            call emit_real(status, y)
        case ('S')
            read(*, *) n, w, p, d, delta, cval, mode
            allocate(x(n))
            if (n > 0) read(*, *) x
            y = savgol_filter(x, w, p, deriv=d, delta=delta, cval=cval, &
                mode=trim(mode), status=status)
            call emit_real(status, y)
            deallocate(x)
        case ('R')
            read(*, *) nb, na, n
            allocate(b(nb), a(na), x(n), zi(max(nb,na)-1))
            read(*, *) b
            read(*, *) a
            if (n > 0) read(*, *) x
            if (size(zi) > 0) read(*, *) zi
            y = lfilter(b, a, x, zi=zi, zf=zf, status=status)
            call emit_real(status, [y, zf])
            deallocate(b, a, x, zi)
        case ('C')
            read(*, *) nb, na, n
            allocate(bc(nb), ac(na), xc(n), zic(max(nb,na)-1))
            read(*, *) bc
            read(*, *) ac
            if (n > 0) read(*, *) xc
            if (size(zic) > 0) read(*, *) zic
            yc = lfilter(bc, ac, xc, zi=zic, zf=zfc, status=status)
            call emit_complex(status, [yc, zfc])
            deallocate(bc, ac, xc, zic)
        case ('Z')
            read(*, *) nb, na
            allocate(b(nb), a(na))
            read(*, *) b
            read(*, *) a
            y = lfilter_zi(b, a, status=status)
            call emit_real(status, y)
            deallocate(b, a)
        case ('Y')
            read(*, *) nb, na
            allocate(bc(nb), ac(na))
            read(*, *) bc
            read(*, *) ac
            yc = lfilter_zi(bc, ac, status=status)
            call emit_complex(status, yc)
            deallocate(bc, ac)
        case default
            error stop 'unknown test command'
        end select
    end do
contains
    subroutine emit_real(code, values)
        integer, intent(in) :: code !! call status
        real(dp), intent(in) :: values(:) !! returned samples
        integer :: j
        write(*, '(i0,1x,i0)', advance='no') code, size(values)
        do j = 1, size(values)
            write(*, '(1x,es26.17e3)', advance='no') values(j)
        end do
        write(*, *)
    end subroutine
    subroutine emit_complex(code, values)
        integer, intent(in) :: code !! call status
        complex(dp), intent(in) :: values(:) !! returned samples
        integer :: j
        write(*, '(i0,1x,i0)', advance='no') code, size(values)
        do j = 1, size(values)
            write(*, '(1x,es26.17e3,1x,es26.17e3)', advance='no') &
                real(values(j), dp), aimag(values(j))
        end do
        write(*, *)
    end subroutine
end program compare_signal_filters
"""


def formatted(array: np.ndarray) -> str:
    if np.iscomplexobj(array):
        return " ".join(f"({z.real:.18e},{z.imag:.18e})" for z in array)
    return " ".join(f"{v:.18e}" for v in array)


def main() -> None:
    rng = np.random.default_rng(20261014)
    cases = []
    tests = []

    def append(tag: str, params: str, arrays: tuple[np.ndarray, ...], expected: np.ndarray):
        record = [tag, params]
        record += [formatted(v) for v in arrays if v.size]
        cases.append("\n".join(record))
        tests.append((tag, np.asarray(expected)))

    # Both coefficient orders and arbitrary (one-based in Fortran) positions.
    for _ in range(440):
        w = int(rng.integers(2, 23))
        p = int(rng.integers(0, min(8, w)))
        d = int(rng.integers(0, min(5, p + 3)))
        delta = float(rng.uniform(0.6, 3.0))
        pos = int(rng.integers(1, w + 1)) if rng.random() < 0.5 else 0
        mode = str(rng.choice(['conv', 'dot']))
        kwargs = dict(deriv=d, delta=delta, use=mode)
        if pos:
            kwargs['pos'] = pos - 1
        expected = signal.savgol_coeffs(w, p, **kwargs)
        append('G', f'{w} {p} {d} {delta:.18e} {pos} {mode}', (), expected)

    # All boundary modes, including even windows and short inputs where valid.
    for _ in range(600):
        n = int(rng.integers(1, 90))
        w = int(rng.integers(1, min(19, n + 9)))
        p = int(rng.integers(0, min(8, w)))
        d = int(rng.integers(0, min(4, p + 2)))
        delta = float(rng.uniform(0.5, 2.5))
        mode = str(rng.choice(['interp', 'nearest', 'mirror', 'wrap', 'constant']))
        if mode == 'interp' and w > n:
            mode = 'mirror'
        cval = float(rng.normal())
        x = rng.normal(size=n)
        expected = signal.savgol_filter(x, w, p, deriv=d, delta=delta, mode=mode, cval=cval)
        append('S', f'{n} {w} {p} {d} {delta:.18e} {cval:.18e} {mode}', (x,), expected)

    for complex_mode in [False, True]:
        for _ in range(520):
            na = int(rng.integers(1, 8))
            nb = int(rng.integers(1, 8))
            n = int(rng.integers(1, 65))
            a = rng.normal(size=na) * 0.12
            a[0] = float(rng.choice([-2., -1., 0.7, 1.5, 3.]))
            b = rng.normal(size=nb) * 0.6
            x = rng.normal(size=n)
            initial = rng.normal(size=max(na, nb) - 1) * 0.2
            if complex_mode:
                a = a + 1j * rng.normal(size=na) * 0.05
                b = b + 1j * rng.normal(size=nb) * 0.2
                x = x + 1j * rng.normal(size=n)
                initial = initial + 1j * rng.normal(size=initial.size) * 0.2
            expected, zf = signal.lfilter(b, a, x, zi=initial)
            tag = 'C' if complex_mode else 'R'
            append(tag, f'{nb} {na} {n}', (b, a, x, initial), np.r_[expected, zf])

        for _ in range(150):
            na = int(rng.integers(2, 7))
            nb = int(rng.integers(1, 7))
            a = rng.normal(size=na) * 0.07
            a[0] = float(rng.choice([0.7, 1., 2.]))
            b = rng.normal(size=nb) * 0.4
            if complex_mode:
                a = a + 1j * rng.normal(size=na) * 0.04
                b = b + 1j * rng.normal(size=nb) * 0.3
            expected = signal.lfilter_zi(b, a)
            append('Y' if complex_mode else 'Z', f'{nb} {na}', (b, a), expected)

    sources = ['kinds.f90', 'signal_convolve.f90', 'signal_convolve2d.f90', 'signal_detrend.f90',
               'signal_peaks.f90', 'signal_savgol.f90', 'signal_lfilter.f90', 'signal_sos.f90', 'signal_filtfilt.f90', 'signal_filtfilt_nd.f90', 'signal_dpss.f90', 'signal_windows_extra.f90', 'signal_spectral.f90', 'signal_fir_advanced.f90', 'signal_spectral_complex.f90', 'signal_fftconvolve.f90', 'signal_analog_design.f90', 'signal_design.f90', 'signal_elliptic.f90', 'signal_analog_sos.f90', 'signal_iir_extra.f90', 'signal_fourier_resample.f90', 'signal_resample.f90', 'signal_analog.f90', 'signal_zpk_conversion.f90', 'signal_state_space.f90', 'signal_cont2discrete.f90', 'signal_discrete_simulation.f90', 'signal_continuous_simulation.f90', 'signal_lti_frequency.f90', 'signal_partial_fraction.f90', 'signal_sequences.f90', 'signal_hilbert.f90', 'signal_waveforms.f90', 'signal_periodic.f90', 'signal_irregular.f90', 'signal_remez.f90', 'signal_czt.f90', 'signal.f90']
    with tempfile.TemporaryDirectory(prefix='scifort-signal-') as temp:
        tmp = Path(temp)
        driver = tmp / 'driver.f90'
        driver.write_text(DRIVER, encoding='ascii')
        cmd = ['gfortran', '-std=f2018', '-O2', '-J', str(tmp), '-I', str(tmp)]
        cmd.extend(str(ROOT / 'src' / source) for source in sources)
        cmd.extend([str(driver), '-o', str(tmp / 'driver')])
        subprocess.run(cmd, cwd=tmp, check=True)
        result = subprocess.run([str(tmp / 'driver')], input='\n'.join(cases) + '\n',
                                text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                cwd=tmp, check=True)
    lines = result.stdout.splitlines()
    if len(lines) != len(tests):
        raise AssertionError(f'expected {len(tests)} results, got {len(lines)}; {result.stderr}')

    per_tag = {}
    maximum = {}
    for case_idx, (line, (tag, ref)) in enumerate(zip(lines, tests)):
        raw = line.split()
        status, n = map(int, raw[:2])
        dtype = complex if tag in ('C', 'Y') else float
        values = np.asarray([float(v) for v in raw[2:]], dtype=float)
        if dtype == complex:
            if len(values) != 2 * n:
                raise AssertionError(f'case {case_idx}: complex length error')
            values = values[::2] + 1j * values[1::2]
        if status != 0 or len(values) != n or len(ref) != n:
            raise AssertionError(f'case {case_idx} ({tag}): invalid Fortran result')
        error = np.max(np.abs(values - ref)) if len(ref) else 0.0
        if not np.allclose(values, ref, rtol=2e-7, atol=2e-8):
            raise AssertionError(f'case {case_idx} ({tag}) failed: max absolute difference={error:.4g}')
        per_tag[tag] = per_tag.get(tag, 0) + 1
        maximum[tag] = max(maximum.get(tag, 0.0), float(error))
    print(f'{len(tests)} SciPy reference cases passed (SciPy {signal.__version__ if hasattr(signal,"__version__") else "installed"}):')
    for tag in 'GSRCZY':
        print(f'  {tag}: {per_tag.get(tag, 0)} cases; worst abs error {maximum.get(tag, 0):.3g}')


if __name__ == '__main__':
    main()
