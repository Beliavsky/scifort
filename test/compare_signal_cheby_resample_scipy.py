#!/usr/bin/env python3
"""Optional numerical cross-check of three signal algorithms against SciPy.

Run: python test/compare_signal_cheby_resample_scipy.py
Requires a GNU Fortran compiler, numpy, and scipy. No sources are downloaded.
"""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[1]
FORTRAN = os.environ.get("FC", "gfortran")
DESIGN_DRIVER = r"""
program design_driver
    use scifort_signal, only : dp, cheby1, butter_result
    implicit none
    type(butter_result) :: r
    integer :: order, kind
    real(dp) :: rp
    character(8) :: names(4) = ['lowpass ', 'highpass', 'bandpass', 'bandstop']
    read(*,*) order,kind,rp
    if (kind <= 2) then
        r=cheby1(order,rp,0.32_dp,btype=trim(names(kind)))
    else
        r=cheby1(order,rp,[0.21_dp,0.64_dp],btype=trim(names(kind)))
    end if
    print *,r%status
    if (r%status/=0) stop
    print *,size(r%b)
    write(*,'(*(es25.16,1x))') r%b
    write(*,'(*(es25.16,1x))') r%a
end program
"""
RESAMPLE_DRIVER = r"""
program resample_driver
    use scifort_signal, only : dp, resample
    implicit none
    integer :: n,m,kind,win,i,st
    real(dp) :: ar,ai
    complex(dp), allocatable :: x(:),y(:)
    real(dp), allocatable :: z(:)
    read(*,*) n,m,kind,win
    allocate(x(n))
    do i=1,n
        read(*,*) ar,ai
        x(i)=cmplx(ar,ai,dp)
    end do
    if (kind==1) then
        if (win==1) then
            z=resample(real(x,dp),m,window='hann',status=st)
        else
            z=resample(real(x,dp),m,status=st)
        end if
        y=cmplx(z,0.0_dp,dp)
    else
        if (win==1) then
            y=resample(x,m,window='hann',status=st)
        else
            y=resample(x,m,status=st)
        end if
    end if
    print *,st,size(y)
    do i=1,size(y)
        write(*,'(2(es25.16,1x))') real(y(i),dp),aimag(y(i))
    end do
end program
"""
DECIMATE_DRIVER = r"""
program decimate_driver
    use scifort_signal, only : dp, decimate
    implicit none
    integer :: n,q,order,phase,i,st
    real(dp), allocatable :: x(:),y(:)
    read(*,*) n,q,order,phase
    allocate(x(n))
    do i=1,n
        read(*,*) x(i)
    end do
    if (order==0) then
        y=decimate(x,q,zero_phase=phase==1,status=st)
    else
        y=decimate(x,q,n=order,zero_phase=phase==1,status=st)
    end if
    print *,st,size(y)
    write(*,'(*(es25.16,1x))') y
end program
"""


def run(exe: Path, request: str) -> list[str]:
    result = subprocess.run([str(exe)], input=request, text=True,
                            capture_output=True, check=True)
    return result.stdout.strip().splitlines()


def build(path: Path) -> tuple[Path, Path, Path]:
    src = ROOT / "src"
    modules = ["kinds", "signal_dpss", "signal_windows_extra", "signal_spectral", "signal_fir_advanced", "signal_convolve", "signal_convolve2d", "signal_fftconvolve",
               "signal_detrend", "signal_peaks", "signal_savgol", "signal_lfilter",
               "signal_sos", "signal_filtfilt", "signal_filtfilt_nd", "signal_spectral_complex",
               "signal_analog_design", "signal_design", "signal_elliptic", "signal_iir_extra", "signal_analog_sos", "signal_fourier_resample", "signal_resample", "signal_analog", "signal_zpk_conversion", "signal_partial_fraction", "signal_sequences", "signal_hilbert", "signal_waveforms", "signal_periodic", "signal_irregular", "signal"]
    cmd = [FORTRAN, "-std=f2018", "-ffree-line-length-none", "-O2", "-c"]
    cmd.extend(str(src / f"{m}.f90") for m in modules)
    subprocess.run(cmd, cwd=path, check=True, capture_output=True)
    built = []
    for name, source in [("design", DESIGN_DRIVER), ("resample", RESAMPLE_DRIVER),
                         ("decimate", DECIMATE_DRIVER)]:
        driver = path / f"{name}.f90"
        driver.write_text(source)
        exe = path / name
        subprocess.run([FORTRAN, "-std=f2018", "-ffree-line-length-none", "-O2",
                        str(driver), *[str(path / f"{m}.o") for m in modules],
                        "-o", str(exe)], cwd=path, check=True, capture_output=True)
        built.append(exe)
    return tuple(built)


def main() -> None:
    rng = np.random.default_rng(20261015)
    checked = 0
    with tempfile.TemporaryDirectory(prefix="scifort_signal_ref_") as directory:
        design, fourier, decimator = build(Path(directory))
        modes = ["lowpass", "highpass", "bandpass", "bandstop"]
        for order in range(1, 9):
            for kind in range(1, 5):
                for ripple in [0.05, 0.5, 1.5]:
                    lines = run(design, f"{order} {kind} {ripple}\n")
                    assert int(lines[0]) == 0
                    b = np.fromstring(lines[2], sep=" ")
                    a = np.fromstring(lines[3], sep=" ")
                    cutoff = 0.32 if kind <= 2 else [0.21, 0.64]
                    ref_b, ref_a = signal.cheby1(order, ripple, cutoff, btype=modes[kind-1])
                    np.testing.assert_allclose(b, ref_b, atol=3e-8, rtol=3e-8)
                    np.testing.assert_allclose(a, ref_a, atol=3e-8, rtol=3e-8)
                    checked += 1
        for _ in range(650):
            n = int(rng.integers(1, 34))
            m = int(rng.integers(1, 40))
            kind = int(rng.integers(1, 3))
            win = int(rng.integers(0, 2))
            x = rng.normal(size=n)
            if kind == 2:
                x = x + 1j*rng.normal(size=n)
            request = f"{n} {m} {kind} {win}\n" + "".join(
                f"{z.real:.17g} {z.imag:.17g}\n" for z in x)
            lines = run(fourier, request)
            status, count = [int(v) for v in lines[0].split()]
            assert status == 0 and count == m
            y = np.array([complex(*[float(v) for v in line.split()])
                          for line in lines[1:]])
            reference = signal.resample(x, m, window="hann" if win else None)
            np.testing.assert_allclose(y, reference, atol=4e-12, rtol=4e-12)
            checked += 1
        for _ in range(180):
            q = int(rng.integers(2, 9))
            order = int(rng.choice([0, 1, 2, 3, 4, 5, 8, 10]))
            phase = int(rng.integers(0, 2))
            n = int(rng.integers(72, 142))
            x = rng.normal(size=n)
            request = f"{n} {q} {order} {phase}\n" + "\n".join(map(str, x)) + "\n"
            lines = run(decimator, request)
            status, count = [int(v) for v in lines[0].split()]
            assert status == 0
            y = np.fromstring(" ".join(lines[1:]), sep=" ")
            assert len(y) == count
            reference = signal.decimate(x, q, n=order or None, zero_phase=bool(phase))
            np.testing.assert_allclose(y, reference, atol=3e-10, rtol=3e-10)
            checked += 1
    print(f"SciFort Chebyshev/IIR/Fourier resampling reference checks passed: {checked}")


if __name__ == "__main__":
    main()
