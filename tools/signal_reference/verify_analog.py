#!/usr/bin/env python3
"""Compare SciFort analog/digital frequency responses with public SciPy APIs.

Requires NumPy, SciPy, and gfortran; no runtime SciFort dependency is added.
"""
from __future__ import annotations

from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[2]
SEED = 20261016
COUNT = 1000


def make_complex_roots(rng, degree, real_coeffs):
    if real_coeffs:
        roots = []
        while len(roots) < degree:
            magnitude = 10.0 ** rng.uniform(-1.0, 1.0)
            if degree - len(roots) >= 2 and rng.random() < 0.6:
                root = complex(-magnitude, magnitude * rng.uniform(0.4, 1.6))
                roots.extend((root, root.conjugate()))
            else:
                roots.append(complex(-magnitude))
        return np.asarray(roots, dtype=complex)
    return rng.uniform(-5, 4, degree) + 1j * rng.uniform(-5, 5, degree)


def format_complex(values):
    return ' '.join(f'{t:.17e}' for z in values for t in (z.real, z.imag)) + '\n'


def main():
    rng = np.random.default_rng(SEED)
    chunks, references, metadata = [], [], []
    # 200 real TF, 200 complex TF, 200 analog ZPK, 200 digital ZPK,
    # and 200 automatic-grid cases (100 analog TF, 100 ZPK).
    for idx in range(COUNT):
        family = idx % 5
        autogrid = idx >= 800
        nw = int(rng.integers(1, 33))
        nz, npoles = (int(rng.integers(0, 5)), int(rng.integers(1, 7)))
        real_coeffs = family == 0 or (autogrid and family % 2 == 0)
        z = make_complex_roots(rng, nz, real_coeffs)
        p = make_complex_roots(rng, npoles, real_coeffs)
        if autogrid and not real_coeffs:
            p[0] = complex(p[0].real,abs(p[0].imag))
        # Avoid exact poles on sampled points; no singularity comparisons here.
        k = float(rng.uniform(0.25, 3.0))
        fs = float(rng.uniform(10.0, 1000.0))
        whole = bool(idx % 2)
        mode = int(autogrid)
        if family in (0, 1) or (autogrid and family in (2, 3)):
            if family in (0, 2):
                command = 'ATF'
                z = make_complex_roots(rng, nz, True)
                p = make_complex_roots(rng, npoles, True)
            else:
                command = 'CTF'
            numerator = np.atleast_1d(k * np.poly(z))
            denominator = np.atleast_1d(np.poly(p))
            # Complex values can be multiplied by complex constants for CTF.
            if command == 'CTF' and idx % 7 == 0:
                numerator = numerator * (1.0 + 0.1j)
            if command == 'ATF':
                numerator = np.asarray(numerator.real, dtype=float)
                denominator = np.asarray(denominator.real, dtype=float)
            if autogrid:
                expected = signal.freqs(numerator, denominator, worN=nw)
            else:
                grid = np.sort(rng.uniform(-12.0, 15.0, nw))
                expected = signal.freqs(numerator, denominator, worN=grid)
        elif family == 4 and not autogrid:
            command = 'DZPK'
            numerator, denominator = z, p
            grid = np.sort(rng.uniform(-fs/3, fs, nw))
            expected = signal.freqz_zpk(z,p,k,worN=grid,fs=fs)
        elif family == 4 and autogrid:
            command = 'AZPK'
            numerator, denominator = z, p
            expected = signal.freqs_zpk(z,p,k,worN=nw)
        elif family == 2 and not autogrid:
            command = 'AZPK'
            numerator, denominator = z, p
            grid = np.sort(rng.uniform(-10.0,15.0,nw))
            expected = signal.freqs_zpk(z,p,k,worN=grid)
        else:
            command = 'DZPK'
            numerator, denominator = z, p
            if autogrid:
                expected = signal.freqz_zpk(z,p,k,worN=nw,fs=fs,whole=whole)
            else:
                grid = np.sort(rng.uniform(-fs/3,fs,nw))
                expected = signal.freqz_zpk(z,p,k,worN=grid,fs=fs)
        # In digital cases, process a mix of auto-grid and caller-supplied grids.
        chunks.append(f'{command} {len(numerator)} {len(denominator)} {nw} {mode} '
                      f'{k:.17e} 0 {fs:.17e} {int(whole)}\n')
        if len(numerator):
            chunks.append(format_complex(np.asarray(numerator,dtype=complex)))
        if len(denominator):
            chunks.append(format_complex(np.asarray(denominator,dtype=complex)))
        if not autogrid:
            chunks.append(' '.join(f'{t:.17e}' for t in grid) + '\n')
        references.append(expected)
        metadata.append(command)

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        exe = tmp / 'analog_reference'
        cmd = ['gfortran', '-std=f2018', '-O2', '-fcheck=all',
               '-Wimplicit-interface', '-Werror=implicit-interface',
               '-J', str(tmp), '-I', str(tmp), '-o', str(exe),
               *(str(ROOT / 'src' / f) for f in ['kinds.f90',
                  'signal_dpss.f90', 'signal_windows_extra.f90', 'signal_spectral.f90', 'signal_analog_design.f90', 'signal_design.f90', 'signal_analog.f90']),
               str(ROOT / 'tools/signal_reference/analog_driver.f90')]
        subprocess.run(cmd,check=True,cwd=tmp)
        run = subprocess.run([str(exe)],input=''.join(chunks),capture_output=True,
                             text=True,cwd=tmp,check=True)

    records = iter(run.stdout.splitlines())
    max_abs, total_bins = 0.0, 0
    for idx, ((want_w,want_h),command) in enumerate(zip(references,metadata)):
        status, n = map(int,next(records).split())
        assert status == 0, (idx,command,'status',status)
        assert n == len(want_w), (idx,command,'size',n,len(want_w))
        actual = np.array([[float(v) for v in next(records).split()]
                           for _ in range(n)]).reshape(n,3)
        w_err = np.abs(actual[:,0] - want_w)
        h_err = np.abs((actual[:,1]+1j*actual[:,2]) - want_h)
        max_abs = max(max_abs,float(np.max(h_err)))
        total_bins += n
        try:
            np.testing.assert_allclose(actual[:,0],want_w,rtol=2e-9,atol=1e-10)
            np.testing.assert_allclose(actual[:,1]+1j*actual[:,2],want_h,
                                       rtol=3e-9,atol=5e-11)
        except AssertionError as exc:
            raise AssertionError(f'case {idx} ({command}): '\
                f'w_err={w_err.max():.3e} h_err={h_err.max():.3e}\n{exc}') from exc
    print(f'PASS: {COUNT} randomized cases, {total_bins} frequency bins, '
          f'max absolute response error {max_abs:.3e}')


if __name__ == '__main__':
    main()
