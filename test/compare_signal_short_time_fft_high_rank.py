#!/usr/bin/env python3
"""Compare rank-6..14 Fortran STFTs with per-channel SciPy ShortTimeFFT.

Usage: python test/compare_signal_short_time_fft_high_rank.py EXECUTABLE
SciPy 1.17 cannot directly evaluate some padded rank>3 arrays, so this
reference explicitly applies the SciPy 1-D operation to each independent
channel and restores the frequency and time axes.
"""
from pathlib import Path
from subprocess import run
from sys import argv
from tempfile import TemporaryDirectory

import numpy as np
from scipy.signal import ShortTimeFFT, get_window

MODES = ('onesided', 'onesided2X', 'twosided', 'centered')
PADDINGS = ('zeros', 'edge', 'even', 'odd')
TRENDS = ('none', 'constant', 'linear')


def all_cases():
    cases = []
    for rank in range(6, 15):
        for axis in sorted({1, (rank + 1) // 2, rank}):
            for j in range(8):
                mode = MODES[(j + rank) % len(MODES)]
                kind = 'complex' if mode in ('twosided', 'centered') and j % 2 else 'real'
                pad = PADDINGS[(j + 2 * axis) % len(PADDINGS)]
                trend = TRENDS[(j + axis) % len(TRENDS)]
                scale = ('none', 'psd', 'magnitude')[j % 3]
                if mode == 'onesided2X' and scale == 'none':
                    scale = 'psd'
                n = 23 if j % 2 else 19
                nw = (7, 8, 9)[(j + axis) % 3]
                hop = 2 + (j % 2)
                nf = nw + 2 + ((j + rank) % 3)
                fs = 21. + j % 3
                phase = j % 2
                plan = ShortTimeFFT(get_window('hann', nw), hop, fs, mfft=nf,
                                    fft_mode=mode, phase_shift=phase,
                                    scale_to=None if scale == 'none' else scale)
                p0, p1 = plan.p_min + (j % 2), plan.p_max(n) - (j % 3)
                off = (j % 5) - 2
                assert p0 < p1
                cases.append((rank, axis, n, nw, hop, nf, fs, mode, scale,
                              phase, pad, kind, p0, p1, off, trend))
    return cases


def independent_scipy_stft(plan, x, axis, opts, trend):
    moved = np.moveaxis(x, axis - 1, -1)
    frames = []
    for idx in np.ndindex(moved.shape[:-1]):
        series = moved[idx]
        spectrum = (plan.stft(series, **opts) if trend == 'none'
                    else plan.stft_detrend(series, trend, **opts))
        frames.append(spectrum)
    spectra = np.stack(frames).reshape(moved.shape[:-1] + frames[0].shape)
    return np.moveaxis(spectra, -2, axis - 1)


def compare(executable):
    cases = all_cases()
    maxerr = 0.
    nvalues = 0
    with TemporaryDirectory() as tmp:
        work = Path(tmp)
        (work / 'stft_high_rank_cases.txt').write_text(
            str(len(cases)) + '\n' + '\n'.join(' '.join(map(str, c)) for c in cases) + '\n')
        p = run([str(Path(executable).resolve())],cwd=work,capture_output=True,text=True)
        if p.returncode:
            raise RuntimeError(f'SciFort high-rank reference failed:\n{p.stdout}\n{p.stderr}')
        with (work / 'stft_high_rank_results.txt').open() as fh:
            for ic, case in enumerate(cases, 1):
                (rank, axis, n, nw, hop, nf, fs, mode, scale, phase, pad,
                 kind, p0, p1, off, trend) = case
                header = list(map(int, fh.readline().split()))
                assert header[:2] == [ic, rank] and header[-1] == p0, (ic, header)
                outshape = tuple(header[2:-1])
                assert len(outshape) == rank + 1
                flat = np.array([complex(*map(float, fh.readline().split()))
                                 for _ in range(int(np.prod(outshape)))])
                got = flat.reshape(outshape, order='F')
                dims = [1] * rank
                dims[axis - 1] = n
                other_axes = [k for k in range(rank) if k != axis - 1]
                dims[other_axes[0]] = 2
                dims[other_axes[-1]] = 2
                i = np.arange(1, int(np.prod(dims)) + 1, dtype=np.float64)
                x = (np.sin(0.071 * i) + 0.2 * np.cos(0.13 * i)).reshape(dims, order='F')
                if kind == 'complex':
                    x = x + 0.15j * np.sin(0.11 * i).reshape(dims, order='F')
                plan = ShortTimeFFT(get_window('hann', nw), hop, fs, mfft=nf,
                                    fft_mode=mode, phase_shift=phase,
                                    scale_to=None if scale == 'none' else scale)
                ref = independent_scipy_stft(plan, x, axis,
                    dict(p0=p0, p1=p1, k_offset=off, padding=pad),trend)
                assert got.shape == ref.shape, (case, got.shape, ref.shape)
                err = float(np.max(np.abs(got-ref)))
                maxerr = max(maxerr, err)
                nvalues += got.size
                if err > 3e-11 * max(1., np.max(abs(ref))):
                    index = np.unravel_index(np.argmax(abs(got-ref)), got.shape)
                    raise AssertionError((case, err, index, got[index], ref[index]))
    print(f'PASS STFT rank6..14: {len(cases)} cases, {nvalues} complex values, '
          f'max absolute error {maxerr:.3e}')


if __name__ == '__main__':
    compare(argv[1])
