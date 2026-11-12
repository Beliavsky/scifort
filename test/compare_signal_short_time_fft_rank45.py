#!/usr/bin/env python3
"""Compare rank-4/5 Fortran STFTs with SciPy, including complex signals.

Usage: python test/compare_signal_short_time_fft_rank45.py PATH_TO_FORTRAN_EXE
"""
import pathlib
import subprocess
import sys
import tempfile

import numpy as np
from scipy.signal import ShortTimeFFT, get_window

W1 = np.array([.13, .07, .11, .05, .03])
W2 = np.array([.17, -.1, .03, .08, .07])
W3 = np.array([.08, .12, -.09, .13, .05])
MODES = ('onesided', 'onesided2X', 'twosided', 'centered')
PADDINGS = ('zeros', 'edge', 'even', 'odd')
TRENDS = ('none', 'constant', 'linear')
SCALES = ('none', 'psd', 'magnitude')


def build_cases():
    cases = []
    for rank in (4, 5):
        for axis in range(1, rank + 1):
            for j in range(12):
                mode = MODES[(j + axis) % 4]
                trend = TRENDS[(j // 2 + axis) % 3]
                pad = PADDINGS[(j // 3 + 2 * axis) % 4]
                scale = SCALES[(j // 4 + axis) % 3]
                if mode == 'onesided2X' and scale == 'none':
                    scale = 'psd'
                kind = ('complex' if mode in ('twosided', 'centered') and j % 2 else 'real')
                n = 23 if j % 2 else 35
                nw = (7, 8, 9)[(j + axis) % 3]
                hop = (2, 3)[j % 2]
                nfft = (13, 16)[(j + axis) % 2]
                fs, phase = 20., (0, 2)[j % 2]
                sft = ShortTimeFFT(get_window('hann', nw), hop, fs, mfft=nfft,
                                   fft_mode=mode,
                                   scale_to=None if scale == 'none' else scale,
                                   phase_shift=phase)
                p0 = sft.p_min + (j % 2)
                p1 = sft.p_max(n) - (j % 3)
                offset = j % 5 - 2
                cases.append((rank, axis, n, nw, hop, nfft, fs, mode, scale,
                              phase, pad, kind, p0, p1, offset, trend))
    return cases


def run(executable):
    cases = build_cases()
    maxerr = 0.
    nvalues = 0
    ndirect = 0
    nfallback = 0
    with tempfile.TemporaryDirectory() as tmp:
        work = pathlib.Path(tmp)
        (work / 'short_fft_rank45_cases.txt').write_text(
            str(len(cases)) + '\n' + '\n'.join(' '.join(map(str, c)) for c in cases) + '\n')
        proc = subprocess.run([str(pathlib.Path(executable).resolve())], cwd=work,
                              capture_output=True, text=True)
        if proc.returncode:
            raise RuntimeError(f'Fortran reference failed: {proc.stderr}\n{proc.stdout}')
        with (work / 'short_fft_rank45_results.txt').open() as fh:
            for ci, case in enumerate(cases, 1):
                rank, axis, n, nw, hop, nfft, fs, mode, scale, phase, pad, kind, p0, p1, offset, trend = case
                header = list(map(int, fh.readline().split()))
                assert header[:2] == [ci, rank] and header[-1] == p0, (ci, header)
                outshape = tuple(header[2:-1])
                assert len(outshape) == rank + 1, (case, header)
                values = np.array([complex(*map(float, fh.readline().split()))
                                   for _ in range(int(np.prod(outshape)))])
                got = values.reshape(outshape, order='F')
                dims = [2] * rank
                dims[axis - 1] = n
                grid = (np.indices(dims, dtype=float) + 1).reshape(rank, -1)
                base = np.sin(W1[:rank] @ grid) + .2 * np.cos(W2[:rank] @ grid)
                imaginary = .17 * np.sin(W3[:rank] @ grid)
                x = (base + 1j * imaginary if kind == 'complex' else base).reshape(dims)
                plan = ShortTimeFFT(get_window('hann', nw), hop, fs, mfft=nfft,
                                    fft_mode=mode, phase_shift=phase,
                                    scale_to=None if scale == 'none' else scale)
                kwargs = dict(p0=p0, p1=p1, k_offset=offset, padding=pad,
                              axis=axis - 1)
                try:
                    ref = (plan.stft(x, **kwargs) if trend == 'none' else
                           plan.stft_detrend(x, trend, **kwargs))
                    ndirect += 1
                except ValueError:
                    # SciPy 1.17 may mishandle high-rank padding. Per-channel
                    # one-dimensional SciPy remains the numerical reference.
                    nfallback += 1
                    moved = np.moveaxis(x, axis - 1, -1)
                    opts = {k: v for k, v in kwargs.items() if k != 'axis'}
                    frames = []
                    for idx in np.ndindex(moved.shape[:-1]):
                        series = moved[idx]
                        frames.append(plan.stft(series, **opts) if trend == 'none'
                                      else plan.stft_detrend(series, trend, **opts))
                    ref = np.moveaxis(
                        np.stack(frames).reshape(moved.shape[:-1] + frames[0].shape),
                        -2, axis - 1)
                assert ref.shape == got.shape, (case, ref.shape, got.shape)
                err = float(np.max(abs(ref - got)))
                maxerr = max(maxerr, err)
                nvalues += got.size
                if err > 3.e-11 * max(1.0, np.max(abs(ref))):
                    pos = np.unravel_index(np.argmax(abs(ref - got)), ref.shape)
                    raise AssertionError(f'case {ci} {case}: error {err} at {pos}, '
                                         f'Fortran={got[pos]} SciPy={ref[pos]}')
            assert not fh.readline(), 'Extra output from Fortran'
    print(f'Rank-4/5 ShortTimeFFT SciPy comparisons: {len(cases)} cases, '
          f'{nvalues} coefficients, maximum error {maxerr:.3e}, '
          f'direct SciPy {ndirect}, per-channel fallback {nfallback}')


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('Usage: compare_signal_short_time_fft_rank45.py EXECUTABLE')
    run(sys.argv[1])
