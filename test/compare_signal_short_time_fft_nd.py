#!/usr/bin/env python3
"""Compare 2D/3D Fortran ShortTimeFFT outputs with scipy.signal.ShortTimeFFT.

Usage: python test/compare_signal_short_time_fft_nd.py /path/to/reference_executable
"""
import itertools
import pathlib
import subprocess
import sys
import tempfile

import numpy as np
from scipy.signal import ShortTimeFFT, get_window


def run(executable):
    cases = []
    for ix, (rank, ax, mode, trend, padding, scale, nfft) in enumerate(
        itertools.product(
            (2, 3),
            (1, 2, 3),
            ('onesided', 'onesided2X', 'twosided', 'centered'),
            ('none', 'constant', 'linear'),
            ('zeros', 'edge', 'even', 'odd'),
            ('none', 'psd', 'magnitude'),
            (13, 16),
        )
    ):
        if ax > rank or (mode == 'onesided2X' and scale == 'none'):
            continue
        if (ix * 31 + rank * 17 + ax * 7) % 29 > 1:
            continue
        n = (21, 35)[ix % 2]
        nw = (7, 8, 9)[ix % 3]
        hop = (2, 3)[ix % 2]
        is_complex = mode in ('centered', 'twosided') and ix % 3 != 0
        phase = (0, 2)[ix % 2]
        fs = 20.0
        sft = ShortTimeFFT(get_window('hann', nw), hop, fs, mfft=nfft,
                           fft_mode=mode, scale_to=scale if scale != 'none' else None,
                           phase_shift=phase)
        start = sft.p_min + (ix % 2)
        end = sft.p_max(n) - (ix % 3)
        if start >= end:
            continue
        offset = ix % 5 - 2
        cases.append((rank,ax,n,nw,hop,nfft,fs,mode,scale,phase,padding,
                      'complex' if is_complex else 'real',start,end,offset,trend))
    if len(cases) < 50:
        raise RuntimeError(f'Insufficient cases: {len(cases)}')
    with tempfile.TemporaryDirectory() as tmp:
        work = pathlib.Path(tmp)
        (work / 'short_fft_nd_cases.txt').write_text(
            f'{len(cases)}\n' + '\n'.join(' '.join(map(str,case)) for case in cases) + '\n'
        )
        p = subprocess.run([str(pathlib.Path(executable).resolve())],
                           cwd=work, capture_output=True, text=True)
        if p.returncode:
            raise RuntimeError(f'Fortran case generator failed: {p.stderr}\n{p.stdout}')
        values_count = 0
        maxerr = 0.0
        scipy_direct = 0
        scipy_workaround = 0
        with (work / 'short_fft_nd_results.txt').open() as fh:
            for ci,case in enumerate(cases,1):
                (rank,ax,n,nw,hop,nfft,fs,mode,scale,phase,padding,kind,p0,p1,
                 offset,trend) = case
                header = list(map(int, fh.readline().split()))
                assert header[0] == ci and header[1] == rank, (ci,header)
                gotshape = header[2:5] if rank == 2 else header[2:6]
                assert header[-1] == p0
                count = int(np.prod(gotshape))
                got = np.array([complex(*map(float, fh.readline().split()))
                                for _ in range(count)]).reshape(gotshape, order='F')
                if rank == 2:
                    dims = (n,3) if ax == 1 else (3,n)
                else:
                    dims = [(n,2,3),(2,n,3),(2,3,n)][ax-1]
                grid = np.indices(dims, dtype=np.float64) + 1
                i,j = grid[:2]
                base = np.sin(.13*i+.07*j)+.2*np.cos(.19*i-.1*j)
                imag = .17*np.sin(.11*i+.03*j)
                if rank == 3:
                    k = grid[2]
                    base = np.sin(.13*i+.07*j+.09*k)+.2*np.cos(.19*i-.1*j+.05*k)
                    imag = .17*np.sin(.11*i+.03*j-.02*k)
                x = base + 1j*imag if kind == 'complex' else base
                sft = ShortTimeFFT(get_window('hann',nw),hop,fs,mfft=nfft,
                                   fft_mode=mode,phase_shift=phase,
                                   scale_to=None if scale=='none' else scale)
                kwargs = dict(p0=p0,p1=p1,k_offset=offset,padding=padding,axis=ax-1)
                # SciPy 1.17's multidimensional padding currently calls np.hstack
                # internally, which fails for several axis/rank combinations when
                # mfft > window length. Per-channel 1-D SciPy remains the oracle.
                try:
                    ref = (sft.stft(x,**kwargs) if trend=='none' else
                           sft.stft_detrend(x,trend,**kwargs))
                    scipy_direct += 1
                except ValueError:
                    scipy_workaround += 1
                    moved = np.moveaxis(x,ax-1,-1)
                    pieces = []
                    opts = {k:v for k,v in kwargs.items() if k != 'axis'}
                    for ind in np.ndindex(moved.shape[:-1]):
                        row = moved[ind]
                        spec = (sft.stft(row,**opts) if trend=='none' else
                                sft.stft_detrend(row,trend,**opts))
                        pieces.append(spec)
                    ref = np.moveaxis(np.array(pieces).reshape(moved.shape[:-1]+spec.shape),-2,ax-1)
                if tuple(ref.shape) != tuple(gotshape):
                    raise AssertionError(f'case {ci}: {ref.shape} != {gotshape}: {case}')
                err = np.max(abs(ref-got))
                maxerr = max(err,maxerr)
                values_count += count
                if err > 3.e-11 * max(1.,np.max(abs(ref))):
                    idx=np.unravel_index(np.argmax(abs(ref-got)),ref.shape)
                    raise AssertionError(f'case {ci}: max err {err}, at {idx}, {got[idx]} != {ref[idx]}: {case}')
            if fh.readline():
                raise AssertionError('Trailing unexpected output')
    print(f'Multidimensional ShortTimeFFT against SciPy: {len(cases)} cases, '
          f'{values_count} complex coefficients, maximum absolute error {maxerr:.3e}; '
          f'{scipy_direct} direct multidimensional SciPy, {scipy_workaround} per-channel SciPy references')


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('Usage: compare_signal_short_time_fft_nd.py EXECUTABLE')
    run(sys.argv[1])
