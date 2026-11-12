#!/usr/bin/env python3
"""Compare ShortTimeFFT coordinates/bounds from Fortran with SciPy.

Usage: python test/compare_signal_short_time_fft_helpers.py EXE
"""
import pathlib
import subprocess
import sys
import tempfile

import numpy as np
from scipy.signal import ShortTimeFFT, get_window


def cases():
    data = []
    for iw, name in enumerate(('hann', 'hamming', 'boxcar')):
        for nw in (7, 8, 11):
            for hop in (1, 2, 3):
                for mode in ('onesided', 'onesided2X', 'twosided', 'centered'):
                    for n in (nw - nw // 2, nw + 13):
                        nf = nw + (iw + hop + n) % 4
                        fs = float(23 + iw)
                        off = 5 - iw - hop
                        data.append((name, nw, hop, nf, n, fs, mode, off))
    return data


def compare(exe):
    data = cases()
    maxerr = 0.
    nvalues = 0
    with tempfile.TemporaryDirectory() as work:
        work = pathlib.Path(work)
        (work / 'stft_helper_cases.txt').write_text(
            str(len(data)) + '\n' + '\n'.join(' '.join(map(str, c)) for c in data) + '\n'
        )
        p = subprocess.run([str(pathlib.Path(exe).resolve())], cwd=work,
                           text=True, capture_output=True)
        if p.returncode:
            raise RuntimeError(f'Fortran helper driver failed:\n{p.stdout}\n{p.stderr}')
        with (work / 'stft_helper_results.txt').open() as fh:
            for ic, case in enumerate(data, 1):
                name, nw, hop, nf, n, fs, mode, off = case
                sft = ShortTimeFFT(get_window(name, nw), hop, fs, mfft=nf,
                                   fft_mode=mode, scale_to='psd')
                values = list(map(int, fh.readline().split()))
                assert values[0] == ic, (case, values)
                expected = [ic, sft.p_min, sft.p_max(n), sft.k_min, sft.k_max(n),
                            *sft.lower_border_end, *sft.upper_border_begin(n),
                            sft.nearest_k_p(-7), sft.nearest_k_p(-7, left=False),
                            sft.f_pts, *sft.p_range(n),
                            1 if mode == 'twosided' else 0]
                if values != expected:
                    raise AssertionError((case, 'indices', values, expected))
                nums = np.array(list(map(float, fh.readline().split())))
                ext = np.zeros(4) if mode == 'twosided' else np.asarray(sft.extent(n, center_bins=True))
                expected_num = np.r_[sft.delta_t, sft.delta_f, ext]
                err = float(np.max(abs(nums-expected_num)))
                if err > 1.e-11:
                    raise AssertionError((case, 'extent and steps', nums, expected_num))
                maxerr = max(maxerr, err)
                nfreq, ntimes = map(int, fh.readline().split())
                freq = np.fromstring(fh.readline(), sep=' ')
                times = np.fromstring(fh.readline(), sep=' ')
                assert len(freq) == nfreq and len(times) == ntimes
                if nfreq != len(sft.f) or ntimes != len(sft.t(n, k_offset=off)):
                    raise AssertionError((case, 'array sizes', nfreq, ntimes))
                for got, ref in ((freq, sft.f), (times, sft.t(n, k_offset=off))):
                    err = float(np.max(abs(got-ref)))
                    if err > 1.e-11:
                        raise AssertionError((case, 'array values', got, ref))
                    maxerr = max(maxerr, err)
                    nvalues += len(got)
    print(f'PASS ShortTimeFFT helpers: {len(data)} cases, {nvalues} coordinates, '
          f'max absolute error {maxerr:.3e}')


if __name__ == '__main__':
    compare(sys.argv[1])
