#!/usr/bin/env python3
"""Compare incremental Fortran inverse STFT against independent scipy.signal.istft.

Usage: python test/compare_signal_istft_stream_scipy.py FORTRAN_DRIVER
"""
import itertools
import subprocess
import sys
import warnings

import numpy as np
from scipy.signal import istft


def main(executable):
    rng = np.random.default_rng(20261025)
    cases = 0
    samples = 0
    worst = 0.0
    for nw, fft_extra, frames, overlap in itertools.product(
        (1, 5, 6, 7, 8, 13), (0, 1, 4), (1, 2, 7), (0, 1, 2)
    ):
        if overlap == 1:
            noverlap = nw // 2
        elif overlap == 2:
            noverlap = nw - 1
        else:
            noverlap = 0
        nfft = nw + fft_extra
        for onesided, real_output in ((True, True), (True, False), (False, False)):
            # Multiple reproducible subsets avoid running thousands of subprocesses.
            selector = nw * 11 + fft_extra * 3 + frames + overlap * 7 + int(onesided) + int(real_output)
            if selector % 4:
                continue
            nfreq = nfft // 2 + 1 if onesided else nfft
            zxx = rng.normal(size=(nfreq, frames)) + 1j * rng.normal(size=(nfreq, frames))
            for window, scaling, boundary in itertools.product(
                ('hann', 'hamming', 'boxcar'), ('spectrum', 'psd'), (True, False)
            ):
                # SciPy 1.17 returns an empty array when nperseg=1 and boundary=True
                # (slice[0:-0]); SciFort correctly retains the single sample.
                if nw == 1 and boundary:
                    continue
                fs = 11.25
                with warnings.catch_warnings():
                    warnings.filterwarnings('ignore', category=UserWarning)
                    reference_times, expected = istft(
                        zxx, fs=fs, window=window, nperseg=nw, noverlap=noverlap,
                        nfft=nfft, boundary=boundary, input_onesided=onesided,
                        scaling=scaling,
                    )
                if real_output:
                    expected = np.real(expected)
                source = f'{nw} {nfft} {noverlap} {frames} {int(onesided)} {int(boundary)} {int(real_output)} {fs}\n'
                source += f'{window} {scaling}\n'
                source += ''.join(f'{value.real:.17g} {value.imag:.17g}\n'
                                  for value in zxx.flatten(order='F'))
                completed = subprocess.run([executable], input=source, capture_output=True,
                                           text=True, check=True)
                lines = completed.stdout.splitlines()
                status, count = (int(x) for x in lines[0].split())
                if status not in (0, 2) or count != expected.size or len(lines) != count+1:
                    raise AssertionError((nw,nfft,noverlap,frames,onesided,real_output,
                                          window,scaling,boundary,status,count,len(lines),completed.stderr))
                actual = np.array([[float(x) for x in line.split()] for line in lines[1:]])
                actual = actual.reshape((-1, 3))
                values = actual[:,1] + 1j * actual[:,2]
                if real_output:
                    values = values.real
                difference = np.max(np.abs(values - expected), initial=0.0)
                worst = max(worst, float(difference))
                if not np.allclose(values,expected,rtol=5e-12,atol=2e-12):
                    raise AssertionError((nw,nfft,noverlap,frames,onesided,real_output,
                                          window,scaling,boundary,difference))
                if not np.allclose(actual[:,0],reference_times,rtol=5e-13,atol=5e-13):
                    raise AssertionError(('time',nw,nfft,noverlap,frames,boundary))
                cases += 1
                samples += expected.size
    print(f'PASS: {cases} independent SciPy streaming-iSTFT cases, '
          f'{samples} samples, max error {worst:.3g}')


if __name__ == '__main__':
    if len(sys.argv)!=2:
        raise SystemExit('usage: compare_signal_istft_stream_scipy.py DRIVER')
    main(sys.argv[1])
