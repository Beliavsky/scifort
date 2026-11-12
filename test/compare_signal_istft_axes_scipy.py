#!/usr/bin/env python3
"""Compare SciFort inverse-STFT axis permutations with SciPy 1.17.

Usage: python test/compare_signal_istft_axes_scipy.py FORTRAN_DRIVER
"""
import itertools
import subprocess
import sys
import warnings

import numpy as np
from scipy.signal import istft


def main(executable):
    rng = np.random.default_rng(20261026)
    ncases = 0
    nsamples = 0
    max_error = 0.0
    for permutation in itertools.permutations((0, 1, 2)):
        fa = permutation.index(0)
        ta = permutation.index(2)
        output_axis = fa - (ta < fa)
        scipy_output_axis = ta - (fa < ta)
        for onesided,real_output in ((True, True),(True,False),(False,False)):
            for nw,no,nfft in ((6,3,6),(7,4,7),(8,6,12)):
                nf = nfft//2+1 if onesided else nfft
                nframes = 5
                nchannels = 2
                base = rng.standard_normal((nf,nchannels,nframes)) + 1j*rng.standard_normal((nf,nchannels,nframes))
                z = np.transpose(base, permutation)
                for window in ('hann','hamming','boxcar'):
                    for scaling in ('spectrum','psd'):
                        for boundary in (True,False):
                            if not boundary and window != 'boxcar':
                                continue
                            with warnings.catch_warnings():
                                warnings.simplefilter('ignore', UserWarning)
                                tref, expect = istft(
                                    z, fs=11.25, window=window, nperseg=nw,
                                    noverlap=no, nfft=nfft, boundary=boundary,
                                    input_onesided=onesided, scaling=scaling,
                                    freq_axis=fa,time_axis=ta)
                            expect = np.moveaxis(expect, scipy_output_axis, output_axis)
                            if real_output:
                                expect = expect.real
                            # SciPy 1.17 incorrectly sizes its reported time
                            # vector when freq_axis/time_axis are reordered;
                            # the reconstructed samples themselves are correct.
                            nout = expect.shape[output_axis]
                            reference_times = np.arange(nout, dtype=float)/11.25
                            if len(tref) == nout and not np.allclose(tref, reference_times):
                                raise AssertionError(('SciPy times',permutation))
                            dims = z.shape
                            for mode in range(3):
                                header = (' '.join(str(v) for v in (
                                    *dims,nout,fa+1,ta+1,nw,no,nfft,
                                    int(real_output),int(boundary),int(onesided),
                                    11.25,mode)) + f'\n{window} {scaling}\n')
                                data = ''.join(f'{v.real:.17g} {v.imag:.17g}\n' for v in z.ravel(order='F'))
                                p = subprocess.run([executable],input=header+data,
                                                   text=True,capture_output=True,check=True)
                                lines = p.stdout.splitlines()
                                if not lines or int(lines[0]) != 0:
                                    raise AssertionError((permutation, mode,onesided,real_output,nw,no,nfft,window,scaling,boundary,dims,nout,lines,p.stderr))
                                if len(lines) != 1+nout+expect.size:
                                    raise AssertionError((permutation,mode,len(lines),expect.size,nout))
                                times = np.asarray([float(v) for v in lines[1:1+nout]])
                                found = np.asarray([[float(v) for v in line.split()]
                                                    for line in lines[1+nout:]])
                                found = (found[:,0]+1j*found[:,1]).reshape(expect.shape,order='F')
                                if not np.allclose(times,reference_times,atol=2e-13,rtol=2e-13):
                                    raise AssertionError(('times',permutation,mode))
                                error = float(np.max(np.abs(found-expect)))
                                max_error=max(max_error,error)
                                if not np.allclose(found,expect,atol=3e-12,rtol=2e-11):
                                    raise AssertionError((permutation,mode,nw,window,
                                                          scaling,boundary,error))
                                ncases += 1
                                nsamples += expect.size
    print(f'PASS: {ncases} independent SciPy inverse-STFT axis cases, '
          f'{nsamples} samples, max absolute error {max_error:.3g}')


if __name__ == '__main__':
    if len(sys.argv)!=2:
        raise SystemExit('usage: compare_signal_istft_axes_scipy.py DRIVER')
    main(sys.argv[1])
