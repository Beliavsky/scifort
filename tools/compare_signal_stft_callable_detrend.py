#!/usr/bin/env python3
"""Compare per-slice ShortTimeFFT callable detrending against SciPy 1.17."""
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
from scipy import signal


def real_detrend(samples: np.ndarray) -> np.ndarray:
    """Apply an independent nonlinear, non-affine real detrend."""
    return samples - samples.mean() + 0.07 * samples**2


def complex_detrend(samples: np.ndarray) -> np.ndarray:
    """Apply an independent nonlinear complex slice transformation."""
    return samples - samples.mean() + 0.07 * samples**2 + 0.04j * samples


def compare(path: Path) -> None:
    """Read independent Fortran output and verify FFT coefficients and grids."""
    modes = {1: 'onesided', 2: 'onesided2X', 3: 'twosided', 4: 'centered'}
    pads = {1: 'zeros', 2: 'edge', 3: 'even', 4: 'odd'}
    scales = {1: None, 2: 'magnitude', 3: 'psd'}
    i = np.arange(1, 49, dtype=np.float64)
    x = np.sin(.15*i) + .2*np.cos(.4*i) + .003*i*i
    z = x + 1j*.1*np.cos(.3*i)
    ncases, ncoeff, max_error = 0, 0, 0.0
    with path.open(encoding='utf-8') as stream:
        while line := stream.readline():
            parts = line.split()
            assert parts[0] == 'S', parts
            (typ, wid, nw, hop, nf, imode, ipad, iscale, phase,
             nfreq, ntime) = map(int, parts[1:])
            p0 = int(stream.readline())
            freqs = np.array([float(stream.readline()) for _ in range(nfreq)])
            times = np.array([float(stream.readline()) for _ in range(ntime)])
            a = np.array([list(map(float, stream.readline().split()))
                          for _ in range(nfreq*ntime)], dtype=np.float64)
            output = (a[:, 0] + 1j*a[:, 1]).reshape(ntime, nfreq).T
            window = signal.get_window('hann' if wid == 1 else 'hamming', nw)
            plan = signal.ShortTimeFFT(window, hop, fs=20., mfft=nf,
                                       fft_mode=modes[imode], scale_to=scales[iscale],
                                       phase_shift=phase)
            data = x if typ == 1 else z
            detrend = real_detrend if typ == 1 else complex_detrend
            reference = plan.stft_detrend(data, detrend, padding=pads[ipad])
            assert reference.shape == output.shape, (reference.shape, output.shape)
            assert p0 == plan.p_min, (p0, plan.p_min)
            np.testing.assert_allclose(freqs, plan.f, atol=3e-14, rtol=3e-14)
            np.testing.assert_allclose(times, plan.t(len(data)), atol=3e-14, rtol=3e-14)
            err = float(np.max(np.abs(reference - output)))
            max_error = max(max_error, err)
            np.testing.assert_allclose(output, reference, atol=3e-12, rtol=3e-12,
                                       err_msg=f'{typ=}, {wid=}, {imode=}, {ipad=}, '
                                               f'{iscale=}, {phase=}')
            ncases += 1
            ncoeff += output.size
    assert ncases == 272, ncases
    print(f'PASS: {ncases} SciPy ShortTimeFFT callable detrending cases, '
          f'{ncoeff} complex coefficients, max abs error {max_error:.3e}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', type=Path)
    compare(parser.parse_args().results)
