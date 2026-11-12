#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent kernel, stationarity and serial-protocol checks for the study."""
from __future__ import annotations

import numpy as np
from scipy.signal.windows import hann

from study_canonical_varp_conditional import (
    BANDS, FREQUENCIES, companion_and_covariance, encode_dataset, models,
    parse_output, population_targets, simulate_windows, window_spectra,
)


def run() -> None:
    rng = np.random.default_rng(884211)
    length, nw = 17, 5
    taper = hann(length, sym=False)
    for model in models():
        f, covariance = companion_and_covariance(model)
        if np.min(np.linalg.eigvalsh(covariance)) <= 0:
            raise AssertionError('non-positive-definite stationary companion covariance')
        windows = simulate_windows(rng, model, nw, length)
        observed = window_spectra(windows, taper)
        assert observed.shape == (3, nw, 5, 5)
        assert np.allclose(observed, observed.swapaxes(-2,-1).conj(), atol=3e-13)
        assert np.min(np.linalg.eigvalsh(observed[:, :, :, :]).real) > -1e-11
        # Check exact analytic band kernels against a completely different
        # 128-node Gauss-Legendre integration of complex tapered Fourier sums.
        nodes, weights = np.polynomial.legendre.leggauss(128)
        time = np.arange(length)
        for j,(left,right) in enumerate(BANDS):
            freqs = (right-left)*nodes/2+(right+left)/2
            phase = np.exp(-2j*np.pi*freqs[:,None]*time[None,:])
            x = windows * taper[None,:,None]
            z = np.einsum('wtc,ft->wfc', x, phase)
            matrix = np.einsum('f,wfi,wfj->wij', weights/2,
                               z, z.conj())/np.dot(taper,taper)
            if not np.allclose(observed[j+1], matrix, rtol=2e-11, atol=2e-11):
                raise AssertionError('finite-window band quadrature disagreement')
        target = population_targets(model, taper)
        if target.shape != (3,3):
            raise AssertionError('unexpected target shape')
        if not np.all((target >= -1e-12) & (target <= 1+1e-12)):
            raise AssertionError('population squared coherence outside [0,1]')
        if not np.all(target[:,0] >= target[:,1]-1e-12):
            raise AssertionError('canonical modes are not descending')
        aggregate = 1-(1-target[:,0])*(1-target[:,1])
        if not np.allclose(aggregate,target[:,2],rtol=2e-12,atol=2e-12):
            raise AssertionError('determinant coherence and canonical roots disagree')
        protocol = encode_dataset(observed,seed=999,block_length=3,resamples=29)
        assert protocol[0] == '5 5 3 3 29 999'
        assert len(protocol) == 17
        assert len(protocol[1].split()) == 3
        assert len(protocol[2].split()) == 50
    # Unusable data must fail parsing, not silently count as invalid coverage.
    try:
        parse_output('1 2 3\n',1)
    except AssertionError:
        pass
    else:
        raise AssertionError('truncated results were accepted')
    print('PASS 4 VAR(p) scenario definitions, stationary initialization, '
          'exact band quadrature, coherent canonical targets and protocol checks')


if __name__ == '__main__':
    run()
