#!/usr/bin/env python3
"""Independent NumPy/SciPy cross-validation of leave-one-out multitaper shrinkage.

Checks 12 configurations (static and moving-window, real and complex,
all three weighting schemes), including every frequency/time bin and their
cross-spectral powers. Does not rely on SciFort's implementation for the
DPSS, FFT, weight updates, or the leave-one-out shrinkage objective.
"""
from __future__ import annotations
import sys
import subprocess

import numpy as np
from scipy.signal.windows import dpss


def samples() -> tuple[np.ndarray, np.ndarray]:
    t = np.arange(1, 97, dtype=float)[:, None]
    ch = np.arange(1, 5, dtype=float)[None, :]
    x = (np.sin(0.08 * t + 0.5 * ch)
         + 0.3 * np.cos(t * (0.11 + 0.015 * ch))
         + 0.18 * np.sin(t * (0.33 + 0.004 * ch * ch)))
    z = x + 1j * (0.26 * np.cos(t * (0.09 + 0.014 * ch))
                  + 0.17 * np.sin(t * 0.23 + ch))
    return x, z


def weights(p: np.ndarray, eig: np.ndarray, variance: float, mode: str) -> np.ndarray:
    """Return per-frequency weights, shape (taper, frequency)."""
    nt, nf = p.shape
    q = np.zeros_like(p)
    for f in range(nf):
        if mode == 'equal':
            q[:, f] = 1.0 / nt
            continue
        if mode == 'eigen' or nt < 3:
            w = np.maximum(eig, 0.0)
            q[:, f] = w / w.sum() if w.sum() > 0 else 1.0 / nt
            continue
        initial = np.maximum(eig[:2], 0.0)
        s = np.dot(initial, p[:2, f]) / initial.sum() if initial.sum() > 0 else 0.0
        cutoff = max(np.finfo(float).tiny * 1000, np.finfo(float).eps * max(s, variance))
        for _ in range(200):
            lambdas = np.clip(eig, 0.0, 1.0)
            denominator = lambdas * s + (1.0 - lambdas) * variance
            w = np.zeros(nt)
            mask = denominator > cutoff
            w[mask] = lambdas[mask] * (s / denominator[mask]) ** 2
            total = w.sum()
            if total <= np.finfo(float).tiny:
                w = np.maximum(eig, 0.0)
                q[:, f] = w / w.sum() if w.sum() > np.finfo(float).tiny else 1.0 / nt
                break
            q[:, f] = w / total
            snext = np.dot(q[:, f], p[:, f])
            delta = abs(snext - s) / max(snext, s, cutoff)
            s = snext
            if delta < 1e-10:
                break
    return q


def independent_reference(x: np.ndarray, nw: float, kmax: int, nfft: int,
                          mode: str) -> tuple[np.ndarray, np.ndarray]:
    """Compute cross spectra and LOO coefficients using a separate vectorized path."""
    n, nc = x.shape
    tapers, eig = dpss(n, nw, Kmax=kmax, sym=False, norm=2, return_ratios=True)
    detrended = x - np.mean(x, axis=0)
    variance = np.sum(abs(detrended) ** 2, axis=0) / n
    fft = np.fft.fft(detrended[:, None, :] * tapers.T[:, :, None], n=nfft, axis=0)
    fft /= np.sqrt(np.sum(tapers**2, axis=1))[None, :, None]
    if np.isrealobj(x):
        fft = fft[:nfft//2+1, :, :]
    nf = len(fft)
    periodograms = abs(fft) ** 2
    full_weights = np.stack([weights(periodograms[:, :, c].T, eig, variance[c], mode).T
                             for c in range(nc)], axis=-1)  # (frequency,taper,channel)
    csd = np.einsum('ftc,ftd,ftc,ftd->fcd', np.sqrt(full_weights),
                    np.sqrt(full_weights), np.conj(fft), fft, optimize=True)
    if np.isrealobj(x):
        if nfft % 2 == 0:
            csd[1:-1] *= 2
        else:
            csd[1:] *= 2
    numerator = np.zeros(nf)
    denominator = np.zeros(nf)
    for leave in range(kmax):
        keep = [k for k in range(kmax) if k != leave]
        loo = fft[:, keep, :]
        loo_weights = np.stack(
            [weights(periodograms[:, keep, c].T, eig[keep], variance[c], mode).T
             for c in range(nc)], axis=-1)
        training = np.einsum('ftc,ftd,ftc,ftd->fcd', np.sqrt(loo_weights),
                             np.sqrt(loo_weights), np.conj(loo), loo, optimize=True)
        heldout = np.conj(fft[:, leave, :, None]) * fft[:, leave, None, :]
        scale = np.max(abs(fft), axis=(1, 2))
        scale = np.maximum(scale, np.finfo(float).tiny)
        train = training / scale[:, None, None]**2
        test = heldout / scale[:, None, None]**2
        tri = np.triu_indices(nc, k=1)
        numerator += np.real(np.sum(np.conj(train[:, tri[0], tri[1]])
                                    * test[:, tri[0], tri[1]], axis=1))
        denominator += np.sum(abs(train[:, tri[0], tri[1]])**2, axis=1)
    alpha = np.ones(nf)
    positive = denominator > np.finfo(float).tiny
    alpha[positive] = np.clip(1-numerator[positive]/denominator[positive],0,1)
    return alpha, csd


def compare(exe: str) -> None:
    data = subprocess.check_output([exe], text=True).strip().splitlines()
    x, z = samples()
    current = 0
    worst_alpha = 0.0
    worst_power = 0.0
    ncases = 0
    nbins = 0
    while current < len(data):
        title, mode, kind, nfreq, nframe = data[current].split()
        current += 1
        nfreq, nframe = int(nfreq), int(nframe)
        raw = np.array([np.fromstring(row, sep=' ') for row in data[current:current+nfreq*nframe]])
        current += nfreq*nframe
        assert raw.shape == (nfreq*nframe, 5), (title, mode, kind, raw.shape)
        raw = raw.reshape(nframe,nfreq,5)
        obs = x if kind == 'real' else z
        expected_alpha = np.empty((nframe, nfreq))
        expected_power = np.empty((nframe, nfreq, 4))
        for frame in range(nframe):
            if title == 'STATIC':
                seg = obs
                nw, nfft = 3.5, (128 if kind == 'real' else 129)
            else:
                seg = obs[frame*16:frame*16+32]
                nw, nfft = 2.5, 48
            alpha, csd = independent_reference(seg, nw, 6, nfft, mode)
            expected_alpha[frame] = alpha
            expected_power[frame] = np.array(
                [np.real(csd[:,0,0]), np.real(csd[:,1,1]),
                 np.real(csd[:,0,1]), np.imag(csd[:,0,1])]).T
        alpha_delta = float(np.max(abs(raw[:,:,0]-expected_alpha)))
        power_delta = float(np.max(abs(raw[:,:,1:]-expected_power)))
        worst_alpha = max(worst_alpha, alpha_delta)
        worst_power = max(worst_power, power_delta)
        assert np.allclose(raw[:,:,0],expected_alpha,atol=1.5e-5,rtol=1.5e-5), (
            title,mode,kind,'alpha',alpha_delta)
        assert np.allclose(raw[:,:,1:],expected_power,atol=2.0e-5,rtol=2.0e-5), (
            title,mode,kind,'spectrum',power_delta)
        ncases += 1
        nbins += nframe*nfreq
    assert ncases==12, ncases
    print(f'PASS {ncases} SciPy/NumPy static/TF real/complex cases, {nbins} bins; '
          f'max alpha error {worst_alpha:.3g}, max CSD error {worst_power:.3g}',flush=True)


if __name__ == '__main__':
    compare(sys.argv[1])
