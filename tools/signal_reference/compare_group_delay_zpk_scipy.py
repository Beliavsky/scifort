#!/usr/bin/env python3
"""Cross-check the Fortran ZPK group delay against independent SciPy algorithms."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import warnings

import numpy as np
from scipy import signal


def root_array(rng: np.random.Generator, n: int, real: bool, pole: bool) -> np.ndarray:
    radii = rng.uniform(0.04, 0.97 if pole else 1.45, n)
    if real:
        return rng.choice((-1.0, 1.0), size=n) * radii
    angles = rng.uniform(-np.pi, np.pi, n)
    return radii * np.exp(1j*angles)


def build_cases() -> list[tuple[np.ndarray, np.ndarray, complex, np.ndarray, float]]:
    rng = np.random.default_rng(20261028)
    cases = []
    sizes = (0, 1, 2, 3, 4, 8, 16, 32, 64, 128, 256)
    for n in sizes:
        for repetition in range(18):
            nz = n if repetition % 4 != 1 else max(0, n-2)
            npoles = n if repetition % 4 != 2 else max(0, n-1)
            if repetition == 3:
                npoles = n + 1
            real_roots = repetition % 3 == 0
            zeros = root_array(rng, nz, real_roots, pole=False)
            poles = root_array(rng, npoles, real_roots, pole=True)
            fs = (2*np.pi, 2.0, 48000.0)[repetition % 3]
            freqs = rng.uniform(0.03*fs, 0.96*fs, 31)
            gain = complex(rng.uniform(0.2, 2),rng.uniform(-2, 2))
            if real_roots:
                gain = complex(gain.real,0)
            cases.append((zeros,poles,gain,freqs,fs))
    # Butterworth cascades of increasing order, but with zeros at Nyquist.
    # Avoid sampling exactly at Nyquist, where group delay is undefined.
    for order in (4, 8, 16, 32, 64, 128):
        for kind in ('lowpass','highpass','bandpass'):
            cutoff = [0.21, 0.57] if kind == 'bandpass' else 0.36
            zeros,poles,gain = signal.butter(order, cutoff, btype=kind, output='zpk')
            fs = 2*np.pi
            freqs = np.linspace(0.017, 0.983, 27)*np.pi
            cases.append((zeros,poles,complex(gain),freqs,fs))
    # All-pass roots very close to the unit circle; neither root is exactly on it.
    for n in (2, 8, 32, 64, 128):
        angles = rng.uniform(0.3,2.8,n)
        poles = 0.999 * np.exp(1j*angles)
        zeros = 1/np.conj(poles)
        fs = 2*np.pi
        freqs = np.linspace(0.05, 0.98, 30)*np.pi
        cases.append((zeros,poles,complex(1),freqs,fs))
    return cases


def scipy_per_root(zeros: np.ndarray, poles: np.ndarray,
                   frequencies: np.ndarray, fs: float) -> np.ndarray:
    """Sum SciPy one-root group delays, with the ZPK net z-power correction."""
    result = np.full(frequencies.shape, float(len(poles)-len(zeros)))
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        for root in zeros:
            _, gd = signal.group_delay((np.array([1.0, -root]),
                                        np.array([1.0])), w=frequencies, fs=fs)
            result += gd
        for root in poles:
            _, gd = signal.group_delay((np.array([1.0]),
                                        np.array([1.0, -root])), w=frequencies, fs=fs)
            result += gd
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable',type=Path)
    args = parser.parse_args()
    cases = build_cases()
    payload = [str(len(cases))]
    for z,p,k,w,fs in cases:
        # kind=0 chooses the real-roots overload in the Fortran driver.
        real_roots = np.all(np.imag(z) == 0) and np.all(np.imag(p) == 0) and k.imag == 0
        payload.append(f'{len(z)} {len(p)} {len(w)} {int(not real_roots)} '
                       f'{fs:.17e} {k.real:.17e} {k.imag:.17e}')
        payload.extend(f'{r.real:.17e} {r.imag:.17e}' for r in z)
        payload.extend(f'{r.real:.17e} {r.imag:.17e}' for r in p)
        payload.append(' '.join(f'{v:.17e}' for v in w))
    run = subprocess.run([str(args.executable)],input='\n'.join(payload)+'\n',
                         text=True,capture_output=True,check=True)
    lines = iter(run.stdout.splitlines())
    biggest_error = 0.0
    npoints = 0
    direct_count = 0
    for ix,(z,p,k,w,fs) in enumerate(cases):
        status_str,n_str = next(lines).split()
        assert int(status_str) == 0, (ix,'status',status_str)
        assert int(n_str) == len(w), (ix,'count',n_str)
        output = np.asarray([[float(v) for v in next(lines).split()] for _ in w])
        np.testing.assert_allclose(output[:,0],w,rtol=1e-15,atol=1e-11)
        independent = scipy_per_root(z,p,w,fs)
        np.testing.assert_allclose(output[:,1],independent,rtol=3e-11,atol=3e-10,
                                   err_msg=f'case {ix} (nz={len(z)}, np={len(p)})')
        biggest_error = max(biggest_error,float(np.max(abs(output[:,1]-independent))))
        # For short *equal-degree* filters, also check SciPy after zpk2tf.
        if len(z)==len(p) and 0 < len(z) <= 4 and ix % 2 == 0:
            b,a = signal.zpk2tf(z,p,k)
            with warnings.catch_warnings():
                warnings.simplefilter('error', RuntimeWarning)
                _, ref = signal.group_delay((b,a),w=w,fs=fs)
            np.testing.assert_allclose(output[:,1],ref,rtol=4e-10,atol=5e-9)
            direct_count += 1
        npoints += len(w)
    assert next(lines,None) is None, 'unread Fortran output'
    print(f'PASS {len(cases)} ZPK group-delay cases, {npoints} frequency values, '
          f'{direct_count} SciPy zpk2tf/group_delay cross-checks; '
          f'max absolute difference {biggest_error:.3e}')


if __name__ == '__main__':
    main()
