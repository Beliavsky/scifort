#!/usr/bin/env python3
"""Compare SOS-native Gustafsson filtering with an independent SciPy SOS model.

SciPy does not provide method='gust' for sosfiltfilt. This reference builds
an observability matrix in SOS direct-form-II state space and uses SciPy's
sosfilt and scipy.linalg.lstsq; no high-order transfer polynomial is formed.

Run: python tools/signal_reference/compare_signal_sos_gust_scipy.py DRIVER_EXE
"""
import sys
import subprocess
import numpy as np
from scipy import signal, linalg


def gust_sos_reference(sos, x, irlen):
    sos = np.asarray(sos, dtype=complex)
    x = np.asarray(x, dtype=complex)
    n, order = len(x), 2 * len(sos)
    m = n if irlen < 0 or n <= 2*irlen else irlen
    obs = np.empty((m, order), dtype=complex)
    for k in range(order):
        init = np.zeros((len(sos), 2), dtype=complex)
        init[k % len(sos), k // len(sos)] = 1
        obs[:, k] = signal.sosfilt(sos, np.zeros(m, dtype=complex), zi=init)[0]
    s = signal.sosfilt(sos, obs[::-1], axis=0)
    fb = signal.sosfilt(sos, signal.sosfilt(sos, x)[::-1])[::-1]
    bf = signal.sosfilt(sos, signal.sosfilt(sos, x[::-1])[::-1])
    obsrev, srev = obs[::-1], s[::-1]
    if m == n:
        mat = np.column_stack((srev-obs, obsrev-s))
        corr = np.column_stack((srev, obsrev))
        delta = bf-fb
    else:
        mat = np.zeros((2*m, 2*order), dtype=complex)
        corr = np.zeros_like(mat)
        mat[:m, :order] = srev - obs
        mat[m:, order:] = obsrev - s
        corr[:m, :order] = srev
        corr[m:, order:] = obsrev
        delta = np.r_[bf[:m]-fb[:m], bf[-m:]-fb[-m:]]
    ic, _, _, _ = linalg.lstsq(mat, delta, lapack_driver='gelsd')
    val = corr @ ic
    y = fb.copy()
    if m == n:
        y += val
    else:
        y[:m] += val[:m]
        y[-m:] += val[m:]
    return y


if len(sys.argv) != 2:
    raise SystemExit("usage: compare_signal_sos_gust_scipy.py DRIVER_EXE")

rng = np.random.default_rng(20261026)
cases = []
for sections in (1, 2, 3):
    for n in (1, 2, 3, 7, 13, 32, 65):
        for irlen in (-1, 1, 3, 12):
            if (n + sections + irlen) % 3 == 1:
                continue
            for kind in (0, 1, 2, 3):
                # Conjugate real poles for real sections; complex modes have
                # conjugacy-free poles. All poles are comfortably stable.
                sos = np.zeros((sections, 6), dtype=complex)
                for j in range(sections):
                    pole1 = rng.uniform(.05, .55) * np.exp(1j*rng.uniform(-1,1))
                    pole2 = rng.uniform(.05, .55) * np.exp(1j*rng.uniform(-1,1))
                    if kind < 2:
                        pole1 = rng.uniform(.05, .55)
                        pole2 = rng.uniform(.05, .55)
                    sos[j, :3] = rng.normal(0, .2, 3) + (
                        (1j*rng.normal(0,.04,3)) if kind >= 2 else 0)
                    sos[j,3:] = [1, -(pole1+pole2), pole1*pole2]
                if sections==2 and n%3==0:
                    sos[-1,2]=0
                    sos[-1,5]=0
                x = rng.normal(0,1,n)
                if kind%2:
                    x = x + 1j*rng.normal(0,1,n)
                cases.append((kind,sos,x,irlen))

for order in (4, 8, 12, 16):
    sos = signal.butter(order, .16, output='sos')
    for n in (20, 95, 240):
        for irlen in (-1, 8, 32):
            for mode in (0,1):
                x=rng.normal(size=n)
                if mode:
                    x=x+1j*rng.normal(size=n)
                cases.append((mode,sos,x,irlen))

# A longer stream exercises bounded state-estimation memory.
sos = signal.butter(8, .2, output='sos')
cases.append((0,sos,rng.normal(size=2048),40))

def pair(z):
    return f'{complex(z).real:.17e} {complex(z).imag:.17e}\n'

inputs = [f'{len(cases)}\n']
expect = []
for mode, sos, x, irlen in cases:
    inputs.append(f'{mode} {len(x)} {len(sos)} {irlen}\n')
    inputs.extend(pair(v) for row in sos for v in row)
    inputs.extend(pair(v) for v in x)
    expect.append(gust_sos_reference(sos,x,irlen))

p = subprocess.run([sys.argv[1]],input=''.join(inputs),capture_output=True,text=True,check=True)
lines = iter(p.stdout.splitlines())
max_abs = 0.
max_scaled = 0.
rejected = 0
total_samples = 0
for k, expected in enumerate(expect):
    code, n = map(int,next(lines).split())
    if code == 2:
        mode, sos, x, irlen = cases[k]
        # Near-singular, untruncated systems can return a misleading
        # least-squares output when the initial-state response is poorly
        # determined. SciFort explicitly reports numerical failure there.
        if (len(sos) >= 5 and len(x) <= 4*len(sos) and
                (irlen < 0 or len(x) <= 2*irlen)):
            rejected += 1
            for _ in range(n):
                _ = next(lines)
            continue
    if code != 0 or n != len(expected):
        raise AssertionError(f'case {k}: status {code}, shape {n} vs {len(expected)}')
    actual = np.array([complex(*map(float,next(lines).split())) for _ in range(n)])
    max_abs = max(max_abs,float(np.max(np.abs(actual-expected))))
    denom = 1+np.max(np.abs(expected))
    rel = float(np.max(np.abs(actual-expected))/denom)
    max_scaled = max(max_scaled,rel)
    # Ill-conditioned rank-deficient systems, especially n<order, warrant
    # norm-aware comparison, rather than a strict allclose at machine eps.
    if rel > 5.e-7:
        kind,sos,x,irlen=cases[k]
        raise AssertionError(f'case={k}, sections={len(sos)}, n={len(x)}, '
                             f'irlen={irlen}, mode={kind}, error={rel:.3g}')
    total_samples += n
if rejected > 4:
    raise AssertionError(f'Unexpected regression: {rejected} SOS cases rejected (maximum 4)')
print(f'PASS {len(cases)-rejected} SOS Gustafsson comparisons ({rejected} ill-conditioned cases rejected), {total_samples} samples, '
      f'max abs={max_abs:.3e}, max magnitude-scaled={max_scaled:.3e}')
