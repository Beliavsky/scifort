"""Compare SciFort Slepian/DPSS tapers with SciPy's tridiagonal implementation.

Usage: python test/compare_signal_dpss_scipy.py [gfortran]
"""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from scipy.signal.windows import dpss

ROOT = Path(__file__).resolve().parents[1]
SOURCES = [ROOT / 'src' / s for s in ('kinds.f90', 'signal_dpss.f90')]
DRIVER = ROOT / 'reference' / 'signal_dpss_reference_driver.f90'


def main() -> None:
    compiler = sys.argv[1] if len(sys.argv) > 1 else os.environ.get('FC', 'gfortran')
    rng = np.random.default_rng(20261020)
    cases: list[tuple[int, float, int, bool, str]] = []
    for m in (2, 3, 4, 5, 8, 15, 16, 25, 32, 63, 64, 97, 128, 192):
        for nw in (0.5, 1.1, 2.5, 3.5, 5.5):
            if nw >= m / 2:
                continue
            for symmetric in (True, False):
                for norm in ('l2', 'approximate'):
                    k = min(m, 5, max(1, int(2 * nw)))
                    cases.append((m, nw, k, symmetric, norm))
    for _ in range(120):
        m = int(rng.integers(5, 100))
        nw = float(rng.uniform(0.1, m / 2 - 0.05))
        k = int(rng.integers(1, min(5, m) + 1))
        cases.append((m, nw, k, bool(rng.integers(2)),
                      'l2' if rng.integers(2) else 'approximate'))
    request = ''.join(f'{m} {nw:.17g} {k} {int(symmetric)} {int(norm == "approximate")}\n'
                      for m, nw, k, symmetric, norm in cases)
    with tempfile.TemporaryDirectory(prefix='scifort-dpss-') as directory:
        tmp = Path(directory)
        exe = tmp / ('dpss.exe' if os.name == 'nt' else 'dpss')
        subprocess.run([compiler, '-std=f2018', '-O2', '-fcheck=all',
                        '-Werror=implicit-interface', *map(str, SOURCES),
                        str(DRIVER), '-o', str(exe)], cwd=tmp, check=True)
        proc = subprocess.run([str(exe)], input=request, capture_output=True,
                              text=True, check=True)
    lines = iter(proc.stdout.splitlines())
    peak_error = ratio_error = 0.0
    for i, (m, nw, k, symmetric, norm) in enumerate(cases):
        header = next(lines).split()
        assert header == ['CASE', str(m), str(k), '0'], (i, header)
        ratio = np.fromstring(next(lines), sep=' ')
        tapers = np.array([np.fromstring(next(lines), sep=' ') for _ in range(k)])
        reference, ref_ratio = dpss(m, nw, k, sym=symmetric,
                                    norm=2 if norm == 'l2' else norm, return_ratios=True)
        if tapers.shape != reference.shape:
            raise AssertionError((i, tapers.shape, reference.shape))
        err = float(np.max(np.abs(tapers - reference)))
        ratio_err = float(np.max(np.abs(ratio - ref_ratio)))
        peak_error = max(peak_error, err)
        ratio_error = max(ratio_error, ratio_err)
        if err > 1e-7 or ratio_err > 1e-8:
            raise AssertionError(f'case {i}, (M,NW,K,sym,norm)={(m,nw,k,symmetric,norm)} '
                                 f'window error {err:.3g}, ratio error {ratio_err:.3g}')
    print(f'SciPy DPSS comparisons: {len(cases)} cases passed, '
          f'peak window error={peak_error:.3g}, ratio error={ratio_error:.3g}')


if __name__ == '__main__':
    main()
