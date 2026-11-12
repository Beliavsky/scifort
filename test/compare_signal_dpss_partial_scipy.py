"""Independent SciPy comparison of selected-eigenpair DPSS and auto dispatch.

Runs GNU Fortran in either checked or optimized mode; SciPy is used only by
this validation script and is not a SciFort dependency.
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
SOURCE = [ROOT / 'src' / name for name in ('kinds.f90', 'signal_dpss.f90')]
DRIVER = ROOT / 'reference' / 'signal_dpss_partial_reference_driver.f90'


def main() -> None:
    compiler = os.environ.get('FC', 'gfortran')
    mode = sys.argv[1] if len(sys.argv) > 1 else 'checked'
    if mode not in ('checked', 'optimized'):
        raise SystemExit('expected checked or optimized')
    rng = np.random.default_rng(20261027)
    # M, NW, K, sym, norm, solver (0=auto, 1=QL, 2=partial).
    cases: list[tuple[int, float, int, int, int, int]] = []
    for m in (2, 3, 4, 9, 31, 32, 63, 64, 97, 128, 255, 256,
              257, 511, 512, 1024, 2048, 4096):
        for nw in (0.25, 0.75, 1.5, 3.5, 5.5, 0.3 * (m-1)):
            if nw >= m/2:
                continue
            k = min(m, 10, max(1, int(2*nw)+1))
            cases.append((m, nw, k, int(rng.integers(2)),
                          int(rng.integers(3)), 2))
    for _ in range(120):
        m = int(rng.integers(3, 1900))
        nw = float(rng.uniform(0.08, max(0.081, m/2-0.05)))
        k = int(rng.integers(1, min(10,m)+1))
        cases.append((m, nw, k, int(rng.integers(2)),
                      int(rng.integers(3)), 2))
    # Long records: the former QL path would allocate n*n doubles.
    cases.extend([(4096, 3.5, 6, 1, 0, 0),
                  (8192, 3.5, 6, 1, 0, 2),
                  (16384, 4.5, 6, 1, 0, 0)])
    # Verify automatic dispatch is equivalent to explicit QL on small inputs.
    cases.extend([(64, 3.5, 5, 1, 0, 0), (128, 4.0, 6, 0, 2, 0)])
    request = ''.join(' '.join(map(str, case))+'\n' for case in cases)
    with tempfile.TemporaryDirectory(prefix='scifort-dpss-partial-') as directory:
        tmp = Path(directory)
        exe = tmp / ('ref.exe' if os.name == 'nt' else 'ref')
        flags = ['-std=f2018', '-Werror=implicit-interface']
        flags += ['-O0', '-g', '-fcheck=all'] if mode == 'checked' else ['-O3']
        subprocess.run([compiler, *flags, *map(str, SOURCE),
                        str(DRIVER), '-o', str(exe)], check=True, cwd=tmp)
        output = subprocess.run([str(exe)], input=request,
                                capture_output=True, text=True, check=True).stdout
    lines = iter(output.splitlines())
    max_taper = max_ratio = 0.0
    count = 0
    for case in cases:
        m, nw, k, sym, norm, solver = case
        header = next(lines).split()
        assert header == ['CASE', str(m), str(k), '0'], (case, header)
        ratios = np.fromstring(next(lines), sep=' ')
        values = np.array([np.fromstring(next(lines), sep=' ') for _ in range(k)])
        reference, ref_ratios = dpss(m, nw, k, sym=bool(sym),
                                      norm=(2, 'approximate', 'subsample')[norm],
                                      return_ratios=True)
        assert values.shape == (k, m), case
        err = float(np.max(np.abs(values-reference)))
        ratio_err = float(np.max(np.abs(ratios-ref_ratios)))
        max_taper = max(err, max_taper)
        max_ratio = max(ratio_err, max_ratio)
        if err > 1e-7 or ratio_err > 2e-9:
            raise AssertionError(f'{case}: taper {err:g}, ratio {ratio_err:g}')
        count += int(values.size)
    assert next(lines, None) is None
    print(f'PASS {mode}: {len(cases)} SciPy DPSS cases, {count} coefficients, '
          f'max taper error={max_taper:.3g}, max ratio error={max_ratio:.3g}')


if __name__ == '__main__':
    main()
