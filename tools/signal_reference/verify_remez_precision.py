#!/usr/bin/env python3
"""Exercise Remez's 80-bit/real64 fallback branches on GNU Fortran.

Forcing an unavailable 128-bit kind in a temporary source copy avoids
requiring a compiler that genuinely lacks binary128. This is a diagnostic
portability test; available real kinds and success counts vary by compiler.
"""
from __future__ import annotations
import tempfile
import subprocess
from pathlib import Path
import numpy as np
from verify_remez import make_cases

ROOT = Path(__file__).resolve().parents[2]


def run_variant(source: str, cases: list, tmp: Path, mode: str) -> tuple[int, int]:
    patched = source.replace('qp = selected_real_kind(30, 400)', 'qp = -1')
    if mode == 'real64':
        patched = patched.replace('ep = selected_real_kind(18, 300)', 'ep = -1')
    assert patched != source
    src = tmp / 'signal_remez.f90'
    src.write_text(patched)
    exe = tmp / 'remez'
    cmd = ['gfortran', '-std=f2018', '-O2', '-fcheck=all', '-Werror=implicit-interface',
           '-J', str(tmp), '-I', str(tmp), str(ROOT/'src/kinds.f90'),str(src),
           str(ROOT/'tools/signal_reference/remez_driver.f90'),'-o',str(exe)]
    subprocess.run(cmd,cwd=tmp,check=True)
    records = [str(len(cases))]
    for n,b,d,w,k,fs in cases:
        records.append(f'{n} {len(d)} {dict(bandpass=1,hilbert=2,differentiator=3)[k]} 16 {fs:.17g}')
        records.extend(' '.join(f'{v:.17g}' for v in arr) for arr in (b,d,w))
    result = subprocess.run([str(exe)],input='\n'.join(records)+'\n',text=True,
                            capture_output=True,check=True,cwd=tmp,timeout=180)
    rows = iter(result.stdout.splitlines())
    successes = 0
    for n,*_ in cases:
        st,_,ntaps,_ = next(rows).split()
        st,ntaps = int(st),int(ntaps)
        if st == 0:
            h = np.fromstring(next(rows),sep=' ')
            if ntaps != n or len(h) != n or not np.all(np.isfinite(h)):
                raise AssertionError(f'{mode}: malformed coefficient result')
            successes += 1
        elif st != 2 or ntaps != 0:
            raise AssertionError(f'{mode}: unexpected status {st}, count {ntaps}')
    if next(rows,None) is not None:
        raise AssertionError('Unexpected trailing results')
    return successes, len(cases)


def main():
    src = (ROOT/'src/signal_remez.f90').read_text()
    cases = make_cases()
    with tempfile.TemporaryDirectory(prefix='scifort-remez-precision-') as path:
        tmp = Path(path)
        eighty = run_variant(src,cases,tmp,'extended-if-present')
        double = run_variant(src,cases,tmp,'real64')
    print(f'Extended-kind fallback: {eighty[0]}/{eighty[1]} designs returned coefficients')
    print(f'Real64-only fallback:   {double[0]}/{double[1]} designs returned coefficients')
    if eighty[0] < double[0]:
        raise AssertionError('Higher-precision fallback unexpectedly lost successful designs')


if __name__ == '__main__':
    main()
