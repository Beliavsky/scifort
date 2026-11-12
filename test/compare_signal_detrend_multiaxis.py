"""Independent SciPy 1.17 comparison for real/complex rank-1/2/3 detrend.

    python test/compare_signal_detrend_multiaxis.py [gfortran]

Fortran axes and breakpoint positions are one-based. SciPy's are zero-based.
"""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from scipy.signal import detrend

ROOT = Path(__file__).resolve().parents[1]
SOURCES = [ROOT / 'src' / f'{s}.f90' for s in ('kinds', 'signal_detrend')]
DRIVER = ROOT / 'reference' / 'signal_detrend_multiaxis_reference_driver.f90'


def main():
    fc = sys.argv[1] if len(sys.argv) > 1 else os.getenv('FC', 'gfortran')
    rng = np.random.default_rng(20261020)
    cases = []
    records = []
    for case_id in range(720):
        rank = case_id % 3 + 1
        shape = tuple(int(rng.integers(2, 19)) for _ in range(rank))
        is_complex = (case_id // 3) % 2 == 1
        axis = int(rng.integers(0, rank))
        trend = ('linear', 'constant')[case_id % 5 == 0]
        data = rng.normal(size=shape)
        if is_complex:
            data = data + 1j * rng.normal(size=shape)
        # Stress narrow or duplicate piecewise breakpoints.
        bp = []
        if trend == 'linear' and case_id % 4 != 0:
            count = int(rng.integers(0, 7))
            bp = [int(rng.integers(0, shape[axis]+1)) for _ in range(count)]
        cases.append((data, axis, trend, bp))
        records.append(' '.join(map(str, [rank, int(is_complex), int(trend == 'constant'),
                                           axis+1, len(bp), *(shape+(1,)*(3-rank))])))
        if bp:
            records.append(' '.join(str(v+1) for v in bp))
        flat = data.ravel(order='F')
        records.append(' '.join(map(str, flat.real)))
        if is_complex:
            records.append(' '.join(map(str, flat.imag)))
    stdin = str(len(cases)) + '\n' + '\n'.join(records) + '\n'
    with tempfile.TemporaryDirectory(prefix='scifort-detrend-') as tmp:
        wd = Path(tmp)
        objects = []
        for src in SOURCES:
            obj = wd / (src.stem+'.o')
            subprocess.run([fc,'-O2','-std=f2018','-I',tmp,'-J',tmp,'-c',str(src),
                            '-o',str(obj)],check=True,cwd=wd)
            objects.append(str(obj))
        exe = wd / 'detrend_ref'
        subprocess.run([fc,'-O2','-std=f2018','-I',tmp,str(DRIVER),*objects,'-o',str(exe)],
                       check=True,cwd=wd)
        result = subprocess.run([str(exe)],input=stdin,text=True,capture_output=True,
                                check=True,cwd=wd)
    lines = result.stdout.splitlines()
    position = 0
    max_error = 0.0
    for i, (data,axis,trend,bp) in enumerate(cases):
        stat,num = map(int,lines[position].split())
        position += 1
        assert stat == 0 and num == data.size,(i,stat,num,data.shape)
        vals = np.array([[float(v) for v in s.split()] for s in lines[position:position+num]])
        position += num
        actual = vals[:,0] + 1j*vals[:,1] if np.iscomplexobj(data) else vals[:,0]
        actual = actual.reshape(data.shape,order='F')
        expected = detrend(data,axis=axis,type=trend,bp=(bp if bp else 0))
        np.testing.assert_allclose(actual,expected,rtol=2e-10,atol=5e-12,err_msg=f'case {i}')
        max_error = max(max_error,float(np.max(abs(actual-expected))))
    assert position == len(lines), (position,len(lines))
    print(f'PASS: {len(cases)} SciPy 1.17 detrend comparisons; max abs error {max_error:.3g}')


if __name__ == '__main__':
    main()
