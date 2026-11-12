#!/usr/bin/env python3
"""Compare arbitrary-rank Fortran Hilbert transforms with SciPy 1.17.

Usage: python test/compare_signal_hilbert_nd_scipy.py /path/to/fortran_driver
"""
import pathlib
import subprocess
import sys
import tempfile

import numpy as np
from scipy.signal import hilbert, hilbert2


def make_cases(seed=20261027):
    rng = np.random.default_rng(seed)
    cases = []
    for rank in (1, 2, 3, 5, 15):
        for idx in range(100):
            mode = 1 if rank == 1 or idx % 2 == 0 else 2
            if rank <= 3:
                dims = rng.integers(2, 9, size=rank).tolist()
            elif rank == 5:
                dims = rng.integers(1, 5, size=rank).tolist()
            else:
                dims = [1] * 15
                for j in rng.choice(rank, size=3, replace=False):
                    dims[j] = int(rng.integers(2, 5))
            selected = rng.choice(rank, size=mode, replace=False).tolist()
            for j in selected:
                dims[j] = int(rng.integers(1, 15))
            axes = [selected[0], selected[-1] if mode == 2 else -1]
            output_dims = dims.copy()
            for j in selected:
                output_dims[j] = int(rng.integers(1, 18))
            signed_axes = [j + 1 if idx % 3 != 0 else j - rank for j in selected]
            if mode == 1:
                signed_axes.append(0)
            x = rng.standard_normal(dims)
            if mode == 1:
                expected = hilbert(x, N=output_dims[selected[0]], axis=selected[0])
            else:
                moved = np.moveaxis(x,selected,[-2,-1])
                expected = np.moveaxis(
                    hilbert2(moved, N=[output_dims[j] for j in selected]),
                    [-2,-1], selected,
                )
            cases.append((rank, mode, dims, output_dims, signed_axes, x, expected))
    return cases


def main():
    driver = pathlib.Path(sys.argv[1]).resolve()
    cases = make_cases()
    with tempfile.TemporaryDirectory() as d:
        input_path = pathlib.Path(d, 'input.txt')
        output_path = pathlib.Path(d, 'output.txt')
        with input_path.open('w') as f:
            f.write(f'{len(cases)}\n')
            for rank, mode, dims, outdims, axes, x, _ in cases:
                meta = [rank, mode, *axes, *dims, *([1]*(15-rank)),
                        *outdims, *([1]*(15-rank))]
                f.write(' '.join(map(str,meta))+'\n')
                f.write(' '.join(format(float(v),'.17g') for v in x.ravel(order='F'))+'\n')
        subprocess.run([str(driver),str(input_path),str(output_path)],check=True)
        lines=output_path.read_text().splitlines()
        assert len(lines) == 2*len(cases), (len(lines),len(cases))
        maximum=0.0
        maximum_scaled=0.0
        compared=0
        for i,(_, _, dims, outdims, _, _, expected) in enumerate(cases):
            status=int(lines[2*i])
            if status != 0:
                raise AssertionError(f'case {i}: Fortran status={status}')
            flat = np.fromstring(lines[2*i+1],sep=' ')
            expected_count=int(np.prod(outdims))
            if flat.size != 2*expected_count:
                raise AssertionError(f'case {i}: unexpected number of outputs {flat.size}')
            y=(flat[::2]+1j*flat[1::2]).reshape(outdims,order='F')
            err=float(np.max(np.abs(y-expected)))
            maximum=max(maximum,err)
            maximum_scaled=max(maximum_scaled,err / max(1.,float(np.max(abs(expected)))))
            compared+=expected_count
            if not np.allclose(y,expected,rtol=3e-11,atol=3e-11):
                loc=np.unravel_index(np.argmax(abs(y-expected)),y.shape)
                raise AssertionError(
                    f'case {i} dims={dims}, result dims={outdims}, max_err={err}, location={loc}, '
                    f'got={y[loc]}, expected={expected[loc]}')
    print(f'{len(cases)} SciPy Hilbert cases passed; {compared} complex values; '
          f'max abs difference {maximum:.3e}, scaled {maximum_scaled:.3e}')


if __name__ == '__main__':
    main()
