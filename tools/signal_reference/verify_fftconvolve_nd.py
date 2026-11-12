#!/usr/bin/env python3
"""Compare rank-2/rank-3 FFT convolution with SciPy on selected axes."""
from __future__ import annotations
import pathlib
import subprocess
import tempfile
import numpy as np
from scipy import signal

ROOT = pathlib.Path(__file__).resolve().parents[2]
RNG = np.random.default_rng(20261018)


def encode_array(arr):
    return [f'{z.real:.17g} {z.imag:.17g}' for z in arr.flatten(order='F')]


def main():
    lines, cases = [], []
    for rank in (2, 3):
        for kind in range(4):
            for sample in range(125):
                axes = [axis for axis in range(rank) if (sample + (axis + 1)*11) % (axis + 3) != 0]
                if not axes:
                    axes = [sample % rank]
                mode = ('full','same','valid')[sample % 3]
                n = list(RNG.integers(1, 7, size=rank))
                m = list(RNG.integers(1, 7, size=rank))
                if mode == 'valid':
                    if sample % 2:
                        for ax in axes:
                            n[ax] = max(n[ax],m[ax])
                    else:
                        for ax in axes:
                            m[ax] = max(n[ax],m[ax])
                for axis in range(rank):
                    if axis not in axes:
                        n[axis],m[axis] = ((1,m[axis]) if sample % 2 else (n[axis],1))
                a = RNG.normal(size=n)
                b = RNG.normal(size=m)
                if kind in (1,3):
                    a = a + 1j*RNG.normal(size=n)
                if kind in (1,2):
                    b = b + 1j*RNG.normal(size=m)
                pad_n = n + [1]*(3-rank)
                pad_m = m + [1]*(3-rank)
                mask = sum(1 << ax for ax in axes)
                lines.append(' '.join(map(str, (rank,kind,('full','same','valid').index(mode)+1,mask,*pad_n,*pad_m))))
                lines.extend(encode_array(a))
                lines.extend(encode_array(b))
                reference = signal.fftconvolve(a,b,mode=mode,axes=axes)
                cases.append((rank,kind,mode,axes,n,m,reference))
    with tempfile.TemporaryDirectory() as tmp:
        work = pathlib.Path(tmp)
        cmd = ['gfortran','-std=f2018','-O2','-Werror=implicit-interface',f'-J{work}',f'-I{work}']
        sources = ['kinds','signal_dpss', 'signal_windows_extra', 'signal_spectral','signal_convolve','signal_fftconvolve']
        cmd += [str(ROOT/'src'/f'{source}.f90') for source in sources]
        cmd += [str(ROOT/'tools/signal_reference/fftconvolve_nd_driver.f90'),'-o',str(work/'driver')]
        subprocess.run(cmd,cwd=work,check=True)
        proc = subprocess.run([str(work/'driver')],input='\n'.join(lines)+'\n',
                              text=True,capture_output=True,check=True,timeout=120)
        result = iter(proc.stdout.splitlines())
        worst = 0.0
        for number,(rank,kind,mode,axes,n,m,reference) in enumerate(cases):
            vals = list(map(int,next(result).split()))
            status, shape = vals[0], tuple(vals[1:])
            if status != 0 or shape != reference.shape:
                raise AssertionError(f'case {number} rank {rank} {mode} axes {axes} shapes {n}, {m}: '
                                     f'status={status}, shape={shape}, expected={reference.shape}')
            values = np.array([complex(*map(float,next(result).split()))
                               for _ in range(reference.size)]).reshape(shape,order='F')
            error = np.max(np.abs(values-reference))
            worst = max(worst,error)
            if not np.allclose(values,reference,atol=1e-10,rtol=2e-11):
                raise AssertionError(f'case {number} max error {error}')
        print(f'PASS: {len(cases)} rank-2/rank-3 FFT convolution cases; max abs error {worst:.3g}')


if __name__ == '__main__':
    main()
