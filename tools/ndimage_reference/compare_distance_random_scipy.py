#!/usr/bin/env python3
"""Stress-test distance transforms using fixed-seed random binary grids."""
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import ndimage as ndi

ROOT=Path(__file__).resolve().parents[2]


def main():
    rng=np.random.default_rng(20261022)
    cases=[]
    for ndim in (1,2,3):
        for n in range(45):
            shape=tuple(int(t) for t in rng.integers(1,9,size=ndim))
            if n%10==0:
                mask=np.ones(shape,dtype=bool,order='F')
            elif n%10==1:
                mask=np.zeros(shape,dtype=bool,order='F')
            else:
                mask=np.asfortranarray(rng.uniform(size=shape)>.20+n%4*.15)
            spacing=10.0**rng.uniform(-2,2,size=ndim)
            cases.append((mask,spacing))
    maxdiff=0.0
    with tempfile.TemporaryDirectory(prefix='scifort-random-distance-') as tmp:
        tmp=Path(tmp)
        with (tmp/'random_distance_input.bin').open('wb') as f:
            f.write(np.array([len(cases)],dtype=np.int32).tobytes())
            for mask,sp in cases:
                d=np.ones(3,dtype=np.int32)
                d[:mask.ndim]=mask.shape
                s=np.ones(3,dtype=np.float64)
                s[:mask.ndim]=sp
                f.write(np.array([mask.ndim],dtype=np.int32).tobytes())
                f.write(d.tobytes())
                f.write(s.tobytes())
                f.write(mask.astype(np.int32).ravel(order='F').tobytes())
        subprocess.run(['gfortran','-std=f2018','-O2','-fcheck=all',
                        '-ffpe-trap=invalid,zero,overflow','-Werror=implicit-interface',
                        '-J',str(tmp),'-I',str(tmp),
                        *(str(ROOT/f) for f in ('src/kinds.f90','src/signal_spline.f90',
                            'src/signal_spline_nd.f90','src/ndimage_morphology.f90',
                            'src/ndimage_distance.f90','src/ndimage_chamfer.f90','src/ndimage_measurements.f90','src/ndimage_watershed.f90',
               'src/ndimage.f90',
                            'tools/ndimage_reference/random_distance_driver.f90')),
                        '-o',str(tmp/'driver')],cwd=tmp,check=True)
        subprocess.run([str(tmp/'driver')],cwd=tmp,check=True)
        with (tmp/'random_distance_output.bin').open('rb') as f:
            for n,(mask,sampling) in enumerate(cases):
                got=np.frombuffer(f.read(mask.size*8),dtype=np.float64).reshape(mask.shape,order='F')
                expected=ndi.distance_transform_edt(mask,sampling=sampling)
                difference=np.max(np.abs(expected-got))
                maxdiff=max(maxdiff,float(difference))
                if difference>2e-11:
                    at=np.unravel_index(np.argmax(np.abs(expected-got)),mask.shape)
                    raise AssertionError((n,mask.shape,sampling,at,got[at],expected[at]))
            if f.read(1):
                raise AssertionError('trailing bytes')
    print(f'PASS: {len(cases)} randomized EDT SciPy comparisons; max error {maxdiff:.3g}')


if __name__=='__main__':
    main()
