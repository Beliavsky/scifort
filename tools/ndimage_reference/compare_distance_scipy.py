#!/usr/bin/env python3
"""Compare exact separable Euclidean distance transforms against SciPy."""
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import ndimage as ndi

ROOT = Path(__file__).resolve().parents[2]


def main():
    m1 = np.asarray([(3*i+2)%5>1 for i in range(1,8)])
    m2 = np.fromfunction(lambda i,j:(7*(i+1)+3*(j+1)+(i+1)*(j+1))%11>2,
                         (6,5),dtype=int)
    m3 = np.fromfunction(lambda i,j,k:(13*(i+1)+5*(j+1)+3*(k+1)+(i+1)*(k+1))%13>1,
                         (4,5,3),dtype=int)
    m15 = np.ones((2,3)+(1,)*13,dtype=bool,order='F')
    m15[(0,1)+(0,)*13] = False
    sp15 = [1.7,.3]+[1.0]*13
    cases = [
        (m1,None), (m1,[2.3]), (np.ones(7,bool),None),
        (np.zeros(7,bool),None),
        (m2,None), (m2,[1.5,.6]),
        (m3,None), (m3,[.9,1.3,2.2]),
        (np.ones((6,5),bool),[1.5,.6]),
        (np.ones((4,5,3),bool),[.9,1.3,2.2]),
        (m15,None), (m15,sp15),
    ]
    largest = 0.0
    with tempfile.TemporaryDirectory(prefix='scifort-distance-') as tmp:
        tmp = Path(tmp)
        subprocess.run(['gfortran','-std=f2018','-O2','-fcheck=all',
                        '-Werror=implicit-interface','-J',str(tmp),'-I',str(tmp),
                        *(str(ROOT/f) for f in ('src/kinds.f90','src/signal_spline.f90',
                            'src/signal_spline_nd.f90','src/ndimage_morphology.f90',
                            'src/ndimage_distance.f90','src/ndimage_chamfer.f90','src/ndimage_measurements.f90','src/ndimage_watershed.f90',
               'src/ndimage.f90',
                            'tools/ndimage_reference/distance_reference_driver.f90')),
                        '-o',str(tmp/'driver')],cwd=tmp,check=True)
        subprocess.run([str(tmp/'driver')],cwd=tmp,check=True)
        with (tmp/'distance_reference_results.bin').open('rb') as f:
            for mask,spacing in cases:
                want = ndi.distance_transform_edt(mask, sampling=spacing)
                got = np.frombuffer(f.read(8*mask.size),dtype=np.float64).reshape(mask.shape,order='F')
                error = np.max(np.abs(got-want))
                largest = max(largest,float(error))
                if error>2e-12:
                    index=np.unravel_index(np.argmax(np.abs(got-want)),mask.shape)
                    raise AssertionError((mask.shape, spacing, error, index, got[index], want[index]))
            if f.read(1):
                raise AssertionError('trailing reference bytes')
    print(f'PASS: {len(cases)} Euclidean distance-transform SciPy cases, max error={largest:.3g}')


if __name__=='__main__':
    main()
