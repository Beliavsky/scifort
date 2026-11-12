#!/usr/bin/env python3
"""Cross-validate morphology composites and geodesic propagation."""
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import ndimage as ndi

ROOT = Path(__file__).resolve().parents[2]


def main():
    m2 = np.fromfunction(lambda i,j: (7*(i+1)+11*(j+1)+(i+1)*(j+1)) % 8 < 3,
                         (5,6),dtype=int)
    m3 = np.fromfunction(lambda i,j,k: (3*(i+1)+7*(j+1)+5*(k+1)+
                           (i+1)*(j+1)+(k+1)*(j+1)) % 9 < 4, (4,5,3), dtype=int)
    ij = np.indices(m2.shape)+1
    ijk = np.indices(m3.shape)+1
    s2 = m2 & ((3*ij[0]+2*ij[1])%4==0)
    a2 = (2*ij[0]+5*ij[1])%7<5
    s3 = m3 & ((ijk[0]+3*ijk[1]+7*ijk[2])%5==0)
    a3 = (2*ijk[0]+ijk[1]+3*ijk[2])%7<5
    h2 = np.ones((7,8),dtype=bool,order='F')
    h2[2:5,2:6] = False
    h2[3,4] = True
    h2[0,3] = h2[1,3] = False
    h3 = np.ones((5,5,4),dtype=bool,order='F')
    h3[1:4,1:4,1:3] = False
    h3[0,0,0] = False
    checks = 0
    with tempfile.TemporaryDirectory(prefix='scifort-ext-morph-') as tmp:
        tmp = Path(tmp)
        subprocess.run(['gfortran','-std=f2018','-O0','-fcheck=all',
                        '-Werror=implicit-interface','-J',str(tmp),'-I',str(tmp),
                        *(str(ROOT/f) for f in ('src/kinds.f90','src/signal_spline.f90',
                            'src/signal_spline_nd.f90','src/ndimage_morphology.f90',
                            'src/ndimage_distance.f90','src/ndimage_chamfer.f90', 'src/ndimage_measurements.f90','src/ndimage_watershed.f90',
               'src/ndimage.f90',
                            'tools/ndimage_reference/extended_morphology_driver.f90')),
                        '-o',str(tmp/'driver')],cwd=tmp,check=True)
        subprocess.run([str(tmp/'driver')],cwd=tmp,check=True)
        with (tmp/'extended_morphology_results.bin').open('rb') as f:
            for x,seed,allowed,holes in ((m2,s2,a2,h2),(m3,s3,a3,h3)):
                for conn in range(1,x.ndim+1):
                    st = ndi.generate_binary_structure(x.ndim,conn)
                    cases = []
                    for it in (1,2):
                        for border in (0,1):
                            cases.extend([
                                ('opening',ndi.binary_opening(x,structure=st,
                                   iterations=it,border_value=border)),
                                ('closing',ndi.binary_closing(x,structure=st,
                                   iterations=it,border_value=border)),
                            ])
                    cases.extend([
                        ('fill_holes',ndi.binary_fill_holes(holes,structure=st)),
                        ('propagate',ndi.binary_propagation(seed,structure=st,mask=allowed)),
                    ])
                    for name,expected in cases:
                        got = np.frombuffer(f.read(4*expected.size),dtype=np.int32).reshape(
                              expected.shape,order='F').astype(bool)
                        if not np.array_equal(got,expected):
                            raise AssertionError((x.shape,conn,name,np.argwhere(got!=expected)))
                        checks += 1
            if f.read(1):
                raise AssertionError('trailing bytes')
    print(f'PASS: {checks} opening/closing/hole-filling/propagation SciPy cases')


if __name__=='__main__':
    main()
