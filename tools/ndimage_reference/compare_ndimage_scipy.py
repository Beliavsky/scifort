#!/usr/bin/env python3
"""Cross-check SciFort ndimage kernels with SciPy (binary64 Fortran ordering)."""
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import ndimage as ndi

ROOT = Path(__file__).resolve().parents[2]
MODES = ('reflect','mirror','nearest','wrap','constant','grid-constant',
         'grid-mirror','grid-wrap')


def main():
    x = np.empty((4, 5, 3), order='F')
    for k in range(1, 4):
        for j in range(1, 6):
            for i in range(1, 5):
                x[i-1,j-1,k-1] = np.sin(.31*i+.42*j)+np.cos(.23*i-.19*j)+.11*k
    z = x + 1j*(.37*x+.21)
    count = 0
    maxima = {}
    with tempfile.TemporaryDirectory(prefix='scifort-ndimage-') as tmp:
        tmp = Path(tmp)
        cmd = ['gfortran', '-std=f2018', '-O2', '-fcheck=all',
               '-Werror=implicit-interface', '-J', str(tmp), '-I', str(tmp)]
        cmd += [str(ROOT / f) for f in ('src/kinds.f90','src/signal_spline.f90',
                 'src/signal_spline_nd.f90','src/ndimage_morphology.f90','src/ndimage_distance.f90','src/ndimage_chamfer.f90','src/ndimage_measurements.f90','src/ndimage_watershed.f90',
               'src/ndimage.f90',
                 'tools/ndimage_reference/ndimage_reference_driver.f90')]
        cmd += ['-o', str(tmp/'driver')]
        subprocess.run(cmd, cwd=tmp, check=True)
        subprocess.run([str(tmp/'driver')], cwd=tmp, check=True, capture_output=True)
        with (tmp/'ndimage_reference_results.bin').open('rb') as f:
            for mode in MODES:
                cases = [
                    ('gaussian', ndi.gaussian_filter(x, [.6,1.3,.9], mode=mode, cval=-1.75)),
                    ('gaussian_complex', ndi.gaussian_filter(z, [.6,1.3,.9], mode=mode,
                                                            cval=-1.75+.6j)),
                    ('derivatives', ndi.gaussian_filter(x, [1.1,.8,1.5], order=[1,0,2],
                                                  mode=mode,cval=-1.75,radius=[3,2,4])),
                    ('one_axis', ndi.gaussian_filter1d(x,1.2,axis=1,order=3,mode=mode,
                                                      cval=-1.75,radius=4)),
                    ('gradient', ndi.gaussian_gradient_magnitude(x,[1.0,.8,1.2],
                                                                 mode=mode,cval=-1.75)),
                    ('gaussian_laplace', ndi.gaussian_laplace(x,[1.0,.8,1.2],
                                                             mode=mode,cval=-1.75)),
                    ('uniform', ndi.uniform_filter(x,size=[2,3,4],mode=mode,cval=-1.75)),
                    ('uniform_complex', ndi.uniform_filter(z,size=[2,3,4],mode=mode,
                                                            cval=-1.75+.6j)),
                    ('minimum', ndi.minimum_filter(x,size=[2,3,4],mode=mode,cval=-1.75)),
                    ('maximum', ndi.maximum_filter(x,size=[2,3,4],mode=mode,cval=-1.75)),
                    ('sobel', ndi.sobel(x,axis=1,mode=mode,cval=-1.75)),
                    ('laplace', ndi.laplace(x,mode=mode,cval=-1.75)),
                    ('gaussian_scalar', ndi.gaussian_filter(x,1.0,mode=mode)),
                    ('uniform_scalar', ndi.uniform_filter(x,3,mode=mode)),
                    ('minimum_scalar', ndi.minimum_filter(x,3,mode=mode)),
                    ('maximum_scalar', ndi.maximum_filter(x,3,mode=mode)),
                ]
                for name,want in cases:
                    dtype=np.complex128 if np.iscomplexobj(want) else np.float64
                    n=x.size*np.dtype(dtype).itemsize
                    got=np.frombuffer(f.read(n),dtype=dtype).reshape(x.shape,order='F')
                    difference=np.max(np.abs(got-want))
                    maxima[name]=max(maxima.get(name,0.0),float(difference))
                    if difference>5e-11:
                        at=np.unravel_index(np.argmax(np.abs(got-want)),x.shape)
                        raise AssertionError((mode,name,difference,at,got[at],want[at]))
                    count+=1
            if f.read(1):
                raise AssertionError('trailing binary reference output')
    print(f'PASS: {count} SciPy ndimage reference comparisons')
    print('Largest absolute errors:',maxima)

if __name__=='__main__':
    main()
