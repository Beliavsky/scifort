#!/usr/bin/env python3
"""Check axis-selective flood fill and binary hole filling against SciPy slices."""
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import ndimage as ndi

ROOT = Path(__file__).resolve().parents[2]
SOURCES = [ROOT / p for p in (
    'src/kinds.f90', 'src/signal_spline.f90', 'src/signal_spline_nd.f90',
    'src/ndimage_morphology.f90', 'src/ndimage_distance.f90',
    'src/ndimage_chamfer.f90', 'src/ndimage_measurements.f90',
    'src/ndimage_watershed.f90', 'src/ndimage.f90',
    'tools/ndimage_reference/propagation_axes_driver.f90',
)]


def main():
    i, j, k = np.mgrid[1:8, 1:5, 1:8]
    seeds = (7*i+11*j+3*k)%23 == 0
    mask = (i+3*j+5*k)%13 != 0
    image = np.ones((7, 4, 7), dtype=bool)
    image[1:6,:,1:6] = False
    image[3,:,3] = True
    image[0,:,2] = False
    footprint = ndi.generate_binary_structure(2, 1)
    with tempfile.TemporaryDirectory(prefix='scifort-selected-binary-') as td:
        tmp = Path(td)
        args = ['gfortran','-std=f2018','-O0','-fcheck=all',
                '-ffpe-trap=invalid,zero,overflow','-Werror=implicit-interface',
                '-J',str(tmp),'-I',str(tmp)]
        objs = []
        for src in SOURCES[:-1]:
            obj = tmp / (src.stem+'.o')
            subprocess.run([*args,'-c',str(src),'-o',str(obj)],
                           cwd=tmp,check=True,capture_output=True,text=True)
            objs.append(str(obj))
        subprocess.run([*args,*objs,str(SOURCES[-1]),'-o',str(tmp/'driver')],
                       cwd=tmp,check=True,capture_output=True,text=True)
        subprocess.run([str(tmp/'driver')],cwd=tmp,check=True,
                       capture_output=True,text=True)
        payload = np.fromfile(tmp/'propagation_axes.bin',dtype=np.int32)
    assert payload.size==6*seeds.size, (payload.size, seeds.size)
    results = [(r.reshape(seeds.shape,order='F')!=0) for r in payload.reshape(6,-1)]
    expected0=np.zeros_like(seeds)
    expected1=np.zeros_like(seeds)
    expected3=np.zeros_like(seeds)
    expected4=np.zeros_like(seeds)
    for j in range(seeds.shape[1]):
        expected0[:,j,:]=ndi.binary_propagation(seeds[:,j,:],structure=footprint,mask=mask[:,j,:])
        expected3[:,j,:]=ndi.binary_fill_holes(image[:,j,:],structure=footprint)
    for i in range(seeds.shape[0]):
        for k in range(seeds.shape[2]):
            expected1[i,:,k]=ndi.binary_propagation(seeds[i,:,k],mask=mask[i,:,k])
            expected4[i,:,k]=ndi.binary_fill_holes(image[i,:,k])
    expected=[expected0,expected1,seeds,expected3,expected4,image]
    for case,(got,ref) in enumerate(zip(results,expected)):
        np.testing.assert_array_equal(got,ref,err_msg=f'segmentation case {case}')
    print(f'PASS: {payload.size} selected-axis binary segmentation values')


if __name__=='__main__':
    main()
