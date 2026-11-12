#!/usr/bin/env python3
"""SciPy comparison of axis-selective rank filters and linear-time moving means."""
from pathlib import Path
from tempfile import TemporaryDirectory
import subprocess
import numpy as np
from scipy import ndimage as ndi

ROOT = Path(__file__).resolve().parents[2]
FILES = ['kinds', 'signal_spline', 'signal_spline_nd', 'ndimage_morphology',
         'ndimage_distance', 'ndimage_chamfer', 'ndimage_measurements',
         'ndimage_watershed', 'ndimage']

def main():
    a, b, c = np.indices((4, 5, 3)) + 1
    x = np.sin(.17*a + .31*b - .19*c) + .05*(a*a - 2*b + c)
    z = x + 1j*(.4*x + .2)
    fp = np.array([[True, True, False], [False, True, True]], order='F')
    order_cases = [
        ndi.median_filter(x,size=(4,3),mode='nearest',axes=(2,0)),
        ndi.median_filter(x,size=5,mode='mirror',origin=-1,axes=(0,)),
        ndi.median_filter(x,size=3,axes=()),
        ndi.rank_filter(x,-2,size=(2,4),mode='reflect',origin=(0,-1),axes=(0,2)),
        ndi.rank_filter(x,0,footprint=fp,origin=(0,0),axes=(2,0),mode='constant',cval=-1.5),
        ndi.percentile_filter(x,35.,size=(4,3),origin=(1,-1),axes=(1,0),mode='wrap'),
        ndi.percentile_filter(x,100.,footprint=fp,origin=(0,0),axes=(2,0),mode='mirror'),
        ndi.percentile_filter(x,80.,size=3,axes=()),
        ndi.median_filter(x,size=(2,3,4),mode='constant',cval=-3.5),
        ndi.rank_filter(x,-1,size=3,mode='nearest'),
        ndi.percentile_filter(x,75.,size=2,mode='mirror'),
        ndi.rank_filter(x,0,size=3,axes=(1,),mode='reflect'),
        ndi.rank_filter(x,1,footprint=fp,origin=(-1,0),axes=(2,0),mode='nearest'),
        ndi.median_filter(x,footprint=fp,origin=(-1,1),axes=(0,2),mode='reflect'),
        ndi.percentile_filter(x,60.,size=(2,4),axes=(2,0),origin=(-1,0),mode='nearest'),
    ]
    real_cases = [
        ndi.uniform_filter1d(x,size=6,axis=0,mode='reflect',origin=-1),
        ndi.uniform_filter1d(x,size=4,axis=2,mode='grid-wrap',origin=1),
        ndi.uniform_filter1d(x,size=50000,axis=1,mode='constant',cval=3.25),
        ndi.uniform_filter1d(x,size=50000,axis=1,mode='nearest',origin=-25000),
        ndi.uniform_filter1d(x,size=3,axis=1,mode='mirror'),
        ndi.uniform_filter(x,size=(4,3),axes=(2,0),mode=('wrap','constant'),cval=-1.25),
        ndi.uniform_filter(x,size=3,mode='wrap'),
        x,
        x,
    ]
    complex_cases = [
        ndi.uniform_filter1d(z,size=6,axis=0,mode='constant',cval=1.25-.7j,origin=-1),
        ndi.uniform_filter1d(z,size=5,axis=2,mode='reflect',origin=1),
        ndi.uniform_filter(z,size=(3,4),axes=(1,0),mode='mirror'),
    ]
    modes = ('reflect','mirror','nearest','wrap','grid-wrap',
             'constant','grid-constant','grid-mirror')
    random_order = []
    random_real = []
    random_complex = []
    for mode in modes:
        random_order.extend([
            ndi.median_filter(x,size=(2,4),mode=mode,axes=(2,0)),
            ndi.rank_filter(x,2,size=(2,4),mode=mode,axes=(2,0)),
            ndi.percentile_filter(x,55.,size=(3,2),mode=mode,axes=(0,2)),
        ])
        for width in range(2,10):
            random_real.append(ndi.uniform_filter1d(x,width,axis=1,mode=mode))
            random_complex.append(ndi.uniform_filter1d(z,width,axis=2,mode=mode))
    worst = 0.
    count = 0
    with TemporaryDirectory() as work:
        work = Path(work)
        cmd = ['gfortran','-std=f2018','-O1','-fcheck=all','-Werror=implicit-interface',
               '-J',str(work),'-I',str(work)]
        cmd += [str(ROOT/'src'/f'{stem}.f90') for stem in FILES]
        cmd += [str(ROOT/'tools/ndimage_reference/order_axes_running_driver.f90'), '-o', str(work/'run')]
        subprocess.run(cmd, cwd=work, check=True)
        subprocess.run([str(work/'run')], cwd=work, check=True)
        for filename,dtype,cases in [
            ('order_axes_output.bin',np.float64,order_cases),
            ('running_uniform_output.bin',np.float64,real_cases),
            ('running_uniform_complex_output.bin',np.complex128,complex_cases),
            ('order_axes_random_output.bin',np.float64,random_order),
            ('running_uniform_random_output.bin',np.float64,random_real),
            ('running_uniform_random_complex.bin',np.complex128,random_complex)]:
            data = np.fromfile(work/filename,dtype=dtype)
            assert data.size == len(cases)*x.size, (filename,data.size,len(cases)*x.size)
            for j,want in enumerate(cases):
                got = data[j*x.size:(j+1)*x.size].reshape(x.shape,order='F')
                err = np.max(np.abs(want-got))
                if err > 8.e-11:
                    raise AssertionError(f'{filename}, case {j}: error={err},\n got={got},\n want={want}')
                worst = max(worst,float(err))
                count += got.size
    print(f'PASS: {count} selected-axis order and running-mean SciPy values (max error {worst:.3g})')

if __name__ == '__main__':
    main()
