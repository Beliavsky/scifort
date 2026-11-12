#!/usr/bin/env python3
"""Randomized SciPy comparisons for city-block/chessboard distances and sites."""
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import ndimage as ndi

ROOT=Path(__file__).resolve().parents[2]
SOURCES=[ROOT/p for p in (
    'src/kinds.f90','src/signal_spline.f90','src/signal_spline_nd.f90',
    'src/ndimage_morphology.f90','src/ndimage_distance.f90','src/ndimage_chamfer.f90',
    'src/ndimage_measurements.f90','src/ndimage_watershed.f90','src/ndimage.f90',
    'tools/ndimage_reference/chamfer_random_driver.f90')]


def main():
    rng=np.random.default_rng(20261022)
    masks=[]
    for n in range(160):
        pair=[]
        for dims in ((5,6),(4,3,2)):
            if n==0:
                x=np.zeros(dims,bool)
            elif n==1:
                x=np.ones(dims,bool)
            elif n==2:
                x=np.ones(dims,bool)
                x.flat[0]=False
            elif n==3:
                x=np.ones(dims,bool)
                x.flat[-1]=False
            else:
                x=rng.random(dims)<[.1,.35,.5,.7,.9][n%5]
            pair.append(x)
        masks.append(pair)
    with tempfile.TemporaryDirectory(prefix='scifort-cdt-') as work:
        work=Path(work)
        with (work/'chamfer_random_inputs.bin').open('wb') as f:
            for pair in masks:
                for x in pair:
                    f.write(np.asarray(x,dtype='<i4').ravel(order='F').tobytes())
        subprocess.run(['gfortran','-std=f2018','-O2','-fcheck=all','-Werror=implicit-interface',
            '-J',str(work),'-I',str(work),*[str(p) for p in SOURCES],'-o',str(work/'driver')],
            check=True,cwd=work)
        subprocess.run([str(work/'driver')],check=True,cwd=work)
        data=np.fromfile(work/'chamfer_random_results.bin',dtype='<i4')
    cursor=0
    comparisons=0
    for pair in masks:
        for metric in ('taxicab','chessboard'):
            for x in pair:
                n=x.size
                nd=x.ndim
                got=data[cursor:cursor+n].reshape(x.shape,order='F')
                cursor+=n
                indices=data[cursor:cursor+nd*n].reshape((nd,n),order='F')
                cursor+=nd*n
                expected=ndi.distance_transform_cdt(x,metric=metric)
                if not np.array_equal(got,expected):
                    raise AssertionError((metric,x,got,expected))
                if np.any(~x):
                    for pos,coordinate in enumerate(np.ndindex(x.shape)):
                        site=indices[:,np.ravel_multi_index(coordinate,x.shape,order='F')]-1
                        if np.any(site<0) or np.any(site>=np.array(x.shape)):
                            raise AssertionError(('out-of-range feature',site))
                        if x[tuple(site)]:
                            raise AssertionError(('foreground feature',site,coordinate))
                        separation=np.abs(np.array(coordinate)-site)
                        measurement=int(np.sum(separation) if metric=='taxicab' else np.max(separation))
                        if measurement!=got[coordinate]:
                            raise AssertionError(('incorrect feature distance',metric,coordinate,site,measurement,got[coordinate]))
                comparisons+=got.size
    assert cursor==data.size,(cursor,data.size)
    print(f'PASS: {comparisons} random SciPy chamfer distance values and nearest-site witnesses')


if __name__=='__main__':
    main()
