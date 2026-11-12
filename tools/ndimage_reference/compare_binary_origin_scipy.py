#!/usr/bin/env python3
"""Check shifted binary dilation, erosion, opening, closing and hit-or-miss."""
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import ndimage as ndi

ROOT=Path(__file__).resolve().parents[2]

def main():
    i,j=np.mgrid[1:7,1:8]
    x=(i*i+2*j*i+3*j)%7<3
    a=np.array([1,0,0,1,1,0,0,1,1],bool).reshape((3,3),order='F')
    b=np.array([1,0,0,1,1,0],bool).reshape((2,3),order='F')
    methods=(ndi.binary_dilation,ndi.binary_erosion,ndi.binary_opening,
             ndi.binary_closing,ndi.binary_hit_or_miss)
    sources=[ROOT/z for z in ('src/kinds.f90','src/signal_spline.f90',
            'src/signal_spline_nd.f90','src/ndimage_morphology.f90',
            'src/ndimage_distance.f90','src/ndimage_chamfer.f90','src/ndimage_measurements.f90','src/ndimage_watershed.f90',
               'src/ndimage.f90',
            'tools/ndimage_reference/binary_origin_driver.f90')]
    with tempfile.TemporaryDirectory() as tmp:
        tmp=Path(tmp)
        subprocess.run(['gfortran','-std=f2018','-O0','-fcheck=all','-Werror=implicit-interface',
                        '-J',str(tmp),'-I',str(tmp),*map(str,sources),
                        '-o',str(tmp/'driver')],cwd=tmp,check=True)
        subprocess.run([str(tmp/'driver')],cwd=tmp,check=True)
        values=np.fromfile(tmp/'binary_origin_results.bin',np.int32)
        count=0
        for orig,fp in [((0,0),a),((-1,1),a),((1,-1),a),((-1,0),b)]:
            for fn in methods:
                kw={'structure1':fp,'origin1':orig} if fn==ndi.binary_hit_or_miss else dict(structure=fp,origin=orig)
                expected=fn(x,**kw)
                actual=values[count*42:(count+1)*42].reshape((6,7),order='F')!=0
                if not np.array_equal(expected,actual):
                    raise AssertionError((fn.__name__,orig,fp.shape,actual.astype(int),expected.astype(int)))
                count+=1
        assert count*42==values.size
        print(f'PASS: binary origins {count} SciPy comparisons')

if __name__=='__main__':
    main()
