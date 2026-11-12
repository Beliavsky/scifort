#!/usr/bin/env python3
"""Compare SciFort assumed-rank interpolation, affine transforms and shifts to SciPy."""
from pathlib import Path
import tempfile
import subprocess
import numpy as np
from scipy import ndimage as ndi

ROOT = Path(__file__).resolve().parents[2]
MODES = ('mirror','constant','reflect','grid-mirror','nearest','wrap','grid-wrap','grid-constant')

def main():
    x = np.empty((4, 5), order='F')
    for j in range(1,6):
        for i in range(1,5):
            x[i-1,j-1] = np.sin(.31*i+.42*j)+np.cos(.23*i-.19*j)
    z = x + 1j*(.37*x+.21)
    coords=np.empty((2,41),order='F')
    for k in range(1,42):
        coords[0,k-1] = .2+.37*(k-11)
        coords[1,k-1] = 1.3+.49*((9*k)%41-12)
    coords[:,0] = [1,1]
    coords[:,1] = [4,5]
    coords[:,2] = [.35,5.4]
    coords[:,3] = [-2,6]
    coords[:,4] = [8,-4]
    coords[:,5] = [.5,5.5]
    coords[:,6] = [.5,1.]
    coords[:,7] = [4.5,1.]
    coords[:,8] = [1.,5.5]
    coords[:,9] = [4.,.5]
    coords[:,10] = [3.5,5.5]
    matrix=np.array([[.83,.12],[-.16,.93]])
    offset=[.3,-.5]
    n=0
    maximum_strict=0.0
    maximum_reflect=0.0
    mode_errors={}
    with tempfile.TemporaryDirectory(prefix='scifort-spline-modes-') as t:
        t=Path(t)
        subprocess.run(['gfortran','-std=f2018','-O2','-fcheck=all',
                        '-Werror=implicit-interface','-J',str(t),'-I',str(t),
                        str(ROOT/'src/kinds.f90'),str(ROOT/'src/signal_spline.f90'),
                        str(ROOT/'src/signal_spline_nd.f90'),
                        str(ROOT/'tools/signal_reference/spline_modes_reference_driver.f90'),
                        '-o',str(t/'driver')],check=True,cwd=t)
        subprocess.run([str(t/'driver')],check=True,cwd=t)
        with (t/'spline_modes_results.bin').open('rb') as f:
            for mode in MODES:
                for order in range(6):
                    for prefilter in (False, True):
                        statuses=np.frombuffer(f.read(24),dtype=np.int32)
                        if any(statuses!=0):
                            raise AssertionError((mode,order,prefilter,statuses))
                        values=[]
                        for count, dtype, shape in ((41,np.float64,(41,)),
                            (41,np.complex128,(41,)),(18,np.float64,(3,6)),
                            (18,np.complex128,(3,6)),(20,np.float64,(4,5)),
                            (20,np.complex128,(4,5))):
                            vec=np.frombuffer(f.read(count*np.dtype(dtype).itemsize),dtype=dtype)
                            values.append(vec.reshape(shape,order='F'))
                        refs=[]
                        for data,fill in ((x,-2.25),(z,-2.25+.5j)):
                            refs.append(ndi.map_coordinates(data,coords-1,order=order,
                                mode=mode,prefilter=prefilter,cval=fill))
                        for data,fill in ((x,-2.25),(z,-2.25+.5j)):
                            refs.append(ndi.affine_transform(data,matrix,offset=offset,
                                output_shape=(3,6),order=order,mode=mode,
                                cval=fill,prefilter=prefilter))
                        for data,fill in ((x,-2.25),(z,-2.25+.5j)):
                            refs.append(ndi.shift(data,shift=[.45,-1.3],order=order,
                                mode=mode,cval=fill,prefilter=prefilter))
                        for j,(got,want) in enumerate(zip(values,refs)):
                            err=float(np.max(abs(got-want)))
                            approximate=mode in ('reflect','grid-mirror') and prefilter and order>1
                            if approximate:
                                maximum_reflect=max(maximum_reflect,err)
                                tolerance=1.5e-3
                            else:
                                maximum_strict=max(maximum_strict,err)
                                tolerance=5.0e-10
                            mode_errors[mode]=max(mode_errors.get(mode,0.0),err)
                            if err>tolerance:
                                ii=np.unravel_index(np.argmax(abs(got-want)),got.shape)
                                raise AssertionError((mode,order,prefilter,j,err,ii,got[ii],want[ii]))
                            n+=1
            if f.read(1): raise AssertionError('extra output bytes')
    print(f'PASS: {n}/576 interpolation/affine/shift comparisons')
    print(f'Max strict-parity error: {maximum_strict:.4g}; '
          f'max short-array reflect-prefilter difference: {maximum_reflect:.4g}')

if __name__=='__main__':
    main()
