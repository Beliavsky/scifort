#!/usr/bin/env python3
"""Check SciFort mixed-type batch callbacks against SciPy and scalar references."""
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import ndimage as ndi

ROOT=Path(__file__).resolve().parents[2]
MODES=('reflect','nearest','mirror','wrap','constant')


def rc(w, *, axis):
    return np.mean(w,axis=axis)+1j*np.mean(w*w,axis=axis)


def cr(w, *, axis):
    return np.mean(np.abs(w)**2,axis=axis)


def rc_second(w, *, axis):
    return np.mean(w*w,axis=axis)+1j*np.mean(w,axis=axis)


def cr_second(w, *, axis):
    return np.mean(w.real,axis=axis)


def main():
    ii,jj,kk=np.indices((5,4,3))
    ii,jj,kk=ii+1,jj+1,kk+1
    x=np.sin(.31*ii+.23*jj)+.17*kk+.09*ii*jj
    z=x+1j*(.37*x+.21)
    cases=[]
    for mode in MODES:
        rs={'size':(2,3,2),'mode':mode}
        cs={'size':(2,3,2),'mode':mode}
        if mode=='constant':
            rs['cval']=-1.75
            cs['cval']=-1.75+.6j
        a=ndi.vectorized_filter(x,rc,**rs)
        b=ndi.vectorized_filter(z,cr,**cs)
        cases.extend([(f'{mode} real->complex',a),(f'{mode} complex->real',b),
                      (f'{mode} real->complex multi',np.stack((a,ndi.vectorized_filter(x,rc_second,**rs)),axis=-1)),
                      (f'{mode} complex->real multi',np.stack((b,ndi.vectorized_filter(z,cr_second,**cs)),axis=-1))])
    rs={'size':(2,3,2),'mode':'valid'}
    cases.extend([('real->complex valid',ndi.vectorized_filter(x,rc,**rs)),
                  ('complex->real valid',ndi.vectorized_filter(z,cr,**rs))])
    cases.append(('real->complex multi valid',np.stack((
        ndi.vectorized_filter(x,rc,**rs),ndi.vectorized_filter(x,rc_second,**rs)),axis=-1)))
    cases.append(('complex->real multi valid',np.stack((
        ndi.vectorized_filter(z,cr,**rs),ndi.vectorized_filter(z,cr_second,**rs)),axis=-1)))
    rs={'size':(2,3),'axes':(2,0),'mode':'nearest'}
    cases.extend([('real->complex selected axes',ndi.vectorized_filter(x,rc,**rs)),
                  ('complex->real selected axes',ndi.vectorized_filter(z,cr,**rs))])
    with tempfile.TemporaryDirectory(prefix='scifort-mixed-batched-') as tmpstr:
        tmp=Path(tmpstr)
        sources=['src/kinds.f90','src/signal_spline.f90','src/signal_spline_nd.f90',
                 'src/ndimage_morphology.f90','src/ndimage_distance.f90',
                 'src/ndimage_chamfer.f90','src/ndimage_measurements.f90',
                 'src/ndimage_watershed.f90','src/ndimage.f90',
                 'tools/ndimage_reference/vectorized_mixed_driver.f90']
        cmd=['gfortran','-std=f2018','-O0','-fcheck=all','-Werror=implicit-interface',
             '-J',str(tmp),'-I',str(tmp)]
        subprocess.run([*cmd,*(str(ROOT/p) for p in sources),'-o',str(tmp/'driver')],
                       cwd=tmp,check=True)
        subprocess.run([str(tmp/'driver')],cwd=tmp,check=True,capture_output=True)
        maximum=0.0
        total=0
        with (tmp/'vectorized_mixed_reference.bin').open('rb') as f:
            for name,expected in cases:
                dtype=np.complex128 if np.iscomplexobj(expected) else np.float64
                data=f.read(expected.size*np.dtype(dtype).itemsize)
                if len(data)!=expected.size*np.dtype(dtype).itemsize:
                    raise AssertionError(f'{name}: truncated results')
                actual=np.frombuffer(data,dtype=dtype).reshape(expected.shape,order='F')
                err=float(np.max(np.abs(actual-expected)))
                if not np.isfinite(err) or err>1e-11:
                    raise AssertionError(f'{name}: max error {err} values {actual.flat[0]} vs {expected.flat[0]}')
                total+=expected.size
                maximum=max(maximum,err)
            if f.read(1):
                raise AssertionError('Unexpected binary output')
    print(f'PASS: {total} mixed real/complex vectorized values across {len(cases)} cases')
    print(f'Maximum absolute error: {maximum:.3g}')


if __name__=='__main__':
    main()
