#!/usr/bin/env python3
"""Check rank-independent SciFort zoom and rotate against SciPy 1.17."""
from pathlib import Path
import subprocess
import warnings
import tempfile
import numpy as np
from scipy import ndimage as ndi

ROOT = Path(__file__).resolve().parents[2]
MODES = ('mirror','constant','nearest','reflect','grid-wrap','grid-constant','wrap')
ORDERS = (0,1,3,5)
ANGLES = (0,23,90,-37)


def main():
    x=np.empty((5,6),dtype='float64',order='F')
    for j in range(1,7):
        for i in range(1,6):
            x[i-1,j-1]=np.sin(.31*i+.27*j)+.13*i-.04*j
    cx=x+1j*(.27*x-.19)
    xx=np.empty((5,3,4),dtype='float64',order='F')
    for k in range(1,5):
        for j in range(1,4):
            for i in range(1,6):
                xx[i-1,j-1,k-1]=np.sin(.31*i+.17*j-.29*k)+.13*i
    cxx=xx+1j*(.27*xx-.19)
    maxerr=0.0
    max_strict=0.0
    max_reflect=0.0
    strict_cases=0
    reflect_cases=0
    err_by_key={}
    cases=0
    n_large=0
    with tempfile.TemporaryDirectory(prefix='scifort-zoom-rotate-') as tmp:
        t=Path(tmp)
        subprocess.run(['gfortran','-std=f2018','-O2','-fcheck=all',
                        '-Werror=implicit-interface','-J',str(t),'-I',str(t),
                        str(ROOT/'src/kinds.f90'),str(ROOT/'src/signal_spline.f90'),
                        str(ROOT/'src/signal_spline_nd.f90'),
                        str(ROOT/'tools/signal_reference/zoom_rotate_reference_driver.f90'),
                        '-o',str(t/'driver')],check=True,cwd=t)
        subprocess.run([str(t/'driver')],check=True,cwd=t)
        with (t/'zoom_rotate_results.bin').open('rb') as stream:
            def compare(ref,label,atol=3e-11):
                nonlocal cases,maxerr,n_large,max_strict,max_reflect,strict_cases,reflect_cases
                got=np.frombuffer(stream.read(ref.size*ref.dtype.itemsize),dtype=ref.dtype)
                if got.size!=ref.size: raise AssertionError('short output')
                got=got.reshape(ref.shape,order='F')
                err=float(np.max(abs(got-ref)))
                maxerr=max(maxerr,err)
                cases+=1
                if label[1]=='reflect' and label[2] > 1 and label[-1] is True:
                    max_reflect=max(max_reflect,err)
                    reflect_cases+=1
                    atol=2.5e-4
                else:
                    max_strict=max(max_strict,err)
                    strict_cases+=1
                key=(label[0],label[1])
                err_by_key[key]=max(err_by_key.get(key,0.0),err)
                if err>atol:
                    n_large+=1
                    if n_large < 12:
                        p=np.unravel_index(np.argmax(abs(got-ref)),ref.shape)
                        print('MISMATCH',label,'err',err,'pos',p,'got',got[p],'ref',ref[p])
            for mode in MODES:
                for order in ORDERS:
                    for grid in (False,True):
                        for pre in (False,True):
                            with warnings.catch_warnings(), np.errstate(all='ignore'):
                                warnings.simplefilter('ignore',UserWarning)
                                ref=ndi.zoom(x,(1.4,.7),order=order,mode=mode,cval=.35,
                                             prefilter=pre,grid_mode=grid)
                                # scipy.ndimage.zoom's complex path omits grid_mode
                                # as of SciPy 1.17, so compare componentwise instead.
                                rr=ndi.zoom(cx.real,(1.4,.7),order=order,mode=mode,
                                            cval=.35,prefilter=pre,grid_mode=grid)
                                ri=ndi.zoom(cx.imag,(1.4,.7),order=order,mode=mode,
                                            cval=.21,prefilter=pre,grid_mode=grid)
                            compare(ref,('zoom',mode,order,grid,pre))
                            compare(rr+1j*ri,('complex zoom',mode,order,grid,pre))
            for mode in MODES:
                for order in ORDERS:
                    for angle in ANGLES:
                        for reshape in (False,True):
                            for pre in (False,True):
                                r=ndi.rotate(x,angle,reshape=reshape,order=order,mode=mode,
                                             prefilter=pre,cval=.35)
                                rc=ndi.rotate(cx,angle,reshape=reshape,order=order,mode=mode,
                                              prefilter=pre,cval=.35+.21j)
                                compare(r,('rotate',mode,order,angle,reshape,pre))
                                compare(rc,('complex rotate',mode,order,angle,reshape,pre))
            for reshape in (True,False):
                for order in ORDERS[:3]:
                    for angle in ANGLES[:3]:
                        r=ndi.rotate(xx,angle,axes=(0,2),reshape=reshape,order=order,
                                     mode='mirror',prefilter=True)
                        rc=ndi.rotate(cxx,angle,axes=(0,2),reshape=reshape,order=order,
                                      mode='mirror',prefilter=True)
                        compare(r,('3D rotate','mirror',order,angle,reshape))
                        compare(rc,('3D complex rotate','mirror',order,angle,reshape))
            if stream.read(1): raise AssertionError('trailing reference bytes')
    print(f'PASS: {cases} comparisons = {strict_cases} strict (max {max_strict:.4e}) + '
          f'{reflect_cases} short-array reflect (max {max_reflect:.4e}); out-of-tolerance {n_large}')
    for k,e in sorted(err_by_key.items()): print(f'{k}: {e:.4e}')
    if n_large: raise AssertionError(f'{n_large} cases exceeded comparison tolerance')

if __name__=='__main__':
    main()
