#!/usr/bin/env python3
"""Direct SciPy checks for compact-axis grayscale morphology and callbacks."""
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
    a,b,c=np.indices((4,5,3))+1
    x=np.sin(.17*a+.31*b-.19*c)+.05*(a*a-2*b+c)
    z=x+1j*(.35*x+.11)
    fp=np.array([[True,True,False], [False,True,True]],order='F')
    s=np.array([[-.3,.2,.4],[.8,-.1,-.2]],order='F')
    modes=('reflect','nearest','mirror','wrap','constant','grid-constant',
           'grid-wrap','grid-mirror')
    # SciPy's nonseparable filter API rejects a sequence of boundary modes.
    # For the separable flat stencil, compose one-dimensional SciPy filters
    # to provide an independent mixed-mode reference.
    gray_mixed=ndi.maximum_filter1d(x,size=2,axis=2,origin=-1,
                                    mode='constant',cval=-1.2)
    gray_mixed=ndi.maximum_filter1d(gray_mixed,size=3,axis=0,mode='nearest')

    def mixed_generic(image,fill,mode_axis2):
        """Coordinate-wise reference, with canonical compact footprint order."""
        output=np.empty_like(image)
        for i,j,k in np.ndindex(image.shape):
            values=[]
            for di in (-1,0):
                for dk in (-1,0,1):
                    pi=i+di
                    pk=k+dk
                    if not 0<=pi<image.shape[0]:
                        values.append(fill)
                    else:
                        # The second selected axis is periodic in this test.
                        assert mode_axis2 in ('wrap','grid-wrap')
                        values.append(image[pi,j,pk%image.shape[2]])
            output[i,j,k]=sum(v*(t+1) for t,v in enumerate(values))
        return output

    def mixed_nonflat(dilate):
        """Independent stencil traversal for mixed-mode asymmetric grayscale filtering."""
        out=np.empty_like(x)
        for i,j,k in np.ndindex(x.shape):
            candidates=[]
            for fi in range(fp.shape[0]):
                for fk in range(fp.shape[1]):
                    if not fp[fi,fk]:
                        continue
                    pi=i-fi if dilate else i+fi-1
                    pk=k+2-fk if dilate else k+fk-1
                    pix=-2.0 if not 0<=pi<x.shape[0] else x[pi,j,pk % x.shape[2]]
                    candidates.append(pix+s[fi,fk] if dilate else pix-s[fi,fk])
            out[i,j,k]=max(candidates) if dilate else min(candidates)
        return out

    gray=[
        ndi.grey_dilation(x,size=(2,3),axes=(2,0)),
        ndi.grey_erosion(x,size=(3,2),axes=(0,2),origin=(-1,0),mode='constant',cval=-3.),
        ndi.grey_opening(x,size=(2,3),axes=(2,0),origin=(0,1)),
        ndi.grey_closing(x,size=(2,3),axes=(2,0),mode='nearest'),
        ndi.morphological_gradient(x,footprint=fp,axes=(2,0),mode='mirror'),
        ndi.morphological_laplace(x,structure=s,axes=(2,0),mode='grid-wrap'),
        ndi.white_tophat(x,footprint=fp,structure=s,axes=(2,0),mode='nearest'),
        ndi.black_tophat(x,footprint=fp,structure=s,axes=(2,0),mode='nearest'),
        ndi.grey_dilation(x,size=(2,3),axes=(2,0),mode='constant',cval=-1.2),
        ndi.grey_erosion(x,size=(),axes=()),
        ndi.grey_dilation(x,size=(2,3,2),mode='reflect'),
        gray_mixed,
        mixed_nonflat(True),
        mixed_nonflat(True)-mixed_nonflat(False),
    ]
    for m in modes:
        gray.append(ndi.grey_erosion(x,structure=s,axes=(2,0),mode=m))
        gray.append(ndi.grey_dilation(x,structure=s,axes=(2,0),mode=m))
        gray.append(ndi.grey_erosion(x,structure=s,axes=(2,0),origin=(-1,0),mode=m))
        gray.append(ndi.grey_dilation(x,structure=s,axes=(2,0),origin=(-1,0),mode=m))

    def weighted(v):
        return float(v @ np.arange(1,len(v)+1))
    generic=[
        ndi.generic_filter(x,weighted,size=(2,3),axes=(2,0)),
        ndi.generic_filter(x,weighted,footprint=fp,axes=(2,0),origin=(1,-1),
                           mode='constant',cval=-1.),
        ndi.generic_filter(x,weighted,size=(3,2),axes=(0,2),origin=(1,-1),mode='nearest'),
        ndi.generic_filter(x,weighted,size=(),axes=()),
        ndi.generic_filter(x,weighted,size=(2,3,2),mode='mirror'),
        mixed_generic(x,-2.0,'wrap'),
    ]
    generic.extend(ndi.generic_filter(x,weighted,footprint=fp,axes=(2,0),mode=m)
                   for m in modes)
    # Real/imaginary channels can be filtered independently for a linear
    # complex callback; scipy.ndimage.generic_filter itself rejects complex.
    complex_cases=[
        ndi.generic_filter(z.real,weighted,size=(2,3),axes=(2,0))+
            1j*ndi.generic_filter(z.imag,weighted,size=(2,3),axes=(2,0)),
        ndi.generic_filter(z.real,weighted,footprint=fp,axes=(2,0),mode='constant',cval=2.)+
            1j*ndi.generic_filter(z.imag,weighted,footprint=fp,axes=(2,0),mode='constant',cval=-.7),
        z.copy(),
        mixed_generic(z,-2.0+.5j,'grid-wrap'),
    ]
    worst=0.
    ncomp=0
    with TemporaryDirectory() as temp:
        wd=Path(temp)
        cmd=['gfortran','-std=f2018','-O1','-fcheck=all','-Werror=implicit-interface',
             '-J',str(wd),'-I',str(wd)]
        cmd += [str(ROOT/'src'/f'{stem}.f90') for stem in FILES]
        cmd += [str(ROOT/'tools/ndimage_reference/gray_generic_axes_driver.f90'),'-o',str(wd/'run')]
        subprocess.run(cmd,cwd=wd,check=True)
        subprocess.run([str(wd/'run')],cwd=wd,check=True)
        for fname,dtype,cases in [('gray_axes_output.bin',np.float64,gray),
            ('generic_axes_output.bin',np.float64,generic),
            ('generic_axes_complex_output.bin',np.complex128,complex_cases)]:
            data=np.fromfile(wd/fname,dtype=dtype)
            assert data.size==x.size*len(cases),(fname,data.size,len(cases))
            for j,expected in enumerate(cases):
                if expected is None: continue
                got=data[j*x.size:(j+1)*x.size].reshape(x.shape,order='F')
                err=float(np.max(np.abs(got-expected)))
                if err>2.e-10:
                    raise AssertionError(f'{fname} case {j}: error {err}, got={got}, expected={expected}')
                worst=max(worst,err)
                ncomp+=got.size
    print(f'PASS: {ncomp} SciPy grayscale/generic selected-axis values (worst={worst:.3g})')

if __name__=='__main__':
    main()
