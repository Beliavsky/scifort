#!/usr/bin/env python3
"""Compare arbitrary footprints, gray morphology and EDT feature indices to SciPy."""
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import ndimage as ndi

ROOT = Path(__file__).resolve().parents[2]
MODES = ('reflect', 'mirror', 'nearest', 'wrap', 'constant',
         'grid-constant', 'grid-mirror', 'grid-wrap')


def main():
    i, j = np.mgrid[1:7, 1:8]
    a = np.sin(.38*i+.67*j)+(7*i+11*j+i*j)%13
    mask = (i*i+3*j+i*j)%7<3
    allowed = (i+2*j)%5!=0
    seeds = (2*i+j)%11==0
    fp3 = np.zeros((3,3), bool, order='F')
    fp3[0,1] = fp3[1,1] = fp3[2,0] = fp3[1,2] = True
    fp23 = np.zeros((2,3), bool, order='F')
    fp23[0,0] = fp23[1,1] = fp23[1,2] = True
    sym3 = np.eye(3,dtype=bool)
    cases = []
    for fp,border in [(fp3,0),(fp23,1)]:
        for iterations in (1,2):
            opts = dict(structure=fp,iterations=iterations,border_value=border)
            for fn in [ndi.binary_dilation,ndi.binary_erosion,ndi.binary_opening,ndi.binary_closing]:
                cases.append((fn.__name__,fn(mask,**opts)))
    for fp in (fp3,fp23):
        cases.append(('binary_propagation',ndi.binary_propagation(seeds,structure=fp,mask=allowed)))
    cases.append(('binary_fill_holes',ndi.binary_fill_holes(mask,structure=fp3)))
    expected_labels,nlabels = ndi.label(mask,structure=sym3)

    for kw in ({},{'structure1':fp3},{'structure1':fp23},
               {'structure1':fp3,'structure2':sym3},{'structure2':sym3}):
        cases.append(('binary_hit_or_miss',ndi.binary_hit_or_miss(mask,**kw)))

    for mode in MODES:
        kw=dict(mode=mode,cval=-2.7)
        for fn in (ndi.grey_dilation,ndi.grey_erosion,ndi.grey_opening,ndi.grey_closing,
                   ndi.morphological_gradient,ndi.morphological_laplace,
                   ndi.white_tophat,ndi.black_tophat):
            cases.append((fn.__name__+'_odd_'+mode,fn(a,footprint=fp3,**kw)))
        for fn in (ndi.grey_dilation,ndi.grey_erosion,ndi.grey_opening,ndi.grey_closing):
            cases.append((fn.__name__+'_even_'+mode,fn(a,footprint=fp23,**kw)))
        for fn in (ndi.grey_dilation,ndi.grey_erosion):
            cases.append((fn.__name__+'_rect_'+mode,fn(a,size=(2,3),**kw)))
    for x,sampling in [(mask,(1.7,.6)),(np.ones((6,7),bool),None)]:
        ds,inds=ndi.distance_transform_edt(x,sampling=sampling,return_indices=True)
        cases.append(('edt_dist',ds))
        cases.append(('edt_indices',inds+1))
    i3,j3,k3=np.mgrid[1:5,1:4,1:4]
    a3=np.cos(.24*i3+.49*j3-.31*k3)+(11*i3+7*j3+3*k3)%9
    m3=(i3*i3+3*j3+5*k3)%6<3
    fp3d=np.zeros((2,1,3),bool,order='F')
    fp3d[0,0,0]=fp3d[0,0,2]=fp3d[1,0,1]=True
    threed=[('binary_dilation3',ndi.binary_dilation(m3,structure=fp3d)),
            ('binary_erosion3',ndi.binary_erosion(m3,structure=fp3d)),
            ('binary_hit_or_miss3',ndi.binary_hit_or_miss(m3,structure1=fp3d))]
    for mode in MODES:
        for fn in (ndi.grey_dilation,ndi.grey_erosion,ndi.morphological_gradient):
            threed.append((fn.__name__+'3_'+mode,
                           fn(a3,footprint=fp3d,mode=mode,cval=-1.5)))

    with tempfile.TemporaryDirectory(prefix='scifort-footprint-feature-') as tmp:
        tmp=Path(tmp)
        files=['src/kinds.f90','src/signal_spline.f90','src/signal_spline_nd.f90',
               'src/ndimage_morphology.f90','src/ndimage_distance.f90','src/ndimage_chamfer.f90','src/ndimage_measurements.f90','src/ndimage_watershed.f90',
               'src/ndimage.f90',
               'tools/ndimage_reference/footprint_feature_driver.f90']
        cmd=['gfortran','-std=f2018','-O2','-g','-fcheck=all','-Werror=implicit-interface',
             '-J',str(tmp),'-I',str(tmp),*(str(ROOT/f) for f in files),'-o',str(tmp/'driver')]
        subprocess.run(cmd, cwd=tmp,check=True)
        subprocess.run([str(tmp/'driver')],cwd=tmp,check=True)
        with (tmp/'footprint_feature_results.bin').open('rb') as f:
            count=0
            maxerr=0.0
            for name,expected in cases[:19]:
                actual=np.frombuffer(f.read(42*4),dtype=np.int32).reshape((6,7),order='F')
                if not np.array_equal(actual!=0,expected):
                    ix=np.unravel_index(np.flatnonzero(actual.ravel(order='F')!=expected.ravel(order='F'))[0],(6,7),order='F')
                    raise AssertionError((name,ix,actual[ix],expected[ix]))
                count+=1
            labels=np.frombuffer(f.read(42*4),dtype=np.int32).reshape((6,7),order='F')
            n=int(np.frombuffer(f.read(4),dtype=np.int32)[0])
            assert n==nlabels,(n,nlabels)
            assert np.array_equal(labels>0,mask)
            for ix in zip(*np.nonzero(mask)):
                if not np.array_equal(labels==labels[ix], expected_labels==expected_labels[ix]):
                    print("label mismatch",ix,"Fortran",labels,"SciPy",expected_labels,sep="\n")
                    raise AssertionError("label partition differs")
            count+=1
            for name,expected in cases[19:24]:
                actual=np.frombuffer(f.read(42*4),dtype=np.int32).reshape((6,7),order='F')
                assert np.array_equal(actual!=0,expected),(name,actual!=0,expected)
                count+=1
            for name,expected in cases[24:]:
                nelems=expected.size
                if name=='edt_indices':
                    actual=np.frombuffer(f.read(4*nelems),dtype=np.int32).reshape(expected.shape,order='F')
                else:
                    actual=np.frombuffer(f.read(8*nelems),dtype=np.float64).reshape(expected.shape,order='F')
                err=np.max(np.abs(actual-expected))
                maxerr=max(maxerr,float(err))
                if err>3e-11:
                    idx=np.unravel_index(np.argmax(np.abs(actual-expected)),expected.shape)
                    raise AssertionError((name,float(err),idx,actual[idx],expected[idx]))
                count+=1
            for name,expected in threed:
                dtype=np.int32 if name.startswith('binary') else np.float64
                got=np.frombuffer(f.read(36*np.dtype(dtype).itemsize),dtype=dtype).reshape(
                    expected.shape,order='F')
                if dtype==np.int32:
                    assert np.array_equal(got!=0,expected),name
                else:
                    err=np.max(abs(got-expected))
                    maxerr=max(maxerr,float(err))
                    assert err<3e-11,(name,err)
                count+=1
            assert f.read(1)==b'', 'unexpected trailing binary data' 
    print(f'PASS: {count} footprint, grayscale and EDT-feature SciPy comparisons; maximum error={maxerr:.3g}')


if __name__=='__main__':
    main()
