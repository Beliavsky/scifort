#!/usr/bin/env python3
"""Validate selected-axis connected-component labeling against SciPy slices."""
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import ndimage as ndi

ROOT=Path(__file__).resolve().parents[2]
SOURCES=[ROOT / p for p in (
    'src/kinds.f90','src/signal_spline.f90','src/signal_spline_nd.f90',
    'src/ndimage_morphology.f90','src/ndimage_distance.f90',
    'src/ndimage_chamfer.f90','src/ndimage_measurements.f90',
    'src/ndimage_watershed.f90','src/ndimage.f90',
    'tools/ndimage_reference/label_axes_driver.f90',
)]


def validate(actual,mask,axes,structure,number):
    """Compare connectivity partitions, not arbitrary numeric component IDs."""
    if not np.array_equal(actual!=0,mask):
        raise AssertionError('foreground mismatch')
    axes=sorted(axes)
    other=[a for a in range(mask.ndim) if a not in axes]
    expected_total=0
    for slicepos in np.ndindex(*(mask.shape[a] for a in other)):
        sl=[slice(None)]*mask.ndim
        for d,p in zip(other,slicepos):
            sl[d]=p
        region=np.asarray(mask[tuple(sl)])
        assigned=np.asarray(actual[tuple(sl)])
        expected,components=ndi.label(region,structure=structure)
        expected_total+=components
        # Each SciPy component must correspond to exactly one SciFort ID,
        # and no SciFort ID may join separate SciPy components.
        mapping={}
        for local_id in range(1,components+1):
            global_ids=np.unique(assigned[expected==local_id])
            if global_ids.size!=1 or global_ids[0]==0:
                raise AssertionError(('split component',slicepos,local_id,global_ids))
            mapping[local_id]=int(global_ids[0])
        if len(set(mapping.values()))!=len(mapping):
            raise AssertionError(('merged components',slicepos,mapping))
    if number!=expected_total:
        raise AssertionError(('component count',number,expected_total))
    if np.unique(actual[actual>0]).size!=number:
        raise AssertionError('reused labels across disconnected slices')
    return actual.size


def main():
    i,j,k=np.mgrid[1:7,1:5,1:8]
    mask3=(2*i+3*j+7*k)%9<4
    i,j,k,l=np.mgrid[1:3,1:5,1:4,1:6]
    mask4=(2*i+3*j+k+5*l)%11<5
    diagonal=np.array([[1,0,1],[0,1,0],[1,0,1]],dtype=bool)
    face=ndi.generate_binary_structure(2,1)
    elements=0
    with tempfile.TemporaryDirectory(prefix='scifort-labelaxes-') as tmp:
        root=Path(tmp)
        subprocess.run(['gfortran','-std=f2018','-O0','-fcheck=all',
                        '-Werror=implicit-interface','-J',str(root),'-I',str(root),
                        *map(str,SOURCES),'-o',str(root/'driver')],cwd=root,check=True)
        subprocess.run([str(root/'driver')],cwd=root,check=True)
        with (root/'selected_axis_labels.bin').open('rb') as stream:
            def read(shape):
                n=np.frombuffer(stream.read(4),np.int32)
                if n.size!=1:
                    raise AssertionError('missing component count')
                count=int(np.prod(shape))
                values=np.frombuffer(stream.read(4*count),np.int32)
                if values.size!=count:
                    raise AssertionError('missing component labels')
                return int(n[0]),values.reshape(shape,order='F')

            n,actual=read(mask3.shape)
            elements+=validate(actual,mask3,[2,0],diagonal,n)
            n,actual=read(mask4.shape)
            elements+=validate(actual,mask4,[3,1],face,n)
            n,actual=read(mask3.shape)
            if n!=np.count_nonzero(mask3):
                raise AssertionError('empty-axis components')
            if np.count_nonzero(actual)!=n or np.unique(actual[actual>0]).size!=n:
                raise AssertionError('empty-axis label uniqueness')
            elements+=actual.size
            if stream.read(1):
                raise AssertionError('unexpected trailing output')
    print(f'PASS: {elements} selected-axis component labels vs SciPy slice partitions')


if __name__=='__main__':
    main()
