#!/usr/bin/env python3
"""Cross-check exact feature indices with randomized anisotropic 2D/3D arrays."""
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import ndimage as ndi

ROOT=Path(__file__).resolve().parents[2]


def main():
    rng=np.random.default_rng(87124)
    masks=[]
    for k in range(64):
        a=rng.uniform(size=(5,6))<rng.uniform(.08,.96)
        b=rng.uniform(size=(4,3,3))<rng.uniform(.08,.96)
        if k%8==0:
            a[:]=True
            a[0,0]=False
        if k%9==0:
            b[:]=True
            b[0,0,0]=False
        if k%11==0:
            a[:]=True
        if k%13==0:
            b[:]=False
        masks.extend([(np.asfortranarray(a),(1.75,.6)),(np.asfortranarray(b),(.85,1.45,2.5))])
    with tempfile.TemporaryDirectory(prefix='scifort-random-features-') as tmp:
        tmp=Path(tmp)
        with (tmp/'feature_random_inputs.bin').open('wb') as f:
            for a,_ in masks:
                f.write(np.asfortranarray(a.astype(np.int32)).tobytes(order='F'))
        sources=['src/kinds.f90','src/ndimage_distance.f90',
                 'tools/ndimage_reference/random_feature_driver.f90']
        flags=['gfortran','-std=f2018','-O2','-g','-fcheck=all',
               '-ffpe-trap=invalid,zero,overflow','-Werror=implicit-interface',
               '-J',str(tmp),'-I',str(tmp)]
        subprocess.run([*flags,*(str(ROOT/f) for f in sources),'-o',str(tmp/'driver')],
                       cwd=tmp,check=True)
        subprocess.run([str(tmp/'driver')],cwd=tmp,check=True)
        exact=0
        count=0
        maxerr=0.0
        with (tmp/'feature_random_results.bin').open('rb') as f:
            for a,spacing in masks:
                ref,feat=ndi.distance_transform_edt(a,sampling=spacing,return_indices=True)
                got=np.frombuffer(f.read(8*a.size),dtype=np.float64).reshape(a.shape,order='F')
                ix=np.frombuffer(f.read(4*a.ndim*a.size),dtype=np.int32).reshape(
                    (a.ndim,a.size),order='F').reshape((a.ndim,*a.shape),order='F')
                err=np.max(abs(ref-got))
                maxerr=max(maxerr,float(err))
                assert err<5e-13,(count,err)
                if np.array_equal(ix,feat+1):
                    exact+=1
                site=ix-1
                if a.all():
                    assert np.all(site[0]==-1)
                    assert all(np.all(site[d]==0) for d in range(1,a.ndim))
                else:
                    assert all(np.all((site[d]>=0)&(site[d]<a.shape[d])) for d in range(a.ndim))
                    assert np.all(~a[tuple(site)])
                xyz=np.indices(a.shape)
                implied=np.sqrt(sum(((xyz[d]-site[d])*spacing[d])**2 for d in range(a.ndim)))
                assert np.allclose(implied,got,rtol=0,atol=1e-12),(count,'nonminimal site')
                count+=1
            assert f.read(1)==b''
    print(f'PASS: {count} randomized EDT feature maps, {exact} exact SciPy index maps, maximum distance error={maxerr:.3g}')

if __name__=='__main__':
    main()
