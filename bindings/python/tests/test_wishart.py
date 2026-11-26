# SPDX-License-Identifier: MIT
"""Wishart and inverse-Wishart shape, buffer and numerical interface contracts."""
import ctypes as C
import math
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

import numpy as np
from scipy.stats import wishart, invwishart

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from build_library import build
from scifort import Scifort
from test_multivariate_t import run

ROOT = Path(__file__).resolve().parents[3]


class WishartTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran'):
            raise unittest.SkipTest('gfortran required')
        cls.builds = []
        for checked in (True, False):
            path, metadata = build(checked=checked)
            cls.builds.append((path, metadata, Scifort(
                path, dll_directories=metadata['dll_directories'])))

    def test_all_twelve_reference_methods(self):
        cases = (([[2.]], [[1.5]], 6.),
                 ([[2.,.4],[.4,1.]], [[1.5,.2],[.2,1.1]], 7.5),
                 ([[2.,.1,.3],[.1,1.,.2],[.3,.2,1.5]],
                  [[1.,.2,.1],[.2,2.,.3],[.1,.3,1.]], 9.))
        for _, _, api in self.builds:
            for name, ref in (('wishart',wishart),('invwishart',invwishart)):
                group=getattr(api,name)
                for scale,x,df in cases:
                    for method in ('pdf','logpdf'):
                        np.testing.assert_allclose(getattr(group,method)(x,df,scale),
                            getattr(ref,method)(x,df,scale),rtol=8e-13,atol=4e-14)
                    for method in ('entropy','mean','mode','var'):
                        np.testing.assert_allclose(getattr(group,method)(df,scale),
                            getattr(ref,method)(df,scale),rtol=8e-13,atol=4e-14)
                    self.assertAlmostEqual(math.log(group.pdf(x,df,scale)),group.logpdf(x,df,scale),places=12)
                    for method in ('mean','mode','var'):
                        out=getattr(group,method)(df,scale)
                        self.assertEqual(out.shape,(len(scale),len(scale)))
                        self.assertTrue(out.flags.f_contiguous)
                        np.testing.assert_allclose(out,out.T,rtol=0,atol=1e-14)
                    self.assertTrue(np.all(group.var(df,scale)>=0))

    def test_df_thresholds_undefined_moments_and_tail(self):
        scale=np.array([[2.,.4],[.4,1.]])
        for _, _, api in self.builds:
            for name,ref in (('wishart',wishart),('invwishart',invwishart)):
                group=getattr(api,name)
                for df in (1.01,1.5,2.,3.,3.01,5.,5.01,20.):
                    np.testing.assert_allclose(group.logpdf(np.eye(2),df,scale),
                        ref.logpdf(np.eye(2),df,scale),rtol=8e-13,atol=5e-14)
                    np.testing.assert_allclose(group.entropy(df,scale),ref.entropy(df,scale),rtol=8e-13,atol=5e-14)
                for df in (1.,0.,-1.,np.nan,np.inf,-np.inf):
                    for method in ('pdf','logpdf'):
                        self.assertTrue(math.isnan(getattr(group,method)(np.eye(2),df,scale)))
                    self.assertTrue(math.isnan(group.entropy(df,scale)))
                    for method in ('mean','mode','var'):
                        self.assertTrue(np.all(np.isnan(getattr(group,method)(df,scale))))
            self.assertTrue(np.all(np.isnan(api.wishart.mode(2.99,scale))))
            np.testing.assert_array_equal(api.wishart.mode(3.,scale),np.zeros((2,2)))
            self.assertTrue(np.all(np.isnan(api.invwishart.mean(3.,scale))))
            self.assertTrue(np.all(np.isnan(api.invwishart.var(5.,scale))))
            self.assertTrue(np.all(np.isfinite(api.invwishart.mean(3.01,scale))))
            self.assertTrue(np.all(np.isfinite(api.invwishart.var(5.01,scale))))
            # Underflow in PDF must not destroy an otherwise finite log density.
            for name,x in (('wishart',np.eye(2)*1e5),('invwishart',np.eye(2)*1e-5)):
                group=getattr(api,name)
                self.assertEqual(group.pdf(x,7.,np.eye(2)),0.)
                actual=group.logpdf(x,7.,np.eye(2))
                self.assertTrue(math.isfinite(actual))
                ref=wishart if name=='wishart' else invwishart
                np.testing.assert_allclose(actual,ref.logpdf(x,7.,np.eye(2)),rtol=4e-14,atol=1e-10)

    def test_invalid_matrices_and_triangle_conventions(self):
        scale=np.array([[2.,.4],[.4,1.]])
        x=np.array([[1.5,.2],[.2,1.1]])
        for _,_,api in self.builds:
            for name in ('wishart','invwishart'):
                group=getattr(api,name)
                for bad in (np.zeros((2,2)),[[1.,2.],[2.,1.]],
                            [[np.nan,0.],[0.,1.]],[[1.,np.inf],[0.,1.]]):
                    self.assertTrue(math.isnan(group.pdf(bad,7.,scale)))
                    self.assertTrue(math.isnan(group.logpdf(x,7.,bad)))
                    self.assertTrue(math.isnan(group.entropy(7.,bad)))
                    for method in ('mean','mode','var'):
                        self.assertTrue(np.all(np.isnan(getattr(group,method)(7.,bad))))
                upper=scale.copy();upper[0,1]=99.
                for method in ('pdf','logpdf'):
                    self.assertEqual(getattr(group,method)(x,7.,upper),getattr(group,method)(x,7.,scale))
                self.assertEqual(group.entropy(7.,upper),group.entropy(7.,scale))
                # The moments deliberately retain the full supplied scale entries.
                self.assertNotEqual(group.mean(7.,upper)[0,1],group.mean(7.,scale)[0,1])
                self.assertNotEqual(group.var(7.,upper)[0,1],group.var(7.,scale)[0,1])
            upper=x.copy();upper[0,1]=9.
            self.assertEqual(api.invwishart.logpdf(upper,7.,scale),api.invwishart.logpdf(x,7.,scale))
            self.assertNotEqual(api.wishart.logpdf(upper,7.,scale),api.wishart.logpdf(x,7.,scale))

    def test_python_shapes_layout_buffers_and_aliases(self):
        scale=np.array([[2.,.4],[.4,1.]])
        for _,_,api in self.builds:
            for name in ('wishart','invwishart'):
                group=getattr(api,name)
                for bad in ([],[1,2],np.empty((0,0)),np.ones((2,3)),np.ones((2,2,1))):
                    with self.assertRaises(ValueError):group.mean(7.,bad)
                with self.assertRaises(ValueError):group.pdf(np.ones((1,1)),7.,scale)
                with self.assertRaises(TypeError):group.mean(7.,scale.astype(complex))
                with self.assertRaises((TypeError,ValueError)):group.mean([7.,8.],scale)
                variants=(scale.tolist(),np.array(scale,order='C'),np.array(scale,dtype='>f8'))
                for value in variants:
                    np.testing.assert_allclose(group.mean(7.,value),group.mean(7.,scale),rtol=0,atol=0)
                # An actual noncontiguous view with the same matrix values.
                storage=np.zeros((4,4));storage[::2,::2]=scale
                np.testing.assert_array_equal(group.mean(7.,storage[::2,::2]),group.mean(7.,scale))
                for method in ('mean','mode','var'):
                    out=np.full((2,2),99.,order='F')
                    self.assertIs(getattr(group,method)(7.,scale,out=out),out)
                    np.testing.assert_array_equal(out,getattr(group,method)(7.,scale))
                    for bad in (np.empty((2,2),order='C'),np.empty((2,2),dtype='float32'),np.empty(4),{}):
                        with self.assertRaises((ValueError,TypeError)):getattr(group,method)(7.,scale,out=bad)
                    out.fill(99.)
                    with self.assertRaises(ValueError):getattr(group,method)(7.,np.ones((2,3)),out=out)
                    np.testing.assert_array_equal(out,99.)
                    out.flags.writeable=False
                    with self.assertRaises(ValueError):getattr(group,method)(7.,scale,out=out)
                a=np.array(scale,order='F')
                with self.assertRaises(ValueError):group.mean(7.,a,out=a)
                with self.assertRaises(ValueError):group.pdf(a,7.,a)

    def test_raw_c_status_and_output_preservation(self):
        for _,_,api in self.builds:
            x=np.eye(2,order='F');scale=np.array([[2.,.4],[.4,1.]],order='F')
            xp=x.ctypes.data_as(C.c_void_p);sp=scale.ctypes.data_as(C.c_void_p)
            for name in ('wishart','invwishart'):
                for method,fn in getattr(api,name)._functions.items():
                    out=np.full(4,99.);op=out.ctypes.data_as(C.POINTER(C.c_double));status=C.c_int(-1)
                    def call(p,df):
                        args=[p,xp,df,sp] if method in ('pdf','logpdf') else [p,df,sp]
                        fn(*args,op,C.byref(status))
                    for p in (0,-1):
                        call(p,7.)
                        self.assertEqual(status.value,1);np.testing.assert_array_equal(out,99.)
                    call(2,1.)
                    self.assertEqual(status.value,0)
                    count=4 if method in ('mean','mode','var') else 1
                    self.assertTrue(np.all(np.isnan(out[:count])))

    def test_native_and_c_clients(self):
        for path,metadata,_ in self.builds:
            with tempfile.TemporaryDirectory() as temp:
                folder=Path(temp)
                env=os.environ.copy();env['PATH']=str(path.parent)+os.pathsep+env.get('PATH','')
                flags=['-O0','-g','-fcheck=all','-fbacktrace'] if '-fcheck=all' in metadata['flags'] else ['-O3']
                exe=folder/'native.exe'
                run([metadata['compiler_path'],*flags,'-I',path.parent,
                     ROOT/'test/test_wishart_c_api.f90',path,'-o',exe]);run([exe],env=env)
                if shutil.which('gcc'):
                    c=folder/'client.c'
                    c.write_text("""#include "scifort.h"
#include <math.h>
int main(void) {
 double s[4]={2,.4,.4,1},x[4]={1.5,.2,.2,1.1},v[4],a,b; int status;
 scifort_wishart_pdf_f64(2,x,7,s,&a,&status); if(status || !(a>0)) return 1;
 scifort_wishart_logpdf_f64(2,x,7,s,&b,&status); if(status || fabs(log(a)-b)>1e-13) return 2;
 scifort_wishart_entropy_f64(2,7,s,&a,&status); if(status || !isfinite(a)) return 3;
 scifort_wishart_mean_f64(2,7,s,v,&status); if(status || fabs(v[0]-14)>1e-13) return 4;
 scifort_wishart_mode_f64(2,7,s,v,&status); if(status || fabs(v[0]-8)>1e-13) return 5;
 scifort_wishart_var_f64(2,7,s,v,&status); if(status || fabs(v[0]-56)>1e-13) return 6;
 scifort_invwishart_pdf_f64(2,x,7,s,&a,&status); if(status || !(a>0)) return 7;
 scifort_invwishart_logpdf_f64(2,x,7,s,&b,&status); if(status || fabs(log(a)-b)>1e-13) return 8;
 scifort_invwishart_entropy_f64(2,7,s,&a,&status); if(status || !isfinite(a)) return 9;
 scifort_invwishart_mean_f64(2,7,s,v,&status); if(status || fabs(v[0]-.5)>1e-13) return 10;
 scifort_invwishart_mode_f64(2,7,s,v,&status); if(status || fabs(v[0]-.2)>1e-13) return 11;
 scifort_invwishart_var_f64(2,7,s,v,&status); return status || fabs(v[0]-.25)>1e-13;
}
""",encoding='ascii')
                    exe=folder/'c.exe'
                    run(['gcc','-std=c11','-Wall','-Wextra','-Werror','-I',ROOT/'include',c,path,'-o',exe]);run([exe],env=env)
                    if shutil.which('g++'):
                        run(['g++','-std=c++17','-Wall','-Wextra','-Werror','-I',ROOT/'include','-fsyntax-only',c])


if __name__ == '__main__':
    unittest.main()
