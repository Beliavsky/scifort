# SPDX-License-Identifier: MIT
"""KDE handle lifecycle, native/SciPy agreement and array contracts."""
import copy
import ctypes as C
import gc
import os
from pathlib import Path
import pickle
import shutil
import sys
import tempfile
import unittest
import weakref
import numpy as np
from scipy.stats import gaussian_kde
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from build_library import build
from scifort import Scifort
from test_multivariate_t import run
ROOT=Path(__file__).resolve().parents[3]
DATA=np.array([[-2.,-.5,.25,1.5,3.,.8],[.5,2.,-1.,.2,2.8,-.4]])

class GaussianKDETests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran'):raise unittest.SkipTest('gfortran required')
        cls.builds=[]
        for checked in (True,False):
            path,meta=build(checked=checked)
            cls.builds.append((path,meta,Scifort(path,dll_directories=meta['dll_directories'])))

    def test_scipy_weighted_bandwidth_and_dimension_references(self):
        for _,_,api in self.builds:
            for data in (DATA[0],DATA):
                for weights in (None,np.array([1.,2.,0.,3.,1.,2.])):
                    for bandwidth in (None,'scott','silverman',.7):
                        reference=gaussian_kde(data,bw_method=bandwidth,weights=weights)
                        with api.gaussian_kde(data,bw_method=bandwidth,weights=weights) as kde:
                            points=np.array([-.8,0.,2.]) if kde.d==1 else np.array([[0.,-1.,2.],[.5,1.5,-.5]])
                            np.testing.assert_allclose(kde(points),reference(points),rtol=8e-14,atol=1e-15)
                            np.testing.assert_allclose(kde.logpdf(points),reference.logpdf(points),rtol=8e-14,atol=2e-14)
                            self.assertEqual((kde.d,kde.n),(reference.d,reference.n))
                            self.assertAlmostEqual(kde.neff,reference.neff,places=14)
                            self.assertAlmostEqual(kde.factor,reference.factor,places=14)
                            self.assertAlmostEqual(kde.scotts_factor(),reference.scotts_factor(),places=14)
                            self.assertAlmostEqual(kde.silverman_factor(),reference.silverman_factor(),places=14)

    def test_layout_empty_batches_tail_and_output_buffers(self):
        for _,_,api in self.builds:
            data=DATA[:,::-1]
            with api.gaussian_kde(data) as kde:
                points=np.array([[0.,1.,2.,3.],[1.,-1.,2.,0.]])[:,::2]
                np.testing.assert_allclose(kde.pdf(points),gaussian_kde(data)(points),rtol=8e-14)
                self.assertEqual(kde.pdf([0.,1.]).shape,(1,))
                self.assertEqual(kde.pdf(np.empty((2,0))).shape,(0,))
                self.assertEqual(kde.logpdf(np.empty((2,0))).shape,(0,))
                out=np.empty(2);self.assertIs(kde.pdf(points,out=out),out)
                readonly=np.empty(2);readonly.flags.writeable=False
                unaligned=np.ndarray((2,),dtype=np.float64,buffer=bytearray(17),offset=1)
                for bad in (np.empty(3),np.empty(2,dtype=np.float32),np.empty(4)[::2],readonly,unaligned):
                    with self.assertRaises((TypeError,ValueError)):kde.pdf(points,out=bad)
                with self.assertRaises(ValueError):kde.pdf(points,out=points[0,:])
                for bad in (np.ones((3,2)),np.array([[np.nan],[0.]]),np.array([[np.inf],[0.]])):
                    with self.assertRaises(ValueError):kde.pdf(bad)
            with api.gaussian_kde(DATA[0]) as kde:
                self.assertEqual(kde.pdf(0.).shape,(1,))
                self.assertEqual(kde.pdf([]).shape,(0,))
                self.assertEqual(kde.pdf([1000.])[0],0.)
                self.assertTrue(np.isfinite(kde.logpdf([1000.])[0]))
                self.assertLess(kde.logpdf([1000.])[0],-1000.)

    def test_lifetime_copying_and_source_independence(self):
        for path,meta,api in self.builds:
            source=DATA.copy();weights=np.array([1.,2.,0.,3.,1.,2.])
            kde=api.gaussian_kde(source,weights=weights)
            expected=kde.pdf([[0.],[1.]])
            source[:]=999;weights[:]=0
            np.testing.assert_array_equal(kde.pdf([[0.],[1.]]),expected)
            for operation in (lambda:copy.copy(kde),lambda:copy.deepcopy(kde),lambda:pickle.dumps(kde)):
                with self.assertRaises(TypeError):operation()
            pointer=kde._handle._pointer;reference=weakref.ref(kde._handle)
            kde.close();kde.close();self.assertTrue(kde.closed);self.assertIsNone(pointer.value)
            for operation in (lambda:kde.pdf([[0.],[1.]]),lambda:kde.factor,lambda:kde.__enter__()):
                with self.assertRaises(ValueError):operation()
            del kde;gc.collect();self.assertIsNone(reference())
            owner=Scifort(path,dll_directories=meta['dll_directories'])
            kde=owner.gaussian_kde(DATA[0]);del owner;gc.collect()
            self.assertTrue(np.all(np.isfinite(kde.pdf([0.,1.]))))
            pointer=kde._handle._pointer;reference=weakref.ref(kde._handle)
            del kde;gc.collect();self.assertIsNone(reference());self.assertIsNone(pointer.value)
            with self.assertRaisesRegex(RuntimeError,'body error'):
                with api.gaussian_kde(DATA[0]) as kde:raise RuntimeError('body error')
            self.assertTrue(kde.closed)
            for _ in range(50):
                with api.gaussian_kde(DATA[0]) as kde:self.assertFalse(kde.closed)

    def test_invalid_initialization_and_handle_ownership(self):
        for path,meta,api in self.builds:
            for data in ([],[1.],[1.,1.,1.],[np.nan,1.,2.],np.ones((3,2)),np.ones((2,5)),np.ones((2,2,2))):
                with self.assertRaises(ValueError):api.gaussian_kde(data)
            for weights in ([1.,2.],[-1.,1.,1.,1.,1.,1.],np.zeros(6),[np.nan]*6,[1.,0.,0.,0.,0.,0.]):
                with self.assertRaises(ValueError):api.gaussian_kde(DATA,weights=weights)
            for bandwidth in (0.,-1.,np.nan,np.inf,'unknown',True,lambda obj:1.):
                with self.assertRaises((TypeError,ValueError)):api.gaussian_kde(DATA,bw_method=bandwidth)
            for data in ([1.+2j,2.,3.],['a','b','c']):
                with self.assertRaises(TypeError):api.gaussian_kde(data)
            owner=Scifort(path,dll_directories=meta['dll_directories'])
            with api.gaussian_kde(DATA) as kde:
                with self.assertRaises(ValueError):owner.kde_bridge.info(kde._handle)
                with self.assertRaises(TypeError):api.kde_bridge.info(C.c_void_p(123))

    def test_raw_c_errors_and_handle_lifecycle(self):
        ptr=lambda a:a.ctypes.data_as(C.c_void_p)
        for _,_,api in self.builds:
            raw=api.kde_bridge._functions;data=np.asfortranarray(DATA);weights=np.ones(6)
            handle=C.c_void_p();status=C.c_int(-1)
            raw['create'](2,6,ptr(data),ptr(weights),0,1.,C.byref(handle),C.byref(status))
            self.assertEqual(status.value,0);self.assertIsNotNone(handle.value)
            try:
                points=np.zeros((2,2),order='F');out=np.full(2,77.)
                for d,m in ((0,2),(2,-1)):
                    raw['pdf'](d,m,handle,ptr(points),ptr(out),C.byref(status))
                    self.assertEqual(status.value,1);np.testing.assert_array_equal(out,77.)
                raw['pdf'](1,2,handle,ptr(points),ptr(out),C.byref(status))
                self.assertEqual(status.value,1);self.assertTrue(np.all(np.isnan(out)))
                raw['pdf'](2,2,None,ptr(points),ptr(out),C.byref(status))
                self.assertEqual(status.value,5);self.assertTrue(np.all(np.isnan(out)))
            finally:raw['destroy'](C.byref(handle),C.byref(status))
            self.assertEqual(status.value,0);self.assertIsNone(handle.value)
            raw['destroy'](C.byref(handle),C.byref(status));self.assertEqual(status.value,0)
            for mode in (-1,3):
                raw['create'](2,6,ptr(data),ptr(weights),mode,1.,C.byref(handle),C.byref(status))
                self.assertEqual(status.value,2);self.assertIsNone(handle.value)
            bad=data.copy(order='F');bad[:]=1.
            raw['create'](2,6,ptr(bad),ptr(weights),0,1.,C.byref(handle),C.byref(status))
            self.assertEqual(status.value,3);self.assertIsNone(handle.value)

    def test_native_regressions_and_c_cpp_clients(self):
        for path,meta,_ in self.builds:
            with tempfile.TemporaryDirectory() as temp:
                folder=Path(temp);env=os.environ.copy();env['PATH']=str(path.parent)+os.pathsep+env.get('PATH','')
                flags=['-O0','-g','-fcheck=all','-fbacktrace'] if '-fcheck=all' in meta['flags'] else ['-O3']
                exe=folder/'native.exe'
                run([meta['compiler_path'],*flags,'-I',path.parent,ROOT/'test/test_gaussian_kde_c_api.f90',path,'-o',exe])
                run([exe],env=env)
                # Exercise the original native suite without pulling every scifort_stats export into the build.
                native=folder/'native_regression.f90'
                native.write_text((ROOT/'test/test_gaussian_kde.f90').read_text(encoding='utf-8').replace(
                    'use scifort_stats, only :','use scifort_gaussian_kde, only :'),encoding='ascii')
                run([meta['compiler_path'],*flags,'-I',path.parent,'-J',folder,ROOT/'test/gaussian_kde_reference.f90',native,path,'-o',exe])
                run([exe],env=env)
                client=folder/'client.c'
                client.write_text("""#include "scifort.h"
#include <math.h>
#include <stddef.h>
int main(void) {
 double data[5]={-2,-.5,.25,1.5,3},weights[5]={1,1,1,1,1},points[3]={-1,0,2};
 double out[3],logs[3],mean=0,var=0; void *handle=NULL; int status;
 for(int i=0;i<5;++i) mean+=data[i]/5;
 for(int i=0;i<5;++i) var+=(data[i]-mean)*(data[i]-mean)/4;
 double sd=sqrt(var)*pow(5.,-.2);
 scifort_gaussian_kde_create_f64(1,5,data,weights,0,1,&handle,&status);
 if(status || !handle) return 1;
 scifort_gaussian_kde_pdf_f64(1,3,handle,points,out,&status); if(status) return 2;
 scifort_gaussian_kde_logpdf_f64(1,3,handle,points,logs,&status); if(status) return 3;
 for(int j=0;j<3;++j) {
  double expected=0;
  for(int i=0;i<5;++i) {
   double z=(points[j]-data[i])/sd;
   expected+=exp(-.5*z*z)/(5*sd*sqrt(2*acos(-1)));
  }
  if(fabs(out[j]-expected)>1e-14 || fabs(logs[j]-log(expected))>1e-13) return 4;
 }
 scifort_gaussian_kde_destroy_f64(&handle,&status); if(status || handle) return 5;
 scifort_gaussian_kde_destroy_f64(&handle,&status); return status || handle!=NULL;
}
""",encoding='ascii')
                for compiler,standard in (('gcc','c11'),('g++','c++17')):
                    if not shutil.which(compiler):continue
                    exe=folder/(compiler+'.exe')
                    run([compiler,'-std='+standard,'-Wall','-Wextra','-Werror','-I',ROOT/'include',client,path,'-o',exe])
                    run([exe],env=env)
