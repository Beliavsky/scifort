# SPDX-License-Identifier: MIT
"""KDE integrals, explicit streams and independently owned marginals."""
import ctypes as C
import gc
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
import weakref
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from scipy.stats import gaussian_kde
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from build_library import build
from scifort import Scifort
from test_multivariate_t import run
ROOT=Path(__file__).resolve().parents[3]
DATA=np.array([[-2.,-.5,.25,1.5,3.,.8],[.5,2.,-1.,.2,2.8,-.4]])
WEIGHTS=np.array([1.,2.,0.,3.,1.,2.])

class GaussianKDEMethodTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran'):raise unittest.SkipTest('gfortran required')
        cls.builds=[]
        for checked in (True,False):
            path,meta=build(checked=checked)
            cls.builds.append((path,meta,Scifort(path,dll_directories=meta['dll_directories'])))

    def test_covariance_bandwidth_updates_and_failed_updates(self):
        for _,_,api in self.builds:
            for data in (DATA[0],DATA):
                ref=gaussian_kde(data,weights=WEIGHTS)
                with api.gaussian_kde(data,weights=WEIGHTS) as kde:
                    for method in ('silverman',.35,None,'scott',.8):
                        kde.set_bandwidth(method);ref.set_bandwidth(method)
                        self.assertAlmostEqual(kde.factor,ref.factor,delta=2e-15)
                        np.testing.assert_allclose(kde.covariance,ref.covariance,rtol=4e-15,atol=1e-16)
                        np.testing.assert_allclose(kde.inv_cov,ref.inv_cov,rtol=5e-15,atol=1e-15)
                        np.testing.assert_allclose(kde.pdf(data),ref.pdf(data),rtol=5e-15,atol=1e-16)
                        np.testing.assert_allclose(kde.covariance@kde.inv_cov,np.eye(kde.d),atol=1e-15)
                    saved=kde.covariance;factor=kde.factor;values=kde.pdf(data)
                    for bad in (0.,-1.,np.nan,np.inf,1e308,1e-300,'bad',True,[.5],1j):
                        with self.assertRaises((ValueError,TypeError)):kde.set_bandwidth(bad)
                        self.assertEqual(kde.factor,factor)
                        np.testing.assert_array_equal(kde.covariance,saved)
                        np.testing.assert_array_equal(kde.pdf(data),values)
                    saved[:]=0.
                    self.assertTrue(np.any(kde.covariance!=0.))
                    for getter in (kde.get_covariance,kde.get_inv_cov):
                        out=np.empty((kde.d,kde.d),order='F')
                        self.assertIs(getter(out=out),out)
                        for badout in (np.zeros((kde.d+1,kde.d)),np.zeros((kde.d,kde.d),dtype=np.float32)):
                            with self.assertRaises((ValueError,TypeError)):getter(out=badout)
                        out.flags.writeable=False
                        with self.assertRaises((ValueError,TypeError)):getter(out=out)
                    for method in ('covariance','inv_cov'):
                        wrong=getattr(api.kde_bridge,method)(kde._handle,d=kde.d+1)
                        self.assertEqual(wrong.status,1)
                        self.assertTrue(np.isnan(wrong.outputs['matrix']).all())
                    self.assertEqual(api.kde_bridge.set_bandwidth(kde._handle,9,1.).status,2)
                for operation in (lambda:kde.covariance,lambda:kde.inv_cov,lambda:kde.set_bandwidth(.5)):
                    with self.assertRaises(ValueError):operation()

    def test_bandwidth_marginal_independence(self):
        for _,_,api in self.builds:
            with api.gaussian_kde(DATA,weights=WEIGHTS) as parent:
                with parent.marginal(0) as child:
                    before=child.covariance;factor=child.factor
                    parent.set_bandwidth(.3)
                    self.assertEqual(child.factor,factor)
                    np.testing.assert_array_equal(child.covariance,before)
                    child.set_bandwidth('silverman')
                    self.assertEqual(parent.factor,.3)
                    parent.close()
                    child.set_bandwidth(.6)
                    self.assertEqual(child.factor,.6)
                    self.assertTrue(np.isfinite(child.inv_cov).all())

    def test_gaussian_and_interval_integrals(self):
        for _,_,api in self.builds:
            for data in (DATA[0],DATA):
                ref=gaussian_kde(data,weights=WEIGHTS)
                with api.gaussian_kde(data,weights=WEIGHTS) as kde:
                    mean=np.arange(kde.d)*.2
                    cov=np.eye(kde.d)*.7
                    self.assertAlmostEqual(kde.integrate_gaussian(mean,cov),ref.integrate_gaussian(mean,cov),delta=2e-14)
                    if kde.d==1:
                        self.assertAlmostEqual(kde.integrate_gaussian(0.,.7),ref.integrate_gaussian([0.],[[.7]]),delta=2e-14)
                        for low,high in ((-np.inf,np.inf),(-1.,2.),(2.,-1.),(0.,0.),(np.inf,np.inf)):
                            self.assertAlmostEqual(kde.integrate_box_1d(low,high),ref.integrate_box_1d(low,high),delta=3e-15)
                        self.assertTrue(np.isnan(kde.integrate_box_1d(np.nan,1.)))
                    else:
                        with self.assertRaises(ValueError):kde.integrate_box_1d(-1.,1.)
                    for badmean,badcov in (([np.nan]*kde.d,cov),(mean,np.full_like(cov,np.nan)),
                                            (mean,-100.*np.eye(kde.d)),(np.zeros(kde.d+1),cov)):
                        with self.assertRaises(ValueError):kde.integrate_gaussian(badmean,badcov)

    def test_box_qmc_references_snapshots_and_signed_limits(self):
        for _,_,api in self.builds:
            with api.gaussian_kde(DATA[0],weights=WEIGHTS) as kde:
                state=api.rng(42);saved=state.state
                self.assertAlmostEqual(kde.integrate_box(-1.,2.,state=state),kde.integrate_box_1d(-1.,2.),delta=3e-15)
                np.testing.assert_array_equal(state.state,saved)
            with api.gaussian_kde(DATA,weights=WEIGHTS) as kde:
                low=np.array([-1.,-.8]);high=np.array([1.5,1.2]);state=api.rng(42);saved=state.state
                value=kde.integrate_box(low,high,state=state,maxpts=16000);final=state.state
                self.assertFalse(np.array_equal(saved,final))
                state.set_state(saved)
                self.assertEqual(kde.integrate_box(low,high,state=state,maxpts=16000),value)
                np.testing.assert_array_equal(state.state,final)
                reference=gaussian_kde(DATA,weights=WEIGHTS).integrate_box(low,high,maxpts=16000)
                self.assertAlmostEqual(value,reference,delta=1e-3)
                reversed_low=np.array([high[0],low[1]]);reversed_high=np.array([low[0],high[1]])
                state.set_state(saved)
                self.assertAlmostEqual(kde.integrate_box(reversed_low,reversed_high,state=state,maxpts=16000),-value,delta=3e-15)
                self.assertAlmostEqual(kde.integrate_box([-np.inf]*2,[np.inf]*2,state=state,maxpts=128),1.,delta=2e-15)
                self.assertEqual(kde.integrate_box(low,low,state=state,maxpts=128),0.)
                saved=state.state
                for kwargs in ({'maxpts':-1},{'maxpts':True},{'maxpts':1.5},{'maxpts':2**31}):
                    with self.assertRaises((TypeError,ValueError)):kde.integrate_box(low,high,state=state,**kwargs)
                for badlow,badhigh in (([0.],high),([np.nan,0.],high),(low,[0.])):
                    with self.assertRaises(ValueError):kde.integrate_box(badlow,badhigh,state=state,maxpts=128)
                np.testing.assert_array_equal(state.state,saved)

    def test_two_handle_products_and_lock_order(self):
        for path,meta,api in self.builds:
            refa=gaussian_kde(DATA,weights=WEIGHTS)
            refb=gaussian_kde(DATA*.7+.2,bw_method=.5)
            with api.gaussian_kde(DATA,weights=WEIGHTS) as a, api.gaussian_kde(DATA*.7+.2,bw_method=.5) as b:
                expected=refa.integrate_kde(refb)
                self.assertAlmostEqual(a.integrate_kde(b),expected,delta=2e-14)
                self.assertAlmostEqual(a.integrate_kde(a),refa.integrate_kde(refa),delta=2e-14)
                with ThreadPoolExecutor(max_workers=4) as pool:
                    futures=[pool.submit(a.integrate_kde,b) if i%2 else pool.submit(b.integrate_kde,a) for i in range(16)]
                    for future in futures:self.assertAlmostEqual(future.result(timeout=10),expected,delta=2e-14)
                with self.assertRaises(TypeError):a.integrate_kde(1.)
                with api.gaussian_kde(DATA[0]) as one:
                    with self.assertRaises(ValueError):a.integrate_kde(one)
                owner=Scifort(path,dll_directories=meta['dll_directories'])
                with owner.gaussian_kde(DATA) as foreign:
                    with self.assertRaises(ValueError):a.integrate_kde(foreign)
                b.close()
                with self.assertRaises(ValueError):a.integrate_kde(b)
            with self.assertRaises(ValueError):a.integrate_kde(a)

    def test_resampling_streams_empty_default_and_buffers(self):
        for path,meta,api in self.builds:
            with api.gaussian_kde(DATA,weights=WEIGHTS) as kde:
                state=api.rng(123);split=api.rng(123);uniforms=api.rng(123)
                whole=kde.resample(state,7)
                first=kde.resample(split,2);rest=kde.resample(split,5)
                np.testing.assert_array_equal(whole,np.concatenate([first,rest],axis=1))
                np.testing.assert_array_equal(state.state,split.state)
                uniforms.uniform(3*7)
                np.testing.assert_array_equal(state.state,uniforms.state)
                default=kde.resample(api.rng(3));self.assertEqual(default.shape,(2,int(kde.neff)))
                saved=state.state
                self.assertEqual(kde.resample(state,0).shape,(2,0))
                np.testing.assert_array_equal(state.state,saved)
                out=np.empty((2,4),order='F');self.assertIs(kde.resample(state,4,out=out),out)
                for size in (-1,True,1.5,2**31):
                    with self.assertRaises((TypeError,ValueError)):kde.resample(state,size)
                readonly=np.empty((2,4),order='F');readonly.flags.writeable=False
                for bad in (np.empty((2,4)),np.empty((2,4),dtype=np.float32,order='F'),
                            np.empty((3,4),order='F'),readonly):
                    with self.assertRaises((TypeError,ValueError)):kde.resample(state,4,out=bad)
                with self.assertRaises(ValueError):
                    kde.resample(state,3,out=state._state.view(np.float64).reshape((2,3),order='F'))
                owner=Scifort(path,dll_directories=meta['dll_directories'])
                with self.assertRaises(TypeError):kde.resample(owner.rng(),1)
                with self.assertRaises(TypeError):kde.resample([1]*6,1)

    def test_resampling_mixture_moments(self):
        ref=gaussian_kde(DATA,weights=WEIGHTS)
        w=WEIGHTS/WEIGHTS.sum();mean=DATA@w;dev=DATA-mean[:,None]
        covariance=(dev*w)@dev.T+ref.covariance
        for _,_,api in self.builds:
            with api.gaussian_kde(DATA,weights=WEIGHTS) as kde:
                draws=kde.resample(api.rng(12345),18000)
                np.testing.assert_allclose(draws.mean(axis=1),mean,rtol=0,atol=.07)
                np.testing.assert_allclose(np.cov(draws),covariance,rtol=0,atol=.14)

    def test_marginal_coordinates_bandwidth_and_independence(self):
        for _,_,api in self.builds:
            reference=gaussian_kde(DATA,weights=WEIGHTS)
            parent=api.gaussian_kde(DATA,weights=WEIGHTS)
            for axes in ([1],[1,0],[-1],[0,-1],0):
                child=parent.marginal(axes);ref=reference.marginal(axes)
                with child:
                    points=np.array([-.5,0.,1.]) if child.d==1 else np.array([[0.,1.],[1.,-.5]])
                    np.testing.assert_allclose(child.pdf(points),ref.pdf(points),rtol=8e-14,atol=1e-15)
                    self.assertEqual(child.factor,parent.factor)
                    self.assertEqual(child.n,parent.n)
                    self.assertAlmostEqual(child.neff,parent.neff,places=14)
                    if child.d==1:self.assertNotEqual(child.factor,child.scotts_factor())
            child=parent.marginal([1]);grandchild=child.marginal([0]);parent_ref=weakref.ref(parent)
            parent.close();del parent;gc.collect();self.assertIsNone(parent_ref())
            np.testing.assert_allclose(child.pdf([0.,1.]),reference.marginal([1]).pdf([0.,1.]),rtol=8e-14)
            child.close()
            self.assertTrue(np.all(np.isfinite(grandchild.pdf([0.,1.]))));grandchild.close()
            with api.gaussian_kde(DATA) as parent:
                for axes in ([],[0,0],[-1,1],[2],[-3],[0.,1.],[True],[2**63],np.array([2**64-1],dtype=np.uint64)):
                    with self.assertRaises((TypeError,ValueError)):parent.marginal(axes)
                self.assertTrue(np.all(np.isfinite(parent.pdf([[0.],[0.]]))))
            with self.assertRaises(ValueError):parent.marginal([0])

    def test_raw_c_error_outputs_state_and_new_ownership(self):
        ptr=lambda a:a.ctypes.data_as(C.c_void_p)
        for _,_,api in self.builds:
            with api.gaussian_kde(DATA,weights=WEIGHTS) as kde:
                raw=api.kde_bridge._functions;h=kde._handle._pointer
                state=api.rng(33).state.copy();saved=state.copy();out=np.full((2,3),77.,order='F');status=C.c_int(-1)
                raw['resample'](2,-1,h,ptr(state),ptr(out),C.byref(status))
                self.assertEqual(status.value,1);np.testing.assert_array_equal(out,77.);np.testing.assert_array_equal(state,saved)
                raw['resample'](1,3,h,ptr(state),ptr(out),C.byref(status))
                self.assertEqual(status.value,1);self.assertTrue(np.all(np.isnan(out.ravel(order='F')[:3])))
                np.testing.assert_array_equal(state,saved)
                out[:]=77.;state[:]=0
                for size in (0,3):
                    raw['resample'](2,size,h,ptr(state),ptr(out),C.byref(status))
                    self.assertEqual(status.value,7);np.testing.assert_array_equal(out,77.);np.testing.assert_array_equal(state,0)
                value=C.c_double(77.);low=np.zeros(2);high=np.ones(2)
                raw['integrate_box'](2,h,ptr(low),ptr(high),128,ptr(state),C.byref(value),C.byref(status))
                self.assertEqual(status.value,7);self.assertEqual(value.value,77.)
                state[:]=saved
                raw['integrate_box'](2,h,ptr(low),ptr(high),-1,ptr(state),C.byref(value),C.byref(status))
                self.assertEqual(status.value,2);self.assertTrue(np.isnan(value.value));np.testing.assert_array_equal(state,saved)
                child=C.c_void_p();axes=np.array([2],dtype=np.intc)
                raw['marginal'](1,h,ptr(axes),C.byref(child),C.byref(status))
                self.assertEqual(status.value,0);self.assertIsNotNone(child.value)
                raw['destroy'](C.byref(child),C.byref(status));self.assertIsNone(child.value)
                for bad in (np.array([0],dtype=np.intc),np.array([1,1],dtype=np.intc)):
                    raw['marginal'](bad.size,h,ptr(bad),C.byref(child),C.byref(status))
                    self.assertEqual(status.value,2);self.assertIsNone(child.value)
                raw['marginal'](1,None,ptr(axes),C.byref(child),C.byref(status))
                self.assertEqual(status.value,5);self.assertIsNone(child.value)
                raw['integrate_kde'](h,None,C.byref(value),C.byref(status))
                self.assertEqual(status.value,5);self.assertTrue(np.isnan(value.value))

    def test_native_and_c_cpp_method_clients(self):
        for path,meta,_ in self.builds:
            with tempfile.TemporaryDirectory() as temp:
                folder=Path(temp);env=os.environ.copy();env['PATH']=str(path.parent)+os.pathsep+env.get('PATH','')
                flags=['-O0','-g','-fcheck=all','-fbacktrace'] if '-fcheck=all' in meta['flags'] else ['-O3']
                exe=folder/'native.exe'
                run([meta['compiler_path'],*flags,'-I',path.parent,ROOT/'test/test_gaussian_kde_methods_c_api.f90',path,'-o',exe])
                run([exe],env=env)
                client=folder/'client.c'
                client.write_text("""#include "scifort.h"
#include <math.h>
#include <stddef.h>
int main(void) {
 double data[5]={-2,-.5,.25,1.5,3},w[5]={1,1,1,1,1},v,mean=0,cov=1,lo=-1,hi=2,out[3];
 void *parent=NULL,*child=NULL; int status,axes[1]={1}; int64_t state[6];
 scifort_gaussian_kde_create_f64(1,5,data,w,0,1,&parent,&status); if(status) return 1;
 scifort_gaussian_kde_integrate_box_1d_f64(parent,-INFINITY,INFINITY,&v,&status);
 if(status || fabs(v-1)>1e-14) return 2;
 scifort_gaussian_kde_integrate_gaussian_f64(1,parent,&mean,&cov,&v,&status);
 if(status || !(v>0)) return 3;
 scifort_gaussian_kde_integrate_kde_f64(parent,parent,&v,&status); if(status || !(v>0)) return 4;
 scifort_rng_seed_i64(123,state,&status); if(status) return 5;
 scifort_gaussian_kde_integrate_box_f64(1,parent,&lo,&hi,128,state,&v,&status);
 if(status || !(v>0 && v<1)) return 6;
 scifort_gaussian_kde_resample_f64(1,3,parent,state,out,&status); if(status) return 7;
 for(int i=0;i<3;++i) if(!isfinite(out[i])) return 8;
 scifort_gaussian_kde_marginal_f64(1,parent,axes,&child,&status); if(status || !child) return 9;
 scifort_gaussian_kde_set_bandwidth_f64(parent,2,.5,&status); if(status) return 12;
 scifort_gaussian_kde_covariance_f64(1,parent,&cov,&status); if(status || !(cov>0)) return 13;
 scifort_gaussian_kde_inv_cov_f64(1,parent,&v,&status);
 if(status || fabs(v*cov-1)>1e-14) return 14;
 scifort_gaussian_kde_covariance_f64(1,NULL,&v,&status);
 if(status!=5 || !isnan(v)) return 15;
 v=77; scifort_gaussian_kde_covariance_f64(0,parent,&v,&status);
 if(status!=1 || v!=77) return 16;
 scifort_gaussian_kde_set_bandwidth_f64(NULL,0,1,&status); if(status!=5) return 17;
 scifort_gaussian_kde_destroy_f64(&parent,&status); if(status || parent) return 10;
 scifort_gaussian_kde_integrate_box_1d_f64(child,-INFINITY,INFINITY,&v,&status);
 if(status || fabs(v-1)>1e-14) return 11;
 scifort_gaussian_kde_destroy_f64(&child,&status); return status || child!=NULL;
}
""",encoding='ascii')
                for compiler,standard in (('gcc','c11'),('g++','c++17')):
                    if not shutil.which(compiler):continue
                    exe=folder/(compiler+'.exe')
                    run([compiler,'-std='+standard,'-Wall','-Wextra','-Werror','-I',ROOT/'include',client,path,'-o',exe])
                    run([exe],env=env)
