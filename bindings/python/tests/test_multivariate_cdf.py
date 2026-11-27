# SPDX-License-Identifier: MIT
"""CDF array/state contracts, analytic identities and randomized reproducibility."""
import ctypes as C
import itertools
import math
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
import numpy as np
from scipy.stats import norm, t
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from build_library import build
from scifort import Scifort
from test_multivariate_t import run
ROOT=Path(__file__).resolve().parents[3]

class MultivariateCDFTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran'):raise unittest.SkipTest('gfortran required')
        cls.builds=[]
        for checked in (True,False):
            path,metadata=build(checked=checked)
            cls.builds.append((path,metadata,Scifort(path,dll_directories=metadata['dll_directories'])))

    def test_univariate_limits_signed_boxes_and_stable_log_tails(self):
        for _,_,api in self.builds:
            for upper,lower in ((.5,-np.inf),(.5,-1.),(-1.,.5),(np.inf,-np.inf),(np.inf,np.inf)):
                for student in (False,True):
                    state=api.rng(8);saved=state.state
                    group=api.multivariate_t if student else api.multivariate_normal
                    kwargs={'df':5} if student else {}
                    actual=group.cdf(state,[upper],lower_limit=[lower],**kwargs)
                    reference=(t.cdf(upper,5)-t.cdf(lower,5) if student else norm.cdf(upper)-norm.cdf(lower))
                    np.testing.assert_allclose(actual,reference,rtol=2e-14,atol=2e-15)
                    np.testing.assert_array_equal(state.state,saved)
            self.assertTrue(math.isnan(api.multivariate_normal.logcdf(api.rng(8),[-1.],lower_limit=[.5])))
            state=api.rng(8);saved=state.state
            logp=api.multivariate_normal.logcdf(state,[-40.,-35.])
            self.assertTrue(math.isfinite(logp));self.assertAlmostEqual(logp,norm.logcdf(-40)+norm.logcdf(-35),places=11)
            np.testing.assert_array_equal(state.state,saved)
            self.assertEqual(api.multivariate_normal.cdf(api.rng(8),[-40.,-35.]),0.)
            self.assertEqual(api.multivariate_normal.logcdf(api.rng(8),[-np.inf]),-np.inf)

    def test_analytic_orthants_independence_monotonicity_and_singular_boxes(self):
        rho=.6;matrix=np.array([[1.,rho],[rho,1.]])
        exact=.25+math.asin(rho)/(2*math.pi)
        for _,_,api in self.builds:
            for student in (False,True):
                group=api.multivariate_t if student else api.multivariate_normal
                kwargs={'df':5} if student else {}
                actual=group.cdf(api.rng(123),[0.,0.],None,matrix,maxpts=16000,**kwargs)
                self.assertAlmostEqual(actual,exact,delta=1e-3 if student else 3e-4)
                low=group.cdf(api.rng(123),[-1.,-1.],None,matrix,maxpts=4000,**kwargs)
                high=group.cdf(api.rng(123),[1.,1.],None,matrix,maxpts=4000,**kwargs)
                self.assertTrue(0<=low<actual<high<=1)
                self.assertTrue(math.isnan(group.cdf(api.rng(123),[0.,0.],None,np.ones((2,2)),**kwargs)))
                singular=group.cdf(api.rng(123),[0.,0.],None,np.ones((2,2)),allow_singular=True,maxpts=4000,**kwargs)
                self.assertAlmostEqual(singular,.5,delta=1e-3)
                self.assertEqual(group.cdf(api.rng(123),[0.],None,[[0.]],allow_singular=True,**kwargs),1.)
            expected=np.prod(norm.cdf([.2,.8])-norm.cdf([-1.,-.5]))
            actual=api.multivariate_normal.cdf(api.rng(4),[.2,.8],lower_limit=[-1.,-.5],maxpts=2000)
            self.assertAlmostEqual(actual,expected,delta=3e-15)

    def test_optional_arrays_snapshots_normal_limit_and_qmc_state(self):
        x=[.5,.3];mean=np.array([.1,-.2]);matrix=np.array([[1.,.4],[.4,2.]])
        lower=np.array([-1.,-.5])
        for _,_,api in self.builds:
            for group,method in ((api.multivariate_normal,'cdf'),(api.multivariate_normal,'logcdf'),(api.multivariate_t,'cdf')):
                fn=getattr(group,method);student=group is api.multivariate_t
                for flags in itertools.product((False,True),repeat=3):
                    location,spread,low=[v if f else None for v,f in zip((mean,matrix,lower),flags)]
                    kwargs={'df':5} if student else {}
                    state=api.rng(777);snapshot=state.state
                    value=fn(state,x,location,spread,lower_limit=low,maxpts=4000,**kwargs)
                    final=state.state;state.set_state(snapshot)
                    self.assertEqual(fn(state,x,location,spread,lower_limit=low,maxpts=4000,**kwargs),value)
                    np.testing.assert_array_equal(state.state,final)
            a=api.rng(19);b=api.rng(19)
            normal=api.multivariate_normal.cdf(a,x,mean,matrix,maxpts=4000)
            student=api.multivariate_t.cdf(b,x,mean,matrix,df=np.inf,maxpts=4000)
            self.assertEqual(normal,student);np.testing.assert_array_equal(a.state,b.state)
            self.assertFalse(np.array_equal(a.state,api.rng(19).state))
            a=api.rng(19);b=api.rng(19)
            p=api.multivariate_normal.cdf(a,x,mean,matrix,maxpts=4000)
            lp=api.multivariate_normal.logcdf(b,x,mean,matrix,maxpts=4000)
            self.assertAlmostEqual(lp,math.log(p),places=14);np.testing.assert_array_equal(a.state,b.state)
            a=api.rng(19);b=api.rng(19)
            p=api.multivariate_normal.cdf(a,x,mean,matrix,maxpts=4000,releps=1e-20)
            q=api.multivariate_normal.cdf(b,x,mean,matrix,maxpts=4000,releps=np.inf)
            self.assertEqual(p,q);np.testing.assert_array_equal(a.state,b.state)

    def test_invalid_parameters_shapes_and_transport_guards(self):
        for _,_,api in self.builds:
            for group,method in ((api.multivariate_normal,'cdf'),(api.multivariate_normal,'logcdf'),(api.multivariate_t,'cdf')):
                fn=getattr(group,method);state=api.rng(7);saved=state.state
                for kwargs in ({'lower_limit':[np.nan,0.]},{'allow_singular':True}):
                    if 'allow_singular' in kwargs:args=([0.,0.],None,[[1.,2.],[2.,1.]])
                    else:args=([0.,0.],)
                    self.assertTrue(math.isnan(fn(state,*args,**kwargs)))
                    np.testing.assert_array_equal(state.state,saved)
                self.assertTrue(math.isnan(fn(state,[np.nan,0.])))
                if group is api.multivariate_normal:
                    for name in ('abseps','releps'):
                        for bad in (0.,-1.,np.nan):
                            self.assertTrue(math.isnan(fn(state,[0.,0.],**{name:bad})))
                    self.assertTrue(math.isnan(fn(state,[0.,0.],abseps=np.inf)))
                else:
                    for bad in (0.,-1.,np.nan,-np.inf):self.assertTrue(math.isnan(fn(state,[0.,0.],df=bad)))
                np.testing.assert_array_equal(state.state,saved)
                for args in (([],),([[1.,2.]],),([1.,2.],[0.]),([1.,2.],None,np.eye(3))):
                    with self.assertRaises((ValueError,TypeError)):fn(state,*args)
                for kwargs in ({'lower_limit':[0.]},{'maxpts':2**31},{'maxpts':1.5},{'allow_singular':1}):
                    with self.assertRaises((ValueError,TypeError)):fn(state,[0.,0.],**kwargs)
                x=np.array([0.,0.])
                with self.assertRaises(ValueError):fn(state,x,lower_limit=x)
                with self.assertRaises(TypeError):fn([1,2,3,4,5,6],[0.,0.])
                np.testing.assert_array_equal(state.state,saved)

    def test_raw_c_results_and_guard_preservation(self):
        ptr=lambda a:a.ctypes.data_as(C.c_void_p)
        for _,_,api in self.builds:
            for group,method in ((api.multivariate_normal,'cdf'),(api.multivariate_normal,'logcdf'),(api.multivariate_t,'cdf')):
                fn=getattr(group,method);student=group is api.multivariate_t
                raw=group._functions[method];rng=api.rng(12);state=rng.state.copy()
                x=np.array([.5,.3]);mean=np.array([.1,-.2]);matrix=np.array([[1.,.4],[.4,2.]],order='F');low=np.array([-1.,-.5])
                value=C.c_double(77.);status=C.c_int(-1)
                def call(n=2,flag=1,singular=0):
                    args=[n,ptr(x),flag,ptr(mean),1,ptr(matrix)]
                    args += [5.,singular,4000] if student else [singular,4000,1e-5,1e-5]
                    raw(*args,1,ptr(low),ptr(state),C.byref(value),C.byref(status))
                call();kwargs={'df':5} if student else {}
                self.assertEqual(status.value,0);self.assertEqual(value.value,fn(rng,x,mean,matrix,maxpts=4000,lower_limit=low,**kwargs))
                np.testing.assert_array_equal(state,rng.state);saved=state.copy();value.value=77
                for kwargs in ({'n':0},{'flag':2},{'singular':2}):
                    call(**kwargs);self.assertEqual(status.value,1);self.assertEqual(value.value,77)
                    np.testing.assert_array_equal(state,saved)
                state[:]=0;call();self.assertEqual(status.value,4);self.assertEqual(value.value,77)
                np.testing.assert_array_equal(state,0)

    def test_native_c_and_cpp_clients(self):
        for path,metadata,_ in self.builds:
            with tempfile.TemporaryDirectory() as temp:
                folder=Path(temp);env=os.environ.copy();env['PATH']=str(path.parent)+os.pathsep+env.get('PATH','')
                flags=['-O0','-g','-fcheck=all','-fbacktrace'] if '-fcheck=all' in metadata['flags'] else ['-O3']
                exe=folder/'native.exe'
                run([metadata['compiler_path'],*flags,'-I',path.parent,ROOT/'test/test_multivariate_cdf_c_api.f90',path,'-o',exe]);run([exe],env=env)
                if shutil.which('gcc'):
                    client=folder/'client.c';client.write_text('''#include "scifort.h"
#include <math.h>
int main(void) {
 int64_t state[6]; double x[1]={0},mu[1]={0},cov[1]={1},low[1]={-INFINITY},value; int status;
 scifort_rng_seed_i64(123,state,&status); if(status) return 1;
 scifort_multivariate_normal_cdf_f64(1,x,0,mu,0,cov,0,0,1e-5,1e-5,0,low,state,&value,&status);
 if(status || fabs(value-.5)>1e-15) return 2;
 scifort_multivariate_normal_logcdf_f64(1,x,0,mu,0,cov,0,0,1e-5,1e-5,0,low,state,&value,&status);
 if(status || fabs(value-log(.5))>1e-15) return 3;
 scifort_multivariate_t_cdf_f64(1,x,0,mu,0,cov,5,0,0,0,low,state,&value,&status);
 return status || fabs(value-.5)>1e-15;
}
''',encoding='ascii')
                    exe=folder/'c.exe';run(['gcc','-std=c11','-Wall','-Wextra','-Werror','-I',ROOT/'include',client,path,'-o',exe]);run([exe],env=env)
                    if shutil.which('g++'):run(['g++','-std=c++17','-Wall','-Wextra','-Werror','-I',ROOT/'include','-fsyntax-only',client])
