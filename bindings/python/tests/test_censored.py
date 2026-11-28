# SPDX-License-Identifier: MIT
"""Censored likelihood adapters: native/C agreement, references and invalid inputs."""
import ctypes as C
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
import numpy as np
from scipy.stats import norm, expon, weibull_min
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from build_library import build
from scifort import Scifort
from test_multivariate_t import run
ROOT = Path(__file__).resolve().parents[3]


class CensoredLikelihoodTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran'):
            raise unittest.SkipTest('gfortran required')
        cls.builds = []
        for checked in (True, False):
            path, meta = build(checked=checked)
            cls.builds.append((path, meta, Scifort(path, dll_directories=meta['dll_directories'])))

    @staticmethod
    def reference(distribution, exact, left, right, interval):
        value = float(np.sum(distribution.logpdf(exact)) + np.sum(distribution.logcdf(left))
                      + np.sum(distribution.logsf(right)))
        for a, b in interval:
            if a == b:
                value += distribution.logpdf(a)
            elif np.isneginf(a):
                value += distribution.logcdf(b)
            elif np.isposinf(b):
                value += distribution.logsf(a)
            else:
                # Stable subtraction of public SciPy tail values, no private routines.
                upper = distribution.logcdf(b)
                lower = distribution.logsf(a)
                if upper <= lower:
                    value += upper + np.log(-np.expm1(distribution.logcdf(a) - upper))
                else:
                    value += lower + np.log(-np.expm1(distribution.logsf(b) - lower))
        return value

    def test_mixed_references_independent_lengths_and_nnlf(self):
        cases = (([.5,1.2,2.5],[.3],[3.5],[[1.1,1.4],[4.,4.4]]),
                 ([],[-1.,.2],[1.,2.,3.],[[.4,.4],[-np.inf,.5],[2.,np.inf]]),
                 ([.2,1.2],[],[],[]), ([],[],[],[[.2,.5]]))
        families = (('normal',norm(loc=.4,scale=1.3),dict(loc=.4,scale=1.3)),
                    ('exponential',expon(scale=2.3),dict(scale=2.3)),
                    ('weibull',weibull_min(1.7,scale=2.3),dict(c=1.7,scale=2.3)))
        for _,_,api in self.builds:
            for family, dist, params in families:
                for exact,left,right,interval in cases:
                    expected = self.reference(dist,exact,left,right,interval)
                    actual = getattr(api.censored,family+'_loglikelihood')(
                        exact,left=left,right=right,interval=interval,**params)
                    np.testing.assert_allclose(actual,expected,rtol=2e-14,atol=2e-13)
                    nll = getattr(api.censored,family+'_nnlf')(
                        exact,left=left,right=right,interval=interval,**params)
                    self.assertEqual(nll,-actual)

    def test_empty_defaults_support_and_weibull_identity(self):
        for _,_,api in self.builds:
            for family in ('normal','exponential','weibull'):
                params = dict(c=1.) if family=='weibull' else {}
                for method in ('loglikelihood','nnlf'):
                    fn = getattr(api.censored,family+'_'+method)
                    self.assertEqual(fn(**params),0.)
                    self.assertEqual(fn([],interval=np.empty((0,2)),**params),0.)
                empty = np.empty(0)
                low = getattr(api.censored_bridge,family+'_loglikelihood')(
                    empty,empty,empty,empty,empty,**params)
                self.assertEqual(low.status,0)
                self.assertEqual(low.outputs['value'],0.)
            exact=[.2,1.,2.]; kwargs=dict(left=[.1],right=[3.],interval=[[.4,.6]])
            self.assertAlmostEqual(api.censored.exponential_loglikelihood(exact,**kwargs),
                                   api.censored.weibull_loglikelihood(exact,c=1.,**kwargs),delta=4e-14)
            for family in ('exponential','weibull'):
                params = dict(c=2.) if family=='weibull' else {}
                self.assertEqual(getattr(api.censored,family+'_loglikelihood')([-1.],**params),-np.inf)
                self.assertEqual(getattr(api.censored,family+'_nnlf')([-1.],**params),np.inf)
                self.assertEqual(getattr(api.censored,family+'_loglikelihood')(right=[-1.],**params),0.)
                a=getattr(api.censored,family+'_loglikelihood')(interval=[[-1.,1.]],**params)
                b=getattr(api.censored,family+'_loglikelihood')(interval=[[0.,1.]],**params)
                self.assertEqual(a,b)
            self.assertAlmostEqual(api.censored.normal_loglikelihood([.2]),float(norm.logpdf(.2)),delta=2e-15)
            # Finite endpoint density may be +inf; opposite infinite terms retain native NaN.
            self.assertEqual(api.censored.weibull_loglikelihood([0.],c=.5),np.inf)
            self.assertTrue(np.isnan(api.censored.weibull_loglikelihood([0.,-1.],c=.5)))

    def test_narrow_and_extreme_intervals_saved_references(self):
        # Existing mpmath references from test/test_censored.f90; exact binary64 endpoints.
        cases=(('normal',{},[0.,1e-20],-46.97064039308559,2e-13),
               ('normal',{},[40.,40.01],-805.7174659453684,3e-11),
               ('exponential',{},[1000.,1000.+1e-10],-1023.0254068566554,3e-12),
               ('exponential',dict(scale=1e300),[0.,1e-300],-1381.5510557964274,2e-12),
               ('weibull',dict(c=2.),[100.,100.+1e-8],-10013.122365005156,3e-10),
               ('weibull',dict(c=.01),[1.,np.nextafter(1.,2.)],-41.64882357510525,3e-13))
        for _,_,api in self.builds:
            for family,params,bounds,expected,tolerance in cases:
                value=getattr(api.censored,family+'_loglikelihood')(interval=[bounds],**params)
                self.assertAlmostEqual(value,expected,delta=tolerance)

    def test_layout_aliases_equal_and_infinite_endpoints(self):
        for _,_,api in self.builds:
            x=np.arange(1.,9.)[::2]; x.flags.writeable=False
            bounds=np.array([[.2,.4],[1.,2.],[3.,3.],[4.,6.]])[::2]
            bounds.flags.writeable=False
            actual=api.censored.normal_loglikelihood(x,left=x,right=x,interval=bounds)
            expected=self.reference(norm(),x,x,x,bounds)
            self.assertAlmostEqual(actual,expected,delta=3e-13)
            np.testing.assert_array_equal(x,[1.,3.,5.,7.])
            a=api.censored.normal_loglikelihood(left=[-2.],right=[3.])
            b=api.censored.normal_loglikelihood(interval=[[-np.inf,-2.],[3.,np.inf]])
            self.assertEqual(a,b)
            a=api.censored.normal_loglikelihood([.2,1.2])
            b=api.censored.normal_loglikelihood(interval=[[.2,.2],[1.2,1.2]])
            self.assertEqual(a,b)
            # One-dimensional views from the same matrix may overlap before conversion.
            matrix=np.array([[.2,.4],[1.,2.]])
            result=api.censored.normal_loglikelihood(matrix[:,0],interval=matrix)
            expected=self.reference(norm(),matrix[:,0],[],[],matrix)
            self.assertAlmostEqual(result,expected,delta=2e-14)

    def test_invalid_data_parameters_and_shapes(self):
        for _,_,api in self.builds:
            for family in ('normal','exponential','weibull'):
                params = dict(c=1.7) if family=='weibull' else {}
                for method in ('loglikelihood','nnlf'):
                    fn=getattr(api.censored,family+'_'+method)
                    for bad in (0.,-1.,np.inf,-np.inf,np.nan):
                        with self.assertRaises(ValueError):fn(scale=bad,**params)
                    for bad in (np.nan,np.inf,-np.inf):
                        with self.assertRaises(ValueError):fn(loc=bad,**params)
                    for group in ('uncensored','left','right'):
                        for bad in ([np.nan],[np.inf],[-np.inf],[[1.]],1.,['a'],[1j]):
                            with self.assertRaises((ValueError,TypeError)):fn(**{group:bad},**params)
                    for bad in ([[2.,1.]],[[np.nan,2.]],[[1.,np.nan]],
                                [[-np.inf,np.inf]],[[np.inf,np.inf]],[[-np.inf,-np.inf]],
                                [1.,2.],[[1.,2.,3.]],np.empty((0,1))):
                        with self.assertRaises((ValueError,TypeError)):fn(interval=bad,**params)
                low=getattr(api.censored_bridge,family+'_loglikelihood')
                for params_bad in (dict(params,scale=0.),dict(params,loc=np.inf)):
                    result=low([],[],[],[],[],**params_bad)
                    self.assertEqual(result.status,2);self.assertTrue(np.isnan(result.outputs['value']))
                result=low([np.inf],[],[],[],[],**params)
                self.assertEqual(result.status,2);self.assertTrue(np.isnan(result.outputs['value']))
                with self.assertRaises(ValueError):low([],[],[],[1.],[2.,3.],**params)
            for bad in (0.,-1.,np.inf,-np.inf,np.nan):
                for method in ('loglikelihood','nnlf'):
                    with self.assertRaises(ValueError):getattr(api.censored,'weibull_'+method)(c=bad)

    def test_raw_c_status_and_output_preservation(self):
        for _,_,api in self.builds:
            for family in ('normal','exponential','weibull'):
                for method in ('loglikelihood','nnlf'):
                    fn=getattr(api.censored_bridge._library,f'scifort_{family}_censored_{method}_f64')
                    count=3 if family=='weibull' else 2
                    fn.argtypes=[C.c_int]*4+[C.c_void_p]*5+[C.c_double]*count+[C.POINTER(C.c_double),C.POINTER(C.c_int)]
                    fn.restype=None
                    exact=np.array([1.]);left=np.array([.2]);right=np.array([2.])
                    lower=np.array([1.]);upper=np.array([2.]);ptrs=[a.ctypes.data for a in (exact,left,right,lower,upper)]
                    params=[1.7,0.,1.] if family=='weibull' else [0.,1.]
                    value=C.c_double(77.);status=C.c_int()
                    for index in range(4):
                        dims=[1]*4;dims[index]=-1
                        fn(*dims,*ptrs,*params,C.byref(value),C.byref(status))
                        self.assertEqual(status.value,1);self.assertEqual(value.value,77.)
                    for a,b in ((2.,1.),(-np.inf,np.inf),(np.nan,2.)):
                        lower[0]=a;upper[0]=b
                        fn(1,1,1,1,*ptrs,*params,C.byref(value),C.byref(status))
                        self.assertEqual(status.value,2);self.assertTrue(np.isnan(value.value))

    def test_native_agreement_and_c_cpp_clients(self):
        for path,meta,_ in self.builds:
            with tempfile.TemporaryDirectory() as temporary:
                folder=Path(temporary);env=os.environ.copy()
                env['PATH']=str(path.parent)+os.pathsep+env.get('PATH','')
                flags=['-O0','-g','-fcheck=all','-fbacktrace'] if '-fcheck=all' in meta['flags'] else ['-O3']
                exe=folder/'native.exe'
                run([meta['compiler_path'],*flags,'-I',path.parent,ROOT/'test/test_censored_c_api.f90',path,'-o',exe])
                run([exe],env=env)
                client=folder/'client.c'
                client.write_text("""#include "scifort.h"
#include <math.h>
int main(void) {
 double x=1,l=.2,r=2,a=1,b=2,value,nll; int status;
 scifort_exponential_censored_loglikelihood_f64(1,0,0,0,&x,&l,&r,&a,&b,0,1,&value,&status);
 if(status || value!=-1) return 1;
 scifort_exponential_censored_nnlf_f64(1,0,0,0,&x,&l,&r,&a,&b,0,1,&nll,&status);
 if(status || nll!=1) return 2;
 scifort_normal_censored_loglikelihood_f64(0,0,0,0,&x,&l,&r,&a,&b,0,1,&value,&status);
 if(status || value!=0) return 3;
 scifort_weibull_censored_loglikelihood_f64(0,0,0,0,&x,&l,&r,&a,&b,0,0,1,&value,&status);
 return status!=2 || !isnan(value);
}
""",encoding='ascii')
                for compiler,standard in (('gcc','c11'),('g++','c++17')):
                    if not shutil.which(compiler):continue
                    exe=folder/(compiler+'.exe')
                    run([compiler,'-std='+standard,'-Wall','-Wextra','-Werror','-I',ROOT/'include',client,path,'-o',exe])
                    run([exe],env=env)


if __name__ == '__main__':
    unittest.main()
