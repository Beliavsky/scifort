# SPDX-License-Identifier: MIT
"""Lognormal/gamma/beta censored interfaces and parameter shape contracts."""
import ctypes as C
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
import numpy as np
from scipy.stats import lognorm, gamma, beta
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from build_library import build
from scifort import Scifort
import test_censored as helpers
from test_multivariate_t import run
ROOT = Path(__file__).resolve().parents[3]
FAMILIES = (('lognormal',dict(s=.8)), ('gamma',dict(a=.6)), ('beta',dict(a=.4,b=.6)))


class CensoredShapeFamilyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran'):
            raise unittest.SkipTest('gfortran required')
        cls.builds=[]
        for checked in (True,False):
            path,meta=build(checked=checked)
            cls.builds.append((path,meta,Scifort(path,dll_directories=meta['dll_directories'])))

    def test_mixed_references_layout_and_nnlf(self):
        cases=(([.1,.3,.6],[.05],[.8],[[.2,.25],[.7,.75]]),
               ([],[.05,.1],[.8,.9,1.],[[.4,.4],[-np.inf,.5],[.8,np.inf]]),
               ([.2,.4],[],[],[]), ([],[],[],[[.2,.5]]))
        families=(('lognormal',dict(s=.8,loc=-.1,scale=2.3),lognorm(.8,loc=-.1,scale=2.3)),
                  ('gamma',dict(a=.6,loc=-.1,scale=2.3),gamma(.6,loc=-.1,scale=2.3)),
                  ('beta',dict(a=.4,b=.6,loc=-.1,scale=2.),beta(.4,.6,loc=-.1,scale=2.)))
        for _,_,api in self.builds:
            for family,params,dist in families:
                for exact,left,right,interval in cases:
                    expected=helpers.CensoredLikelihoodTests.reference(dist,exact,left,right,interval)
                    fn=getattr(api.censored,family+'_loglikelihood')
                    actual=fn(exact,left=left,right=right,interval=interval,**params)
                    np.testing.assert_allclose(actual,expected,rtol=3e-13,atol=3e-12)
                    nll=getattr(api.censored,family+'_nnlf')(exact,left=left,right=right,interval=interval,**params)
                    self.assertEqual(nll,-actual)
                x=np.array([.1,99.,.3,99.,.6])[::2];x.flags.writeable=False
                bounds=np.array([[.2,.25],[99.,99.],[.7,.75]])[::2];bounds.flags.writeable=False
                expected=helpers.CensoredLikelihoodTests.reference(dist,x,x,x,bounds)
                actual=fn(x,left=x,right=x,interval=bounds,**params)
                self.assertAlmostEqual(actual,expected,delta=4e-12)
                np.testing.assert_array_equal(x,[.1,.3,.6])

    def test_empty_support_endpoints_and_identities(self):
        for _,_,api in self.builds:
            for family,params in FAMILIES:
                for method in ('loglikelihood','nnlf'):
                    fn=getattr(api.censored,family+'_'+method)
                    self.assertEqual(fn(**params),0.)
                    self.assertEqual(fn([],interval=[],**params),0.)
                    self.assertEqual(fn([-1.],**params),-np.inf if method=='loglikelihood' else np.inf)
                a=getattr(api.censored,family+'_loglikelihood')(left=[.2],right=[.8],**params)
                b=getattr(api.censored,family+'_loglikelihood')(interval=[[-np.inf,.2],[.8,np.inf]],**params)
                self.assertEqual(a,b)
                a=getattr(api.censored,family+'_loglikelihood')([.2,.4],**params)
                b=getattr(api.censored,family+'_loglikelihood')(interval=[[.2,.2],[.4,.4]],**params)
                self.assertEqual(a,b)
            kwargs=dict(left=[.1],right=[3.],interval=[[.4,.6]],loc=-.2,scale=2.3)
            self.assertEqual(api.censored.gamma_loglikelihood([.2,1.],a=1.,**kwargs),
                             api.censored.exponential_loglikelihood([.2,1.],**kwargs))
            x=np.array([-.3,.2,.7]);left=np.array([-.5]);right=np.array([1.])
            bounds=np.array([[-.1,.4]])
            n=api.censored.normal_loglikelihood(x,left=left,right=right,interval=bounds,loc=.3,scale=.8)
            ln=api.censored.lognormal_loglikelihood(np.exp(x),s=.8,left=np.exp(left),right=np.exp(right),
                                                   interval=np.exp(bounds),scale=np.exp(.3))
            self.assertAlmostEqual(ln,n-float(np.sum(x)),delta=3e-14)
            a=api.censored.beta_loglikelihood([.2,.6],a=.4,b=.6,left=[.1],right=[.8],interval=[[.25,.4]])
            b=api.censored.beta_loglikelihood([.8,.4],a=.6,b=.4,left=[.2],right=[.9],interval=[[.6,.75]])
            self.assertAlmostEqual(a,b,delta=3e-13)
            self.assertEqual(api.censored.lognormal_loglikelihood([0.],s=1.),-np.inf)
            self.assertEqual(api.censored.gamma_loglikelihood([0.],a=.4),np.inf)
            self.assertEqual(api.censored.gamma_loglikelihood([0.],a=2.),-np.inf)
            self.assertEqual(api.censored.beta_loglikelihood([0.],a=.4,b=2.),np.inf)
            self.assertEqual(api.censored.beta_loglikelihood([1.],a=.4,b=2.),-np.inf)
            self.assertTrue(np.isnan(api.censored.beta_loglikelihood([0.,1.],a=.4,b=2.)))
            self.assertEqual(api.censored.beta_loglikelihood(interval=[[-2.,2.]],a=.4,b=.6),0.)
            # This valid parameter regime retains the native unsupported-expansion NaN.
            low=api.censored_bridge.beta_loglikelihood([], [1e-300], [], [], [], a=.4,b=1e308,scale=1e24)
            self.assertEqual(low.status,0);self.assertTrue(np.isnan(low.outputs['value']))

    def test_saved_narrow_tail_references(self):
        # Existing mpmath references from test/test_censored.f90, binary64 endpoints.
        cases=(('lognormal',dict(s=.8,scale=1e100),[1e100,np.nextafter(1e100,np.inf)],-36.87309372451035,3e-13),
               ('lognormal',dict(s=1.),[np.exp(40.),np.exp(40.01)],-805.7174659453684,4e-11),
               ('gamma',dict(a=.4,scale=1e300),[0.,1e-300],-552.5008094043986,3e-12),
               ('gamma',dict(a=.4),[1.,np.nextafter(1.,np.inf)],-37.84033120681894,3e-13),
               ('gamma',dict(a=3.5),[1000.,1000.+1e-10],-1006.956992261547,5e-11),
               ('beta',dict(a=1000.,b=1000.),[.5,np.nextafter(.5,np.inf)],-33.16226569254558,4e-12))
        for _,_,api in self.builds:
            for family,params,bounds,expected,tolerance in cases:
                actual=getattr(api.censored,family+'_loglikelihood')(interval=[bounds],**params)
                self.assertAlmostEqual(actual,expected,delta=tolerance)
            self.assertAlmostEqual(api.censored.beta_loglikelihood(right=[1e-300],a=1e-20,b=.6,scale=1e300),
                                   -38.82004266769975,delta=4e-13)

    def test_invalid_parameters_data_and_shapes(self):
        for _,_,api in self.builds:
            for family,params in FAMILIES:
                for method in ('loglikelihood','nnlf'):
                    fn=getattr(api.censored,family+'_'+method)
                    for key in (*params,'scale'):
                        for bad in (0.,-1.,np.nan,np.inf,-np.inf):
                            with self.assertRaises(ValueError):fn(**dict(params,**{key:bad}))
                    for bad in (np.nan,np.inf,-np.inf):
                        with self.assertRaises(ValueError):fn(**params,loc=bad)
                    for key in ('uncensored','left','right'):
                        for bad in ([np.nan],[np.inf],[-np.inf],[[1.]],1.,[1j],['a']):
                            with self.assertRaises((ValueError,TypeError)):fn(**params,**{key:bad})
                    for bad in ([[2.,1.]],[[np.nan,2.]],[[1.,np.nan]],
                                [[-np.inf,np.inf]],[[np.inf,np.inf]],[[-np.inf,-np.inf]],
                                [1.,2.],[[1.,2.,3.]],np.empty((0,1))):
                        with self.assertRaises((ValueError,TypeError)):fn(**params,interval=bad)
                    if family=='beta':
                        with self.assertRaises(ValueError):fn(**params,loc=1e308,scale=1e308)
                low=getattr(api.censored_bridge,family+'_loglikelihood')
                with self.assertRaises(ValueError):low([],[],[],[.2],[.4,.6],**params)
                for bad in (dict(params,scale=0.),dict(params,loc=np.inf),dict(params,**{next(iter(params)):0.})):
                    result=low([],[],[],[],[],**bad)
                    self.assertEqual(result.status,2);self.assertTrue(np.isnan(result.outputs['value']))
                result=low([np.nan],[],[],[],[],**params)
                self.assertEqual(result.status,2);self.assertTrue(np.isnan(result.outputs['value']))
            result=api.censored_bridge.beta_loglikelihood([],[],[],[],[],a=.4,b=.6,loc=1e308,scale=1e308)
            self.assertEqual(result.status,2);self.assertTrue(np.isnan(result.outputs['value']))

    def test_raw_c_status_and_output_preservation(self):
        for _,_,api in self.builds:
            for family,shape_params in FAMILIES:
                for method in ('loglikelihood','nnlf'):
                    fn=getattr(api.censored_bridge._library,f'scifort_{family}_censored_{method}_f64')
                    parameters=[*shape_params.values(),0.,1.]
                    fn.argtypes=[C.c_int]*4+[C.c_void_p]*5+[C.c_double]*len(parameters)+[C.POINTER(C.c_double),C.POINTER(C.c_int)]
                    fn.restype=None
                    arrays=[np.array([v]) for v in (.3,.1,.8,.2,.4)]
                    ptrs=[a.ctypes.data for a in arrays]
                    value=C.c_double(77.);status=C.c_int()
                    for index in range(4):
                        dims=[1]*4;dims[index]=-1
                        fn(*dims,*ptrs,*parameters,C.byref(value),C.byref(status))
                        self.assertEqual(status.value,1);self.assertEqual(value.value,77.)
                    for index in range(len(shape_params)):
                        for invalid in (0.,-1.,np.nan,np.inf):
                            bad=parameters.copy();bad[index]=invalid
                            fn(1,1,1,1,*ptrs,*bad,C.byref(value),C.byref(status))
                            self.assertEqual(status.value,2);self.assertTrue(np.isnan(value.value))
                    arrays[3][0]=.5;arrays[4][0]=.2
                    fn(1,1,1,1,*ptrs,*parameters,C.byref(value),C.byref(status))
                    self.assertEqual(status.value,2);self.assertTrue(np.isnan(value.value))

    def test_c_cpp_new_symbols(self):
        for path,meta,_ in self.builds:
            with tempfile.TemporaryDirectory() as temporary:
                folder=Path(temporary);env=os.environ.copy()
                env['PATH']=str(path.parent)+os.pathsep+env.get('PATH','')
                client=folder/'client.c'
                client.write_text("""#include "scifort.h"
#include <math.h>
int main(void) {
 double x=1,l=.1,r=.8,a=.2,b=.4,value,nll; int status;
 scifort_lognormal_censored_loglikelihood_f64(1,0,0,0,&x,&l,&r,&a,&b,1,0,1,&value,&status);
 if(status || fabs(value+.9189385332046727)>2e-15) return 1;
 scifort_lognormal_censored_nnlf_f64(1,0,0,0,&x,&l,&r,&a,&b,1,0,1,&nll,&status);
 if(status || fabs(value+nll)>2e-15) return 2;
 scifort_gamma_censored_loglikelihood_f64(1,0,0,0,&x,&l,&r,&a,&b,1,0,1,&value,&status);
 if(status || value!=-1) return 3;
 scifort_gamma_censored_nnlf_f64(1,0,0,0,&x,&l,&r,&a,&b,1,0,1,&nll,&status);
 if(status || nll!=1) return 4;
 x=.5;
 scifort_beta_censored_loglikelihood_f64(1,0,0,0,&x,&l,&r,&a,&b,1,1,0,1,&value,&status);
 if(status || fabs(value)>2e-15) return 5;
 scifort_beta_censored_nnlf_f64(1,0,0,0,&x,&l,&r,&a,&b,1,1,0,1,&nll,&status);
 return status || fabs(nll)>2e-15;
}
""",encoding='ascii')
                for compiler,standard in (('gcc','c11'),('g++','c++17')):
                    if not shutil.which(compiler):continue
                    exe=folder/(compiler+'.exe')
                    run([compiler,'-std='+standard,'-Wall','-Wextra','-Werror','-I',ROOT/'include',client,path,'-o',exe])
                    run([exe],env=env)


if __name__=='__main__':
    unittest.main()
