# SPDX-License-Identifier: MIT
"""Location-scale censored adapters, numerical identities and flat C contracts."""
import ctypes as C
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
import numpy as np
from scipy.stats import logistic, laplace, cauchy, rayleigh
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from build_library import build
from scifort import Scifort
import test_censored as helpers
from test_multivariate_t import run
ROOT = Path(__file__).resolve().parents[3]
FAMILIES = dict(logistic=logistic, laplace=laplace, cauchy=cauchy, rayleigh=rayleigh)


class CensoredLocationScaleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran'):
            raise unittest.SkipTest('gfortran required')
        cls.builds=[]
        for checked in (True,False):
            path,meta=build(checked=checked)
            cls.builds.append((path,meta,Scifort(path,dll_directories=meta['dll_directories'])))

    def test_mixed_references_layout_and_nnlf(self):
        cases=(([.5,1.2,2.5],[.3],[3.5],[[1.1,1.4],[4.,4.4]]),
               ([],[.05,.1],[.8,.9,1.],[[.4,.4],[-np.inf,.5],[.8,np.inf]]),
               ([.2,.4],[],[],[]), ([],[],[],[[.2,.5]]))
        for _,_,api in self.builds:
            for family,distribution in FAMILIES.items():
                dist=distribution(loc=-.4,scale=1.3)
                fn=getattr(api.censored,family+'_loglikelihood')
                for exact,left,right,interval in cases:
                    expected=helpers.CensoredLikelihoodTests.reference(dist,exact,left,right,interval)
                    actual=fn(exact,left=left,right=right,interval=interval,loc=-.4,scale=1.3)
                    np.testing.assert_allclose(actual,expected,rtol=2e-13,atol=2e-12)
                    nll=getattr(api.censored,family+'_nnlf')(
                        exact,left=left,right=right,interval=interval,loc=-.4,scale=1.3)
                    self.assertEqual(nll,-actual)
                x=np.array([.1,99.,.3,99.,.6])[::2];x.flags.writeable=False
                bounds=np.array([[.2,.25],[99.,99.],[.7,.75]])[::2];bounds.flags.writeable=False
                expected=helpers.CensoredLikelihoodTests.reference(dist,x,x,x,bounds)
                actual=fn(x,left=x,right=x,interval=bounds,loc=-.4,scale=1.3)
                self.assertAlmostEqual(actual,expected,delta=4e-12)
                np.testing.assert_array_equal(x,[.1,.3,.6])

    def test_empty_endpoints_reflection_and_rayleigh_identity(self):
        for _,_,api in self.builds:
            for family in FAMILIES:
                fn=getattr(api.censored,family+'_loglikelihood')
                self.assertEqual(fn(),0.)
                self.assertEqual(getattr(api.censored,family+'_nnlf')(),0.)
                self.assertEqual(fn([],interval=[]),0.)
                a=fn(left=[.2],right=[.8])
                b=fn(interval=[[-np.inf,.2],[.8,np.inf]])
                self.assertEqual(a,b)
                self.assertEqual(fn([.2,.4]),fn(interval=[[.2,.2],[.4,.4]]))
                if family!='rayleigh':
                    a=fn([-.5,.2,1.],left=[-.8],right=[1.2,2.],interval=[[-.2,.4]],loc=.3,scale=1.7)
                    b=fn([.5,-.2,-1.],left=[-1.2,-2.],right=[.8],interval=[[-.4,.2]],loc=-.3,scale=1.7)
                    self.assertAlmostEqual(a,b,delta=3e-13)
            kwargs=dict(left=[.3],right=[3.],interval=[[.4,.6]],loc=-.2)
            a=api.censored.rayleigh_loglikelihood([.2,1.],scale=1.3,**kwargs)
            b=api.censored.weibull_loglikelihood([.2,1.],c=2.,scale=np.sqrt(2.)*1.3,**kwargs)
            self.assertAlmostEqual(a,b,delta=4e-14)
            fn=api.censored.rayleigh_loglikelihood
            self.assertEqual(fn([-1.]),-np.inf)
            self.assertEqual(fn([0.]),-np.inf)
            self.assertEqual(fn(left=[0.]),-np.inf)
            self.assertEqual(fn(right=[-1.,0.]),0.)
            self.assertEqual(fn(interval=[[-1.,1.]]),fn(interval=[[0.,1.]]))
            self.assertEqual(api.censored.rayleigh_nnlf([-1.]),np.inf)

    def test_narrow_extreme_tail_and_underflow_references(self):
        # Formulas and saved references already used in the corresponding native tests.
        a=1000.;b=np.nextafter(a,np.inf);width=b-a
        tail_cases=(('logistic',-1029.8053287640778,3e-13),
                    ('laplace',-a-np.log(2.)+np.log(width),3e-13),
                    ('cauchy',-np.log(np.pi)-np.log1p(a*a)+np.log(width),3e-13))
        for _,_,api in self.builds:
            for family,expected,tolerance in tail_cases:
                fn=getattr(api.censored,family+'_loglikelihood')
                for bounds in ([a,b],[-b,-a]):
                    self.assertAlmostEqual(fn(interval=[bounds]),expected,delta=tolerance)
            ray_a=40.;ray_b=np.nextafter(ray_a,np.inf)
            expected=np.log(ray_a)-.5*ray_a*ray_a+np.log(ray_b-ray_a)
            self.assertAlmostEqual(api.censored.rayleigh_loglikelihood(
                interval=[[ray_a,ray_b]]),expected,delta=1e-12)
            span=np.log(1e-300)-np.log(1e300)
            for family,expected in (('logistic',span-np.log(4.)),('laplace',span-np.log(2.)),
                                    ('cauchy',span-np.log(np.pi))):
                actual=getattr(api.censored,family+'_loglikelihood')(
                    interval=[[1e-300,2e-300]],scale=1e300)
                self.assertAlmostEqual(actual,expected,delta=5e-13)
            actual=api.censored.rayleigh_loglikelihood(interval=[[1e-300,2e-300]],scale=1e300)
            self.assertAlmostEqual(actual,2*span+np.log(1.5),delta=1e-12)
            # Subnormal interval width must not become an exact observation.
            tiny=np.nextafter(0.,1.)
            actual=api.censored.logistic_loglikelihood(interval=[[0.,tiny]])
            self.assertAlmostEqual(actual,np.log(tiny)-np.log(4.),delta=2e-13)

    def test_invalid_data_parameters_and_shapes(self):
        for _,_,api in self.builds:
            for family in FAMILIES:
                for method in ('loglikelihood','nnlf'):
                    fn=getattr(api.censored,family+'_'+method)
                    for bad in (0.,-1.,np.nan,np.inf,-np.inf):
                        with self.assertRaises(ValueError):fn(scale=bad)
                    for bad in (np.nan,np.inf,-np.inf):
                        with self.assertRaises(ValueError):fn(loc=bad)
                    for key in ('uncensored','left','right'):
                        for bad in ([np.nan],[np.inf],[-np.inf],[[1.]],1.,[1j],['a']):
                            with self.assertRaises((ValueError,TypeError)):fn(**{key:bad})
                    for bad in ([[2.,1.]],[[np.nan,2.]],[[1.,np.nan]],
                                [[-np.inf,np.inf]],[[np.inf,np.inf]],[[-np.inf,-np.inf]],
                                [1.,2.],[[1.,2.,3.]],np.empty((0,1))):
                        with self.assertRaises((ValueError,TypeError)):fn(interval=bad)
                low=getattr(api.censored_bridge,family+'_loglikelihood')
                with self.assertRaises(ValueError):low([],[],[],[.2],[.4,.6])
                for bad in (dict(scale=0.),dict(loc=np.inf)):
                    result=low([],[],[],[],[],**bad)
                    self.assertEqual(result.status,2);self.assertTrue(np.isnan(result.outputs['value']))
                result=low([np.nan],[],[],[],[])
                self.assertEqual(result.status,2);self.assertTrue(np.isnan(result.outputs['value']))

    def test_raw_c_status_and_output_preservation(self):
        for _,_,api in self.builds:
            for family in FAMILIES:
                for method in ('loglikelihood','nnlf'):
                    fn=getattr(api.censored_bridge._library,f'scifort_{family}_censored_{method}_f64')
                    fn.argtypes=[C.c_int]*4+[C.c_void_p]*5+[C.c_double]*2+[C.POINTER(C.c_double),C.POINTER(C.c_int)]
                    fn.restype=None
                    arrays=[np.array([v]) for v in (.3,.1,.8,.2,.4)]
                    ptrs=[a.ctypes.data for a in arrays]
                    value=C.c_double(77.);status=C.c_int()
                    for index in range(4):
                        dims=[1]*4;dims[index]=-1
                        fn(*dims,*ptrs,0.,1.,C.byref(value),C.byref(status))
                        self.assertEqual(status.value,1);self.assertEqual(value.value,77.)
                    for loc,scale in ((np.inf,1.),(0.,0.),(0.,np.nan)):
                        fn(1,1,1,1,*ptrs,loc,scale,C.byref(value),C.byref(status))
                        self.assertEqual(status.value,2);self.assertTrue(np.isnan(value.value))
                    arrays[3][0]=.5;arrays[4][0]=.2
                    fn(1,1,1,1,*ptrs,0.,1.,C.byref(value),C.byref(status))
                    self.assertEqual(status.value,2);self.assertTrue(np.isnan(value.value))

    def test_c_cpp_new_symbols(self):
        for path,meta,_ in self.builds:
            with tempfile.TemporaryDirectory() as temporary:
                folder=Path(temporary);env=os.environ.copy()
                env['PATH']=str(path.parent)+os.pathsep+env.get('PATH','')
                calls=[]
                for family,expected in (('logistic','-log(4.)'),('laplace','-log(2.)'),
                                        ('cauchy','-log(acos(-1.))'),('rayleigh','-.5')):
                    calls.append(f'x={1 if family=="rayleigh" else 0};')
                    calls.append(f'scifort_{family}_censored_loglikelihood_f64(1,0,0,0,&x,&l,&r,&a,&b,0,1,&value,&status);')
                    calls.append(f'if(status || fabs(value-({expected}))>3e-15) return 1;')
                    calls.append(f'scifort_{family}_censored_nnlf_f64(1,0,0,0,&x,&l,&r,&a,&b,0,1,&nll,&status);')
                    calls.append('if(status || fabs(value+nll)>3e-15) return 2;')
                client=folder/'client.c'
                client.write_text('#include "scifort.h"\n#include <math.h>\nint main(void) {\n'
                                  'double x=0,l=.1,r=.8,a=.2,b=.4,value,nll; int status;\n'
                                  +'\n'.join(calls)+'\nreturn 0;\n}\n',encoding='ascii')
                for compiler,standard in (('gcc','c11'),('g++','c++17')):
                    if not shutil.which(compiler):continue
                    exe=folder/(compiler+'.exe')
                    run([compiler,'-std='+standard,'-Wall','-Wextra','-Werror','-I',ROOT/'include',client,path,'-o',exe])
                    run([exe],env=env)


if __name__=='__main__':
    unittest.main()
