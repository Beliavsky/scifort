# SPDX-License-Identifier: MIT
"""Reviewed ECDF ownership, array contracts and native/SciPy agreement."""
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
import warnings
import numpy as np
from scipy.stats import ecdf, CensoredData
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from build_library import build
from scifort import Scifort
from test_multivariate_t import run
ROOT=Path(__file__).resolve().parents[3]

class ECDFTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran'):raise unittest.SkipTest('gfortran required')
        cls.builds=[]
        for checked in (True,False):
            path,meta=build(checked=checked)
            cls.builds.append((path,meta,Scifort(path,dll_directories=meta['dll_directories'])))

    def test_scipy_steps_ties_empty_and_endpoints(self):
        cases=([3.,1.,2.,2.],[]),([1.,2.,2.,4.],[2.,5.]),([],[1.,2.]),([],[])
        x=np.array([-np.inf,0.,1.,1.5,2.,4.,6.,np.inf,np.nan])
        for _,_,api in self.builds:
            for events,right in cases:
                reference=ecdf(CensoredData(uncensored=events,right=right))
                with api.ecdf(events,right=right) as result:
                    for name in ('cdf','sf'):
                        f=getattr(result,name);ref=getattr(reference,name)
                        np.testing.assert_array_equal(f.quantiles,ref.quantiles)
                        np.testing.assert_allclose(f.probabilities,ref.probabilities,atol=2e-16,rtol=2e-15)
                        np.testing.assert_allclose(f.evaluate(x),ref.evaluate(x),atol=2e-16,rtol=2e-15)
                        self.assertEqual(f(1.5),float(ref.evaluate(1.5)))
                        self.assertEqual(f.evaluate([]).shape,(0,))
                        self.assertTrue(np.all(np.diff(f.quantiles)>0))
                    np.testing.assert_allclose(result.cdf(x)+result.sf(x),1.,atol=2e-16)

    def test_confidence_intervals(self):
        for _,_,api in self.builds:
            for events,right in (([1.,2.,2.,4.],[2.,5.]),([1.,2.,3.],[]),([],[1.,2.]),([],[])):
                reference=ecdf(CensoredData(uncensored=events,right=right))
                with api.ecdf(events,right=right) as result:
                    for name in ('cdf','sf'):
                        for method in ('linear','log-log'):
                            for level in (0.,.9,1.):
                                with warnings.catch_warnings():
                                    warnings.simplefilter('ignore')
                                    ref=getattr(reference,name).confidence_interval(level,method=method)
                                with getattr(result,name).confidence_interval(level,method) as ci:
                                    for bound in ('low','high'):
                                        actual=getattr(ci,bound);expected=getattr(ref,bound)
                                        np.testing.assert_array_equal(actual.quantiles,expected.quantiles)
                                        np.testing.assert_allclose(actual.probabilities,expected.probabilities,
                                                                   rtol=2e-13,atol=3e-15,equal_nan=True)
                                    with self.assertRaises(ValueError):ci.low.confidence_interval()
                    for bad in (-.1,1.1,np.nan,np.inf):
                        with self.assertRaises(ValueError):result.cdf.confidence_interval(bad)
                    with self.assertRaises(ValueError):result.cdf.confidence_interval(method='bad')
                    rejected=api.ecdf_bridge.interval(result.cdf._handle,method=9)
                    self.assertEqual(rejected.status,2)
                    self.assertTrue(all(value is None for value in rejected.outputs.values()))

    def test_buffers_and_shape_errors(self):
        for _,_,api in self.builds:
            with api.ecdf([3.,1.,2.,2.]) as result:
                f=result.cdf
                q=np.empty(3);p=np.empty(3)
                qout,pout=f.get_values(out={'quantiles':q,'probabilities':p})
                self.assertIs(qout,q);self.assertIs(pout,p)
                np.testing.assert_allclose(p,[.25,.75,1.])
                q[:]=-99;p[:]=-99
                np.testing.assert_array_equal(f.quantiles,[1.,2.,3.])
                x=np.arange(8.)[::2];out=np.empty(x.size)
                self.assertIs(f.evaluate(x,out=out),out)
                np.testing.assert_array_equal(out,[0.,.75,1.,1.])
                for bad in (np.empty(3),np.empty(4,dtype=np.float32),np.empty(8)[::2]):
                    with self.assertRaises((TypeError,ValueError)):f.evaluate(x,out=bad)
                out.flags.writeable=False
                with self.assertRaises((TypeError,ValueError)):f.evaluate(x,out=out)
                with self.assertRaises(ValueError):f.evaluate([[1.,2.]])
                with self.assertRaises(ValueError):f.get_values(out={'quantiles':q,'probabilities':q})
                bad=api.ecdf_bridge.values(f._handle,m=4)
                self.assertEqual(bad.status,1)
                self.assertTrue(all(np.isnan(a).all() for a in bad.outputs.values()))
            for sample,right in (([np.nan],[]),([np.inf],[]),([],[np.inf]),([[1.,2.]],[])):
                with self.assertRaises(ValueError):api.ecdf(sample,right=right)
            shared=np.array([1.,2.,3.])
            with api.ecdf(shared,right=shared) as result:self.assertEqual(result.cdf.quantiles.size,3)

    def test_independent_ownership_and_cleanup(self):
        for path,meta,api in self.builds:
            values=np.array([1.,2.,3.]);pair=api.ecdf(values);f=pair.cdf;sf=pair.sf
            ci=f.confidence_interval();low=ci.low
            values[:]=-100
            sf.close();sf.close()
            del pair;gc.collect()
            np.testing.assert_array_equal(f.quantiles,[1.,2.,3.])
            for obj in (f,ci):
                for operation in (copy.copy,copy.deepcopy,pickle.dumps):
                    with self.assertRaises(TypeError):operation(obj)
            f.close()
            self.assertTrue(np.isfinite(low.quantiles).all())
            for operation in (lambda:f.quantiles,lambda:f.evaluate([]),lambda:f.confidence_interval()):
                with self.assertRaises(ValueError):operation()
            other=Scifort(path,dll_directories=meta['dll_directories'])
            with self.assertRaises(ValueError):other.ecdf_bridge.info(low._handle)
            ci.high.close();del ci;gc.collect()
            self.assertEqual(low.quantiles.size,3)
            low.close()

    def test_raw_c_guards(self):
        for _,_,api in self.builds:
            # The group owns the exact CDLL; obtain functions from a live owner's library.
            with api.ecdf([1.,2.]) as pair:
                lib=pair.cdf._handle._library
                fn=lib.scifort_ecdf_evaluate_f64
                fn.argtypes=[C.c_int,C.c_void_p,C.c_void_p,C.c_void_p,C.POINTER(C.c_int)];fn.restype=None
                x=np.array([0.,2.]);out=np.full(2,77.);status=C.c_int()
                fn(2,None,x.ctypes.data,out.ctypes.data,C.byref(status))
                self.assertEqual(status.value,5);self.assertTrue(np.isnan(out).all())
                out[:]=77.
                fn(-1,None,x.ctypes.data,out.ctypes.data,C.byref(status))
                self.assertEqual(status.value,1);np.testing.assert_array_equal(out,77.)
                create=lib.scifort_ecdf_create_f64
                create.argtypes=[C.c_int,C.c_int,C.c_void_p,C.c_void_p,C.POINTER(C.c_void_p),
                                 C.POINTER(C.c_void_p),C.POINTER(C.c_int)];create.restype=None
                a=C.c_void_p();b=C.c_void_p();x[0]=np.nan
                create(2,0,x.ctypes.data,out.ctypes.data,C.byref(a),C.byref(b),C.byref(status))
                self.assertEqual(status.value,2);self.assertIsNone(a.value);self.assertIsNone(b.value)

    def test_native_regression_and_c_cpp_clients(self):
        for path,meta,_ in self.builds:
            with tempfile.TemporaryDirectory() as temporary:
                folder=Path(temporary);env=os.environ.copy();env['PATH']=str(path.parent)+os.pathsep+env.get('PATH','')
                flags=['-O0','-g','-fcheck=all','-fbacktrace'] if '-fcheck=all' in meta['flags'] else ['-O3']
                program=folder/'regression.f90'
                original=(ROOT/'test/test_empirical_distribution.f90').read_text(encoding='utf-8')
                old='    use scifort_stats, only : confidence_interval, ecdf, ecdf_confidence_interval, &\n        ecdf_result, normal_cdf, normal_pdf, normal_ppf, order_statistic, &\n        order_statistic_distribution, quantile_test, quantile_test_result'
                new='    use scifort_contingency_meta, only : confidence_interval\n    use scifort_normal, only : normal_cdf, normal_pdf, normal_ppf\n    use scifort_empirical_distribution, only : ecdf, ecdf_confidence_interval, &\n        ecdf_result, order_statistic, order_statistic_distribution, quantile_test, quantile_test_result'
                self.assertIn(old,original)
                program.write_text(original.replace(old,new),encoding='utf-8')
                for source in (program,ROOT/'test/test_ecdf_c_api.f90'):
                    exe=folder/(source.stem+'.exe')
                    run([meta['compiler_path'],*flags,'-I',path.parent,source,path,'-o',exe])
                    run([exe],env=env)
                client=folder/'client.c'
                client.write_text("""#include "scifort.h"
#include <math.h>
#include <stddef.h>
int main(void) {
 double x[4]={1,2,2,4},right[2]={2,5},q[4],p[4],v[4];
 void *cdf=NULL,*sf=NULL,*lo=NULL,*hi=NULL; int status,n,variance,survival;
 scifort_ecdf_create_f64(4,2,x,right,&cdf,&sf,&status); if(status || !cdf || !sf) return 1;
 scifort_ecdf_info_f64(cdf,&n,&variance,&survival,&status);
 if(status || n!=4 || variance!=1 || survival!=0) return 2;
 scifort_ecdf_values_f64(n,cdf,q,p,&status); if(status || fabs(p[2]-.75)>1e-15) return 3;
 scifort_ecdf_evaluate_f64(4,cdf,q,v,&status); if(status || v[2]!=p[2]) return 4;
 scifort_ecdf_interval_f64(cdf,.9,1,&lo,&hi,&status); if(status || !lo || !hi) return 5;
 scifort_ecdf_destroy_f64(&cdf,&status); if(status || cdf) return 6;
 scifort_ecdf_evaluate_f64(4,lo,q,v,&status); if(status) return 7;
 scifort_ecdf_destroy_f64(&lo,&status); scifort_ecdf_destroy_f64(&hi,&status);
 scifort_ecdf_destroy_f64(&sf,&status); scifort_ecdf_destroy_f64(&sf,&status);
 return status || sf!=NULL;
}
""",encoding='ascii')
                for compiler,standard in (('gcc','c11'),('g++','c++17')):
                    if not shutil.which(compiler):continue
                    exe=folder/(compiler+'.exe')
                    run([compiler,'-std='+standard,'-Wall','-Wextra','-Werror','-I',ROOT/'include',client,path,'-o',exe])
                    run([exe],env=env)
