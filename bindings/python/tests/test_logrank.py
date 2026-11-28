# SPDX-License-Identifier: MIT
"""Flat survival interfaces, native agreement and lightweight type dependencies."""
import ctypes as C
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
import warnings
import numpy as np
from scipy.stats import CensoredData, logrank
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from build_library import build
from scifort import Scifort, LogrankResult
from test_multivariate_t import run
ROOT=Path(__file__).resolve().parents[3]

class LogrankTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran'):raise unittest.SkipTest('gfortran required')
        cls.builds=[]
        for checked in (True,False):
            path,meta=build(checked=checked)
            cls.builds.append((path,meta,Scifort(path,dll_directories=meta['dll_directories'])))

    def test_scipy_references_independent_lengths_and_ties(self):
        cases=(([1.,2.],[],[2.,3.],[]),([1.,2.,2.],[2.,5.],[2.,3.,4.,6.],[4.]),
               ([1.],[7.,8.,9.],[2.,4.,5.,6.],[]),([],[1.,2.,3.],[1.,2.],[3.]),
               ([1.,2.,3.],[3.],[],[1.,2.,4.,5.]))
        for _,_,api in self.builds:
            for x,xr,y,yr in cases:
                for alternative in ('two-sided','less','greater'):
                    ref=logrank(CensoredData(uncensored=x,right=xr),CensoredData(uncensored=y,right=yr),
                                alternative=alternative)
                    actual=api.logrank(x,y,x_right=xr,y_right=yr,alternative=alternative)
                    self.assertIsInstance(actual,LogrankResult)
                    self.assertAlmostEqual(actual.statistic,ref.statistic,delta=2e-14)
                    self.assertAlmostEqual(actual.pvalue,ref.pvalue,delta=3e-15)
                    self.assertGreaterEqual(actual.pvalue,0.);self.assertLessEqual(actual.pvalue,1.)

    def test_group_swap_identities_and_tails(self):
        for _,_,api in self.builds:
            x=np.array([1.,2.,2.,5.]);xr=np.array([3.,6.]);y=np.array([2.,4.,7.]);yr=np.array([4.])
            a=api.logrank(x,y,x_right=xr,y_right=yr)
            b=api.logrank(y,x,x_right=yr,y_right=xr)
            self.assertAlmostEqual(a.statistic,-b.statistic,delta=3e-15)
            self.assertAlmostEqual(a.pvalue,b.pvalue,delta=3e-15)
            less=api.logrank(x,y,x_right=xr,y_right=yr,alternative='less')
            greater=api.logrank(y,x,x_right=yr,y_right=xr,alternative='greater')
            self.assertAlmostEqual(less.pvalue,greater.pvalue,delta=3e-15)
            x=np.arange(1.,81.);y=np.arange(81.,161.)
            a=api.logrank(x,y,alternative='greater');ref=logrank(x,y,alternative='greater')
            self.assertLess(a.pvalue,1e-20)
            np.testing.assert_allclose(a.statistic,ref.statistic,rtol=5e-14,atol=0.)
            np.testing.assert_allclose(a.pvalue,ref.pvalue,rtol=2e-11,atol=0.)

    def test_empty_all_censored_and_zero_variance(self):
        cases=(([],[],[],[]),([],[1.,2.],[],[3.]),([1.],[],[1.],[]),([1.],[],[],[]))
        for _,_,api in self.builds:
            for x,xr,y,yr in cases:
                for alternative in ('two-sided','less','greater'):
                    with warnings.catch_warnings():
                        warnings.simplefilter('ignore')
                        ref=logrank(CensoredData(uncensored=x,right=xr),CensoredData(uncensored=y,right=yr),
                                    alternative=alternative)
                    actual=api.logrank(x,y,x_right=xr,y_right=yr,alternative=alternative)
                    np.testing.assert_allclose(actual,(ref.statistic,ref.pvalue),equal_nan=True,rtol=0.,atol=0.)
                    self.assertTrue(np.isnan(actual.statistic));self.assertTrue(np.isnan(actual.pvalue))

    def test_input_validation_aliases_layout_and_low_level(self):
        for _,_,api in self.builds:
            x=np.arange(1.,9.)[::2];y=x.copy();saved=x.copy();x.flags.writeable=False
            same=api.logrank(x,x)
            self.assertEqual(same.statistic,0.);self.assertEqual(same.pvalue,1.)
            np.testing.assert_array_equal(x,saved)
            expected=api.logrank(x,y)
            low=api.survival.logrank(x,[],y,[])
            self.assertEqual(low.status,0)
            self.assertEqual(expected,LogrankResult(low.outputs['statistic'],low.outputs['pvalue']))
            for name in ('x','y','x_right','y_right'):
                for bad in ([np.nan],[np.inf],[-np.inf],[[1.,2.]],1.):
                    kwargs={'x':[1.],'y':[2.],'x_right':[],'y_right':[]};kwargs[name]=bad
                    with self.assertRaises((TypeError,ValueError)):api.logrank(**kwargs)
            for bad in ('bad',None,0,True):
                with self.assertRaises(ValueError):api.logrank([1.],[2.],alternative=bad)
            for bad in (-1,3):
                result=api.survival.logrank([1.],[],[2.],[],alternative=bad)
                self.assertEqual(result.status,2)
                self.assertTrue(all(np.isnan(v) for v in result.outputs.values()))
            with self.assertRaises(TypeError):api.survival.logrank([1.],[],[2.],[],alternative=True)

    def test_raw_c_invalid_extents_and_parameters(self):
        for _,_,api in self.builds:
            fn=api.survival._library.scifort_logrank_f64
            fn.argtypes=[C.c_int]*4+[C.c_void_p]*4+[C.c_int]+[C.POINTER(C.c_double)]*2+[C.POINTER(C.c_int)]
            fn.restype=None
            x=np.array([1.]);y=np.array([2.]);empty=np.empty(0)
            z=C.c_double(77.);p=C.c_double(88.);status=C.c_int()
            args=(x.ctypes.data,empty.ctypes.data,y.ctypes.data,empty.ctypes.data)
            fn(-1,0,1,0,*args,0,C.byref(z),C.byref(p),C.byref(status))
            self.assertEqual(status.value,1);self.assertEqual(z.value,77.);self.assertEqual(p.value,88.)
            for alternative in (-1,3):
                fn(1,0,1,0,*args,alternative,C.byref(z),C.byref(p),C.byref(status))
                self.assertEqual(status.value,2);self.assertTrue(np.isnan(z.value));self.assertTrue(np.isnan(p.value))
            x[0]=np.nan
            fn(1,0,1,0,*args,0,C.byref(z),C.byref(p),C.byref(status))
            self.assertEqual(status.value,2);self.assertTrue(np.isnan(z.value));self.assertTrue(np.isnan(p.value))

    def test_native_regression_and_c_cpp_clients(self):
        for path,meta,_ in self.builds:
            with tempfile.TemporaryDirectory() as temporary:
                folder=Path(temporary);env=os.environ.copy();env['PATH']=str(path.parent)+os.pathsep+env.get('PATH','')
                flags=['-O0','-g','-fcheck=all','-fbacktrace'] if '-fcheck=all' in meta['flags'] else ['-O3']
                program=folder/'native_logrank.f90'
                text=(ROOT/'test/test_logrank.f90').read_text(encoding='utf-8')
                text=text.replace('use scifort_stats, only :','use scifort_empirical_distribution, only :')
                program.write_text(text,encoding='utf-8')
                for source in (program,ROOT/'test/test_logrank_c_api.f90'):
                    exe=folder/(source.stem+'.exe')
                    run([meta['compiler_path'],*flags,'-I',path.parent,source,path,'-o',exe])
                    run([exe],env=env)
                client=folder/'client.c'
                client.write_text("""#include "scifort.h"
#include <math.h>
int main(void) {
 double x=1,y=2,empty=0,z,p; int status;
 scifort_logrank_f64(1,0,1,0,&x,&empty,&y,&empty,0,&z,&p,&status);
 if(status || fabs(z-1)>1e-15 || fabs(p-.31731050786291415)>2e-15) return 1;
 scifort_logrank_f64(0,0,0,0,&x,&empty,&y,&empty,0,&z,&p,&status);
 if(status || !isnan(z) || !isnan(p)) return 2;
 scifort_logrank_f64(1,0,1,0,&x,&empty,&y,&empty,3,&z,&p,&status);
 return status!=2 || !isnan(z) || !isnan(p);
}
""",encoding='ascii')
                for compiler,standard in (('gcc','c11'),('g++','c++17')):
                    if not shutil.which(compiler):continue
                    exe=folder/(compiler+'.exe')
                    run([compiler,'-std='+standard,'-Wall','-Wextra','-Werror','-I',ROOT/'include',client,path,'-o',exe])
                    run([exe],env=env)

    def test_lightweight_dependency_plan(self):
        plan=json.loads((ROOT/'bindings/python/build_sources.json').read_text(encoding='utf-8'))
        sources=set()
        pending=['src/empirical_distribution.f90']
        while pending:
            source=pending.pop()
            if source in sources:continue
            sources.add(source)
            pending.extend(plan['dependencies'][source])
        self.assertIn('src/confidence_interval.f90',sources)
        self.assertNotIn('src/contingency_meta.f90',sources)
        self.assertNotIn('src/hypothesis_extended.f90',sources)
        self.assertFalse(any('/descriptive' in name for name in sources))
        self.assertEqual(plan['dependencies']['src/empirical_distribution.f90'].count('src/confidence_interval.f90'),1)
