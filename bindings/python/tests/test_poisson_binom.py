# SPDX-License-Identifier: MIT
"""Poisson-binomial scalar C/Python interfaces and conversion-range regressions."""
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
from scipy.stats import poisson_binom, binom
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from build_library import build
from scifort import Scifort
from test_multivariate_t import run
ROOT=Path(__file__).resolve().parents[3]
METHODS=('pmf','logpmf','cdf','sf','logcdf','logsf','ppf','isf')


class PoissonBinomTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran'):raise unittest.SkipTest('gfortran required')
        cls.builds=[]
        for checked in (True,False):
            path,metadata=build(checked=checked)
            cls.builds.append((path,metadata,Scifort(path,dll_directories=metadata['dll_directories'])))

    def test_all_methods_scipy_and_binomial_references(self):
        for _,_,api in self.builds:
            group=api.poisson_binom
            for p in ([.2,.6,.8],[.3],[0.,1.,.2,.7],[0.,0.,1.,1.]):
                for loc in (0.,2.5):
                    for method in METHODS:
                        values=(.01,.25,.6,.9,.99) if method in ('ppf','isf') else (-1.,0.,.2,1.,2.,3.,5.)
                        for x in values:
                            arg=x if method in ('ppf','isf') else x+loc
                            actual=getattr(group,method)(arg,p,loc)
                            expected=getattr(poisson_binom,method)(arg,np.array(p),loc=loc)
                            np.testing.assert_allclose(actual,expected,rtol=8e-13,atol=4e-14)
                            self.assertIsInstance(actual,float)
            for method in METHODS:
                arg=.8 if method in ('ppf','isf') else 3.
                np.testing.assert_allclose(getattr(group,method)(arg,[.3]*8),getattr(binom,method)(arg,8,.3),rtol=5e-13,atol=2e-14)

    def test_independent_enumeration_and_quantile_inequalities(self):
        p=[.1,.35,.65,.9]
        mass=np.zeros(5)
        for outcomes in itertools.product((0,1),repeat=4):
            mass[sum(outcomes)]+=math.prod(q if bit else 1-q for q,bit in zip(p,outcomes))
        for _,_,api in self.builds:
            group=api.poisson_binom
            for k in range(5):
                np.testing.assert_allclose(group.pmf(k,p),mass[k],rtol=3e-14,atol=1e-16)
                np.testing.assert_allclose(group.cdf(k,p),sum(mass[:k+1]),rtol=3e-14,atol=1e-16)
                np.testing.assert_allclose(group.sf(k,p),sum(mass[k+1:]),rtol=3e-14,atol=1e-16)
            for q in (.01,.2,.5,.8,.99):
                k=group.ppf(q,p)
                self.assertGreaterEqual(group.cdf(k,p),q)
                self.assertLess(group.cdf(k-1,p),q)
                k=group.isf(q,p)
                self.assertLessEqual(group.sf(k,p),q)
                self.assertGreater(group.sf(k-1,p),q)

    def test_endpoints_nonfinite_and_large_observations(self):
        for _,_,api in self.builds:
            group=api.poisson_binom;p=[.2,.6,.8]
            for k in (-1e100,-np.inf):
                self.assertEqual(group.cdf(k,p),0.);self.assertEqual(group.sf(k,p),1.)
                self.assertEqual(group.logcdf(k,p),-np.inf);self.assertEqual(group.logsf(k,p),0.)
            for k in (1e100,2**31+5.,np.inf):
                self.assertEqual(group.cdf(k,p),1.);self.assertEqual(group.sf(k,p),0.)
                self.assertEqual(group.logcdf(k,p),0.);self.assertEqual(group.logsf(k,p),-np.inf)
            self.assertEqual(group.cdf(-.2,p),0.)
            self.assertEqual(group.cdf(1.9,p),group.cdf(1,p))
            self.assertEqual(group.pmf(1.5,p),0.);self.assertEqual(group.logpmf(1.5,p),-np.inf)
            for method in METHODS:
                arg=.5 if method in ('ppf','isf') else 1.
                self.assertTrue(math.isnan(getattr(group,method)(np.nan,p)))
                for bad in ([np.nan,.5],[np.inf,.5],[-.1,.5],[1.1,.5]):
                    self.assertTrue(math.isnan(getattr(group,method)(arg,bad)))
                self.assertEqual(getattr(group,method)(arg,p),getattr(group,method)(arg,p,0.))
            for q in (-.1,1.1,np.inf,-np.inf):
                self.assertTrue(math.isnan(group.ppf(q,p)));self.assertTrue(math.isnan(group.isf(q,p)))
            for loc in (0.,2.5):
                self.assertEqual(group.ppf(0.,p,loc),loc-1)
                self.assertEqual(group.ppf(1.,p,loc),loc+3)
                self.assertEqual(group.isf(1.,p,loc),loc-1)
                self.assertEqual(group.isf(0.,p,loc),loc+3)

    def test_log_tail_underflow(self):
        for _,_,api in self.builds:
            group=api.poisson_binom
            p=[1e-5]*80
            self.assertEqual(group.sf(79,p),0.)
            np.testing.assert_allclose(group.logsf(79,p),80*math.log(1e-5),rtol=3e-14,atol=1e-12)
            p=[1.-1e-5]*80
            self.assertEqual(group.cdf(0,p),0.)
            np.testing.assert_allclose(group.logcdf(0,p),80*math.log1p(-p[0]),rtol=3e-14,atol=1e-12)

    def test_shapes_layout_and_raw_c_preservation(self):
        for _,_,api in self.builds:
            group=api.poisson_binom
            for method in METHODS:
                fn=getattr(group,method);arg=.5 if method in ('ppf','isf') else 1.
                for bad in ([],[[.2,.8]],np.ones((2,2,1))):
                    with self.assertRaises(ValueError):fn(arg,bad)
                with self.assertRaises(TypeError):fn(arg,[.2+1j,.8])
                with self.assertRaises((ValueError,TypeError)):fn([arg,arg],[.2,.8])
                p=np.array([.2,.8]);readonly=p.copy();readonly.flags.writeable=False
                for variant in (p.astype('>f8'),np.array([.2,9.,.8])[::2],readonly):
                    self.assertEqual(fn(arg,variant),fn(arg,p))
                raw=group._functions[method];value=C.c_double(99.);status=C.c_int(-1)
                for n in (0,-1):
                    raw(n,arg,p.ctypes.data_as(C.c_void_p),0.,C.byref(value),C.byref(status))
                    self.assertEqual(status.value,1);self.assertEqual(value.value,99.)
            huge=np.broadcast_to(np.ones(1),(np.iinfo(np.intc).max+1,))
            with self.assertRaises(ValueError):group.pmf(1.,huge)

    def test_native_c_clients_and_trapped_conversion_regressions(self):
        for path,metadata,_ in self.builds:
            with tempfile.TemporaryDirectory() as temp:
                folder=Path(temp);env=os.environ.copy();env['PATH']=str(path.parent)+os.pathsep+env.get('PATH','')
                flags=['-O0','-g','-fcheck=all','-fbacktrace'] if '-fcheck=all' in metadata['flags'] else ['-O3']
                source=ROOT/'test/test_poisson_binom_c_api.f90';exe=folder/'native.exe'
                run([metadata['compiler_path'],*flags,'-I',path.parent,source,path,'-o',exe]);run([exe],env=env)
                exe=folder/'trapped.exe'
                run([metadata['compiler_path'],*flags,'-ftrapv','-ffpe-trap=invalid','-I',path.parent,'-J',folder,
                     ROOT/'src/stats/poisson_binom.f90',source,path,'-o',exe]);run([exe],env=env)
                if shutil.which('gcc'):
                    c=folder/'client.c'
                    c.write_text("""#include "scifort.h"
#include <math.h>
int main(void) {
 double p[2]={.2,.8},a; int status;
 scifort_poisson_binom_pmf_explicit_f64(2,1,p,0,&a,&status); if(status || fabs(a-.68)>1e-14) return 1;
 scifort_poisson_binom_logpmf_explicit_f64(2,1,p,0,&a,&status); if(status || fabs(a-log(.68))>1e-14) return 2;
 scifort_poisson_binom_cdf_explicit_f64(2,1,p,0,&a,&status); if(status || fabs(a-.84)>1e-14) return 3;
 scifort_poisson_binom_sf_explicit_f64(2,1,p,0,&a,&status); if(status || fabs(a-.16)>1e-14) return 4;
 scifort_poisson_binom_logcdf_explicit_f64(2,1,p,0,&a,&status); if(status || fabs(a-log(.84))>1e-14) return 5;
 scifort_poisson_binom_logsf_explicit_f64(2,1,p,0,&a,&status); if(status || fabs(a-log(.16))>1e-14) return 6;
 scifort_poisson_binom_ppf_explicit_f64(2,.6,p,0,&a,&status); if(status || a!=1) return 7;
 scifort_poisson_binom_isf_explicit_f64(2,.4,p,0,&a,&status); return status || a!=1;
}
""",encoding='ascii')
                    exe=folder/'c.exe'
                    run(['gcc','-std=c11','-Wall','-Wextra','-Werror','-I',ROOT/'include',c,path,'-o',exe]);run([exe],env=env)
                    if shutil.which('g++'):run(['g++','-std=c++17','-Wall','-Wextra','-Werror','-I',ROOT/'include','-fsyntax-only',c])
