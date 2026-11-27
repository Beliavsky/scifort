# SPDX-License-Identifier: MIT
"""Paired normal-inverse-gamma results, boundary behavior and explicit streams."""
import ctypes as C
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from build_library import build
from scifort import Scifort
from test_multivariate_t import run
ROOT=Path(__file__).resolve().parents[3]

class NormalInverseGammaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran'): raise unittest.SkipTest('gfortran required')
        cls.builds=[]
        for checked in (True,False):
            path,meta=build(checked=checked)
            cls.builds.append((path,meta,Scifort(path,dll_directories=meta['dll_directories'])))

    def test_density_defaults_references_and_endpoints(self):
        from scipy.stats import invgamma, norm
        for _,_,api in self.builds:
            family=api.normal_inverse_gamma
            for mu,lam,a,b in ((0.,1.,1.,1.),(-2.,3.,4.,2.)):
                for x,s2 in ((0.,1.),(2.,.2),(-4.,20.),(40.,.1)):
                    expected=invgamma.logpdf(s2,a,scale=b)+norm.logpdf(x,mu,np.sqrt(s2/lam))
                    actual=family.logpdf(x,s2,mu,lam,a,b)
                    self.assertAlmostEqual(actual,expected,delta=2e-12*max(1.,abs(expected)))
                    np.testing.assert_allclose(family.pdf(x,s2,mu,lam,a,b),np.exp(expected),rtol=3e-13,atol=0)
            self.assertEqual(family.pdf(0.,1.),family.pdf(0.,1.,0.,1.,1.,1.))
            for x,s2 in ((0.,0.),(0.,-1.),(np.inf,1.),(-np.inf,1.),(0.,np.inf)):
                self.assertEqual(family.pdf(x,s2),0.)
                self.assertEqual(family.logpdf(x,s2),-np.inf)
            for x,s2 in ((np.nan,1.),(0.,np.nan)):
                self.assertTrue(np.isnan(family.pdf(x,s2)))
            for key in ('mu','lmbda','a','b'):
                for value in ((np.nan,np.inf,-np.inf) if key=='mu' else (0.,-1.,np.nan,np.inf)):
                    self.assertTrue(np.isnan(family.pdf(0.,1.,**{key:value})))
            with self.assertRaises((TypeError,ValueError)): family.pdf([0.,1.],1.)

    def test_paired_moments_and_boundaries(self):
        for _,_,api in self.builds:
            f=api.normal_inverse_gamma
            for a in (.5,1.,1.5,2.,3.,np.nan):
                means=f.mean(mu=2.,lmbda=4.,a=a,b=3.)
                var=f.var(mu=2.,lmbda=4.,a=a,b=3.)
                self.assertEqual(means.status,0); self.assertEqual(var.status,0)
                for key,expected in dict(mean_x=2. if a>1 else np.nan,
                                         mean_s2=3./(a-1) if a>1 else np.nan).items():
                    np.testing.assert_equal(means.outputs[key],expected)
                np.testing.assert_equal(var.outputs['var_x'],3./((a-1)*4.) if a>1 else np.nan)
                np.testing.assert_equal(var.outputs['var_s2'],9./((a-1)**2*(a-2)) if a>2 else np.nan)
            self.assertTrue(all(np.isnan(v) for v in f.mean().outputs.values()))
            self.assertTrue(all(np.isnan(v) for v in f.var(lmbda=0.).outputs.values()))

    def test_single_batch_split_streams_and_empty_calls(self):
        for _,_,api in self.builds:
            f=api.normal_inverse_gamma; state=api.rng(2**40+123); repeated=api.rng(2**40+123)
            params=dict(mu=2.,lmbda=3.,a=4.,b=2.)
            batch=f.rvs_array(state,**params,n_samples=16)
            self.assertEqual(batch.status,0)
            singles=[f.rvs(repeated,**params).outputs for _ in range(16)]
            for key in ('x','s2'):
                np.testing.assert_array_equal(batch.outputs[key],[s[key] for s in singles])
            np.testing.assert_array_equal(state.state,repeated.state)
            state=api.rng(123); split=api.rng(123)
            whole=f.rvs_array(state,n_samples=7).outputs
            first=f.rvs_array(split,n_samples=2).outputs; rest=f.rvs_array(split,n_samples=5).outputs
            for key in ('x','s2'): np.testing.assert_array_equal(whole[key],np.r_[first[key],rest[key]])
            np.testing.assert_array_equal(state.state,split.state)
            saved=state.state
            empty=f.rvs_array(state,a=-1.,n_samples=0)
            self.assertEqual(empty.status,0); self.assertEqual(empty.outputs['x'].shape,(0,))
            for result in (f.rvs(state,a=-1.),f.rvs_array(state,a=-1.,n_samples=3)):
                self.assertEqual(result.status,0)
                self.assertTrue(all(np.all(np.isnan(v)) for v in result.outputs.values()))
            np.testing.assert_array_equal(state.state,saved)

    def test_buffers_and_raw_c_guards(self):
        ptr=lambda a:a.ctypes.data_as(C.c_void_p)
        for _,_,api in self.builds:
            f=api.normal_inverse_gamma; state=api.rng(12); saved=state.state
            for count in (-1,True,1.5,2**31):
                with self.assertRaises((TypeError,ValueError)): f.rvs_array(state,n_samples=count)
            for out in ({'x':np.empty(2),'s2':np.empty(3)},
                        {'x':np.empty(3,dtype=np.float32)},):
                with self.assertRaises((TypeError,ValueError)): f.rvs_array(state,n_samples=3,out=out)
            same=np.empty(3)
            with self.assertRaises(ValueError): f.rvs_array(state,n_samples=3,out={'x':same,'s2':same})
            np.testing.assert_array_equal(state.state,saved)
            out={'x':np.empty(3),'s2':np.empty(3)}
            result=f.rvs_array(state,n_samples=3,out=out)
            for k in out:self.assertIs(out[k],result.outputs[k])
            raw=f._functions['rvs_array']; serialized=state.state.copy(); before=serialized.copy()
            x=np.full(3,77.); s2=np.full(3,88.); status=C.c_int(-1)
            raw(-1,ptr(serialized),ptr(x),ptr(s2),0.,1.,1.,1.,C.byref(status))
            self.assertEqual(status.value,1)
            np.testing.assert_array_equal(serialized,before)
            np.testing.assert_array_equal(x,77.); np.testing.assert_array_equal(s2,88.)
            serialized[:]=0
            for n in (0,3):
                raw(n,ptr(serialized),ptr(x),ptr(s2),0.,1.,1.,1.,C.byref(status))
                self.assertEqual(status.value,4)
                np.testing.assert_array_equal(x,77.); np.testing.assert_array_equal(s2,88.)
                np.testing.assert_array_equal(serialized,0)

    def test_native_and_c_cpp_clients(self):
        for path,meta,_ in self.builds:
            with tempfile.TemporaryDirectory() as temp:
                folder=Path(temp); env=os.environ.copy()
                env['PATH']=str(path.parent)+os.pathsep+env.get('PATH','')
                flags=['-O0','-g','-fcheck=all','-fbacktrace'] if '-fcheck=all' in meta['flags'] else ['-O3']
                exe=folder/'native.exe'
                run([meta['compiler_path'],*flags,'-I',path.parent,ROOT/'test/test_normal_inverse_gamma_c_api.f90',path,'-o',exe])
                run([exe],env=env)
                client=folder/'client.c'
                client.write_text("""#include "scifort.h"
#include <math.h>
int main(void) {
 int64_t state[6]; int status; double x,s2,p;
 scifort_normal_inverse_gamma_pdf_f64(0,1,0,1,1,1,&p,&status);
 if(status || fabs(p-exp(-1)/sqrt(2*acos(-1)))>1e-14) return 1;
 scifort_normal_inverse_gamma_mean_f64(2,3,4,2,&x,&s2,&status);
 if(status || x!=2 || fabs(s2-2./3)>1e-14) return 2;
 scifort_rng_seed_i64(123,state,&status); if(status) return 3;
 scifort_normal_inverse_gamma_rvs_f64(state,&x,&s2,0,1,1,1,&status);
 return status || !isfinite(x) || !(s2>0);
}
""",encoding='ascii')
                for compiler,standard in (('gcc','c11'),('g++','c++17')):
                    if not shutil.which(compiler):continue
                    exe=folder/(compiler+'.exe')
                    run([compiler,'-std='+standard,'-Wall','-Wextra','-Werror','-I',ROOT/'include',client,path,'-o',exe])
                    run([exe],env=env)
