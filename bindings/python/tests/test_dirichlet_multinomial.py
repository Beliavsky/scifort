# SPDX-License-Identifier: MIT
"""Dirichlet-multinomial C/Python contracts and count conversion regressions."""
import ctypes as C
import math
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

import numpy as np
from scipy.stats import dirichlet_multinomial as reference

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from build_library import build
from scifort import Scifort
from test_multivariate_t import run

ROOT = Path(__file__).resolve().parents[3]


class DirichletMultinomialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran'):
            raise unittest.SkipTest('gfortran required')
        cls.builds = []
        for checked in (True, False):
            path, metadata = build(checked=checked)
            cls.builds.append((path, metadata, Scifort(
                path, dll_directories=metadata['dll_directories'])))

    def test_references_and_moment_identities(self):
        for _, _, api in self.builds:
            group = api.dirichlet_multinomial
            for alpha, n, x in (([2.,3.,4.],6,[1,2,3]), ([.2,.7],8,[2,6]),
                                 ([1.],5,[5])):
                for method in ('pmf','logpmf'):
                    np.testing.assert_allclose(getattr(group,method)(x,alpha,n),
                        getattr(reference,method)(x,alpha,n),rtol=4e-13,atol=3e-14)
                for method in ('mean','var','cov'):
                    np.testing.assert_allclose(getattr(group,method)(alpha,n),
                        getattr(reference,method)(alpha,n),rtol=4e-13,atol=3e-14)
                self.assertAlmostEqual(group.mean(alpha,n).sum(),n,places=13)
                np.testing.assert_allclose(np.diag(group.cov(alpha,n)),group.var(alpha,n),
                                           rtol=2e-15,atol=1e-16)
                np.testing.assert_allclose(group.cov(alpha,n).sum(axis=0),0,rtol=0,atol=1e-14)

    def test_native_invalid_counts_support_and_parameters(self):
        for _, _, api in self.builds:
            group = api.dirichlet_multinomial
            for bad in (-1.,.5,np.nan,np.inf,-np.inf,-1e100):
                for method in ('pmf','logpmf'):
                    self.assertTrue(math.isnan(getattr(group,method)([bad,0],[2.,3.],0)))
            for x in ([0,0],[2**31,0],[1e100,0]):
                self.assertEqual(group.pmf(x,[2.,3.],1),0.)
                self.assertEqual(group.logpmf(x,[2.,3.],1),-np.inf)
            for alpha in ([0.,3.],[-1.,3.],[np.nan,3.],[np.inf,3.]):
                self.assertTrue(math.isnan(group.pmf([1,2],alpha,3)))
                for method in ('mean','var','cov'):
                    self.assertTrue(np.all(np.isnan(getattr(group,method)(alpha,3))))
            self.assertTrue(math.isnan(group.logpmf([0,0],[2.,3.],-1)))
            np.testing.assert_array_equal(group.mean([2.,3.],0),[0.,0.])
            np.testing.assert_array_equal(group.cov([2.,3.],0),np.zeros((2,2)))
            self.assertEqual(group.pmf([0,0],[2.,3.],0),1.)
            self.assertEqual(group.pmf([1],[1.],1),1.)

    def test_python_validation_and_buffers(self):
        for _, _, api in self.builds:
            group = api.dirichlet_multinomial
            for method in ('mean','var','cov'):
                shape = (2,2) if method=='cov' else (2,)
                out = np.full(shape,99.,order='F')
                self.assertIs(getattr(group,method)([2.,3.],4,out=out),out)
                np.testing.assert_array_equal(out,getattr(group,method)([2.,3.],4))
                for n in (2147483647,2**31,-2**31-1,True,4.,'4'):
                    out.fill(99.)
                    with self.assertRaises((ValueError,TypeError)):
                        getattr(group,method)([2.,3.],n,out=out)
                    np.testing.assert_array_equal(out,99.)
                for invalid in (np.zeros(shape,dtype='float32'),np.zeros(3),{}):
                    with self.assertRaises((ValueError,TypeError)):
                        getattr(group,method)([2.,3.],4,out=invalid)
                out.flags.writeable=False
                with self.assertRaises(ValueError):getattr(group,method)([2.,3.],4,out=out)
            for alpha in ([],[[2.,3.]]):
                with self.assertRaises(ValueError):group.mean(alpha,4)
            with self.assertRaises(ValueError):group.pmf([1],[2.,3.],4)
            with self.assertRaises(TypeError):group.mean([2+1j,3],4)
            alpha=np.array([2.,3.])
            with self.assertRaises(ValueError):group.mean(alpha,4,out=alpha)
            for alpha in (np.array([2.,99.,3.])[::2],np.array([2.,3.],dtype='>f8')):
                np.testing.assert_allclose(group.mean(alpha,5),[2.,3.],rtol=0,atol=1e-15)

    def test_raw_c_guards_preserve_outputs(self):
        for _, _, api in self.builds:
            group=api.dirichlet_multinomial
            alpha=np.array([2.,3.]);x=np.array([1.,3.])
            ap=alpha.ctypes.data_as(C.c_void_p);xp=x.ctypes.data_as(C.c_void_p)
            for method,fn in group._functions.items():
                out=np.full(4,99.);op=out.ctypes.data_as(C.POINTER(C.c_double))
                status=C.c_int(-1)
                def call(k,n):
                    args=[k,xp,ap,n] if method in ('pmf','logpmf') else [k,ap,n]
                    fn(*args,op,C.byref(status))
                for k,n in ((2,2147483647),(0,4),(-1,4)):
                    call(k,n)
                    self.assertEqual(status.value,1);np.testing.assert_array_equal(out,99.)
                call(2,-1)
                self.assertEqual(status.value,0)
                count=4 if method=='cov' else (2 if method in ('mean','var') else 1)
                self.assertTrue(np.all(np.isnan(out[:count])))

    def test_native_regressions_and_c_clients(self):
        for path,metadata,_ in self.builds:
            with tempfile.TemporaryDirectory() as temp:
                folder=Path(temp)
                env=os.environ.copy();env['PATH']=str(path.parent)+os.pathsep+env.get('PATH','')
                flags=['-O0','-g','-fcheck=all','-fbacktrace','-ftrapv'] if '-fcheck=all' in metadata['flags'] else ['-O3']
                exe=folder/'counts.exe'
                run([metadata['compiler_path'],*flags,'-I',path.parent,
                     ROOT/'test/test_dirichlet_multinomial_counts.f90',path,'-o',exe])
                run([exe],env=env)
                exe=folder/'agreement.exe'
                run([metadata['compiler_path'],*flags,'-I',path.parent,
                     ROOT/'test/test_dirichlet_multinomial_c_api.f90',path,'-o',exe])
                run([exe],env=env)
                # Also compile the changed numerical module with traps, rather than
                # relying on flags of the shared library used by the other checks.
                exe=folder/'trapped.exe'
                run([metadata['compiler_path'],*flags,'-ffpe-trap=invalid','-I',path.parent,
                     '-J',folder,ROOT/'src/stats/dirichlet_multinomial.f90',
                     ROOT/'test/test_dirichlet_multinomial_counts.f90',path,'-o',exe])
                run([exe],env=env)
                if shutil.which('gcc'):
                    c=folder/'client.c'
                    c.write_text("""#include "scifort.h"
#include <math.h>
int main(void) {
 double x[2]={1,3},a[2]={2,3},v[2],mat[4],value; int status;
 scifort_dirichlet_multinomial_pmf_f64(2,x,a,4,&value,&status);
 if (status || fabs(value-2.0/7)>1e-13) return 1;
 scifort_dirichlet_multinomial_logpmf_f64(2,x,a,4,&value,&status);
 if (status || fabs(value-log(2.0/7))>1e-13) return 2;
 scifort_dirichlet_multinomial_mean_f64(2,a,4,v,&status);
 if (status || fabs(v[0]-1.6)>1e-14) return 3;
 scifort_dirichlet_multinomial_var_f64(2,a,4,v,&status);
 if (status || fabs(v[0]-1.44)>1e-14) return 4;
 scifort_dirichlet_multinomial_cov_f64(2,a,4,mat,&status);
 return status || fabs(mat[0]-1.44)>1e-14 || fabs(mat[1]+1.44)>1e-14;
}
""",encoding='ascii')
                    exe=folder/'c.exe'
                    run(['gcc','-std=c11','-Wall','-Wextra','-Werror','-I',ROOT/'include',c,path,'-o',exe])
                    run([exe],env=env)
                    if shutil.which('g++'):
                        run(['g++','-std=c++17','-Wall','-Wextra','-Werror','-I',ROOT/'include','-fsyntax-only',c])


if __name__ == '__main__':
    unittest.main()
