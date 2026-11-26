# SPDX-License-Identifier: MIT
"""Integer-input multinomial and multivariate hypergeometric integration contracts."""
import ctypes as C
import math
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

import numpy as np
from scipy.stats import multinomial, multivariate_hypergeom

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from build_library import build
from scifort import Scifort
from test_multivariate_t import run

ROOT=Path(__file__).resolve().parents[3]


class DiscreteMultivariateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran'):raise unittest.SkipTest('gfortran required')
        cls.builds=[]
        for checked in (True,False):
            path,metadata=build(checked=checked)
            cls.builds.append((path,metadata,Scifort(path,dll_directories=metadata['dll_directories'])))

    def test_all_ten_reference_methods_and_identities(self):
        for _,_,api in self.builds:
            a=api.multinomial;b=api.multivariate_hypergeom
            for n,p,x in ((5,[.2,.3,.5],[1,1,3]),(7,[.4,.6],[2,5]),(0,[.1,.9],[0,0])):
                for method in ('pmf','logpmf'):
                    np.testing.assert_allclose(getattr(a,method)(x,n,p),getattr(multinomial,method)(x,n,p),rtol=3e-13,atol=2e-14)
                for method in ('mean','cov','entropy'):
                    np.testing.assert_allclose(getattr(a,method)(n,p),getattr(multinomial,method)(n,p),rtol=3e-13,atol=2e-14)
                self.assertAlmostEqual(a.mean(n,p).sum(),n,places=13)
                np.testing.assert_allclose(a.cov(n,p).sum(axis=0),np.zeros(len(p)),rtol=0,atol=3e-15)
            for m,n,x in (([3,5,7],4,[1,1,2]),([0,5,2],4,[0,3,1]),([2,3],5,[2,3])):
                for method in ('pmf','logpmf'):
                    np.testing.assert_allclose(getattr(b,method)(x,m,n),getattr(multivariate_hypergeom,method)(x,m,n),rtol=3e-13,atol=2e-14)
                for method in ('mean','var','cov'):
                    np.testing.assert_allclose(getattr(b,method)(m,n),getattr(multivariate_hypergeom,method)(m,n),rtol=3e-13,atol=2e-14)
                np.testing.assert_allclose(np.diag(b.cov(m,n)),b.var(m,n),rtol=3e-15,atol=3e-17)
                np.testing.assert_allclose(b.cov(m,n).sum(axis=0),np.zeros(len(m)),rtol=0,atol=3e-15)
            for group,args in ((a,(5,[.2,.8])),(b,([3,5],4))):
                self.assertTrue(group.mean(*args).flags.f_contiguous)
                self.assertTrue(group.cov(*args).flags.f_contiguous)
                self.assertIsInstance(group.mean(*args),np.ndarray)

    def test_native_parameter_and_support_behavior(self):
        for _,_,api in self.builds:
            a=api.multinomial;b=api.multivariate_hypergeom
            # Native final-entry adjustment happens before the final range check.
            for p in ([.2,.3,.1],[.2,.3,99.],[.2,.3,-1.]):
                np.testing.assert_allclose(a.mean(10,p),[2,3,5],rtol=0,atol=1e-14)
                np.testing.assert_allclose(a.pmf([2,3,5],10,p),a.pmf([2,3,5],10,[.2,.3,.5]),rtol=0,atol=0)
            for p in ([-.2,1.2],[.7,.7,0.],[np.nan,.5],[np.inf,.5]):
                self.assertTrue(math.isnan(a.pmf([0]*len(p),0,p)))
                self.assertTrue(np.all(np.isnan(a.mean(1,p))))
            for group,params in ((a,(5,[.2,.8])),(b,([3,5],4))):
                for x in ([-1,2],[0,0]):
                    self.assertEqual(group.pmf(x,*params),0.)
                    self.assertEqual(group.logpmf(x,*params),-np.inf)
            self.assertEqual(a.pmf([1.5,3.5],5,[.2,.8]),0.)
            self.assertEqual(a.pmf([np.nan,5],5,[.2,.8]),0.)
            self.assertEqual(a.pmf([5],5,[.2]),1.)
            self.assertEqual(a.entropy(0,[.2,.8]),0.)
            self.assertTrue(math.isnan(a.entropy(-1,[.2,.8])))
            self.assertTrue(math.isnan(a.pmf([0,0],-1,[.2,.8])))
            for m,n in (([-1,5],2),([3,5],9),([3,5],-1)):
                self.assertTrue(math.isnan(b.pmf([0,0],m,n)))
                self.assertTrue(np.all(np.isnan(b.cov(m,n))))
            self.assertEqual(b.pmf([0,0],[0,0],0),1.)
            np.testing.assert_array_equal(b.mean([0,0],0),[0.,0.])
            np.testing.assert_array_equal(b.var([0,0],0),[0.,0.])
            np.testing.assert_array_equal(b.cov([1],1),[[0.]])
            limit=np.iinfo(np.intc).max-1
            np.testing.assert_array_equal(a.mean(limit,[1.]),[float(limit)])
            np.testing.assert_array_equal(b.mean([limit,0],0),[0.,0.])

    def test_integer_dtype_layout_shape_and_range_validation(self):
        for _,_,api in self.builds:
            a=api.multinomial;b=api.multivariate_hypergeom
            reference=b.pmf([1,3],[3,5],4)
            unaligned=np.ndarray((2,),dtype=np.intc,buffer=bytearray(2*np.dtype(np.intc).itemsize+1),offset=1)
            unaligned[:]=[3,5]
            readonly=np.array([3,5]);readonly.flags.writeable=False
            variants=([3,5],np.array([3,5],dtype='>i8'),np.array([3,5],dtype='uint64'),
                      np.array([3,99,5,99],dtype='int16')[::2],unaligned,readonly)
            for m in variants:
                self.assertEqual(b.pmf(np.array([1,3],dtype='>i4'),m,np.array(4,dtype='uint32')),reference)
            for n in (True,3.,'3',np.array([3])):
                with self.assertRaises((ValueError,TypeError)):a.mean(n,[.2,.8])
                with self.assertRaises((ValueError,TypeError)):b.mean([3,5],n)
            for m in ([3.,5.],[True,False],['3','5'],[3+1j,5],np.array([3,5],dtype=object)):
                with self.assertRaises(TypeError):b.mean(m,4)
            with self.assertRaises(TypeError):b.pmf([1.,3.],[3,5],4)
            for n in (2**31,2**31-1,-2**31-1,np.uint64(2**64-1)):
                with self.assertRaises(ValueError):a.mean(n,[.2,.8])
                with self.assertRaises(ValueError):b.mean([3,5],n)
            for m in ([2**31,0],[2**31-2,1],[-2**31,0],np.array([2**64-1,0],dtype='uint64')):
                with self.assertRaises(ValueError):b.mean(m,0)
            with self.assertRaises(ValueError):b.pmf([2**31-2,1],[3,5],4)
            for call in (lambda:a.pmf([1],4,[.2,.8]),lambda:b.pmf([1],[3,5],4),
                         lambda:a.mean(4,[]),lambda:b.mean(np.array([],dtype=int),0),
                         lambda:b.mean([[3,5]],4)):
                with self.assertRaises(ValueError):call()
            big=np.broadcast_to(np.ones(1,dtype=np.intc),(np.iinfo(np.intc).max+1,))
            with self.assertRaises(ValueError):b.mean(big,0)
            m=np.array([3,5])
            with self.assertRaises(ValueError):b.pmf(m,m,8)

    def test_caller_buffers_and_validation_before_writes(self):
        for _,_,api in self.builds:
            a=api.multinomial;b=api.multivariate_hypergeom
            for group,args in ((a,(4,[.2,.8])),(b,([3,5],4))):
                for method in ('mean','cov'):
                    shape=(2,) if method=='mean' else (2,2)
                    out=np.full(shape,99.,order='F')
                    self.assertIs(getattr(group,method)(*args,out=out),out)
                    np.testing.assert_allclose(out,getattr(group,method)(*args),rtol=0,atol=0)
                    for invalid in (np.empty(shape,dtype='float32'),np.empty((3,)),{}):
                        with self.assertRaises((ValueError,TypeError)):getattr(group,method)(*args,out=invalid)
                    out.flags.writeable=False
                    with self.assertRaises(ValueError):getattr(group,method)(*args,out=out)
            out=np.full(2,99.)
            for call in (lambda:a.mean(2**31-1,[.2,.8],out=out),lambda:b.mean([2**31-2,1],0,out=out)):
                with self.assertRaises(ValueError):call()
                np.testing.assert_array_equal(out,99.)
            self.assertIs(a.mean(-1,[.2,.8],out=out),out)
            self.assertTrue(np.all(np.isnan(out)))
            p=np.array([.2,.8])
            with self.assertRaises(ValueError):a.mean(4,p,out=p)
            storage=np.zeros(8,dtype=np.intc)
            with self.assertRaises(ValueError):b.mean(storage[:2],0,out=storage.view(np.float64)[:2])

    def test_c_integer_guards_and_output_preservation(self):
        for _,_,api in self.builds:
            for group in (api.multinomial,api.multivariate_hypergeom):
                hyper=group is api.multivariate_hypergeom
                x=np.array([1,3],dtype=np.intc) if hyper else np.array([1.,3.])
                params=np.array([3,5],dtype=np.intc) if hyper else np.array([.2,.8])
                xp=x.ctypes.data_as(C.c_void_p);pp=params.ctypes.data_as(C.c_void_p)
                for method in (('pmf','logpmf','mean','var','cov') if hyper
                               else ('pmf','logpmf','entropy','mean','cov')):
                    fn=group._functions[method]
                    out=np.full(4,99.);op=out.ctypes.data_as(C.POINTER(C.c_double));status=C.c_int(-1)
                    def call(n):
                        args=([2,xp,pp,n] if hyper else [2,xp,n,pp]) if method in ('pmf','logpmf') else ([2,pp,n] if hyper else [2,n,pp])
                        fn(*args,op,C.byref(status))
                    call(2147483647)
                    self.assertEqual(status.value,1);np.testing.assert_array_equal(out,99.)
                    call(-1)
                    self.assertEqual(status.value,0)
                    count=4 if method=='cov' else (2 if method in ('mean','var') else 1)
                    self.assertTrue(np.all(np.isnan(out[:count])))
                    if hyper:
                        params[:]=[2147483646,1];out.fill(99.);call(0)
                        self.assertEqual(status.value,1);np.testing.assert_array_equal(out,99.)
                        params[:]=[3,5]
                        if method in ('pmf','logpmf'):
                            x[:]=[-2147483648,0];out.fill(99.);call(0)
                            self.assertEqual(status.value,1);np.testing.assert_array_equal(out,99.)
                            x[:]=[1,3]

    def test_native_and_c_clients_all_ten_methods(self):
        for path,metadata,api in self.builds:
            with tempfile.TemporaryDirectory() as temp:
                folder=Path(temp);native=folder/'native.f90'
                native.write_text("""program native_client
use iso_c_binding, only: c_double,c_int
use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
use scifort_multinomial
use scifort_multivariate_hypergeom
use scifort_c_api_generated_multinomial
use scifort_c_api_generated_multivariate_hypergeom
implicit none
real(c_double) :: xr(2),p(2),actual,expected,v(2),mat(2,2)
integer(c_int) :: x(2),m(2),n,status
integer :: i
xr=[1.0_c_double,3.0_c_double];p=[0.2_c_double,0.8_c_double]
x=[1,3];m=[3,5]
do i=1,3
    n=4
    if (i==2) n=-1
    if (i==3) p=[0.2_c_double,0.1_c_double]
    expected=multinomial_pmf(xr,n,p)
    call scifort_multinomial_pmf_f64(2_c_int,xr,n,p,actual,status)
    call check()
    expected=multinomial_logpmf(xr,n,p)
    call scifort_multinomial_logpmf_f64(2_c_int,xr,n,p,actual,status)
    call check()
    expected=multinomial_entropy(n,p)
    call scifort_multinomial_entropy_f64(2_c_int,n,p,actual,status)
    call check()
    call scifort_multinomial_mean_f64(2_c_int,n,p,v,status)
    if (status/=0 .or. .not. all(equal(v,multinomial_mean(n,p)))) stop 1
    call scifort_multinomial_cov_f64(2_c_int,n,p,mat,status)
    if (status/=0 .or. .not. all(equal(mat,multinomial_cov(n,p)))) stop 2
    expected=multivariate_hypergeom_pmf(x,m,n)
    call scifort_multivariate_hypergeom_pmf_f64(2_c_int,x,m,n,actual,status)
    call check()
    expected=multivariate_hypergeom_logpmf(x,m,n)
    call scifort_multivariate_hypergeom_logpmf_f64(2_c_int,x,m,n,actual,status)
    call check()
    call scifort_multivariate_hypergeom_mean_f64(2_c_int,m,n,v,status)
    if (status/=0 .or. .not. all(equal(v,multivariate_hypergeom_mean(m,n)))) stop 3
    call scifort_multivariate_hypergeom_var_f64(2_c_int,m,n,v,status)
    if (status/=0 .or. .not. all(equal(v,multivariate_hypergeom_var(m,n)))) stop 4
    call scifort_multivariate_hypergeom_cov_f64(2_c_int,m,n,mat,status)
    if (status/=0 .or. .not. all(equal(mat,multivariate_hypergeom_cov(m,n)))) stop 5
end do
contains
pure elemental function equal(a,b) result(ok)
    implicit none
    real(c_double), intent(in) :: a !! actual
    real(c_double), intent(in) :: b !! expected
    logical :: ok
    ok=(a==b) .or. (ieee_is_nan(a).and.ieee_is_nan(b))
end function
subroutine check()
    implicit none
    if (status/=0 .or. .not. equal(actual,expected)) stop 6
end subroutine
end program
""",encoding='ascii')
                env=os.environ.copy();env['PATH']=str(path.parent)+os.pathsep+env.get('PATH','')
                flags=['-O0','-g','-fcheck=all','-fbacktrace'] if '-fcheck=all' in metadata['flags'] else ['-O3']
                exe=folder/'native.exe';run([metadata['compiler_path'],*flags,'-I',path.parent,native,path,'-o',exe]);run([exe],env=env)
                if shutil.which('gcc'):
                    c=folder/'client.c'
                    c.write_text("""#include "scifort.h"
#include <math.h>
int main(void) {
 double p[2]={.2,.8},xr[2]={1,3},v[2],mat[4],value; int x[2]={1,3},m[2]={3,5},status;
 scifort_multinomial_pmf_f64(2,xr,4,p,&value,&status);
 if (status || fabs(value-.4096)>1e-13) return 1;
 scifort_multinomial_logpmf_f64(2,xr,4,p,&value,&status);
 if (status || fabs(value-log(.4096))>1e-13) return 2;
 scifort_multinomial_entropy_f64(2,4,p,&value,&status); if (status || !isfinite(value)) return 3;
 scifort_multinomial_mean_f64(2,4,p,v,&status); if (status || fabs(v[0]-.8)>1e-14) return 4;
 scifort_multinomial_cov_f64(2,4,p,mat,&status); if (status || fabs(mat[0]-.64)>1e-14) return 5;
 scifort_multivariate_hypergeom_pmf_f64(2,x,m,4,&value,&status);
 if (status || fabs(value-3.0/7)>1e-13) return 6;
 scifort_multivariate_hypergeom_logpmf_f64(2,x,m,4,&value,&status);
 if (status || fabs(value-log(3.0/7))>1e-13) return 7;
 scifort_multivariate_hypergeom_mean_f64(2,m,4,v,&status); if (status || v[0]!=1.5) return 8;
 scifort_multivariate_hypergeom_var_f64(2,m,4,v,&status); if (status || fabs(v[0]-15.0/28)>1e-14) return 9;
 scifort_multivariate_hypergeom_cov_f64(2,m,4,mat,&status);
 return status || fabs(mat[0]-15.0/28)>1e-14 || fabs(mat[1]-mat[2])>1e-15;
}
""",encoding='ascii')
                    exe=folder/'c.exe';run(['gcc','-std=c11','-Wall','-Wextra','-Werror','-I',ROOT/'include',c,path,'-o',exe]);run([exe],env=env)
                    if shutil.which('g++'):run(['g++','-std=c++17','-Wall','-Wextra','-Werror','-I',ROOT/'include','-fsyntax-only',c])
