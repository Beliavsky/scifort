# SPDX-License-Identifier: MIT
"""Dirichlet scalar/array results, simplex length choices and C/native agreement."""
import ctypes as C
import math
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

import numpy as np
from scipy.stats import dirichlet

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from build_library import build
from scifort import Scifort
from test_multivariate_t import run

ROOT = Path(__file__).resolve().parents[3]
METHODS = ('pdf','logpdf','entropy','mean','var','cov')


class DirichletTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran'):
            raise unittest.SkipTest('gfortran required')
        cls.builds = []
        for checked in (True,False):
            path,metadata = build(checked=checked)
            cls.builds.append((path,metadata,Scifort(path,dll_directories=metadata['dll_directories'])))

    def test_reference_values_full_and_short_points_and_moments(self):
        for _,_,api in self.builds:
            group = api.dirichlet
            for alpha,x in (([2.,3.],[.25,.75]),([.3,2.,4.],[.2,.3,.5]),([1.,1.,1.,1.],[.1,.2,.3,.4])):
                for method in METHODS:
                    call = getattr(group,method)
                    actual = call(x,alpha) if method in ('pdf','logpdf') else call(alpha)
                    reference = getattr(dirichlet,method)(x,alpha) if method in ('pdf','logpdf') else getattr(dirichlet,method)(alpha)
                    np.testing.assert_allclose(actual,reference,rtol=3e-13,atol=3e-14)
                    if method in ('pdf','logpdf'):
                        np.testing.assert_allclose(call(x[:-1],alpha),reference,rtol=3e-13,atol=3e-14)
                    elif method in ('mean','var','cov'):
                        self.assertIsInstance(actual,np.ndarray)
                        self.assertTrue(actual.flags.f_contiguous)
                        self.assertEqual(actual.dtype,np.float64)
                    else:
                        self.assertIsInstance(actual,float)
                cov = group.cov(alpha)
                np.testing.assert_allclose(np.diag(cov),group.var(alpha),rtol=3e-15,atol=3e-17)
                np.testing.assert_allclose(cov.sum(axis=0),np.zeros(len(alpha)),rtol=0,atol=5e-17)
                self.assertAlmostEqual(group.mean(alpha).sum(),1.,places=14)
            self.assertEqual(group.pdf([1.],[.5]),1.)
            self.assertEqual(group.logpdf([1.],[.5]),0.)
            self.assertEqual(group.entropy([.5]),0.)
            np.testing.assert_array_equal(group.mean([.5]),[1.])
            np.testing.assert_array_equal(group.var([.5]),[0.])
            np.testing.assert_array_equal(group.cov([.5]),[[0.]])

    def test_simplex_boundaries_and_invalid_numerical_inputs(self):
        for _,_,api in self.builds:
            group = api.dirichlet
            for alpha in ([1.,1.],[2.,1.]):
                for method in ('pdf','logpdf'):
                    call = getattr(group,method)
                    np.testing.assert_allclose(call([0.,1.],alpha),getattr(dirichlet,method)([0.,1.],alpha))
            # Native behavior returns NaN at zero coordinates with alpha<1.
            self.assertTrue(math.isnan(group.pdf([0.,1.],[.5,2.])))
            self.assertTrue(math.isnan(group.logpdf([0.,1.],[.5,2.])))
            for method in METHODS:
                call = getattr(group,method)
                for alpha in ([0.,1.],[-1.,2.],[np.nan,1.],[np.inf,1.],[-np.inf,1.]):
                    actual = call([.2,.8],alpha) if method in ('pdf','logpdf') else call(alpha)
                    self.assertTrue(np.all(np.isnan(actual)))
            for point in ([.2,.2],[-.1,1.1],[np.nan,.5],[np.inf,.5],[1.1]):
                for method in ('pdf','logpdf'):
                    self.assertTrue(math.isnan(getattr(group,method)(point,[2.,3.])))
            # The unchanged native simplex sum tolerance is 1e-9.
            self.assertTrue(math.isfinite(group.logpdf([.2,.8+1e-10],[2.,3.])))
            self.assertTrue(math.isnan(group.logpdf([.2,.8+1e-8],[2.,3.])))

    def test_conversions_shape_choices_and_type_errors(self):
        for _,_,api in self.builds:
            group = api.dirichlet
            alpha=np.array([2.,3.])
            unaligned=np.ndarray((2,),dtype='float64',buffer=bytearray(17),offset=1)
            unaligned[:]=alpha
            readonly=alpha.copy();readonly.flags.writeable=False
            for aa in (alpha.tolist(),alpha.astype('float32'),alpha.astype('>f8'),np.array([2.,9.,3.,9.])[::2],unaligned,readonly):
                np.testing.assert_allclose(group.mean(aa),[.4,.6],rtol=2e-15)
                np.testing.assert_allclose(group.logpdf(np.array([.25],dtype='>f8'),aa),group.logpdf([.25],alpha),rtol=2e-15)
            for method in METHODS:
                call=getattr(group,method)
                for bad in ([],2.,[[2.,3.]],[2+1j,3],['2','3']):
                    with self.assertRaises((ValueError,TypeError)):
                        call([.2,.8],bad) if method in ('pdf','logpdf') else call(bad)
            for method in ('pdf','logpdf'):
                call=getattr(group,method)
                for point,aa in (([],[2.]),([.2],[2.,3.,4.]),([.2,.3,.5],[2.,3.]),([[.2,.8]],[2.,3.]),([1+1j],[2.,3.])):
                    with self.assertRaises((ValueError,TypeError)):
                        call(point,aa)
                with self.assertRaises(ValueError):
                    call(alpha,alpha)
                with self.assertRaises(TypeError):
                    call([.2,.8],alpha,out=np.zeros(1))
            with self.assertRaises(ValueError):
                group.mean(np.broadcast_to(np.ones(1),(np.iinfo(np.intc).max+1,)))
            np.testing.assert_array_equal(alpha,[2.,3.])

    def test_caller_owned_array_results_and_output_safety(self):
        for _,_,api in self.builds:
            group=api.dirichlet
            for method in ('mean','var','cov'):
                call=getattr(group,method)
                shape=(2,2) if method=='cov' else (2,)
                out=np.full(shape,99.,order='F')
                result=call([2.,3.],out=out)
                self.assertIs(result,out)
                np.testing.assert_allclose(out,getattr(dirichlet,method)([2.,3.]),rtol=3e-14,atol=3e-17)
                self.assertIs(call([0.,3.],out=out),out)
                self.assertTrue(np.all(np.isnan(out)))
                out.fill(99.)
                with self.assertRaises(ValueError):
                    call([],out=out)
                np.testing.assert_array_equal(out,99.)
                for bad in (np.empty(shape,dtype='float32'),np.empty(shape,dtype='>f8'),np.empty((3,)),{},[0.,0.]):
                    with self.assertRaises((ValueError,TypeError)):
                        call([2.,3.],out=bad)
                readonly=np.empty(shape,order='F');readonly.flags.writeable=False
                with self.assertRaises(ValueError):
                    call([2.,3.],out=readonly)
                unaligned=np.ndarray(shape,dtype='float64',buffer=bytearray(8*np.prod(shape)+1),offset=1,order='F')
                with self.assertRaises(ValueError):
                    call([2.,3.],out=unaligned)
            with self.assertRaises(ValueError):
                group.cov([2.,3.],out=np.empty((2,2),order='C'))
            storage=np.full(4,2.)
            for method in ('mean','var'):
                with self.assertRaises(ValueError):
                    getattr(group,method)(storage[:2],out=storage[:2])
            with self.assertRaises(ValueError):
                group.cov(storage[:2],out=storage.reshape((2,2),order='F'))
            np.testing.assert_array_equal(storage,2.)

    def test_c_shape_status_preservation_and_numerical_nan(self):
        for _,_,api in self.builds:
            alpha=np.array([2.,3.,4.]);x=np.array([.2,.8,.1])
            ap=alpha.ctypes.data_as(C.c_void_p);xp=x.ctypes.data_as(C.c_void_p)
            for method in METHODS:
                fn=api.dirichlet._functions[method]
                out=np.full(9,99.)
                op=out.ctypes.data_as(C.POINTER(C.c_double))
                status=C.c_int(-1)
                for k in (0,-1):
                    args=[k,1,xp,ap] if method in ('pdf','logpdf') else [k,ap]
                    fn(*args,op,C.byref(status))
                    self.assertEqual(status.value,1)
                    np.testing.assert_array_equal(out,99.)
                if method in ('pdf','logpdf'):
                    for k,m in ((1,0),(2,3),(3,1)):
                        fn(k,m,xp,ap,op,C.byref(status))
                        self.assertEqual(status.value,1)
                        np.testing.assert_array_equal(out,99.)
                alpha[0]=0.
                args=[2,2,xp,ap] if method in ('pdf','logpdf') else [2,ap]
                fn(*args,op,C.byref(status))
                self.assertEqual(status.value,0)
                count=4 if method=='cov' else (2 if method in ('mean','var') else 1)
                self.assertTrue(np.all(np.isnan(out[:count])))
                alpha[0]=2.

    def test_native_and_c_clients_all_six_methods(self):
        for path,metadata,api in self.builds:
            with tempfile.TemporaryDirectory() as temp:
                folder=Path(temp);client=folder/'native.f90'
                client.write_text("""program native_client
use iso_c_binding, only: c_double,c_int
use, intrinsic :: ieee_arithmetic, only: ieee_is_nan,ieee_value,ieee_quiet_nan,ieee_positive_inf
use scifort_dirichlet
use scifort_c_api_generated_dirichlet
implicit none
real(c_double) :: alpha(3),x(3),actual(3),matrix(3,3),scalar,expected
integer(c_int) :: status
integer :: i,m
do i=1,5
    alpha=[2.0_c_double,3.0_c_double,4.0_c_double]
    x=[0.2_c_double,0.3_c_double,0.5_c_double]
    if (i==2) alpha(1)=0
    if (i==3) alpha(1)=ieee_value(0.0_c_double,ieee_quiet_nan)
    if (i==4) alpha(1)=ieee_value(0.0_c_double,ieee_positive_inf)
    if (i==5) x(1)=-1
    do m=2,3
        expected=dirichlet_pdf(x(:m),alpha)
        call scifort_dirichlet_pdf_f64(3_c_int,int(m,c_int),x,alpha,scalar,status)
        call check_scalar()
        expected=dirichlet_logpdf(x(:m),alpha)
        call scifort_dirichlet_logpdf_f64(3_c_int,int(m,c_int),x,alpha,scalar,status)
        call check_scalar()
    end do
    expected=dirichlet_entropy(alpha)
    call scifort_dirichlet_entropy_f64(3_c_int,alpha,scalar,status)
    call check_scalar()
    call scifort_dirichlet_mean_f64(3_c_int,alpha,actual,status)
    if (status/=0) stop 1
    if (i<2 .or. i==5) then
        if (any(actual/=dirichlet_mean(alpha))) stop 2
    else
        if (.not. all(ieee_is_nan(actual))) stop 3
    end if
    call scifort_dirichlet_var_f64(3_c_int,alpha,actual,status)
    if (status/=0) stop 4
    if (i<2 .or. i==5) then
        if (any(actual/=dirichlet_var(alpha))) stop 5
    else
        if (.not. all(ieee_is_nan(actual))) stop 6
    end if
    call scifort_dirichlet_cov_f64(3_c_int,alpha,matrix,status)
    if (status/=0) stop 7
    if (i<2 .or. i==5) then
        if (any(matrix/=dirichlet_cov(alpha))) stop 8
    else
        if (.not. all(ieee_is_nan(matrix))) stop 9
    end if
end do
contains
subroutine check_scalar()
    implicit none
    if (status/=0) stop 10
    if (ieee_is_nan(expected)) then
        if (.not. ieee_is_nan(scalar)) stop 11
    else if (expected/=scalar) then
        stop 12
    end if
end subroutine
end program
""",encoding='ascii')
                env=os.environ.copy();env['PATH']=str(path.parent)+os.pathsep+env.get('PATH','')
                flags=['-O0','-g','-fcheck=all','-fbacktrace'] if '-fcheck=all' in metadata['flags'] else ['-O3']
                exe=folder/'native.exe'
                run([metadata['compiler_path'],*flags,'-I',path.parent,client,path,'-o',exe])
                run([exe],env=env)
                if shutil.which('gcc'):
                    c=folder/'client.c'
                    c.write_text("""#include "scifort.h"
#include <math.h>
int main(void) {
 double alpha[2]={2,3},x[2]={0.25,0.75},value,vector[2],matrix[4]; int status;
 scifort_dirichlet_pdf_f64(2,2,x,alpha,&value,&status);
 if (status || fabs(value-1.6875)>1e-13) return 1;
 scifort_dirichlet_logpdf_f64(2,1,x,alpha,&value,&status);
 if (status || fabs(value-log(1.6875))>1e-13) return 2;
 scifort_dirichlet_entropy_f64(2,alpha,&value,&status);
 if (status || !isfinite(value)) return 3;
 scifort_dirichlet_mean_f64(2,alpha,vector,&status);
 if (status || fabs(vector[0]-.4)>1e-15 || fabs(vector[1]-.6)>1e-15) return 4;
 scifort_dirichlet_var_f64(2,alpha,vector,&status);
 if (status || fabs(vector[0]-.04)>1e-15) return 5;
 scifort_dirichlet_cov_f64(2,alpha,matrix,&status);
 return status || fabs(matrix[0]-.04)>1e-15 || fabs(matrix[1]+.04)>1e-15 || matrix[1]!=matrix[2];
}
""",encoding='ascii')
                    exe=folder/'c.exe'
                    run(['gcc','-std=c11','-Wall','-Wextra','-Werror','-I',ROOT/'include',c,path,'-o',exe])
                    run([exe],env=env)
                    if shutil.which('g++'):
                        run(['g++','-std=c++17','-Wall','-Wextra','-Werror','-I',ROOT/'include','-fsyntax-only',c])
