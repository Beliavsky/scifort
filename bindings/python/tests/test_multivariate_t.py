# SPDX-License-Identifier: MIT
"""Student t density adapters: shapes, defaults, singularity and native agreement."""
import ctypes as C
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np
from scipy.stats import multivariate_t

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from build_library import build
from scifort import Scifort

ROOT = Path(__file__).resolve().parents[3]


def run(command, **kwargs):
    result = subprocess.run([str(s) for s in command],capture_output=True,text=True,**kwargs)
    if result.returncode:
        raise AssertionError(f'{command}\n{result.stdout}\n{result.stderr}')


class MultivariateTTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran'):
            raise unittest.SkipTest('gfortran required')
        cls.builds = []
        for checked in (True,False):
            path,metadata = build(checked=checked)
            cls.builds.append((path,metadata,Scifort(path,dll_directories=metadata['dll_directories'])))

    def test_reference_values_defaults_and_normal_limit(self):
        for _,_,api in self.builds:
            group = api.multivariate_t
            for x,loc,shape in (([.3],[.1],[[2.]]),
                                 ([.2,-.4],[.5,-.1],[[2.,.3],[.3,1.]]),
                                 ([1.,2.,3.],[0.,1.,-1.],[[4.,1.,2.],[1.,3.,.5],[2.,.5,5.]])):
                for df in (.5,1.,5.,30.,100.):
                    for method in ('pdf','logpdf'):
                        value = getattr(group,method)(x,loc,shape,df)
                        self.assertIsInstance(value,float)
                        reference = getattr(multivariate_t,method)(x,loc,shape,df)
                        np.testing.assert_allclose(value,reference,rtol=3e-13,atol=2e-14)
                for method in ('pdf','logpdf'):
                    value = getattr(group,method)(x,loc,shape,float('inf'))
                    reference = getattr(api.multivariate_normal,method)(x,loc,shape)
                    self.assertEqual(value,reference)
            for method in ('pdf','logpdf'):
                call = getattr(group,method)
                for kwargs in ({},{'loc':[.2,-.3]},{'shape':[[2.,.1],[.1,3.]]},{'df':5.}):
                    reference_args = dict(loc=np.zeros(2),shape=np.eye(2),df=1.)
                    reference_args.update(kwargs)
                    np.testing.assert_allclose(call([.5,-.7],**kwargs),
                        getattr(multivariate_t,method)([.5,-.7],**reference_args),rtol=3e-14,atol=2e-14)
                self.assertEqual(call([.5,-.7]),call([.5,-.7],None,None,1.))
            self.assertEqual(group.pdf([1e100,-1e100],df=5.),0.)
            self.assertTrue(math.isfinite(group.logpdf([1e100,-1e100],df=5.)))

    def test_singular_policies_invalid_parameters_and_shape_semantics(self):
        for _,_,api in self.builds:
            group = api.multivariate_t
            shape = np.diag([4.,0.])
            for x in ([1.,0.],[1.,2.]):
                self.assertTrue(math.isnan(group.pdf(x,shape=shape)))
                pdf = group.pdf(x,shape=shape,df=5.,allow_singular=np.bool_(True))
                logpdf = group.logpdf(x,shape=shape,df=5.)
                np.testing.assert_allclose(pdf,multivariate_t.pdf(x,shape=shape,df=5.,allow_singular=True),rtol=3e-14)
                np.testing.assert_allclose(logpdf,multivariate_t.logpdf(x,shape=shape,df=5.),rtol=3e-14)
                np.testing.assert_allclose(pdf,math.exp(logpdf),rtol=3e-15)
            # Finite-df native/SciPy formulas project onto the nonzero eigenspace.
            self.assertEqual(group.logpdf([1.,0.],shape=shape,df=5.),
                             group.logpdf([1.,2.],shape=shape,df=5.))
            # At infinite df, the native normal branch applies its support mask.
            self.assertEqual(group.logpdf([1.,2.],shape=shape,df=np.inf),-np.inf)
            self.assertEqual(group.pdf([1.,2.],shape=shape,df=np.inf,allow_singular=True),0.)
            self.assertTrue(math.isnan(group.pdf([1.,0.],shape=shape,df=np.inf)))
            for method in ('pdf','logpdf'):
                call = getattr(group,method)
                for df in (0.,-1.,-np.inf,np.nan):
                    self.assertTrue(math.isnan(call([0.,0.],df=df)))
                for bad in (np.diag([-1.,2.]),np.diag([np.nan,1.]),np.diag([np.inf,1.]),
                            np.array([[2.,np.nan],[.2,3.]])):
                    self.assertTrue(math.isnan(call([0.,0.],shape=bad)))
                self.assertTrue(math.isnan(call([np.nan,0.])))
                self.assertTrue(math.isnan(call([0.,0.],loc=[np.inf,0.])))
                self.assertEqual(call([.3,.2],shape=[[2.,99.],[.2,3.]]),
                                 call([.3,.2],shape=[[2.,.2],[.2,3.]]))
            with self.assertRaises(TypeError):
                group.logpdf([0.,0.],allow_singular=True)

    def test_layout_conversions_and_validation(self):
        for _,_,api in self.builds:
            group = api.multivariate_t
            x = np.array([.25,-.5])
            shape = np.array([[2.,.125],[.125,3.]])
            unaligned = np.ndarray((2,),dtype='float64',buffer=bytearray(17),offset=1)
            unaligned[:] = x
            readonly = shape.copy()
            readonly.flags.writeable = False
            for xx in (x.tolist(),x.astype('float32'),x.astype('>f8'),np.array([.25,9.,-.5,9.])[::2],unaligned):
                for ss in (shape,shape.astype('float32'),shape.astype('>f8'),np.asfortranarray(shape),readonly):
                    np.testing.assert_allclose(group.logpdf(xx,shape=ss,df=np.array(5.)),
                        group.logpdf(x,shape=shape,df=5.),rtol=2e-14,atol=0)
            for method in ('pdf','logpdf'):
                call = getattr(group,method)
                for xx in ([],.5,[[1.,2.]],np.ones((2,2)),[1+2j],['x']):
                    with self.assertRaises((ValueError,TypeError)):
                        call(xx)
                for kwargs in ({'loc':[0.]},{'shape':np.eye(3)},{'shape':[1.,1.]},
                               {'df':[5.]},{'df':1+2j},{'df':'5'}):
                    with self.assertRaises((ValueError,TypeError)):
                        call([0.,0.],**kwargs)
                with self.assertRaises(ValueError):
                    call(x,loc=x)
                with self.assertRaises(ValueError):
                    call(np.broadcast_to(np.zeros(1),(np.iinfo(np.intc).max+1,)))
                with self.assertRaises(TypeError):
                    call([0.,0.],out={})
            for flag in (1,np.array(True)):
                with self.assertRaises(TypeError):
                    group.pdf([0.,0.],allow_singular=flag)
            np.testing.assert_array_equal(x,[.25,-.5])
            np.testing.assert_array_equal(shape,[[2.,.125],[.125,3.]])

    def test_c_call_status_and_output_preservation(self):
        for _,_,api in self.builds:
            arrays = [np.array([.2,-.1]),np.zeros(2),np.eye(2,order='F')]
            pointers = [a.ctypes.data_as(C.c_void_p) for a in arrays]
            for method in ('pdf','logpdf'):
                function = api.multivariate_t._functions[method]
                for n in (0,-1):
                    output,status = C.c_double(123.),C.c_int(-1)
                    flags = [0] if method=='pdf' else []
                    function(n,*pointers,5.,*flags,C.byref(output),C.byref(status))
                    self.assertEqual(status.value,1)
                    self.assertEqual(output.value,123.)
                output,status = C.c_double(123.),C.c_int(-1)
                if method=='pdf':
                    function(2,*pointers,5.,7,C.byref(output),C.byref(status))
                    self.assertEqual(status.value,1)
                    self.assertEqual(output.value,123.)
                function(2,*pointers,-1.,*flags,C.byref(output),C.byref(status))
                self.assertEqual(status.value,0)
                self.assertTrue(math.isnan(output.value))

    def test_native_fortran_and_c_clients(self):
        for path,metadata,api in self.builds:
            with tempfile.TemporaryDirectory() as temp:
                folder = Path(temp)
                client = folder/'native_client.f90'
                client.write_text("""program native_client
use iso_c_binding, only: c_double, c_int
use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan, ieee_positive_inf
use scifort_multivariate_t, only: multivariate_t_pdf, multivariate_t_logpdf
use scifort_c_api_generated_multivariate_t
implicit none
real(c_double) :: x(2), loc(2), shape(2,2), df, actual, expected
integer(c_int) :: status
integer :: i
logical :: allow
loc=0
do i=1,9
    x=[0.5_c_double,-0.25_c_double]
    shape=reshape([4.0_c_double,0.0_c_double,0.0_c_double,1.0_c_double],[2,2])
    df=5
    allow=.false.
    if (i==2 .or. i==3) then
        shape(2,2)=0
        allow=i==3
    end if
    if (i==4) x(2)=ieee_value(0.0_c_double,ieee_quiet_nan)
    if (i==5) x(1)=ieee_value(0.0_c_double,ieee_positive_inf)
    if (i==6) shape(1,1)=-1
    if (i==7) df=-1
    if (i==8) df=ieee_value(0.0_c_double,ieee_positive_inf)
    if (i==9) df=ieee_value(0.0_c_double,ieee_quiet_nan)
    expected=multivariate_t_pdf(x,loc,shape,df,allow)
    call scifort_multivariate_t_pdf_f64(2_c_int,x,loc,shape,df, &
        merge(1_c_int,0_c_int,allow),actual,status)
    call check_result()
    expected=multivariate_t_logpdf(x,loc,shape,df)
    call scifort_multivariate_t_logpdf_f64(2_c_int,x,loc,shape,df,actual,status)
    call check_result()
end do
contains
subroutine check_result()
    implicit none
    if (status/=0) stop 1
    if (ieee_is_nan(expected)) then
        if (.not. ieee_is_nan(actual)) stop 2
    else if (expected/=actual) then
        stop 3
    end if
end subroutine
end program
""",encoding='ascii')
                env = os.environ.copy()
                env['PATH'] = str(path.parent)+os.pathsep+env.get('PATH','')
                flags = ['-O0','-g','-fcheck=all','-fbacktrace'] if '-fcheck=all' in metadata['flags'] else ['-O3']
                exe = folder/'native.exe'
                run([metadata['compiler_path'],*flags,'-I',path.parent,client,path,'-o',exe])
                run([exe],env=env)
                if shutil.which('gcc'):
                    c = folder/'client.c'
                    c.write_text("""#include "scifort.h"
#include <math.h>
int main(void) {
 double x[2]={0,0}, loc[2]={0,0}, shape[4]={1,0,0,1}, value; int status;
 scifort_multivariate_t_pdf_f64(2,x,loc,shape,5,0,&value,&status);
 if (status || fabs(value-0.15915494309189535)>1e-15) return 1;
 scifort_multivariate_t_logpdf_f64(2,x,loc,shape,5,&value,&status);
 return status || fabs(value+1.8378770664093453)>1e-14;
}
""",encoding='ascii')
                    exe = folder/'c.exe'
                    run(['gcc','-std=c11','-Wall','-Wextra','-Werror','-I',ROOT/'include',c,path,'-o',exe])
                    run([exe],env=env)
                    if shutil.which('g++'):
                        run(['g++','-std=c++17','-Wall','-Wextra','-Werror','-I',ROOT/'include','-fsyntax-only',c])
