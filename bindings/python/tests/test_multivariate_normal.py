# SPDX-License-Identifier: MIT
"""Multivariate density shapes, native/C agreement and Python conversion contracts."""
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
from scipy.stats import multivariate_normal

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from build_library import build
from scifort import Scifort

ROOT = Path(__file__).resolve().parents[3]


class MultivariateNormalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran'):
            raise unittest.SkipTest('gfortran required')
        cls.builds = []
        for checked in (True,False):
            path,metadata = build(checked=checked)
            api = Scifort(path,dll_directories=metadata['dll_directories'])
            cls.builds.append((path,metadata,api))

    def test_reference_values_defaults_and_one_dimension(self):
        for _,_,api in self.builds:
            group = api.multivariate_normal
            for x,mean,cov in (([.3],[.1],[[2.]]),
                               ([.2,-.4],[.5,-.1],[[2.,.3],[.3,1.]]),
                               ([1.,2.,3.],[0.,1.,-1.],[[4.,1.,2.],[1.,3.,.5],[2.,.5,5.]])):
                for method in ('pdf','logpdf'):
                    value = getattr(group,method)(x,mean,cov)
                    self.assertIsInstance(value,float)
                    reference = getattr(multivariate_normal,method)(x,mean,cov)
                    np.testing.assert_allclose(value,reference,rtol=2e-13,atol=2e-14)
            for method in ('pdf','logpdf'):
                call = getattr(group,method)
                for kwargs in ({},{'mean':[.2,-.3]},{'cov':[[2.,.1],[.1,3.]]}):
                    value = call([.5,-.7],**kwargs)
                    reference_args = dict(mean=np.zeros(2),cov=np.eye(2))
                    reference_args.update(kwargs)
                    reference = getattr(multivariate_normal,method)([.5,-.7],**reference_args)
                    np.testing.assert_allclose(value,reference,rtol=2e-14,atol=2e-14)
                self.assertEqual(call([.5,-.7]),call([.5,-.7],None,None,False))
            self.assertEqual(group.pdf([1000.,-1000.]),0.)
            self.assertTrue(math.isfinite(group.logpdf([1000.,-1000.])))

    def test_singular_invalid_parameters_and_nonfinite_values(self):
        for _,_,api in self.builds:
            group = api.multivariate_normal
            cov = np.diag([4.,0.])
            for method in ('pdf','logpdf'):
                call = getattr(group,method)
                self.assertTrue(math.isnan(call([1.,0.],cov=cov)))
                for x in ([1.,0.],[1.,1.]):
                    value = call(x,cov=cov,allow_singular=np.bool_(True))
                    reference = getattr(multivariate_normal,method)(x,cov=cov,allow_singular=True)
                    np.testing.assert_allclose(value,reference,rtol=2e-14,atol=2e-14)
                for bad in (np.diag([-1.,2.]),np.diag([np.nan,1.]),np.diag([np.inf,1.])):
                    self.assertTrue(math.isnan(call([0.,0.],cov=bad)))
                self.assertTrue(math.isnan(call([np.nan,0.])))
                self.assertTrue(math.isnan(call([0.,0.],mean=[np.inf,0.])))
                # All covariance entries must be finite, even the unused upper triangle.
                bad = np.array([[2.,np.nan],[.2,3.]])
                self.assertTrue(math.isnan(call([0.,0.],cov=bad)))
                lower = np.array([[2.,99.],[.2,3.]])
                symmetric = np.array([[2.,.2],[.2,3.]])
                self.assertEqual(call([.3,.2],cov=lower),call([.3,.2],cov=symmetric))

    def test_layout_conversions_and_shape_errors(self):
        for _,_,api in self.builds:
            group = api.multivariate_normal
            x = np.array([.25,-.5])
            cov = np.array([[2.,.125],[.125,3.]])
            unaligned = np.ndarray((2,),dtype='float64',buffer=bytearray(17),offset=1)
            unaligned[:] = x
            readonly = cov.copy(); readonly.flags.writeable = False
            for xx in (x.tolist(),x.astype('float32'),x.astype('>f8'),np.array([.25,9.,-.5,9.])[::2],unaligned):
                for cc in (cov,cov.astype('float32'),cov.astype('>f8'),np.asfortranarray(cov),readonly):
                    np.testing.assert_allclose(group.logpdf(xx,cov=cc),group.logpdf(x,cov=cov),rtol=2e-14,atol=0)
            for method in ('pdf','logpdf'):
                call = getattr(group,method)
                for xx in ([],.5,[[1.,2.]],np.ones((2,2)),[1+2j],['x']):
                    with self.assertRaises((ValueError,TypeError)):
                        call(xx)
                for kwargs in ({'mean':[0.]},{'cov':np.eye(3)},{'cov':[1.,1.]},
                               {'allow_singular':1},{'allow_singular':np.array(True)}):
                    with self.assertRaises((ValueError,TypeError)):
                        call([0.,0.],**kwargs)
                with self.assertRaises(ValueError):
                    call(x,mean=x)
                with self.assertRaises(ValueError):
                    call(np.broadcast_to(np.zeros(1),(np.iinfo(np.intc).max+1,)))
                with self.assertRaises(TypeError):
                    call([0.,0.],out={})
            np.testing.assert_array_equal(x,[.25,-.5])
            np.testing.assert_array_equal(cov,[[2.,.125],[.125,3.]])

    def test_c_call_status_and_output_preservation(self):
        for _,_,api in self.builds:
            arrays = [np.array([.2,-.1]),np.zeros(2),np.eye(2,order='F')]
            pointers = [a.ctypes.data_as(C.c_void_p) for a in arrays]
            for method in ('pdf','logpdf'):
                function = api.multivariate_normal._functions[method]
                for n,flag in ((0,0),(-1,0),(2,7)):
                    output,status = C.c_double(123.),C.c_int(-1)
                    function(n,*pointers,flag,C.byref(output),C.byref(status))
                    self.assertEqual(status.value,1)
                    self.assertEqual(output.value,123.)
                arrays[2][0,0] = -1.
                output,status = C.c_double(123.),C.c_int(-1)
                function(2,*pointers,0,C.byref(output),C.byref(status))
                self.assertEqual(status.value,0)
                self.assertTrue(math.isnan(output.value))
                arrays[2][0,0] = 1.

    def test_native_fortran_and_c_clients(self):
        for path,metadata,api in self.builds:
            with tempfile.TemporaryDirectory() as temp:
                folder = Path(temp)
                client = folder/'native_client.f90'
                client.write_text("""program native_client
use iso_c_binding, only: c_double, c_int
use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan, ieee_positive_inf
use scifort_multivariate_normal, only: multivariate_normal_pdf, multivariate_normal_logpdf
use scifort_c_api_generated_multivariate_normal
implicit none
real(c_double) :: x(2), mean(2), cov(2,2), actual, expected
integer(c_int) :: status
integer :: i, j
logical :: allow
mean=0
cov=reshape([4.0_c_double,0.0_c_double,0.0_c_double,1.0_c_double],[2,2])
do i=1,6
    x=[0.5_c_double,-0.25_c_double]
    allow=.false.
    if (i==2 .or. i==3) then
        cov(2,2)=0
        x(2)=0
        allow=i==3
    end if
    if (i==4) x(2)=ieee_value(0.0_c_double,ieee_quiet_nan)
    if (i==5) x(1)=ieee_value(0.0_c_double,ieee_positive_inf)
    if (i==6) cov(1,1)=-1
    do j=1,2
        if (j==1) then
            expected=multivariate_normal_pdf(x,mean,cov,allow)
            call scifort_multivariate_normal_pdf_f64(2_c_int,x,mean,cov, &
                merge(1_c_int,0_c_int,allow),actual,status)
        else
            expected=multivariate_normal_logpdf(x,mean,cov,allow)
            call scifort_multivariate_normal_logpdf_f64(2_c_int,x,mean,cov, &
                merge(1_c_int,0_c_int,allow),actual,status)
        end if
        if (status/=0) stop 1
        if (ieee_is_nan(expected)) then
            if (.not. ieee_is_nan(actual)) stop 2
        else if (expected/=actual) then
            stop 3
        end if
    end do
end do
end program
""",encoding='ascii')
                env = os.environ.copy()
                env['PATH'] = str(path.parent)+os.pathsep+env.get('PATH','')
                exe = folder/'native.exe'
                flags = ['-O0','-g','-fcheck=all','-fbacktrace'] if '-fcheck=all' in metadata['flags'] else ['-O3']
                subprocess.run([metadata['compiler_path'],*flags,'-I',str(path.parent),str(client),str(path),'-o',str(exe)],check=True,capture_output=True,text=True)
                subprocess.run([str(exe)],env=env,check=True,capture_output=True,text=True)
                if shutil.which('gcc'):
                    c = folder/'client.c'
                    c.write_text("""#include "scifort.h"
#include <math.h>
int main(void) {
 double x[2]={0,0}, mean[2]={0,0}, cov[4]={1,0,0,1}, value; int status;
 scifort_multivariate_normal_pdf_f64(2,x,mean,cov,0,&value,&status);
 if (status || fabs(value-0.15915494309189535)>1e-15) return 1;
 scifort_multivariate_normal_logpdf_f64(2,x,mean,cov,0,&value,&status);
 return status || fabs(value+1.8378770664093453)>1e-14;
}
""",encoding='ascii')
                    exe = folder/'c.exe'
                    subprocess.run(['gcc','-std=c11','-Wall','-Wextra','-Werror','-I',str(ROOT/'include'),str(c),str(path),'-o',str(exe)],check=True,capture_output=True,text=True)
                    subprocess.run([str(exe)],env=env,check=True,capture_output=True,text=True)
                    if shutil.which('g++'):
                        subprocess.run(['g++','-std=c++17','-Wall','-Wextra','-Werror','-I',str(ROOT/'include'),'-fsyntax-only',str(c)],check=True,capture_output=True,text=True)
