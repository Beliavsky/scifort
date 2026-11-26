# SPDX-License-Identifier: MIT
"""Rectangular matrix distributions, optional arrays and entropy dimension defaults."""
import ctypes as C
import math
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

import numpy as np
from scipy.stats import matrix_normal, t

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from build_library import build
from scifort import Scifort
from test_multivariate_t import run
ROOT=Path(__file__).resolve().parents[3]


class MatrixDensityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran'):raise unittest.SkipTest('gfortran required')
        cls.builds=[]
        for checked in (True,False):
            path,metadata=build(checked=checked)
            cls.builds.append((path,metadata,Scifort(path,dll_directories=metadata['dll_directories'])))

    def test_matrix_normal_references_and_rectangular_inputs(self):
        for _,_,api in self.builds:
            for m,n in ((1,1),(2,3),(3,2),(1,3),(3,1)):
                mean=np.arange(m*n).reshape(m,n)/5
                x=mean+np.arange(m*n).reshape(m,n)/7+.2
                row=np.eye(m)*2+.2;col=np.eye(n)*1.5+.1
                for method in ('pdf','logpdf'):
                    actual=getattr(api.matrix_normal,method)(x,mean,row,col)
                    np.testing.assert_allclose(actual,getattr(matrix_normal,method)(x,mean,row,col),rtol=8e-13,atol=4e-14)
                    self.assertIsInstance(actual,float)
                np.testing.assert_allclose(api.matrix_normal.entropy(row,col),matrix_normal.entropy(rowcov=row,colcov=col),rtol=8e-13,atol=4e-14)
                self.assertAlmostEqual(api.matrix_normal.logpdf(x,mean,row,col),api.matrix_normal.logpdf(x.T,mean.T,col,row),places=12)

    def test_matrix_t_saved_reference_reductions_and_transpose(self):
        for _,_,api in self.builds:
            group=api.matrix_t
            # Existing SciPy 1.17 reference from test_matrix_distributions_reference.f90.
            mean=np.arange(6.).reshape(3,2);x=mean+.1
            np.testing.assert_allclose(group.logpdf(x,mean,np.diag([1.,2.,3.]),np.eye(2)*.3,3.),
                                       -.00271565604466317279,rtol=0,atol=5e-14)
            np.testing.assert_allclose(group.pdf(x,mean,np.diag([1.,2.,3.]),np.eye(2)*.3,3.),
                                       .997288028013580408,rtol=0,atol=5e-14)
            # The native spread parameterization has no df factor in its quadratic form.
            for df in (.5,1.,5.,30.):
                for value in (-3.,0.,.7,20.):
                    actual=group.logpdf([[value]],[[.2]],[[2.]],[[3.]],df)
                    expected=t.logpdf(value,df,loc=.2,scale=math.sqrt(6/df))
                    np.testing.assert_allclose(actual,expected,rtol=5e-13,atol=3e-14)
            row=np.array([[2.,.3],[.3,1.]]);col=np.eye(3)*1.5+.1
            x=np.arange(6.).reshape(2,3)/3;mean=np.ones((2,3))*.2
            for df in (.2,1.,3.5,20.):
                a=group.logpdf(x,mean,row,col,df)
                b=group.logpdf(x.T,mean.T,col,row,df)
                np.testing.assert_allclose(a,b,rtol=3e-13,atol=5e-14)
                np.testing.assert_allclose(group.pdf(x,mean,row,col,df),math.exp(a),rtol=2e-15,atol=0)

    def test_optional_arrays_and_entropy_independent_dimension_fallbacks(self):
        for _,_,api in self.builds:
            x=np.arange(6.).reshape(2,3)/5
            for group,keys in ((api.matrix_normal,('rowcov','colcov')),(api.matrix_t,('row_spread','col_spread'))):
                for method in ('pdf','logpdf'):
                    fn=getattr(group,method)
                    self.assertEqual(fn(x),fn(x,np.zeros((2,3)),np.eye(2),np.eye(3)))
                    for kwargs in ({'mean':np.ones((2,3))*.2},{keys[0]:np.eye(2)*2},{keys[1]:np.eye(3)*3}):
                        expected=dict(mean=np.zeros((2,3)),**{keys[0]:np.eye(2),keys[1]:np.eye(3)})
                        expected.update(kwargs)
                        self.assertEqual(fn(x,**kwargs),fn(x,**expected))
            entropy=api.matrix_normal.entropy
            for row,col in ((None,None),(np.eye(2)*2,None),(None,np.eye(3)*3),(np.eye(2)*2,np.eye(3)*3)):
                explicit_row=np.eye(1) if row is None else row
                explicit_col=np.eye(1) if col is None else col
                self.assertEqual(entropy(row,col),entropy(explicit_row,explicit_col))
                np.testing.assert_allclose(entropy(row,col),matrix_normal.entropy(rowcov=explicit_row,colcov=explicit_col),rtol=5e-14,atol=1e-14)
            self.assertEqual(entropy(),entropy(None,None))

    def test_invalid_parameters_layouts_and_shapes(self):
        for _,_,api in self.builds:
            x=np.arange(6.).reshape(2,3)/5
            for group,keys in ((api.matrix_normal,('rowcov','colcov')),(api.matrix_t,('row_spread','col_spread'))):
                for method in ('pdf','logpdf'):
                    fn=getattr(group,method)
                    for bad in ([],[1.,2.],np.zeros((0,3)),np.ones((2,3,1))):
                        with self.assertRaises(ValueError):fn(bad)
                    for kwargs in ({'mean':np.zeros((3,2))},{keys[0]:np.eye(3)},{keys[1]:np.ones((3,2))}):
                        with self.assertRaises(ValueError):fn(x,**kwargs)
                    with self.assertRaises(TypeError):fn(x.astype(complex))
                    mean=np.zeros((2,3));readonly=x.copy();readonly.flags.writeable=False
                    storage=np.zeros((4,6));storage[::2,::2]=x
                    for value in (x.tolist(),np.array(x,order='F'),x.astype('>f8'),storage[::2,::2],readonly):
                        self.assertEqual(fn(value),fn(x))
                    with self.assertRaises(ValueError):fn(x,mean=x)
                    for bad in (np.zeros((2,2)),[[1.,2.],[2.,1.]],[[np.nan,0.],[0.,1.]],[[1.,np.inf],[0.,1.]]):
                        self.assertTrue(math.isnan(fn(x,**{keys[0]:bad})))
                    self.assertTrue(math.isnan(fn(np.full((2,3),np.inf))))
                    self.assertTrue(math.isnan(fn(x,mean=np.full((2,3),np.nan))))
                row=np.array([[2.,.3],[.3,1.]]);upper=row.copy();upper[0,1]=99.
                self.assertEqual(group.logpdf(x,**{keys[0]:row}),group.logpdf(x,**{keys[0]:upper}))
            for df in (0.,-1.,np.nan,np.inf,-np.inf):
                self.assertTrue(math.isnan(api.matrix_t.logpdf(x,df=df)))
            for row in ([],np.zeros((0,0)),np.ones((2,3))):
                with self.assertRaises(ValueError):api.matrix_normal.entropy(row)
            huge=np.broadcast_to(np.ones(1),(np.iinfo(np.intc).max+1,1))
            with self.assertRaises(ValueError):api.matrix_normal.entropy(huge)
            self.assertTrue(math.isnan(api.matrix_normal.entropy(np.zeros((2,2)))))
            self.assertEqual(api.matrix_normal.pdf(np.full((2,3),1e3)),0.)
            self.assertTrue(math.isfinite(api.matrix_normal.logpdf(np.full((2,3),1e3))))

    def test_raw_c_status_and_output_preservation(self):
        for _,_,api in self.builds:
            x=np.zeros((2,3),order='F');mean=np.zeros((2,3),order='F')
            row=np.eye(2,order='F');col=np.eye(3,order='F')
            xp,mp,rp,cp=[a.ctypes.data_as(C.c_void_p) for a in (x,mean,row,col)]
            for group in (api.matrix_normal,api.matrix_t):
                for method,fn in group._functions.items():
                    value=C.c_double(99.);status=C.c_int(-1)
                    def call(m,n):
                        args=[rp,cp] if method=='entropy' else [xp,mp,rp,cp]
                        if group is api.matrix_t:args.append(0.)
                        fn(m,n,*args,C.byref(value),C.byref(status))
                    for m,n in ((0,3),(2,0),(-1,3),(2,-1)):
                        call(m,n);self.assertEqual(status.value,1);self.assertEqual(value.value,99.)
                    if group is api.matrix_t:
                        call(2,3);self.assertEqual(status.value,0);self.assertTrue(math.isnan(value.value))

    def test_native_and_c_clients(self):
        for path,metadata,_ in self.builds:
            with tempfile.TemporaryDirectory() as temp:
                folder=Path(temp);env=os.environ.copy();env['PATH']=str(path.parent)+os.pathsep+env.get('PATH','')
                flags=['-O0','-g','-fcheck=all','-fbacktrace'] if '-fcheck=all' in metadata['flags'] else ['-O3']
                exe=folder/'native.exe'
                run([metadata['compiler_path'],*flags,'-I',path.parent,'-J',folder,
                     ROOT/'test/test_matrix_distributions_reference.f90',ROOT/'test/test_matrix_densities_c_api.f90',path,'-o',exe])
                run([exe],env=env)
                if shutil.which('gcc'):
                    c=folder/'client.c'
                    c.write_text("""#include "scifort.h"
#include <math.h>
int main(void) {
 double x[1]={0},mu[1]={0},r[1]={1},c[1]={1},value; int status;
 scifort_matrix_normal_pdf_f64(1,1,x,mu,r,c,&value,&status);
 if(status || fabs(value-1/sqrt(2*acos(-1.0)))>1e-14) return 1;
 scifort_matrix_normal_logpdf_f64(1,1,x,mu,r,c,&value,&status);
 if(status || fabs(value+log(2*acos(-1.0))/2)>1e-14) return 2;
 scifort_matrix_normal_entropy_f64(1,1,r,c,&value,&status);
 if(status || fabs(value-(1+log(2*acos(-1.0)))/2)>1e-14) return 3;
 scifort_matrix_t_pdf_f64(1,1,x,mu,r,c,1,&value,&status);
 if(status || fabs(value-1/acos(-1.0))>1e-14) return 4;
 scifort_matrix_t_logpdf_f64(1,1,x,mu,r,c,1,&value,&status);
 return status || fabs(value+log(acos(-1.0)))>1e-14;
}
""",encoding='ascii')
                    exe=folder/'c.exe'
                    run(['gcc','-std=c11','-Wall','-Wextra','-Werror','-I',ROOT/'include',c,path,'-o',exe]);run([exe],env=env)
                    if shutil.which('g++'):run(['g++','-std=c++17','-Wall','-Wextra','-Werror','-I',ROOT/'include','-fsyntax-only',c])
