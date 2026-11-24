# SPDX-License-Identifier: MIT
"""Checked/optimized ABI agreement with the actual native matrix-helper module."""
import ctypes as C
import math
import os
from pathlib import Path
import shutil
import tempfile
import unittest
from test_generated_linalg_decomposition import run

ROOT = Path(__file__).resolve().parents[3]
METHODS = ('spd_factor','invert_lower','solve_lower_matrix','solve_spd_matrix')


class MatrixHelperTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        fc, cc = shutil.which('gfortran'),shutil.which('gcc')
        if not fc or not cc:
            raise unittest.SkipTest('gfortran and gcc required')
        cls.temp = tempfile.TemporaryDirectory()
        cls.folder = Path(cls.temp.name)
        cls.libs = []
        cls.handle = os.add_dll_directory(str(Path(fc).parent)) if os.name=='nt' else None
        sources = [ROOT/p for p in ('src/kinds.f90','src/constants.f90','src/math.f90','src/linalg.f90',
            'src/stats/normal.f90','src/special/digamma.f90','src/special/log_gamma.f90',
            'src/special/elementary.f90','src/stats/matrix_distribution_helpers.f90',
            'src/c_api/generated/matrix_helpers_generated.f90')]
        oracle = ['module matrix_oracles','use iso_c_binding','use scifort_matrix_distribution_helpers',
                  'implicit none','contains']
        for name in METHODS:
            solve = name.startswith('solve_')
            args = 'n,k,l,b,x,status' if solve else ('n,a,l,logdet,status' if name=='spd_factor'
                                                   else 'n,l,linv,status')
            oracle += [f'subroutine native_{name}({args}) bind(c)','implicit none',
                       'integer(c_int), value :: '+('n,k' if solve else 'n'),
                       'integer(c_int) :: status']
            if solve:
                oracle += ['real(c_double) :: l(n,n),b(n,k),x(n,k)',f'call {name}(l,b,x,status)']
            elif name=='spd_factor':
                oracle += ['real(c_double) :: a(n,n),l(n,n),logdet','call spd_factor(a,l,logdet,status)']
            else:
                oracle += ['real(c_double) :: l(n,n),linv(n,n)','call invert_lower(l,linv,status)']
            oracle += ['end subroutine']
        oracle += ['end module']
        oracle_path = cls.folder/'oracle.f90'
        oracle_path.write_text('\n'.join(oracle),encoding='ascii')
        client = cls.folder/'client.c'
        client.write_text('''#include "scifort.h"
#include <math.h>
int main(void) {
 double a[4]={4,2,2,10}, l[4], inverse[4], logdet;
 double b[6]={2,3,-1,4,7,-2}, x[6]; int status;
 scifort_spd_factor_f64(2,a,l,&logdet,&status);
 if(status || l[0]!=2 || l[1]!=1 || l[2]!=0 || l[3]!=3) return 1;
 if(fabs(logdet-log(36.0))>1e-14) return 2;
 scifort_invert_lower_f64(2,l,inverse,&status);
 if(status || inverse[0]!=0.5 || inverse[2]!=0) return 3;
 scifort_solve_lower_matrix_f64(2,3,l,b,x,&status);
 if(status) return 4;
 for(int j=0;j<3;++j) {
  if(fabs(2*x[2*j]-b[2*j])>1e-14) return 5;
  if(fabs(x[2*j]+3*x[2*j+1]-b[2*j+1])>1e-14) return 6;
 }
 scifort_solve_spd_matrix_f64(2,3,l,b,x,&status);
 if(status) return 7;
 for(int j=0;j<3;++j) {
  if(fabs(4*x[2*j]+2*x[2*j+1]-b[2*j])>1e-14) return 8;
  if(fabs(2*x[2*j]+10*x[2*j+1]-b[2*j+1])>1e-14) return 9;
 }
 return 0;
}
''',encoding='ascii')
        try:
            for profile,flags in (('checked',['-O0','-g','-fcheck=all','-fbacktrace']),
                                   ('optimized',['-O3'])):
                folder = cls.folder/profile
                folder.mkdir()
                common = [fc,*flags,'-std=f2018','-fimplicit-none','-Werror=implicit-interface',
                          '-J',folder,'-I',folder]
                lib = folder/('helpers.dll' if os.name=='nt' else 'helpers.so')
                run([*common,'-shared','-fPIC',*sources,oracle_path,'-o',lib],folder)
                cls.libs.append((profile,C.CDLL(str(lib))))
                obj = folder/'client.o'
                run([cc,'-std=c11','-Wall','-Wextra','-Werror','-I',ROOT/'include',
                     '-c',client,'-o',obj],folder)
                if shutil.which('g++'):
                    run(['g++','-std=c++17','-Wall','-Wextra','-Werror','-fsyntax-only',
                         '-I',ROOT/'include',client],folder)
                exe = folder/('client.exe' if os.name=='nt' else 'client')
                run([*common,*sources,obj,'-o',exe],folder)
                run([exe],folder)
        except Exception:
            cls.tearDownClass()
            raise

    @classmethod
    def tearDownClass(cls):
        import _ctypes
        for _,lib in cls.libs:
            (_ctypes.FreeLibrary if os.name=='nt' else _ctypes.dlclose)(lib._handle)
        cls.libs.clear()
        if cls.handle:
            cls.handle.close()
        cls.temp.cleanup()

    @staticmethod
    def call(lib,name,a,n=2,k=3,native=False):
        func = getattr(lib,('native_' if native else 'scifort_')+name+('' if native else '_f64'))
        func.restype = None
        inp = (C.c_double*4)(*a)
        out = (C.c_double*6)(*[99.]*6)
        status = C.c_int(-1)
        logdet = C.c_double(99)
        if name.startswith('solve_'):
            func.argtypes = [C.c_int,C.c_int,C.c_void_p,C.c_void_p,C.c_void_p,C.POINTER(C.c_int)]
            rhs = (C.c_double*6)(2,3,-1,4,7,-2)
            func(n,k,inp,rhs,out,C.byref(status))
        elif name=='spd_factor':
            func.argtypes = [C.c_int,C.c_void_p,C.c_void_p,C.POINTER(C.c_double),C.POINTER(C.c_int)]
            func(n,inp,out,C.byref(logdet),C.byref(status))
        else:
            func.argtypes = [C.c_int,C.c_void_p,C.c_void_p,C.POINTER(C.c_int)]
            func(n,inp,out,C.byref(status))
        return status.value,[*out,logdet.value]

    def test_native_success_and_failure_output_agreement(self):
        for profile,lib in self.libs:
            for name in METHODS:
                for a in ([4,2,2,10],[1,1,1,1],[-1,0,0,2],[0,0,0,1],
                          [math.nan,0,0,1],[math.inf,0,0,1],[4,2,math.nan,10]):
                    for n in (1,2):
                        with self.subTest(profile=profile,name=name,matrix=a,n=n):
                            actual = self.call(lib,name,a,n=n)
                            expected = self.call(lib,name,a,n=n,native=True)
                            self.assertEqual(actual[0],expected[0])
                            for x,y in zip(actual[1],expected[1]):
                                self.assertTrue(math.isnan(x) if math.isnan(y) else x==y)
            # Independent expected status distinctions: nonfinite versus non-SPD.
            self.assertEqual(self.call(lib,'spd_factor',[4,2,math.nan,10])[0],2)
            self.assertEqual(self.call(lib,'spd_factor',[1,1,1,1])[0],3)

    def test_invalid_extents_preserve_output(self):
        for _,lib in self.libs:
            for name in METHODS:
                cases = [(0,3),(-1,3)] + ([(2,-1)] if name.startswith('solve_') else [])
                for n,k in cases:
                    status,out = self.call(lib,name,[4,2,2,10],n=n,k=k)
                    self.assertEqual(status,1)
                    self.assertEqual(out,[99.]*7)

    def test_empty_rhs_skips_diagonals(self):
        for _,lib in self.libs:
            for name in ('solve_lower_matrix','solve_spd_matrix'):
                # Valid non-NULL scratch buffers are supplied for the zero-size arrays.
                for a in ([0,0,0,0],[math.nan,0,0,math.inf]):
                    self.assertEqual(self.call(lib,name,a,k=0),(0,[99.]*7))
                    self.assertEqual(self.call(lib,name,a,k=0,native=True),(0,[99.]*7))
