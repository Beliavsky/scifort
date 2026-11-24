# SPDX-License-Identifier: MIT
"""Compare generated dense-matrix C adapters with independent native calls."""
import ctypes as C
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[3]


def run(command, cwd):
    result = subprocess.run([str(v) for v in command], cwd=cwd, capture_output=True, text=True)
    if result.returncode:
        raise AssertionError(result.stdout + result.stderr)


class LinalgAbiTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        fc, cc = shutil.which('gfortran'), shutil.which('gcc')
        if not fc or not cc:
            raise unittest.SkipTest('gfortran and gcc required')
        cls.temp = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temp.name)
        cls.libs = []
        cls.handle = os.add_dll_directory(str(Path(fc).parent)) if os.name == 'nt' else None
        oracle = cls.root / 'oracle.f90'
        oracle.write_text('''module native_oracles
use iso_c_binding
use scifort_linalg
implicit none
contains
subroutine native_cholesky(n,a,l,status) bind(c)
implicit none
integer(c_int), value :: n
real(c_double) :: a(n,n),l(n,n)
integer(c_int) :: status
call cholesky_lower(a,l,status)
end subroutine
subroutine native_eigen(n,a,w,v,status) bind(c)
implicit none
integer(c_int), value :: n
real(c_double) :: a(n,n),w(n),v(n,n)
integer(c_int) :: status
call symmetric_eigen_jacobi(a,w,v,status)
end subroutine
subroutine native_psd(n,a,w,v,rank,logdet,cutoff,status,flag) bind(c)
implicit none
integer(c_int), value :: n,flag
real(c_double) :: a(n,n),w(n),v(n,n),logdet,cutoff
integer(c_int) :: status,rank
call psd_decompose(a,w,v,rank,logdet,cutoff,status,flag == 1)
end subroutine
end module
''', encoding='ascii')
        client = cls.root / 'client.c'
        client.write_text('''#include "scifort.h"
#include <math.h>
int main(void) {
 double a[4]={4,1,99,3}, l[4], w[2], v[4], logdet, cutoff;
 double b[2]={2,3}, x[2]; int status, rank;
 scifort_cholesky_lower_f64(2,a,l,&status);
 if(status || l[0]!=2 || l[1]!=0.5 || l[2]!=0) return 1;
 scifort_symmetric_eigen_jacobi_f64(2,a,w,v,&status);
 if(status || w[0]>w[1]) return 2;
 scifort_psd_decompose_f64(2,a,w,v,&rank,&logdet,&cutoff,&status,1);
 if(status || rank!=2 || fabs(logdet-log(11.0))>1e-14) return 3;
 scifort_solve_lower_f64(2,l,b,x,&status);
 if(status || fabs(l[1]*x[0]+l[3]*x[1]-3)>1e-14) return 4;
 return 0;
}
''', encoding='ascii')
        sources = [ROOT / name for name in ('src/kinds.f90','src/linalg.f90',
                                            'src/c_api/generated/linalg_generated.f90')]
        try:
            for profile, flags in (('checked',['-O0','-g','-fcheck=all','-fbacktrace']),
                                   ('optimized',['-O3'])):
                folder = cls.root / profile
                folder.mkdir()
                common = [fc,*flags,'-std=f2018','-fimplicit-none','-Werror=implicit-interface',
                          '-J',folder,'-I',folder]
                library = folder / ('linalg.dll' if os.name == 'nt' else 'linalg.so')
                run([*common,'-shared','-fPIC',*sources,oracle,'-o',library],folder)
                cls.libs.append((profile,C.CDLL(str(library))))
                obj = folder / 'client.o'
                run([cc,'-std=c11','-Wall','-Wextra','-Werror','-I',ROOT/'include',
                     '-c',client,'-o',obj],folder)
                cxx = shutil.which('g++')
                if cxx:
                    run([cxx,'-std=c++17','-Wall','-Wextra','-Werror','-fsyntax-only',
                         '-I',ROOT/'include',client],folder)
                executable = folder / ('client.exe' if os.name == 'nt' else 'client')
                run([*common,*sources,obj,'-o',executable],folder)
                run([executable],folder)
        except Exception:
            cls.tearDownClass()
            raise

    @classmethod
    def tearDownClass(cls):
        import _ctypes
        for _, lib in cls.libs:
            (_ctypes.FreeLibrary if os.name == 'nt' else _ctypes.dlclose)(lib._handle)
        cls.libs.clear()
        if cls.handle:
            cls.handle.close()
        cls.temp.cleanup()

    @staticmethod
    def call(lib, symbol, a, n=2, flag=1):
        matrix = (C.c_double*4)(*a)
        w, v = (C.c_double*2)(99,99), (C.c_double*4)(99,99,99,99)
        rank, status = C.c_int(99), C.c_int(-1)
        logdet, cutoff = C.c_double(99), C.c_double(99)
        func = getattr(lib,symbol)
        func.restype = None
        args = [C.c_int(n),matrix]
        if 'cholesky' in symbol:
            args += [v,C.byref(status)]
        elif 'psd' in symbol:
            args += [w,v,C.byref(rank),C.byref(logdet),C.byref(cutoff),C.byref(status),C.c_int(flag)]
        else:
            args += [w,v,C.byref(status)]
        func(*args)
        return status.value,[*w,*v,rank.value,logdet.value,cutoff.value]

    def test_native_agreement_and_failure_outputs(self):
        matrices = ([4,1,99,3],[1,1,1,1],[0,0,0,0],[-1,0,0,2],
                    [math.nan,0,0,1],[math.inf,0,0,1],[2,0,0,2],[4,1,math.nan,3])
        for profile,lib in self.libs:
            for generated,native in (('cholesky_lower','cholesky'),
                                     ('symmetric_eigen_jacobi','eigen'),('psd_decompose','psd')):
                for a in matrices:
                    for flag in ((0,1) if native == 'psd' else (1,)):
                        with self.subTest(profile=profile,method=generated,matrix=a,flag=flag):
                            actual = self.call(lib,'scifort_'+generated+'_f64',a,flag=flag)
                            expected = self.call(lib,'native_'+native,a,flag=flag)
                            self.assertEqual(actual[0],expected[0])
                            for x,y in zip(actual[1],expected[1]):
                                self.assertTrue(math.isnan(x) if math.isnan(y) else x == y)
            # Minimum positive order and independently known one-dimensional results.
            status,outputs = self.call(lib,'scifort_cholesky_lower_f64',[4,0,0,0],n=1)
            self.assertEqual(status,0)
            self.assertEqual(outputs[2],2.)
            status,outputs = self.call(lib,'scifort_symmetric_eigen_jacobi_f64',[4,0,0,0],n=1)
            self.assertEqual(status,0)
            self.assertEqual(outputs[0],4.)
            self.assertEqual(abs(outputs[2]),1.)
            status,outputs = self.call(lib,'scifort_psd_decompose_f64',[4,0,0,0],n=1)
            self.assertEqual(status,0)
            self.assertEqual(outputs[-3],1)
            self.assertAlmostEqual(outputs[-2],math.log(4.),places=14)
            for flag,expected in ((0,5),(1,0)):
                result = self.call(lib,'scifort_psd_decompose_f64',[1,1,1,1],flag=flag)
                self.assertEqual(result[0],expected)
                self.assertEqual(result[1][-3],1)

    def test_adapter_errors_preserve_all_outputs(self):
        for _,lib in self.libs:
            for name in ('cholesky_lower','symmetric_eigen_jacobi','psd_decompose'):
                for n in (0,-1):
                    status,outputs = self.call(lib,'scifort_'+name+'_f64',[4,1,1,3],n=n)
                    self.assertEqual(status,1)
                    self.assertTrue(all(v == 99 for v in outputs))
            for flag in (-1,2,2147483647):
                status,outputs = self.call(lib,'scifort_psd_decompose_f64',[4,1,1,3],flag=flag)
                self.assertEqual(status,1)
                self.assertTrue(all(v == 99 for v in outputs))
