# SPDX-License-Identifier: MIT
"""Compile all seven benchmark families and compare bulk, scalar and native calls."""
import ctypes as C
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / 'tools/cabi_generator'))
import generate
from test_generated_normal import run

FAMILIES = ('normal', 'exponential', 'logistic', 'laplace', 'cauchy', 'weibull', 'lognormal')
METHODS = ('pdf', 'logpdf', 'cdf', 'sf', 'logcdf', 'logsf', 'ppf', 'isf')


class DistributionAbiTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fc = shutil.which(os.environ.get('FC', 'gfortran'))
        cls.cc = shutil.which(os.environ.get('CC', 'gcc'))
        if not cls.fc or not cls.cc:
            raise unittest.SkipTest('gfortran and gcc required')
        cls.temp = tempfile.TemporaryDirectory()
        cls.directory = Path(cls.temp.name)
        cls.libs = []
        cls.dll_dir = os.add_dll_directory(str(Path(cls.fc).parent)) if os.name == 'nt' else None
        try:
            # Compile unchanged production scalar wrappers, extracted verbatim to avoid
            # dependencies on the other hundred distribution families in api.f90.
            source = (ROOT / 'src/c_api/api.f90').read_text(encoding='utf-8')
            manual = ['module production_subset', 'use iso_c_binding',
                      'use ieee_arithmetic, only: ieee_is_finite', 'use scifort_kinds, only: dp']
            manual += [f'use scifort_{family}' for family in FAMILIES]
            manual += ['implicit none', 'integer(c_int), parameter :: status_ok=0, status_invalid_argument=1',
                       'contains']
            for family in FAMILIES:
                for method in ('pdf', 'cdf', 'ppf'):
                    name = f'scifort_{family}_{method}_f64'
                    match = re.search(rf'^    function {name}\(.*?^    end function {name}',
                                      source, re.M | re.S)
                    if not match:
                        raise AssertionError(f'missing production scalar {name}')
                    manual.append(match.group())
            manual.append('end module production_subset')
            manual_path = cls.directory / 'production_subset.f90'
            manual_path.write_text('\n'.join(manual), encoding='ascii')
            # Independent direct native calls do not use generated conversion logic.
            oracle = ['module native_reference', 'use iso_c_binding',
                      'use scifort_kinds, only: dp']
            oracle += [f'use scifort_{family}' for family in FAMILIES]
            oracle += ['implicit none', 'contains']
            for family in FAMILIES:
                shape = family in ('weibull', 'lognormal')
                args = 'v,shape,loc,scale' if shape else 'v,loc,scale'
                for method in METHODS:
                    name = f'native_{family}_{method}'
                    oracle += [f'function {name}({args}) result(y) bind(c)', 'implicit none',
                               f'real(c_double), value :: {args}', 'real(c_double) :: y',
                               f'y={family}_{method}({args})', f'end function {name}']
            oracle += ['end module native_reference']
            oracle_path = cls.directory / 'native_reference.f90'
            oracle_path.write_text('\n'.join(oracle), encoding='ascii')
            # C clients link every newly exposed symbol, including bulk-only entries.
            client = ['#include "scifort.h"', '#include <math.h>', 'int main(void) {',
                      'double x[1]={0.3}, y[1]; int status=-1;']
            for family in FAMILIES[1:]:
                args = '1.7,0.25,1.5' if family in ('weibull', 'lognormal') else '0.25,1.5'
                for method in METHODS:
                    name = f'scifort_{family}_{method}'
                    client += [f'{name}_vec_f64(1,x,{args},y,&status);',
                               f'if(status || y[0]!={name}_f64(x[0],{args})) return 1;']
            client += ['return 0;', '}']
            client_path = cls.directory / 'client.c'
            client_path.write_text('\n'.join(client), encoding='ascii')
            artifacts = generate.collect(ROOT / 'bindings/cabi_config/manifests.json', ROOT, True)
            for path, text in artifacts.items():
                if path.read_text(encoding='utf-8') != text:
                    raise AssertionError(f'stale artifact: {path}')
            sources = [ROOT / p for p in ('src/kinds.f90', 'src/constants.f90', 'src/math.f90')]
            sources += [ROOT / f'src/stats/{f}.f90' for f in FAMILIES]
            sources += [ROOT / f'src/c_api/generated/{f}_generated.f90' for f in FAMILIES]
            sources += [manual_path, oracle_path]
            for profile, flags in (('checked', ['-O0', '-g', '-fcheck=all', '-fbacktrace']),
                                   ('optimized', ['-O3'])):
                directory = cls.directory / profile
                directory.mkdir()
                common = [cls.fc, *flags, '-std=f2018', '-fimplicit-none', '-Werror=implicit-interface',
                          '-J', directory, '-I', directory]
                library = directory / ('families.dll' if os.name == 'nt' else 'families.so')
                run([*common, '-shared', '-fPIC', *sources, '-o', library], directory)
                cls.libs.append((profile, C.CDLL(str(library))))
                obj = directory / 'client.o'
                run([cls.cc, '-std=c11', '-Wall', '-Wextra', '-Werror', '-I', ROOT / 'include',
                     '-c', client_path, '-o', obj], directory)
                cxx = shutil.which('g++')
                if cxx:
                    run([cxx, '-std=c++17', '-Wall', '-Wextra', '-Werror', '-fsyntax-only',
                         '-I', ROOT / 'include', client_path], directory)
                exe = directory / ('client.exe' if os.name == 'nt' else 'client')
                run([*common, *sources, obj, '-o', exe], directory)
                run([exe], directory)
        except Exception:
            cls.tearDownClass()
            raise

    @classmethod
    def tearDownClass(cls):
        import _ctypes
        for _, lib in cls.libs:
            (_ctypes.FreeLibrary if os.name == 'nt' else _ctypes.dlclose)(lib._handle)
        cls.libs.clear()
        if cls.dll_dir:
            cls.dll_dir.close()
        cls.temp.cleanup()

    def test_native_scalar_and_bulk_agreement(self):
        for profile, lib in self.libs:
            for family in FAMILIES:
                params = [1.7, .25, 1.5] if family in ('weibull', 'lognormal') else [.25, 1.5]
                for method in METHODS:
                    # Original normal PDF/CDF/PPF bulk wrappers are covered elsewhere;
                    # this shared library contains only their production scalar functions.
                    if family == 'normal' and method in ('pdf', 'cdf', 'ppf'):
                        continue
                    with self.subTest(profile=profile, family=family, method=method):
                        scalar = getattr(lib, f'scifort_{family}_{method}_f64')
                        native = getattr(lib, f'native_{family}_{method}')
                        for func in (scalar, native):
                            func.argtypes = [C.c_double] * (1+len(params))
                            func.restype = C.c_double
                        bulk = getattr(lib, f'scifort_{family}_{method}_vec_f64')
                        bulk.argtypes = [C.c_size_t, C.c_void_p, *([C.c_double]*len(params)),
                                         C.c_void_p, C.POINTER(C.c_int)]
                        bulk.restype = None
                        values = ([0., math.nextafter(0., 1.), 1e-300, 1e-12, .1, .5, .9,
                                   math.nextafter(1., 0.), 1., -.1, 1.1, math.nan]
                                  if method in ('ppf', 'isf') else
                                  [-math.inf, -1000., -8., 0., .25, 1., 8., 1000., math.inf, math.nan])
                        x = (C.c_double * len(values))(*values)
                        y = (C.c_double * len(values))()
                        status = C.c_int(-1)
                        bulk(len(values), x, *params, y, C.byref(status))
                        self.assertEqual(status.value, 0)
                        for value, actual in zip(values, y):
                            expected = native(value, *params)
                            for result in (actual, scalar(value, *params)):
                                if math.isnan(expected):
                                    self.assertTrue(math.isnan(result))
                                else:
                                    self.assertEqual(result, expected)
                        # All uniform parameters, including shapes, have explicit policies.
                        for index in range(len(params)):
                            bad_values = [math.inf, -math.inf, math.nan]
                            if index != len(params)-2:  # location need not be positive
                                bad_values += [0., -1.]
                            for bad in bad_values:
                                invalid = list(params)
                                invalid[index] = bad
                                for i in range(len(y)):
                                    y[i] = 99.
                                bulk(len(values), x, *invalid, y, C.byref(status))
                                self.assertEqual(status.value, 1)
                                self.assertEqual(list(y), [99.]*len(y))
                                bulk(0, x, *invalid, y, C.byref(status))
                                self.assertEqual(status.value, 1)
                        bulk(0, x, *params, y, C.byref(status))
                        self.assertEqual(status.value, 0)
                        bulk(C.c_size_t(-1).value, x, *params, y, C.byref(status))
                        self.assertEqual(status.value, 1)
                        self.assertEqual(list(y), [99.]*len(y))


if __name__ == '__main__':
    unittest.main()
