# SPDX-License-Identifier: MIT
"""SciFort-specific integration checks for generated normal C interfaces."""
import ctypes as C
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "tools/cabi_generator"))
import generate as generator


def run(command, cwd):
    result = subprocess.run([str(x) for x in command], cwd=cwd, capture_output=True, text=True)
    if result.returncode:
        raise AssertionError(f"{command}\n{result.stdout}\n{result.stderr}")
    return result


class GeneratedNormalTests(unittest.TestCase):
    """Compile checked and optimized shared libraries without building the full project."""
    @classmethod
    def setUpClass(cls):
        cls.fc = shutil.which(os.environ.get("FC", "gfortran"))
        cls.cc = shutil.which(os.environ.get("CC", "gcc"))
        if not cls.fc or not cls.cc:
            raise unittest.SkipTest("gfortran and gcc required for compiled ABI checks")
        cls.temp = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temp.name)
        cls.libs = []
        cls.dll_directory = os.add_dll_directory(str(Path(cls.fc).parent)) if os.name == "nt" else None
        try:
            config = json.loads((ROOT / "bindings/cabi_config/normal.json").read_text())
            for path, value in generator.generate(config, ROOT).items():
                if not path.exists() or path.read_text(encoding="utf-8") != value:
                    raise AssertionError(f"stale generated integration artifact: {path}")
            for profile, flags in (("checked", ["-O0", "-g", "-fcheck=all", "-fbacktrace"]),
                                   ("optimized", ["-O3"])):
                directory = cls.root / profile
                directory.mkdir()
                common = [cls.fc, *flags, "-std=f2018", "-fimplicit-none", "-Werror=implicit-interface",
                          "-J", directory, "-I", directory]
                extension = ".dll" if os.name == "nt" else ".so"
                normal = directory / ("normal" + extension)
                sources = [ROOT / path for path in ("src/kinds.f90", "src/constants.f90", "src/math.f90",
                                                     "src/stats/normal.f90",
                                                     "src/c_api/generated/normal_generated.f90")]
                run([*common, "-shared", "-fPIC", *sources, HERE / "fixtures/normal_oracle.f90",
                     "-o", normal], directory)
                normal_lib = C.CDLL(str(normal))
                cls.libs.append((profile, normal_lib))
                # A real C client tests the installed umbrella header and link symbols.
                client = directory / "client.c"
                client.write_text('''#include "scifort.h"
#include <math.h>
int main(void) {
    double x[3] = {-1.0, 0.0, 1.0}, y[3] = {99.0, 99.0, 99.0};
    int status = -1;
    scifort_normal_sf_vec_f64(3, x, 0.0, 1.0, y, &status);
    if (status != 0 || y[1] != 0.5) return 1;
    if (fabs(scifort_normal_logpdf_f64(0.0, 0.0, 1.0) + 0.9189385332046727) > 1e-15) return 2;
    scifort_normal_sf_vec_f64(0, NULL, 0.0, 1.0, NULL, &status);
    return status;
}
''')
                obj = directory / "client.o"
                run([cls.cc, "-std=c11", "-Wall", "-Wextra", "-Werror", "-I", ROOT / "include",
                     "-c", client, "-o", obj], directory)
                cxx = shutil.which(os.environ.get("CXX", "g++"))
                if cxx:
                    run([cxx, "-std=c++17", "-Wall", "-Wextra", "-Werror", "-fsyntax-only",
                         "-I", ROOT / "include", client], directory)
                exe = directory / ("client.exe" if os.name == "nt" else "client")
                run([*common, *sources, obj, "-o", exe], directory)
                run([exe], directory)
        except Exception:
            cls.tearDownClass()
            raise

    @classmethod
    def tearDownClass(cls):
        import _ctypes
        for _, normal in cls.libs:
            for lib in (normal,):
                if os.name == "nt":
                    _ctypes.FreeLibrary(lib._handle)
                else:
                    _ctypes.dlclose(lib._handle)
        cls.libs.clear()
        if cls.dll_directory:
            cls.dll_directory.close()
        cls.temp.cleanup()

    @staticmethod
    def bind_scalar(lib, name):
        function = getattr(lib, name)
        function.argtypes = [C.c_double] * 3
        function.restype = C.c_double
        return function

    @staticmethod
    def bind_vector(lib, name):
        function = getattr(lib, name)
        function.argtypes = [C.c_size_t, C.c_void_p, C.c_double, C.c_double,
                             C.c_void_p, C.POINTER(C.c_int)]
        function.restype = None
        return function

    def same_number(self, actual, expected):
        if math.isnan(expected):
            self.assertTrue(math.isnan(actual))
        else:
            self.assertEqual(actual, expected)

    def test_normal_scalar_bulk_native_and_reference_agreement(self):
        from statistics import NormalDist
        methods = ("isf", "logcdf", "logpdf", "logsf", "sf")
        for profile, lib in self.libs:
            oracle = lib.normal_oracle
            oracle.argtypes = [C.c_int, C.c_double, C.c_double, C.c_double]
            oracle.restype = C.c_double
            for index, method in enumerate(methods):
                with self.subTest(profile=profile, method=method):
                    scalar = self.bind_scalar(lib, f"scifort_normal_{method}_f64")
                    vector = self.bind_vector(lib, f"scifort_normal_{method}_vec_f64")
                    inputs = ([0., math.nextafter(0., 1.), 1e-300, 1e-12, .1, .5, .9,
                               1-1e-12, math.nextafter(1., 0.), 1., -.1, 1.1,
                               -math.inf, math.inf, float("nan")]
                              if method == "isf" else
                              [-math.inf, -40., -8., -1., 0., 1., 8., 40., math.inf, float("nan")])
                    x = (C.c_double * len(inputs))(*inputs)
                    y = (C.c_double * len(inputs))()
                    status = C.c_int(-1)
                    vector(len(inputs), x, .25, 1.5, y, C.byref(status))
                    self.assertEqual(status.value, 0)
                    for value, bulk in zip(inputs, y):
                        expected = oracle(index, value, .25, 1.5)
                        self.same_number(scalar(value, .25, 1.5), expected)
                        self.same_number(bulk, expected)
                    for bad_scale in (0., -1., math.inf, math.nan):
                        self.assertTrue(math.isnan(scalar(.5, .25, bad_scale)))
                        y[0] = 99.
                        vector(len(inputs), x, .25, bad_scale, y, C.byref(status))
                        self.assertEqual(status.value, 1)
                        self.assertEqual(y[0], 99.)
                    for bad_loc in (math.inf, -math.inf, math.nan):
                        self.assertTrue(math.isnan(scalar(.5, bad_loc, 1.)))
                        y[0] = 99.
                        vector(len(inputs), x, bad_loc, 1., y, C.byref(status))
                        self.assertEqual(status.value, 1)
                        self.assertEqual(y[0], 99.)
                    value = .1 if method == "isf" else .75
                    z = (value-.25)/1.5
                    expected = {"isf": lambda: .25-1.5*NormalDist().inv_cdf(value),
                                "logpdf": lambda: -.5*z*z-math.log(1.5*math.sqrt(2*math.pi)),
                                "logcdf": lambda: math.log(.5*math.erfc(-z/math.sqrt(2))),
                                "logsf": lambda: math.log(.5*math.erfc(z/math.sqrt(2))),
                                "sf": lambda: .5*math.erfc(z/math.sqrt(2))}[method]()
                    self.assertTrue(math.isclose(scalar(value, .25, 1.5), expected,
                                                 rel_tol=3e-14, abs_tol=3e-15))


if __name__ == "__main__":
    unittest.main()
