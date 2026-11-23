# SPDX-License-Identifier: MIT
"""C and Fortran integration checks for the generated explicit-shape solver."""
import ctypes as C
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools/cabi_generator"))
import generate as generator


def run(command, cwd):
    result = subprocess.run([str(x) for x in command], cwd=cwd, capture_output=True, text=True)
    if result.returncode:
        raise AssertionError(f"{command}\n{result.stdout}\n{result.stderr}")


class GeneratedLinalgTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        fc = shutil.which(os.environ.get("FC", "gfortran"))
        cc = shutil.which(os.environ.get("CC", "gcc"))
        if not fc or not cc:
            raise unittest.SkipTest("gfortran and gcc required for compiled ABI checks")
        config = json.loads((ROOT / "bindings/cabi_config/linalg.json").read_text())
        for path, value in generator.generate(config, ROOT).items():
            if not path.exists() or path.read_text(encoding="utf-8") != value:
                raise AssertionError(f"stale generated integration artifact: {path}")
        cls.temp = tempfile.TemporaryDirectory()
        cls.libs = []
        cls.dll_directory = os.add_dll_directory(str(Path(fc).parent)) if os.name == "nt" else None
        try:
            for profile, flags in (("checked", ["-O0", "-g", "-fcheck=all", "-fbacktrace"]),
                                   ("optimized", ["-O3"])):
                directory = Path(cls.temp.name) / profile
                directory.mkdir()
                common = [fc, *flags, "-std=f2018", "-fimplicit-none", "-Werror=implicit-interface",
                          "-J", directory, "-I", directory]
                sources = [ROOT / name for name in ("src/kinds.f90", "src/linalg.f90",
                                                     "src/c_api/generated/linalg_generated.f90")]
                library = directory / ("solver.dll" if os.name == "nt" else "solver.so")
                run([*common, "-shared", "-fPIC", *sources, "-o", library], directory)
                cls.libs.append((profile, C.CDLL(str(library))))
                client = directory / "client.c"
                client.write_text('''#include "scifort.h"
int main(void) {
    double l[9] = {2,1,-1, 99,3,2, -55,8,4};
    double b[3] = {2,-5,7}, x[3] = {99,99,99};
    int status = -1;
    scifort_solve_lower_f64(3, l, b, x, &status);
    if (status != 0 || x[0] != 1 || x[1] != -2 || x[2] != 3) return 1;
    scifort_solve_lower_f64(0, l, b, x, &status);
    if (status != 1 || x[0] != 1 || x[2] != 3) return 2;
    return 0;
}
''')
                obj = directory / "client.o"
                run([cc, "-std=c11", "-Wall", "-Wextra", "-Werror", "-I", ROOT / "include",
                     "-c", client, "-o", obj], directory)
                cxx = shutil.which(os.environ.get("CXX", "g++"))
                if cxx:
                    run([cxx, "-std=c++17", "-Wall", "-Wextra", "-Werror", "-fsyntax-only",
                         "-I", ROOT / "include", client], directory)
                exe = directory / ("client.exe" if os.name == "nt" else "client")
                run([*common, *sources, obj, "-o", exe], directory)
                run([exe], directory)
                fclient = directory / "native_agreement.f90"
                fclient.write_text('''program native_agreement
    use, intrinsic :: iso_c_binding, only : c_int
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_linalg, only : solve_lower
    use scifort_c_api_generated_linalg, only : scifort_solve_lower_f64
    implicit none
    real(dp) :: l(3,3), b(3), native(3), wrapped(3)
    integer :: native_status, scenario, i
    integer(c_int) :: status
    do scenario=1,3
        l=reshape([2.0_dp,1.0_dp,-1.0_dp,99.0_dp,3.0_dp,2.0_dp, &
            -55.0_dp,8.0_dp,4.0_dp],[3,3])
        b=[2.0_dp,-5.0_dp,7.0_dp]
        if (scenario == 2) l(3,3)=0.0_dp
        if (scenario == 3) b(2)=ieee_value(0.0_dp,ieee_quiet_nan)
        call solve_lower(l,b,native,native_status)
        call scifort_solve_lower_f64(3_c_int,l,b,wrapped,status)
        if (status /= native_status) error stop 1
        do i=1,3
            if (ieee_is_nan(native(i))) then
                if (.not. ieee_is_nan(wrapped(i))) error stop 2
            else
                if (native(i) /= wrapped(i)) error stop 3
            end if
        end do
    end do
end program native_agreement
''')
                fexe = directory / ("native_agreement.exe" if os.name == "nt" else "native_agreement")
                run([*common, *sources, fclient, "-o", fexe], directory)
                run([fexe], directory)
        except Exception:
            cls.tearDownClass()
            raise

    @classmethod
    def tearDownClass(cls):
        import _ctypes
        for _, lib in cls.libs:
            if os.name == "nt":
                _ctypes.FreeLibrary(lib._handle)
            else:
                _ctypes.dlclose(lib._handle)
        cls.libs.clear()
        if cls.dll_directory:
            cls.dll_directory.close()
        cls.temp.cleanup()

    def test_success_singular_and_invalid_orders(self):
        for profile, lib in self.libs:
            with self.subTest(profile=profile):
                solve = lib.scifort_solve_lower_f64
                solve.argtypes = [C.c_int, C.c_void_p, C.c_void_p, C.c_void_p, C.POINTER(C.c_int)]
                solve.restype = None
                matrix = (C.c_double * 9)(2, 1, -1, 99, 3, 2, -55, 8, 4)
                rhs = (C.c_double * 3)(2, -5, 7)
                result = (C.c_double * 3)(99, 99, 99)
                status = C.c_int(-1)
                solve(3, matrix, rhs, result, C.byref(status))
                self.assertEqual(status.value, 0)
                self.assertEqual(list(result), [1., -2., 3.])
                for n in (0, -1):
                    solve(n, matrix, rhs, result, C.byref(status))
                    self.assertEqual(status.value, 1)
                    self.assertEqual(list(result), [1., -2., 3.])
                matrix[8] = 0.
                solve(3, matrix, rhs, result, C.byref(status))
                self.assertEqual(status.value, 5)
                self.assertEqual(list(result), [1., -2., 0.])


if __name__ == "__main__":
    unittest.main()
