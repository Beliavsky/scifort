# SPDX-License-Identifier: MIT
"""Project-independent generator and compiled arithmetic ABI checks."""
import copy
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
sys.path.insert(0, str(HERE.parent))
import generate as generator


def fixture_config():
    return dict(schema_version=1, source="arithmetic.f90", module="example_arithmetic",
                kind_module="iso_c_binding", kind_name="c_double", wrapper_module="example_c_api",
                header_guard="EXAMPLE_C_API_H", outputs=dict(fortran="api.f90", header="api.h",
                                                              coverage="coverage.json"),
                argument_docs=dict(t="evaluation coordinate", shift="additive offset", gain="multiplier"),
                functions=dict(affine=dict(symbol="example_affine", vector_symbol="example_affine_array",
                                           array_argument="t", optional_policy="required",
                                           finite=["shift", "gain"], positive=["gain"])),
                excluded=dict(total="Array arguments need a later interface pattern."))


def run(command, cwd):
    result = subprocess.run([str(x) for x in command], cwd=cwd, capture_output=True, text=True)
    if result.returncode:
        raise AssertionError(f"{command}\n{result.stdout}\n{result.stderr}")
    return result


class GeneratorTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        shutil.copyfile(HERE / "fixtures/arithmetic.f90", self.root / "arithmetic.f90")
        self.config = fixture_config()

    def test_independent_discovery_and_deterministic_coverage(self):
        artifacts = generator.generate(self.config, self.root)
        self.assertEqual(artifacts, generator.generate(self.config, self.root))
        report = json.loads(artifacts[self.root / "coverage.json"])
        self.assertEqual(report["public_api_count"], 2)
        self.assertEqual(report["generated_symbol_count"], 2)
        self.assertEqual(report["procedures"][0]["optional_promoted_to_required"], ["shift", "gain"])
        self.assertEqual(report["procedures"][1]["state"], "excluded")
        self.assertNotIn("scifort", artifacts[self.root / "api.f90"].lower())

    def test_unaccounted_public_api_is_an_error(self):
        self.config["excluded"].clear()
        with self.assertRaisesRegex(ValueError, "unaccounted=.*total"):
            generator.generate(self.config, self.root)

    def test_line_endings_do_not_change_generated_artifacts(self):
        source = self.root / "arithmetic.f90"
        original = generator.generate(self.config, self.root)
        source.write_bytes(source.read_bytes().replace(b"\r\n", b"\n").replace(b"\n", b"\r\n"))
        self.assertEqual(original, generator.generate(self.config, self.root))

    def test_public_type_and_generic_require_explicit_exclusions(self):
        source = self.root / "arithmetic.f90"
        source.write_text(source.read_text().replace("contains\n", """type, public :: box
        real(c_double) :: value
    end type box
    interface evaluate
        module procedure affine
    end interface evaluate
    public :: evaluate
contains
""", 1))
        with self.assertRaisesRegex(ValueError, "unaccounted=.*box.*evaluate"):
            generator.generate(self.config, self.root)
        self.config["excluded"].update(box="Opaque handles not implemented.",
                                       evaluate="Wrap the specific procedure instead.")
        report = generator.generate(self.config, self.root)[self.root / "coverage.json"]
        self.assertEqual(json.loads(report)["public_api_count"], 4)

    def test_default_public_is_rejected_instead_of_missing_reexports(self):
        source = self.root / "arithmetic.f90"
        source.write_text(source.read_text().replace("    private\n", ""))
        with self.assertRaisesRegex(ValueError, "private-by-default"):
            generator.generate(self.config, self.root)

    def test_optionals_need_an_explicit_policy(self):
        del self.config["functions"]["affine"]["optional_policy"]
        with self.assertRaisesRegex(ValueError, "optional_policy"):
            generator.generate(self.config, self.root)

    def test_unsupported_array_cannot_generate_a_scalar_wrapper(self):
        self.config["excluded"].clear()
        self.config["functions"]["total"] = dict(symbol="example_total")
        with self.assertRaisesRegex(ValueError, "scalar intent"):
            generator.generate(self.config, self.root)

    def test_invalid_policy_and_duplicate_symbols_are_rejected(self):
        for change, message in ((dict(finte=["gain"]), "policy fields"),
                                (dict(vector_symbol="example_affine"), "duplicate C symbol"),
                                (dict(positive=["t"]), "uniform parameter"),
                                (dict(symbol="invalid-name"), "identifier")):
            with self.subTest(change=change):
                config = copy.deepcopy(self.config)
                config["functions"]["affine"].update(change)
                with self.assertRaisesRegex(ValueError, message):
                    generator.generate(config, self.root)

    def test_cli_check_detects_stale_files_without_rewriting(self):
        config = self.root / "config.json"
        config.write_text(json.dumps(self.config))
        command = [sys.executable, HERE.parent / "generate.py", config, "--root", self.root]
        run(command, self.root)
        run([*command, "--check"], self.root)
        output = self.root / "api.h"
        output.write_text(output.read_text() + "/* deliberate drift */\n")
        result = subprocess.run([str(x) for x in [*command, "--check"]], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("stale", result.stderr)
        self.assertIn("deliberate drift", output.read_text())
        # A bad configuration is rejected before any artifact is overwritten.
        self.config["functions"]["affine"]["array_argument"] = "unknown"
        config.write_text(json.dumps(self.config))
        result = subprocess.run([str(x) for x in command], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("deliberate drift", output.read_text())


class CompiledTests(unittest.TestCase):
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
            config = fixture_config()
            shutil.copyfile(HERE / "fixtures/arithmetic.f90", cls.root / "arithmetic.f90")
            for path, value in generator.generate(config, cls.root).items():
                path.write_text(value)
            for profile, flags in (("checked", ["-O0", "-g", "-fcheck=all", "-fbacktrace"]),
                                   ("optimized", ["-O3"])):
                directory = cls.root / profile
                directory.mkdir()
                common = [cls.fc, *flags, "-std=f2018", "-fimplicit-none", "-Werror=implicit-interface",
                          "-J", directory, "-I", directory]
                extension = ".dll" if os.name == "nt" else ".so"
                fixture = directory / ("fixture" + extension)
                run([*common, "-shared", "-fPIC", cls.root / "arithmetic.f90", cls.root / "api.f90",
                     "-o", fixture], directory)
                cls.libs.append((profile, C.CDLL(str(fixture))))
                # A real C client tests the generated header and link symbols.
                client = directory / "client.c"
                client.write_text('''#include "api.h"
int main(void) {
    double x[3] = {-1.0, 0.0, 1.0}, y[3] = {99.0, 99.0, 99.0};
    int status = -1;
    example_affine_array(3, x, 2.0, 3.0, y, &status);
    if (status != 0 || y[0] != -1.0 || y[1] != 2.0 || y[2] != 5.0) return 1;
    if (example_affine(5.0, 3.0, 2.0) != 13.0) return 2;
    example_affine_array(0, NULL, 0.0, 1.0, NULL, &status);
    return status;
}
''')
                obj = directory / "client.o"
                run([cls.cc, "-std=c11", "-Wall", "-Wextra", "-Werror", "-I", cls.root,
                     "-c", client, "-o", obj], directory)
                cxx = shutil.which(os.environ.get("CXX", "g++"))
                if cxx:
                    run([cxx, "-std=c++17", "-Wall", "-Wextra", "-Werror", "-fsyntax-only",
                         "-I", cls.root, client], directory)
                exe = directory / ("client.exe" if os.name == "nt" else "client")
                run([*common, cls.root / "arithmetic.f90", cls.root / "api.f90", obj, "-o", exe], directory)
                run([exe], directory)
        except Exception:
            cls.tearDownClass()
            raise

    @classmethod
    def tearDownClass(cls):
        import _ctypes
        for _, fixture in cls.libs:
            for lib in (fixture,):
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

    def test_independent_fixture_scalar_and_bulk(self):
        for profile, lib in self.libs:
            with self.subTest(profile=profile):
                scalar = self.bind_scalar(lib, "example_affine")
                vector = self.bind_vector(lib, "example_affine_array")
                self.assertEqual(scalar(5., 3., 2.), 13.)
                x = (C.c_double * 4)(-2., 0., 7., float("nan"))
                y = (C.c_double * 4)()
                status = C.c_int(-1)
                vector(4, x, 3., 2., y, C.byref(status))
                self.assertEqual(status.value, 0)
                for actual, expected in zip(y, [-1., 3., 17., float("nan")]):
                    self.same_number(actual, expected)

    def test_invalid_bulk_calls_preserve_output_and_empty_calls_allow_null(self):
        for profile, lib in self.libs:
            vector = self.bind_vector(lib, "example_affine_array")
            for n, x, shift, gain, out, expected_status in (
                    (0, None, 0., 1., None, 0), (0, None, float("nan"), 1., None, 1),
                    (1, None, 0., 1., True, 1), (1, True, 0., 1., None, 1),
                    (1, True, 0., 0., True, 1), (1, True, 0., float("inf"), True, 1),
                    (2**63-1, True, 0., 1., True, 1)):
                with self.subTest(profile=profile, n=n, shift=shift, gain=gain):
                    values = (C.c_double * 1)(2.)
                    output = (C.c_double * 1)(99.)
                    status = C.c_int(-1)
                    vector(n, values if x else None, shift, gain,
                           output if out else None, C.byref(status))
                    self.assertEqual(status.value, expected_status)
                    self.assertEqual(output[0], 99.)


if __name__ == "__main__":
    unittest.main()
