# SPDX-License-Identifier: MIT
"""Project-independent generator and compiled arithmetic ABI checks."""
import copy
from contextlib import redirect_stderr
import io
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
                                           array_argument="t", array_style="pointer", optional_policy="required",
                                           finite=["shift", "gain"], positive=["gain"])),
                excluded=dict(total="Array arguments need a later interface pattern."))


def array_fixture_config():
    config = fixture_config()
    config.update(schema_version=2, source="array_ops.f90", module="example_arrays", wrapper_module="example_array_api",
                  header_guard="EXAMPLE_ARRAY_API_H", functions={}, excluded={},
                  outputs=dict(fortran="array_api.f90", header="array_api.h", coverage="array_coverage.json"),
                  argument_docs=dict(a="input matrix", gain="row-sum multiplier", y="scaled row sums",
                                     status="0 success; 1 invalid dimensions", nrows="row count, >= 1",
                                     ncols="column count, >= 1"),
                  subroutines=dict(row_sums=dict(symbol="example_row_sums",
                                                dimensions=[dict(name="nrows", minimum=1),
                                                            dict(name="ncols", minimum=1)],
                                                shapes=dict(a=["nrows", "ncols"], y=["nrows"]),
                                                status_argument="status", invalid_status=1,
                                                description="Scale the row sums of a real matrix.")))
    return config


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

    def test_elementwise_policy_validation_and_coverage(self):
        entry = self.config['functions']['affine']
        entry['elementwise_symbol'] = 'example_elements'
        outputs = generator.generate(self.config, self.root)
        report = json.loads(outputs[self.root / 'coverage.json'])
        self.assertEqual(report['generated_symbol_count'], 3)
        self.assertIn('example_elements', report['procedures'][0]['symbols'])
        entry['elementwise_symbol'] = entry['symbol']
        with self.assertRaisesRegex(ValueError, 'duplicate C symbol'):
            generator.generate(self.config, self.root)
        entry['elementwise_symbol'] = 'example_elements'
        source = self.root / 'arithmetic.f90'
        source.write_text(source.read_text().replace('elemental ', ''))
        with self.assertRaisesRegex(ValueError, 'elemental'):
            generator.generate(self.config, self.root)

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

    def test_bulk_only_records_existing_scalar_without_emitting_it(self):
        policy = self.config["functions"]["affine"]
        policy["existing_scalar_symbol"] = policy.pop("symbol")
        artifacts = generator.generate(self.config, self.root)
        self.assertNotIn("function example_affine(", artifacts[self.root / "api.f90"])
        self.assertNotIn("double example_affine(", artifacts[self.root / "api.h"])
        report = json.loads(artifacts[self.root / "coverage.json"])
        self.assertEqual(report["generated_symbol_count"], 1)
        row = report["procedures"][0]
        self.assertEqual(row["existing_scalar_symbol"], "example_affine")
        self.assertEqual(row["symbols"], ["example_affine_array"])
        policy["symbol"] = "other_scalar"
        with self.assertRaisesRegex(ValueError, "mutually exclusive"):
            generator.generate(self.config, self.root)
        del policy["symbol"]
        policy["existing_scalar_symbol"] = policy["vector_symbol"]
        with self.assertRaisesRegex(ValueError, "cannot also be generated"):
            generator.generate(self.config, self.root)

    def test_batch_validates_all_manifests_and_collisions(self):
        config = self.root / "config.json"
        config.write_text(json.dumps(self.config))
        batch = self.root / "batch.json"
        batch.write_text(json.dumps(["config.json"]))
        self.assertEqual(generator.collect(batch, self.root, True),
                         generator.generate(self.config, self.root))
        batch.write_text(json.dumps(["config.json", "config.json"]))
        with self.assertRaisesRegex(ValueError, "duplicate manifest"):
            generator.collect(batch, self.root, True)
        second = copy.deepcopy(self.config)
        second["outputs"] = dict(fortran="second.f90", header="second.h", coverage="second.json")
        (self.root / "second_config.json").write_text(json.dumps(second))
        batch.write_text(json.dumps(["config.json", "second_config.json"]))
        with self.assertRaisesRegex(ValueError, "duplicate generated C symbol"):
            generator.collect(batch, self.root, True)
        second["functions"]["affine"]["symbol"] = "second_affine"
        second["functions"]["affine"]["vector_symbol"] = "second_array"
        second["outputs"]["header"] = "arithmetic.f90"
        (self.root / "second_config.json").write_text(json.dumps(second))
        with self.assertRaisesRegex(ValueError, "overwrite"):
            generator.collect(batch, self.root, True)
        batch.write_text('{}')
        with self.assertRaisesRegex(ValueError, "nonempty JSON list"):
            generator.collect(batch, self.root, True)

    def test_batch_cli_is_read_only_on_drift_and_later_invalid_manifest(self):
        config = self.root / "config.json"
        config.write_text(json.dumps(self.config))
        batch = self.root / "batch.json"
        batch.write_text(json.dumps(["config.json"]))
        command = [sys.executable, HERE.parent / "generate.py", batch,
                   "--batch", "--root", self.root]
        run(command, self.root)
        run([*command, "--check"], self.root)
        header = self.root / "api.h"
        header.write_text("sentinel")
        result = subprocess.run([str(x) for x in [*command, "--check"]],
                                capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("stale", result.stderr)
        self.assertEqual(header.read_text(), "sentinel")
        (self.root / "invalid.json").write_text('{}')
        batch.write_text(json.dumps(["config.json", "invalid.json"]))
        result = subprocess.run([str(x) for x in command], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(header.read_text(), "sentinel")

    def test_explicit_shape_is_the_default_for_new_bulk_wrappers(self):
        del self.config["functions"]["affine"]["array_style"]
        legacy = generator.generate(self.config, self.root)[self.root / "api.f90"]
        self.assertIn("c_f_pointer", legacy)
        self.config["schema_version"] = 2
        text = generator.generate(self.config, self.root)[self.root / "api.f90"]
        self.assertIn("real(c_double), intent(in) :: cg_arg_0(cg_count)", text)
        self.assertIn("real(c_double), intent(inout) :: cg_output(cg_count)", text)
        self.assertNotIn("c_f_pointer", text)

    def test_assumed_shape_subroutine_gets_an_explicit_shape_adapter(self):
        shutil.copyfile(HERE / "fixtures/array_ops.f90", self.root / "array_ops.f90")
        artifacts = generator.generate(array_fixture_config(), self.root)
        text = artifacts[self.root / "array_api.f90"]
        self.assertIn("real(c_double), intent(in) :: a(nrows,ncols)", text)
        self.assertIn("real(c_double), intent(inout) :: y(nrows)", text)
        self.assertNotIn("c_f_pointer", text)
        report = json.loads(artifacts[self.root / "array_coverage.json"])
        self.assertEqual(report["generated_symbol_count"], 1)
        self.assertEqual(report["procedures"][0]["interface"], "explicit_shape_subroutine")
        self.assertEqual(report["procedures"][0]["dimension_order"], ["nrows", "ncols"])
        # JSON object key reordering must never change the C argument order.
        reordered = json.loads(json.dumps(array_fixture_config(), sort_keys=True))
        self.assertEqual(artifacts, generator.generate(reordered, self.root))

    def test_native_array_shape_comments_are_reviewed_without_rewriting(self):
        source = self.root / 'array_ops.f90'
        shutil.copyfile(HERE / 'fixtures/array_ops.f90', source)
        config = array_fixture_config()
        original = source.read_text()
        def review():
            artifacts = generator.generate(config, self.root)
            report = json.loads(artifacts[self.root / 'array_coverage.json'])
            return artifacts, report['procedures'][0]['native_shape_documentation_review']
        artifacts, issues = review()
        self.assertEqual(issues, [])
        # Bare tuples, uppercase names, spaces and a continued declaration are accepted.
        source.write_text(original.replace(':: a(:,:) !! shape (nrows,ncols);',
            ':: a( &\n            :,:) !! ( NROWS, NCOLS );'))
        self.assertEqual(review()[1], [])
        source.write_text(original.replace('shape (nrows,ncols); ', ''))
        before = source.read_bytes()
        artifacts, issues = review()
        self.assertEqual(len(issues), 1)
        self.assertEqual(issues[0]['argument'], 'a')
        self.assertIsNone(issues[0]['documented_shape'])
        self.assertEqual(issues[0]['required_shape'], ['nrows', 'ncols'])
        self.assertEqual(source.read_bytes(), before)
        output = io.StringIO()
        with redirect_stderr(output):
            generator.print_documentation_reviews(artifacts)
        self.assertIn('array_ops.f90:row_sums:a: expected shape (nrows,ncols)', output.getvalue())
        source.write_text(original.replace('shape (nrows,ncols)', 'shape (ncols,nrows)'))
        issues = review()[1]
        self.assertEqual(issues[0]['documented_shape'], ['ncols', 'nrows'])
        self.assertIn('differs', issues[0]['reason'])
        source.write_text(original.replace('shape (nrows)', 'shape (nrows,ncols)'))
        self.assertEqual(review()[1][0]['argument'], 'y')

    def test_missing_shape_comments_are_advisory_in_cli_and_check(self):
        source = self.root / 'array_ops.f90'
        shutil.copyfile(HERE / 'fixtures/array_ops.f90', source)
        source.write_text(source.read_text().replace('shape (nrows,ncols); ', ''))
        config = self.root / 'array.json'
        config.write_text(json.dumps(array_fixture_config()))
        command = [sys.executable, HERE.parent / 'generate.py', config, '--root', self.root]
        generated = run(command, self.root)
        self.assertIn('documentation review', generated.stderr)
        checked = run([*command, '--check'], self.root)
        self.assertIn('documentation review', checked.stderr)
        self.assertIn('missing or unsupported', checked.stderr)

    def test_invalid_array_shape_and_dimension_policies_are_rejected(self):
        shutil.copyfile(HERE / "fixtures/array_ops.f90", self.root / "array_ops.f90")
        for change, message in ((dict(shapes=dict(a=["nrows"], y=["nrows"])), "rank mismatch"),
                                (dict(shapes=dict(a=["nrows", "unknown"], y=["nrows"])), "dimension names"),
                                (dict(dimensions=[dict(name="a", minimum=1)]), "conflicts"),
                                (dict(status_argument="gain"), "status must"),
                                (dict(invalid_status=0), "invalid_status")):
            with self.subTest(change=change):
                config = array_fixture_config()
                config["subroutines"]["row_sums"].update(change)
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
            config["functions"]["affine"]["elementwise_symbol"] = "example_affine_elements"
            shutil.copyfile(HERE / "fixtures/arithmetic.f90", cls.root / "arithmetic.f90")
            for path, value in generator.generate(config, cls.root).items():
                path.write_text(value)
            explicit = copy.deepcopy(config)
            explicit.update(schema_version=2, wrapper_module="example_explicit_api", header_guard="EXAMPLE_EXPLICIT_API_H",
                            outputs=dict(fortran="explicit.f90", header="explicit.h", coverage="explicit.json"))
            explicit["functions"]["affine"].update(symbol="example_explicit_affine",
                                                   vector_symbol="example_explicit_affine_array")
            del explicit["functions"]["affine"]["array_style"]
            del explicit["functions"]["affine"]["elementwise_symbol"]
            for path, value in generator.generate(explicit, cls.root).items():
                path.write_text(value)
            shutil.copyfile(HERE / "fixtures/array_ops.f90", cls.root / "array_ops.f90")
            for path, value in generator.generate(array_fixture_config(), cls.root).items():
                path.write_text(value)
            sources = [cls.root / name for name in ("arithmetic.f90", "api.f90", "explicit.f90",
                                                    "array_ops.f90", "array_api.f90")]
            for profile, flags in (("checked", ["-O0", "-g", "-fcheck=all", "-fbacktrace"]),
                                   ("optimized", ["-O3"])):
                directory = cls.root / profile
                directory.mkdir()
                common = [cls.fc, *flags, "-std=f2018", "-fimplicit-none", "-Werror=implicit-interface",
                          "-J", directory, "-I", directory]
                extension = ".dll" if os.name == "nt" else ".so"
                fixture = directory / ("fixture" + extension)
                run([*common, "-shared", "-fPIC", *sources,
                     "-o", fixture], directory)
                cls.libs.append((profile, C.CDLL(str(fixture))))
                # A real C client tests the generated header and link symbols.
                client = directory / "client.c"
                client.write_text('''#include "api.h"
#include "explicit.h"
#include "array_api.h"
int main(void) {
    double x[3] = {-1.0, 0.0, 1.0}, y[3] = {99.0, 99.0, 99.0};
    int status = -1;
    example_affine_array(3, x, 2.0, 3.0, y, &status);
    if (status != 0 || y[0] != -1.0 || y[1] != 2.0 || y[2] != 5.0) return 1;
    if (example_affine(5.0, 3.0, 2.0) != 13.0) return 2;
    example_affine_array(0, NULL, 0.0, 1.0, NULL, &status);
    if (status != 0) return 3;
    example_explicit_affine_array(3, x, 2.0, 3.0, y, &status);
    if (status != 0 || y[0] != -1.0 || y[2] != 5.0) return 4;
    double a[6] = {1,2,3,4,5,6}, sums[2] = {99,99};
    example_row_sums(2, 3, a, 2.0, sums, &status);
    if (status != 0 || sums[0] != 18 || sums[1] != 24) return 5;
    double shifts[3] = {1,2,3}, gains[3] = {2,3,-1};
    example_affine_elements(3, x, shifts, gains, y, &status);
    if (status != 0 || y[0] != -1 || y[1] != 2 || y[2] != 2) return 6;
    example_affine_elements(0, x, shifts, gains, y, &status);
    if (status != 0) return 7;
    example_affine_elements((size_t)2147483648ULL, x, shifts, gains, y, &status);
    if (status != 1 || y[0] != -1) return 8;
    return 0;
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
                run([*common, *sources, obj, "-o", exe], directory)
                run([exe], directory)
                fclient = directory / "fortran_client.f90"
                fclient.write_text('''program check_adapter
    use, intrinsic :: iso_c_binding, only : c_double, c_int
    use example_arrays, only : row_sums
    use example_array_api, only : example_row_sums
    implicit none
    real(c_double) :: a(2,3), y(2), native(2)
    integer(c_int) :: status
    integer :: native_status
    a=reshape([1.0_c_double,2.0_c_double,3.0_c_double,4.0_c_double, &
        5.0_c_double,6.0_c_double],[2,3])
    call row_sums(a,2.0_c_double,native,native_status)
    call example_row_sums(2_c_int,3_c_int,a,2.0_c_double,y,status)
    if (status /= 0 .or. native_status /= 0) error stop 1
    if (any(y /= native)) error stop 2
end program check_adapter
''')
                fexe = directory / ("fortran_client.exe" if os.name == "nt" else "fortran_client")
                run([*common, *sources, fclient, "-o", fexe], directory)
                run([fexe], directory)
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

    def test_explicit_bulk_and_rectangular_subroutine(self):
        for profile, lib in self.libs:
            with self.subTest(profile=profile):
                vector = self.bind_vector(lib, "example_explicit_affine_array")
                values, output = (C.c_double * 2)(2., 5.), (C.c_double * 2)(99., 99.)
                status = C.c_int(-1)
                vector(2, values, 3., 2., output, C.byref(status))
                self.assertEqual(list(output), [7., 13.])
                self.assertEqual(status.value, 0)
                vector(0, values, 0., 1., output, C.byref(status))
                self.assertEqual(status.value, 0)
                vector(2, values, 0., 0., output, C.byref(status))
                self.assertEqual(status.value, 1)
                self.assertEqual(list(output), [7., 13.])
                row_sums = lib.example_row_sums
                row_sums.argtypes = [C.c_int, C.c_int, C.c_void_p, C.c_double, C.c_void_p,
                                     C.POINTER(C.c_int)]
                row_sums.restype = None
                a = (C.c_double * 6)(1., 2., 3., 4., 5., 6.)
                row_sums(2, 3, a, 2., output, C.byref(status))
                self.assertEqual(status.value, 0)
                self.assertEqual(list(output), [18., 24.])
                for n, m in ((0, 3), (2, 0), (-1, 3)):
                    row_sums(n, m, a, 2., output, C.byref(status))
                    self.assertEqual(status.value, 1)
                    self.assertEqual(list(output), [18., 24.])


if __name__ == "__main__":
    unittest.main()
