#!/usr/bin/env python3
"""Build user-supplied perazz/FITPACK and SciFort optional adapter OFFLINE.

Requires Python >=3.11, GNU Fortran, ar, and (unless --no-cpp) g++.
For independent numerical verification, SciPy and NumPy are also needed.

Example:
    python tools/interpolate_reference/check_optional_fitpack_local.py fitpack-main.zip

This does not modify SciFort's FPM manifest, vendor FITPACK, or fetch anything.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import tomllib
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[2]
CORE_SOURCES = (
    'fitpack_core.F90', 'fitpack_fitters.f90', 'fitpack_curves.f90',
    'fitpack_parametric.f90', 'fitpack_gridded_splines.f90',
    'fitpack_gridded_surfaces.f90', 'fitpack_surfaces.f90',
    'fitpack_parametric_surfaces.f90', 'fitpack_convex_curves.f90',
    'fitpack_polar.f90', 'fitpack_spheres.f90', 'fitpack_gridded_polar.f90',
    'fitpack_gridded_sphere.f90', 'fitpack.f90',
)
TEST_MODULES = (
    'fitpack_test_data.f90', 'fitpack_curve_tests.f90',
    'fitpack_grid_nd_tests.f90', 'fitpack_tests.f90', 'fitpack_cpp_tests.f90',
)
SCIFORT_SOURCES = (
    ROOT / 'src/kinds.f90', ROOT / 'src/interpolate_bspline.f90',
    ROOT / 'src/interpolate_fitpack_bivariate.f90',
    ROOT / 'src/interpolate_curfit_core.f90', ROOT / 'src/interpolate_fitpack_compat.f90',
    ROOT / 'optional/fitpack_backend/src/scifort_optional_fitpack.f90',
)


def invoke(command: list[str | Path], cwd: Path, name: str) -> None:
    args = [str(a) for a in command]
    p = subprocess.run(args, cwd=cwd, capture_output=True, text=True)
    (cwd / f'{name}.log').write_text('COMMAND: ' + ' '.join(args) + '\n\n' + p.stdout + p.stderr)
    if p.returncode:
        raise RuntimeError(f'{name} failed with exit code {p.returncode}:\n'
                           f'{p.stdout[-1400:]}\n{p.stderr[-3600:]}')


def unpack_or_resolve(source: Path, build: Path) -> Path:
    if source.is_dir():
        candidates = [source, *(f.parent for f in source.rglob('fpm.toml'))]
    else:
        extracted = build / 'supplied_fitpack'
        extracted.mkdir(parents=True, exist_ok=True)
        with ZipFile(source) as z:
            for info in z.infolist():
                destination = (extracted / info.filename).resolve()
                if not destination.is_relative_to(extracted.resolve()):
                    raise ValueError('ZIP path escapes working directory')
                if (info.external_attr >> 16) & 0o170000 == 0o120000:
                    raise ValueError('ZIP symlinks are unsupported')
            z.extractall(extracted)
        candidates = [f.parent for f in extracted.rglob('fpm.toml')]
    valid = [p for p in candidates if (p / 'src/fitpack_core.F90').is_file()]
    valid = sorted(set(valid))
    if len(valid) != 1:
        raise ValueError(f'Expected one FITPACK source tree; found {valid}')
    return valid[0]


def capi_order(source: Path) -> list[Path]:
    files = list((source / 'src/capi').glob('*.f90'))
    def modules(file: Path) -> set[str]:
        return set(re.findall(r'^\s*module\s+(?!procedure\b|function\b|subroutine\b)(\w+)',
                              file.read_text(), re.I | re.M))
    def uses(file: Path) -> set[str]:
        return set(re.findall(r'^\s*use\s*(?:,\s*(?:non_)?intrinsic\s*::\s*)?(\w+)',
                              file.read_text(), re.I | re.M))
    provider = {m.lower(): file for file in files for m in modules(file)}
    ordered: list[Path] = []
    while len(ordered) < len(files):
        ready = [f for f in files if f not in ordered and
                 all(provider.get(m.lower()) in (None, *ordered) for m in uses(f))]
        if not ready:
            raise RuntimeError('C API module dependency cycle or missing parser support')
        ordered.extend(ready)
    return ordered


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fitpack', type=Path, help='Original FITPACK ZIP or extracted source directory')
    parser.add_argument('--work', type=Path, help='Directory to retain build products and logs')
    parser.add_argument('--compiler', default='gfortran')
    parser.add_argument('--cxx', default='g++')
    parser.add_argument('--cases', type=int, default=140)
    parser.add_argument('--optimized', action='store_true', help='Build with -O3 instead of checked -O0')
    parser.add_argument('--no-cpp', action='store_true', help='Skip upstream C++ integration tests')
    parser.add_argument('--no-numerical', action='store_true', help='Skip independent SciPy comparison')
    args = parser.parse_args()
    source_input = args.fitpack.resolve()
    if not source_input.exists():
        parser.error(f'Not found: {source_input}')
    if args.work:
        work = args.work.resolve()
    else:
        work = ROOT / 'build' / 'fitpack-offline'
    work.mkdir(parents=True, exist_ok=True)
    fitpack_root = unpack_or_resolve(source_input, work)
    manifest = tomllib.loads((fitpack_root / 'fpm.toml').read_text())
    print(f'FITPACK {manifest.get("version", "unknown")} source: {fitpack_root}', flush=True)
    suffix = '.exe' if os.name == 'nt' else ''
    libdir = work / 'upstream'
    appdir = work / 'adapter'
    libdir.mkdir(exist_ok=True)
    appdir.mkdir(exist_ok=True)
    flags = ['-O3' if args.optimized else '-O0', '-g', '-cpp', '-ffree-line-length-none']
    if not args.optimized:
        flags += ['-fcheck=all', '-fbacktrace']
    # Force a clean rebuild to avoid stale .mod files from other compilers.
    for folder in (libdir, appdir):
        for item in folder.iterdir():
            if item.is_file() and item.suffix.lower() in ('.mod', '.o', '.obj'):
                item.unlink()
    objs: list[Path] = []
    def build_source(file: Path, out: Path, extra_include: Path | None = None) -> Path:
        obj = out / f'{file.stem}.o'
        command = [args.compiler, *flags, '-c', '-J', out, '-I', out]
        if extra_include is not None:
            command += ['-I', extra_include]
        command += [file, '-o', obj]
        invoke(command, out, 'compile-' + file.name)
        return obj
    for name in CORE_SOURCES:
        objs.append(build_source(fitpack_root / 'src' / name, libdir))
    capi = capi_order(fitpack_root)
    for file in capi:
        objs.append(build_source(file, libdir))
    invoke(['ar', 'rcs', libdir / 'libfitpack.a', *objs], libdir, 'archive-fitpack')
    print(f'Compiled FITPACK: {len(CORE_SOURCES)} core + {len(capi)} C-API modules', flush=True)

    if not args.no_cpp:
        test_objs = [build_source(fitpack_root / 'test' / name, libdir) for name in TEST_MODULES]
        cpp_obj = libdir / 'test_curve_cpp.o'
        invoke([args.cxx, '-std=c++17', '-O0', '-c', fitpack_root / 'test/test_curve.cpp',
                '-o', cpp_obj], libdir, 'compile-upstream-cpp-test')
        exe = libdir / ('test_fitpack_all' + suffix)
        invoke([args.compiler, *flags, '-I', libdir, fitpack_root / 'test/test.f90',
                *test_objs, cpp_obj, libdir / 'libfitpack.a', '-lstdc++', '-o', exe],
               libdir, 'link-upstream-tests')
        invoke([exe], libdir, 'run-upstream-tests')
        runlog = (libdir / 'run-upstream-tests.log').read_text()
        assert '[fitpack] SUCCESS!' in runlog, 'Upstream test log does not report success'
        print('Upstream FITPACK complete test executable: PASS', flush=True)

    adapt_objs = [build_source(file, appdir, libdir) for file in SCIFORT_SOURCES]
    test_exe = appdir / ('test_fitpack_smoothing' + suffix)
    invoke([args.compiler, *flags, '-I', appdir, '-I', libdir,
            ROOT / 'optional/fitpack_backend/test/test_fitpack_smoothing.f90',
            *adapt_objs, libdir / 'libfitpack.a', '-o', test_exe],
           appdir, 'link-scifort-adapter-test')
    invoke([test_exe], appdir, 'run-scifort-adapter-test')
    print('SciFort optional adaptive-smoothing adapter: PASS', flush=True)
    periodic_test = appdir / ('test_fitpack_periodic_parametric' + suffix)
    invoke([args.compiler, *flags, '-I', appdir, '-I', libdir,
            ROOT / 'optional/fitpack_backend/test/test_fitpack_periodic_parametric.f90',
            *adapt_objs, libdir / 'libfitpack.a', '-o', periodic_test],
           appdir, 'link-fitpack-periodic-parametric-test')
    invoke([periodic_test], appdir, 'run-fitpack-periodic-parametric-test')
    print('SciFort periodic/parametric adapter native tests: PASS', flush=True)

    bivariate_test = appdir / ('test_fitpack_bivariate' + suffix)
    invoke([args.compiler, *flags, '-I', appdir, '-I', libdir,
            ROOT / 'optional/fitpack_backend/test/test_fitpack_bivariate.f90',
            *adapt_objs, libdir / 'libfitpack.a', '-o', bivariate_test],
           appdir, 'link-fitpack-bivariate-test')
    invoke([bivariate_test], appdir, 'run-fitpack-bivariate-test')
    print('SciFort bivariate FITPACK adapter: PASS', flush=True)

    bivariate_driver = appdir / ('fitpack_bivariate_driver' + suffix)
    invoke([args.compiler, *flags, '-I', appdir, '-I', libdir,
            ROOT / 'optional/fitpack_backend/app/fitpack_bivariate_driver.f90',
            *adapt_objs, libdir / 'libfitpack.a', '-o', bivariate_driver],
           appdir, 'link-fitpack-bivariate-driver')
    if not args.no_numerical:
        bivariate_compare = ROOT / 'tools/interpolate_reference/compare_optional_fitpack_bivariate_scipy.py'
        invoke([sys.executable, bivariate_compare, '--executable', bivariate_driver,
                '--cases', str(min(args.cases,80))], appdir, 'run-fitpack-bivariate-scipy')
        print((appdir / 'run-fitpack-bivariate-scipy.log').read_text()
              .split('\n\n', 1)[-1].strip(), flush=True)

    periodic_driver = appdir / ('fitpack_periodic_parametric_driver' + suffix)
    invoke([args.compiler, *flags, '-I', appdir, '-I', libdir,
            ROOT / 'optional/fitpack_backend/app/fitpack_periodic_parametric_driver.f90',
            *adapt_objs, libdir / 'libfitpack.a', '-o', periodic_driver],
           appdir, 'link-fitpack-periodic-parametric-driver')
    if not args.no_numerical:
        periodic_compare = ROOT / 'tools/interpolate_reference/compare_optional_fitpack_periodic_parametric_scipy.py'
        invoke([sys.executable, periodic_compare, '--executable', periodic_driver,
                '--cases', str(args.cases)], appdir, 'run-fitpack-periodic-parametric-scipy')
        report = (appdir / 'run-fitpack-periodic-parametric-scipy.log').read_text()
        print(report.split('\n\n', 1)[-1].strip(), flush=True)

    driver = appdir / ('fitpack_reference_driver' + suffix)
    invoke([args.compiler, *flags, '-I', appdir, '-I', libdir,
            ROOT / 'optional/fitpack_backend/app/fitpack_reference_driver.f90',
            *adapt_objs, libdir / 'libfitpack.a', '-o', driver],
           appdir, 'link-fitpack-reference-driver')
    if not args.no_numerical:
        compare = ROOT / 'tools/interpolate_reference/compare_optional_fitpack_scipy.py'
        invoke([sys.executable, compare, '--executable', driver,
                '--cases', str(args.cases)], appdir, 'run-fitpack-scipy-comparison')
        print((appdir / 'run-fitpack-scipy-comparison.log').read_text().split('COMMAND:', 1)[-1]
              .split('\n\n', 1)[-1].strip(), flush=True)
    print(f'SUCCESS. Logs, objects, binaries: {work}', flush=True)


if __name__ == '__main__':
    main()
