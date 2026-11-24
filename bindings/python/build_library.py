#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Build the distribution-only development shared library, without FPM or parser dependencies."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
FAMILIES = ('normal', 'exponential', 'logistic', 'laplace', 'cauchy', 'weibull', 'lognormal',
            'gamma', 'beta', 't')
FAMILY_SOURCES = {'t': 'student_t'}
SOURCES = ['src/kinds.f90', 'src/constants.f90', 'src/math.f90']
SOURCES += ['src/stats/normal.f90', 'src/special/log_gamma.f90',
            'src/special/incomplete_gamma.f90', 'src/special/incomplete_beta.f90']
SOURCES += [f'src/stats/{FAMILY_SOURCES.get(f, f)}.f90' for f in FAMILIES if f != 'normal']
SOURCES += [f'src/c_api/generated/{FAMILY_SOURCES.get(f, f)}_generated.f90' for f in FAMILIES]
SOURCES += [f'bindings/python/native/{f}_compat.f90' for f in FAMILIES]


def build(compiler='gfortran', checked=False, directory=None):
    compiler = shutil.which(compiler)
    if not compiler:
        raise RuntimeError('Fortran compiler not found')
    version = subprocess.check_output([compiler, '--version'], text=True).strip()
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if checked else ['-O3']
    flags += ['-std=f2018', '-fimplicit-none', '-Werror=implicit-interface', '-shared', '-fPIC']
    digest = hashlib.sha256((compiler + version + repr(flags)).encode())
    for source in SOURCES:
        digest.update((ROOT / source).read_bytes())
    directory = Path(directory) if directory else ROOT / 'build/python-bindings' / digest.hexdigest()[:16]
    directory = directory.resolve()
    directory.mkdir(parents=True, exist_ok=True)
    suffix = '.dll' if os.name == 'nt' else ('.dylib' if sys.platform == 'darwin' else '.so')
    library = directory / ('scifort' + suffix)
    metadata = dict(compiler=version, compiler_path=compiler, flags=flags,
                    source_sha256=digest.hexdigest(), library=str(library),
                    dll_directories=[str(Path(compiler).parent)] if os.name == 'nt' else [])
    cache = directory / 'build.json'
    if not library.exists() or not cache.exists() or json.loads(cache.read_text()) != metadata:
        # Compile/link to a temporary name so failed builds cannot leave a valid-looking cache.
        temporary = library.with_name('pending' + suffix)
        try:
            subprocess.run([compiler, *flags, '-J', str(directory), '-I', str(directory),
                            *[str(ROOT / p) for p in SOURCES], '-o', str(temporary)],
                           cwd=directory, check=True)
            temporary.replace(library)
            cache.write_text(json.dumps(metadata, indent=2)+'\n', encoding='utf-8')
        finally:
            temporary.unlink(missing_ok=True)
    return library, metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default=os.environ.get('FC', 'gfortran'))
    parser.add_argument('--checked', action='store_true')
    parser.add_argument('--directory', type=Path)
    args = parser.parse_args()
    library, _ = build(args.compiler, args.checked, args.directory)
    print(library)


if __name__ == '__main__':
    main()
