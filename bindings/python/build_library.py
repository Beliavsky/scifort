#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Build the reviewed development shared library, without FPM or parser dependencies."""
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
PLAN_PATH = Path(__file__).with_name('build_sources.json')
try:
    PLAN = json.loads(PLAN_PATH.read_text(encoding='utf-8'))
except (OSError, ValueError) as error:
    raise RuntimeError('Missing/invalid generated build plan; regenerate the reviewed interface pipeline') from error
if (not isinstance(PLAN, dict) or PLAN.get('schema_version') != 1 or not isinstance(PLAN.get('sources'), list)
        or not isinstance(PLAN.get('families'), list)
        or set(PLAN.get('source_sha256', {})) != set(PLAN['sources'])):
    raise RuntimeError('Invalid generated build plan; regenerate the reviewed interface pipeline')
FAMILIES = tuple(PLAN['families'])
SOURCES = tuple(PLAN['sources'])


def build(compiler='gfortran', checked=False, directory=None):
    compiler = shutil.which(compiler)
    if not compiler:
        raise RuntimeError('Fortran compiler not found')
    version = subprocess.check_output([compiler, '--version'], text=True).strip()
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if checked else ['-O3']
    flags += ['-std=f2018', '-fimplicit-none', '-Werror=implicit-interface', '-shared', '-fPIC']
    digest = hashlib.sha256((compiler + version + repr(flags)).encode())
    digest.update(PLAN_PATH.read_bytes())
    for source in SOURCES:
        path = ROOT / source
        normalized = path.read_text(encoding='utf-8').encode('utf-8')
        if hashlib.sha256(normalized).hexdigest() != PLAN['source_sha256'][source]:
            raise RuntimeError(f'{source} changed; regenerate the reviewed interface pipeline before building')
        digest.update(path.read_bytes())
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
