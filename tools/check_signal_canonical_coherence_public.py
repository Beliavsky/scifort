#!/usr/bin/env python3
"""Compile the full scifort_signal closure and test canonical spectral public API."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

from check_signal_mixed_public import closure

ROOT = Path(__file__).resolve().parents[1]
TESTS = (
    'test_signal_canonical_coherence_public.f90',
    'test_signal_group_coherence_public.f90',
    'test_signal_band_coherence_public.f90',
    'test_signal_multiple_coherence_public.f90',
    'test_signal_partial_shrinkage_cv_public.f90',
    'test_signal_partial_shrinkage_auto_public.f90',
    'test_signal_partial_coherency_public.f90',
    'test_signal_partial_coherence_public.f90',
    'test_signal_tf_partial_coherence_public.f90',
)


def check(folder: Path) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    order = closure()
    flags = ['gfortran', '-std=f2018', '-fimplicit-none',
             '-ffree-line-length-none', '-O0', '-J.', '-I.']
    print(f'Public-module dependency closure: {len(order)} sources', flush=True)
    objects = []
    for index, source in enumerate(order, 1):
        obj = folder / (source.stem + '.o')
        subprocess.run([*flags, '-c', str(source), '-o', str(obj)],
                       cwd=folder, check=True)
        objects.append(str(obj))
        if index % 15 == 0 or index == len(order):
            print(f'Compiled {index}/{len(order)}', flush=True)
    for name in TESTS:
        executable = folder / Path(name).stem
        subprocess.run([*flags, str(ROOT / 'test' / name), *objects,
                        '-o', str(executable)], cwd=folder, check=True)
        subprocess.run([str(executable)], cwd=folder, check=True)
    print(f'PASS {len(TESTS)} public canonical coherence regressions', flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir', type=Path)
    args = parser.parse_args()
    if args.build_dir is not None:
        check(args.build_dir.resolve())
    else:
        with tempfile.TemporaryDirectory(prefix='scifort-canonical-public-') as tmp:
            check(Path(tmp))


if __name__ == '__main__':
    main()
