#!/usr/bin/env python3
"""Benchmark native 3D Delaunay construction against a predecessor archive.

Compile both variants with GNU Fortran -O3. The benchmark's input coordinates
are deterministic. Reports the median of repeated construction CPU timings;
verifies tetrahedron count, checksum, fill count, and vertex interpolation.
Results are machine- and workload-specific, not guaranteed speedup factors.

Example:
    python tools/interpolate_reference/benchmark_delaunay_3d_cached.py \\
      --baseline-zip ../scifort-interpolate-tetra-bvh-2026-10-31.zip
"""
from __future__ import annotations

import argparse
from pathlib import Path
import statistics
import subprocess
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[2]
SOURCE_FILES = ('kinds', 'interpolate_kdtree', 'interpolate_scattered')
BENCHMARK = ROOT / 'benchmarks' / 'bench_interpolate_tetra_bvh.f90'


def execute(command: list[str], directory: Path) -> str:
    result = subprocess.run(command, cwd=directory, text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError(f"Command failed: {' '.join(command)}\n{result.stderr[-3000:]}")
    return result.stdout


def build(root: Path, scratch: Path) -> Path:
    scratch.mkdir(parents=True)
    compiler = ['gfortran', '-O3', '-std=f2018', '-ffree-line-length-none']
    for filename in SOURCE_FILES:
        execute([*compiler, '-I'+str(scratch), '-J'+str(scratch), '-c',
                 str(root / 'src' / (filename + '.f90')),
                 '-o', str(scratch / (filename + '.o'))], scratch)
    exe = scratch / 'benchmark'
    execute([*compiler, '-I'+str(scratch), str(BENCHMARK),
             *(str(scratch / (name + '.o')) for name in reversed(SOURCE_FILES)),
             '-o', str(exe)], scratch)
    return exe


def measure(exe: Path, scratch: Path, n: int, repetitions: int) -> dict[str, float]:
    observations = []
    for _ in range(repetitions):
        text = execute([str(exe), str(n), '20'], scratch)
        observations.append({k: float(v) for k, v in
                             (line.split('=', 1) for line in text.splitlines() if '=' in line)})
    first = observations[0]
    for other in observations[1:]:
        for key in ('tetrahedra', 'filled', 'checksum'):
            if abs(other[key]-first[key]) > 1e-9*max(1, abs(first[key])):
                raise AssertionError(f'Non-deterministic {key} for {n} sites')
    first['construct_s'] = statistics.median(item['construct_s'] for item in observations)
    return first


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--baseline-zip', type=Path, required=True)
    parser.add_argument('--sizes', type=int, nargs='+', default=[160, 320, 440])
    parser.add_argument('--repetitions', type=int, default=4)
    args = parser.parse_args()
    if args.repetitions < 1 or any(n < 4 for n in args.sizes):
        parser.error('repetitions must be positive and sites must be at least 4')
    with tempfile.TemporaryDirectory() as temporary:
        scratch = Path(temporary)
        with zipfile.ZipFile(args.baseline_zip) as archive:
            archive.extractall(scratch / 'baseline-source')
        baseline = scratch / 'baseline-source' / 'scifort-main'
        if not baseline.is_dir():
            raise ValueError('Expected scifort-main/ at the root of the baseline archive')
        old_exe = build(baseline, scratch / 'baseline-build')
        new_exe = build(ROOT, scratch / 'new-build')
        print('sites tetrahedra baseline_s cached_s construction_speedup', flush=True)
        for sites in args.sizes:
            original = measure(old_exe, old_exe.parent, sites, args.repetitions)
            updated = measure(new_exe, new_exe.parent, sites, args.repetitions)
            for key in ('tetrahedra', 'filled'):
                if original[key] != updated[key]:
                    raise AssertionError(f'{key} mismatch for {sites} sites')
            if abs(original['checksum']-updated['checksum']) > 1e-9*max(1,abs(original['checksum'])):
                raise AssertionError(f'Interpolation checksum mismatch for {sites} sites')
            ratio = original['construct_s']/updated['construct_s']
            print(f"{sites:>5d} {int(updated['tetrahedra']):>11d} "
                  f"{original['construct_s']:.6f} {updated['construct_s']:.6f} {ratio:.2f}x",flush=True)
        print('Tetrahedron counts, checksums, and fill counts agree: PASS', flush=True)


if __name__ == '__main__':
    main()
