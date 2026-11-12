#!/usr/bin/env python3
"""Benchmark indexed 3-D tetra location; optionally compare the predecessor ZIP.

Usage:
    python tools/interpolate_reference/benchmark_tetra_bvh.py
    python tools/interpolate_reference/benchmark_tetra_bvh.py --baseline-zip path/to/previous.zip

Compiles from source with GNU Fortran, uses a deterministic benchmark sample,
and verifies baseline/new interpolation checksums when the old ZIP is provided.
No SciPy, BLAS, LAPACK, Qhull, or network access required.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / 'benchmarks' / 'bench_interpolate_tetra_bvh.f90'
SOURCES = ('kinds', 'interpolate_kdtree', 'interpolate_scattered')


def command(args: list[str], cwd: Path) -> str:
    proc = subprocess.run(args, cwd=cwd, text=True, capture_output=True)
    if proc.returncode:
        raise RuntimeError(f"{' '.join(args)} FAILED\n{proc.stdout}\n{proc.stderr[-3000:]}")
    return proc.stdout


def run(source: Path, scratch: Path, n: int, nq: int) -> dict[str, float]:
    scratch.mkdir()
    flags = ['gfortran', '-std=f2018', '-O3', '-ffree-line-length-none']
    for module in SOURCES:
        command([*flags, '-J'+str(scratch), '-I'+str(scratch), '-c',
                 str(source / 'src' / (module + '.f90')),
                 '-o', str(scratch / (module + '.o'))], scratch)
    exe = scratch / 'benchmark'
    command([*flags, '-I'+str(scratch), str(BENCH),
             *map(str, [scratch / (module + '.o') for module in reversed(SOURCES)]),
             '-o', str(exe)], scratch)
    result = command([str(exe), str(n), str(nq)], scratch)
    parsed = dict(line.strip().split('=', 1) for line in result.splitlines())
    print(source.name, '\n' + result.strip(), flush=True)
    return {k: float(v) for k, v in parsed.items()}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--baseline-zip', type=Path)
    parser.add_argument('--sites', type=int, default=160)
    parser.add_argument('--queries', type=int, default=12000)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as temp:
        work = Path(temp)
        new = run(ROOT, work / 'new', args.sites, args.queries)
        if args.baseline_zip is not None:
            with zipfile.ZipFile(args.baseline_zip) as zf:
                zf.extractall(work / 'baseline_source')
            original = work / 'baseline_source' / 'scifort-main'
            if not original.is_dir():
                raise ValueError('Expected scifort-main/ in baseline archive')
            old = run(original, work / 'old', args.sites, args.queries)
            if old['tetrahedra'] != new['tetrahedra'] or old['filled'] != new['filled']:
                raise AssertionError('mesh / out-of-hull result differs')
            if abs(old['checksum'] - new['checksum']) > 1e-9 * max(1, abs(old['checksum'])):
                raise AssertionError('indexed and linear-scan checksums differ')
            print(f"evaluation speedup: {old['evaluate_s']/new['evaluate_s']:.2f}x", flush=True)
            print('baseline/new checksums and out-of-hull counts agree', flush=True)


if __name__ == '__main__':
    main()
