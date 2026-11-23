#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Time Python calling SciFort versus Python calling public scipy.stats methods."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import time
import sys

# compare_scipy sets thread environment before importing NumPy/SciPy.
import compare_scipy as common
import numpy as np
import scipy
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'bindings/python'))
from build_library import build
from scifort import Scifort

CASES = [case for case in common.CASES if case[0].endswith('_array')
         and common.category(case[1]) == 'distribution']
SCIPY_NAMES = dict(normal='norm', exponential='expon', weibull='weibull_min', lognormal='lognorm')


def prepare(api, mode, n):
    _, inputs, _ = common.prepare(mode, n)
    family, method, _ = mode.split('_')
    shape = (1.7,) if family == 'weibull' else ((.8,) if family == 'lognormal' else ())
    fort = getattr(getattr(api, family), method)
    python = getattr(getattr(stats, SCIPY_NAMES.get(family, family)), method)
    # Both public APIs allocate outputs; input construction and library loading are excluded.
    return (lambda k: fort(inputs[:, k], *shape),
            lambda k: python(inputs[:, k], *shape), inputs)


def markdown(report):
    lines = ['# Python SciPy / Python SciFort timing comparison', '',
             f"Started: {report['started_utc']}", '',
             'Both backends are called from Python through their public distribution methods.',
             'Times include Python validation, dispatch, allocation and FFI overhead.',
             'Inputs and library loading/building are outside timings. Both backends allocate output;',
             'SciFort does not use its optional reusable output buffer. No direct scipy.special kernels.',
             'Each batch alternates two preconstructed inputs and consumes endpoint checksums.',
             'Backend order alternates by sample. Ratios above one favor SciFort.', '',
             f"Quick smoke workloads: {report['quick']}. Validation-only: {report['check_only']}.", '']
    if not report['check_only']:
        lines += common.ratio_summary_lines(report['cases'])
        lines += ['', '| Workload | N | Python SciFort ms | Python SciPy ms | SciPy/SciFort |',
                  '| --- | ---: | ---: | ---: | ---: |']
        for row in report['cases']:
            lines.append(f"| {row['case']} | {row['n']} | {row['scifort_median_seconds']*1e3:.6f} | "
                         f"{row['scipy_median_seconds']*1e3:.6f} | {row['scipy_over_scifort']:.3f} |")
    else:
        lines += [f"Validated {len(report['cases'])} workloads; no timings collected."]
    lines += ['', 'Per-call timing samples in milliseconds (sampling order):', '']
    if not report['check_only']:
        lines += ['| Workload | SciFort samples | SciPy samples |', '| --- | --- | --- |']
        for row in report['cases']:
            lines.append('| '+row['case']+' | '+', '.join(f'{s*1e3:.6f}' for s in row['scifort_seconds'])
                         +' | '+', '.join(f'{s*1e3:.6f}' for s in row['scipy_seconds'])+' |')
    lines += ['', 'Compiler, source hashes, versions, validation errors and batch sizes are in the sibling JSON.', '']
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default=os.environ.get('FC', 'gfortran'))
    parser.add_argument('--case', action='append', choices=[c[0] for c in CASES])
    parser.add_argument('--quick', action='store_true')
    parser.add_argument('--check-only', action='store_true')
    parser.add_argument('--list', action='store_true')
    parser.add_argument('--repeats', type=int, default=7)
    parser.add_argument('--target-seconds', type=float, default=.1)
    parser.add_argument('--output', type=Path, default=ROOT / 'build/python-scipy-speed/results.json')
    args = parser.parse_args()
    if args.repeats < 3 or not np.isfinite(args.target_seconds) or args.target_seconds <= 0:
        parser.error('require repeats >= 3 and finite target-seconds > 0')
    cases = [c for c in CASES if not args.case or c[0] in args.case]
    if args.list:
        for name, _, n, quick in cases:
            print(f'{name:32} n={n} quick={quick}')
        return
    path, metadata = build(args.compiler, args.check_only)
    api = Scifort(path, dll_directories=metadata['dll_directories'])
    revision = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, capture_output=True, text=True)
    report = dict(git_revision=revision.stdout.strip(), python_clock=vars(time.get_clock_info('perf_counter')),
                  shared_harness_sha256=hashlib.sha256(Path(common.__file__).read_bytes()).hexdigest(),
                  started_utc=datetime.now(timezone.utc).isoformat(), build=metadata,
                  python=sys.version, numpy=np.__version__, scipy=scipy.__version__,
                  platform=platform.platform(), quick=args.quick, check_only=args.check_only,
                  requested_repeats=args.repeats, target_batch_seconds=args.target_seconds,
                  thread_environment={v: os.environ[v] for v in common.THREAD_ENV},
                  harness_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  binding_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in sorted((ROOT / 'bindings/python/scifort').glob('*.py'))},
                  cases=[])
    for label, mode, normal_size, quick_size in cases:
        n = quick_size if args.quick else normal_size
        fort, python, inputs = prepare(api, mode, n)
        actual = np.column_stack([fort(k) for k in range(2)])
        accuracy, expected = common.validate(actual, python, mode, inputs, None)
        row = dict(case=label, category='distribution', n=n, accuracy=accuracy)
        if not args.check_only:
            fb = lambda count: common.python_batch(fort, count)
            pb = lambda count: common.python_batch(python, count)
            fb(2)
            pb(2)
            nf, npy = common.calibrate(fb, args.target_seconds), common.calibrate(pb, args.target_seconds)
            ftimes, ptimes = [], []
            for sample in range(args.repeats):
                if sample % 2:
                    pt, pc = pb(npy)
                    ft, fc = fb(nf)
                else:
                    ft, fc = fb(nf)
                    pt, pc = pb(npy)
                if ft <= 0 or pt <= 0:
                    raise RuntimeError('clock resolution too coarse')
                common.check_checksum(fc, expected, nf)
                common.check_checksum(pc, expected, npy)
                ftimes.append(ft/nf)
                ptimes.append(pt/npy)
            fm, pm = statistics.median(ftimes), statistics.median(ptimes)
            row.update(scifort_seconds=ftimes, scipy_seconds=ptimes,
                       scifort_batch_size=nf, scipy_batch_size=npy,
                       scifort_median_seconds=fm, scipy_median_seconds=pm,
                       scipy_over_scifort=pm/fm)
            print(f'{label:32} SciFort {fm*1e3:9.5f} ms  SciPy {pm*1e3:9.5f} ms  ratio {pm/fm:.3f}', flush=True)
        else:
            print(f'{label}: agreement passed', flush=True)
        report['cases'].append(row)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    args.output.with_suffix('.md').write_text(markdown(report), encoding='utf-8')
    print(f'Reports: {args.output} and {args.output.with_suffix(".md")}')


if __name__ == '__main__':
    main()
