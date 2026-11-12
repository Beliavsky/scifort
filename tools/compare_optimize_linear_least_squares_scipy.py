"""Compare SciFort NNLS/BVLS with SciPy using random reproducible dense cases.

Build/usage, from the project root (requires GNU Fortran, Python SciPy):
  python tools/compare_optimize_linear_least_squares_scipy.py --compiler gfortran

The driver is built from only the modules needed for this new algorithm;
FPM is not required for this independent reference comparison.
"""
from __future__ import annotations

import argparse
import pathlib
import subprocess
import tempfile

import numpy as np
from scipy.optimize import lsq_linear, nnls

ROOT = pathlib.Path(__file__).resolve().parents[1]
SOURCES = [
    ROOT / "src/kinds.f90",
    ROOT / "src/linalg.f90",
    ROOT / "src/linalg_svd.f90",
    ROOT / "src/optimize/linear_least_squares.f90",
    ROOT / "tools/optimize_linear_least_squares_compare.f90",
]


def make_cases(seed: int = 20261110, count: int = 480):
    rng = np.random.default_rng(seed)
    cases = []
    for case in range(count):
        n = int(rng.integers(1, 9))
        m = int(rng.integers(1, 14))
        a = rng.standard_normal((m, n))
        if case % 11 == 0:
            a[:, -1] = a[:, 0]  # rank deficient
        if case % 17 == 0:
            a[:, :] = 0.0  # entirely zero operator
        x_true = rng.standard_normal(n)
        b = a @ x_true + 0.08 * rng.standard_normal(m)
        if case % 9 == 0:
            b = rng.standard_normal(m) * 2.0
        # Avoid coefficient scaling that causes problems in SciPy itself.
        scale = 10.0 ** int(rng.integers(-3, 4))
        a *= scale
        b *= scale
        mode = 1 if case % 3 == 0 else 0
        if mode == 1:
            lower = np.zeros(n)
            upper = np.full(n, 1e100)
        else:
            lower = rng.uniform(-2.5, 0.2, n)
            upper = lower + rng.uniform(0.2, 3.0, n)
            if case % 5 == 0:
                lower[::2] = -1e100
            if case % 7 == 0:
                upper[1::2] = 1e100
        cases.append((mode, a, b, lower, upper))
    return cases


def input_text(cases):
    lines = [str(len(cases))]
    for mode, a, b, lower, upper in cases:
        m, n = a.shape
        lines.append(f"{mode} {m} {n}")
        lines += [" ".join(f"{v:.17e}" for v in row) for row in a]
        lines.append(" ".join(f"{v:.17e}" for v in b))
        lines.append(" ".join(f"{v:.17e}" for v in lower))
        lines.append(" ".join(f"{v:.17e}" for v in upper))
    return "\n".join(lines) + "\n"


def check_results(cases, output):
    rows = output.strip().splitlines()
    assert len(rows) == 2 * len(cases), (len(rows), len(cases), rows[:10])
    worst_cost_gap = 0.0
    worst_kkt = 0.0
    failures = []
    ref_failures = 0
    for index, (mode, a, b, lower, upper) in enumerate(cases):
        headline = rows[2*index].split()
        xi = np.array([float(v) for v in rows[2*index+1].split()])
        number, status, cost = int(headline[0]), int(headline[1]), float(headline[2])
        assert number == index + 1
        if xi.size != a.shape[1] or not np.all(np.isfinite(xi)):
            failures.append((index, 'nonfinite or wrong size', status))
            continue
        if status < 0:
            failures.append((index, 'status', status))
            continue
        if np.any(xi < lower - 1e-7) or np.any(xi > upper + 1e-7):
            failures.append((index, 'bounds', status))
            continue
        actual_cost = 0.5 * np.linalg.norm(a @ xi - b)**2
        cost_gap = abs(cost-actual_cost) / max(1.0, actual_cost)
        if cost_gap > 1e-7:
            failures.append((index, 'reported cost mismatch', cost_gap))
            continue
        try:
            if mode == 1:
                x_ref, _ = nnls(a, b, maxiter=1000)
            else:
                ref = lsq_linear(a, b, (lower, upper), method='bvls', tol=1e-11, max_iter=300)
                x_ref = ref.x
        except (ValueError, RuntimeError, np.linalg.LinAlgError):
            ref_failures += 1
            continue
        reference_cost = 0.5 * np.linalg.norm(a @ x_ref - b)**2
        gap = max(0.0, actual_cost-reference_cost)/max(1.0,reference_cost)
        worst_cost_gap = max(worst_cost_gap, gap)
        if gap > 2e-7:
            failures.append((index, 'larger cost than SciPy', gap))
            continue
        grad = a.T @ (a @ xi - b)
        violation = 0.0
        for j, gj in enumerate(grad):
            if xi[j] <= lower[j] + 1e-8:
                v = max(0., -gj)
            elif xi[j] >= upper[j] - 1e-8:
                v = max(0., gj)
            else:
                v = abs(gj)
            violation = max(violation, v)
        relative_kkt = violation/max(1.0, np.max(abs(a.T@b)))
        worst_kkt = max(worst_kkt,relative_kkt)
        if relative_kkt>2e-7:
            failures.append((index,'KKT',relative_kkt))
    print(f"Compared {len(cases)} SciFort cases; SciPy reference errors: {ref_failures}; "
          f"worst relative objective gap: {worst_cost_gap:.3e}; "
          f"worst scaled KKT residual: {worst_kkt:.3e}")
    if failures:
        print('FAILURES:', failures[:20])
        raise AssertionError(f'{len(failures)} cases failed')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--compiler', default='gfortran')
    parser.add_argument('--count',type=int,default=480)
    parser.add_argument('--opt',default='-O0')
    args=parser.parse_args()
    cases = make_cases(count=args.count)
    with tempfile.TemporaryDirectory(prefix='scifort_lsq_') as tmp:
        exe = pathlib.Path(tmp) / 'compare'
        cmd = [args.compiler,'-std=f2018','-fcheck=all','-fbacktrace',args.opt,'-o',str(exe)]
        cmd += [str(s) for s in SOURCES]
        subprocess.run(cmd,cwd=tmp,check=True)
        proc=subprocess.run([str(exe)],input=input_text(cases),capture_output=True,text=True,check=True)
        check_results(cases,proc.stdout)

if __name__ == '__main__':
    main()
