"""Compare new SciFort signal utilities against SciPy and independent WH solves.

Usage:
  python test/compare_signal_missing_extra_scipy.py path/to/test_signal_missing_extra_reference
"""
import sys
import warnings
from math import comb
from pathlib import Path
from subprocess import run

import numpy as np
from scipy import signal
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve


def row(v):
    return " ".join(f"{float(a):.17g}" for a in np.asarray(v).ravel()) + "\n"


def pair_row(v):
    return " ".join(f"({z.real:.17g},{z.imag:.17g})" for z in np.asarray(v).ravel()) + "\n"


def wh_reference(y, w, order, lam):
    n = len(y)
    coeff = [(-1) ** (order - j) * comb(order, j) for j in range(order + 1)]
    diff = diags([np.full(n - order, c, dtype=np.float64) for c in coeff],
                 list(range(order + 1)), shape=(n - order, n), format="csr")
    matrix = diags(w, format="csc") + lam * (diff.T @ diff)
    return spsolve(matrix, w * np.where(w > 0, y, 0))


def main(binary):
    rng = np.random.default_rng(20261024)
    cases = []
    records = []
    # Whittaker-Henderson: smoothness, weighted gaps and several orders.
    for n in [7, 15, 38, 72, 161]:
        for order in [1, 2, 3, 4]:
            if order >= n:
                continue
            for lam in [0., .04, 1., 40., 1000.]:
                y = rng.standard_normal(n)
                w = rng.uniform(.1, 2., size=n)
                if lam > 0:
                    w[rng.choice(n, size=n // 7, replace=False)] = 0
                if lam == 0:
                    w += .1
                exp = wh_reference(y, w, order, lam)
                records.extend(["WH\n", f"{n} {order} {lam}\n", row(y), row(w)])
                cases.append(("WH", exp))
    # 1-D and multi-output transfer-function normalization.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", signal.BadCoefficients)
        for _ in range(95):
            n, m = int(rng.integers(1, 14)), int(rng.integers(1, 12))
            b, a = rng.normal(size=n), rng.normal(size=m)
            b[:int(rng.integers(0, n))] = 0
            a[:int(rng.integers(0, m))] = 0
            if not np.any(a):
                a[-1] = 1
            exp_b, exp_a = signal.normalize(b, a)
            records.extend(["NOR\n", f"{n} {m}\n", row(b), row(a)])
            cases.append(("NOR", (np.asarray(exp_b).ravel(), exp_a)))
        for _ in range(60):
            nr, nc, na = int(rng.integers(2, 6)), int(rng.integers(2, 12)), int(rng.integers(1, 10))
            b, a = rng.normal(size=(nr, nc)), rng.normal(size=na)
            b[:, :int(rng.integers(0, nc))] = 0
            a[:int(rng.integers(0, na))] = 0
            if not np.any(a):
                a[-1] = 1
            exp_b, exp_a = signal.normalize(b, a)
            records.extend(["NOM\n", f"{nr} {nc} {na}\n", *[row(v) for v in b], row(a)])
            cases.append(("NOM", (np.atleast_2d(exp_b), exp_a)))
    for complex_flag in [False, True]:
        for correlate_flag in [False, True]:
            for j in range(90):
                n = int(rng.choice([1, 2, 4, 15, 39, 130, 201, 512]))
                m = int(rng.choice([1, 2, 7, 43, 94, 203]))
                x, y = rng.normal(size=n), rng.normal(size=m)
                if complex_flag:
                    x = x + 1j * rng.normal(size=n)
                    y = y + 1j * rng.normal(size=m)
                mode = ["full", "same", "valid"][j % 3]
                method = ["auto", "fft", "direct"][j % 3]
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    expect = (signal.correlate if correlate_flag else signal.convolve)(x, y, mode=mode,
                                                                                     method=method)
                tag = ("COR" if correlate_flag else "CON") + ("C" if complex_flag else "R")
                records.extend([tag + "\n", f"{n} {m} {mode} {method}\n",
                                pair_row(x) if complex_flag else row(x),
                                pair_row(y) if complex_flag else row(y)])
                cases.append((tag, expect))

    cp = run([str(binary)], input="".join(records), text=True, capture_output=True, check=True)
    lines = iter(cp.stdout.splitlines())
    maximum_absolute = 0.
    maximum_scaled = 0.
    values = 0
    for i, (tag, expectation) in enumerate(cases):
        h = [int(t) for t in next(lines).split()]
        if h[0] != 0:
            raise AssertionError(f"{i}: {tag} returned status {h[0]}")
        if tag == "WH":
            size = h[1]
            actual = np.fromstring(next(lines), sep=" ")
            assert size == len(expectation) == len(actual), (i, tag, size)
            targets = [expectation]
            results = [actual]
        elif tag == "NOR":
            actual_b = np.fromstring(next(lines), sep=" ")
            actual_a = np.fromstring(next(lines), sep=" ")
            assert h[2] == actual_b.size and h[3] == actual_a.size
            targets = list(expectation)
            results = [actual_b, actual_a]
        elif tag == "NOM":
            actual_b = np.array([np.fromstring(next(lines), sep=" ") for _ in range(h[2])])
            actual_a = np.fromstring(next(lines), sep=" ")
            assert actual_b.shape == expectation[0].shape, (i, actual_b.shape, expectation[0].shape)
            targets = list(expectation)
            results = [actual_b, actual_a]
        else:
            actual = np.fromstring(next(lines), sep=" ")
            if tag.endswith("C"):
                actual = actual + 1j * np.fromstring(next(lines), sep=" ")
            assert h[1] == np.asarray(expectation).size == actual.size, (i, tag)
            targets = [expectation]
            results = [actual]
        for ref, got in zip(targets, results):
            diff = np.abs(np.asarray(ref) - got)
            if diff.size:
                err = float(np.max(diff))
                scale = max(1., float(np.max(np.abs(ref))))
                maximum_absolute = max(maximum_absolute, err)
                maximum_scaled = max(maximum_scaled, err / scale)
                if err > 8e-9 * scale:
                    raise AssertionError(f"{i}: {tag} max error {err:.3g}, scale {scale:.3g}")
                values += diff.size
    print(f"Passed {len(cases)} cases / {values} numerical values; max abs error "
          f"{maximum_absolute:.3g}, max scaled error {maximum_scaled:.3g}")


if __name__ == "__main__":
    path = Path(sys.argv[1] if len(sys.argv) > 1 else "build_signal_extra/test_signal_missing_extra_reference")
    main(path)
