#!/usr/bin/env python3
"""Compare SciFort CountSketch to independently generated SciPy sketches.

For each input, create SciPy's sparse CountSketch with a seeded Generator,
extract its signed row mapping, and feed that mapping into the Fortran driver.
Compare the result to the *public* scipy.linalg.clarkson_woodruff_transform,
using a fresh Generator with the same seed. This avoids requiring identical
Fortran and NumPy random number generator algorithms.

Usage: python test/compare_linalg_clarkson_woodruff_scipy.py ./cwt_driver
"""

from __future__ import annotations

import subprocess
import sys

import numpy as np
import scipy
from scipy import sparse
from scipy.linalg import clarkson_woodruff_transform
from scipy.linalg._sketches import cwt_matrix


def create_cases():
    rng = np.random.default_rng(20261109)
    cases = []
    for kind in range(6):  # real, complex, real CSR, complex CSR, real batch, complex batch
        for n in (0, 1, 2, 3, 5, 9, 17, 37):
            for d in (0, 1, 2, 4, 11):
                for k in (0, 1, 2, 4, 13):
                    if k == 0 and n > 0:
                        continue
                    if kind >= 4:
                        batch = 2 if (n + d + k) % 2 else 3
                        a = rng.normal(size=(batch, n, d))
                    else:
                        batch = 1
                        a = rng.normal(size=(n, d))
                    if kind % 2:
                        a = a + 1j * rng.normal(size=a.shape)
                    if kind in (2, 3):
                        mask = rng.random((n, d)) < 0.28
                        a = sparse.csr_array(a * mask)
                    seed = int(rng.integers(0, 2**32))
                    cases.append((kind, n, d, k, batch, a, seed))
    return cases


def make_input(cases):
    data = [str(len(cases))]
    references = []
    max_input_rows = 0
    for kind, n, d, k, batch, a, seed in cases:
        rng = np.random.default_rng(seed)
        s = cwt_matrix(k, n, rng=rng)
        if s.shape != (k, n) or s.nnz != n:
            raise AssertionError(("bad SciPy CountSketch matrix", s.shape, s.nnz))
        buckets = s.indices + 1
        signs = s.data
        data.append(f"{kind} {n} {d} {k} {batch} {a.nnz if sparse.issparse(a) else 0}")
        data.extend(f"{int(row)} {int(sign)}" for row, sign in zip(buckets, signs))
        if sparse.issparse(a):
            data.extend(str(int(v)+1) for v in a.indptr)
            data.extend(str(int(v)+1) for v in a.indices)
            values = a.data
        else:
            values = a.flat
        for value in values:
            data.append(f"{float(np.real(value)):.18e} {float(np.imag(value)):.18e}")
        reference = clarkson_woodruff_transform(a, k, rng=np.random.default_rng(seed))
        if sparse.issparse(reference):
            reference = reference.toarray()
        references.append(np.asarray(reference))
        max_input_rows = max(max_input_rows, n)
    return "\n".join(data) + "\n", references


def run(driver):
    cases = create_cases()
    stdin, expected = make_input(cases)
    completed = subprocess.run(
        [driver], input=stdin, text=True, capture_output=True,
        timeout=180, check=True,
    )
    lines = iter(completed.stdout.splitlines())
    worst_error = 0.0
    kinds = [0] * 6
    count_outputs = 0
    for index, (case, reference) in enumerate(zip(cases, expected), 1):
        kind, n, d, k, batch, a, seed = case
        header = next(lines).split()
        if header != ["CASE", str(index), "0"]:
            raise AssertionError((index, a.shape, seed, "bad Fortran status", header))
        output = np.empty((batch, k, d), dtype=np.complex128)
        for b in range(batch):
            for i in range(k):
                for j in range(d):
                    re, im = map(float, next(lines).split())
                    output[b, i, j] = re + 1j*im
        if kind < 4:
            output = output[0]
        if output.shape != reference.shape:
            raise AssertionError((index, "wrong output shape", output.shape, reference.shape))
        difference = abs(output-reference)
        error = float(np.max(difference)) if difference.size else 0.0
        scale = max(1.0, float(np.max(abs(reference))) if reference.size else 0.0)
        if error / scale > 2e-13:
            raise AssertionError((index, kind, a.shape, seed, error/scale))
        worst_error = max(worst_error, error/scale)
        kinds[kind] += 1
        count_outputs += reference.size
    if next(lines, None) is not None:
        raise AssertionError("Unexpected trailing Fortran output")
    print(f"SciPy {scipy.__version__}: {len(cases)} independent CountSketch cases passed")
    print(f"Cases by type (real, complex, CSR real, CSR complex, batches real, batches complex): {kinds}")
    print(f"Output matrix entries compared: {count_outputs}")
    print(f"Maximum normalized difference from public SciPy CountSketch: {worst_error:.3e}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: compare_linalg_clarkson_woodruff_scipy.py PATH_TO_DRIVER")
    run(sys.argv[1])
