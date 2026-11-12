#!/usr/bin/env python3
"""Compare SciFort batched rotations with scipy.spatial.transform.Rotation.

All quaternion comparisons are modulo sign; power, inversion, composition,
matrix and vector application use independent SciPy reference calculations.

Usage: python tools/spatial_reference/compare_rotation_stack_scipy.py --mode both
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parents[2]


def compile_program(tmp: Path, mode: str) -> Path:
    flags = ["-O0", "-g", "-fcheck=all", "-fbacktrace"] if mode == "checked" else ["-O3"]
    fc = ["gfortran", "-std=f2018", "-ffree-line-length-none", *flags]
    objects = []
    for name in ("kinds", "spatial_transform", "spatial_rotation_stack"):
        obj = tmp / (name + ".o")
        subprocess.run([*fc, "-I", str(tmp), "-J", str(tmp), "-c",
                        str(ROOT / "src" / f"{name}.f90"), "-o", str(obj)],
                       cwd=tmp, check=True, capture_output=True)
        objects.append(str(obj))
    exe = tmp / "rotation_stack_driver"
    subprocess.run([*fc, "-I", str(tmp),
                    str(ROOT / "tools" / "spatial_reference" / "rotation_stack_driver.f90"),
                    *objects, "-o", str(exe)], cwd=tmp, check=True, capture_output=True)
    return exe


def make_cases(seed: int, count: int):
    rng = np.random.default_rng(seed)
    cases = []
    for j in range(count):
        n = int(rng.integers(1, 10))
        m = n if j % 3 == 0 else 1
        k = n if j % 4 != 0 else 1
        angles = (rng.uniform(-np.pi, np.pi, n)
                  if j % 5 else np.full(n, np.pi))
        axes = rng.normal(size=(n, 3))
        axes /= np.linalg.norm(axes, axis=1)[:, None]
        a = Rotation.from_rotvec(angles[:, None] * axes)
        b = Rotation.random(k, random_state=rng)
        qa = a.as_quat()
        qb = b.as_quat()
        qa[::2] *= -1
        qb[::2] *= -1
        # Scaling raw input quaternions is legitimate; compare to normalized
        # quaternions because SciPy 1.17 can overflow on raw huge components.
        if j % 8 == 0:
            qa *= 1e150
        elif j % 8 == 1:
            qb *= 1e-150
        elif j % 8 == 2:
            qa *= 1e-150
        elif j % 8 == 3:
            qb *= 1e150
        vectors = rng.normal(size=(m, 3)) * 10 ** rng.uniform(-2, 2)
        weights = rng.uniform(.3, 5, size=n)
        cases.append((qa, qb, vectors, weights))
    return cases


def quat_errors(got: np.ndarray, exp: np.ndarray) -> np.ndarray:
    return np.minimum(np.max(np.abs(got - exp), axis=-1),
                      np.max(np.abs(got + exp), axis=-1))


def compare(exe: Path, seed: int, count: int) -> tuple[float, int]:
    cases = make_cases(seed, count)
    lines = [str(len(cases))]
    for qa, qb, vectors, weights in cases:
        n, m, k = len(qa), len(vectors), len(qb)
        lines.append(f"{n} {m} {k}")
        for row in (*qa, *qb, *vectors):
            lines.append(" ".join(f"{v:.17e}" for v in row))
        lines.append(" ".join(f"{v:.17e}" for v in weights))
    proc = subprocess.run([str(exe)], input="\n".join(lines) + "\n",
                          text=True, capture_output=True, check=True)
    output = iter(proc.stdout.splitlines())
    worst = 0.0
    values_checked = 0
    for i, (qa, qb, vectors, weights) in enumerate(cases):
        n, k, m = len(qa), len(qb), len(vectors)
        status = int(next(output))
        assert status == 0, (i, status)

        def get_rows(rows, columns):
            nonlocal values_checked
            arr = np.array([[float(v) for v in next(output).split()] for _ in range(rows)])
            assert arr.shape == (rows, columns), (i, arr.shape, rows, columns)
            values_checked += arr.size
            return arr

        q = get_rows(n, 4)
        rv = get_rows(n, 3)
        mrp = get_rows(n, 3)
        magnitudes = np.array([float(v) for v in next(output).split()])
        values_checked += n
        mats = get_rows(n, 9)
        got_applied = get_rows(max(n, m), 3)
        got_inverse = get_rows(max(n, m), 3)
        qc = get_rows(max(n, k), 4)
        qp = get_rows(n, 4)
        qm = np.array([float(v) for v in next(output).split()])
        assert qm.shape == (4,)
        values_checked += 4

        a = Rotation.from_quat(qa / np.max(np.abs(qa), axis=1, keepdims=True))
        b = Rotation.from_quat(qb / np.max(np.abs(qb), axis=1, keepdims=True))
        comparisons = (
            float(np.max(np.abs(q - a.as_quat(canonical=True)))),
            float(np.max(np.abs(rv - a.as_rotvec()))),
            float(np.max(np.abs(mrp - a.as_mrp()))),
            float(np.max(np.abs(magnitudes - a.magnitude()))),
            float(np.max(np.abs(np.array([v.reshape((3, 3), order="F") for v in mats])
                                - a.as_matrix()))),
            float(np.max(np.abs(got_applied - a.apply(vectors)))) /
                   max(1.0, float(np.max(np.abs(vectors)))),
            float(np.max(np.abs(got_inverse - a.apply(vectors, inverse=True)))) /
                   max(1.0, float(np.max(np.abs(vectors)))),
            float(np.max(quat_errors(qc, (a * b).as_quat()))),
            float(np.max(quat_errors(qp, (a ** 0.7).as_quat()))),
            float(np.max(quat_errors(qm, a.mean(weights=weights).as_quat()))),
        )
        value = max(comparisons)
        worst = max(value, worst)
        assert value < 1e-8, (i, comparisons, value)
    assert next(output, None) is None, "Unexpected extra Fortran output"
    return worst, values_checked


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("checked", "optimized", "both"), default="both")
    parser.add_argument("--count", type=int, default=260)
    args = parser.parse_args()
    modes = ["checked", "optimized"] if args.mode == "both" else [args.mode]
    for mode in modes:
        with tempfile.TemporaryDirectory(prefix="scifort-rotation-stack-") as td:
            exe = compile_program(Path(td), mode)
            total, nvalues = 0.0, 0
            for seed in (20261102, 20261103):
                error, values = compare(exe, seed, args.count)
                total = max(total, error)
                nvalues += values
                print(f"{mode}: seed={seed}, cases={args.count}, maximum error={error:.3e}", flush=True)
            print(f"{mode}: {2 * args.count} cases, {nvalues} values, max error={total:.3e}: PASS", flush=True)


if __name__ == "__main__":
    main()
