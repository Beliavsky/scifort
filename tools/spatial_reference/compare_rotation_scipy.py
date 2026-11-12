#!/usr/bin/env python3
"""Independent SciPy reference for SciFort 3D rotations and piecewise Slerp.

Checks scalar-last quaternions, rotation vectors, matrix conversions,
composition order, inverse vector application, magnitudes, and interpolated
orientations. Quaternion comparisons are performed on SO(3), not on a
particular choice of q or -q.

Usage: python tools/spatial_reference/compare_rotation_scipy.py --mode both
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.spatial.transform import Rotation, Slerp

ROOT = Path(__file__).resolve().parents[2]


def compile_program(tmp: Path, mode: str, driver: str) -> Path:
    flags = ["-O0", "-g", "-fcheck=all", "-fbacktrace"] if mode == "checked" else ["-O3"]
    cmd = ["gfortran", "-std=f2018", *flags]
    objects = []
    for name in ("kinds", "spatial_transform"):
        obj = tmp / (name + ".o")
        subprocess.run([*cmd, "-I", str(tmp), "-J", str(tmp), "-c",
                        str(ROOT / "src" / f"{name}.f90"), "-o", str(obj)],
                       cwd=tmp, check=True, capture_output=True)
        objects.append(str(obj))
    exe = tmp / driver
    subprocess.run([*cmd, "-I", str(tmp),
                    str(ROOT / "tools" / "spatial_reference" / (driver + ".f90")),
                    *objects, "-o", str(exe)], cwd=tmp, check=True, capture_output=True)
    return exe


def close_rotation(got: np.ndarray, expected: np.ndarray) -> float:
    got = np.asarray(got, dtype=float)
    expected = np.asarray(expected, dtype=float)
    # q and -q are equivalent. This is also reliable for rotations of pi.
    return float(min(np.max(np.abs(got - expected)),
                     np.max(np.abs(got + expected))))


def make_cases(seed: int, n: int) -> list[tuple[np.ndarray, np.ndarray, np.ndarray, float]]:
    rng = np.random.default_rng(seed)
    cases = []
    for i in range(n):
        a = Rotation.random(random_state=rng)
        b = Rotation.random(random_state=rng)
        mode = i % 17
        if mode == 1:
            b = a
        elif mode == 2:
            b = Rotation.from_quat(-a.as_quat())
        elif mode == 3:
            b = a * Rotation.from_rotvec(1e-9 * rng.normal(size=3))
        elif mode == 4:
            b = a * Rotation.from_rotvec(1e-14 * rng.normal(size=3))
        elif mode == 5:
            axis = rng.normal(size=3)
            axis /= np.linalg.norm(axis)
            a = Rotation.from_rotvec((np.pi - 1e-8) * axis)
        elif mode == 6:
            axis = rng.normal(size=3)
            axis /= np.linalg.norm(axis)
            a = Rotation.from_rotvec((np.pi + 1e-8) * axis)
        elif mode == 7:
            a = Rotation.from_rotvec(1e-12 * rng.normal(size=3))
        elif mode == 8:
            a = Rotation.identity()
        elif mode == 9:
            b = Rotation.identity()
        elif mode == 10:
            a = Rotation.from_rotvec(np.pi * np.array([1.0, 0.0, 0.0]))
        elif mode == 11:
            a = Rotation.from_rotvec(np.pi * np.array([0.0, 1.0, 0.0]))
        elif mode == 12:
            a = Rotation.from_rotvec(np.pi * np.array([0.0, 0.0, 1.0]))
        qa = a.as_quat().copy()
        qb = b.as_quat().copy()
        if i % 4 == 0:
            qa *= -1.0
        if i % 3 == 0:
            qb *= -1.0
        if i % 15 == 0:
            qa *= 1e200
        elif i % 15 == 1:
            qb *= 1e-200
        v = rng.normal(size=3) * float(10.0 ** rng.uniform(-2.0, 2.0))
        t = float((0.0, 1.0, .5)[i % 3]) if i % 9 == 0 else float(rng.uniform(0, 1))
        cases.append((qa, qb, v, t))
    return cases


def check_rotation(exe: Path, seed: int, count: int) -> tuple[float, float]:
    cases = make_cases(seed, count)
    lines = [str(len(cases))]
    for qa, qb, v, fraction in cases:
        lines.append(" ".join(f"{x:.17e}" for x in (*qa, *qb, *v, fraction)))
    proc = subprocess.run([str(exe)], input="\n".join(lines) + "\n",
                          text=True, capture_output=True, check=True)
    output = iter(proc.stdout.splitlines())
    max_error = 0.0
    max_slerp = 0.0
    for i, (qa, qb, v, fraction) in enumerate(cases):
        code = int(next(output))
        assert code == 0, (i, code)
        data = np.array([float(x) for x in next(output).split()])
        assert data.size == 39, (i, data.size)
        q, mat, rv, applied, inverted = data[:4], data[4:13], data[13:16], data[16:19], data[19:22]
        magnitude = data[22]
        qm, qr, qc, qs = data[23:27], data[27:31], data[31:35], data[35:39]
        # SciPy 1.17 normalizes raw 1e200 inputs with an overflowing squared norm.
        # Normalize for the oracle; the Fortran driver receives the raw input.
        a = Rotation.from_quat(qa / np.max(np.abs(qa)))
        b = Rotation.from_quat(qb / np.max(np.abs(qb)))
        expected_q = a.as_quat(canonical=True)
        expected_mat = a.as_matrix()
        got_mat = mat.reshape((3,3),order="F")
        checks = [
            np.max(np.abs(q - expected_q)),
            np.max(np.abs(got_mat - expected_mat)),
            np.max(np.abs(rv - a.as_rotvec())),
            np.max(np.abs(applied - a.apply(v))) / max(1.0,np.max(abs(v))),
            np.max(np.abs(inverted - a.apply(v, inverse=True))) / max(1.0,np.max(abs(v))),
            abs(magnitude - a.magnitude()),
            close_rotation(qm,a.as_quat()),
            close_rotation(qr,a.as_quat()),
            close_rotation(qc,(a*b).as_quat()),
        ]
        if max(checks) >= 3e-12:
            raise AssertionError((i, 'rotation', max(checks), checks))
        expected_slerp = Slerp([0,1],Rotation.from_quat([a.as_quat(),b.as_quat()]))(
            [fraction]).as_quat()[0]
        slerp_err = close_rotation(qs,expected_slerp)
        if slerp_err >= 6e-12:
            raise AssertionError((i,'slerp',slerp_err, fraction))
        max_error = max(max_error, max(checks))
        max_slerp = max(max_slerp,slerp_err)
    assert not list(output), 'Unexpected rotation driver output'
    return max_error, max_slerp


def check_paths(exe: Path, seed: int, count: int) -> float:
    rng = np.random.default_rng(seed)
    maximum = 0.0
    for i in range(count):
        n = 2 + i % 9
        nq = 1 + i % 20
        # Strictly increasing, nonuniform times, including negative times.
        times = np.cumsum(rng.uniform(.01,4.0,n)) - 7.0
        quats = Rotation.random(n,random_state=rng).as_quat()
        if i % 3 == 0:
            quats[1::2] *= -1
        if i % 7 == 0:
            quats[-1] = -quats[-2]
        if i % 11 == 0:
            quats[1:] = quats[0]
        query = rng.uniform(times[0],times[-1],nq)
        if nq >= 2:
            query[0] = times[0]
            query[1] = times[-1]
        lines = [f"{n} {nq}"," ".join(f"{t:.17e}" for t in times)]
        lines.extend(" ".join(f"{x:.17e}" for x in q) for q in quats)
        lines.append(" ".join(f"{t:.17e}" for t in query))
        proc = subprocess.run([str(exe)],input="\n".join(lines)+"\n",
                              text=True,capture_output=True,check=True)
        output = proc.stdout.splitlines()
        assert int(output[0]) == 0, (i,output[:4],proc.stderr)
        got = np.array([[float(z) for z in line.split()] for line in output[1:]])
        ref = Slerp(times,Rotation.from_quat(quats))(query).as_quat()
        assert got.shape == ref.shape
        delta = max(close_rotation(a,b) for a,b in zip(got,ref))
        assert delta < 8e-12, (i,delta)
        maximum = max(maximum,delta)
    return maximum


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=("checked","optimized","both"), default="both")
    ap.add_argument("--seed", type=int, default=17231)
    ap.add_argument("--count", type=int, default=450)
    ap.add_argument("--paths", type=int, default=110)
    args = ap.parse_args()
    modes = ("checked","optimized") if args.mode == "both" else (args.mode,)
    for mode in modes:
        with tempfile.TemporaryDirectory(prefix="scifort-rotation-scipy-") as td:
            tmp = Path(td)
            driver = compile_program(tmp,mode,"rotation_driver")
            other = compile_program(tmp,mode,"rotation_slerp_path_driver")
            rotation_err, slerp_err = check_rotation(driver,args.seed,args.count)
            path_err = check_paths(other,args.seed+100,args.paths)
            print(f"{mode}: {args.count} rotations, {args.paths} Slerp paths PASS; "
                  f"max rotation error={rotation_err:.3g}, "
                  f"max scalar Slerp error={slerp_err:.3g}, "
                  f"max piecewise Slerp error={path_err:.3g}",flush=True)


if __name__ == "__main__":
    main()
