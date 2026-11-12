#!/usr/bin/env python3
"""Compare native SciFort SO(3) cubic splines with scipy RotationSpline.

Each run creates reproducible random 2-160 keyframe configurations. Evaluates
orientations, angular velocities, angular accelerations at keyframes, segment
interiors, and extrapolated points. No SciPy source is used by the Fortran
program; this file supplies reference outputs separately.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.spatial.transform import Rotation, RotationSpline

ROOT = Path(__file__).resolve().parents[2]


def build(tmp: Path, flags: list[str]) -> Path:
    fc = ["gfortran", "-std=f2018", "-ffree-line-length-none", *flags]
    objs = []
    for name in ("kinds", "spatial_transform", "spatial_rotation_spline"):
        obj = tmp / (name + ".o")
        subprocess.run([*fc, "-J", str(tmp), "-I", str(tmp), "-c",
                        str(ROOT / "src" / f"{name}.f90"), "-o", str(obj)],
                       cwd=tmp, check=True, capture_output=True, text=True)
        objs.append(str(obj))
    exe = tmp / "rotation_spline_driver"
    subprocess.run([*fc, "-I", str(tmp),
                    str(ROOT / "tools/spatial_reference/rotation_spline_driver.f90"),
                    *objs, "-o", str(exe)], cwd=tmp, check=True, capture_output=True, text=True)
    return exe


def generate(seed: int, cases: int):
    rng = np.random.default_rng(seed)
    for i in range(cases):
        n = 160 if i % 61 == 60 else int(rng.integers(2, 13))
        scale = float(10.0 ** rng.uniform(-2, 2))
        dt = rng.uniform(0.1, 2.5, n-1)*scale
        times = np.r_[0.0, np.cumsum(dt)]
        # Four regimes: ordinary, tiny increments, near-pi jumps, fixed axis,
        # and identical orientations. Quaternions may change sign.
        regime = i % 5
        if regime == 0:
            rotations = Rotation.random(n, random_state=rng)
        elif regime == 1:
            orientations = np.cumsum(rng.normal(scale=1.0e-5,size=(n,3)),axis=0)
            rotations = Rotation.from_rotvec(orientations)
        elif regime == 2:
            axis = rng.normal(size=(n-1,3))
            axis /= np.linalg.norm(axis,axis=1)[:,None]
            jumps = Rotation.from_rotvec(axis*(np.pi-rng.uniform(1e-3,0.3,(n-1,1))))
            rr=[Rotation.random(random_state=rng)]
            for jump in jumps:
                rr.append(rr[-1]*jump)
            rotations=Rotation.from_quat([r.as_quat() for r in rr])
        elif regime == 3:
            axis = rng.normal(size=3)
            axis /= np.linalg.norm(axis)
            angle = np.cumsum(rng.uniform(-1.2,1.2,n))
            rotations = Rotation.from_rotvec(angle[:,None]*axis[None,:])
        else:
            q = Rotation.random(random_state=rng).as_quat()
            rotations = Rotation.from_quat(np.tile(q,(n,1)))
        sample = np.r_[times, times[:-1] + 0.45*dt,
                       rng.uniform(-dt[0],times[-1]+dt[-1],10)]
        rng.shuffle(sample)
        yield times, rotations, sample


def check(exe: Path, seed: int, cases: int) -> tuple[int,float,float,float,float]:
    data = []
    scenarios = []
    for times,rots,query in generate(seed,cases):
        data.append(f"{len(times)} {len(query)}")
        data.append(" ".join(f"{v:.17e}" for v in times))
        for quat in rots.as_quat():
            data.append(" ".join(f"{v:.17e}" for v in quat))
        data.append(" ".join(f"{v:.17e}" for v in query))
        scenarios.append((times,rots,query))
    proc = subprocess.run([str(exe)], input="\n".join(data)+"\n", text=True,
                          capture_output=True, timeout=150)
    if proc.returncode:
        raise RuntimeError(f"native driver failed: {proc.stderr[-2500:]}")
    rows = iter(proc.stdout.splitlines())
    max_r = max_v = max_a = max_absolute = 0.0
    values = 0
    rejected = 0
    for case,(times,rots,query) in enumerate(scenarios):
        status = int(next(rows))
        ref = RotationSpline(times,rots)
        if status == 2:
            # SciPy's fixed-iteration nonlinear solve also has a divergent
            # regime for some near-180-degree trajectories. When its knot
            # rates are millions of times the angular rate implied by a
            # segment's principal rotation, a stable reference orientation
            # can no longer be inferred from double-precision coefficients.
            scale = np.max(np.abs(ref(times,order=1))) * np.min(np.diff(times)) / np.pi
            if scale < 1.0e7:
                raise AssertionError(f"case {case}: unexpectedly rejected a regular path, scale={scale}")
            rejected += 1
            continue
        if status:
            raise AssertionError(f"case {case}: unexpected Fortran status {status}")
        got = np.array([np.fromstring(next(rows),sep=" ") for _ in query])
        if got.shape != (len(query),10) or not np.all(np.isfinite(got)):
            raise AssertionError(f"case {case}: malformed values")
        error_rotation = (Rotation.from_quat(got[:,:4]).inv()*ref(query)).magnitude()
        rate = ref(query,order=1)
        acc = ref(query,order=2)
        e_v = np.linalg.norm(got[:,4:7]-rate,axis=1)
        e_a = np.linalg.norm(got[:,7:10]-acc,axis=1)
        norm_v = np.linalg.norm(rate,axis=1)
        norm_a = np.linalg.norm(acc,axis=1)
        norm_factor = max(1.,(1./np.min(np.diff(times))))
        err_v = e_v / (norm_factor + norm_v)
        err_a = e_a / (norm_factor**2+norm_a)
        max_r = max(max_r,float(error_rotation.max()))
        max_v = max(max_v,float(err_v.max()))
        max_a = max(max_a,float(err_a.max()))
        max_absolute = max(max_absolute,float(e_a.max()))
        if error_rotation.max() > 1e-8 or err_v.max()>1e-8 or err_a.max()>1e-7:
            worst = int(np.argmax(err_a))
            raise AssertionError(f"case={case}, n={len(times)}, times={times}, "
                                 f"r={error_rotation.max():.2g}, v={err_v.max():.2g}, "
                                 f"a={err_a.max():.2g}, worst t={query[worst]:.4g}")
        values += len(query)
    return values,max_r,max_v,max_a,max_absolute,rejected


def main() -> None:
    p=argparse.ArgumentParser()
    p.add_argument("--mode",choices=("checked","optimized","both"),default="both")
    p.add_argument("--seed",type=int,default=41)
    p.add_argument("--cases",type=int,default=240)
    a=p.parse_args()
    modes=[]
    if a.mode!="optimized":modes.append(("checked",["-O0","-g","-fcheck=all","-fbacktrace"]))
    if a.mode!="checked":modes.append(("optimized",["-O3"]))
    for mode,flags in modes:
        with tempfile.TemporaryDirectory(prefix="scifort-rot-spline-") as td:
            exe=build(Path(td),flags)
            n,r,v,acc,absolute,rejected=check(exe,a.seed,a.cases)
            print(f"{mode}: seed={a.seed} cases={a.cases}, compared={a.cases-rejected}, "
                  f"rejected-divergent={rejected}, samples={n}, "
                  f"max rotation={r:.3e} rad, normalized rate={v:.3e}, "
                  f"normalized acceleration={acc:.3e}, absolute acceleration={absolute:.3e} PASS",flush=True)


if __name__=="__main__":
    main()
