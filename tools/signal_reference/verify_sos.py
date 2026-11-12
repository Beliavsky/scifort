#!/usr/bin/env python3
"""Randomized SciPy regression for 1D SOS and forward/backward filters.

Requires gfortran, NumPy, SciPy. Compiles a standalone driver in a temporary
folder, without adding project dependencies.
"""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy import signal
import scipy

ROOT = Path(__file__).resolve().parents[2]
SOURCES = [
    "kinds.f90", "signal_lfilter.f90", "signal_sos.f90", "signal_filtfilt.f90",
]


def formatted(x: np.ndarray) -> str:
    array = np.asarray(x).ravel()
    if np.iscomplexobj(array):
        return " ".join(f"({v.real:.18e},{v.imag:.18e})" for v in array)
    return " ".join(f"{float(v):.18e}" for v in array)


def main() -> None:
    rng = np.random.default_rng(20261015)
    lines: list[str] = []
    expected: list[tuple[str, np.ndarray]] = []

    def add(tag: str, parameters: str, arrays: list[np.ndarray], reference: np.ndarray) -> None:
        lines.append("\n".join([tag, parameters] + [formatted(a) for a in arrays]))
        expected.append((tag, np.asarray(reference).reshape(-1)))

    for is_complex in (False, True):
        for _ in range(310):
            order = int(rng.integers(1, 9))
            sos = signal.butter(order, float(rng.uniform(0.08, 0.75)), output="sos")
            if is_complex:
                sos = sos.astype(complex)
                sos[:, :3] *= 1 + 1j * rng.uniform(-0.04, 0.04, size=(len(sos), 1))
                sos[:, 4:] += 1j * rng.uniform(-0.025, 0.025, size=(len(sos), 2))
            n = int(rng.integers(1, 150))
            x = rng.normal(size=n)
            zi = rng.normal(size=(len(sos), 2)) * 0.03
            if is_complex:
                x = x + 1j*rng.normal(size=n)
                zi = zi + 1j*rng.normal(size=zi.shape)*0.03
            y, zf = signal.sosfilt(sos,x,zi=zi)
            add("C" if is_complex else "R",f"{len(sos)} {n}",
                [*sos, x, *zi],np.r_[y,zf.ravel(order="F")])

        for _ in range(180):
            order = int(rng.integers(1, 10))
            sos = signal.butter(order, float(rng.uniform(0.1, 0.8)), output="sos")
            if is_complex:
                sos = sos.astype(complex)
                sos[:, :3] *= 1 + 1j * rng.uniform(-0.04, 0.04, size=(len(sos), 1))
                sos[:,4:] += 1j * rng.uniform(-0.02, 0.02, size=(len(sos),2))
            zi = signal.sosfilt_zi(sos)
            add("Y" if is_complex else "Z",f"{len(sos)}",list(sos),zi.ravel(order="F"))

        for _ in range(280):
            order = int(rng.integers(1, 6))
            b, a = signal.butter(order, float(rng.uniform(0.15, 0.8)))
            if is_complex:
                b = b.astype(complex) * (1+0.03j)
                a = a.astype(complex)
                a[1:] += 1j * rng.normal(size=order)*0.004
            n = int(rng.integers(35, 165))
            x = rng.normal(size=n)
            if is_complex:
                x = x + 1j*rng.normal(size=n)
            padtype = str(rng.choice(["odd", "even", "constant", "none"]))
            pad = int(rng.integers(0, 12)) if rng.random()<0.35 else -999
            y = signal.filtfilt(b,a,x,padtype=None if padtype=="none" else padtype,
                **({} if pad==-999 else {"padlen":pad}))
            add("G" if is_complex else "F",f"{len(b)} {len(a)} {n} {padtype} {pad}",
                [b,a,x],y)

        for _ in range(280):
            order = int(rng.integers(1, 9))
            sos = signal.butter(order,float(rng.uniform(0.12,0.78)),output="sos")
            if is_complex:
                sos = sos.astype(complex)
                sos[:, :3] *= 1 + 1j * rng.uniform(-0.03, 0.03, size=(len(sos), 1))
                sos[:,4:] += 1j * rng.uniform(-0.01, 0.01, size=(len(sos),2))
            n = int(rng.integers(50,155))
            x = rng.normal(size=n)
            if is_complex:
                x = x + 1j*rng.normal(size=n)
            padtype = str(rng.choice(["odd", "even", "constant", "none"]))
            pad = int(rng.integers(0, 18)) if rng.random()<0.35 else -999
            y = signal.sosfiltfilt(sos,x,padtype=None if padtype=="none" else padtype,
                **({} if pad==-999 else {"padlen":pad}))
            add("T" if is_complex else "S",f"{len(sos)} {n} {padtype} {pad}",list(sos)+[x],y)

    with tempfile.TemporaryDirectory(prefix="scifort-sos-") as temp:
        tmp = Path(temp)
        cmd = ["gfortran", "-std=f2018", "-O2", "-J",str(tmp),"-I",str(tmp)]
        cmd.extend(str(ROOT/"src"/s) for s in SOURCES)
        cmd.extend([str(ROOT/"tools/signal_reference/sos_driver.f90"),"-o",str(tmp/"driver")])
        subprocess.run(cmd,cwd=tmp,check=True)
        result=subprocess.run([str(tmp/"driver")],input="\n".join(lines)+"\n",
            text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,cwd=tmp,check=True)
    observed=result.stdout.splitlines()
    if len(observed)!=len(expected):
        raise AssertionError(f"expected {len(expected)} outputs; got {len(observed)}; {result.stderr}")
    worst=0.0
    counts:dict[str,int]={}
    for j,(text,(tag,ref)) in enumerate(zip(observed,expected)):
        fields=text.split()
        status, n=map(int,fields[:2])
        nums=np.asarray([float(v) for v in fields[2:]])
        if tag in ("C","Y","G","T"):
            actual=nums[::2]+1j*nums[1::2]
        else:
            actual=nums
        if status!=0 or len(actual)!=n or n!=len(ref):
            raise AssertionError(f"case {j}, {tag}: status={status} size={n}, expected={len(ref)}")
        err=np.max(np.abs(actual-ref)) if n else 0.0
        scale=np.max(np.abs(ref)) if n else 0.0
        tolerance=1e-10+3e-9*scale
        if err>tolerance:
            raise AssertionError(f"case {j}, {tag}: maximum difference={err:g} tolerance={tolerance:g}")
        worst=max(worst,float(err))
        counts[tag]=counts.get(tag,0)+1
    print(f"SciPy {scipy.__version__}: "
          f"{len(expected)} SOS/filtfilt comparisons passed; max abs diff {worst:.3g}")
    print("by case:",counts)


if __name__ == "__main__":
    main()
