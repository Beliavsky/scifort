# SPDX-License-Identifier: MIT
"""Check shifts/2D FFTs against SciPy and fast lengths against a heap oracle."""
from pathlib import Path
import argparse
import bisect
import heapq
import subprocess

import numpy as np
import scipy
from scipy import fft

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--module-dir", type=Path, required=True)
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--compiler", default="gfortran")
    args = parser.parse_args()
    rng = np.random.default_rng(20261112)
    inputs, records = [], []
    names = ["fft2", "ifft2", "rfft2", "irfft2", "hfft2", "ihfft2", "fftshift", "ifftshift"]
    norms = ["backward", "forward", "ortho"]

    def add(op, x, axes=None, s=None, mode=0, real_input=False):
        if real_input:
            x = x.real
        function = getattr(fft, names[op-11])
        kwargs = dict(axes=axes)
        if op < 17:
            if axes is None:
                kwargs.pop("axes")
            kwargs.update(s=s, norm=norms[mode])
        expected = function(x, **kwargs)
        label = f"{names[op-11]} shape={x.shape}, axes={axes}, s={s}, norm={norms[mode]}"
        scale = max(np.max(np.abs(x), initial=0), np.max(np.abs(expected), initial=0),
                    np.finfo(float).tiny)
        records.append((label, expected.ravel(order="F"), scale, op >= 17))
        inputs.append(" ".join(map(str, [op, int(real_input), 2, mode, -1,
                                       -1 if s is None else len(s),
                                       -1 if axes is None else len(axes),
                                       *x.shape, *expected.shape])))
        if s:
            inputs.append(" ".join(map(str, s)))
        if axes:
            inputs.append(" ".join(str(a+1) for a in axes))
        inputs.extend(f"{complex(z).real:.17g} {complex(z).imag:.17g}"
                      for z in x.ravel(order="F"))

    for shape in [(2, 3, 5, 4), (1, 2, 7, 3), (2, 1, 1, 1)]:
        x = rng.normal(size=shape)+1j*rng.normal(size=shape)
        for op in range(11, 17):
            for axes in [None, (0, 2), (3, 1)]:
                selected = (2, 3) if axes is None else axes
                for s in [None, (3, 8), (1, 1)]:
                    if op in (14, 15) and s is None and shape[selected[-1]] == 1:
                        continue
                    for mode in range(3):
                        for real_input in ([True] if op in (13, 16) else [False, True]):
                            add(op, x, axes, s, mode, real_input)
        for op in [17, 18]:
            for axes in [None, (), (0,), (3, 1), (3, 2, 1, 0)]:
                for real_input in [False, True]:
                    add(op, x, axes, real_input=real_input)
    for shape in [(0, 2, 3, 4), (2, 0, 3, 4), (1, 1, 1, 1)]:
        x = np.zeros(shape, dtype=complex)
        for op in [17, 18]:
            add(op, x)
    for amplitude in [1e-100, 1e100]:
        x = amplitude*(rng.normal(size=(2, 3, 7, 5))+1j*rng.normal(size=(2, 3, 7, 5)))
        for op in range(11, 17):
            for mode in range(3):
                add(op, x, axes=(2, 0), mode=mode, real_input=op in (13, 16))
    # Generate the complete representable 7-smooth set using a priority queue,
    # independently of the Fortran nested prime-power search.
    maximum = np.iinfo(np.int32).max
    pending, seen, smooth = [1], {1}, []
    while pending:
        value = heapq.heappop(pending)
        smooth.append(value)
        for prime in [2, 3, 5, 7]:
            candidate = value*prime
            if candidate <= maximum and candidate not in seen:
                seen.add(candidate)
                heapq.heappush(pending, candidate)
    targets = {-maximum, -1, 0, maximum, maximum-1, smooth[-1], smooth[-1]+1}
    targets.update(range(1, 2049))
    targets.update(map(int, rng.integers(1, maximum, size=500)))
    for prime in [2, 3, 5, 7]:
        value = 1
        while value <= maximum:
            targets.update(t for t in [value-1, value, value+1] if t <= maximum)
            value *= prime
    for target in sorted(targets):
        if target < 1:
            above = below = 0
        else:
            index = bisect.bisect_left(smooth, target)
            above = smooth[index] if index < len(smooth) else 0
            below = smooth[bisect.bisect_right(smooth, target)-1]
        expected = np.array([complex(above, int(above == 0)), complex(below, int(below == 0))]*2)
        records.append((f"fast lengths target={target}", expected, 1., True))
        inputs.append(f"19 0 {target} 0 -1 -1 -1 0 1 1 1 4 1 1 1")
    stdin = "\n".join(inputs)+"\n"
    print(f"SciPy {scipy.__version__}; {len(records)} helper/2D cases; "
          f"{len(smooth)} representable 7-smooth lengths", flush=True)
    for name, flags in [("checked", ["-O0", "-g", "-fcheck=all",
                                   "-ffpe-trap=invalid,zero,overflow", "-fbacktrace"]),
                        ("optimized", ["-O3", "-fcheck=bounds"])]:
        directory = ROOT/"build/fft_helpers_reference"/name
        directory.mkdir(parents=True, exist_ok=True)
        exe = directory/"driver.exe"
        subprocess.run([args.compiler, *flags, "-Werror=implicit-interface",
                        "-J"+str(directory), "-I"+str(directory),
                        "-I"+str(args.module_dir.resolve()), str(ROOT/"src/fft.f90"),
                        str(ROOT/"validation/fft_nd_reference_driver.f90"),
                        str(args.library.resolve()), "-o", str(exe)], check=True)
        completed = subprocess.run([str(exe)], input=stdin, text=True,
                                   capture_output=True, check=True)
        lines = iter(completed.stdout.splitlines())
        worst = 0.
        for label, expected, scale, exact in records:
            code, count = map(int, next(lines).split())
            if code != 0 or count != expected.size:
                raise AssertionError((label, code, count))
            actual = np.array([complex(*map(float, next(lines).split())) for _ in range(count)])
            if exact:
                if not np.array_equal(actual, expected):
                    raise AssertionError((label, actual, expected))
            else:
                error = np.max(np.abs(actual-expected), initial=0)/scale
                if not np.isfinite(error) or error > 3e-12:
                    raise AssertionError((label, error))
                worst = max(worst, error)
        if next(lines, None) is not None:
            raise AssertionError("unexpected output")
        print(f"{name}: {len(records)} passed; max transform scaled error {worst:.3e}",
              flush=True)


if __name__ == "__main__":
    main()
