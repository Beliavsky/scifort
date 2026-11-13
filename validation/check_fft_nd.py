# SPDX-License-Identifier: MIT
"""Compare native multidimensional transforms with public SciPy definitions."""
from pathlib import Path
import argparse
import subprocess

import numpy as np
import scipy
from scipy import fft

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", default="gfortran")
    parser.add_argument("--module-dir", type=Path, required=True)
    parser.add_argument("--library", type=Path, required=True)
    args = parser.parse_args()
    rng = np.random.default_rng(20261112)
    names = ["fftn", "ifftn", "rfftn", "irfftn", "hfftn", "ihfftn",
             "dctn", "idctn", "dstn", "idstn"]
    norms = ["backward", "forward", "ortho"]
    records, inputs = [], []

    def add(op, x, s=None, axes=None, mode=0, kind=2, orth=-1, real_input=False):
        if real_input:
            x = x.real
        kwargs = dict(s=s, axes=axes, norm=norms[mode])
        if op >= 7:
            kwargs["type"] = kind
            if orth >= 0:
                kwargs["orthogonalize"] = bool(orth)
        function = getattr(fft, names[op-1])
        # SciPy 1.15.3 omits orthogonalize in its complex DCT/DST path.
        if op >= 7 and not real_input:
            expected = function(x.real, **kwargs)+1j*function(x.imag, **kwargs)
        else:
            expected = function(x, **kwargs)
        label = (f"{names[op-1]} shape={x.shape}, s={s}, axes={axes}, "
                 f"norm={norms[mode]}, type={kind}, orth={orth}, real={real_input}")
        scale = max(np.max(np.abs(x), initial=0),
                    np.max(np.abs(expected), initial=0), np.finfo(float).tiny)
        records.append((label, expected.ravel(order="F"), scale))
        header = [op, int(real_input), kind, mode, orth,
                  -1 if s is None else len(s), -1 if axes is None else len(axes),
                  *x.shape, *expected.shape]
        inputs.append(" ".join(map(str, header)))
        if s:
            inputs.append(" ".join(map(str, s)))
        if axes:
            inputs.append(" ".join(str(a+1) for a in axes))
        inputs.extend(f"{complex(z).real:.17g} {complex(z).imag:.17g}"
                      for z in x.ravel(order="F"))

    for shape in [(3, 4, 5, 2), (2, 3, 8, 1), (1, 1, 17, 7)]:
        x = rng.normal(size=shape)+1j*rng.normal(size=shape)
        for axes in [None, (0,), (3, 1), (2, 0, 1), (3, 2, 1, 0)]:
            selected = tuple(range(4)) if axes is None else axes
            for resize in [False, True]:
                s = tuple(max(2, shape[a]+(-1 if j % 2 else 2))
                          for j, a in enumerate(selected)) if resize else None
                for mode in range(3):
                    for op in range(1, 7):
                        for real_input in ([True] if op in (3, 6) else [False, True]):
                            if op in (4, 5) and s is None and shape[selected[-1]] == 1:
                                continue
                            add(op, x, s, axes, mode, real_input=real_input)
                    for kind in range(1, 5):
                        for orth in [-1, 0, 1]:
                            for op in range(7, 11):
                                if kind == 1 and op in (7, 8):
                                    lengths = s if s is not None else [shape[a] for a in selected]
                                    if min(lengths) < 2:
                                        continue
                                for real_input in [False, True]:
                                    add(op, x, s, axes, mode, kind, orth, real_input)
    # Inferred last axes from s; empty axes identity; empty cropped/padded input;
    # large/small amplitudes; half-spectrum defaults and explicit singleton length.
    for shape, s, axes in [((2, 3, 4, 5), (7, 3), None),
                          ((2, 3, 4, 5), (), ()),
                          ((0, 3, 4, 2), (3,), (0,)),
                          ((2, 3, 4, 0), (3,), (3,)),
                          ((1, 1, 1, 1), (1,), (3,))]:
        x = rng.normal(size=shape)+1j*rng.normal(size=shape)
        for op in range(1, 11):
            if op in (3, 4, 5, 6) and axes == ():
                continue
            add(op, x, s, axes, real_input=op in (3, 6))
    for amplitude in [1e-100, 1e100]:
        x = amplitude*(rng.normal(size=(2, 3, 17, 5))+1j*rng.normal(size=(2, 3, 17, 5)))
        for op in range(1, 11):
            for mode in range(3):
                add(op, x, axes=(2, 0), mode=mode, real_input=op in (3, 6))
    stdin = "\n".join(inputs)+"\n"
    print(f"SciPy {scipy.__version__}; {len(records)} multidimensional cases", flush=True)
    for name, flags in [("checked", ["-O0", "-g", "-fcheck=all",
                                   "-ffpe-trap=invalid,zero,overflow", "-fbacktrace"]),
                        ("optimized", ["-O3", "-fcheck=bounds"])]:
        directory = ROOT/"build/fft_nd_reference"/name
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
        worst, worst_label = 0., ""
        for label, expected, scale in records:
            code, count = map(int, next(lines).split())
            if code != 0 or count != expected.size:
                raise AssertionError(f"{label}: status={code}, size={count}")
            actual = np.array([complex(*map(float, next(lines).split())) for _ in range(count)])
            error = np.max(np.abs(actual-expected), initial=0)/scale
            if not np.isfinite(error) or error > 3e-12:
                raise AssertionError(f"{label}: scaled error {error:.3e}")
            if error > worst:
                worst, worst_label = error, label
        if next(lines, None) is not None:
            raise AssertionError("unexpected trailing output")
        print(f"{name}: {len(records)} passed; max scaled error {worst:.3e}; "
              f"{worst_label}", flush=True)


if __name__ == "__main__":
    main()
