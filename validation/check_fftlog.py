# SPDX-License-Identifier: MIT AND BSD-3-Clause
# Coefficient-magnitude calculation adapts SciPy v1.15.3 _fftlog_backend.py.
# Copyright (c) 2001-2002 Enthought, Inc. 2003-2024, SciPy Developers.
# Full terms: LICENSES/scipy-fftlog-BSD-3-Clause.txt.
"""Compare FFTLog with SciPy, upstream fixtures, and independent Mellin sums."""
from pathlib import Path
import argparse
import ast
import subprocess
import warnings

import mpmath as mp
import numpy as np
import scipy
from scipy import fft, special

ROOT = Path(__file__).resolve().parents[1]


def mellin_sum(x, dln, mu, shift, bias, inverse):
    """Direct O(n**2) Fourier sums at 90 digits, without an FFT backend."""
    with mp.workdps(90):
        n = len(x)
        h, order, offset, q = map(mp.mpf, (dln, mu, shift, bias))
        xp, xm = (order+1+q)/2, (order+1-q)/2
        center = mp.mpf(n-1)/2
        pole = lambda z: z <= 0 and z == mp.floor(z)
        coefficients = [mp.mpc(0) for _ in range(n)]
        for m in range(n//2+1):
            omega = 2*mp.pi*m/(n*h)
            if m == 0 and (pole(xp) or pole(xm)):
                if pole(xp) and pole(xm):
                    u = 2**q*(-1)**int(xp-xm)*mp.factorial(-xm)/mp.factorial(-xp)
                    if inverse:
                        u = 1/u
                else:
                    u = 0
            else:
                u = mp.exp(q*mp.log(2)+1j*omega*(mp.log(2)-offset))
                u *= mp.gamma(xp+1j*omega/2)/mp.gamma(xm-1j*omega/2)
                if n % 2 == 0 and m == n//2:
                    u = mp.re(u)
                if inverse:
                    u = 1/mp.conj(u)
            coefficients[m] = u
            if m > 0 and m != n-m:
                coefficients[n-m] = mp.conj(u)
        a = [mp.mpf(float(z))*mp.exp((1 if inverse else -1)*q*(j-center)*h
                                   + (q*offset if inverse else 0)) for j, z in enumerate(x)]
        spectrum = [sum(a[j]*mp.exp(-2j*mp.pi*m*j/n) for j in range(n))/n
                    for m in range(n)]
        result = []
        for j in range(n):
            value = sum(spectrum[m]*coefficients[m]*mp.exp(2j*mp.pi*m*(n-1-j)/n)
                        for m in range(n))
            value *= mp.exp((1 if inverse else -1)*q*(j-center)*h
                            - (q*offset if not inverse else 0))
            result.append(float(mp.re(value)))
        return np.array(result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", default="gfortran")
    parser.add_argument("--module-dir", type=Path, required=True)
    parser.add_argument("--library", type=Path, required=True)
    args = parser.parse_args()
    rng = np.random.default_rng(20261112)
    inputs, records = [], []
    counts = {"scipy": 0, "mpmath": 0, "upstream": 0}

    def add(op, x, dln, mu, shift=0., bias=0., axis=-1, reference=None, category="scipy"):
        x = np.asarray(x)
        if op in (4, 5):
            dims = x.shape
            selected = axis % 3
            samples = np.moveaxis(x, selected, -1)
        else:
            dims, selected, samples = (len(x), 1, 1), 0, x
        if reference is None:
            if op == 3:
                reference = np.array([fft.fhtoffset(dln, mu, shift, bias)])
            else:
                function = fft.ifht if op in (2, 5) else fft.fht
                with warnings.catch_warnings():
                    warnings.simplefilter("error")
                    reference = function(samples, dln, mu, shift, bias)
                if op in (4, 5):
                    reference = np.moveaxis(reference, -1, selected)
        reference = np.asarray(reference)
        scale = max(np.max(np.abs(x), initial=0), np.max(np.abs(reference), initial=0),
                    np.finfo(float).tiny)
        label = f"op={op},shape={x.shape},dln={dln},mu={mu},offset={shift},bias={bias}"
        comparison_scale = scale
        if op == 3:
            comparison_scale = max(scale, 1., abs(dln))
        else:
            length = samples.shape[-1]
            if abs(bias)*(length-1)*abs(dln) > 8 and x.size:
                function = fft.ifht if op in (2, 5) else fft.fht
                # Bound rounding amplification through the individual bias
                # diagonals and FFTs, not the already-cancelled dense matrix.
                omega = 2*np.pi*np.arange(length//2+1)/(length*dln)
                xp, xm = (mu+1+bias)/2, (mu+1-bias)/2
                with np.errstate(invalid="ignore"):
                    plus = special.loggamma(xp+0.5j*omega)
                    minus = special.loggamma(xm+0.5j*omega)
                    u = np.exp(bias*np.log(2)+plus.real-minus.real
                               + 1j*(plus.imag+minus.imag+omega*(np.log(2)-shift)))
                if not np.isfinite(u[0]):
                    u[0] = 2**bias*special.poch(xm, xp-xm)
                if length % 2 == 0:
                    u[-1] = u[-1].real
                amplitudes = np.abs(u)
                if op in (2, 5):
                    amplitudes = np.divide(1., amplitudes, out=np.zeros_like(amplitudes),
                                           where=amplitudes != 0)
                coordinate = (np.arange(length)-(length-1)/2)*dln
                pre = -bias*coordinate
                post = -bias*(coordinate+shift)
                if op in (2, 5):
                    pre, post = -post, -pre
                weighted_input = np.max(np.abs(samples)*np.exp(pre), initial=0)
                bound = np.sqrt(length)*np.max(amplitudes)*weighted_input*np.max(np.exp(post))
                comparison_scale = max(scale, bound)
        records.append((label, reference.ravel(order="F"), scale, comparison_scale, category))
        counts[category] += 1
        inputs.append(f"{op} {x.size} {dln:.17g} {mu:.17g} {shift:.17g} {bias:.17g} "
                      + " ".join(map(str, dims))+f" {selected+1}")
        if x.size:
            inputs.append(" ".join(f"{z:.17g}" for z in x.ravel(order="F")))

    for n in [1, 2, 3, 7, 8, 16, 31, 64, 97, 128, 257]:
        x = rng.normal(size=n)
        for mu in [-4., -2.5, -1.2, -.5, 0., .3, 2., 10., 100.]:
            for bias in [-.8, 0., .8]:
                for dln in [-.2, .05, .7]:
                    for optimal in [False, True]:
                        shift = fft.fhtoffset(dln, mu, .37, bias) if optimal else .37
                        for op in [1, 2]:
                            xp, xm = (mu+1+bias)/2, (mu+1-bias)/2
                            pp = xp <= 0 and xp == round(xp)
                            pm = xm <= 0 and xm == round(xm)
                            if pp != pm and ((op == 1 and pp) or (op == 2 and pm)):
                                continue
                            add(op, x, dln, mu, shift, bias)
    for dln in [-10., -.01, .001, .1, 1., 10., 1000.]:
        for mu in [-100., -4., -.5, 0., 2., 100., 1000.]:
            for bias in [-2., 0., .5, 2.]:
                for initial in [-100., -.5, 0., .5, 100.]:
                    add(3, [], dln, mu, initial, bias)
    for shape in [(2, 3, 7), (3, 8, 2), (1, 1, 1), (0, 2, 3)]:
        x = rng.normal(size=shape)
        for axis in range(3):
            if shape[axis] == 0:
                continue
            for op in [4, 5]:
                for bias in [-.3, 0., .3]:
                    add(op, x, .17, -.5, .1, bias, axis)
    for amplitude in [1e-200, 1e200]:
        x = amplitude*rng.normal(size=33)
        for op in [1, 2]:
            for bias in [-.3, 0., .3]:
                add(op, x, .2, .5, -.3, bias)
    for mu, bias in [(-4., 1.), (-2.5, .5), (-3.5, .5), (-3., 0.)]:
        x = rng.normal(size=17)
        for op in [1, 2]:
            xp, xm = (mu+1+bias)/2, (mu+1-bias)/2
            pp, pm = xp <= 0 and xp == round(xp), xm <= 0 and xm == round(xm)
            if pp != pm and ((op == 1 and pp) or (op == 2 and pm)):
                continue
            add(op, x, .2, mu, .1, bias)
    for n in [1, 5, 8, 9]:
        for mu, bias in [(-4., 1.), (-.5, -.3), (0., 0.), (2., .3), (100., 0.)]:
            for amplitude in [1e-200, 1., 1e200]:
                x = amplitude*rng.normal(size=n)
                for op in [1, 2]:
                    reference = mellin_sum(x, .17, mu, .13, bias, op == 2)
                    add(op, x, .17, mu, .13, bias, reference=reference, category="mpmath")
    rng_sensitive = np.random.default_rng(20261112)
    for n in [1, 2, 3, 7, 8, 16, 31]:
        x = rng_sensitive.normal(size=n)
    for op, bias in [(1, -.8), (2, .8)]:
        shift = fft.fhtoffset(.7, 100., .37, bias)
        reference = mellin_sum(x, .7, 100., shift, bias, op == 2)
        add(op, x, .7, 100., shift, bias, reference=reference, category="mpmath")
    # Inspect installed upstream fixtures rather than vendor their source/data.
    path = Path(scipy.__file__).parent/"fft/tests/test_fftlog.py"
    if path.exists():
        tree = ast.parse(path.read_text())
        test = next(n for n in tree.body if isinstance(n, ast.FunctionDef)
                    and n.name == "test_fht_agrees_with_fftlog")
        fixtures = [ast.literal_eval(n.value) for n in test.body if isinstance(n, ast.Assign)
                    and isinstance(n.value, ast.List)]
        r = np.logspace(-4, 4, 16)
        x = r**1.3*np.exp(-r*r/2)
        dln = np.log(r[1]/r[0])
        for (bias, optimal), fixture in zip([(0., False), (0., True), (.8, True), (-.8, True)],
                                          fixtures):
            shift = fft.fhtoffset(dln, .3, bias=bias) if optimal else 0.
            add(1, x, dln, .3, shift, bias, reference=fixture, category="upstream")
    print(f"SciPy {scipy.__version__}; {len(records)} FFTLog cases; {counts}", flush=True)
    stdin = "\n".join(inputs)+"\n"
    for name, flags in [("checked", ["-O0", "-g", "-fcheck=all",
                                   "-ffpe-trap=invalid,zero,overflow", "-fbacktrace"]),
                        ("optimized", ["-O3", "-fcheck=bounds"])]:
        directory = ROOT/"build/fftlog_reference"/name
        directory.mkdir(parents=True, exist_ok=True)
        exe = directory/"driver.exe"
        subprocess.run([args.compiler, *flags, "-Werror=implicit-interface",
                        "-J"+str(directory), "-I"+str(directory),
                        "-I"+str(args.module_dir.resolve()), str(ROOT/"src/fftlog.f90"),
                        str(ROOT/"src/fft.f90"), str(ROOT/"validation/fftlog_reference_driver.f90"),
                        str(args.library.resolve()), "-o", str(exe)], check=True)
        process = subprocess.run([str(exe)], input=stdin, text=True,
                                 capture_output=True, check=True)
        lines = iter(process.stdout.splitlines())
        worst, worst_label, worst_bound, worst_direct = 0., "", 0., 0.
        for label, expected, scale, comparison_scale, category in records:
            code, count = map(int, next(lines).split())
            if code != 0 or count != expected.size:
                raise AssertionError((label, code, count))
            actual = np.array([float(next(lines)) for _ in range(count)])
            difference = np.max(np.abs(actual-expected), initial=0)
            error = difference/scale
            bounded_error = difference/comparison_scale
            if not np.isfinite(bounded_error) or bounded_error > 3e-12:
                raise AssertionError((category, label, error, bounded_error))
            worst_bound = max(worst_bound, bounded_error)
            if category == "mpmath":
                worst_direct = max(worst_direct, error)
            if error > worst:
                worst, worst_label = error, label
        if next(lines, None) is not None:
            raise AssertionError("unexpected output")
        print(f"{name}: {len(records)} passed; max raw scaled error {worst:.3e}; "
              f"max conditioning-aware error {worst_bound:.3e}; "
              f"max mpmath raw error {worst_direct:.3e}; {worst_label}",
              flush=True)


if __name__ == "__main__":
    main()
