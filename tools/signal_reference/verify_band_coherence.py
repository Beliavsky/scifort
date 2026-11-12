#!/usr/bin/env python3
"""Independent NumPy oracle for piecewise-linear integrated CSD/coherence.

Generate Hermitian positive-definite spectral matrices over nonuniform and
wrapped two-sided FFT frequency grids. Independently integrate interpolated
matrices with numpy.trapezoid and calculate Schur/precision quantities with
numpy.linalg.inv. Compare static and time-frequency Fortran outputs.

Usage: python tools/signal_reference/verify_band_coherence.py path/to/driver
"""
from __future__ import annotations

import argparse
import subprocess

import numpy as np


def oracle_integral(freq, matrices, lo, hi, average):
    """Integrate each matrix entry using a separate interpolation/trapezoid path."""
    perm = np.argsort(freq)
    f = freq[perm]
    s = matrices[perm]
    knots = np.concatenate(([lo], f[(f > lo) & (f < hi)], [hi]))
    interpolated = np.empty((knots.size,) + s.shape[1:], dtype=complex)
    for i in range(s.shape[1]):
        for j in range(s.shape[2]):
            real = np.interp(knots, f, s[:, i, j].real)
            imag = np.interp(knots, f, s[:, i, j].imag)
            interpolated[:, i, j] = real + 1j * imag
    result = np.trapezoid(interpolated, knots, axis=0)
    return result / (hi - lo) if average else result


def oracle_coherence(spectrum, ridge, shrinkage):
    """Separate NumPy precision and conditional Schur-complement reference."""
    n = spectrum.shape[0]
    conditioned = spectrum.copy()
    conditioned *= 1.0 - shrinkage
    np.fill_diagonal(conditioned, np.diag(spectrum))
    conditioned += np.eye(n) * (ridge * np.trace(spectrum).real / n)
    precision = np.linalg.inv(conditioned)
    diag = np.diag(conditioned).real
    residual = 1.0 / np.diag(precision).real
    r2 = 1.0 - residual / diag
    partial_complex = -precision / np.sqrt(np.outer(np.diag(precision).real, np.diag(precision).real))
    np.fill_diagonal(partial_complex, 1.0 + 0.0j)
    partial_sq = np.abs(partial_complex) ** 2
    return r2, residual, diag - residual, partial_sq, partial_complex


def cases(seed=3826, count=180):
    rng = np.random.default_rng(seed)
    for case_id in range(count):
        nf = int(rng.integers(3, 15))
        nc = int(rng.integers(2, 6))
        nb = int(rng.integers(1, 5))
        nt = int(rng.integers(1, 4))
        f = np.cumsum(rng.uniform(0.035, 0.75, nf))
        f -= rng.uniform(0.3, 0.7) * f[-1]
        if case_id % 3 == 0:
            # Mimic complex FFT wrapping: nonnegative followed by negative f.
            p = int(np.searchsorted(f, 0.0))
            order = np.r_[np.arange(p, nf), np.arange(p)]
        else:
            order = rng.permutation(nf)
        freq = f[order]
        bands = []
        for ib in range(nb):
            lo = float(rng.uniform(f[0], f[-2]))
            hi = float(rng.uniform(max(lo + 1e-8, f[1]), f[-1]))
            bands.append((lo, hi))
        if case_id % 13 == 0:
            bands[0] = (float(f[0]), float(f[-1]))
        ridge = float(rng.uniform(0.0, 0.12))
        shrinkage = float(rng.uniform(0.0, 0.55))
        average = bool(case_id % 2)

        def matrices():
            raw = rng.standard_normal((nf, nc, nc)) + 1j * rng.standard_normal((nf, nc, nc))
            return (raw @ raw.conj().transpose(0, 2, 1) / nc +
                    np.eye(nc) * 0.25)[order]

        spectra = matrices()
        frames = np.stack([matrices() for _ in range(nt)], axis=0)
        yield (nf, nc, nb, nt, freq, bands, ridge, shrinkage, average, spectra, frames)


def check(actual, reference, name, max_error):
    if not np.allclose(actual, reference, rtol=8e-11, atol=2e-11):
        deviation = np.max(np.abs(actual - reference))
        raise AssertionError(f"{name}: max error {deviation:.3e} vs reference")
    return max(max_error, float(np.max(np.abs(actual - reference))))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("driver", help="compiled native Fortran driver")
    parser.add_argument("--count", type=int, default=180)
    args = parser.parse_args()
    scenarios = list(cases(count=args.count))
    content = [str(len(scenarios))]
    for nf, nc, nb, nt, f, bands, ridge, shrinkage, average, spectra, frames in scenarios:
        content.append(f"{nf} {nc} {nb} {nt} {ridge:.17g} {shrinkage:.17g} {int(average)}")
        content.append(" ".join(f"{v:.17g}" for v in f))
        content.extend(f"{lo:.17g} {hi:.17g}" for lo, hi in bands)
        content.append(" ".join(str(k) for k in range(nt)))
        for matrix in spectra:
            for j in range(nc):
                for i in range(nc):
                    z = matrix[i, j]
                    content.append(f"{z.real:.17g} {z.imag:.17g}")
        for frame in frames:
            for matrix in frame:
                for j in range(nc):
                    for i in range(nc):
                        z = matrix[i, j]
                        content.append(f"{z.real:.17g} {z.imag:.17g}")
    result = subprocess.run([args.driver], input="\n".join(content) + "\n",
                            capture_output=True, text=True, check=True)
    lines = iter(result.stdout.splitlines())

    max_error = 0.0
    comparisons = 0

    def read_row(length):
        v = np.fromstring(next(lines), sep=" ")
        if v.size != length:
            raise ValueError(f"expected {length} floats, got {v}")
        return v

    def read_complex():
        pair = read_row(2)
        return complex(pair[0], pair[1])

    for idx, case in enumerate(scenarios, start=1):
        nf, nc, nb, nt, freq, bands, ridge, shrinkage, average, spectra, frames = case
        header = next(lines)
        if header != f"CASE {idx}":
            raise AssertionError(f"wrong header {header!r}, case {idx}")
        integrated = np.empty((nb, nc, nc), dtype=complex)
        frame_integrated = np.empty((nt, nb, nc, nc), dtype=complex)
        for ib in range(nb):
            for j in range(nc):
                for i in range(nc):
                    integrated[ib, i, j] = read_complex()
        for t in range(nt):
            for ib in range(nb):
                for j in range(nc):
                    for i in range(nc):
                        frame_integrated[t, ib, i, j] = read_complex()
        mc = np.empty((nb, nc, 3))
        pc_sq = np.empty((nb, nc, nc))
        pc_complex = np.empty((nb, nc, nc), dtype=complex)
        for ib in range(nb):
            for i in range(nc):
                mc[ib, i, :] = read_row(3)
            for j in range(nc):
                for i in range(nc):
                    pc_sq[ib, i, j] = read_row(1)[0]
                    pc_complex[ib, i, j] = read_complex()
        tmc = np.empty((nt, nb, nc, 3))
        tpc_sq = np.empty((nt, nb, nc, nc))
        tpc_complex = np.empty((nt, nb, nc, nc), dtype=complex)
        for t in range(nt):
            for ib in range(nb):
                for i in range(nc):
                    tmc[t, ib, i, :] = read_row(3)
                for j in range(nc):
                    for i in range(nc):
                        tpc_sq[t, ib, i, j] = read_row(1)[0]
                        tpc_complex[t, ib, i, j] = read_complex()
        for ib, (lo, hi) in enumerate(bands):
            reference = oracle_integral(freq, spectra, lo, hi, average)
            max_error = check(integrated[ib], reference, f"case {idx} band {ib} integrated", max_error)
            mc_ref = oracle_coherence(reference, ridge, shrinkage)
            for component in range(3):
                max_error = check(mc[ib, :, component], mc_ref[component],
                                  f"case {idx} static multiple component {component}", max_error)
            max_error = check(pc_sq[ib], mc_ref[3], f"case {idx} static partial^2", max_error)
            max_error = check(pc_complex[ib], mc_ref[4], f"case {idx} static coherency", max_error)
            comparisons += 1
            for t in range(nt):
                reference = oracle_integral(freq, frames[t], lo, hi, average)
                max_error = check(frame_integrated[t, ib], reference, f"case {idx} frame {t}", max_error)
                mc_ref = oracle_coherence(reference, ridge, shrinkage)
                for component in range(3):
                    max_error = check(tmc[t, ib, :, component], mc_ref[component],
                                      f"case {idx} TF component {component}", max_error)
                max_error = check(tpc_sq[t, ib], mc_ref[3], f"case {idx} TF partial^2", max_error)
                max_error = check(tpc_complex[t, ib], mc_ref[4], f"case {idx} TF coherency", max_error)
                comparisons += 1
    print(f"PASS independent NumPy band integration + coherence: {args.count} scenarios, "
          f"{comparisons} band/time cells, maximum absolute error {max_error:.3e}")


if __name__ == "__main__":
    main()
