#!/usr/bin/env python3
"""Compare SciFort rank-1..15 spline coefficients and samples with scipy.ndimage.

Run from anywhere with a gfortran executable on PATH. The test builds a
standalone Fortran driver and uses binary streams for precise comparisons.
"""
from pathlib import Path
import shutil
import struct
import subprocess
import tempfile
import numpy as np
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[2]
SHAPES = [(7,), (4, 6), (3, 5, 4), (3, 2, 4, 2),
          (2, 3, 2, 2, 3), (2, 3) + (1,) * 13]


def make_cases():
    rng = np.random.default_rng(107128)
    for shape in SHAPES:
        rank = len(shape)
        for degree in range(6):
            for complex_data in (False, True):
                for constant_mode in (False, True):
                    for trial in range(2):
                        x = rng.standard_normal(shape)
                        if complex_data:
                            x = x + 1j * rng.standard_normal(shape)
                        # Include fractional and out-of-bounds coordinates.
                        coords = np.vstack([
                            rng.uniform(-4.0, n + 4.0, 16) for n in shape
                        ])
                        # Add four integer coordinates: exact lattice reconstruction.
                        knots = np.vstack([rng.integers(0, n, 4) for n in shape])
                        coords = np.hstack([coords, knots.astype(np.float64)])
                        yield rank, degree, complex_data, constant_mode, np.asfortranarray(x), coords


def main():
    cases = list(make_cases())
    compiler = shutil.which('gfortran')
    if compiler is None:
        raise RuntimeError('gfortran required for this comparison')
    max_coefficient_error = 0.0
    max_interpolation_error = 0.0
    n_values = 0
    with tempfile.TemporaryDirectory(prefix='scifort-spline-nd-') as tmp:
        work = Path(tmp)
        cmd = [compiler, '-std=f2018', '-O2', '-fcheck=all',
               '-Werror=implicit-interface', '-ffree-line-length-none',
               str(ROOT/'src/kinds.f90'),
               str(ROOT/'src/signal_spline.f90'),
               str(ROOT/'src/signal_spline_nd.f90'),
               str(ROOT/'tools/signal_reference/spline_nd_reference_driver.f90'),
               '-o', str(work/'spline_nd_driver')]
        subprocess.run(cmd, cwd=work, check=True, capture_output=True, text=True)
        with (work/'spline_nd_reference_input.bin').open('wb') as handle:
            for rank, order, complex_data, constant_mode, x, coords in cases:
                handle.write(struct.pack('=iiiii', rank, order,
                    int(complex_data), int(constant_mode), coords.shape[1]))
                handle.write(np.asarray(x.shape, dtype=np.int32).tobytes())
                handle.write(np.asfortranarray(x).ravel(order='F').tobytes())
                handle.write(np.asfortranarray(coords + 1.0).ravel(order='F').tobytes())
        subprocess.run([str(work/'spline_nd_driver')], cwd=work, check=True,
                       capture_output=True, text=True)
        with (work/'spline_nd_reference_output.bin').open('rb') as handle:
            for rank, order, complex_data, constant_mode, x, coords in cases:
                statuses = struct.unpack('=iii', handle.read(12))
                if statuses != (0, 0, 0):
                    raise AssertionError(('Fortran status',rank,order,statuses))
                dtype = np.complex128 if complex_data else np.float64
                item_count = x.size
                coeff = np.frombuffer(handle.read(item_count * np.dtype(dtype).itemsize),
                                      dtype=dtype).reshape(x.shape, order='F')
                npoints = coords.shape[1]
                values = np.frombuffer(handle.read(npoints * np.dtype(dtype).itemsize), dtype=dtype)
                direct = np.frombuffer(handle.read(npoints * np.dtype(dtype).itemsize), dtype=dtype)
                mode = 'constant' if constant_mode else 'mirror'
                cval = (-2.25 + 0.125j) if complex_data else -2.25
                target_coeff = (ndimage.spline_filter(x, order=order, mode=mode)
                                if order > 1 else x)
                # For rank 15, SciPy's tensor enumerator can require enormous
                # temporary storage, even for singleton axes. Collapse the 13
                # unit-length axes in the independent reference instead.
                if rank == 15:
                    base_shape = x.shape[:2]
                    reduced = x.reshape(base_shape)
                    reduced_coeff = target_coeff.reshape(base_shape)
                    target = ndimage.map_coordinates(reduced_coeff,coords[:2],
                              order=order, mode=mode, cval=cval, prefilter=False)
                    target_direct = ndimage.map_coordinates(reduced,coords[:2],
                              order=order, mode=mode, cval=cval, prefilter=True)
                    if constant_mode:
                        outside = np.any(coords[2:] != 0.0, axis=0)
                        target[outside] = cval
                        target_direct[outside] = cval
                else:
                    target = ndimage.map_coordinates(target_coeff, coords, order=order,
                             mode=mode, cval=cval, prefilter=False)
                    target_direct = ndimage.map_coordinates(x, coords, order=order,
                             mode=mode, cval=cval, prefilter=True)
                cerror = float(np.max(np.abs(coeff-target_coeff)))
                verror = float(max(np.max(np.abs(values-target)),
                                   np.max(np.abs(direct-target_direct))))
                if cerror > 8.e-11 or verror > 8.e-11:
                    raise AssertionError((rank,order,complex_data,constant_mode,
                                          cerror,verror, np.argmax(np.abs(values-target)),
                                          coords[:,np.argmax(np.abs(values-target))],
                                          values[np.argmax(np.abs(values-target))],
                                          target[np.argmax(np.abs(values-target))],
                                          np.argmax(np.abs(direct-target_direct)),
                                          np.max(np.abs(values-target)),
                                          np.max(np.abs(direct-target_direct)),
                                          coeff.ravel(order='F')[:6],
                                          target_coeff.ravel(order='F')[:6],values[:4],target[:4]))
                max_coefficient_error = max(max_coefficient_error,cerror)
                max_interpolation_error = max(max_interpolation_error,verror)
                n_values += npoints
    print(f'PASS {len(cases)} cases, {n_values} coordinate evaluations')
    print(f'Maximum coefficient difference: {max_coefficient_error:.3g}')
    print(f'Maximum interpolation difference: {max_interpolation_error:.3g}')


if __name__ == '__main__':
    main()
