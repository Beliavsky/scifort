#!/usr/bin/env python3
"""Build and run the standalone scifort_sparse native Fortran regressions."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=('checked', 'optimized', 'both'), default='both')
    parser.add_argument('--start-at', default=None,
                        help='Resume at a named test (per requested mode)')
    parser.add_argument('--fast-tests', action='store_true',
                        help='Compile library objects with selected mode but use -O0 for test drivers')
    args = parser.parse_args()
    for mode in (('checked', 'optimized') if args.mode == 'both' else (args.mode,)):
        with tempfile.TemporaryDirectory(prefix='scifort-sparse-native-') as temp:
            build = Path(temp)
            flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
            cc = ['gfortran', '-std=f2018', *flags, '-J', str(build), '-I', str(build)]
            test_cc = (['gfortran', '-std=f2018', '-O0', '-J', str(build), '-I', str(build)]
                       if args.fast_tests else cc)
            for name, file in (('kinds', 'src/kinds.f90'),
                               ('random', 'src/random.f90'),
                               ('sparse', 'src/sparse.f90'),
                               ('sparse_coo', 'src/sparse_coo.f90'),
                               ('sparse_coo_nd', 'src/sparse_coo_nd.f90'),
                               ('sparse_coo_nd_ops', 'src/sparse_coo_nd_ops.f90'),
                               ('sparse_coo_nd_complex', 'src/sparse_coo_nd_complex.f90'),
                               ('sparse_complex', 'src/sparse_complex.f90'),
                               ('sparse_coo_nd_complex_ops', 'src/sparse_coo_nd_complex_ops.f90'),
                               ('sparse_coo_nd_index', 'src/sparse_coo_nd_index.f90'),
                               ('sparse_coo_nd_edit', 'src/sparse_coo_nd_edit.f90'),
                               ('sparse_coo_nd_filter', 'src/sparse_coo_nd_filter.f90'),
                               ('sparse_coo_nd_bool', 'src/sparse_coo_nd_bool.f90'),
                               ('sparse_coo_nd_bool_reduce', 'src/sparse_coo_nd_bool_reduce.f90'),
                               ('sparse_dok_lil', 'src/sparse_dok_lil.f90'),
                               ('sparse_laplacian_nd', 'src/sparse_laplacian_nd.f90'),
                               ('sparse_constructors', 'src/sparse_constructors.f90'),
                               ('lsmr', 'src/optimize/lsmr.f90'),
                               ('sparse_linalg', 'src/sparse_linalg.f90'),
                               ('sparse_complex_ilu', 'src/sparse_complex_ilu.f90'),
                               ('sparse_linalg_complex', 'src/sparse_linalg_complex.f90'),
                               ('sparse_linalg_complex_direct', 'src/sparse_linalg_complex_direct.f90'),
                               ('sparse_linalg_complex_triangular', 'src/sparse_linalg_complex_triangular.f90'),
                               ('sparse_linalg_complex_eigen', 'src/sparse_linalg_complex_eigen.f90'),
                               ('sparse_linalg_complex_lsqr', 'src/sparse_linalg_complex_lsqr.f90'),
                               ('sparse_linalg_complex_lsmr', 'src/sparse_linalg_complex_lsmr.f90'),
                               ('sparse_linalg_complex_svds', 'src/sparse_linalg_complex_svds.f90'),
                               ('sparse_linalg_complex_eigs', 'src/sparse_linalg_complex_eigs.f90'),
                               ('sparse_linalg_complex_expm', 'src/sparse_linalg_complex_expm.f90'),
                               ('sparse_linalg_complex_funm', 'src/sparse_linalg_complex_funm.f90')):
                subprocess.run([*cc, '-c', str(ROOT / file),
                                '-o', str(build / f'{name}.o')], cwd=build, check=True)
            resumed = False
            for test in ('test_sparse_complex_funm', 'test_sparse_complex_expm', 'test_sparse_complex_eigs_operator', 'test_sparse_complex_eigs_shift', 'test_sparse_complex_eigs', 'test_sparse_complex_svds', 'test_sparse_complex_lsmr', 'test_sparse_complex_lsqr', 'test_sparse_complex_triangular_eigen', 'test_sparse_complex_direct', 'test_sparse_complex_ilu_bicg_qmr', 'test_sparse_complex', 'test_sparse_coo_nd_complex', 'test_sparse_coo_nd_bool_reduce', 'test_sparse_coo_nd_bool', 'test_sparse_coo_nd_filter', 'test_sparse_coo_nd_edit', 'test_sparse_coo_nd_index', 'test_sparse_coo_nd_ops', 'test_sparse_coo_nd', 'test_sparse_dok_lil', 'test_sparse_coo', 'test_sparse_funm_krylov', 'test_sparse_laplacian_nd', 'test_sparse_recycling', 'test_sparse_qmr_tfqmr', 'test_sparse_bicg_cgs', 'test_sparse_kronsum', 'test_sparse_constructors', 'test_sparse_deletion_reshape', 'test_sparse_block_assignment', 'test_sparse_assembly', 'test_sparse_matrix', 'test_sparse_csc_diagonals',
                         'test_sparse_linalg', 'test_sparse_lsqr_ilu0',
                         'test_sparse_minres', 'test_sparse_eigsh', 'test_sparse_eigs',
                         'test_sparse_eigs_shift', 'test_sparse_arnoldi',
                         'test_sparse_lanczos', 'test_sparse_lanczos_locked',
                         'test_sparse_restarted_shift', 'test_sparse_arnoldi_schur',
                         'test_sparse_arnoldi_block', 'test_sparse_svds',
                         'test_sparse_expm_multiply', 'test_sparse_onenorm_expm', 'test_sparse_direct',
                         'test_sparse_spilu', 'test_sparse_lobpcg',
                         'test_sparse_matrix_functions', 'test_sparse_formats_norm'):
                if args.start_at is not None and test != args.start_at and not resumed:
                    continue
                resumed = True
                exe = build / test
                subprocess.run([*test_cc, str(ROOT / 'test' / (test + '.f90')),
                                str(build / 'kinds.o'), str(build / 'random.o'),
                                str(build / 'sparse.o'), str(build / 'sparse_coo.o'),
                                str(build / 'sparse_coo_nd.o'), str(build / 'sparse_coo_nd_ops.o'),
                                str(build / 'sparse_coo_nd_complex.o'), str(build / 'sparse_complex.o'),
                                str(build / 'sparse_coo_nd_complex_ops.o'),
                                str(build / 'sparse_coo_nd_index.o'), str(build / 'sparse_coo_nd_edit.o'),
                                str(build / 'sparse_coo_nd_filter.o'), str(build / 'sparse_coo_nd_bool.o'),
                                str(build / 'sparse_coo_nd_bool_reduce.o'),
                                str(build / 'sparse_dok_lil.o'), str(build / 'sparse_laplacian_nd.o'),
                                str(build / 'sparse_constructors.o'),
                                str(build / 'lsmr.o'), str(build / 'sparse_linalg.o'),
                                str(build / 'sparse_complex_ilu.o'), str(build / 'sparse_linalg_complex.o'),
                                str(build / 'sparse_linalg_complex_direct.o'),
                                str(build / 'sparse_linalg_complex_triangular.o'),
                                str(build / 'sparse_linalg_complex_eigen.o'),
                                str(build / 'sparse_linalg_complex_lsqr.o'),
                                str(build / 'sparse_linalg_complex_lsmr.o'),
                                str(build / 'sparse_linalg_complex_svds.o'),
                                str(build / 'sparse_linalg_complex_eigs.o'),
                                str(build / 'sparse_linalg_complex_expm.o'),
                                str(build / 'sparse_linalg_complex_funm.o'),
                                '-o', str(exe)], cwd=build, check=True)
                result = subprocess.run([str(exe)], cwd=build, capture_output=True,
                                        text=True, check=True)
                print(f'{mode}: {result.stdout.strip()}', flush=True)


if __name__ == '__main__':
    main()
