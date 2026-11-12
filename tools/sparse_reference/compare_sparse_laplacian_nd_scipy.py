#!/usr/bin/env python3
"""Randomized SciPy comparisons for LaplacianNd and sparse structural queries.

Compiles the source tree using gfortran and checks matrix-free products,
canonical CSR matrix data, grid eigenvalues, orthonormal eigenvectors, and
structural bandwidth/triangularity against scipy.sparse.linalg. The known
SciPy 1.17 ``LaplacianNd.tosparse`` bug for singleton axes is isolated:
those cases are compared to the matrix-free ``matvec`` instead.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import sparse
from scipy.sparse import linalg as spla

ROOT = Path(__file__).resolve().parents[2]


def build_drivers(build: Path, mode: str) -> tuple[Path, Path]:
    flags = (['-O0', '-g', '-fcheck=all', '-fbacktrace',
              '-ffpe-trap=invalid,zero,overflow'] if mode == 'checked' else ['-O3'])
    compiler = ['gfortran', '-std=f2018', *flags, '-J', str(build), '-I', str(build)]
    sources = [('kinds', 'src/kinds.f90'), ('sparse', 'src/sparse.f90'),
               ('laplacian', 'src/sparse_laplacian_nd.f90'),
               ('lsmr', 'src/optimize/lsmr.f90'),
               ('sparse_linalg', 'src/sparse_linalg.f90')]
    objects = []
    for name, source in sources:
        obj = build / (name + '.o')
        subprocess.run([*compiler, '-c', str(ROOT/source), '-o', str(obj)],
                       cwd=build,check=True,capture_output=True)
        objects.append(obj)
    drivers = []
    for name in ('sparse_laplacian_nd_reference_driver', 'sparse_structure_reference_driver'):
        exe = build / name
        subprocess.run([*compiler, str(ROOT/'tools/sparse_reference'/f'{name}.f90'),
                        *(str(obj) for obj in objects), '-o', str(exe)],
                       cwd=build,check=True,capture_output=True)
        drivers.append(exe)
    return tuple(drivers)


def numbers(line: str, dtype=float) -> np.ndarray:
    return np.fromstring(line, dtype=dtype, sep=' ')


def assert_close(lhs: np.ndarray, rhs: np.ndarray, context: object, atol=2e-11) -> None:
    assert lhs.shape == rhs.shape, (context, lhs.shape, rhs.shape)
    assert np.allclose(lhs, rhs, rtol=2e-11, atol=atol), (context, np.max(abs(lhs-rhs)))


def check_laplacians(exe: Path, count: int, seed: int) -> None:
    rng = np.random.default_rng(seed)
    cases = []
    lines = [str(count)]
    for k in range(count):
        d = int(rng.integers(1,5))
        shape = tuple(int(x) for x in rng.integers(1 if k % 9 == 0 else 2,6,size=d))
        while np.prod(shape) > 125:
            shape = tuple(int(x) for x in rng.integers(2,5,size=d))
        n = int(np.prod(shape))
        bc_code = k % 3
        bc = ('dirichlet','neumann','periodic')[bc_code]
        m = int(rng.integers(1,min(n,11)+1)) if k % 3 else n
        x = rng.normal(size=n)
        lines.extend((f'{d} {bc_code+1} {m}', ' '.join(map(str,shape)),
                      ' '.join(f'{v:.17e}' for v in x)))
        cases.append((shape, bc, m, x))
    result = subprocess.run([str(exe)], input='\n'.join(lines)+'\n',
                            text=True,capture_output=True,check=True)
    output = iter(result.stdout.splitlines())
    for k,(shape,bc,m,x) in enumerate(cases):
        label = (k,shape,bc,m)
        n = int(np.prod(shape))
        h = next(output).split()
        assert h[0]=='CASE' and int(h[1])==n and int(h[3])==m, (label,h)
        nnz = int(h[2])
        y = numbers(next(output))
        eig = numbers(next(output))
        vec = numbers(next(output)).reshape((n,m),order='F')
        ptr = numbers(next(output),int)
        cols = numbers(next(output),int)
        data = numbers(next(output))
        stats = numbers(next(output),int)
        assert len(ptr)==n+1 and len(cols)==nnz and len(data)==nnz
        assert ptr[0]==1 and ptr[-1]==nnz+1
        assert np.all(np.diff(ptr)>=0)
        actual = sparse.csr_array((data,cols-1,ptr-1),shape=(n,n))
        ref = spla.LaplacianNd(shape,boundary_conditions=bc,dtype=np.float64)
        assert_close(y, np.asarray(ref.matvec(x)).ravel(), (label,'matvec'))
        assert_close(actual @ x,y,(label,'CSR matvec'))
        if min(shape)>1:
            assert_close(actual.toarray(),ref.tosparse().toarray(),(label,'CSR dense'))
        # SciPy 1.17's m-optimization uses lexicographic tuple min,
        # occasionally introducing nonexistent modes (e.g. (5,2,2), m=4).
        # Compare against the full SciPy eigenspectrum's exact top-m slice.
        assert_close(eig,ref.eigenvalues()[-m:],(label,'eigenvalues'))
        assert_close(actual @ vec, vec*eig[np.newaxis,:],(label,'eigenpair'))
        assert_close(vec.T @ vec,np.eye(m),(label,'orthonormal'))
        expected_stats = (*spla.spbandwidth(actual),
                          *(int(z) for z in spla.is_sptriangular(actual)))
        assert tuple(stats[:4]) == expected_stats, (label,'CSR stats',stats,expected_stats)
        assert tuple(stats[4:]) == expected_stats, (label,'CSC stats',stats,expected_stats)
    print(f'LaplacianNd: {count}/{count} SciPy comparison cases pass')


def check_structure(exe: Path, count: int, seed: int) -> None:
    rng = np.random.default_rng(seed)
    lines = [str(count)]
    cases = []
    for k in range(count):
        n,m = (int(q) for q in rng.integers(1,46,size=2))
        nnz = int(rng.integers(0,160))
        if k % 13 == 0:
            nnz = 0
        rr = rng.integers(0,n,size=nnz)
        cc = rng.integers(0,m,size=nnz)
        if k % 7 == 0:
            cc = np.minimum(cc,rr)
        elif k % 7 == 1:
            rr = np.minimum(rr,cc)
        vv = rng.normal(size=nnz)
        vv[rng.random(nnz)<0.3] = 0.0
        lines.append(f'{n} {m} {nnz}')
        lines.extend(f'{r+1} {c+1} {v:.17e}' for r,c,v in zip(rr,cc,vv))
        cases.append(sparse.coo_array((vv,(rr,cc)),shape=(n,m)).tocsr())
    result = subprocess.run([str(exe)],input='\n'.join(lines)+'\n',
                            text=True,capture_output=True,check=True)
    output = result.stdout.splitlines()
    assert len(output) == count, (len(output), count)
    for k,(a,line) in enumerate(zip(cases,output)):
        got = tuple(int(q) for q in line.split())
        assert len(got)==8
        if a.nnz:
            bandwidth = spla.spbandwidth(a)
        else:
            # SciPy 1.17 raises on empty CSR; the Fortran API returns (0,0).
            bandwidth = (0,0)
        expect = (*bandwidth,*(int(t) for t in spla.is_sptriangular(a)))
        assert got[:4] == expect and got[4:] == expect, (k,got,expect)
    print(f'Structural queries: {count}/{count} SciPy comparison cases pass')


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode',choices=('checked','optimized','both'),default='both')
    parser.add_argument('--laplacian-cases',type=int,default=160)
    parser.add_argument('--structure-cases',type=int,default=240)
    parser.add_argument('--seed',type=int,default=20261105)
    args = parser.parse_args()
    for mode in (('checked','optimized') if args.mode=='both' else (args.mode,)):
        with tempfile.TemporaryDirectory(prefix='scifort-laplacian-') as folder:
            laplace_exe,structure_exe = build_drivers(Path(folder),mode)
            print(f'Build: {mode}',flush=True)
            check_laplacians(laplace_exe,args.laplacian_cases,args.seed)
            check_structure(structure_exe,args.structure_cases,args.seed+1)


if __name__=='__main__':
    main()
