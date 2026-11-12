#!/usr/bin/env python3
"""Compile-check optional rpoly adapter with a signature-equivalent test stub.

This does NOT validate the genuine upstream Jenkins--Traub implementation.
Run `cd optional/polyroots_backend && fpm test` with network access and system
BLAS/LAPACK to validate the actual pinned polyroots-fortran dependency.
"""
from pathlib import Path
import subprocess
import tempfile

ROOT=Path(__file__).resolve().parents[2]
SOURCE=["src/kinds.f90", "src/interpolate_piecewise.f90",
        "src/interpolate_bspline.f90", "src/interpolate_regular_grid.f90",
        "src/interpolate_barycentric.f90", "src/interpolate_global.f90", "src/interpolate_kdtree.f90", "src/interpolate_rbf.f90", "src/interpolate_legacy.f90", "src/interpolate_scattered.f90", "src/interpolate_clough_tocher.f90", "src/interpolate_fitpack_bivariate.f90", "src/interpolate_fitpack_compat.f90", "src/interpolate.f90",
        "tools/interpolate_reference/polyroots_rpoly_test_stub.f90",
        "optional/polyroots_backend/src/scifort_optional_polyroots.f90",
        "optional/polyroots_backend/test/test_optional_polyroots.f90"]
for name,flags in [('checked',['-O0','-g','-fcheck=all']),('optimized',['-O3'])]:
    with tempfile.TemporaryDirectory() as tmp:
        cmd=['gfortran','-std=f2018','-ffree-line-length-none',*flags,
             '-J'+tmp,'-I'+tmp,*[str(ROOT/s) for s in SOURCE],'-o',tmp+'/test']
        p=subprocess.run(cmd,cwd=tmp,capture_output=True,text=True)
        if p.returncode: raise RuntimeError(p.stderr[-5000:])
        p=subprocess.run([tmp+'/test'],cwd=tmp,capture_output=True,text=True)
        if p.returncode: raise RuntimeError(p.stderr[-5000:])
        print(name,': API adapter stub PASS; upstream implementation NOT tested')
