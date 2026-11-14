#!/usr/bin/env python3
"""Build and test SciFort's optional adapter against an uploaded polyroots ZIP.

This script does not use FPM or fetch source files from the network. The
upstream source is extracted into a temporary build directory and is never
copied into SciFort. GNU Fortran and (optionally) system LAPACK/BLAS are needed.

Examples
--------
python tools/interpolate_reference/check_optional_polyroots_local.py \
    /path/to/polyroots-fortran-master.zip
python tools/interpolate_reference/check_optional_polyroots_local.py \
    /path/to/polyroots-fortran-master.zip --link-mode lapack

The default 'sections' link mode uses -ffunction-sections and --gc-sections to
link only the Jenkins--Traub path, which does not require LAPACK/BLAS.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import shlex
import subprocess
import tempfile
from zipfile import ZipFile

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SOURCE_FILES = (
    "src/kinds.f90",
    "src/interpolate_piecewise.f90",
    "src/interpolate_bspline.f90",
    "src/interpolate_regular_grid.f90",
    "src/interpolate_barycentric.f90", "src/interpolate_global.f90",
    "src/interpolate_kdtree.f90", "src/interpolate_rbf.f90", "src/interpolate_legacy.f90", "src/interpolate_scattered.f90", "src/interpolate_clough_tocher.f90", "src/interpolate_fitpack_bivariate.f90", "src/interpolate_curfit_core.f90", "src/interpolate_fitpack_compat.f90", "src/interpolate.f90",
    "optional/polyroots_backend/src/scifort_optional_polyroots.f90",
)
TESTS = (
    "optional/polyroots_backend/test/test_optional_polyroots.f90",
    "optional/polyroots_backend/test/test_polyroots_real_stress.f90",
)
BENCHMARK = "optional/polyroots_backend/app/benchmark_polyroots.f90"


def run(command: list[str], cwd: Path) -> str:
    """Run a build or test command and include its output in error reports."""
    completed = subprocess.run(
        command, cwd=cwd, text=True, capture_output=True, check=False
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"Command failed (exit {completed.returncode}):\n"
            + shlex.join(command)
            + "\nstdout:\n"
            + completed.stdout[-7000:]
            + "\nstderr:\n"
            + completed.stderr[-7000:]
        )
    return completed.stdout


def unpack_upstream(archive: Path, destination: Path) -> Path:
    """Extract only the upstream module source from the supplied archive."""
    with ZipFile(archive) as zf:
        members = [
            name
            for name in zf.namelist()
            if name.endswith("/src/polyroots_module.F90")
            or name == "src/polyroots_module.F90"
        ]
        if len(members) != 1:
            raise ValueError("Expected exactly one src/polyroots_module.F90 in ZIP")
        destination.write_bytes(zf.read(members[0]))
    return destination


def build(compiler: str, upstream: Path, build_dir: Path,
          variant: str, link_mode: str) -> None:
    """Compile the actual upstream module and all SciFort adapter sources."""
    build_dir.mkdir(parents=True, exist_ok=True)
    flags = ["-O0", "-g", "-fcheck=all", "-fbacktrace"] if variant == "checked" else ["-O3"]
    sections = ["-ffunction-sections", "-fdata-sections"] if link_mode == "sections" else []
    common = ["-cpp", "-ffree-line-length-none", *flags, *sections]
    objects = []
    upstream_obj = build_dir / "polyroots_module.o"
    # Upstream currently uses a nonstandard COMPLEX*16 declaration. Avoid
    # -std=f2018 for this one source, without modifying third-party code.
    run([compiler, *common, "-J", str(build_dir), "-c", str(upstream),
         "-o", str(upstream_obj)], build_dir)
    objects.append(upstream_obj)
    for relative in SOURCE_FILES:
        source = PROJECT_ROOT / relative
        obj = build_dir / (source.stem + ".o")
        run([compiler, *common, "-std=f2018", "-I", str(build_dir),
             "-J", str(build_dir), "-c", str(source), "-o", str(obj)], build_dir)
        objects.append(obj)
    link_flags = ["-Wl,--gc-sections"] if link_mode == "sections" else ["-llapack", "-lblas"]
    for relative in TESTS + ((BENCHMARK,) if variant == "optimized" else ()):
        source = PROJECT_ROOT / relative
        output = build_dir / source.stem
        run([compiler, *common, "-std=f2018", "-I", str(build_dir),
             "-J", str(build_dir), *[str(obj) for obj in objects],
             str(source), *link_flags, "-o", str(output)], build_dir)
        result = run([str(output)], build_dir)
        print(f"{variant}: {source.name}\n{result.rstrip()}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("upstream_zip", type=Path,
                        help="Uploaded polyroots-fortran source ZIP")
    parser.add_argument("--compiler", default="gfortran")
    parser.add_argument("--link-mode", choices=("sections", "lapack"),
                        default="sections",
                        help="No LAPACK via linker GC (default), or system LAPACK")
    args = parser.parse_args()
    if not args.upstream_zip.is_file():
        parser.error(f"ZIP does not exist: {args.upstream_zip}")
    with tempfile.TemporaryDirectory(prefix="scifort-polyroots-real-") as tmp:
        tmpdir = Path(tmp)
        source = unpack_upstream(args.upstream_zip, tmpdir / "polyroots_module.F90")
        for variant in ("checked", "optimized"):
            build(args.compiler, source, tmpdir / variant, variant, args.link_mode)
    print("All genuine upstream adapter tests PASS.")


if __name__ == "__main__":
    main()
