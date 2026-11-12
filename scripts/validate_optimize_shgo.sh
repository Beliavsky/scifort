#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
# Standalone checked/optimized SHGO core tests and independent SciPy comparisons.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEMP="$(mktemp -d)"
trap 'rm -rf "$TEMP"' EXIT
for mode in checked optimized; do
    mkdir -p "$TEMP/$mode"
    cd "$TEMP/$mode"
    if [[ "$mode" == checked ]]; then
        flags=(-std=f2008 -O0 -g -Wall -Wextra -fcheck=all \
            -ffpe-trap=invalid,zero,overflow)
    else
        flags=(-std=f2008 -O3 -Wall -Wextra)
    fi
    gfortran "${flags[@]}" -c \
        "$ROOT/src/kinds.f90" "$ROOT/src/optimize/shgo.f90"
    gfortran "${flags[@]}" "$ROOT/test/optimize_shgo_reference.f90" \
        kinds.o shgo.o -o native_shgo
    ./native_shgo
    gfortran "${flags[@]}" "$ROOT/test/optimize_shgo_adaptive_reference.f90" \
        kinds.o shgo.o -o native_shgo_adaptive
    ./native_shgo_adaptive
    gfortran "${flags[@]}" "$ROOT/test/optimize_shgo_hgr_reference.f90" \
        kinds.o shgo.o -o native_shgo_hgr
    ./native_shgo_hgr
    gfortran "${flags[@]}" "$ROOT/test/optimize_shgo_conforming_reference.f90" \
        kinds.o shgo.o -o native_shgo_conforming
    ./native_shgo_conforming
    gfortran "${flags[@]}" "$ROOT/validation/optimize_shgo_driver.f90" \
        kinds.o shgo.o -o driver_shgo
    python3 "$ROOT/validation/compare_optimize_shgo_scipy.py" ./driver_shgo
    gfortran "${flags[@]}" "$ROOT/validation/optimize_shgo_adaptive_driver.f90" \
        kinds.o shgo.o -o driver_shgo_adaptive
    python3 "$ROOT/validation/compare_optimize_shgo_adaptive_scipy.py" ./driver_shgo_adaptive
    gfortran "${flags[@]}" "$ROOT/validation/optimize_shgo_hgr_driver.f90" \
        kinds.o shgo.o -o driver_shgo_hgr
    python3 "$ROOT/validation/compare_optimize_shgo_hgr_scipy.py" ./driver_shgo_hgr
    gfortran "${flags[@]}" "$ROOT/validation/optimize_shgo_conforming_driver.f90" \
        kinds.o shgo.o -o driver_shgo_conforming
    python3 "$ROOT/validation/compare_optimize_shgo_conforming_scipy.py" ./driver_shgo_conforming
    gfortran "${flags[@]}" "$ROOT/test/optimize_shgo_conforming_3d_reference.f90" \
        kinds.o shgo.o -o native_shgo_conforming_3d
    ./native_shgo_conforming_3d
    gfortran "${flags[@]}" "$ROOT/validation/optimize_shgo_conforming_3d_driver.f90" \
        kinds.o shgo.o -o driver_shgo_conforming_3d
    python3 "$ROOT/validation/compare_optimize_shgo_conforming_3d_scipy.py" ./driver_shgo_conforming_3d
    gfortran "${flags[@]}" "$ROOT/example/optimize_shgo_conforming_3d_example.f90" \
        kinds.o shgo.o -o example_shgo_conforming_3d
    ./example_shgo_conforming_3d
    gfortran "${flags[@]}" "$ROOT/test/optimize_shgo_constraints_reference.f90" \
        kinds.o shgo.o -o native_shgo_constraints
    ./native_shgo_constraints
    gfortran "${flags[@]}" "$ROOT/validation/optimize_shgo_constraints_driver.f90" \
        kinds.o shgo.o -o driver_shgo_constraints
    python3 "$ROOT/validation/compare_optimize_shgo_constraints_scipy.py" ./driver_shgo_constraints
    gfortran "${flags[@]}" "$ROOT/test/optimize_shgo_restoration_reference.f90" \
        kinds.o shgo.o -o native_shgo_restoration
    ./native_shgo_restoration
    gfortran "${flags[@]}" "$ROOT/validation/optimize_shgo_restoration_driver.f90" \
        kinds.o shgo.o -o driver_shgo_restoration
    python3 "$ROOT/validation/compare_optimize_shgo_restoration_scipy.py" ./driver_shgo_restoration
    gfortran "${flags[@]}" "$ROOT/test/optimize_shgo_multiconstraints_reference.f90" \
        kinds.o shgo.o -o native_shgo_multiconstraints
    ./native_shgo_multiconstraints
    gfortran "${flags[@]}" "$ROOT/validation/optimize_shgo_multiconstraints_driver.f90" \
        kinds.o shgo.o -o driver_shgo_multiconstraints
    python3 "$ROOT/validation/compare_optimize_shgo_multiconstraints_scipy.py" ./driver_shgo_multiconstraints
    # Use a test-only stand-in for SciFort's FPM-dependent optimizer module.
    gfortran "${flags[@]}" -c \
        "$ROOT/validation/optimize_shgo_test_minimize_stub.f90" \
        "$ROOT/src/optimize/shgo_minimize.f90"
    gfortran "${flags[@]}" "$ROOT/validation/optimize_shgo_adapter_smoke.f90" \
        kinds.o shgo.o optimize_shgo_test_minimize_stub.o shgo_minimize.o \
        -o shgo_adapter_smoke
    ./shgo_adapter_smoke
    gfortran "${flags[@]}" "$ROOT/example/optimize_shgo_constraints_example.f90" \
        kinds.o shgo.o optimize_shgo_test_minimize_stub.o shgo_minimize.o \
        -o shgo_constrained_example
    ./shgo_constrained_example
    gfortran "${flags[@]}" "$ROOT/example/optimize_shgo_restoration_example.f90" \
        kinds.o shgo.o optimize_shgo_test_minimize_stub.o shgo_minimize.o \
        -o shgo_restoration_example
    ./shgo_restoration_example
    gfortran "${flags[@]}" "$ROOT/example/optimize_shgo_multiconstraints_example.f90" \
        kinds.o shgo.o optimize_shgo_test_minimize_stub.o shgo_minimize.o \
        -o example_shgo_multiconstraints
    ./example_shgo_multiconstraints
    echo "$mode: SHGO core verification passed"
done
# To validate public adapter and full library dependencies, run: fpm test
