#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
# Standalone numerical test for new basin-hopping engine, not the full FPM suite.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEMP="$(mktemp -d)"
trap 'rm -rf "$TEMP"' EXIT
for mode in checked optimized; do
    mkdir -p "$TEMP/$mode"
    cd "$TEMP/$mode"
    if [[ "$mode" == checked ]]; then
        flags=(-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow)
    else
        flags=(-O3 -march=native)
    fi
    gfortran -std=f2008 "${flags[@]}" -c \
        "$ROOT/src/kinds.f90" "$ROOT/src/optimize/basinhopping.f90"
    gfortran -std=f2008 "${flags[@]}" \
        "$ROOT/test/optimize_basinhopping_reference.f90" kinds.o basinhopping.o -o native
    ./native
    gfortran -std=f2008 "${flags[@]}" \
        "$ROOT/validation/optimize_basinhopping_driver.f90" kinds.o basinhopping.o -o driver
    python3 "$ROOT/validation/compare_optimize_basinhopping_scipy.py" ./driver
    # This substitute tests adapter signatures, NOT SciFort's real local algorithms.
    gfortran -std=f2008 "${flags[@]}" -c \
        "$ROOT/validation/optimize_basinhopping_test_minimize.f90" \
        "$ROOT/src/optimize/basinhopping_minimize.f90"
    gfortran -std=f2008 "${flags[@]}" \
        "$ROOT/test/optimize_basinhopping_public_reference.f90" \
        kinds.o basinhopping.o optimize_basinhopping_test_minimize.o \
        basinhopping_minimize.o -o public_contract
    ./public_contract
    printf '%s validation passed\n' "$mode"
done
