#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
# Reproduces standalone dual-annealing checks. Not the full SciFort/FPM suite.
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
        "$ROOT/src/kinds.f90" "$ROOT/src/optimize/dual_annealing.f90"
    gfortran -std=f2008 "${flags[@]}" \
        "$ROOT/test/optimize_dual_annealing_reference.f90" kinds.o dual_annealing.o \
        -o native
    ./native
    gfortran -std=f2008 "${flags[@]}" \
        "$ROOT/validation/optimize_dual_annealing_driver.f90" kinds.o dual_annealing.o \
        -o driver
    python3 "$ROOT/validation/compare_optimize_dual_annealing_scipy.py" ./driver
    gfortran -std=f2008 "${flags[@]}" -c \
        "$ROOT/validation/optimize_dual_annealing_test_minimize.f90" \
        "$ROOT/src/optimize/dual_annealing_minimize.f90"
    gfortran -std=f2008 "${flags[@]}" \
        "$ROOT/test/optimize_dual_annealing_public_reference.f90" kinds.o dual_annealing.o \
        optimize_dual_annealing_test_minimize.o dual_annealing_minimize.o \
        -o public_contract
    ./public_contract
    gfortran -std=f2008 "${flags[@]}" \
        "$ROOT/example/optimize_dual_annealing_example.f90" kinds.o dual_annealing.o \
        -o example
    ./example
    echo "${mode} dual annealing validation passed"
done
