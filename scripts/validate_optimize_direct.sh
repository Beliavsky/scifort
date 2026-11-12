#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
# Standalone checked/optimized native and SciPy validation for DIRECT.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEMP="$(mktemp -d)"
trap 'rm -rf "$TEMP"' EXIT
for mode in checked optimized; do
  mkdir -p "$TEMP/$mode"
  cd "$TEMP/$mode"
  if [[ "$mode" == checked ]]; then
    flags=(-std=f2008 -O0 -g -Wall -Wextra -fcheck=all -ffpe-trap=invalid,zero,overflow)
  else
    flags=(-std=f2008 -O3 -Wall -Wextra)
  fi
  gfortran "${flags[@]}" -c \
      "$ROOT/src/kinds.f90" "$ROOT/src/optimize/direct.f90"
  gfortran "${flags[@]}" "$ROOT/test/optimize_direct_reference.f90" \
      kinds.o direct.o -o native_direct
  ./native_direct
  gfortran "${flags[@]}" "$ROOT/validation/optimize_direct_driver.f90" \
      kinds.o direct.o -o driver_direct
  python3 "$ROOT/validation/compare_optimize_direct_scipy.py" ./driver_direct
  gfortran "${flags[@]}" "$ROOT/example/optimize_direct_example.f90" \
      kinds.o direct.o -o example_direct
  ./example_direct
  echo "$mode: DIRECT validation passed"
done
