/* SPDX-License-Identifier: MIT */
/* Copyright (c) 2026 SciFort contributors */
#ifndef SCIFORT_RNG_H
#define SCIFORT_RNG_H
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
/* MRG32k3a state is six caller-owned components in native get/set order.
 * First triple: [0,4294967087); second: [0,4294944443); neither all zero.
 * All pointers must be valid, aligned, disjoint and fully sized, even for n=0.
 * Seed accepts every int64_t value. Uniform draws match the native stream.
 * Uniform status: 0 success, 1 negative n, 3 invalid state. Errors preserve
 * state/output. Empty calls validate state and do not advance valid streams.
 * Independent states are independent objects; do not concurrently mutate one. */
void scifort_rng_seed_i64(int64_t seed, int64_t values[6], int *status);
void scifort_rng_uniform_f64(int n, int64_t values[6], double *output, int *status);
#ifdef __cplusplus
}
#endif
#endif
