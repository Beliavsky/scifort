/* SPDX-License-Identifier: MIT */
/* Copyright (c) 2026 SciFort contributors */

#ifndef SCIFORT_H
#define SCIFORT_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

enum {
    SCIFORT_STATUS_OK = 0,
    SCIFORT_STATUS_INVALID_ARGUMENT = 1
};

int scifort_version_major(void);
int scifort_version_minor(void);
int scifort_version_patch(void);

double scifort_normal_pdf_f64(double x, double loc, double scale);
double scifort_normal_cdf_f64(double x, double loc, double scale);
double scifort_normal_ppf_f64(double p, double loc, double scale);

void scifort_normal_pdf_vec_f64(
    size_t n,
    const double *x,
    double loc,
    double scale,
    double *y,
    int *status
);

void scifort_normal_cdf_vec_f64(
    size_t n,
    const double *x,
    double loc,
    double scale,
    double *y,
    int *status
);

void scifort_normal_ppf_vec_f64(
    size_t n,
    const double *p,
    double loc,
    double scale,
    double *x,
    int *status
);

#ifdef __cplusplus
}
#endif

#endif
