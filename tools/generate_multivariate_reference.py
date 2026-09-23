#!/usr/bin/env python3
"""Generate SciPy 1.17 references for SciFort multivariate distributions."""
from pathlib import Path
import warnings
import numpy as np
import scipy
from scipy import stats

OUT = Path(__file__).resolve().parents[1] / "test" / "multivariate_reference.f90"


def fmt(v):
    v = float(v)
    if np.isnan(v):
        return "transfer(9221120237041090560_int64, 0.0_dp)"
    if np.isposinf(v):
        return "transfer(9218868437227405312_int64, 0.0_dp)"
    if np.isneginf(v):
        return "transfer(-4503599627370496_int64, 0.0_dp)"
    return f"{v:.17e}_dp"


def rscalar(name, value):
    return f"    real(dp), parameter :: {name} = {fmt(value)}\n"


def rvec(name, values):
    a = np.asarray(values, dtype=float).reshape(-1)
    lines = [f"    real(dp), parameter :: {name}({a.size}) = [ &\n"]
    for i, v in enumerate(a):
        suffix = ", &" if i < a.size - 1 else " ]"
        lines.append(f"        {fmt(v)}{suffix}\n")
    return "".join(lines)


def rmat(name, values):
    a = np.asarray(values, dtype=float)
    flat = a.flatten(order="F")
    lines = [f"    real(dp), parameter :: {name}({a.shape[0]},{a.shape[1]}) = reshape([ &\n"]
    for i, v in enumerate(flat):
        suffix = ", &" if i < flat.size - 1 else " ], &"
        lines.append(f"        {fmt(v)}{suffix}\n")
    lines.append(f"        [{a.shape[0]},{a.shape[1]}])\n")
    return "".join(lines)


# Multivariate normal references.
mvn1_x = np.array([0.2, -0.4])
mvn1_mean = np.array([0.5, -1.0])
mvn1_cov = np.array([[2.0, 0.4], [0.4, 1.5]])

mvn2_x = np.array([1.2, 0.1, -0.7])
mvn2_mean = np.array([0.0, 0.3, -1.0])
mvn2_cov = np.array([[1.0, 0.2, -0.1], [0.2, 2.0, 0.5], [-0.1, 0.5, 1.2]])

mvn_lower_cov = np.array([[2.0, 99.0], [0.4, 1.5]])
mvn_singular_cov = np.array([[1.0, 1.0], [1.0, 1.0]])
mvn_singular_on = np.array([0.5, 0.5])
mvn_singular_off = np.array([0.5, 0.4])

fit_data = np.array([
    [0.2, -1.1],
    [1.4, 0.3],
    [-0.7, 0.8],
    [2.0, -0.2],
    [0.1, 1.7],
])
fit_mean, fit_cov = stats.multivariate_normal.fit(fit_data)
fixed_mean = np.array([0.25, -0.5])
fit_fixed_mean, fit_cov_fixed_mean = stats.multivariate_normal.fit(fit_data, fix_mean=fixed_mean)

mvn_box_lower = np.array([-1.0, -2.0])
mvn1_cdf = stats.multivariate_normal.cdf(
    mvn1_x, mvn1_mean, mvn1_cov, rng=np.random.default_rng(123),
    abseps=1e-10, maxpts=1_000_000)
mvn1_logcdf = stats.multivariate_normal.logcdf(
    mvn1_x, mvn1_mean, mvn1_cov, rng=np.random.default_rng(123),
    abseps=1e-10, maxpts=1_000_000)
mvn1_box_cdf = stats.multivariate_normal.cdf(
    mvn1_x, mvn1_mean, mvn1_cov, lower_limit=mvn_box_lower,
    rng=np.random.default_rng(123), abseps=1e-10, maxpts=1_000_000)
mvn2_cdf = stats.multivariate_normal.cdf(
    mvn2_x, mvn2_mean, mvn2_cov, rng=np.random.default_rng(123),
    abseps=1e-9, maxpts=3_000_000)
mvn_singular_on_cdf = stats.multivariate_normal.cdf(
    mvn_singular_on, cov=mvn_singular_cov, allow_singular=True,
    rng=np.random.default_rng(123), abseps=1e-10, maxpts=1_000_000)
mvn_singular_off_cdf = stats.multivariate_normal.cdf(
    mvn_singular_off, cov=mvn_singular_cov, allow_singular=True,
    rng=np.random.default_rng(123), abseps=1e-10, maxpts=1_000_000)
mvn_marg_dims = np.array([2, 0], dtype=int)
mvn_marg_mean = mvn2_mean[mvn_marg_dims]
mvn_marg_cov = mvn2_cov[np.ix_(mvn_marg_dims, mvn_marg_dims)]

# Multivariate Student t references.
mvt_x = np.array([0.4, 5.0])
mvt_loc = np.array([0.0, 1.0])
mvt_shape = np.array([[1.0, 0.1], [0.1, 1.0]])
mvt_df = 7.0
mvt_lower = np.array([-1.0, 0.0])
mvt_cdf = stats.multivariate_t.cdf(
    mvt_x, mvt_loc, mvt_shape, mvt_df, maxpts=1_000_000,
    random_state=np.random.default_rng(123))
mvt_box_cdf = stats.multivariate_t.cdf(
    mvt_x, mvt_loc, mvt_shape, mvt_df, lower_limit=mvt_lower, maxpts=1_000_000,
    random_state=np.random.default_rng(123))
mvt_singular_shape = np.array([[1.0, 1.0], [1.0, 1.0]])
mvt_singular_x = np.array([0.5, 0.4])
mvt_marg_dims = np.array([1, 0], dtype=int)
mvt_marg_loc = mvt_loc[mvt_marg_dims]
mvt_marg_shape = mvt_shape[np.ix_(mvt_marg_dims, mvt_marg_dims)]

# Four-dimensional CDF cases exercise the full permuted-Cholesky pivot path.
mv4_x = np.array([0.2, -0.3, 0.7, 0.1])
mv4_loc = np.array([0.1, -0.2, 0.0, 0.3])
mv4_cov = np.array([
    [1.0, 0.5, -0.2, 0.1],
    [0.5, 2.0, 0.3, -0.4],
    [-0.2, 0.3, 1.5, 0.25],
    [0.1, -0.4, 0.25, 0.8],
])
mv4_lower = np.array([-1.2, -1.0, -0.8, -1.5])
mv4_cdf = stats.multivariate_normal.cdf(
    mv4_x, mv4_loc, mv4_cov, rng=np.random.default_rng(321),
    abseps=5e-8, maxpts=5_000_000)
mv4_box_cdf = stats.multivariate_normal.cdf(
    mv4_x, mv4_loc, mv4_cov, lower_limit=mv4_lower,
    rng=np.random.default_rng(321), abseps=5e-8, maxpts=5_000_000)
mvt4_cdf = stats.multivariate_t.cdf(
    mv4_x, mv4_loc, mv4_cov, 6.5, maxpts=2_000_000,
    random_state=np.random.default_rng(321))
mvt4_box_cdf = stats.multivariate_t.cdf(
    mv4_x, mv4_loc, mv4_cov, 6.5, lower_limit=mv4_lower, maxpts=2_000_000,
    random_state=np.random.default_rng(321))

# Dirichlet references.
dir_alpha = np.array([0.4, 5.0, 15.0])
dir_x = np.array([0.2, 0.2, 0.6])
dir_x_short = dir_x[:-1]

# Multinomial references.
mn_n = 8
mn_p = np.array([0.3, 0.2, 0.5])
mn_x = np.array([1, 3, 4])
mn_adjust_n = 6
mn_adjust_p = np.array([0.2, 0.2, 0.2])
mn_adjust_x = np.array([1, 1, 4])
with warnings.catch_warnings():
    warnings.simplefilter("ignore", FutureWarning)
    mn_adjust_logpmf = stats.multinomial.logpmf(mn_adjust_x, mn_adjust_n, mn_adjust_p)
    mn_adjust_pmf = stats.multinomial.pmf(mn_adjust_x, mn_adjust_n, mn_adjust_p)
    mn_adjust_mean = stats.multinomial.mean(mn_adjust_n, mn_adjust_p)

# Dirichlet multinomial references.
dm_alpha = np.array([3.0, 4.0, 5.0])
dm_n = 6
dm_x = np.array([1, 2, 3])

parts = [
    "! SPDX-License-Identifier: MIT\n",
    "! Generated by tools/generate_multivariate_reference.py\n",
    f"! SciPy version: {scipy.__version__}\n",
    "module test_multivariate_reference\n",
    "    use, intrinsic :: iso_fortran_env, only : int64\n",
    "    use scifort_kinds, only : dp\n",
    "    implicit none\n",
    rvec("mvn1_x_ref", mvn1_x), rvec("mvn1_mean_ref", mvn1_mean), rmat("mvn1_cov_ref", mvn1_cov),
    rscalar("mvn1_pdf_ref", stats.multivariate_normal.pdf(mvn1_x, mvn1_mean, mvn1_cov)),
    rscalar("mvn1_logpdf_ref", stats.multivariate_normal.logpdf(mvn1_x, mvn1_mean, mvn1_cov)),
    rscalar("mvn1_entropy_ref", stats.multivariate_normal(mean=mvn1_mean, cov=mvn1_cov).entropy()),
    rvec("mvn2_x_ref", mvn2_x), rvec("mvn2_mean_ref", mvn2_mean), rmat("mvn2_cov_ref", mvn2_cov),
    rscalar("mvn2_pdf_ref", stats.multivariate_normal.pdf(mvn2_x, mvn2_mean, mvn2_cov)),
    rscalar("mvn2_logpdf_ref", stats.multivariate_normal.logpdf(mvn2_x, mvn2_mean, mvn2_cov)),
    rscalar("mvn2_entropy_ref", stats.multivariate_normal(mean=mvn2_mean, cov=mvn2_cov).entropy()),
    rmat("mvn_lower_cov_ref", mvn_lower_cov),
    rscalar("mvn_lower_pdf_ref", stats.multivariate_normal.pdf(mvn1_x, mvn1_mean, mvn_lower_cov)),
    rmat("mvn_singular_cov_ref", mvn_singular_cov), rvec("mvn_singular_on_ref", mvn_singular_on),
    rvec("mvn_singular_off_ref", mvn_singular_off),
    rscalar("mvn_singular_on_pdf_ref", stats.multivariate_normal.pdf(mvn_singular_on, cov=mvn_singular_cov, allow_singular=True)),
    rscalar(
        "mvn_singular_on_logpdf_ref",
        stats.multivariate_normal.logpdf(
            mvn_singular_on, cov=mvn_singular_cov, allow_singular=True
        ),
    ),
    rscalar("mvn_singular_off_pdf_ref", stats.multivariate_normal.pdf(mvn_singular_off, cov=mvn_singular_cov, allow_singular=True)),
    rmat("mvn_fit_data_ref", fit_data.T), rvec("mvn_fit_mean_ref", fit_mean), rmat("mvn_fit_cov_ref", fit_cov),
    rvec("mvn_fixed_mean_ref", fixed_mean), rmat("mvn_fit_cov_fixed_mean_ref", fit_cov_fixed_mean),
    rscalar("mvn1_cdf_ref", mvn1_cdf), rscalar("mvn1_logcdf_ref", mvn1_logcdf),
    rvec("mvn_box_lower_ref", mvn_box_lower), rscalar("mvn1_box_cdf_ref", mvn1_box_cdf),
    rscalar("mvn2_cdf_ref", mvn2_cdf),
    rscalar("mvn_singular_on_cdf_ref", mvn_singular_on_cdf),
    rscalar("mvn_singular_off_cdf_ref", mvn_singular_off_cdf),
    rvec("mvn_marg_mean_ref", mvn_marg_mean), rmat("mvn_marg_cov_ref", mvn_marg_cov),
    rvec("mvt_x_ref", mvt_x), rvec("mvt_loc_ref", mvt_loc), rmat("mvt_shape_ref", mvt_shape),
    rscalar("mvt_df_ref", mvt_df),
    rscalar("mvt_pdf_ref", stats.multivariate_t.pdf(mvt_x, mvt_loc, mvt_shape, mvt_df)),
    rscalar("mvt_logpdf_ref", stats.multivariate_t.logpdf(mvt_x, mvt_loc, mvt_shape, mvt_df)),
    rscalar("mvt_entropy_ref", stats.multivariate_t.entropy(mvt_loc, mvt_shape, mvt_df)),
    rscalar("mvt_cdf_ref", mvt_cdf), rvec("mvt_lower_ref", mvt_lower),
    rscalar("mvt_box_cdf_ref", mvt_box_cdf),
    rmat("mvt_singular_shape_ref", mvt_singular_shape), rvec("mvt_singular_x_ref", mvt_singular_x),
    rscalar("mvt_singular_pdf_ref", stats.multivariate_t.pdf(mvt_singular_x, shape=mvt_singular_shape, df=4, allow_singular=True)),
    rscalar("mvt_singular_logpdf_ref", stats.multivariate_t.logpdf(mvt_singular_x, shape=mvt_singular_shape, df=4)),
    rvec("mvt_marg_loc_ref", mvt_marg_loc), rmat("mvt_marg_shape_ref", mvt_marg_shape),
    rvec("mv4_x_ref", mv4_x), rvec("mv4_loc_ref", mv4_loc), rmat("mv4_cov_ref", mv4_cov),
    rvec("mv4_lower_ref", mv4_lower), rscalar("mv4_cdf_ref", mv4_cdf),
    rscalar("mv4_box_cdf_ref", mv4_box_cdf), rscalar("mvt4_cdf_ref", mvt4_cdf),
    rscalar("mvt4_box_cdf_ref", mvt4_box_cdf),
    rvec("dir_alpha_ref", dir_alpha), rvec("dir_x_ref", dir_x), rvec("dir_x_short_ref", dir_x_short),
    rscalar("dir_pdf_ref", stats.dirichlet.pdf(dir_x, dir_alpha)),
    rscalar("dir_logpdf_ref", stats.dirichlet.logpdf(dir_x, dir_alpha)),
    rvec("dir_mean_ref", stats.dirichlet.mean(dir_alpha)),
    rvec("dir_var_ref", stats.dirichlet.var(dir_alpha)),
    rmat("dir_cov_ref", stats.dirichlet.cov(dir_alpha)),
    rscalar("dir_entropy_ref", stats.dirichlet.entropy(dir_alpha)),
    rvec("mn_p_ref", mn_p), rvec("mn_x_ref", mn_x),
    rscalar("mn_pmf_ref", stats.multinomial.pmf(mn_x, mn_n, mn_p)),
    rscalar("mn_logpmf_ref", stats.multinomial.logpmf(mn_x, mn_n, mn_p)),
    rvec("mn_mean_ref", stats.multinomial.mean(mn_n, mn_p)),
    rmat("mn_cov_ref", stats.multinomial.cov(mn_n, mn_p)),
    rscalar("mn_entropy_ref", stats.multinomial.entropy(mn_n, mn_p)),
    rvec("mn_adjust_p_ref", mn_adjust_p), rvec("mn_adjust_x_ref", mn_adjust_x),
    rscalar("mn_adjust_pmf_ref", mn_adjust_pmf), rscalar("mn_adjust_logpmf_ref", mn_adjust_logpmf),
    rvec("mn_adjust_mean_ref", mn_adjust_mean),
    rvec("dm_alpha_ref", dm_alpha), rvec("dm_x_ref", dm_x),
    rscalar("dm_pmf_ref", stats.dirichlet_multinomial.pmf(dm_x, dm_alpha, dm_n)),
    rscalar("dm_logpmf_ref", stats.dirichlet_multinomial.logpmf(dm_x, dm_alpha, dm_n)),
    rvec("dm_mean_ref", stats.dirichlet_multinomial.mean(dm_alpha, dm_n)),
    rvec("dm_var_ref", stats.dirichlet_multinomial.var(dm_alpha, dm_n)),
    rmat("dm_cov_ref", stats.dirichlet_multinomial.cov(dm_alpha, dm_n)),
    "end module test_multivariate_reference\n",
]
OUT.write_text("".join(parts), encoding="ascii")
print(OUT)
