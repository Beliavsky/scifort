# SPDX-License-Identifier: MIT
# Copyright (c) 2026 SciFort contributors
"""Print reference values used by descriptive and hypothesis tests.

This is a developer validation helper, not a runtime dependency. It records
NumPy/SciPy versions and reproduces the reference values embedded in
``test/test_descriptive.f90`` and ``test/test_hypothesis.f90``.
It also reproduces the SciPy-compatibility values in
``test/test_descriptive_scipy_compat.f90``,
``test/test_descriptive_axis_keepdims.f90``, and
``test/test_descriptive_nan_policy.f90``,
``test/test_descriptive_rank3.f90``, and
``test/test_descriptive_rank3_weights.f90``, and
``test/test_descriptive_rank4.f90``,
``test/test_descriptive_highrank_means.f90``,
``test/test_descriptive_highrank_weights.f90``, and
``test/test_descriptive_highrank_stats.f90``, and
``test/test_descriptive_named_scale.f90``, and
``test/test_descriptive_standardize.f90``,
``test/test_descriptive_trimmed.f90``, and
``test/test_descriptive_proportion_trim.f90``, and
``test/test_descriptive_kstat.f90``, and
``test/test_descriptive_circular.f90``, and
``test/test_descriptive_frequency.f90``, and
``test/test_descriptive_summary.f90``, and
``test/test_descriptive_expectile_sigmaclip.f90``, and
``test/test_descriptive_lmoment.f90``,
``test/test_power_transform.f90``, and
``test/test_probability_plot.f90``, and
``test/test_hypothesis_extended.f90``.
"""

import numpy as np
import scipy
from scipy import stats


def main():
    x = np.array([1.0, 2.0, 2.0, 4.0, 8.0, 16.0])
    y = np.array([3.0, -1.0, 5.0, 4.0, 9.0, 12.0])
    weights = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    q = np.array([0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0])

    print(f"numpy {np.__version__}")
    print(f"scipy {scipy.__version__}")
    print("mean", np.mean(x))
    print("variance_ddof0", np.var(x, ddof=0))
    print("variance_ddof1", np.var(x, ddof=1))
    print("standard_deviation_ddof1", np.std(x, ddof=1))
    for order in range(5):
        print(f"central_moment_{order}", stats.moment(x, order=order))
    print("quantile_linear", np.quantile(x, q, method="linear"))
    print("covariance_ddof0", np.cov(x, y, ddof=0)[0, 1])
    print("covariance_ddof1", np.cov(x, y, ddof=1)[0, 1])
    print("pearson_correlation", np.corrcoef(x, y)[0, 1])
    for method in ("average", "min", "max", "dense", "ordinal"):
        print(f"rankdata_{method}", stats.rankdata(x, method=method))

    print("gmean", stats.gmean(x))
    print("gmean_weighted", stats.gmean(x, weights=weights))
    print("hmean", stats.hmean(x))
    print("hmean_weighted", stats.hmean(x, weights=weights))
    for power in (-2.0, 0.0, 0.5, 1.0, 2.0, 3.5):
        print(f"pmean_{power}", stats.pmean(x, power))
    print("pmean_weighted_2", stats.pmean(x, 2.0, weights=weights))
    for ddof in (0, 1, 2):
        print(f"sem_ddof{ddof}", stats.sem(x, ddof=ddof))
        print(f"variation_ddof{ddof}", stats.variation(x, ddof=ddof))
    for percentile_range in ((25, 75), (10, 90), (75, 25), (0, 100)):
        print(
            f"iqr_{percentile_range[0]}_{percentile_range[1]}",
            stats.iqr(x, rng=percentile_range),
        )
    print("iqr_scaled", stats.iqr(x, scale=1.3489795003921634))
    print("iqr_scale_normal", stats.iqr(x, scale="normal"))
    print("mad_scale_normal", stats.median_abs_deviation(x, scale="normal"))
    for scale in (1.0, 2.0, 1.482602218505602):
        print(f"median_abs_deviation_{scale}", stats.median_abs_deviation(x, scale=scale))
    for order in range(6):
        print(f"moment_{order}", stats.moment(x, order=order))
    for center in (0.0, 3.0, 5.5):
        print(f"moment_2_center_{center}", stats.moment(x, order=2, center=center))
    print("skew_biased", stats.skew(x, bias=True))
    print("skew_unbiased", stats.skew(x, bias=False))
    print("kurtosis_fisher_biased", stats.kurtosis(x, fisher=True, bias=True))
    print("kurtosis_fisher_unbiased", stats.kurtosis(x, fisher=True, bias=False))
    print("kurtosis_pearson_biased", stats.kurtosis(x, fisher=False, bias=True))
    print("kurtosis_pearson_unbiased", stats.kurtosis(x, fisher=False, bias=False))

    matrix = np.array(
        [[1.0, 2.0, 4.0, 8.0], [3.0, 5.0, 7.0, 9.0], [2.0, 6.0, 10.0, 14.0]]
    )
    weights4 = np.array([1.0, 2.0, 3.0, 4.0])
    weights31 = np.array([[1.0], [2.0], [4.0]])
    for axis in (0, 1):
        print(f"matrix_gmean_axis{axis}", stats.gmean(matrix, axis=axis))
        print(f"matrix_hmean_axis{axis}", stats.hmean(matrix, axis=axis))
        print(f"matrix_pmean_axis{axis}", stats.pmean(matrix, 2.0, axis=axis))
        print(f"matrix_sem_axis{axis}", stats.sem(matrix, axis=axis))
        print(f"matrix_variation_axis{axis}", stats.variation(matrix, axis=axis))
        print(f"matrix_iqr_axis{axis}", stats.iqr(matrix, axis=axis))
        print(
            f"matrix_mad_axis{axis}",
            stats.median_abs_deviation(matrix, axis=axis),
        )
        print(f"matrix_moment_axis{axis}", stats.moment(matrix, order=2, axis=axis))
        print(f"matrix_skew_axis{axis}", stats.skew(matrix, axis=axis))
        print(f"matrix_kurtosis_axis{axis}", stats.kurtosis(matrix, axis=axis))
    print("matrix_iqr_all", stats.iqr(matrix))
    print("matrix_gmean_axis0_weights4", stats.gmean(matrix, axis=0, weights=weights4))
    print("matrix_gmean_axis1_weights4", stats.gmean(matrix, axis=1, weights=weights4))
    print("matrix_gmean_axis0_weights31", stats.gmean(matrix, axis=0, weights=weights31))
    print("matrix_gmean_keepdims", stats.gmean(matrix, axis=0, keepdims=True))
    print("matrix_iqr_keepdims_all", stats.iqr(matrix, keepdims=True))

    standardize_matrix = np.array([1.0, 2.0, 4.0, 8.0, 16.0, 32.0]).reshape(
        (3, 2), order="F"
    )
    for axis in (0, 1):
        print(f"standardize_zscore_axis{axis}", stats.zscore(standardize_matrix, axis=axis))
        print(f"standardize_gzscore_axis{axis}", stats.gzscore(standardize_matrix, axis=axis))
        print(f"standardize_gstd_axis{axis}", stats.gstd(standardize_matrix, axis=axis))
    print("standardize_zscore_all", stats.zscore(standardize_matrix, axis=None))
    print("standardize_gzscore_all", stats.gzscore(standardize_matrix, axis=None))
    print("standardize_gstd_all", stats.gstd(standardize_matrix, axis=None))
    print("standardize_gstd_all_ddof0", stats.gstd(standardize_matrix, axis=None, ddof=0))
    zmap_scores = np.array([0.5, 2.0, 2.5, 3.0])
    zmap_compare = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    print("standardize_zmap_vector", stats.zmap(zmap_scores, zmap_compare))
    print("standardize_zmap_ddof1", stats.zmap(zmap_scores, zmap_compare, ddof=1))
    zmap_scores_2d = np.array([[1.0], [2.0], [3.0]])
    zmap_compare_2d = np.array(
        [[0.0, 10.0], [2.0, 12.0], [4.0, 14.0], [6.0, 16.0]]
    )
    print(
        "standardize_zmap_broadcast",
        stats.zmap(zmap_scores_2d, zmap_compare_2d, axis=0),
    )
    print(
        "standardize_zmap_all",
        stats.zmap(zmap_scores_2d, zmap_compare_2d, axis=None),
    )
    print("standardize_zscore_ddof1", stats.zscore(standardize_matrix, axis=0, ddof=1))
    print("standardize_gstd_ddof0", stats.gstd(standardize_matrix, axis=0, ddof=0))
    print(
        "standardize_gstd_keepdims_axis1",
        stats.gstd(standardize_matrix, axis=1, keepdims=True),
    )
    standardize_cube = np.arange(1.0, 25.0).reshape((2, 3, 4), order="F")
    print(
        "standardize_gstd_tuple_02_keepdims",
        stats.gstd(standardize_cube, axis=(0, 2), keepdims=True),
    )
    print(
        "standardize_gstd_tuple_all_keepdims",
        stats.gstd(standardize_cube, axis=(0, 1, 2), keepdims=True),
    )
    standardize_nan = standardize_matrix.copy()
    standardize_nan[1, 0] = np.nan
    standardize_cube_nan = standardize_cube.copy()
    standardize_cube_nan[1, 1, 2] = np.nan
    zmap_compare_nan = zmap_compare_2d.copy()
    zmap_compare_nan[2, 0] = np.nan
    for policy in ("propagate", "omit"):
        print(
            f"standardize_gstd_tuple_02_nan_{policy}",
            stats.gstd(
                standardize_cube_nan,
                axis=(0, 2),
                keepdims=True,
                nan_policy=policy,
            ),
        )
        print(
            f"standardize_zscore_all_nan_{policy}",
            stats.zscore(standardize_nan, axis=None, nan_policy=policy),
        )
        print(
            f"standardize_gzscore_all_nan_{policy}",
            stats.gzscore(standardize_nan, axis=None, nan_policy=policy),
        )
        print(
            f"standardize_gstd_all_nan_{policy}",
            stats.gstd(standardize_nan, axis=None, nan_policy=policy),
        )
        print(
            f"standardize_zmap_all_nan_{policy}",
            stats.zmap(zmap_scores_2d, zmap_compare_nan, axis=None, nan_policy=policy),
        )
        print(
            f"standardize_zscore_nan_{policy}",
            stats.zscore(standardize_nan, axis=0, nan_policy=policy),
        )
        print(
            f"standardize_gzscore_nan_{policy}",
            stats.gzscore(standardize_nan, axis=0, nan_policy=policy),
        )
        print(
            f"standardize_gstd_nan_{policy}",
            stats.gstd(standardize_nan, axis=0, nan_policy=policy),
        )
    print(
        "standardize_gzscore_negative_omit",
        stats.gzscore(np.array([1.0, -2.0, 4.0, 8.0]), nan_policy="omit"),
    )
    print("standardize_zscore_ddof_minus1", stats.zscore(np.array([1.0, 2.0, 4.0]), ddof=-1))
    print("standardize_gstd_ddof_minus1", stats.gstd(np.array([1.0, 2.0, 4.0]), ddof=-1))

    trimmed_x = np.arange(20.0)
    print("trimmed_tmean", stats.tmean(trimmed_x, limits=(3.0, 17.0)))
    print("trimmed_tvar", stats.tvar(trimmed_x, limits=(3.0, 17.0)))
    print("trimmed_tstd", stats.tstd(trimmed_x, limits=(3.0, 17.0)))
    print("trimmed_tsem", stats.tsem(trimmed_x, limits=(3.0, 17.0)))
    print(
        "trimmed_tmean_upper_exclusive",
        stats.tmean(trimmed_x, limits=(3.0, 17.0), inclusive=(True, False)),
    )
    print("trimmed_tmin", stats.tmin(trimmed_x, lowerlimit=13.0))
    print(
        "trimmed_tmin_exclusive",
        stats.tmin(trimmed_x, lowerlimit=13.0, inclusive=False),
    )
    print("trimmed_tmax", stats.tmax(trimmed_x, upperlimit=13.0))
    print(
        "trimmed_tmax_exclusive",
        stats.tmax(trimmed_x, upperlimit=13.0, inclusive=False),
    )
    trimmed_matrix = np.arange(1.0, 13.0).reshape((3, 4), order="F")
    print(
        "trimmed_tmean_axis0",
        stats.tmean(trimmed_matrix, limits=(3.0, 10.0), axis=0),
    )
    print(
        "trimmed_tmean_axis1",
        stats.tmean(trimmed_matrix, limits=(3.0, 10.0), axis=1),
    )
    print(
        "trimmed_tvar_axis0",
        stats.tvar(trimmed_matrix, limits=(3.0, 10.0), axis=0),
    )
    print(
        "trimmed_tmean_all",
        stats.tmean(trimmed_matrix, limits=(3.0, 10.0), axis=None),
    )
    print(
        "trimmed_tvar_all",
        stats.tvar(trimmed_matrix, limits=(3.0, 10.0), axis=None),
    )

    proportion_matrix = np.array(
        [[1.0, 2.0, 3.0, 5.0], [10.0, 20.0, 30.0, 50.0]]
    )
    print("proportion_trim_mean", stats.trim_mean(np.arange(20.0), 0.1))
    print("proportion_trim_mean_axis0", stats.trim_mean(proportion_matrix, 0.25))
    print(
        "proportion_trim_mean_axis1",
        stats.trim_mean(proportion_matrix, 0.25, axis=1),
    )
    print(
        "proportion_trim_mean_all",
        stats.trim_mean(proportion_matrix, 0.25, axis=None),
    )
    print(
        "proportion_trimboth",
        np.sort(stats.trimboth(np.arange(10.0), 0.1)),
    )
    print(
        "proportion_trim1_left",
        np.sort(stats.trim1(np.arange(10.0), 0.2, tail="left")),
    )
    print(
        "proportion_trimboth_all",
        np.sort(stats.trimboth(proportion_matrix, 0.25, axis=None)),
    )
    proportion_nan = np.array([1.0, np.nan, 3.0, 5.0])
    print(
        "proportion_trim_mean_nan_omit",
        stats.trim_mean(proportion_nan, 0.25, nan_policy="omit"),
    )

    kstat_x = np.arange(1.0, 7.0)
    for order in range(1, 5):
        print(f"kstat_n{order}", stats.kstat(kstat_x, n=order))
    for order in (1, 2):
        print(f"kstatvar_n{order}", stats.kstatvar(kstat_x, n=order))
    kstat_matrix = np.arange(1.0, 13.0).reshape((3, 4), order="F")
    for axis in (0, 1):
        print(f"kstat_matrix_n2_axis{axis}", stats.kstat(kstat_matrix, n=2, axis=axis))
        print(f"kstat_matrix_n4_axis{axis}", stats.kstat(kstat_matrix, n=4, axis=axis))
        print(f"kstatvar_matrix_n1_axis{axis}", stats.kstatvar(kstat_matrix, n=1, axis=axis))
        print(f"kstatvar_matrix_n2_axis{axis}", stats.kstatvar(kstat_matrix, n=2, axis=axis))
    kstat_nan = np.array([1.0, 2.0, np.nan, 4.0, 5.0])
    for order in range(1, 5):
        print(
            f"kstat_nan_omit_n{order}",
            stats.kstat(kstat_nan, n=order, nan_policy="omit"),
        )
    for order in (1, 2):
        print(
            f"kstatvar_nan_omit_n{order}",
            stats.kstatvar(kstat_nan, n=order, nan_policy="omit"),
        )

    circular_angles = np.deg2rad(np.array([20.0, 30.0, 330.0]))
    print("circular_mean", stats.circmean(circular_angles))
    print("circular_var", stats.circvar(circular_angles))
    print("circular_std", stats.circstd(circular_angles))
    circular_degrees = np.array([-180.0, -170.0, 170.0])
    print(
        "circular_mean_degrees",
        stats.circmean(circular_degrees, high=180.0, low=-180.0),
    )
    print(
        "circular_std_degrees",
        stats.circstd(circular_degrees, high=180.0, low=-180.0),
    )
    print(
        "circular_std_degrees_normalized",
        stats.circstd(circular_degrees, high=180.0, low=-180.0, normalize=True),
    )
    circular_matrix = np.array(
        [[0.0, np.pi / 2.0, np.pi], [np.pi / 4.0, np.pi, 7.0 * np.pi / 4.0]]
    )
    for axis in (0, 1):
        print(f"circular_mean_axis{axis}", stats.circmean(circular_matrix, axis=axis))
        print(f"circular_var_axis{axis}", stats.circvar(circular_matrix, axis=axis))
        print(f"circular_std_axis{axis}", stats.circstd(circular_matrix, axis=axis))
    print("circular_mean_all", stats.circmean(circular_matrix, axis=None))
    circular_nan = np.array([0.0, np.nan, np.pi / 2.0])
    print("circular_mean_nan_omit", stats.circmean(circular_nan, nan_policy="omit"))
    print("circular_var_nan_omit", stats.circvar(circular_nan, nan_policy="omit"))
    print("circular_std_nan_omit", stats.circstd(circular_nan, nan_policy="omit"))

    rank_matrix = np.array([[3.0, 1.0, 1.0], [2.0, np.nan, 2.0]], order="F")
    rank_flat_f = rank_matrix.reshape(-1, order="F")
    for method in ("average", "min", "max", "dense", "ordinal"):
        print(
            f"rankdata_fortran_flat_{method}",
            stats.rankdata(rank_flat_f, method=method, nan_policy="omit"),
        )
    print("rankdata_axis0_omit", stats.rankdata(rank_matrix, axis=0, nan_policy="omit"))
    print("rankdata_axis1_omit", stats.rankdata(rank_matrix, axis=1, nan_policy="omit"))
    print("tiecorrect_example", stats.tiecorrect([1.0, 2.5, 2.5, 4.0]))
    print("tiecorrect_nan", stats.tiecorrect([1.0, np.nan, 1.0]))

    summary_mode = np.array(
        [[1.0, 2.0, 5.0, np.nan], [3.0, 2.0, 5.0, np.nan], [1.0, 4.0, 6.0, 7.0]],
        order="F",
    )
    print("summary_mode_axis0", stats.mode(summary_mode, axis=0, nan_policy="propagate"))
    print("summary_mode_axis0_omit", stats.mode(summary_mode, axis=0, nan_policy="omit"))
    print("summary_mode_all", stats.mode(summary_mode, axis=None, nan_policy="propagate"))
    summary_vector = np.arange(4.0)
    print("summary_describe", stats.describe(summary_vector))
    print("summary_describe_unbiased", stats.describe(summary_vector, bias=False))
    print("summary_describe_ddof4", stats.describe(summary_vector, ddof=4))
    print("summary_describe_ddof5", stats.describe(summary_vector, ddof=5))
    summary_nan = np.array([[1.0, np.nan, np.inf], [3.0, 5.0, 7.0]], order="F")
    print("summary_describe_nan_omit_axis0", stats.describe(summary_nan, axis=0, nan_policy="omit"))

    entropy_pk = np.array([[1.0, 3.0, 5.0], [2.0, 4.0, 6.0]], order="F")
    entropy_q = np.array([1.0, 2.0, 4.0])
    print("entropy_axis0", stats.entropy(entropy_pk, axis=0))
    print("entropy_axis1", stats.entropy(entropy_pk, axis=1))
    print("entropy_axis1_q", stats.entropy(entropy_pk, entropy_q, axis=1))
    print("entropy_axis1_q_base2", stats.entropy(entropy_pk, entropy_q, axis=1, base=2))
    print("entropy_all", stats.entropy(entropy_pk, axis=None))
    entropy_nan = entropy_pk.copy()
    entropy_nan[0, 1] = np.nan
    entropy_q_nan = np.broadcast_to(entropy_q, entropy_pk.shape).copy()
    entropy_q_nan[1, 2] = np.nan
    print("entropy_nan_omit", stats.entropy(entropy_nan, axis=1, nan_policy="omit"))
    print(
        "entropy_paired_nan_omit",
        stats.entropy(entropy_nan, entropy_q_nan, axis=1, nan_policy="omit"),
    )

    differential_x = np.array(
        [0.2, 0.5, 1.1, 1.7, 2.0, 3.4, 4.1, 5.6, 7.0, 8.5, 10.0]
    )
    for method in ("vasicek", "van es", "ebrahimi", "correa", "auto"):
        print(
            f"differential_entropy_{method.replace(' ', '_')}",
            stats.differential_entropy(differential_x, window_length=2, method=method),
        )
    print("differential_entropy_auto_default", stats.differential_entropy(differential_x))
    print(
        "differential_entropy_base2",
        stats.differential_entropy(
            differential_x, window_length=2, method="ebrahimi", base=2
        ),
    )
    differential_matrix = np.column_stack([differential_x, differential_x**1.1])
    print(
        "differential_entropy_axis0",
        stats.differential_entropy(
            differential_matrix, axis=0, window_length=2, method="ebrahimi"
        ),
    )
    differential_nan = differential_matrix.copy()
    differential_nan[1, 0] = np.nan
    differential_nan[5, 1] = np.nan
    print(
        "differential_entropy_nan_omit",
        stats.differential_entropy(
            differential_nan, axis=0, window_length=2, method="ebrahimi",
            nan_policy="omit",
        ),
    )

    frequency_a = np.array([1.0, 2.0, 3.0, 3.0, 4.0])
    for kind in ("rank", "strict", "weak", "mean"):
        print(
            f"frequency_percentileofscore_{kind}",
            stats.percentileofscore(frequency_a, 3.0, kind=kind),
        )
    print(
        "frequency_percentileofscore_vector",
        stats.percentileofscore(frequency_a, np.array([2.0, 3.0])),
    )
    frequency_matrix = np.arange(1.0, 13.0).reshape((3, 4), order="F")
    print(
        "frequency_scoreatpercentile_all",
        stats.scoreatpercentile(frequency_matrix, 50.0),
    )
    print(
        "frequency_scoreatpercentile_axis0",
        stats.scoreatpercentile(frequency_matrix, 50.0, axis=0),
    )
    print(
        "frequency_scoreatpercentile_axis1",
        stats.scoreatpercentile(frequency_matrix, 50.0, axis=1),
    )
    print(
        "frequency_scoreatpercentile_vector_axis0",
        stats.scoreatpercentile(frequency_matrix, [25.0, 50.0, 75.0], axis=0),
    )
    frequency_cum = stats.cumfreq(
        [1.0, 4.0, 2.0, 1.0, 3.0, 1.0],
        numbins=4,
        defaultreallimits=(1.5, 5.0),
    )
    print("frequency_cumfreq", frequency_cum)
    frequency_rel = stats.relfreq([2.0, 4.0, 1.0, 2.0, 3.0, 2.0], numbins=4)
    print("frequency_relfreq", frequency_rel)
    print(
        "frequency_relfreq_multidim",
        stats.relfreq(
            np.array([[1.0, 2.0], [3.0, 4.0]]),
            numbins=2,
            defaultreallimits=(1.0, 4.0),
        ),
    )

    named_matrix = np.array(
        [[1.0, 4.0, 9.0, 16.0], [2.0, 8.0, 18.0, 32.0], [3.0, 12.0, 27.0, 48.0]]
    )
    print("named_iqr_pair", stats.iqr(np.array([1.0, 3.0]), scale="normal"))
    print(
        "named_mad_pair",
        stats.median_abs_deviation(np.array([1.0, 3.0]), scale="normal"),
    )
    print(
        "named_iqr_custom_range",
        stats.iqr(np.array([1.0, 2.0, 4.0, 8.0, 16.0]), rng=(10, 90), scale="normal"),
    )
    print("named_matrix_iqr_axis0", stats.iqr(named_matrix, axis=0, scale="normal"))
    print(
        "named_matrix_mad_axis1",
        stats.median_abs_deviation(named_matrix, axis=1, scale="normal"),
    )
    print("named_matrix_iqr_all", stats.iqr(named_matrix, scale="normal"))
    print(
        "named_matrix_iqr_keepdims_axis1",
        stats.iqr(named_matrix, axis=1, scale="normal", keepdims=True),
    )

    nan_vector = np.array([1.0, np.nan, 4.0, 9.0])
    nan_weights = np.array([1.0, np.nan, 3.0, 4.0])
    print("nan_gmean_weighted", stats.gmean(nan_vector, weights=weights[:4], nan_policy="omit"))
    print("nan_gmean_weight_nan", stats.gmean(np.array([1.0, 2.0, 4.0, 9.0]), weights=nan_weights, nan_policy="omit"))
    print("nan_hmean_weighted", stats.hmean(nan_vector, weights=weights[:4], nan_policy="omit"))
    print("nan_pmean_weighted_2", stats.pmean(nan_vector, 2.0, weights=weights[:4], nan_policy="omit"))
    print("nan_sem", stats.sem(nan_vector, nan_policy="omit"))
    print("nan_variation", stats.variation(nan_vector, nan_policy="omit"))
    print("nan_iqr", stats.iqr(nan_vector, nan_policy="omit"))
    print("nan_mad", stats.median_abs_deviation(nan_vector, nan_policy="omit"))
    print("nan_moment_2", stats.moment(nan_vector, order=2, nan_policy="omit"))
    print("nan_skew", stats.skew(nan_vector, nan_policy="omit"))
    print("nan_kurtosis", stats.kurtosis(nan_vector, nan_policy="omit"))

    nan_matrix = np.array(
        [[1.0, np.nan, 4.0, 8.0], [3.0, 5.0, np.nan, 9.0], [np.nan, 6.0, 10.0, 14.0]]
    )
    nan_weights31 = np.array([[1.0], [np.nan], [4.0]])
    for axis in (0, 1):
        print(f"nan_matrix_gmean_axis{axis}", stats.gmean(nan_matrix, axis=axis, nan_policy="omit"))
        print(f"nan_matrix_hmean_axis{axis}", stats.hmean(nan_matrix, axis=axis, nan_policy="omit"))
        print(f"nan_matrix_pmean_axis{axis}", stats.pmean(nan_matrix, 2.0, axis=axis, nan_policy="omit"))
        print(f"nan_matrix_sem_axis{axis}", stats.sem(nan_matrix, axis=axis, nan_policy="omit"))
        print(f"nan_matrix_variation_axis{axis}", stats.variation(nan_matrix, axis=axis, nan_policy="omit"))
        print(f"nan_matrix_iqr_axis{axis}", stats.iqr(nan_matrix, axis=axis, nan_policy="omit"))
        print(f"nan_matrix_mad_axis{axis}", stats.median_abs_deviation(nan_matrix, axis=axis, nan_policy="omit"))
        print(f"nan_matrix_moment_axis{axis}", stats.moment(nan_matrix, order=2, axis=axis, nan_policy="omit"))
        print(f"nan_matrix_skew_axis{axis}", stats.skew(nan_matrix, axis=axis, nan_policy="omit"))
        print(f"nan_matrix_kurtosis_axis{axis}", stats.kurtosis(nan_matrix, axis=axis, nan_policy="omit"))
    print("nan_matrix_iqr_all", stats.iqr(nan_matrix, nan_policy="omit"))
    print(
        "nan_matrix_gmean_weights31_nan",
        stats.gmean(nan_matrix, axis=0, weights=nan_weights31, nan_policy="omit"),
    )
    print(
        "nan_matrix_gmean_keepdims",
        stats.gmean(nan_matrix, axis=0, keepdims=True, nan_policy="omit"),
    )
    print("nan_matrix_iqr_keepdims_all", stats.iqr(nan_matrix, keepdims=True, nan_policy="omit"))

    tensor_flat = np.array(
        [
            1.0, 2.0, 4.0, 9.0, 3.0, 5.0, 8.0, 15.0,
            2.0, 6.0, 11.0, 20.0, 4.0, 7.0, 13.0, 24.0,
            5.0, 10.0, 18.0, 31.0, 6.0, 12.0, 22.0, 40.0,
        ]
    )
    tensor = tensor_flat.reshape((4, 3, 2), order="F")
    for axis in (0, 1, 2):
        print(f"tensor_gmean_axis{axis}", stats.gmean(tensor, axis=axis))
        print(f"tensor_hmean_axis{axis}", stats.hmean(tensor, axis=axis))
        print(f"tensor_pmean_axis{axis}", stats.pmean(tensor, 2.0, axis=axis))
        print(f"tensor_sem_axis{axis}", stats.sem(tensor, axis=axis))
        print(f"tensor_variation_axis{axis}", stats.variation(tensor, axis=axis))
        print(f"tensor_iqr_axis{axis}", stats.iqr(tensor, axis=axis))
        print(f"tensor_mad_axis{axis}", stats.median_abs_deviation(tensor, axis=axis))
        print(f"tensor_moment_axis{axis}", stats.moment(tensor, order=2, axis=axis))
        print(f"tensor_skew_axis{axis}", stats.skew(tensor, axis=axis, bias=False))
        print(f"tensor_kurtosis_axis{axis}", stats.kurtosis(tensor, axis=axis, bias=False))
    print("tensor_iqr_all", stats.iqr(tensor))
    print("tensor_gmean_keepdims_axis1", stats.gmean(tensor, axis=1, keepdims=True))
    print("tensor_iqr_keepdims_all", stats.iqr(tensor, keepdims=True))

    tensor_weights1 = np.array([1.0, 3.0])
    tensor_weights2 = np.array([[1.0, 4.0], [2.0, 5.0], [3.0, 6.0]])
    tensor_weights3 = np.arange(1.0, 25.0).reshape((4, 3, 2), order="F")
    tensor_weights31 = np.array([1.0, 2.0, 3.0, 4.0]).reshape((4, 1, 1))
    print(
        "tensor_pmean_axis2_weights1",
        stats.pmean(tensor, 2.0, axis=2, weights=tensor_weights1),
    )
    print(
        "tensor_hmean_axis1_weights2",
        stats.hmean(tensor, axis=1, weights=tensor_weights2),
    )
    print(
        "tensor_gmean_axis0_weights3",
        stats.gmean(tensor, axis=0, weights=tensor_weights3),
    )
    print(
        "tensor_gmean_axis0_weights31",
        stats.gmean(tensor, axis=0, weights=tensor_weights31),
    )
    print(
        "tensor_gmean_keepdims_axis1_weights2",
        stats.gmean(tensor, axis=1, weights=tensor_weights2, keepdims=True),
    )
    singleton_tensor = tensor[:, :1, :]
    print(
        "tensor_expand_gmean_axis0_weights2",
        stats.gmean(singleton_tensor, axis=0, weights=tensor_weights2),
    )
    print(
        "tensor_expand_hmean_axis1_weights2",
        stats.hmean(singleton_tensor, axis=1, weights=tensor_weights2),
    )
    print(
        "tensor_expand_pmean_axis2_weights2",
        stats.pmean(singleton_tensor, 2.0, axis=2, weights=tensor_weights2),
    )

    nan_tensor = tensor.copy()
    nan_tensor[0, 0, 0] = np.nan
    nan_tensor[2, 1, 0] = np.nan
    nan_tensor[1, 2, 1] = np.nan
    nan_tensor[3, 0, 1] = np.nan
    nan_tensor_weights1 = np.array([1.0, np.nan])
    print(
        "nan_tensor_gmean_axis2_weights1",
        stats.gmean(
            nan_tensor, axis=2, weights=nan_tensor_weights1, nan_policy="omit"
        ),
    )
    print(
        "nan_tensor_pmean_axis2_weights1",
        stats.pmean(
            nan_tensor, 2.0, axis=2, weights=nan_tensor_weights1, nan_policy="omit"
        ),
    )
    print("nan_tensor_gmean_axis0", stats.gmean(nan_tensor, axis=0, nan_policy="omit"))
    print("nan_tensor_hmean_axis1", stats.hmean(nan_tensor, axis=1, nan_policy="omit"))
    print("nan_tensor_pmean_axis2", stats.pmean(nan_tensor, 2.0, axis=2, nan_policy="omit"))
    print("nan_tensor_sem_axis0", stats.sem(nan_tensor, axis=0, nan_policy="omit"))
    print("nan_tensor_variation_axis1", stats.variation(nan_tensor, axis=1, nan_policy="omit"))
    print("nan_tensor_iqr_axis2", stats.iqr(nan_tensor, axis=2, nan_policy="omit"))
    print(
        "nan_tensor_mad_axis0",
        stats.median_abs_deviation(nan_tensor, axis=0, nan_policy="omit"),
    )
    print(
        "nan_tensor_moment_axis1",
        stats.moment(nan_tensor, order=2, axis=1, nan_policy="omit"),
    )
    print(
        "nan_tensor_skew_axis0",
        stats.skew(nan_tensor, axis=0, bias=False, nan_policy="omit"),
    )
    print(
        "nan_tensor_kurtosis_axis1",
        stats.kurtosis(nan_tensor, axis=1, bias=False, nan_policy="omit"),
    )
    print("nan_tensor_iqr_all", stats.iqr(nan_tensor, nan_policy="omit"))
    print(
        "nan_tensor_sem_keepdims_axis2",
        stats.sem(nan_tensor, axis=2, keepdims=True, nan_policy="omit"),
    )


    rank4_n = np.arange(1.0, 49.0)
    rank4 = (0.25 * rank4_n**2 + 0.5 * rank4_n + 1.0).reshape(
        (4, 3, 2, 2), order="F"
    )
    print("rank4_gmean_axis0", stats.gmean(rank4, axis=0))
    print("rank4_hmean_axis1", stats.hmean(rank4, axis=1))
    print("rank4_pmean_axis2", stats.pmean(rank4, 2.0, axis=2))
    print("rank4_sem_axis3", stats.sem(rank4, axis=3))
    print("rank4_variation_axis1", stats.variation(rank4, axis=1))
    print("rank4_iqr_all", stats.iqr(rank4))
    print("rank4_iqr_axis3", stats.iqr(rank4, axis=3))
    print("rank4_mad_axis0", stats.median_abs_deviation(rank4, axis=0))
    print("rank4_moment_axis1", stats.moment(rank4, order=2, axis=1))
    print("rank4_skew_axis0", stats.skew(rank4, axis=0, bias=False))
    print("rank4_kurtosis_axis1", stats.kurtosis(rank4, axis=1, bias=False))
    print("rank4_gmean_keepdims_axis1", stats.gmean(rank4, axis=1, keepdims=True))
    print("rank4_iqr_keepdims_all", stats.iqr(rank4, keepdims=True))

    rank4_weight_data = np.arange(1.0, 7.0).reshape((2, 1, 1, 3), order="F")
    rank4_weights2 = np.arange(1.0, 13.0).reshape((4, 3), order="F")
    print(
        "rank4_weighted_gmean_joint",
        stats.gmean(rank4_weight_data, axis=0, weights=rank4_weights2),
    )
    print(
        "rank4_weighted_hmean_joint",
        stats.hmean(rank4_weight_data, axis=0, weights=rank4_weights2),
    )
    print(
        "rank4_weighted_pmean_joint",
        stats.pmean(rank4_weight_data, 2.0, axis=0, weights=rank4_weights2),
    )

    nan_rank4 = rank4.copy()
    nan_rank4[0, 0, 0, 0] = np.nan
    nan_rank4[2, 1, 0, 0] = np.nan
    nan_rank4[1, 2, 1, 1] = np.nan
    nan_rank4[3, 0, 1, 1] = np.nan
    nan_rank4[0, 2, 0, 1] = np.nan
    print(
        "nan_rank4_gmean_axis0",
        stats.gmean(nan_rank4, axis=0, nan_policy="omit"),
    )
    print(
        "nan_rank4_sem_axis3",
        stats.sem(nan_rank4, axis=3, nan_policy="omit"),
    )
    print("nan_rank4_iqr_all", stats.iqr(nan_rank4, nan_policy="omit"))
    print(
        "nan_rank4_skew_axis1",
        stats.skew(nan_rank4, axis=1, bias=False, nan_policy="omit"),
    )
    print(
        "nan_rank4_variation_axis0_propagate",
        stats.variation(nan_rank4, axis=0, nan_policy="propagate"),
    )

    rank15 = np.arange(1.0, 7.0).reshape((2, 3) + (1,) * 13, order="F")
    print("rank15_gmean_axis1", stats.gmean(rank15, axis=1))
    print("rank15_hmean_axis1", stats.hmean(rank15, axis=1))
    print("rank15_pmean_axis1", stats.pmean(rank15, 2.0, axis=1))
    print("rank15_gmean_keepdims_axis1", stats.gmean(rank15, axis=1, keepdims=True))
    rank15_weights1 = np.array([1.0, 2.0, 3.0, 4.0])
    rank15_weights15 = np.arange(1.0, 7.0).reshape(rank15.shape, order="F")
    print(
        "rank15_weighted_gmean_trailing_expand",
        stats.gmean(rank15, axis=1, weights=rank15_weights1),
    )
    print(
        "rank15_weighted_gmean_full",
        stats.gmean(rank15, axis=1, weights=rank15_weights15),
    )
    print(
        "rank15_weighted_hmean_full",
        stats.hmean(rank15, axis=1, weights=rank15_weights15),
    )
    print(
        "rank15_weighted_pmean_full",
        stats.pmean(rank15, 2.0, axis=1, weights=rank15_weights15),
    )
    nan_rank15 = rank15.copy()
    nan_rank15[(0, 1) + (0,) * 13] = np.nan
    print(
        "nan_rank15_gmean_axis1_omit",
        stats.gmean(nan_rank15, axis=1, nan_policy="omit"),
    )
    print(
        "nan_rank15_hmean_axis1_omit",
        stats.hmean(nan_rank15, axis=1, nan_policy="omit"),
    )
    print(
        "nan_rank15_pmean_axis1_omit",
        stats.pmean(nan_rank15, 2.0, axis=1, nan_policy="omit"),
    )
    nan_rank15_weights = rank15_weights15.copy()
    nan_rank15_weights[(0, 1) + (0,) * 13] = np.nan
    print(
        "nan_rank15_weighted_gmean_axis1_omit",
        stats.gmean(
            rank15, axis=1, weights=nan_rank15_weights, nan_policy="omit"
        ),
    )

    highrank_stats = np.arange(1.0, 49.0).reshape(
        (2, 3, 2, 2, 2) + (1,) * 10, order="F"
    )
    print("rank15_sem_axis4", stats.sem(highrank_stats, axis=4))
    print("rank15_variation_axis4", stats.variation(highrank_stats, axis=4))
    print("rank15_iqr_axis4", stats.iqr(highrank_stats, axis=4))
    print(
        "rank15_mad_axis4",
        stats.median_abs_deviation(highrank_stats, axis=4),
    )
    print("rank15_moment2_axis4", stats.moment(highrank_stats, moment=2, axis=4))
    print(
        "rank15_skew_axis4_unbiased",
        stats.skew(highrank_stats, axis=4, bias=False),
    )
    print(
        "rank15_kurtosis_axis4_unbiased",
        stats.kurtosis(highrank_stats, axis=4, fisher=True, bias=False),
    )
    print("rank15_iqr_all", stats.iqr(highrank_stats))
    print(
        "rank15_sem_keepdims_axis4",
        stats.sem(highrank_stats, axis=4, keepdims=True),
    )
    highrank_stats_sq = highrank_stats**2
    print(
        "rank15_moment3_axis1_nonlinear",
        stats.moment(highrank_stats_sq, moment=3, axis=1),
    )
    print(
        "rank15_skew_axis1_nonlinear_unbiased",
        stats.skew(highrank_stats_sq, axis=1, bias=False),
    )
    print(
        "rank15_sem_axis1_nonlinear",
        stats.sem(highrank_stats_sq, axis=1),
    )
    nan_highrank_stats = highrank_stats.copy()
    nan_highrank_stats[(0,) * 15] = np.nan
    for policy in ("propagate", "omit"):
        print(
            f"rank15_sem_axis4_{policy}",
            stats.sem(nan_highrank_stats, axis=4, nan_policy=policy),
        )
    print(
        "rank15_variation_axis4_omit",
        stats.variation(nan_highrank_stats, axis=4, nan_policy="omit"),
    )
    print(
        "rank15_iqr_axis4_omit",
        stats.iqr(nan_highrank_stats, axis=4, nan_policy="omit"),
    )
    print(
        "rank15_mad_axis4_omit",
        stats.median_abs_deviation(
            nan_highrank_stats, axis=4, nan_policy="omit"
        ),
    )
    print(
        "rank15_moment2_axis4_omit",
        stats.moment(nan_highrank_stats, moment=2, axis=4, nan_policy="omit"),
    )
    print(
        "rank15_skew_axis4_omit",
        stats.skew(nan_highrank_stats, axis=4, nan_policy="omit"),
    )
    print(
        "rank15_kurtosis_axis4_omit",
        stats.kurtosis(nan_highrank_stats, axis=4, nan_policy="omit"),
    )

    expectile_x = np.array([1.0, 4.0, 2.0, -1.0])
    for alpha in (0.0, 0.2, 0.5, 0.8, 1.0):
        print(f"expectile_alpha_{alpha}", stats.expectile(expectile_x, alpha))
    print(
        "expectile_weighted",
        stats.expectile(expectile_x, 0.8, weights=[1.0, 3.0, 1.0, 1.0]),
    )
    expectile_matrix = np.arange(1.0, 13.0).reshape((3, 4), order="F")
    print("expectile_matrix", stats.expectile(expectile_matrix, 0.3))
    print(
        "expectile_matrix_weights4",
        stats.expectile(expectile_matrix, 0.3, weights=[1.0, 2.0, 3.0, 4.0]),
    )
    print(
        "expectile_matrix_weights31",
        stats.expectile(
            expectile_matrix, 0.3, weights=np.array([[1.0], [2.0], [4.0]])
        ),
    )
    print(
        "expectile_negative_weights",
        stats.expectile([1.0, 2.0, 3.0], 0.2, weights=[1.0, -1.0, 1.0]),
    )
    clip_data = np.concatenate((np.linspace(9.5, 10.5, 31), np.linspace(0, 20, 5)))
    clipped, clip_lower, clip_upper = stats.sigmaclip(clip_data, 1.5, 1.5)
    print("sigmaclip_example", clipped)
    print("sigmaclip_example_limits", clip_lower, clip_upper)

    lmoment_x = np.array([1.0, 2.0, 2.0, 4.0, 8.0, 16.0])
    print("lmoment_default", stats.lmoment(lmoment_x))
    print("lmoment_raw", stats.lmoment(lmoment_x, standardize=False))
    print("lmoment_orders_1_3_5", stats.lmoment(lmoment_x, [1, 3, 5]))
    print("lmoment_sorted_true", stats.lmoment([3.0, 1.0, 2.0], sorted=True))
    lmoment_matrix = np.arange(1.0, 13.0).reshape((3, 4), order="F")
    print("lmoment_matrix_axis0", stats.lmoment(lmoment_matrix, axis=0))
    print("lmoment_matrix_axis1", stats.lmoment(lmoment_matrix, axis=1))
    print("lmoment_matrix_order2_axis1", stats.lmoment(lmoment_matrix, 2, axis=1))
    print("lmoment_matrix_orders_1_3_axis1", stats.lmoment(lmoment_matrix, [1, 3], axis=1))
    print("lmoment_matrix_all", stats.lmoment(lmoment_matrix, axis=None))
    lmoment_nan = np.array([1.0, np.nan, 3.0])
    print("lmoment_nan_omit", stats.lmoment(lmoment_nan, nan_policy="omit"))

    power_x = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    power_y = np.array([-2.0, -1.0, 0.0, 1.0, 2.0, 4.0])
    power_a = power_x.reshape((3, 2), order="F")
    power_b = power_y.reshape((3, 2), order="F")
    print("boxcox_transform_0p5", stats.boxcox(power_x, 0.5))
    print("boxcox_llf_0p5", stats.boxcox_llf(0.5, power_x))
    print("boxcox_llf_axis0", stats.boxcox_llf(0.5, power_a, axis=0))
    print("boxcox_normmax_pearsonr", stats.boxcox_normmax(power_x, method="pearsonr"))
    print("boxcox_normmax_mle", stats.boxcox_normmax(power_x, method="mle"))
    print("boxcox_normmax_all", stats.boxcox_normmax(power_x, method="all"))
    print("yeojohnson_transform_0p5", stats.yeojohnson(power_y, 0.5))
    print("yeojohnson_llf_0p5", stats.yeojohnson_llf(0.5, power_y))
    print("yeojohnson_llf_axis0", stats.yeojohnson_llf(0.5, power_b, axis=0))
    print("yeojohnson_normmax", stats.yeojohnson_normmax(power_y))


    prob_x = np.array([1.0, 2.0, 4.0, 8.0, 16.0, 32.0])
    (osm, osr), fit = stats.probplot(prob_x)
    print("probplot_norm_osm", osm)
    print("probplot_norm_osr", osr)
    print("probplot_norm_fit", fit)
    (_, _), fit_t = stats.probplot(prob_x, sparams=(5.0,), dist="t")
    print("probplot_t_fit", fit_t)
    print("ppcc_plot_tukey", stats.ppcc_plot(prob_x, -1.0, 1.0, N=5))
    print("ppcc_max_tukey_default", stats.ppcc_max(prob_x))
    print("ppcc_max_tukey_brack", stats.ppcc_max(prob_x, brack=(-1.0, 1.0)))
    ppcc_wx = np.array([0.5, 1.0, 1.7, 2.8, 4.4, 6.3, 8.9, 12.4])
    print(
        "ppcc_max_weibull",
        stats.ppcc_max(ppcc_wx, brack=(0.5, 4.0), dist="weibull_min"),
    )
    print("boxcox_normplot", stats.boxcox_normplot(prob_x, -1.0, 1.0, N=5))
    prob_y = np.array([-3.0, -1.0, 0.0, 1.0, 2.0, 5.0])
    print("yeojohnson_normplot", stats.yeojohnson_normplot(prob_y, -1.0, 2.0, N=5))

    distance_u = np.array([3.4, 3.9, 7.5, 7.8])
    distance_v = np.array([4.5, 1.4])
    distance_uw = np.array([1.4, 0.9, 3.1, 7.2])
    distance_vw = np.array([3.2, 3.5])
    print(
        "wasserstein_distance_weighted",
        stats.wasserstein_distance(distance_u, distance_v, distance_uw, distance_vw),
    )
    print(
        "energy_distance_weighted",
        stats.energy_distance(distance_u, distance_v, distance_uw, distance_vw),
    )
    print(
        "wasserstein_distance_repeated",
        stats.wasserstein_distance([0.0, 0.0, 1.0], [0.0, 2.0]),
    )
    print(
        "energy_distance_repeated",
        stats.energy_distance([0.0, 0.0, 1.0], [0.0, 2.0]),
    )

    distance_nd_u = np.array([[0.0, 2.0, 3.0], [1.0, 2.0, 5.0]])
    distance_nd_v = np.array([[3.0, 2.0, 3.0], [4.0, 2.0, 5.0]])
    print(
        "wasserstein_distance_nd_3d",
        stats.wasserstein_distance_nd(distance_nd_u, distance_nd_v),
    )
    distance_nd_uw = np.array([0.4, 5.2, 0.114])
    distance_nd_vw = np.array([0.8, 1.5])
    distance_nd_u2 = np.array([[0.0, 2.75], [2.0, 209.3], [0.0, 0.0]])
    distance_nd_v2 = np.array([[0.2, 0.322], [4.5, 25.1808]])
    print(
        "wasserstein_distance_nd_weighted",
        stats.wasserstein_distance_nd(
            distance_nd_u2, distance_nd_v2, distance_nd_uw, distance_nd_vw
        ),
    )

    ecdf_sample = np.array([6.23, 5.58, 7.06, 6.42, 5.20])
    ecdf_result = stats.ecdf(ecdf_sample)
    print("ecdf_uncensored_quantiles", ecdf_result.cdf.quantiles)
    print("ecdf_uncensored_cdf", ecdf_result.cdf.probabilities)
    ecdf_ci = ecdf_result.cdf.confidence_interval(0.9, method="linear")
    print("ecdf_uncensored_linear_low", ecdf_ci.low.probabilities)
    ecdf_ci = ecdf_result.cdf.confidence_interval(0.9, method="log-log")
    print("ecdf_uncensored_loglog_low", ecdf_ci.low.probabilities)

    broken = np.array([77.0, 47.0, 81.0, 56.0, 80.0])
    unbroken = np.array([62.0, 60.0, 43.0, 71.0, 37.0])
    censored = stats.CensoredData(uncensored=broken, right=unbroken)
    ecdf_censored = stats.ecdf(censored)
    print("ecdf_right_censored_quantiles", ecdf_censored.sf.quantiles)
    print("ecdf_right_censored_sf", ecdf_censored.sf.probabilities)

    quantile_x = np.array([0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
    for alternative in ("two-sided", "less", "greater"):
        result = stats.quantile_test(quantile_x, q=1.0, p=0.5, alternative=alternative)
        print(f"quantile_test_{alternative}", result)
        print(f"quantile_test_{alternative}_ci90", result.confidence_interval(0.9))

    order_dist = stats.order_statistic(stats.Normal(), r=2, n=5)
    print("order_statistic_normal_pdf", order_dist.pdf([-1.0, 0.0, 1.0]))
    print("order_statistic_normal_cdf", order_dist.cdf([-1.0, 0.0, 1.0]))
    print("order_statistic_normal_ppf", order_dist.icdf([0.1, 0.5, 0.9]))

    logrank_x = stats.CensoredData(
        uncensored=[6, 13, 21, 30, 37, 38, 49, 50, 63, 79, 86, 98, 202, 219],
        right=[31, 47, 80, 82, 82, 149],
    )
    logrank_y = stats.CensoredData(
        uncensored=[10, 10, 12, 13, 14, 15, 16, 17, 18, 20, 24, 24, 25, 28,
                    30, 33, 35, 37, 40, 40, 46, 48, 76, 81, 82, 91, 112, 181],
        right=[34, 40, 70],
    )
    for alternative in ("two-sided", "less", "greater"):
        print(f"logrank_{alternative}", stats.logrank(logrank_x, logrank_y, alternative=alternative))

    for alternative in ("two-sided", "less", "greater"):
        print(
            f"ttest_1samp_{alternative}",
            stats.ttest_1samp(x, 3.0, alternative=alternative),
        )
        print(
            f"pearsonr_{alternative}",
            stats.pearsonr(x, y, alternative=alternative),
        )
        print(
            f"spearmanr_{alternative}",
            stats.spearmanr(x, y, alternative=alternative),
        )
        print(
            f"mannwhitneyu_{alternative}",
            stats.mannwhitneyu(
                x,
                y,
                alternative=alternative,
                method="asymptotic",
                use_continuity=True,
            ),
        )
    print("ttest_ind_equal", stats.ttest_ind(x, y, equal_var=True))
    print("ttest_ind_welch", stats.ttest_ind(x, y, equal_var=False))
    print(
        "ttest_ind_from_stats_equal",
        stats.ttest_ind_from_stats(15.0, np.sqrt(87.5), 13, 12.0, np.sqrt(39.0), 11),
    )
    print(
        "ttest_ind_from_stats_welch",
        stats.ttest_ind_from_stats(
            15.0, np.sqrt(87.5), 13, 12.0, np.sqrt(39.0), 11, equal_var=False
        ),
    )
    for alternative in ("less", "greater"):
        print(
            f"ttest_ind_from_stats_{alternative}",
            stats.ttest_ind_from_stats(
                15.0,
                np.sqrt(87.5),
                13,
                12.0,
                np.sqrt(39.0),
                11,
                alternative=alternative,
            ),
        )
    print(
        "ttest_ind_from_stats_negative_std",
        stats.ttest_ind_from_stats(1.0, -2.0, 5, 2.0, -3.0, 7),
    )
    print(
        "ttest_ind_from_stats_single_equal",
        stats.ttest_ind_from_stats(1.0, 2.0, 1, 2.0, 3.0, 5),
    )
    print(
        "ttest_ind_from_stats_single_welch_zero_std",
        stats.ttest_ind_from_stats(1.0, 0.0, 1, 2.0, 3.0, 5, equal_var=False),
    )
    print(
        "ttest_ind_from_stats_vectorized",
        stats.ttest_ind_from_stats(
            np.array([15.0, 0.2]),
            np.array([np.sqrt(87.5), np.sqrt(0.161073)]),
            np.array([13, 150]),
            np.array([12.0, 0.225]),
            np.array([np.sqrt(39.0), np.sqrt(0.175251)]),
            np.array([11, 200]),
        ),
    )
    print("ttest_rel", stats.ttest_rel(x, y))
    print(
        "mannwhitneyu_no_continuity",
        stats.mannwhitneyu(
            x,
            y,
            method="asymptotic",
            use_continuity=False,
        ),
    )

    obrien_x1 = np.array([0.0, 2.0, 4.0])
    obrien_x2 = np.array([0.0, 3.0, 6.0, 9.0])
    obrien_doc_x = np.array([10.0, 11.0, 13.0, 9.0, 7.0, 12.0, 12.0, 9.0, 10.0])
    obrien_doc_y = np.array([13.0, 21.0, 5.0, 10.0, 8.0, 14.0, 10.0, 12.0, 7.0, 15.0])
    print("obrientransform_n3", stats.obrientransform(obrien_x1))
    print("obrientransform_n4", stats.obrientransform(obrien_x2))
    obrien_tx, obrien_ty = stats.obrientransform(obrien_doc_x, obrien_doc_y)
    print("obrientransform_doc_x", obrien_tx)
    print("obrientransform_doc_y", obrien_ty)
    print("obrientransform_doc_anova", stats.f_oneway(obrien_tx, obrien_ty))
    print("obrientransform_constant", stats.obrientransform([2.0, 2.0, 2.0, 2.0]))


if __name__ == "__main__":
    main()
