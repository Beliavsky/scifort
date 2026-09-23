# SPDX-License-Identifier: MIT
# Copyright (c) 2026 SciFort contributors
"""Print reference values used by descriptive and hypothesis tests.

This is a developer validation helper, not a runtime dependency. It records
NumPy/SciPy versions and reproduces the reference values embedded in
``test/test_descriptive.f90`` and ``test/test_hypothesis.f90``.
"""

import numpy as np
import scipy
from scipy import stats


def main():
    x = np.array([1.0, 2.0, 2.0, 4.0, 8.0, 16.0])
    y = np.array([3.0, -1.0, 5.0, 4.0, 9.0, 12.0])
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


if __name__ == "__main__":
    main()
