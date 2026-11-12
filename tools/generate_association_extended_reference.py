#!/usr/bin/env python3
"""Generate SciPy 1.17.0 references for association and trend procedures."""
from pathlib import Path
import numpy as np
from scipy import stats

OUT = Path(__file__).resolve().parents[1] / "test" / "association_extended_reference.f90"


def fmt(x):
    return f"{float(x):.17e}_dp"


def arr(name, vals):
    vals = np.asarray(vals).ravel()
    chunks = []
    for i in range(0, len(vals), 2):
        chunks.append(", ".join(fmt(v) for v in vals[i:i+2]))
    if len(chunks) == 1:
        return f"    real(dp), parameter, public :: {name}({len(vals)}) = [{chunks[0]}]\n"
    body = ", &\n        ".join(chunks)
    return f"    real(dp), parameter, public :: {name}({len(vals)}) = [{body}]\n"

# Kendall tau: tied, exact untied, asymptotic untied, and tau-c.
kx_tie = np.array([12., 2., 1., 12., 2.])
ky_tie = np.array([1., 4., 7., 1., 0.])
kendall_tie = []
for alt in ("two-sided", "less", "greater"):
    r = stats.kendalltau(kx_tie, ky_tie, alternative=alt)
    kendall_tie.extend([r.statistic, r.pvalue])
r = stats.kendalltau(kx_tie, ky_tie, variant="c")
kendall_tie.extend([r.statistic, r.pvalue])

kx = np.arange(1., 6.)
ky = np.array([5., 1., 4., 2., 3.])
kendall_exact = []
for alt in ("two-sided", "less", "greater"):
    r = stats.kendalltau(kx, ky, method="exact", alternative=alt)
    kendall_exact.extend([r.statistic, r.pvalue])
for alt in ("two-sided", "less", "greater"):
    r = stats.kendalltau(kx, ky, method="asymptotic", alternative=alt)
    kendall_exact.extend([r.statistic, r.pvalue])

kx6 = np.arange(1., 7.)
ky6 = np.array([1., 3., 2., 6., 4., 5.])
kendall_6 = []
for alt in ("two-sided", "less", "greater"):
    r = stats.kendalltau(kx6, ky6, method="exact", alternative=alt)
    kendall_6.extend([r.statistic, r.pvalue])

pbx = np.array([0., 0., 0., 1., 1., 1., 1.])
pby = np.arange(7.)
r = stats.pointbiserialr(pbx, pby)
point_ref = [r.statistic, r.pvalue]

# Spearman rho: SciPy 1.17 replacement API, asymptotic and exact permutation inference.
srx = np.array([1., 2., 3., 4., 5.])
sry = np.array([5., 6., 7., 8., 7.])
spearmanrho_ref = []
for alt in ("two-sided", "less", "greater"):
    r = stats.spearmanrho(srx, sry, alternative=alt)
    spearmanrho_ref.extend([r.statistic, r.pvalue])
method = stats.PermutationMethod(n_resamples=9999, rng=np.random.default_rng(24680))
for alt in ("two-sided", "less", "greater"):
    r = stats.spearmanrho(srx, sry, alternative=alt, method=method)
    spearmanrho_ref.extend([r.statistic, r.pvalue])

# Weighted Kendall tau: default/multiplicative weights, rank modes, custom ranks, and weighers.
wx = np.array([12., 2., 1., 12., 2.])
wy = np.array([1., 4., 7., 1., 0.])
weightedtau_ref = [
    stats.weightedtau(wx, wy).statistic,
    stats.weightedtau(wx, wy, additive=False).statistic,
    stats.weightedtau(wx, wy, weigher=lambda r: 1.0).statistic,
    stats.weightedtau(wx, wy, rank=None).statistic,
    stats.weightedtau(wx, wy, rank=False).statistic,
    stats.weightedtau(wx, wy, rank=np.array([0, 3, 1, 4, 2])).statistic,
    stats.weightedtau(wx, wy, rank=np.array([0., 30., 10., 40., 20.])).statistic,
    stats.weightedtau(wx, wy, weigher=lambda r: 1.0 / (r + 1.0)**2).statistic,
    stats.weightedtau(wx, wy, rank=np.array([0, 3, 1, 4, 2]),
                      weigher=lambda r: 1.0 / (r + 1.0)**2, additive=False).statistic,
    stats.weightedtau(wx, np.array([1., 4., 7., 1., np.nan])).statistic,
    stats.weightedtau(wx, wy, rank=np.array([5, 1, 5, 2, 1])).statistic,
]

# Somers D: canonical contingency table, transpose asymmetry, and tied samples.
somers_table = np.array([[27, 25, 14, 7, 0],
                         [7, 14, 18, 35, 12],
                         [1, 3, 2, 7, 17]])
somers_table_ref = []
for alt in ("two-sided", "less", "greater"):
    r = stats.somersd(somers_table, alternative=alt)
    somers_table_ref.extend([r.statistic, r.pvalue])
r = stats.somersd(somers_table.T)
somers_table_ref.extend([r.statistic, r.pvalue])

sx = np.array([1., 1., 2., 2., 3., 3., 3.])
sy = np.array([1., 2., 2., 3., 2., 4., 4.])
somers_sample_ref = []
for alt in ("two-sided", "less", "greater"):
    r = stats.somersd(sx, sy, alternative=alt)
    somers_sample_ref.extend([r.statistic, r.pvalue])

# Chatterjee xi: tied-y, continuous-y assumption, negative statistic, and exact permutation.
cx = np.arange(-3., 4.)
cy = np.array([9., 4., 1., 0., 1., 4., 9.])
chatterjee_ref = []
for continuous in (False, True):
    r = stats.chatterjeexi(cx, cy, y_continuous=continuous)
    chatterjee_ref.extend([r.statistic, r.pvalue])

cnx = np.arange(10.)
cny = (cnx % 2 == 0).astype(float)
for continuous in (False, True):
    r = stats.chatterjeexi(cnx, cny, y_continuous=continuous)
    chatterjee_ref.extend([r.statistic, r.pvalue])

cpx = np.arange(5.)
cpy = np.array([4., 1., 3., 0., 2.])
method = stats.PermutationMethod(n_resamples=9999, rng=np.random.default_rng(12345))
r = stats.chatterjeexi(cpx, cpy, method=method)
chatterjee_ref.extend([r.statistic, r.pvalue])

lx = np.array([0.2, 1.1, 2.7, 3.0, 5.2, 8.1])
ly = np.array([1.0, 1.7, 4.2, 4.0, 8.9, 12.3])
lin_ref = []
for alt in ("two-sided", "less", "greater"):
    r = stats.linregress(lx, ly, alternative=alt)
    lin_ref.extend([r.slope, r.intercept, r.rvalue, r.pvalue, r.stderr, r.intercept_stderr])

rx = np.arange(7.)
ry = np.array([1., 2.2, 2.7, 10., 5.1, 6.2, 7.1])
theil_ref = []
for method in ("separate", "joint"):
    r = stats.theilslopes(ry, rx, alpha=.9, method=method)
    theil_ref.extend([r.slope, r.intercept, r.low_slope, r.high_slope])
siegel_ref = []
for method in ("hierarchical", "separate"):
    r = stats.siegelslopes(ry, rx, method=method)
    siegel_ref.extend([r.slope, r.intercept])

rx2 = np.array([0., 0., 1., 2., 3., 4.])
ry2 = np.array([1., 1.3, 2.1, 3.2, 4.5, 5.])
r = stats.theilslopes(ry2, rx2, alpha=.95, method="separate")
theil_repeat_ref = [r.slope, r.intercept, r.low_slope, r.high_slope]
r = stats.siegelslopes(ry2, rx2, method="hierarchical")
siegel_repeat_ref = [r.slope, r.intercept]

bx = np.array([1,2,1,1,1,1,1,1,1,1,2,4,1,1], dtype=float)
by = np.array([3,3,4,3,1,2,3,1,1,5,4], dtype=float)
bm_ref = []
for dist in ("t", "normal"):
    for alt in ("two-sided", "less", "greater"):
        r = stats.brunnermunzel(bx, by, alternative=alt, distribution=dist)
        bm_ref.extend([r.statistic, r.pvalue])

page_data = np.array([
    [3., 4., 3.],
    [2., 2., 4.],
    [3., 3., 5.],
    [1., 3., 2.],
    [2., 3., 2.],
    [2., 4., 5.],
    [1., 2., 4.],
    [3., 4., 4.],
    [2., 4., 5.],
    [1., 3., 4.],
])
page_ref = []
for method in ("exact", "asymptotic"):
    r = stats.page_trend_test(page_data, method=method)
    page_ref.extend([r.statistic, r.pvalue])
# Pre-ranked/reordered prediction gives the same physical trend.
ranks = np.vstack([stats.rankdata(row) for row in page_data])
r = stats.page_trend_test(ranks[:, [2, 0, 1]], ranked=True,
                          predicted_ranks=[3, 1, 2], method="exact")
page_ref.extend([r.statistic, r.pvalue])

page4 = np.array([[4.,1.,3.,2.], [2.,1.,4.,3.], [1.,2.,3.,4.]])
page4_ref = []
for method in ("exact", "asymptotic"):
    r = stats.page_trend_test(page4, method=method)
    page4_ref.extend([r.statistic, r.pvalue])

lines = [
    "! Generated by tools/generate_association_extended_reference.py with SciPy 1.17.0.\n",
    "module association_extended_reference\n",
    "    use scifort_kinds, only : dp\n",
    "    implicit none\n",
    "    private\n",
    arr("kendall_tie_ref", kendall_tie),
    arr("kendall_exact_ref", kendall_exact),
    arr("kendall_6_ref", kendall_6),
    arr("pointbiserial_ref", point_ref),
    arr("spearmanrho_ref", spearmanrho_ref),
    arr("weightedtau_ref", weightedtau_ref),
    arr("somers_table_ref", somers_table_ref),
    arr("somers_sample_ref", somers_sample_ref),
    arr("chatterjee_ref", chatterjee_ref),
    arr("linregress_ref", lin_ref),
    arr("theil_ref", theil_ref),
    arr("siegel_ref", siegel_ref),
    arr("theil_repeat_ref", theil_repeat_ref),
    arr("siegel_repeat_ref", siegel_repeat_ref),
    arr("brunnermunzel_ref", bm_ref),
    arr("page_ref", page_ref),
    arr("page4_ref", page4_ref),
    "end module association_extended_reference\n",
]
OUT.write_text("".join(lines), encoding="ascii")
