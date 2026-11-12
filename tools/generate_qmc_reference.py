#!/usr/bin/env python3
"""Generate deterministic SciPy 1.17 reference values for QMC tests."""

from pathlib import Path
import numpy as np
import scipy
from scipy.stats import qmc, sobol_indices
from scipy.spatial import Voronoi
from scipy.stats._qmc import van_der_corput

assert scipy.__version__ == "1.17.0", scipy.__version__

space_raw = np.array([[1, 3], [2, 6], [3, 2], [4, 5], [5, 1], [6, 4]], dtype=float)
space = qmc.scale(space_raw, [0.5, 0.5], [6.5, 6.5], reverse=True)
methods = ["CD", "WD", "MD", "L2-star"]
disc = [qmc.discrepancy(space, method=m) for m in methods]
disc_iter = qmc.discrepancy(space[:-1], iterative=True, method="CD")
disc_update = qmc.update_discrepancy(space[-1], space[:-1], disc_iter)
geo_min = qmc.geometric_discrepancy(space, method="mindist")
geo_mst = qmc.geometric_discrepancy(space, method="mst")
geo_city = qmc.geometric_discrepancy(space, method="mindist", metric="cityblock")

# Geometry reference for SciFort's dependency-free Lloyd fallback.  SciPy's
# public QMC optimizer uses Qhull with QJ joggle; SciFort deliberately omits
# that joggle, so validate one mathematical Voronoi-relaxation step against
# the same Qhull geometry without QJ.
lloyd_input = np.array([
    [0.4338946866958654, 0.39327530366669794],
    [0.3424063572208076, 0.3637226514133418],
    [0.22580773595182554, 0.4099602705749197],
    [0.6536383363505450, 0.6980264897828301],
    [0.7504540230597282, 0.5743578379353608],
    [0.3009291948166416, 0.8688066670720066],
    [0.8197004358425417, 0.2606038141121997],
    [0.6613257075833527, 0.21651468999055032],
    [0.2556896234297943, 0.4860364686385993],
    [0.10496270977626461, 0.24000260422988468],
], dtype=float)
lloyd_reference = lloyd_input.copy()
vor = Voronoi(lloyd_input, qhull_options="Qbb Qc Qz")
for i, region_index in enumerate(vor.point_region):
    region = [vertex for vertex in vor.regions[region_index] if vertex != -1]
    if not region:
        continue
    centroid = vor.vertices[region].mean(axis=0)
    candidate = lloyd_input[i] + 1.9 * (centroid - lloyd_input[i])
    if np.all((candidate >= 0.0) & (candidate <= 1.0)):
        lloyd_reference[i] = candidate


halton = qmc.Halton(d=3, scramble=False)
halton_first = halton.random(8)
halton.fast_forward(5)
halton_later = halton.random(4)
vdc = van_der_corput(8, base=3, start_index=2, scramble=False)

scale_in = np.array([[0.5, 0.75], [0.5, 0.5], [0.75, 0.25]])
scale_out = qmc.scale(scale_in, [-2, 0], [6, 5])

sobol = qmc.Sobol(d=3, scramble=False, bits=30)
sobol_first = sobol.random_base2(3)
sobol.reset()
sobol.fast_forward(4)
sobol_later = sobol.random_base2(2)
sobol_hi = qmc.Sobol(d=21201, scramble=False, bits=64).random(16)[15, [0, 99, 999, 21200]]
sobol_deep_engine = qmc.Sobol(d=7, scramble=False, bits=20)
sobol_deep_engine.fast_forward(12345)
sobol_deep = sobol_deep_engine.random(6)
sobol_bits4 = qmc.Sobol(d=2, scramble=False, bits=4)
sobol_bits4.random_base2(3)
sobol_second8 = sobol_bits4.random_base2(3)

multinomial_qmc = qmc.MultinomialQMC(
    pvals=[0.2, 0.4, 0.4], n_trials=8,
    engine=qmc.Sobol(d=1, scramble=False, bits=30),
)
multinomial_sample = multinomial_qmc.random(4).astype(int)

mvn_qmc = qmc.MultivariateNormalQMC(
    mean=[1.0, -2.0], cov=[[2.0, 0.5], [0.5, 1.0]],
    engine=qmc.Sobol(d=2, scramble=False, bits=30),
)
mvn_sample = mvn_qmc.random(8)
mvn_singular_qmc = qmc.MultivariateNormalQMC(
    mean=[0.5, -1.0], cov=[[1.0, 1.0], [1.0, 1.0]],
    engine=qmc.Sobol(d=2, scramble=False, bits=30),
)
mvn_singular_sample = mvn_singular_qmc.random(4)


# scipy.stats.sobol_indices references from precomputed model evaluations.
sobol_fa = np.array([[1., 2., 3., 4., 5., 6., 7., 8.]])
sobol_fb = np.array([[8., 7., 6., 5., 4., 3., 2., 1.]])
sobol_fab = np.empty((3, 1, 8))
sobol_fab[0, 0] = [2., 3., 4., 5., 6., 7., 8., 1.]
sobol_fab[1, 0] = [1., 3., 2., 4., 6., 5., 8., 7.]
sobol_fab[2, 0] = [8., 2., 3., 4., 5., 6., 7., 1.]
sobol_sensitivity = sobol_indices(
    func={'f_A': sobol_fa, 'f_B': sobol_fb, 'f_AB': sobol_fab}, n=8
)

sobol_fa_multi = np.array([[1., 2., 4., 8.], [0., 1., 0., 1.]])
sobol_fb_multi = np.array([[3., 5., 7., 9.], [1., 1., 0., 0.]])
sobol_fab_multi = np.empty((2, 2, 4))
sobol_fab_multi[0] = np.array([[2., 3., 5., 7.], [1., 1., 0., 1.]])
sobol_fab_multi[1] = np.array([[1., 4., 6., 8.], [0., 0., 1., 1.]])
sobol_sensitivity_multi = sobol_indices(
    func={'f_A': sobol_fa_multi, 'f_B': sobol_fb_multi, 'f_AB': sobol_fab_multi}, n=4
)


def fr(x):
    text = f"{float(x):.17g}"
    if "e" not in text and "." not in text:
        text += ".0"
    return text + "_dp"


def write_vec(name, arr):
    values = [fr(v) for v in np.asarray(arr).ravel()]
    lines = [f"    real(dp), parameter :: {name}({len(values)}) = [ &"]
    for i in range(0, len(values), 3):
        chunk = ", ".join(values[i:i+3])
        suffix = ", &" if i + 3 < len(values) else " &"
        lines.append(f"        {chunk}{suffix}")
    lines.append("    ]")
    return "\n".join(lines) + "\n"


def write_mat(name, arr):
    arr = np.asarray(arr)
    # Fortran reshape consumes column-major values; flatten in Fortran order.
    vals = [fr(v) for v in arr.flatten(order="F")]
    lines = [f"    real(dp), parameter :: {name}({arr.shape[0]},{arr.shape[1]}) = reshape([ &"]
    for i in range(0, len(vals), 4):
        chunk = ", ".join(vals[i:i+4])
        suffix = ", &" if i + 4 < len(vals) else " &"
        lines.append(f"        {chunk}{suffix}")
    lines.append(f"    ], [{arr.shape[0]},{arr.shape[1]}])")
    return "\n".join(lines) + "\n"


def write_int_mat(name, arr):
    arr = np.asarray(arr, dtype=int)
    vals = [str(int(v)) for v in arr.flatten(order="F")]
    lines = [f"    integer, parameter :: {name}({arr.shape[0]},{arr.shape[1]}) = reshape([ &"]
    for i in range(0, len(vals), 8):
        chunk = ", ".join(vals[i:i+8])
        suffix = ", &" if i + 8 < len(vals) else " &"
        lines.append(f"        {chunk}{suffix}")
    lines.append(f"    ], [{arr.shape[0]},{arr.shape[1]}])")
    return "\n".join(lines) + "\n"

out = [
    "! Generated by tools/generate_qmc_reference.py with SciPy 1.17.0.\n",
    "module qmc_reference\n",
    "    use scifort_kinds, only : dp\n",
    "    implicit none\n",
    "    private\n",
    "    public :: qmc_ref_space, qmc_ref_discrepancy, qmc_ref_iterative, qmc_ref_update\n",
    "    public :: qmc_ref_geo, qmc_ref_halton_first, qmc_ref_halton_later, qmc_ref_vdc\n",
    "    public :: qmc_ref_scale, qmc_ref_sobol_first, qmc_ref_sobol_later, qmc_ref_sobol_hi\n",
    "    public :: qmc_ref_sobol_deep, qmc_ref_sobol_second8\n",
    "    public :: qmc_ref_multinomial, qmc_ref_mvn, qmc_ref_mvn_singular\n",
    "    public :: qmc_ref_lloyd_input, qmc_ref_lloyd_once\n",
    "    public :: qmc_ref_sensitivity_first, qmc_ref_sensitivity_total\n",
    "    public :: qmc_ref_sensitivity_multi_first, qmc_ref_sensitivity_multi_total\n\n",
]
out.append(write_mat("qmc_ref_space", space))
out.append(write_vec("qmc_ref_discrepancy", disc))
out.append(f"    real(dp), parameter :: qmc_ref_iterative = {fr(disc_iter)}\n")
out.append(f"    real(dp), parameter :: qmc_ref_update = {fr(disc_update)}\n")
out.append(write_vec("qmc_ref_geo", [geo_min, geo_mst, geo_city]))
out.append(write_mat("qmc_ref_lloyd_input", lloyd_input))
out.append(write_mat("qmc_ref_lloyd_once", lloyd_reference))
out.append(write_mat("qmc_ref_halton_first", halton_first))
out.append(write_mat("qmc_ref_halton_later", halton_later))
out.append(write_vec("qmc_ref_vdc", vdc))
out.append(write_mat("qmc_ref_scale", scale_out))
out.append(write_mat("qmc_ref_sobol_first", sobol_first))
out.append(write_mat("qmc_ref_sobol_later", sobol_later))
out.append(write_vec("qmc_ref_sobol_hi", sobol_hi))
out.append(write_mat("qmc_ref_sobol_deep", sobol_deep))
out.append(write_mat("qmc_ref_sobol_second8", sobol_second8))
out.append(write_int_mat("qmc_ref_multinomial", multinomial_sample))
out.append(write_mat("qmc_ref_mvn", mvn_sample))
out.append(write_mat("qmc_ref_mvn_singular", mvn_singular_sample))
out.append(write_vec("qmc_ref_sensitivity_first", sobol_sensitivity.first_order))
out.append(write_vec("qmc_ref_sensitivity_total", sobol_sensitivity.total_order))
out.append(write_mat("qmc_ref_sensitivity_multi_first", sobol_sensitivity_multi.first_order))
out.append(write_mat("qmc_ref_sensitivity_multi_total", sobol_sensitivity_multi.total_order))
out.append("end module qmc_reference\n")
Path("test/qmc_reference.f90").write_text("".join(out), encoding="ascii")
