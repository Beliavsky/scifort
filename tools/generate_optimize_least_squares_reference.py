#!/usr/bin/env python3
"""Generate SciPy 1.17 references for nonlinear least squares tests."""
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import curve_fit, least_squares, leastsq
from scipy.sparse.linalg import LinearOperator

OUT = Path(__file__).resolve().parents[1] / "test" / "optimize_least_squares_reference.f90"


def residual(x):
    return np.array([10.0 * (x[1] - x[0] ** 2), 1.0 - x[0]])


def residual_jac(x):
    return np.array([[-20.0 * x[0], 10.0], [-1.0, 0.0]])


def residual_jac_operator(x):
    jac = residual_jac(x)
    return LinearOperator(
        jac.shape, matvec=lambda v: jac @ v, rmatvec=lambda w: jac.T @ w
    )



SPARSE_PATTERN = np.array([
    [1, 1, 0, 0],
    [1, 1, 0, 0],
    [0, 0, 1, 1],
    [0, 0, 1, 1],
], dtype=bool)
SPARSE_X0 = np.array([0.5, 1.5, 2.5, 3.5])
sparse_calls = 0


def sparse_residual(x):
    global sparse_calls
    sparse_calls += 1
    return np.array([
        x[0] + x[1] - 3.0,
        x[0] - x[1] + 1.0,
        x[2] + x[3] - 7.0,
        x[2] - x[3] + 1.0,
    ])


def robust_residual(x):
    data = np.array([1.0, 1.2, 0.8, 10.0])
    return x[0] - data


def robust_jac(x):
    return np.ones((4, 1))



def stop_after_two(intermediate_result):
    if intermediate_result.nit >= 2:
        raise StopIteration

def model(x, a, b):
    return a * x + b


def model_jac(x, a, b):
    return np.column_stack([x, np.ones_like(x)])


def freal(x):
    return f"{float(x):.17e}_dp"


def vec(name, values):
    vals = np.asarray(values).ravel()
    lines = [f"    real(dp), parameter :: {name}({vals.size}) = [ &"]
    for i, value in enumerate(vals):
        end = "]" if i == vals.size - 1 else ", &"
        lines.append(f"        {freal(value)}{end}")
    return lines


def mat(name, values):
    if hasattr(values, "toarray"):
        values = values.toarray()
    a = np.asarray(values)
    flat = a.flatten(order="F")
    lines = [f"    real(dp), parameter :: {name}({a.shape[0]},{a.shape[1]}) = reshape([ &"]
    for i, value in enumerate(flat):
        end = "], [" if i == flat.size - 1 else ", &"
        lines.append(f"        {freal(value)}{end}")
    lines[-1] += f"{a.shape[0]},{a.shape[1]}])"
    return lines


def scalar(name, value):
    return [f"    real(dp), parameter :: {name} = {freal(value)}"]


def integer(name, value):
    return [f"    integer, parameter :: {name} = {int(value)}"]


def ivec(name, values):
    vals = np.asarray(values, dtype=int).ravel()
    body = ", ".join(str(int(v)) for v in vals)
    return [f"    integer, parameter :: {name}({vals.size}) = [{body}]"]


ls_fd = least_squares(residual, [-2.0, 1.0], method="lm", jac="2-point")
ls_3 = least_squares(residual, [-2.0, 1.0], method="lm", jac="3-point")
ls_cs = least_squares(residual, [-2.0, 1.0], method="lm", jac="cs")
ls_jac = least_squares(residual, [-2.0, 1.0], method="lm", jac=residual_jac)
ls_limit = least_squares(residual, [-2.0, 1.0], method="lm", jac=residual_jac, max_nfev=2)
ls_scale = least_squares(
    residual, [-2.0, 1.0], method="lm", jac=residual_jac, x_scale=np.array([2.0, 0.5])
)

trf_fd = least_squares(residual, [-2.0, 1.0], method="trf", jac="2-point")
trf_3 = least_squares(residual, [-2.0, 1.0], method="trf", jac="3-point")
trf_cs = least_squares(residual, [-2.0, 1.0], method="trf", jac="cs")
trf_jac = least_squares(residual, [-2.0, 1.0], method="trf", jac=residual_jac)
trf_lsmr_fd = least_squares(
    residual, [-2.0, 1.0], method="trf", jac="2-point", tr_solver="lsmr"
)
trf_lsmr_jac = least_squares(
    residual, [-2.0, 1.0], method="trf", jac=residual_jac, tr_solver="lsmr"
)
trf_lsmr_noreg = least_squares(
    residual, [-2.0, 1.0], method="trf", jac=residual_jac, tr_solver="lsmr",
    tr_options={"regularize": False}
)
trf_operator = least_squares(
    residual, [-2.0, 1.0], method="trf", jac=residual_jac_operator
)

sparse_calls = 0
trf_sparse_2 = least_squares(
    sparse_residual, SPARSE_X0, method="trf", jac="2-point",
    jac_sparsity=SPARSE_PATTERN
)
trf_sparse_2_calls = sparse_calls
sparse_calls = 0
trf_sparse_3 = least_squares(
    sparse_residual, SPARSE_X0, method="trf", jac="3-point",
    jac_sparsity=SPARSE_PATTERN
)
trf_sparse_3_calls = sparse_calls
sparse_calls = 0
trf_sparse_cs = least_squares(
    sparse_residual, SPARSE_X0, method="trf", jac="cs",
    jac_sparsity=SPARSE_PATTERN
)
trf_sparse_cs_calls = sparse_calls
sparse_calls = 0
dog_sparse_2 = least_squares(
    sparse_residual, SPARSE_X0, method="dogbox", jac="2-point",
    jac_sparsity=SPARSE_PATTERN
)
dog_sparse_2_calls = sparse_calls
sparse_calls = 0
dog_sparse_3 = least_squares(
    sparse_residual, SPARSE_X0, method="dogbox", jac="3-point",
    jac_sparsity=SPARSE_PATTERN
)
dog_sparse_3_calls = sparse_calls
trf_bound = least_squares(
    residual, [-2.0, 1.0], method="trf", jac=residual_jac,
    bounds=([-3.0, -1.0], [0.5, 2.0])
)
trf_lsmr_bound = least_squares(
    residual, [-2.0, 1.0], method="trf", jac=residual_jac,
    bounds=([-3.0, -1.0], [0.5, 2.0]), tr_solver="lsmr"
)
trf_operator_bound = least_squares(
    residual, [-2.0, 1.0], method="trf", jac=residual_jac_operator,
    bounds=([-3.0, -1.0], [0.5, 2.0])
)
trf_bound_3 = least_squares(
    residual, [-2.0, 1.0], method="trf", jac="3-point",
    bounds=([-3.0, -1.0], [0.5, 2.0])
)
trf_bound_cs = least_squares(
    residual, [-2.0, 1.0], method="trf", jac="cs",
    bounds=([-3.0, -1.0], [0.5, 2.0])
)
trf_scale = least_squares(
    residual, [-2.0, 1.0], method="trf", jac=residual_jac,
    x_scale=np.array([2.0, 0.5])
)
trf_jac_scale = least_squares(
    residual, [-2.0, 1.0], method="trf", jac=residual_jac, x_scale="jac"
)
trf_bound_jac_scale = least_squares(
    residual, [-2.0, 1.0], method="trf", jac=residual_jac, x_scale="jac",
    bounds=([-3.0, -1.0], [0.5, 2.0])
)
trf_limit = least_squares(
    residual, [-2.0, 1.0], method="trf", jac=residual_jac, max_nfev=2
)

trf_callback = least_squares(
    residual, [-2.0, 1.0], method="trf", jac=residual_jac, callback=stop_after_two
)

dog_fd = least_squares(residual, [-2.0, 1.0], method="dogbox", jac="2-point")
dog_3 = least_squares(residual, [-2.0, 1.0], method="dogbox", jac="3-point")
dog_cs = least_squares(residual, [-2.0, 1.0], method="dogbox", jac="cs")
dog_jac = least_squares(residual, [-2.0, 1.0], method="dogbox", jac=residual_jac)
dog_lsmr_fd = least_squares(
    residual, [-2.0, 1.0], method="dogbox", jac="2-point", tr_solver="lsmr"
)
dog_lsmr_jac = least_squares(
    residual, [-2.0, 1.0], method="dogbox", jac=residual_jac, tr_solver="lsmr"
)
dog_bound = least_squares(
    residual, [-2.0, 1.0], method="dogbox", jac=residual_jac,
    bounds=([-3.0, -1.0], [0.5, 2.0])
)
dog_lsmr_bound = least_squares(
    residual, [-2.0, 1.0], method="dogbox", jac=residual_jac,
    bounds=([-3.0, -1.0], [0.5, 2.0]), tr_solver="lsmr"
)
dog_bound_3 = least_squares(
    residual, [-2.0, 1.0], method="dogbox", jac="3-point",
    bounds=([-3.0, -1.0], [0.5, 2.0])
)
dog_bound_cs = least_squares(
    residual, [-2.0, 1.0], method="dogbox", jac="cs",
    bounds=([-3.0, -1.0], [0.5, 2.0])
)
dog_scale = least_squares(
    residual, [-2.0, 1.0], method="dogbox", jac=residual_jac,
    x_scale=np.array([2.0, 0.5])
)
dog_jac_scale = least_squares(
    residual, [-2.0, 1.0], method="dogbox", jac=residual_jac, x_scale="jac"
)
dog_limit = least_squares(
    residual, [-2.0, 1.0], method="dogbox", jac=residual_jac, max_nfev=2
)

dog_callback = least_squares(
    residual, [-2.0, 1.0], method="dogbox", jac=residual_jac, callback=stop_after_two
)

robust_trf_soft = least_squares(
    robust_residual, [0.0], method="trf", jac=robust_jac,
    loss="soft_l1", f_scale=0.5
)
robust_trf_cauchy = least_squares(
    robust_residual, [0.0], method="trf", jac=robust_jac,
    loss="cauchy", f_scale=0.5
)
robust_dog_huber = least_squares(
    robust_residual, [0.0], method="dogbox", jac=robust_jac,
    loss="huber", f_scale=0.5
)
robust_dog_arctan = least_squares(
    robust_residual, [0.0], method="dogbox", jac=robust_jac,
    loss="arctan", f_scale=0.5
)
robust_trf_bound = least_squares(
    robust_residual, [0.2], method="trf", jac=robust_jac,
    loss="soft_l1", f_scale=0.5, bounds=([0.0], [1.0])
)

lq_fd = leastsq(residual, [-2.0, 1.0], full_output=True)
lq_jac = leastsq(residual, [-2.0, 1.0], Dfun=residual_jac, full_output=True)

xdata = np.array([0.0, 1.0, 2.0, 3.0])
ydata = np.array([1.1, 2.9, 5.2, 6.8])
sigma = np.array([1.0, 2.0, 1.0, 2.0])
cf_plain = curve_fit(model, xdata, ydata, p0=[2.0, 1.0], method="lm", full_output=True)
cf_rel = curve_fit(
    model, xdata, ydata, p0=[2.0, 1.0], sigma=sigma, method="lm", full_output=True
)
cf_abs = curve_fit(
    model, xdata, ydata, p0=[2.0, 1.0], sigma=sigma, absolute_sigma=True,
    method="lm", full_output=True
)
cf_jac = curve_fit(
    model, xdata, ydata, p0=[2.0, 1.0], sigma=sigma, absolute_sigma=True,
    method="lm", jac=model_jac, full_output=True
)

cf_bound = curve_fit(
    model, xdata, ydata, p0=[1.0, 1.0], bounds=([0.0, 0.0], [1.5, 3.0]),
    full_output=True
)
cf_bound_jac = curve_fit(
    model, xdata, ydata, p0=[1.0, 1.0], bounds=([0.0, 0.0], [1.5, 3.0]),
    method="trf", jac=model_jac, full_output=True
)
cf_bound_dog = curve_fit(
    model, xdata, ydata, p0=[1.0, 1.0], bounds=([0.0, 0.0], [1.5, 3.0]),
    method="dogbox", jac=model_jac, full_output=True
)
cf_trf = curve_fit(
    model, xdata, ydata, p0=[1.0, 1.0], method="trf", jac=model_jac,
    full_output=True
)
cf_trf_3 = curve_fit(
    model, xdata, ydata, p0=[1.0, 1.0], method="trf", jac="3-point",
    full_output=True
)
cf_trf_cs = curve_fit(
    model, xdata, ydata, p0=[1.0, 1.0], method="trf", jac="cs",
    full_output=True
)
cf_bound_dog_3 = curve_fit(
    model, xdata, ydata, p0=[1.0, 1.0], bounds=([0.0, 0.0], [1.5, 3.0]),
    method="dogbox", jac="3-point", full_output=True
)
cf_bound_dog_cs = curve_fit(
    model, xdata, ydata, p0=[1.0, 1.0], bounds=([0.0, 0.0], [1.5, 3.0]),
    method="dogbox", jac="cs", full_output=True
)

ydata_outlier = np.array([1.1, 2.9, 5.2, 20.0])
cf_robust_trf = curve_fit(
    model, xdata, ydata_outlier, p0=[2.0, 1.0], method="trf", jac=model_jac,
    loss="soft_l1", f_scale=0.5, full_output=True
)
cf_robust_dog = curve_fit(
    model, xdata, ydata_outlier, p0=[2.0, 1.0], method="dogbox", jac=model_jac,
    loss="soft_l1", f_scale=0.5, full_output=True
)

lines = [
    "! SPDX-License-Identifier: MIT",
    "! Generated by tools/generate_optimize_least_squares_reference.py",
    f"! SciPy version: {scipy.__version__}",
    "module test_optimize_least_squares_reference",
    "    use scifort_kinds, only : dp",
    "    implicit none",
]

for prefix, r in [
    ("ls_fd", ls_fd), ("ls_3", ls_3), ("ls_cs", ls_cs), ("ls_jac", ls_jac),
    ("ls_limit", ls_limit), ("ls_scale", ls_scale),
]:
    lines += vec(prefix + "_x_ref", r.x)
    lines += scalar(prefix + "_cost_ref", r.cost)
    lines += vec(prefix + "_fun_ref", r.fun)
    lines += mat(prefix + "_jac_ref", r.jac)
    lines += vec(prefix + "_grad_ref", r.grad)
    lines += scalar(prefix + "_optimality_ref", r.optimality)
    lines += integer(prefix + "_nfev_ref", r.nfev)
    lines += integer(prefix + "_njev_ref", -1 if r.njev is None else r.njev)
    lines += integer(prefix + "_status_ref", r.status)

for prefix, r in [
    ("trf_fd", trf_fd), ("trf_3", trf_3), ("trf_cs", trf_cs), ("trf_jac", trf_jac),
    ("trf_lsmr_fd", trf_lsmr_fd), ("trf_lsmr_jac", trf_lsmr_jac),
    ("trf_lsmr_noreg", trf_lsmr_noreg), ("trf_lsmr_bound", trf_lsmr_bound),
    ("trf_bound", trf_bound), ("trf_bound_3", trf_bound_3),
    ("trf_bound_cs", trf_bound_cs),
    ("trf_scale", trf_scale), ("trf_jac_scale", trf_jac_scale),
    ("trf_bound_jac_scale", trf_bound_jac_scale), ("trf_limit", trf_limit),
    ("trf_callback", trf_callback),
    ("trf_sparse_2", trf_sparse_2), ("trf_sparse_3", trf_sparse_3),
    ("trf_sparse_cs", trf_sparse_cs),
    ("dog_sparse_2", dog_sparse_2), ("dog_sparse_3", dog_sparse_3),
    ("dog_fd", dog_fd), ("dog_3", dog_3), ("dog_cs", dog_cs), ("dog_jac", dog_jac),
    ("dog_lsmr_fd", dog_lsmr_fd), ("dog_lsmr_jac", dog_lsmr_jac),
    ("dog_bound", dog_bound), ("dog_lsmr_bound", dog_lsmr_bound), ("dog_bound_3", dog_bound_3),
    ("dog_bound_cs", dog_bound_cs),
    ("dog_scale", dog_scale), ("dog_jac_scale", dog_jac_scale),
    ("dog_limit", dog_limit), ("dog_callback", dog_callback),
    ("robust_trf_soft", robust_trf_soft), ("robust_trf_cauchy", robust_trf_cauchy),
    ("robust_dog_huber", robust_dog_huber), ("robust_dog_arctan", robust_dog_arctan),
    ("robust_trf_bound", robust_trf_bound),
]:
    lines += vec(prefix + "_x_ref", r.x)
    lines += scalar(prefix + "_cost_ref", r.cost)
    lines += vec(prefix + "_fun_ref", r.fun)
    lines += mat(prefix + "_jac_ref", r.jac)
    lines += vec(prefix + "_grad_ref", r.grad)
    lines += scalar(prefix + "_optimality_ref", r.optimality)
    lines += ivec(prefix + "_active_ref", r.active_mask)
    lines += integer(prefix + "_nfev_ref", r.nfev)
    lines += integer(prefix + "_njev_ref", -1 if r.njev is None else r.njev)
    lines += integer(prefix + "_status_ref", r.status)

lines += integer("trf_sparse_2_calls_ref", trf_sparse_2_calls)
lines += integer("trf_sparse_3_calls_ref", trf_sparse_3_calls)
lines += integer("trf_sparse_cs_calls_ref", trf_sparse_cs_calls)
lines += integer("dog_sparse_2_calls_ref", dog_sparse_2_calls)
lines += integer("dog_sparse_3_calls_ref", dog_sparse_3_calls)

for prefix, r in [("trf_operator", trf_operator), ("trf_operator_bound", trf_operator_bound)]:
    lines += vec(prefix + "_x_ref", r.x)
    lines += scalar(prefix + "_cost_ref", r.cost)
    lines += vec(prefix + "_fun_ref", r.fun)
    lines += vec(prefix + "_grad_ref", r.grad)
    lines += scalar(prefix + "_optimality_ref", r.optimality)
    lines += ivec(prefix + "_active_ref", r.active_mask)
    lines += integer(prefix + "_nfev_ref", r.nfev)
    lines += integer(prefix + "_njev_ref", -1 if r.njev is None else r.njev)
    lines += integer(prefix + "_status_ref", r.status)

for prefix, r in [("lq_fd", lq_fd), ("lq_jac", lq_jac)]:
    x, cov, info, _msg, ier = r
    lines += vec(prefix + "_x_ref", x)
    lines += mat(prefix + "_cov_ref", cov)
    lines += vec(prefix + "_fvec_ref", info["fvec"])
    lines += integer(prefix + "_nfev_ref", info["nfev"])
    lines += integer(prefix + "_njev_ref", info.get("njev", -1))
    lines += integer(prefix + "_status_ref", ier)

for prefix, r in [
    ("cf_plain", cf_plain), ("cf_rel", cf_rel), ("cf_abs", cf_abs),
    ("cf_jac", cf_jac), ("cf_bound", cf_bound), ("cf_bound_jac", cf_bound_jac),
    ("cf_bound_dog", cf_bound_dog), ("cf_trf", cf_trf),
    ("cf_trf_3", cf_trf_3), ("cf_trf_cs", cf_trf_cs),
    ("cf_bound_dog_3", cf_bound_dog_3), ("cf_bound_dog_cs", cf_bound_dog_cs),
    ("cf_robust_trf", cf_robust_trf), ("cf_robust_dog", cf_robust_dog),
]:
    popt, pcov, info, _msg, ier = r
    lines += vec(prefix + "_popt_ref", popt)
    lines += mat(prefix + "_pcov_ref", pcov)
    lines += vec(prefix + "_fvec_ref", info["fvec"])
    lines += integer(prefix + "_nfev_ref", info["nfev"])
    lines += integer(prefix + "_status_ref", ier)

lines.append("end module test_optimize_least_squares_reference")
OUT.write_text("\n".join(lines) + "\n", encoding="ascii")
