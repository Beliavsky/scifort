#!/usr/bin/env python3
"""Generate SciPy 1.17 reference values for trust-constr regressions."""

from pathlib import Path

import numpy as np
import scipy
from scipy import optimize


def fmt(value):
    return f"{float(value):.17e}_dp"


def emit_vector(parts, name, values):
    values = np.asarray(values, dtype=float).ravel()
    parts.append(f"    real(dp), parameter :: {name}({values.size}) = [ &\n")
    for i, value in enumerate(values):
        suffix = ", &" if i + 1 < values.size else "]"
        parts.append(f"        {fmt(value)}{suffix}\n")


def emit_matrix(parts, name, values):
    values = np.asarray(values, dtype=float)
    rows, cols = values.shape
    flat = values.ravel(order="F")
    parts.append(
        f"    real(dp), parameter :: {name}({rows}, {cols}) = reshape([ &\n"
    )
    for i, value in enumerate(flat):
        suffix = ", &" if i + 1 < flat.size else f"], [{rows}, {cols}])"
        parts.append(f"        {fmt(value)}{suffix}\n")


def emit_int_vector(parts, name, values):
    values = np.asarray(values, dtype=int).ravel()
    joined = ", ".join(str(int(value)) for value in values)
    parts.append(f"    integer, parameter :: {name}({values.size}) = [{joined}]\n")


def emit_constraint_diagnostics(parts, prefix, result):
    parts.append(
        f"    character(len={len(result.method)}), parameter :: {prefix}_method_ref = "
        f"'{result.method}'\n"
    )
    emit_vector(parts, f"{prefix}_constr_ref", np.concatenate(result.constr))
    jac_blocks = [jac.toarray() if hasattr(jac, "toarray") else np.asarray(jac) for jac in result.jac]
    emit_matrix(parts, f"{prefix}_constr_jac_ref", np.vstack(jac_blocks))
    emit_int_vector(parts, f"{prefix}_constr_nfev_ref", result.constr_nfev)
    emit_int_vector(parts, f"{prefix}_constr_njev_ref", result.constr_njev)
    emit_int_vector(parts, f"{prefix}_constr_nhev_ref", result.constr_nhev)


def emit_result(parts, prefix, result):
    emit_vector(parts, f"{prefix}_x_ref", result.x)
    parts.append(f"    real(dp), parameter :: {prefix}_fun_ref = {fmt(result.fun)}\n")
    parts.append(
        f"    real(dp), parameter :: {prefix}_maxcv_ref = "
        f"{fmt(result.constr_violation)}\n"
    )
    parts.append(f"    integer, parameter :: {prefix}_status_ref = {int(result.status)}\n")
    parts.append(
        f"    logical, parameter :: {prefix}_success_ref = "
        f"{'.true.' if result.success else '.false.'}\n"
    )


def linear_objective(x):
    return (x[0] - 1.0) ** 2 + (x[1] - 2.5) ** 2


def linear_gradient(x):
    return np.array([2.0 * (x[0] - 1.0), 2.0 * (x[1] - 2.5)])


def linear_hessian(x):
    del x
    return 2.0 * np.eye(2)


def linear_constraint(x):
    return np.array([x[0] + x[1] - 2.0])


def linear_constraint_jacobian(x):
    del x
    return np.array([[1.0, 1.0]])


def circle_objective(x):
    return (x[0] - 2.0) ** 2 + (x[1] - 1.0) ** 2


def circle_gradient(x):
    return np.array([2.0 * (x[0] - 2.0), 2.0 * (x[1] - 1.0)])


def circle_hessian(x):
    del x
    return 2.0 * np.eye(2)


def circle_constraint(x):
    return np.array([x[0] ** 2 + x[1] ** 2 - 1.0])


def circle_constraint_jacobian(x):
    return np.array([[2.0 * x[0], 2.0 * x[1]]])


def circle_constraint_hessian(x, v):
    del x
    return 2.0 * v[0] * np.eye(2)


def maratos_objective(x):
    return 2.0 * (x[0] ** 2 + x[1] ** 2 - 1.0) - x[0]


def maratos_gradient(x):
    return np.array([4.0 * x[0] - 1.0, 4.0 * x[1]])


def maratos_hessian(x):
    del x
    return 4.0 * np.eye(2)


def maratos_constraint(x):
    return np.array([x[0] ** 2 + x[1] ** 2 - 1.0])


def maratos_constraint_jacobian(x):
    return np.array([[2.0 * x[0], 2.0 * x[1]]])


def maratos_constraint_hessian(x, v):
    del x
    return 2.0 * v[0] * np.eye(2)


def active_inequality(x):
    return np.array([1.0 - x[0] - x[1]])


def active_inequality_jacobian(x):
    del x
    return np.array([[-1.0, -1.0]])


def circle_inequality(x):
    return np.array([1.0 - x[0] ** 2 - x[1] ** 2])


def circle_inequality_jacobian(x):
    return np.array([[-2.0 * x[0], -2.0 * x[1]]])


def circle_inequality_hessian(x, v):
    del x
    return -2.0 * v[0] * np.eye(2)


def bound_objective(x):
    return (x[0] + 1.0) ** 2 + (x[1] - 3.0) ** 2


def bound_gradient(x):
    return np.array([2.0 * (x[0] + 1.0), 2.0 * (x[1] - 3.0)])


def mixed_constraints(x):
    return np.array([x[0] + x[1] - 2.0, x[0] - 0.5])


def mixed_constraint_jacobian(x):
    del x
    return np.array([[1.0, 1.0], [1.0, 0.0]])


def bounded_constraints(x):
    return np.array([x[0] + x[1], x[0], x[1], x[0] - x[1]])


def bounded_constraint_jacobian(x):
    del x
    return np.array([[1.0, 1.0], [1.0, 0.0], [0.0, 1.0], [1.0, -1.0]])


def keep_feasible_objective(x):
    return (x[0] + 2.0) ** 2


def keep_feasible_gradient(x):
    return np.array([2.0 * (x[0] + 2.0)])


def keep_feasible_hessian(x):
    del x
    return np.array([[2.0]])


def keep_feasible_constraint(x):
    return np.array([x[0]])


def keep_feasible_constraint_jacobian(x):
    del x
    return np.array([[1.0]])


def fd_clip_objective(x):
    return (x[0] - 2.0) ** 2


def fd_clip_constraint(x):
    return np.array([x[0] ** 2])


parts = [
    "! SPDX-License-Identifier: MIT\n",
    "! Generated by tools/generate_optimize_trust_constr_reference.py\n",
    f"! SciPy version: {scipy.__version__}\n",
    "module test_optimize_trust_constr_reference\n",
    "    use scifort_kinds, only : dp\n",
    "    implicit none\n",
]

linear = optimize.NonlinearConstraint(
    linear_constraint,
    0.0,
    0.0,
    jac=linear_constraint_jacobian,
    hess=lambda x, v: np.zeros((2, 2)),
)
emit_result(
    parts,
    "linear",
    optimize.minimize(
        linear_objective,
        np.array([0.5, 1.5]),
        method="trust-constr",
        jac=linear_gradient,
        hess=linear_hessian,
        constraints=[linear],
        options={"gtol": 1.0e-10, "xtol": 1.0e-10},
    ),
)

penalty_result = optimize.minimize(
    linear_objective,
    np.array([0.5, 1.5]),
    method="trust-constr",
    jac=linear_gradient,
    hess=linear_hessian,
    constraints=[linear],
    options={
        "gtol": 1.0e-10,
        "xtol": 1.0e-10,
        "initial_constr_penalty": 7.5,
    },
)
parts.append(
    "    real(dp), parameter :: initial_penalty_ref = "
    f"{fmt(penalty_result.constr_penalty)}\n"
)

infeasible_xtol_result = optimize.minimize(
    linear_objective,
    np.array([0.0, 0.0]),
    method="trust-constr",
    jac=linear_gradient,
    hess=linear_hessian,
    constraints=[linear],
    options={"gtol": 1.0e-8, "xtol": 2.0, "maxiter": 20},
)
parts.append(
    "    integer, parameter :: infeasible_xtol_status_ref = "
    f"{int(infeasible_xtol_result.status)}\n"
)
parts.append(
    "    logical, parameter :: infeasible_xtol_success_ref = "
    f"{'.true.' if infeasible_xtol_result.success else '.false.'}\n"
)
parts.append(
    "    real(dp), parameter :: infeasible_xtol_maxcv_ref = "
    f"{fmt(infeasible_xtol_result.constr_violation)}\n"
)
parts.append(
    "    integer, parameter :: infeasible_xtol_nit_ref = "
    f"{int(infeasible_xtol_result.nit)}\n"
)

circle = optimize.NonlinearConstraint(
    circle_constraint,
    0.0,
    0.0,
    jac=circle_constraint_jacobian,
)
emit_result(
    parts,
    "circle",
    optimize.minimize(
        circle_objective,
        np.array([1.0, 0.2]),
        method="trust-constr",
        jac=circle_gradient,
        constraints=[circle],
        options={"gtol": 1.0e-9, "xtol": 1.0e-10, "maxiter": 200},
    ),
)

circle_exact = optimize.NonlinearConstraint(
    circle_constraint,
    0.0,
    0.0,
    jac=circle_constraint_jacobian,
    hess=circle_constraint_hessian,
)
circle_hess_result = optimize.minimize(
    circle_objective,
    np.array([1.0, 0.2]),
    method="trust-constr",
    jac=circle_gradient,
    hess=circle_hessian,
    constraints=[circle_exact],
    options={"gtol": 1.0e-10, "xtol": 1.0e-10, "maxiter": 200},
)
emit_result(parts, "circle_hess", circle_hess_result)

maratos = optimize.NonlinearConstraint(
    maratos_constraint,
    0.0,
    0.0,
    jac=maratos_constraint_jacobian,
    hess=maratos_constraint_hessian,
)
maratos_result = optimize.minimize(
    maratos_objective,
    np.array([0.0, 1.0]),
    method="trust-constr",
    jac=maratos_gradient,
    hess=maratos_hessian,
    constraints=[maratos],
    options={"gtol": 1.0e-10, "xtol": 1.0e-10, "maxiter": 200},
)
emit_result(parts, "maratos", maratos_result)
emit_constraint_diagnostics(parts, "maratos", maratos_result)
parts.append(f"    integer, parameter :: maratos_nfev_ref = {int(maratos_result.nfev)}\n")
emit_vector(parts, "maratos_lagrangian_grad_ref", maratos_result.lagrangian_grad)
parts.append(
    "    integer, parameter :: maratos_cg_niter_ref = "
    f"{int(maratos_result.cg_niter)}\n"
)
parts.append(
    "    integer, parameter :: maratos_cg_stop_cond_ref = "
    f"{int(maratos_result.cg_stop_cond)}\n"
)

barrier_maratos = optimize.NonlinearConstraint(
    maratos_constraint,
    0.0,
    np.inf,
    jac=maratos_constraint_jacobian,
    hess=maratos_constraint_hessian,
)
barrier_maratos_result = optimize.minimize(
    maratos_objective,
    np.array([0.0, 1.0]),
    method="trust-constr",
    jac=maratos_gradient,
    hess=maratos_hessian,
    constraints=[barrier_maratos],
    options={"gtol": 1.0e-10, "xtol": 1.0e-10, "maxiter": 500},
)
emit_result(parts, "barrier_maratos", barrier_maratos_result)
emit_constraint_diagnostics(parts, "barrier_maratos", barrier_maratos_result)
parts.append(
    "    integer, parameter :: barrier_maratos_nit_ref = "
    f"{int(barrier_maratos_result.nit)}\n"
)
parts.append(
    "    integer, parameter :: barrier_maratos_nfev_ref = "
    f"{int(barrier_maratos_result.nfev)}\n"
)
parts.append(
    "    integer, parameter :: barrier_maratos_njev_ref = "
    f"{int(barrier_maratos_result.njev)}\n"
)
parts.append(
    "    integer, parameter :: barrier_maratos_nhev_ref = "
    f"{int(barrier_maratos_result.nhev)}\n"
)
emit_vector(
    parts,
    "barrier_maratos_lagrangian_grad_ref",
    barrier_maratos_result.lagrangian_grad,
)
parts.append(
    "    integer, parameter :: barrier_maratos_cg_niter_ref = "
    f"{int(barrier_maratos_result.cg_niter)}\n"
)
parts.append(
    "    integer, parameter :: barrier_maratos_cg_stop_cond_ref = "
    f"{int(barrier_maratos_result.cg_stop_cond)}\n"
)
parts.append(
    "    real(dp), parameter :: barrier_maratos_optimality_ref = "
    f"{fmt(barrier_maratos_result.optimality)}\n"
)
parts.append(
    "    real(dp), parameter :: barrier_maratos_barrier_parameter_ref = "
    f"{fmt(barrier_maratos_result.barrier_parameter)}\n"
)
parts.append(
    "    real(dp), parameter :: barrier_maratos_barrier_tolerance_ref = "
    f"{fmt(barrier_maratos_result.barrier_tolerance)}\n"
)

ineq = optimize.NonlinearConstraint(
    active_inequality,
    0.0,
    np.inf,
    jac=active_inequality_jacobian,
)
emit_result(
    parts,
    "ineq",
    optimize.minimize(
        linear_objective,
        np.array([0.0, 0.0]),
        method="trust-constr",
        jac=linear_gradient,
        hess=linear_hessian,
        constraints=[ineq],
        options={"gtol": 1.0e-8, "xtol": 1.0e-8},
    ),
)

emit_result(
    parts,
    "barrier_options",
    optimize.minimize(
        linear_objective,
        np.array([0.0, 0.0]),
        method="trust-constr",
        jac=linear_gradient,
        hess=linear_hessian,
        constraints=[ineq],
        options={
            "gtol": 1.0e-8,
            "xtol": 1.0e-8,
            "barrier_tol": 1.0e-6,
            "initial_barrier_parameter": 5.0e-2,
            "initial_barrier_tolerance": 2.5e-2,
        },
    ),
)

circle_ineq = optimize.NonlinearConstraint(
    circle_inequality,
    0.0,
    np.inf,
    jac=circle_inequality_jacobian,
)
emit_result(
    parts,
    "circle_ineq",
    optimize.minimize(
        circle_objective,
        np.array([0.2, 0.1]),
        method="trust-constr",
        jac=circle_gradient,
        constraints=[circle_ineq],
        options={"gtol": 1.0e-8, "xtol": 1.0e-8, "maxiter": 400},
    ),
)

circle_ineq_exact = optimize.NonlinearConstraint(
    circle_inequality,
    0.0,
    np.inf,
    jac=circle_inequality_jacobian,
    hess=circle_inequality_hessian,
)
emit_result(
    parts,
    "circle_ineq_hess",
    optimize.minimize(
        circle_objective,
        np.array([0.2, 0.1]),
        method="trust-constr",
        jac=circle_gradient,
        hess=circle_hessian,
        constraints=[circle_ineq_exact],
        options={"gtol": 1.0e-8, "xtol": 1.0e-8, "maxiter": 400},
    ),
)

bounds_result = optimize.minimize(
    bound_objective,
    np.array([0.5, 1.0]),
    method="trust-constr",
    jac=bound_gradient,
    hess=linear_hessian,
    bounds=optimize.Bounds([0.0, -np.inf], [np.inf, 2.0]),
    options={"gtol": 1.0e-8, "xtol": 1.0e-8},
)
emit_result(parts, "bounds", bounds_result)
emit_constraint_diagnostics(parts, "bounds", bounds_result)

mixed_combined = optimize.NonlinearConstraint(
    mixed_constraints,
    np.array([0.0, 0.0]),
    np.array([0.0, np.inf]),
    jac=mixed_constraint_jacobian,
)
mixed_bounds_result = optimize.minimize(
    linear_objective,
    np.array([0.8, 1.2]),
    method="trust-constr",
    jac=linear_gradient,
    hess=linear_hessian,
    constraints=[mixed_combined],
    bounds=optimize.Bounds([0.0, -np.inf], [np.inf, 2.0]),
    options={"gtol": 1.0e-8, "xtol": 1.0e-8, "maxiter": 300},
)
emit_result(parts, "mixed_bounds", mixed_bounds_result)
emit_constraint_diagnostics(parts, "mixed_bounds", mixed_bounds_result)
parts.append(
    "    integer, parameter :: mixed_bounds_nit_ref = "
    f"{int(mixed_bounds_result.nit)}\n"
)

bounded_nonlinear = optimize.NonlinearConstraint(
    bounded_constraints,
    np.array([2.0, 0.5, -np.inf, -1.2]),
    np.array([2.0, np.inf, 1.6, -0.6]),
    jac=bounded_constraint_jacobian,
)
bounded_nonlinear_result = optimize.minimize(
    linear_objective,
    np.array([0.8, 1.2]),
    method="trust-constr",
    jac=linear_gradient,
    hess=linear_hessian,
    constraints=[bounded_nonlinear],
    options={"gtol": 1.0e-8, "xtol": 1.0e-8, "maxiter": 300},
)
emit_result(parts, "bounded_nonlinear", bounded_nonlinear_result)
emit_constraint_diagnostics(parts, "bounded_nonlinear", bounded_nonlinear_result)
parts.append(
    "    integer, parameter :: bounded_nonlinear_nit_ref = "
    f"{int(bounded_nonlinear_result.nit)}\n"
)
parts.append(
    "    integer, parameter :: bounded_nonlinear_nfev_ref = "
    f"{int(bounded_nonlinear_result.nfev)}\n"
)
parts.append(
    "    integer, parameter :: bounded_nonlinear_cg_niter_ref = "
    f"{int(bounded_nonlinear_result.cg_niter)}\n"
)

keep_feasible_nonlinear = optimize.NonlinearConstraint(
    keep_feasible_constraint,
    np.array([0.5]),
    np.array([np.inf]),
    jac=keep_feasible_constraint_jacobian,
    keep_feasible=np.array([True]),
)
keep_feasible_nonlinear_result = optimize.minimize(
    keep_feasible_objective,
    np.array([2.0]),
    method="trust-constr",
    jac=keep_feasible_gradient,
    hess=keep_feasible_hessian,
    constraints=[keep_feasible_nonlinear],
    options={"gtol": 1.0e-8, "xtol": 1.0e-8, "maxiter": 100},
)
emit_result(parts, "keep_feasible_nonlinear", keep_feasible_nonlinear_result)
emit_constraint_diagnostics(parts, "keep_feasible_nonlinear", keep_feasible_nonlinear_result)
for field in ("nit", "nfev", "njev", "nhev", "cg_niter"):
    parts.append(
        f"    integer, parameter :: keep_feasible_nonlinear_{field}_ref = "
        f"{int(getattr(keep_feasible_nonlinear_result, field))}\n"
    )

keep_feasible_bounds_result = optimize.minimize(
    keep_feasible_objective,
    np.array([0.6]),
    method="trust-constr",
    jac=keep_feasible_gradient,
    hess=keep_feasible_hessian,
    bounds=optimize.Bounds(
        np.array([0.5]), np.array([np.inf]), keep_feasible=np.array([True])
    ),
    options={"gtol": 1.0e-8, "xtol": 1.0e-8, "barrier_tol": 1.0e-3},
)
emit_result(parts, "keep_feasible_bounds", keep_feasible_bounds_result)
emit_constraint_diagnostics(parts, "keep_feasible_bounds", keep_feasible_bounds_result)
for field in ("nit", "nfev", "njev", "nhev", "cg_niter"):
    parts.append(
        f"    integer, parameter :: keep_feasible_bounds_{field}_ref = "
        f"{int(getattr(keep_feasible_bounds_result, field))}\n"
    )


fd_clip_constraint_obj = optimize.NonlinearConstraint(
    fd_clip_constraint,
    np.array([0.2]),
    np.array([np.inf]),
    jac="2-point",
)
fd_clip_result = optimize.minimize(
    fd_clip_objective,
    np.array([0.99999999]),
    method="trust-constr",
    jac="2-point",
    constraints=[fd_clip_constraint_obj],
    bounds=optimize.Bounds(
        np.array([0.0]), np.array([1.0]), keep_feasible=np.array([True])
    ),
    options={"gtol": 1.0e-8, "xtol": 1.0e-8, "maxiter": 100},
)
emit_result(parts, "fd_clip", fd_clip_result)
emit_constraint_diagnostics(parts, "fd_clip", fd_clip_result)
for field in ("nit", "nfev", "njev", "nhev", "cg_niter"):
    parts.append(
        f"    integer, parameter :: fd_clip_{field}_ref = "
        f"{int(getattr(fd_clip_result, field))}\n"
    )

mixed_eq = optimize.NonlinearConstraint(
    lambda x: mixed_constraints(x)[:1],
    0.0,
    0.0,
    jac=lambda x: mixed_constraint_jacobian(x)[:1, :],
)
mixed_ineq = optimize.NonlinearConstraint(
    lambda x: mixed_constraints(x)[1:],
    0.0,
    np.inf,
    jac=lambda x: mixed_constraint_jacobian(x)[1:, :],
)
emit_result(
    parts,
    "mixed",
    optimize.minimize(
        linear_objective,
        np.array([0.8, 1.2]),
        method="trust-constr",
        jac=linear_gradient,
        hess=linear_hessian,
        constraints=[mixed_eq, mixed_ineq],
        options={"gtol": 1.0e-8, "xtol": 1.0e-8},
    ),
)

emit_result(
    parts,
    "fixed_bound",
    optimize.minimize(
        lambda x: (x[0] - 1.0) ** 2 + (x[1] - 3.0) ** 2,
        np.array([0.7, 1.0]),
        method="trust-constr",
        jac=lambda x: np.array([2.0 * (x[0] - 1.0), 2.0 * (x[1] - 3.0)]),
        hess=linear_hessian,
        bounds=optimize.Bounds([0.5, -np.inf], [0.5, 2.0]),
        options={"gtol": 1.0e-8, "xtol": 1.0e-8},
    ),
)

parts.append("end module test_optimize_trust_constr_reference\n")
out = Path(__file__).resolve().parents[1] / "test" / "optimize_trust_constr_reference.f90"
out.write_text("".join(parts), encoding="ascii")
print(out)
