! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright (c) 2021 Jacob Williams
! Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.

module scifort_optimize
    use scifort_kinds, only : dp
    use scifort_optimize_shgo, only : shgo_core, shgo_result, &
        shgo_objective, shgo_local_minimizer, shgo_callback, shgo_constraint_function, &
        shgo_constraint_jacobian, shgo_nonlinear_constraint, &
        shgo_constraint_blocks_valid, &
        shgo_completed, shgo_maxfun, shgo_target, shgo_callback_stop, &
        shgo_mesh_limit, shgo_hgr_stable, shgo_infeasible, shgo_invalid, shgo_numerical
    use scifort_optimize_shgo_minimize, only : shgo
    use scifort_optimize_direct, only : direct, direct_result, &
        direct_objective, direct_callback, direct_maxfun, direct_maxiter, &
        direct_vol_tol, direct_len_tol, direct_target, direct_callback_stop, &
        direct_invalid, direct_numerical
    use scifort_optimize_dual_annealing, only : dual_annealing_core, &
        dual_annealing_result, da_objective, da_local_minimizer, da_callback, &
        da_maxiter, da_maxfun, da_callback_stop, da_invalid, da_no_initial
    use scifort_optimize_dual_annealing_minimize, only : dual_annealing
    use scifort_optimize_basinhopping, only : basinhopping_result, basinhopping_core, &
        bh_objective, bh_local_minimizer, bh_take_step, bh_accept_test, bh_callback, &
        bh_completed, bh_no_improvement, bh_callback_stop, bh_invalid, bh_local_failure
    use scifort_optimize_basinhopping_minimize, only : basinhopping
    use scifort_optimize_differential_evolution, only : differential_evolution, &
        differential_evolution_result, de_objective, de_callback, &
        de_constraint_function, de_polisher, &
        de_status_converged, de_status_maxiter, de_status_callback, &
        de_status_maxfun, de_status_invalid, de_status_infeasible
    use scifort_optimize_differential_evolution_polish, only : de_lbfgsb_polisher
    use scifort_optimize_milp, only : milp, milp_result, &
        milp_optimal, milp_node_limit, milp_infeasible, milp_unbounded, &
        milp_numerical, milp_invalid
    use scifort_optimize_linprog, only : linprog, linprog_result, &
        linprog_optimal, linprog_iteration_limit, linprog_infeasible, &
        linprog_unbounded, linprog_numerical, linprog_invalid
    use scifort_optimize_linear_least_squares, only : nnls, lsq_linear, linear_lsq_result, &
        linear_lsq_success, linear_lsq_invalid, linear_lsq_maxiter, linear_lsq_numerical
    use scifort_optimize_linear_least_squares_slsqp, only : nnls_slsqp, lsq_linear_slsqp
    use scifort_optimize_assignment, only : assignment_status_success, &
        assignment_status_invalid_input, assignment_status_infeasible, linear_sum_assignment
    use scifort_optimize_quadratic_assignment, only : quadratic_assignment, &
        quadratic_assignment_score, qap_result, qap_success, qap_invalid_input, &
        qap_assignment_failure, qap_numerical_failure
    use scifort_optimize_isotonic, only : isotonic_regression, isotonic_result, &
        isotonic_status_success, isotonic_status_invalid_input, isotonic_status_numerical_failure
    use scifort_optimize_elementwise_batch, only : elementwise_argument, &
        elementwise_argument_from, elementwise_typed_argument, &
        elementwise_typed_argument_from
    use scifort_optimize_least_squares, only : curve_fit, curve_fit_result, &
        curve_model_complex_function, curve_model_function, curve_model_jacobian, &
        least_squares, least_squares_callback, &
        least_squares_complex_residual_function, least_squares_jacobian_function, &
        least_squares_jacobian_matvec_function, least_squares_jacobian_rmatvec_function, &
        least_squares_residual_function, &
        least_squares_result, leastsq, leastsq_result, robust_loss_function
    use scifort_optimize_minimize, only : bfgs, cg, cobyla, dogleg, lbfgsb, linear_constraint, minimize, &
        multivariate_constraint_function, multivariate_constraint_hessian, &
        multivariate_constraint_hessian_product, &
        multivariate_constraint_jacobian, multivariate_gradient_function, &
        multivariate_hessian_function, &
        multivariate_hessian_product_function, multivariate_minimize_result, &
        multivariate_objective_function, nelder_mead, newton_cg, nonlinear_constraint, &
        powell, slsqp, tnc, trust_exact, &
        trust_constr, trust_constr_iteration_callback, trust_krylov, trust_ncg
    use scifort_optimize_minimize_scalar, only : bounded, bracket, bracket_minimum, &
        bracket_minimum_vectorized, bracket_minimum_vectorized_into, &
        bracket_minimum_vectorized_args_into, bracket_minimum_vectorized_indexed_args_into, &
        bracket_minimum_result, bracket_minimum_status_limit, bracket_result, brent, &
        elementwise_minimum_vector_function, elementwise_minimum_args_function, &
        elementwise_minimum_typed_args_function, &
        elementwise_minimum_indexed_args_function, &
        find_minimum, find_minimum_vectorized, find_minimum_vectorized_args_into, &
        find_minimum_vectorized_into, find_minimum_vectorized_indexed_args_into, &
        find_minimum_callback, find_minimum_result, &
        find_minimum_status_callback, golden, &
        minimize_result, minimize_scalar, minimize_status_bracket_error, &
        minimize_status_convergence_error, minimize_status_invalid_input, &
        minimize_status_success, minimize_status_value_error, scalar_objective_function
    use scifort_optimize_root, only : anderson, broyden1, broyden2, diagbroyden, &
        excitingmixing, fsolve, linearmixing, newton_krylov, root, &
        vector_root_function, &
        vector_root_jacobian, vector_root_norm_function, vector_root_result, &
        vector_root_linear_operator, vector_root_preconditioner_update, &
        vector_root_inner_solver, vector_root_jacobian_product
    use scifort_optimize_root_scalar, only : bisect, bracket_root, bracket_root_result, &
        bracket_root_vectorized, bracket_root_vectorized_into, &
        bracket_root_vectorized_args_into, bracket_root_vectorized_indexed_args_into, &
        bracket_root_status_limit, brenth, brentq, find_root, find_root_callback, &
        find_root_result, find_root_status_callback, find_root_vectorized, &
        find_root_vectorized_into, find_root_vectorized_args_into, &
        find_root_vectorized_indexed_args_into, &
        elementwise_root_vector_function, elementwise_root_args_function, &
        elementwise_root_typed_args_function, &
        elementwise_root_indexed_args_function, &
        newton, ridder, &
        root_result, root_scalar, root_status_convergence_error, root_status_invalid_input, &
        root_status_sign_error, root_status_success, root_status_value_error, &
        scalar_root_function, toms748
    implicit none
    private

    public :: shgo, shgo_core, shgo_result, shgo_objective, shgo_local_minimizer
    public :: shgo_callback, shgo_constraint_function, shgo_constraint_jacobian
    public :: shgo_nonlinear_constraint, shgo_constraint_blocks_valid
    public :: shgo_completed, shgo_maxfun, shgo_target
    public :: shgo_callback_stop, shgo_mesh_limit, shgo_hgr_stable, shgo_infeasible, shgo_invalid, shgo_numerical
    public :: direct, direct_result, direct_objective, direct_callback
    public :: direct_maxfun, direct_maxiter, direct_vol_tol, direct_len_tol
    public :: direct_target, direct_callback_stop, direct_invalid, direct_numerical
    public :: dual_annealing, dual_annealing_core, dual_annealing_result
    public :: da_objective, da_local_minimizer, da_callback
    public :: da_maxiter, da_maxfun, da_callback_stop, da_invalid, da_no_initial
    public :: basinhopping, basinhopping_core, basinhopping_result
    public :: bh_objective, bh_local_minimizer, bh_take_step, bh_accept_test, bh_callback
    public :: bh_completed, bh_no_improvement, bh_callback_stop, bh_invalid, bh_local_failure
    public :: differential_evolution, differential_evolution_result
    public :: de_objective, de_callback, de_status_converged
    public :: de_status_maxiter, de_status_callback, de_status_maxfun, de_status_invalid
    public :: de_status_infeasible, de_constraint_function, de_polisher
    public :: de_lbfgsb_polisher
    public :: milp, milp_result
    public :: milp_optimal, milp_node_limit, milp_infeasible, milp_unbounded, milp_numerical, milp_invalid
    public :: linprog, linprog_result
    public :: linprog_optimal, linprog_iteration_limit, linprog_infeasible
    public :: linprog_unbounded, linprog_numerical, linprog_invalid
    public :: nnls, lsq_linear, linear_lsq_result
    public :: linear_lsq_success, linear_lsq_invalid, linear_lsq_maxiter, linear_lsq_numerical
    public :: nnls_slsqp, lsq_linear_slsqp
    public :: quadratic_assignment, quadratic_assignment_score, qap_result
    public :: qap_success, qap_invalid_input, qap_assignment_failure, qap_numerical_failure
    public :: linear_sum_assignment
    public :: assignment_status_success, assignment_status_invalid_input, assignment_status_infeasible
    public :: isotonic_regression, isotonic_result
    public :: isotonic_status_success, isotonic_status_invalid_input
    public :: isotonic_status_numerical_failure
    public :: elementwise_argument
    public :: elementwise_argument_from
    public :: elementwise_typed_argument
    public :: elementwise_typed_argument_from
    public :: curve_fit
    public :: curve_fit_result
    public :: curve_model_complex_function
    public :: curve_model_function
    public :: curve_model_jacobian
    public :: least_squares
    public :: least_squares_callback
    public :: least_squares_complex_residual_function
    public :: least_squares_jacobian_function
    public :: least_squares_jacobian_matvec_function
    public :: least_squares_jacobian_rmatvec_function
    public :: least_squares_residual_function
    public :: least_squares_result
    public :: robust_loss_function
    public :: leastsq
    public :: leastsq_result
    public :: bfgs
    public :: cg
    public :: cobyla
    public :: dogleg
    public :: lbfgsb
    public :: linear_constraint
    public :: minimize
    public :: multivariate_constraint_function
    public :: multivariate_constraint_hessian
    public :: multivariate_constraint_hessian_product
    public :: multivariate_constraint_jacobian
    public :: multivariate_gradient_function
    public :: multivariate_hessian_function
    public :: multivariate_hessian_product_function
    public :: multivariate_minimize_result
    public :: multivariate_objective_function
    public :: nelder_mead
    public :: nonlinear_constraint
    public :: newton_cg
    public :: trust_constr
    public :: trust_constr_iteration_callback
    public :: trust_exact
    public :: trust_krylov
    public :: trust_ncg
    public :: powell
    public :: slsqp
    public :: tnc
    public :: anderson
    public :: broyden1
    public :: broyden2
    public :: diagbroyden
    public :: excitingmixing
    public :: fsolve
    public :: linearmixing
    public :: newton_krylov
    public :: root
    public :: vector_root_function
    public :: vector_root_jacobian
    public :: vector_root_norm_function
    public :: vector_root_linear_operator
    public :: vector_root_preconditioner_update
    public :: vector_root_inner_solver
    public :: vector_root_jacobian_product
    public :: vector_root_result
    public :: bounded
    public :: bracket
    public :: bracket_minimum
    public :: bracket_minimum_vectorized
    public :: bracket_minimum_vectorized_args_into
    public :: bracket_minimum_vectorized_indexed_args_into
    public :: bracket_minimum_vectorized_into
    public :: bracket_minimum_result
    public :: bracket_minimum_status_limit
    public :: bracket_result
    public :: find_minimum
    public :: find_minimum_vectorized
    public :: find_minimum_vectorized_args_into
    public :: find_minimum_vectorized_indexed_args_into
    public :: find_minimum_vectorized_into
    public :: elementwise_minimum_args_function
    public :: elementwise_minimum_typed_args_function
    public :: elementwise_minimum_indexed_args_function
    public :: elementwise_minimum_vector_function
    public :: find_minimum_callback
    public :: find_minimum_result
    public :: find_minimum_status_callback
    public :: brent
    public :: golden
    public :: minimize_result
    public :: minimize_scalar
    public :: minimize_status_bracket_error
    public :: minimize_status_convergence_error
    public :: minimize_status_invalid_input
    public :: minimize_status_success
    public :: minimize_status_value_error
    public :: scalar_objective_function
    public :: bisect
    public :: bracket_root
    public :: bracket_root_vectorized
    public :: bracket_root_vectorized_args_into
    public :: bracket_root_vectorized_indexed_args_into
    public :: bracket_root_vectorized_into
    public :: bracket_root_result
    public :: bracket_root_status_limit
    public :: find_root
    public :: find_root_vectorized
    public :: find_root_vectorized_args_into
    public :: find_root_vectorized_indexed_args_into
    public :: find_root_vectorized_into
    public :: elementwise_root_args_function
    public :: elementwise_root_typed_args_function
    public :: elementwise_root_indexed_args_function
    public :: elementwise_root_vector_function
    public :: find_root_callback
    public :: find_root_result
    public :: find_root_status_callback
    public :: brenth
    public :: brentq
    public :: dp
    public :: newton
    public :: ridder
    public :: root_result
    public :: root_scalar
    public :: root_status_convergence_error
    public :: root_status_invalid_input
    public :: root_status_sign_error
    public :: root_status_success
    public :: root_status_value_error
    public :: scalar_root_function
    public :: toms748

end module scifort_optimize
