program test_optimize_trust_constr
    use, intrinsic :: ieee_arithmetic, only : ieee_negative_inf, ieee_positive_inf, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_minimize, only : minimize, multivariate_minimize_result, trust_constr
    use test_optimize_trust_constr_reference
    implicit none

    type(multivariate_minimize_result) :: r
    type(multivariate_minimize_result) :: r_alias
    type(multivariate_minimize_result) :: r_legacy
    integer :: failures
    integer :: constraint_hess_calls
    real(dp) :: fd_clip_objective_max
    real(dp) :: fd_clip_constraint_max

    failures = 0
    constraint_hess_calls = 0

    r = trust_constr(linear_objective, [0.5_dp, 1.5_dp], jac=linear_gradient, &
        hess=linear_hessian, constraints=linear_constraint, &
        constraint_jac=linear_constraint_jacobian, meq=1, &
        gtol=1.0e-10_dp, xtol=1.0e-10_dp)
    call check_result('linear equality analytic', r, linear_x_ref, linear_fun_ref, &
        linear_maxcv_ref, 5.0e-11_dp, 5.0e-12_dp, 5.0e-12_dp, failures)
    call check_int('linear equality status', r%status, linear_status_ref, failures)
    if (r%success .neqv. linear_success_ref) then
        write(*, '(a)') 'linear equality success mismatch'
        failures = failures + 1
    end if

    r = trust_constr(linear_objective, [0.5_dp, 1.5_dp], jac=linear_gradient, &
        hess=linear_hessian, constraints=linear_constraint, &
        constraint_jac=linear_constraint_jacobian, &
        constraint_lower=[0.0_dp], constraint_upper=[0.0_dp], &
        gtol=1.0e-10_dp, xtol=1.0e-10_dp)
    call check_result('bounded-form equality', r, linear_x_ref, linear_fun_ref, &
        linear_maxcv_ref, 5.0e-11_dp, 5.0e-12_dp, 5.0e-12_dp, failures)
    if (trim(r%method) /= 'equality_constrained_sqp') then
        write(*, '(a,a)') 'bounded-form equality method mismatch: ', trim(r%method)
        failures = failures + 1
    end if

    r = trust_constr(linear_objective, [0.5_dp, 1.5_dp], jac=linear_gradient, &
        hess=linear_hessian, constraints=linear_constraint, &
        constraint_jac=linear_constraint_jacobian, meq=1, &
        initial_constr_penalty=7.5_dp, gtol=1.0e-10_dp, xtol=1.0e-10_dp)
    call check_result('initial constraint penalty', r, linear_x_ref, linear_fun_ref, &
        linear_maxcv_ref, 5.0e-11_dp, 5.0e-12_dp, 5.0e-12_dp, failures)
    call check_real('initial constraint penalty result', r%constr_penalty, &
        initial_penalty_ref, 5.0e-14_dp, failures)

    r = trust_constr(linear_objective, [0.0_dp, 0.0_dp], jac=linear_gradient, &
        hess=linear_hessian, constraints=linear_constraint, &
        constraint_jac=linear_constraint_jacobian, meq=1, &
        gtol=1.0e-8_dp, xtol=2.0_dp, maxiter=20)
    call check_int('infeasible xtol status', r%status, infeasible_xtol_status_ref, failures)
    call check_int('infeasible xtol iterations', r%nit, infeasible_xtol_nit_ref, failures)
    call check_real('infeasible xtol constraint violation', r%constr_violation, &
        infeasible_xtol_maxcv_ref, 5.0e-14_dp, failures)
    if (r%success .neqv. infeasible_xtol_success_ref) then
        write(*, '(a)') 'infeasible xtol success mismatch'
        failures = failures + 1
    end if

    r = minimize(linear_objective, [0.5_dp, 1.5_dp], method='trust-constr', &
        jac=linear_gradient, hess=linear_hessian, constraints=linear_constraint, &
        constraint_jac=linear_constraint_jacobian, meq=1, &
        gtol=1.0e-10_dp, xtol=1.0e-10_dp)
    call check_result('minimize dispatch', r, linear_x_ref, linear_fun_ref, &
        linear_maxcv_ref, 5.0e-11_dp, 5.0e-12_dp, 5.0e-12_dp, failures)
    if (trim(r%method) /= trim(maratos_method_ref)) then
        write(*, '(a,a)') 'minimize dispatch method mismatch: ', trim(r%method)
        failures = failures + 1
    end if

    r_legacy = trust_constr(linear_objective, [0.5_dp, 1.5_dp], &
        jac=linear_gradient, hess=linear_hessian, constraints=linear_constraint, &
        constraint_jac=linear_constraint_jacobian, meq=1, &
        initial_trust_radius=0.2_dp, gtol=1.0e-12_dp, xtol=1.0e-12_dp, maxiter=3)
    r_alias = trust_constr(linear_objective, [0.5_dp, 1.5_dp], &
        jac=linear_gradient, hess=linear_hessian, constraints=linear_constraint, &
        constraint_jac=linear_constraint_jacobian, meq=1, &
        initial_tr_radius=0.2_dp, gtol=1.0e-12_dp, xtol=1.0e-12_dp, maxiter=3)
    call check_result('initial_tr_radius alias', r_alias, r_legacy%x, r_legacy%fun, &
        r_legacy%maxcv, 5.0e-14_dp, 5.0e-14_dp, 5.0e-14_dp, failures)
    call check_int('initial_tr_radius alias status', r_alias%status, r_legacy%status, failures)
    call check_int('initial_tr_radius alias nit', r_alias%nit, r_legacy%nit, failures)
    call check_int('initial_tr_radius alias nfev', r_alias%nfev, r_legacy%nfev, failures)

    r = minimize(linear_objective, [0.5_dp, 1.5_dp], method='trust-constr', &
        jac=linear_gradient, hess=linear_hessian, constraints=linear_constraint, &
        constraint_jac=linear_constraint_jacobian, meq=1, &
        initial_tr_radius=0.2_dp, gtol=1.0e-12_dp, xtol=1.0e-12_dp, maxiter=3)
    call check_result('initial_tr_radius minimize alias', r, r_legacy%x, r_legacy%fun, &
        r_legacy%maxcv, 5.0e-14_dp, 5.0e-14_dp, 5.0e-14_dp, failures)
    call check_int('initial_tr_radius minimize nit', r%nit, r_legacy%nit, failures)

    r = trust_constr(linear_objective, [0.5_dp, 1.5_dp], &
        constraints=linear_constraint, meq=1, initial_trust_radius=0.5_dp, &
        initial_tr_radius=0.25_dp)
    call check_int('conflicting trust radius spellings', r%status, -1, failures)
    if (r%success) then
        write(*, '(a)') 'conflicting trust radius spellings unexpectedly succeeded'
        failures = failures + 1
    end if

    r = trust_constr(linear_objective, [0.5_dp, 1.5_dp], &
        constraints=linear_constraint, meq=1, initial_tr_radius=-1.0_dp)
    call check_int('invalid initial_tr_radius status', r%status, -1, failures)
    if (r%success) then
        write(*, '(a)') 'invalid initial_tr_radius unexpectedly succeeded'
        failures = failures + 1
    end if

    r = minimize(linear_objective, [0.5_dp, 1.5_dp], method='trust-constr', &
        jac=linear_gradient, hess=linear_hessian, constraints=linear_constraint, &
        constraint_jac=linear_constraint_jacobian, meq=1, &
        initial_constr_penalty=7.5_dp, gtol=1.0e-10_dp, xtol=1.0e-10_dp)
    call check_real('initial penalty minimize dispatch', r%constr_penalty, &
        initial_penalty_ref, 5.0e-14_dp, failures)

    r = trust_constr(linear_objective, [2.0_dp, -1.0_dp], &
        constraints=linear_constraint, meq=1, gtol=1.0e-8_dp, xtol=1.0e-8_dp)
    call check_result('finite differences and BFGS', r, linear_x_ref, linear_fun_ref, &
        linear_maxcv_ref, 2.0e-7_dp, 2.0e-12_dp, 2.0e-10_dp, failures)

    r = trust_constr(circle_objective, [1.0_dp, 0.2_dp], jac=circle_gradient, &
        constraints=circle_constraint, constraint_jac=circle_constraint_jacobian, &
        meq=1, gtol=1.0e-9_dp, xtol=1.0e-10_dp, maxiter=200)
    call check_result('nonlinear equality BFGS', r, circle_x_ref, circle_fun_ref, &
        circle_maxcv_ref, 5.0e-8_dp, 5.0e-10_dp, 5.0e-9_dp, failures)

    constraint_hess_calls = 0
    r = trust_constr(circle_objective, [1.0_dp, 0.2_dp], jac=circle_gradient, &
        hess=circle_hessian, constraints=circle_constraint, &
        constraint_jac=circle_constraint_jacobian, constraint_hess=circle_constraint_hessian, &
        meq=1, gtol=1.0e-10_dp, xtol=1.0e-10_dp, maxiter=200)
    call check_result('analytic constraint Hessian', r, circle_hess_x_ref, &
        circle_hess_fun_ref, circle_hess_maxcv_ref, 5.0e-9_dp, 5.0e-10_dp, &
        5.0e-10_dp, failures)
    call check_int('analytic constraint Hessian status', r%status, &
        circle_hess_status_ref, failures)
    if (constraint_hess_calls <= 0) then
        write(*, '(a)') 'analytic constraint Hessian callback was not used'
        failures = failures + 1
    end if

    r = trust_constr(maratos_objective, [0.0_dp, 1.0_dp], jac=maratos_gradient, &
        hess=maratos_hessian, constraints=maratos_constraint, &
        constraint_jac=maratos_constraint_jacobian, &
        constraint_hess=maratos_constraint_hessian, meq=1, &
        gtol=1.0e-10_dp, xtol=1.0e-10_dp, maxiter=200)
    call check_result('Maratos SOC', r, maratos_x_ref, maratos_fun_ref, &
        maratos_maxcv_ref, 5.0e-9_dp, 5.0e-12_dp, 5.0e-12_dp, failures)
    call check_int('Maratos SOC objective evaluations', r%nfev, maratos_nfev_ref, failures)
    call check_vector('Maratos SOC Lagrangian gradient', r%lagrangian_grad, &
        maratos_lagrangian_grad_ref, 5.0e-12_dp, failures)
    call check_int('Maratos SOC CG iterations', r%cg_niter, maratos_cg_niter_ref, failures)
    call check_int('Maratos SOC CG stop condition', r%cg_stop_cond, &
        maratos_cg_stop_cond_ref, failures)
    if (trim(r%method) /= trim(maratos_method_ref)) then
        write(*, '(a,2a)') 'Maratos SOC method mismatch: ', trim(r%method), &
            ' /= '//trim(maratos_method_ref)
        failures = failures + 1
    end if
    call check_vector('Maratos SOC final constraints', r%constr, maratos_constr_ref, &
        5.0e-12_dp, failures)
    call check_matrix('Maratos SOC final constraint Jacobian', r%constr_jac, &
        maratos_constr_jac_ref, 5.0e-12_dp, failures)
    call check_int_vector('Maratos SOC constraint function evaluations', r%constr_nfev, &
        maratos_constr_nfev_ref, failures)
    call check_int_vector('Maratos SOC constraint Jacobian evaluations', r%constr_njev, &
        maratos_constr_njev_ref, failures)
    call check_int_vector('Maratos SOC constraint Hessian evaluations', r%constr_nhev, &
        maratos_constr_nhev_ref, failures)

    r = trust_constr(maratos_objective, [0.0_dp, 1.0_dp], jac=maratos_gradient, &
        hess=maratos_hessian, constraints=maratos_constraint, &
        constraint_jac=maratos_constraint_jacobian, &
        constraint_hess=maratos_constraint_hessian, mineq=1, &
        gtol=1.0e-10_dp, xtol=1.0e-10_dp, maxiter=500)
    call check_result('barrier Maratos SOC', r, barrier_maratos_x_ref, &
        barrier_maratos_fun_ref, barrier_maratos_maxcv_ref, 5.0e-8_dp, &
        5.0e-9_dp, 5.0e-12_dp, failures)
    call check_int('barrier Maratos SOC status', r%status, barrier_maratos_status_ref, failures)
    call check_int('barrier Maratos SOC iterations', r%nit, barrier_maratos_nit_ref, failures)
    call check_int('barrier Maratos SOC objective evaluations', r%nfev, &
        barrier_maratos_nfev_ref, failures)
    call check_int('barrier Maratos SOC gradient evaluations', r%njev, &
        barrier_maratos_njev_ref, failures)
    call check_int('barrier Maratos SOC Hessian evaluations', r%nhev, &
        barrier_maratos_nhev_ref, failures)
    call check_real('barrier Maratos SOC constraint violation', r%constr_violation, &
        barrier_maratos_maxcv_ref, 5.0e-12_dp, failures)
    call check_real('barrier Maratos SOC optimality', r%optimality, &
        barrier_maratos_optimality_ref, 5.0e-9_dp, failures)
    call check_vector('barrier Maratos SOC Lagrangian gradient', r%lagrangian_grad, &
        barrier_maratos_lagrangian_grad_ref, 5.0e-9_dp, failures)
    call check_int('barrier Maratos SOC CG iterations', r%cg_niter, &
        barrier_maratos_cg_niter_ref, failures)
    call check_int('barrier Maratos SOC CG stop condition', r%cg_stop_cond, &
        barrier_maratos_cg_stop_cond_ref, failures)
    call check_real('barrier Maratos SOC barrier parameter', r%barrier_parameter, &
        barrier_maratos_barrier_parameter_ref, 1.0e-15_dp, failures)
    call check_real('barrier Maratos SOC barrier tolerance', r%barrier_tolerance, &
        barrier_maratos_barrier_tolerance_ref, 1.0e-15_dp, failures)
    if (trim(r%method) /= trim(barrier_maratos_method_ref)) then
        write(*, '(a,2a)') 'barrier Maratos SOC method mismatch: ', trim(r%method), &
            ' /= '//trim(barrier_maratos_method_ref)
        failures = failures + 1
    end if
    call check_vector('barrier Maratos SOC final constraints', r%constr, &
        barrier_maratos_constr_ref, 5.0e-12_dp, failures)
    call check_matrix('barrier Maratos SOC final constraint Jacobian', r%constr_jac, &
        barrier_maratos_constr_jac_ref, 5.0e-12_dp, failures)
    call check_int_vector('barrier Maratos SOC constraint function evaluations', &
        r%constr_nfev, barrier_maratos_constr_nfev_ref, failures)
    call check_int_vector('barrier Maratos SOC constraint Jacobian evaluations', &
        r%constr_njev, barrier_maratos_constr_njev_ref, failures)
    call check_int_vector('barrier Maratos SOC constraint Hessian evaluations', &
        r%constr_nhev, barrier_maratos_constr_nhev_ref, failures)

    constraint_hess_calls = 0
    r = minimize(circle_objective, [1.0_dp, 0.2_dp], method='trust-constr', &
        jac=circle_gradient, hess=circle_hessian, constraints=circle_constraint, &
        constraint_jac=circle_constraint_jacobian, constraint_hess=circle_constraint_hessian, &
        meq=1, gtol=1.0e-10_dp, xtol=1.0e-10_dp, maxiter=200)
    call check_result('constraint Hessian minimize dispatch', r, circle_hess_x_ref, &
        circle_hess_fun_ref, circle_hess_maxcv_ref, 5.0e-9_dp, 5.0e-10_dp, &
        5.0e-10_dp, failures)
    if (constraint_hess_calls <= 0) then
        write(*, '(a)') 'minimize did not forward the constraint Hessian callback'
        failures = failures + 1
    end if

    r = trust_constr(circle_objective, [1.0_dp, 0.2_dp], jac=circle_gradient, &
        constraints=circle_constraint, constraint_jac=circle_constraint_jacobian, &
        constraint_hess=circle_constraint_hessian, meq=1, gtol=1.0e-9_dp)
    call check_result('constraint Hessian with BFGS objective', r, circle_hess_x_ref, &
        circle_hess_fun_ref, circle_hess_maxcv_ref, 5.0e-6_dp, 1.0e-8_dp, &
        1.0e-8_dp, failures)
    if (.not. allocated(r%constr_nhev)) then
        write(*, '(a)') 'missing BFGS/object constraint Hessian count'
        failures = failures + 1
    else if (r%constr_nhev(1) < 1) then
        write(*, '(a)') 'BFGS/object constraint Hessian was not evaluated'
        failures = failures + 1
    end if

    r = trust_constr(linear_objective, [0.0_dp, 0.0_dp], jac=linear_gradient, &
        hess=linear_hessian, constraints=active_inequality, &
        constraint_jac=active_inequality_jacobian, mineq=1, &
        gtol=1.0e-8_dp, xtol=1.0e-8_dp, maxiter=300)
    call check_result('active inequality', r, ineq_x_ref, ineq_fun_ref, &
        ineq_maxcv_ref, 5.0e-6_dp, 2.0e-6_dp, 1.0e-10_dp, failures)
    call check_int('active inequality status', r%status, ineq_status_ref, failures)

    r = trust_constr(linear_objective, [0.0_dp, 0.0_dp], jac=linear_gradient, &
        hess=linear_hessian, constraints=active_inequality, &
        constraint_jac=active_inequality_jacobian, mineq=1, &
        initial_constr_penalty=7.5_dp, gtol=1.0e-8_dp, xtol=1.0e-8_dp, maxiter=300)
    call check_result('initial penalty interior point', r, ineq_x_ref, ineq_fun_ref, &
        ineq_maxcv_ref, 5.0e-6_dp, 2.0e-6_dp, 1.0e-10_dp, failures)
    call check_real('interior-point constraint penalty result', r%constr_penalty, &
        initial_penalty_ref, 5.0e-12_dp, failures)

    r = trust_constr(linear_objective, [0.0_dp, 0.0_dp], jac=linear_gradient, &
        hess=linear_hessian, constraints=active_inequality, &
        constraint_jac=active_inequality_jacobian, mineq=1, &
        gtol=1.0e-8_dp, xtol=1.0e-8_dp, barrier_tol=1.0e-6_dp, &
        initial_barrier_parameter=5.0e-2_dp, &
        initial_barrier_tolerance=2.5e-2_dp, maxiter=300)
    call check_result('custom barrier controls', r, barrier_options_x_ref, &
        barrier_options_fun_ref, barrier_options_maxcv_ref, 5.0e-6_dp, &
        3.0e-5_dp, 1.0e-10_dp, failures)
    call check_int('custom barrier controls status', r%status, &
        barrier_options_status_ref, failures)

    r = minimize(linear_objective, [0.0_dp, 0.0_dp], method='trust-constr', &
        jac=linear_gradient, hess=linear_hessian, constraints=active_inequality, &
        constraint_jac=active_inequality_jacobian, mineq=1, &
        gtol=1.0e-8_dp, xtol=1.0e-8_dp, barrier_tol=1.0e-6_dp, &
        initial_barrier_parameter=5.0e-2_dp, &
        initial_barrier_tolerance=2.5e-2_dp, maxiter=300)
    call check_result('custom barrier minimize dispatch', r, barrier_options_x_ref, &
        barrier_options_fun_ref, barrier_options_maxcv_ref, 5.0e-6_dp, &
        3.0e-5_dp, 1.0e-10_dp, failures)

    r = trust_constr(linear_objective, [0.0_dp, 0.0_dp], &
        constraints=active_inequality, mineq=1, gtol=1.0e-6_dp, &
        xtol=1.0e-8_dp, maxiter=300)
    call check_result('active inequality finite differences', r, ineq_x_ref, &
        ineq_fun_ref, ineq_maxcv_ref, 8.0e-6_dp, 8.0e-6_dp, 1.0e-8_dp, failures)

    r = trust_constr(circle_objective, [0.2_dp, 0.1_dp], jac=circle_gradient, &
        constraints=circle_inequality, constraint_jac=circle_inequality_jacobian, &
        mineq=1, gtol=1.0e-8_dp, xtol=1.0e-8_dp, maxiter=400)
    call check_result('nonlinear inequality BFGS', r, circle_ineq_x_ref, &
        circle_ineq_fun_ref, circle_ineq_maxcv_ref, 5.0e-4_dp, 1.0e-3_dp, &
        1.0e-8_dp, failures)

    r = minimize(linear_objective, [0.0_dp, 0.0_dp], method='trust-constr', &
        jac=linear_gradient, hess=linear_hessian, constraints=active_inequality, &
        constraint_jac=active_inequality_jacobian, mineq=1, &
        gtol=1.0e-8_dp, xtol=1.0e-8_dp, maxiter=300)
    call check_result('inequality minimize dispatch', r, ineq_x_ref, ineq_fun_ref, &
        ineq_maxcv_ref, 5.0e-6_dp, 2.0e-6_dp, 1.0e-10_dp, failures)

    r = trust_constr(bound_objective, [0.5_dp, 1.0_dp], jac=bound_gradient, &
        hess=linear_hessian, &
        lower=[0.0_dp, ieee_value(0.0_dp, ieee_negative_inf)], &
        upper=[ieee_value(0.0_dp, ieee_positive_inf), 2.0_dp], &
        gtol=1.0e-8_dp, xtol=1.0e-8_dp, maxiter=300)
    call check_result('one-sided bounds', r, bounds_x_ref, bounds_fun_ref, &
        bounds_maxcv_ref, 2.0e-4_dp, 4.0e-4_dp, 1.0e-10_dp, failures)
    if (trim(r%method) /= trim(bounds_method_ref)) then
        write(*, '(a,2a)') 'one-sided bounds method mismatch: ', trim(r%method), &
            ' /= '//trim(bounds_method_ref)
        failures = failures + 1
    end if
    call check_vector('one-sided bounds final constraints', r%constr, bounds_constr_ref, &
        2.0e-4_dp, failures)
    call check_matrix('one-sided bounds final constraint Jacobian', r%constr_jac, &
        bounds_constr_jac_ref, 5.0e-14_dp, failures)
    call check_int_vector('one-sided bounds constraint function evaluations', r%constr_nfev, &
        bounds_constr_nfev_ref, failures)
    call check_int_vector('one-sided bounds constraint Jacobian evaluations', r%constr_njev, &
        bounds_constr_njev_ref, failures)
    call check_int_vector('one-sided bounds constraint Hessian evaluations', r%constr_nhev, &
        bounds_constr_nhev_ref, failures)

    r = trust_constr(linear_objective, [0.8_dp, 1.2_dp], jac=linear_gradient, &
        hess=linear_hessian, constraints=mixed_constraint, &
        constraint_jac=mixed_constraint_jacobian, meq=1, mineq=1, &
        gtol=1.0e-8_dp, xtol=1.0e-8_dp, maxiter=300)
    call check_result('mixed equality inequality', r, mixed_x_ref, mixed_fun_ref, &
        mixed_maxcv_ref, 2.0e-4_dp, 2.0e-4_dp, 1.0e-10_dp, failures)

    r = trust_constr(linear_objective, [0.8_dp, 1.2_dp], jac=linear_gradient, &
        hess=linear_hessian, constraints=mixed_constraint, &
        constraint_jac=mixed_constraint_jacobian, meq=1, mineq=1, &
        lower=[0.0_dp, ieee_value(0.0_dp, ieee_negative_inf)], &
        upper=[ieee_value(0.0_dp, ieee_positive_inf), 2.0_dp], &
        gtol=1.0e-8_dp, xtol=1.0e-8_dp, maxiter=300)
    call check_result('mixed constraints and bounds diagnostics', r, mixed_bounds_x_ref, &
        mixed_bounds_fun_ref, mixed_bounds_maxcv_ref, 5.0e-10_dp, 5.0e-12_dp, &
        1.0e-10_dp, failures)
    call check_int('mixed constraints and bounds iterations', r%nit, mixed_bounds_nit_ref, failures)
    if (trim(r%method) /= trim(mixed_bounds_method_ref)) then
        write(*, '(a,2a)') 'mixed constraints and bounds method mismatch: ', trim(r%method), &
            ' /= '//trim(mixed_bounds_method_ref)
        failures = failures + 1
    end if
    call check_vector('mixed constraints and bounds final constraints', r%constr, &
        mixed_bounds_constr_ref, 5.0e-10_dp, failures)
    call check_matrix('mixed constraints and bounds final constraint Jacobian', r%constr_jac, &
        mixed_bounds_constr_jac_ref, 5.0e-14_dp, failures)
    call check_int_vector('mixed constraints and bounds function evaluations', r%constr_nfev, &
        mixed_bounds_constr_nfev_ref, failures)
    call check_int_vector('mixed constraints and bounds Jacobian evaluations', r%constr_njev, &
        mixed_bounds_constr_njev_ref, failures)
    call check_int_vector('mixed constraints and bounds Hessian evaluations', r%constr_nhev, &
        mixed_bounds_constr_nhev_ref, failures)

    r = trust_constr(linear_objective, [0.8_dp, 1.2_dp], jac=linear_gradient, &
        hess=linear_hessian, constraints=bounded_constraint, &
        constraint_jac=bounded_constraint_jacobian, &
        constraint_lower=[2.0_dp, 0.5_dp, ieee_value(0.0_dp, ieee_negative_inf), -1.2_dp], &
        constraint_upper=[2.0_dp, ieee_value(0.0_dp, ieee_positive_inf), 1.6_dp, -0.6_dp], &
        gtol=1.0e-8_dp, xtol=1.0e-8_dp, maxiter=300)
    call check_result('bounded nonlinear components', r, bounded_nonlinear_x_ref, &
        bounded_nonlinear_fun_ref, bounded_nonlinear_maxcv_ref, 5.0e-10_dp, &
        5.0e-11_dp, 1.0e-12_dp, failures)
    call check_int('bounded nonlinear status', r%status, bounded_nonlinear_status_ref, failures)
    call check_int('bounded nonlinear iterations', r%nit, bounded_nonlinear_nit_ref, failures)
    call check_int('bounded nonlinear objective evaluations', r%nfev, &
        bounded_nonlinear_nfev_ref, failures)
    call check_int('bounded nonlinear CG iterations', r%cg_niter, &
        bounded_nonlinear_cg_niter_ref, failures)
    if (trim(r%method) /= trim(bounded_nonlinear_method_ref)) then
        write(*, '(a,2a)') 'bounded nonlinear method mismatch: ', trim(r%method), &
            ' /= '//trim(bounded_nonlinear_method_ref)
        failures = failures + 1
    end if
    call check_vector('bounded nonlinear final constraints', r%constr, &
        bounded_nonlinear_constr_ref, 5.0e-10_dp, failures)
    call check_matrix('bounded nonlinear final constraint Jacobian', r%constr_jac, &
        bounded_nonlinear_constr_jac_ref, 5.0e-14_dp, failures)
    call check_int_vector('bounded nonlinear constraint function evaluations', &
        r%constr_nfev, bounded_nonlinear_constr_nfev_ref, failures)
    call check_int_vector('bounded nonlinear constraint Jacobian evaluations', &
        r%constr_njev, bounded_nonlinear_constr_njev_ref, failures)
    call check_int_vector('bounded nonlinear constraint Hessian evaluations', &
        r%constr_nhev, bounded_nonlinear_constr_nhev_ref, failures)

    constraint_hess_calls = 0
    r = trust_constr(circle_objective, [0.2_dp, 0.1_dp], jac=circle_gradient, &
        hess=circle_hessian, constraints=circle_value_constraint, &
        constraint_jac=circle_value_constraint_jacobian, &
        constraint_hess=circle_value_constraint_hessian, &
        constraint_upper=[1.0_dp], gtol=1.0e-8_dp, xtol=1.0e-8_dp, maxiter=400)
    call check_result('upper-bounded nonlinear constraint Hessian', r, &
        circle_ineq_hess_x_ref, circle_ineq_hess_fun_ref, circle_ineq_hess_maxcv_ref, &
        5.0e-4_dp, 1.0e-3_dp, 1.0e-8_dp, failures)
    if (constraint_hess_calls <= 0) then
        write(*, '(a)') 'upper-bounded constraint Hessian callback was not used'
        failures = failures + 1
    end if

    r = minimize(linear_objective, [0.8_dp, 1.2_dp], method='trust-constr', &
        jac=linear_gradient, hess=linear_hessian, constraints=bounded_constraint, &
        constraint_jac=bounded_constraint_jacobian, &
        constraint_lower=[2.0_dp, 0.5_dp, ieee_value(0.0_dp, ieee_negative_inf), -1.2_dp], &
        constraint_upper=[2.0_dp, ieee_value(0.0_dp, ieee_positive_inf), 1.6_dp, -0.6_dp], &
        gtol=1.0e-8_dp, xtol=1.0e-8_dp, maxiter=300)
    call check_result('bounded nonlinear minimize dispatch', r, bounded_nonlinear_x_ref, &
        bounded_nonlinear_fun_ref, bounded_nonlinear_maxcv_ref, 5.0e-10_dp, &
        5.0e-11_dp, 1.0e-12_dp, failures)

    r = trust_constr(keep_feasible_objective, [2.0_dp], jac=keep_feasible_gradient, &
        hess=keep_feasible_hessian, constraints=keep_feasible_constraint, &
        constraint_jac=keep_feasible_constraint_jacobian, constraint_lower=[0.5_dp], &
        constraint_upper=[ieee_value(0.0_dp, ieee_positive_inf)], &
        constraint_keep_feasible=[.true.], gtol=1.0e-8_dp, xtol=1.0e-8_dp, maxiter=100)
    call check_result('keep-feasible nonlinear constraint', r, keep_feasible_nonlinear_x_ref, &
        keep_feasible_nonlinear_fun_ref, keep_feasible_nonlinear_maxcv_ref, 5.0e-13_dp, &
        5.0e-13_dp, 5.0e-14_dp, failures)
    call check_int('keep-feasible nonlinear status', r%status, &
        keep_feasible_nonlinear_status_ref, failures)
    call check_int('keep-feasible nonlinear iterations', r%nit, &
        keep_feasible_nonlinear_nit_ref, failures)
    call check_int('keep-feasible nonlinear objective evaluations', r%nfev, &
        keep_feasible_nonlinear_nfev_ref, failures)
    call check_int('keep-feasible nonlinear gradient evaluations', r%njev, &
        keep_feasible_nonlinear_njev_ref, failures)
    call check_int('keep-feasible nonlinear Hessian evaluations', r%nhev, &
        keep_feasible_nonlinear_nhev_ref, failures)
    call check_int('keep-feasible nonlinear CG iterations', r%cg_niter, &
        keep_feasible_nonlinear_cg_niter_ref, failures)
    if (r%success .neqv. keep_feasible_nonlinear_success_ref) then
        write(*, '(a)') 'keep-feasible nonlinear success mismatch'
        failures = failures + 1
    end if
    if (trim(r%method) /= trim(keep_feasible_nonlinear_method_ref)) then
        write(*, '(a,2a)') 'keep-feasible nonlinear method mismatch: ', trim(r%method), &
            ' /= '//trim(keep_feasible_nonlinear_method_ref)
        failures = failures + 1
    end if
    call check_vector('keep-feasible nonlinear final constraints', r%constr, &
        keep_feasible_nonlinear_constr_ref, 5.0e-13_dp, failures)
    call check_matrix('keep-feasible nonlinear final constraint Jacobian', r%constr_jac, &
        keep_feasible_nonlinear_constr_jac_ref, 5.0e-14_dp, failures)
    call check_int_vector('keep-feasible nonlinear constraint function evaluations', &
        r%constr_nfev, keep_feasible_nonlinear_constr_nfev_ref, failures)
    call check_int_vector('keep-feasible nonlinear constraint Jacobian evaluations', &
        r%constr_njev, keep_feasible_nonlinear_constr_njev_ref, failures)
    call check_int_vector('keep-feasible nonlinear constraint Hessian evaluations', &
        r%constr_nhev, keep_feasible_nonlinear_constr_nhev_ref, failures)

    r = minimize(keep_feasible_objective, [0.6_dp], method='trust-constr', &
        jac=keep_feasible_gradient, hess=keep_feasible_hessian, lower=[0.5_dp], &
        upper=[ieee_value(0.0_dp, ieee_positive_inf)], bounds_keep_feasible=[.true.], &
        gtol=1.0e-8_dp, xtol=1.0e-8_dp, barrier_tol=1.0e-3_dp, maxiter=100)
    call check_result('keep-feasible variable bound', r, keep_feasible_bounds_x_ref, &
        keep_feasible_bounds_fun_ref, keep_feasible_bounds_maxcv_ref, 5.0e-13_dp, &
        5.0e-13_dp, 5.0e-14_dp, failures)
    call check_int('keep-feasible bound status', r%status, keep_feasible_bounds_status_ref, failures)
    call check_int('keep-feasible bound iterations', r%nit, keep_feasible_bounds_nit_ref, failures)
    call check_int('keep-feasible bound objective evaluations', r%nfev, &
        keep_feasible_bounds_nfev_ref, failures)
    call check_int('keep-feasible bound gradient evaluations', r%njev, &
        keep_feasible_bounds_njev_ref, failures)
    call check_int('keep-feasible bound Hessian evaluations', r%nhev, &
        keep_feasible_bounds_nhev_ref, failures)
    call check_int('keep-feasible bound CG iterations', r%cg_niter, &
        keep_feasible_bounds_cg_niter_ref, failures)
    if (r%success .neqv. keep_feasible_bounds_success_ref) then
        write(*, '(a)') 'keep-feasible bound success mismatch'
        failures = failures + 1
    end if
    if (trim(r%method) /= trim(keep_feasible_bounds_method_ref)) then
        write(*, '(a,2a)') 'keep-feasible bound method mismatch: ', trim(r%method), &
            ' /= '//trim(keep_feasible_bounds_method_ref)
        failures = failures + 1
    end if
    call check_vector('keep-feasible bound final constraints', r%constr, &
        keep_feasible_bounds_constr_ref, 5.0e-13_dp, failures)
    call check_matrix('keep-feasible bound final constraint Jacobian', r%constr_jac, &
        keep_feasible_bounds_constr_jac_ref, 5.0e-14_dp, failures)
    call check_int_vector('keep-feasible bound constraint function evaluations', &
        r%constr_nfev, keep_feasible_bounds_constr_nfev_ref, failures)
    call check_int_vector('keep-feasible bound constraint Jacobian evaluations', &
        r%constr_njev, keep_feasible_bounds_constr_njev_ref, failures)
    call check_int_vector('keep-feasible bound constraint Hessian evaluations', &
        r%constr_nhev, keep_feasible_bounds_constr_nhev_ref, failures)


    fd_clip_objective_max = -huge(1.0_dp)
    fd_clip_constraint_max = -huge(1.0_dp)
    r = trust_constr(fd_clip_objective, [0.99999999_dp], &
        constraints=fd_clip_constraint, constraint_lower=[0.2_dp], &
        constraint_upper=[ieee_value(0.0_dp, ieee_positive_inf)], &
        lower=[0.0_dp], upper=[1.0_dp], bounds_keep_feasible=[.true.], &
        gtol=1.0e-8_dp, xtol=1.0e-8_dp, maxiter=100)
    call check_result('keep-feasible finite-difference clipping', r, fd_clip_x_ref, &
        fd_clip_fun_ref, fd_clip_maxcv_ref, 2.0e-8_dp, 5.0e-8_dp, 1.0e-12_dp, failures)
    call check_int('finite-difference clipping status', r%status, fd_clip_status_ref, failures)
    if (fd_clip_objective_max > 1.0_dp) then
        write(*, '(a,es16.8)') 'objective finite-difference probe crossed kept upper bound: ', &
            fd_clip_objective_max
        failures = failures + 1
    end if
    if (fd_clip_constraint_max > 1.0_dp) then
        write(*, '(a,es16.8)') 'constraint finite-difference probe crossed kept upper bound: ', &
            fd_clip_constraint_max
        failures = failures + 1
    end if

    r = trust_constr(keep_feasible_objective, [0.4_dp], jac=keep_feasible_gradient, &
        hess=keep_feasible_hessian, constraints=keep_feasible_constraint, &
        constraint_jac=keep_feasible_constraint_jacobian, constraint_lower=[0.5_dp], &
        constraint_upper=[ieee_value(0.0_dp, ieee_positive_inf)], &
        constraint_keep_feasible=[.true.])
    call check_int('keep-feasible nonlinear infeasible start', r%status, -1, failures)
    if (r%success) then
        write(*, '(a)') 'keep-feasible nonlinear infeasible start unexpectedly succeeded'
        failures = failures + 1
    end if

    r = trust_constr(keep_feasible_objective, [0.4_dp], jac=keep_feasible_gradient, &
        hess=keep_feasible_hessian, lower=[0.5_dp], &
        upper=[ieee_value(0.0_dp, ieee_positive_inf)], bounds_keep_feasible=[.true.])
    call check_int('keep-feasible bound infeasible start', r%status, -1, failures)
    if (r%success) then
        write(*, '(a)') 'keep-feasible bound infeasible start unexpectedly succeeded'
        failures = failures + 1
    end if

    r = trust_constr(linear_objective, [0.8_dp, 1.2_dp], constraints=bounded_constraint, &
        meq=1, constraint_lower=[2.0_dp, 0.5_dp, -1.0_dp, -1.2_dp], &
        constraint_upper=[2.0_dp, 2.0_dp, 1.6_dp, -0.6_dp])
    call check_int('bounded nonlinear legacy-count conflict', r%status, -1, failures)

    r = trust_constr(fixed_bound_objective, [0.7_dp, 1.0_dp], &
        jac=fixed_bound_gradient, hess=linear_hessian, &
        lower=[0.5_dp, ieee_value(0.0_dp, ieee_negative_inf)], &
        upper=[0.5_dp, 2.0_dp], gtol=1.0e-8_dp, xtol=1.0e-8_dp, maxiter=300)
    call check_result('fixed and upper bounds', r, fixed_bound_x_ref, &
        fixed_bound_fun_ref, fixed_bound_maxcv_ref, 2.0e-6_dp, 3.0e-6_dp, &
        2.0e-8_dp, failures)

    r = trust_constr(linear_objective, [0.0_dp, 0.0_dp], &
        lower=[1.0_dp, 0.0_dp], upper=[0.0_dp, 1.0_dp])
    call check_int('invalid bounds status', r%status, -1, failures)
    if (r%success) then
        write(*, '(a)') 'invalid bounds unexpectedly succeeded'
        failures = failures + 1
    end if

    constraint_hess_calls = 0
    r = trust_constr(circle_objective, [0.2_dp, 0.1_dp], jac=circle_gradient, &
        hess=circle_hessian, constraints=circle_inequality, &
        constraint_jac=circle_inequality_jacobian, &
        constraint_hess=circle_inequality_hessian, mineq=1, &
        gtol=1.0e-8_dp, xtol=1.0e-8_dp, maxiter=400)
    call check_result('nonlinear inequality analytic constraint Hessian', r, &
        circle_ineq_hess_x_ref, circle_ineq_hess_fun_ref, circle_ineq_hess_maxcv_ref, &
        5.0e-4_dp, 1.0e-3_dp, 1.0e-8_dp, failures)
    call check_int('nonlinear inequality analytic constraint Hessian status', r%status, &
        circle_ineq_hess_status_ref, failures)
    if (constraint_hess_calls <= 0) then
        write(*, '(a)') 'inequality constraint Hessian callback was not used'
        failures = failures + 1
    end if

    constraint_hess_calls = 0
    r = minimize(circle_objective, [0.2_dp, 0.1_dp], method='trust-constr', &
        jac=circle_gradient, hess=circle_hessian, constraints=circle_inequality, &
        constraint_jac=circle_inequality_jacobian, &
        constraint_hess=circle_inequality_hessian, mineq=1, &
        gtol=1.0e-8_dp, xtol=1.0e-8_dp, maxiter=400)
    call check_result('inequality constraint Hessian minimize dispatch', r, &
        circle_ineq_hess_x_ref, circle_ineq_hess_fun_ref, circle_ineq_hess_maxcv_ref, &
        5.0e-4_dp, 1.0e-3_dp, 1.0e-8_dp, failures)
    if (constraint_hess_calls <= 0) then
        write(*, '(a)') 'minimize did not forward the inequality constraint Hessian callback'
        failures = failures + 1
    end if

    r = trust_constr(linear_objective, [0.0_dp, 0.0_dp], &
        constraints=active_inequality, mineq=1, initial_barrier_parameter=0.0_dp)
    call check_int('invalid initial barrier status', r%status, -1, failures)
    if (r%success) then
        write(*, '(a)') 'invalid initial barrier unexpectedly succeeded'
        failures = failures + 1
    end if

    if (failures /= 0) then
        write(*, '(a,i0)') 'test_optimize_trust_constr: FAIL ', failures
        stop 1
    end if
    write(*, '(a)') 'test_optimize_trust_constr: PASS'

contains

    function linear_objective(x) result(fx)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp) :: fx

        fx = (x(1) - 1.0_dp)**2 + (x(2) - 2.5_dp)**2
    end function linear_objective

    subroutine linear_gradient(x, grad)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: grad(:) !! objective gradient

        grad = [2.0_dp * (x(1) - 1.0_dp), 2.0_dp * (x(2) - 2.5_dp)]
    end subroutine linear_gradient

    subroutine linear_hessian(x, hmat)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: hmat(:, :) !! objective Hessian

        if (size(x) < 0) hmat = 0.0_dp
        hmat = 0.0_dp
        hmat(1, 1) = 2.0_dp
        hmat(2, 2) = 2.0_dp
    end subroutine linear_hessian

    subroutine linear_constraint(x, c)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: c(:) !! equality residual

        c(1) = x(1) + x(2) - 2.0_dp
    end subroutine linear_constraint

    subroutine linear_constraint_jacobian(x, amat)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: amat(:, :) !! constraint Jacobian

        if (size(x) < 0) amat = 0.0_dp
        amat(1, :) = [1.0_dp, 1.0_dp]
    end subroutine linear_constraint_jacobian

    function circle_objective(x) result(fx)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp) :: fx

        fx = (x(1) - 2.0_dp)**2 + (x(2) - 1.0_dp)**2
    end function circle_objective

    subroutine circle_gradient(x, grad)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: grad(:) !! objective gradient

        grad = [2.0_dp * (x(1) - 2.0_dp), 2.0_dp * (x(2) - 1.0_dp)]
    end subroutine circle_gradient

    subroutine circle_constraint(x, c)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: c(:) !! equality residual

        c(1) = x(1)**2 + x(2)**2 - 1.0_dp
    end subroutine circle_constraint

    subroutine circle_constraint_jacobian(x, amat)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: amat(:, :) !! constraint Jacobian

        amat(1, :) = [2.0_dp * x(1), 2.0_dp * x(2)]
    end subroutine circle_constraint_jacobian

    subroutine circle_hessian(x, hmat)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: hmat(:, :) !! objective Hessian

        if (size(x) < 0) hmat = 0.0_dp
        hmat = 0.0_dp
        hmat(1, 1) = 2.0_dp
        hmat(2, 2) = 2.0_dp
    end subroutine circle_hessian

    subroutine circle_constraint_hessian(x, v, hmat)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(in) :: v(:) !! equality multiplier
        real(dp), intent(out) :: hmat(:, :) !! Hessian of v(1) * circle_constraint(x)

        if (size(x) < 0) hmat = 0.0_dp
        constraint_hess_calls = constraint_hess_calls + 1
        hmat = 0.0_dp
        hmat(1, 1) = 2.0_dp * v(1)
        hmat(2, 2) = 2.0_dp * v(1)
    end subroutine circle_constraint_hessian

    function maratos_objective(x) result(fx)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp) :: fx

        fx = 2.0_dp * (x(1)**2 + x(2)**2 - 1.0_dp) - x(1)
    end function maratos_objective

    subroutine maratos_gradient(x, grad)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: grad(:) !! objective gradient

        grad = [4.0_dp * x(1) - 1.0_dp, 4.0_dp * x(2)]
    end subroutine maratos_gradient

    subroutine maratos_hessian(x, hmat)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: hmat(:, :) !! objective Hessian

        if (size(x) < 0) hmat = 0.0_dp
        hmat = 0.0_dp
        hmat(1, 1) = 4.0_dp
        hmat(2, 2) = 4.0_dp
    end subroutine maratos_hessian

    subroutine maratos_constraint(x, c)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: c(:) !! equality residual

        c(1) = x(1)**2 + x(2)**2 - 1.0_dp
    end subroutine maratos_constraint

    subroutine maratos_constraint_jacobian(x, amat)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: amat(:, :) !! equality Jacobian

        amat(1, :) = [2.0_dp * x(1), 2.0_dp * x(2)]
    end subroutine maratos_constraint_jacobian

    subroutine maratos_constraint_hessian(x, v, hmat)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(in) :: v(:) !! equality multiplier
        real(dp), intent(out) :: hmat(:, :) !! Hessian of v(1) * constraint(x)

        if (size(x) < 0) hmat = 0.0_dp
        hmat = 0.0_dp
        hmat(1, 1) = 2.0_dp * v(1)
        hmat(2, 2) = 2.0_dp * v(1)
    end subroutine maratos_constraint_hessian

    subroutine active_inequality(x, c)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: c(:) !! inequality residual, feasible when nonnegative

        c(1) = 1.0_dp - x(1) - x(2)
    end subroutine active_inequality

    subroutine active_inequality_jacobian(x, amat)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: amat(:, :) !! inequality Jacobian

        if (size(x) < 0) amat = 0.0_dp
        amat(1, :) = [-1.0_dp, -1.0_dp]
    end subroutine active_inequality_jacobian

    subroutine circle_inequality(x, c)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: c(:) !! feasible inside the unit circle

        c(1) = 1.0_dp - x(1)**2 - x(2)**2
    end subroutine circle_inequality

    subroutine circle_inequality_jacobian(x, amat)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: amat(:, :) !! inequality Jacobian

        amat(1, :) = [-2.0_dp * x(1), -2.0_dp * x(2)]
    end subroutine circle_inequality_jacobian

    subroutine circle_inequality_hessian(x, v, hmat)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(in) :: v(:) !! inequality multiplier
        real(dp), intent(out) :: hmat(:, :) !! Hessian of v(1) * circle_inequality(x)

        if (size(x) < 0) hmat = 0.0_dp
        constraint_hess_calls = constraint_hess_calls + 1
        hmat = 0.0_dp
        hmat(1, 1) = -2.0_dp * v(1)
        hmat(2, 2) = -2.0_dp * v(1)
    end subroutine circle_inequality_hessian

    function bound_objective(x) result(fx)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp) :: fx

        fx = (x(1) + 1.0_dp)**2 + (x(2) - 3.0_dp)**2
    end function bound_objective

    subroutine bound_gradient(x, grad)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: grad(:) !! objective gradient

        grad = [2.0_dp * (x(1) + 1.0_dp), 2.0_dp * (x(2) - 3.0_dp)]
    end subroutine bound_gradient

    function keep_feasible_objective(x) result(fx)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp) :: fx

        fx = (x(1) + 2.0_dp)**2
    end function keep_feasible_objective

    subroutine keep_feasible_gradient(x, grad)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: grad(:) !! objective gradient

        grad(1) = 2.0_dp * (x(1) + 2.0_dp)
    end subroutine keep_feasible_gradient

    subroutine keep_feasible_hessian(x, hmat)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: hmat(:, :) !! objective Hessian

        if (size(x) < 0) hmat = 0.0_dp
        hmat = 0.0_dp
        hmat(1, 1) = 2.0_dp
    end subroutine keep_feasible_hessian

    subroutine keep_feasible_constraint(x, c)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: c(:) !! raw nonlinear constraint value

        c(1) = x(1)
    end subroutine keep_feasible_constraint

    subroutine keep_feasible_constraint_jacobian(x, amat)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: amat(:, :) !! raw nonlinear constraint Jacobian

        if (size(x) < 0) amat = 0.0_dp
        amat(1, 1) = 1.0_dp
    end subroutine keep_feasible_constraint_jacobian


    function fd_clip_objective(x) result(fx)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp) :: fx !! objective value

        fd_clip_objective_max = max(fd_clip_objective_max, x(1))
        fx = (x(1) - 2.0_dp)**2
    end function fd_clip_objective

    subroutine fd_clip_constraint(x, c)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: c(:) !! raw nonlinear constraint value

        fd_clip_constraint_max = max(fd_clip_constraint_max, x(1))
        c(1) = x(1)**2
    end subroutine fd_clip_constraint

    subroutine mixed_constraint(x, c)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: c(:) !! equality first, then inequality >= 0

        c = [x(1) + x(2) - 2.0_dp, x(1) - 0.5_dp]
    end subroutine mixed_constraint

    subroutine mixed_constraint_jacobian(x, amat)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: amat(:, :) !! combined constraint Jacobian

        if (size(x) < 0) amat = 0.0_dp
        amat(1, :) = [1.0_dp, 1.0_dp]
        amat(2, :) = [1.0_dp, 0.0_dp]
    end subroutine mixed_constraint_jacobian

    subroutine circle_value_constraint(x, c)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: c(:) !! raw circle value

        c(1) = x(1)**2 + x(2)**2
    end subroutine circle_value_constraint

    subroutine circle_value_constraint_jacobian(x, amat)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: amat(:, :) !! raw circle Jacobian

        amat(1, :) = [2.0_dp * x(1), 2.0_dp * x(2)]
    end subroutine circle_value_constraint_jacobian

    subroutine circle_value_constraint_hessian(x, v, hmat)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(in) :: v(:) !! raw constraint multiplier
        real(dp), intent(out) :: hmat(:, :) !! weighted raw circle Hessian

        if (size(x) < 0) hmat = 0.0_dp
        hmat = 0.0_dp
        hmat(1, 1) = 2.0_dp * v(1)
        hmat(2, 2) = 2.0_dp * v(1)
        constraint_hess_calls = constraint_hess_calls + 1
    end subroutine circle_value_constraint_hessian

    subroutine bounded_constraint(x, c)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: c(:) !! raw nonlinear constraint values

        c = [x(1) + x(2), x(1), x(2), x(1) - x(2)]
    end subroutine bounded_constraint

    subroutine bounded_constraint_jacobian(x, amat)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: amat(:, :) !! raw nonlinear constraint Jacobian

        if (size(x) < 0) amat = 0.0_dp
        amat(1, :) = [1.0_dp, 1.0_dp]
        amat(2, :) = [1.0_dp, 0.0_dp]
        amat(3, :) = [0.0_dp, 1.0_dp]
        amat(4, :) = [1.0_dp, -1.0_dp]
    end subroutine bounded_constraint_jacobian

    function fixed_bound_objective(x) result(fx)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp) :: fx

        fx = (x(1) - 1.0_dp)**2 + (x(2) - 3.0_dp)**2
    end function fixed_bound_objective

    subroutine fixed_bound_gradient(x, grad)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: grad(:) !! objective gradient

        grad = [2.0_dp * (x(1) - 1.0_dp), 2.0_dp * (x(2) - 3.0_dp)]
    end subroutine fixed_bound_gradient

    subroutine check_result(name, actual, expected_x, expected_fun, expected_maxcv, &
            xtol, ftol, cvtol, nfail)
        character(len=*), intent(in) :: name !! check label
        type(multivariate_minimize_result), intent(in) :: actual !! actual result
        real(dp), intent(in) :: expected_x(:) !! SciPy reference point
        real(dp), intent(in) :: expected_fun !! SciPy reference objective
        real(dp), intent(in) :: expected_maxcv !! SciPy reference constraint violation
        real(dp), intent(in) :: xtol !! point tolerance
        real(dp), intent(in) :: ftol !! objective tolerance
        real(dp), intent(in) :: cvtol !! constraint-violation tolerance
        integer, intent(inout) :: nfail !! accumulated failures

        if (size(actual%x) /= size(expected_x) .or. &
                maxval(abs(actual%x - expected_x)) > xtol) then
            write(*, '(a,a)') trim(name), ': x mismatch'
            nfail = nfail + 1
        end if
        if (abs(actual%fun - expected_fun) > ftol) then
            write(*, '(a,a,2es16.8)') trim(name), ': objective mismatch ', &
                actual%fun, expected_fun
            nfail = nfail + 1
        end if
        if (abs(actual%maxcv - expected_maxcv) > cvtol) then
            write(*, '(a,a,2es16.8)') trim(name), ': maxcv mismatch ', &
                actual%maxcv, expected_maxcv
            nfail = nfail + 1
        end if
        if (.not. actual%success) then
            write(*, '(a,1x,i0,2a)') trim(name)//' failed with status', actual%status, &
                ': ', trim(actual%message)
            nfail = nfail + 1
        end if
    end subroutine check_result

    subroutine check_vector(name, actual, expected, tol, nfail)
        character(len=*), intent(in) :: name !! check label
        real(dp), intent(in) :: actual(:) !! actual vector
        real(dp), intent(in) :: expected(:) !! expected vector
        real(dp), intent(in) :: tol !! maximum absolute error tolerance
        integer, intent(inout) :: nfail !! accumulated failures

        if (size(actual) /= size(expected)) then
            write(*, '(a,a,2i8)') trim(name), ': size mismatch ', size(actual), size(expected)
            nfail = nfail + 1
        else if (maxval(abs(actual - expected)) > tol) then
            write(*, '(a,a,es16.8)') trim(name), ': maximum error ', &
                maxval(abs(actual - expected))
            nfail = nfail + 1
        end if
    end subroutine check_vector

    subroutine check_matrix(name, actual, expected, tol, nfail)
        character(len=*), intent(in) :: name !! check label
        real(dp), intent(in) :: actual(:, :) !! actual matrix
        real(dp), intent(in) :: expected(:, :) !! expected matrix
        real(dp), intent(in) :: tol !! maximum absolute error tolerance
        integer, intent(inout) :: nfail !! accumulated failures

        if (any(shape(actual) /= shape(expected))) then
            write(*, '(a,a,4i8)') trim(name), ': shape mismatch ', shape(actual), shape(expected)
            nfail = nfail + 1
        else if (maxval(abs(actual - expected)) > tol) then
            write(*, '(a,a,es16.8)') trim(name), ': maximum error ', &
                maxval(abs(actual - expected))
            nfail = nfail + 1
        end if
    end subroutine check_matrix

    subroutine check_int_vector(name, actual, expected, nfail)
        character(len=*), intent(in) :: name !! check label
        integer, intent(in) :: actual(:) !! actual vector
        integer, intent(in) :: expected(:) !! expected vector
        integer, intent(inout) :: nfail !! accumulated failures

        if (size(actual) /= size(expected)) then
            write(*, '(a,a,2i8)') trim(name), ': size mismatch ', size(actual), size(expected)
            nfail = nfail + 1
        else if (any(actual /= expected)) then
            write(*, '(a)') trim(name)//': values mismatch'
            nfail = nfail + 1
        end if
    end subroutine check_int_vector

    subroutine check_real(name, actual, expected, tol, nfail)
        character(len=*), intent(in) :: name !! check label
        real(dp), intent(in) :: actual !! actual value
        real(dp), intent(in) :: expected !! expected value
        real(dp), intent(in) :: tol !! absolute tolerance
        integer, intent(inout) :: nfail !! accumulated failures

        if (abs(actual - expected) > tol) then
            write(*, '(a,a,2es16.8)') trim(name), ': mismatch ', actual, expected
            nfail = nfail + 1
        end if
    end subroutine check_real

    subroutine check_int(name, actual, expected, nfail)
        character(len=*), intent(in) :: name !! check label
        integer, intent(in) :: actual !! actual value
        integer, intent(in) :: expected !! expected value
        integer, intent(inout) :: nfail !! accumulated failures

        if (actual /= expected) then
            write(*, '(a,a,2i8)') trim(name), ': mismatch ', actual, expected
            nfail = nfail + 1
        end if
    end subroutine check_int

end program test_optimize_trust_constr
