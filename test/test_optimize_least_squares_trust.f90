program test_optimize_least_squares_trust
    use scifort_kinds, only : dp
    use scifort_optimize, only : curve_fit, curve_fit_result, least_squares, &
        least_squares_result
    use test_optimize_least_squares_reference
    implicit none

    real(dp), parameter :: lb(2) = [-3.0_dp, -1.0_dp]
    real(dp), parameter :: ub(2) = [0.5_dp, 2.0_dp]
    real(dp), parameter :: xdata(4) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
    real(dp), parameter :: ydata(4) = [1.1_dp, 2.9_dp, 5.2_dp, 6.8_dp]
    real(dp), parameter :: fit_lb(2) = [0.0_dp, 0.0_dp]
    real(dp), parameter :: fit_ub(2) = [1.5_dp, 3.0_dp]
    real(dp), parameter :: robust_lb(1) = [0.0_dp]
    real(dp), parameter :: robust_ub(1) = [1.0_dp]
    logical, parameter :: sparse_pattern(4,4) = reshape([ &
        .true., .true., .false., .false., &
        .true., .true., .false., .false., &
        .false., .false., .true., .true., &
        .false., .false., .true., .true.], [4,4])
    real(dp), parameter :: ydata_outlier(4) = [1.1_dp, 2.9_dp, 5.2_dp, 20.0_dp]
    type(least_squares_result) :: res
    type(curve_fit_result) :: cf
    integer :: callback_calls, callback_last_nfev, callback_last_nit, sparse_fun_calls
    integer :: failures
    integer :: bounded_model_calls
    real(dp) :: callback_last_cost, callback_last_fun(2), callback_last_x(2)

    failures = 0

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2)
    call check_result('default TRF finite difference', res, trf_fd_x_ref, &
        trf_fd_cost_ref, trf_fd_fun_ref, trf_fd_jac_ref, trf_fd_grad_ref, &
        trf_fd_optimality_ref, trf_fd_active_ref, trf_fd_nfev_ref, &
        trf_fd_njev_ref, trf_fd_status_ref, 5.0e-7_dp, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        method='trf', diff_scheme='3-point')
    call check_result('TRF three-point', res, trf_3_x_ref, trf_3_cost_ref, &
        trf_3_fun_ref, trf_3_jac_ref, trf_3_grad_ref, trf_3_optimality_ref, &
        trf_3_active_ref, trf_3_nfev_ref, trf_3_njev_ref, trf_3_status_ref, &
        5.0e-9_dp, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        method='trf', diff_scheme='cs', complex_fun=rosen_residual_complex)
    call check_result('TRF complex-step', res, trf_cs_x_ref, trf_cs_cost_ref, &
        trf_cs_fun_ref, trf_cs_jac_ref, trf_cs_grad_ref, trf_cs_optimality_ref, &
        trf_cs_active_ref, trf_cs_nfev_ref, trf_cs_njev_ref, trf_cs_status_ref, &
        5.0e-12_dp, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        jac=rosen_jac, method='trf')
    call check_result('TRF analytic', res, trf_jac_x_ref, trf_jac_cost_ref, &
        trf_jac_fun_ref, trf_jac_jac_ref, trf_jac_grad_ref, &
        trf_jac_optimality_ref, trf_jac_active_ref, trf_jac_nfev_ref, &
        trf_jac_njev_ref, trf_jac_status_ref, 5.0e-10_dp, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        method='trf', tr_solver='lsmr')
    call check_result('TRF LSMR finite difference', res, trf_lsmr_fd_x_ref, &
        trf_lsmr_fd_cost_ref, trf_lsmr_fd_fun_ref, trf_lsmr_fd_jac_ref, &
        trf_lsmr_fd_grad_ref, trf_lsmr_fd_optimality_ref, trf_lsmr_fd_active_ref, &
        trf_lsmr_fd_nfev_ref, trf_lsmr_fd_njev_ref, trf_lsmr_fd_status_ref, &
        5.0e-7_dp, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        jac=rosen_jac, method='trf', tr_solver='lsmr')
    call check_result('TRF LSMR analytic', res, trf_lsmr_jac_x_ref, &
        trf_lsmr_jac_cost_ref, trf_lsmr_jac_fun_ref, trf_lsmr_jac_jac_ref, &
        trf_lsmr_jac_grad_ref, trf_lsmr_jac_optimality_ref, &
        trf_lsmr_jac_active_ref, trf_lsmr_jac_nfev_ref, trf_lsmr_jac_njev_ref, &
        trf_lsmr_jac_status_ref, 5.0e-10_dp, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        jac_matvec=rosen_jac_matvec, jac_rmatvec=rosen_jac_rmatvec, method='trf')
    call check_operator_result('TRF matrix-free LSMR', res, trf_operator_x_ref, &
        trf_operator_cost_ref, trf_operator_fun_ref, trf_operator_grad_ref, &
        trf_operator_optimality_ref, trf_operator_active_ref, trf_operator_nfev_ref, &
        trf_operator_njev_ref, trf_operator_status_ref, 5.0e-10_dp, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        jac=rosen_jac, method='trf', tr_solver='lsmr', tr_regularize=.false.)
    call check_result('TRF LSMR no regularization', res, trf_lsmr_noreg_x_ref, &
        trf_lsmr_noreg_cost_ref, trf_lsmr_noreg_fun_ref, trf_lsmr_noreg_jac_ref, &
        trf_lsmr_noreg_grad_ref, trf_lsmr_noreg_optimality_ref, &
        trf_lsmr_noreg_active_ref, trf_lsmr_noreg_nfev_ref, &
        trf_lsmr_noreg_njev_ref, trf_lsmr_noreg_status_ref, 5.0e-10_dp, failures)

    sparse_fun_calls = 0
    res = least_squares(sparse_residual, [0.5_dp, 1.5_dp, 2.5_dp, 3.5_dp], 4, &
        method='trf', jac_sparsity=sparse_pattern)
    call check_result('TRF jac_sparsity two-point', res, trf_sparse_2_x_ref, &
        trf_sparse_2_cost_ref, trf_sparse_2_fun_ref, trf_sparse_2_jac_ref, &
        trf_sparse_2_grad_ref, trf_sparse_2_optimality_ref, trf_sparse_2_active_ref, &
        trf_sparse_2_nfev_ref, trf_sparse_2_njev_ref, trf_sparse_2_status_ref, &
        5.0e-8_dp, failures)
    call check_int('TRF jac_sparsity two-point raw calls', sparse_fun_calls, &
        trf_sparse_2_calls_ref, failures)

    sparse_fun_calls = 0
    res = least_squares(sparse_residual, [0.5_dp, 1.5_dp, 2.5_dp, 3.5_dp], 4, &
        method='trf', jac_sparsity=sparse_pattern, diff_scheme='3-point')
    call check_result('TRF jac_sparsity three-point', res, trf_sparse_3_x_ref, &
        trf_sparse_3_cost_ref, trf_sparse_3_fun_ref, trf_sparse_3_jac_ref, &
        trf_sparse_3_grad_ref, trf_sparse_3_optimality_ref, trf_sparse_3_active_ref, &
        trf_sparse_3_nfev_ref, trf_sparse_3_njev_ref, trf_sparse_3_status_ref, &
        5.0e-8_dp, failures)
    call check_int('TRF jac_sparsity three-point raw calls', sparse_fun_calls, &
        trf_sparse_3_calls_ref, failures)

    sparse_fun_calls = 0
    res = least_squares(sparse_residual, [0.5_dp, 1.5_dp, 2.5_dp, 3.5_dp], 4, &
        method='trf', jac_sparsity=sparse_pattern, diff_scheme='cs', &
        complex_fun=sparse_residual_complex)
    call check_result('TRF jac_sparsity complex-step', res, trf_sparse_cs_x_ref, &
        trf_sparse_cs_cost_ref, trf_sparse_cs_fun_ref, trf_sparse_cs_jac_ref, &
        trf_sparse_cs_grad_ref, trf_sparse_cs_optimality_ref, trf_sparse_cs_active_ref, &
        trf_sparse_cs_nfev_ref, trf_sparse_cs_njev_ref, trf_sparse_cs_status_ref, &
        5.0e-12_dp, failures)
    call check_int('TRF jac_sparsity complex-step raw calls', sparse_fun_calls, &
        trf_sparse_cs_calls_ref, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        jac=rosen_jac, method='trf', lower=lb, upper=ub)
    call check_result('TRF bounded', res, trf_bound_x_ref, trf_bound_cost_ref, &
        trf_bound_fun_ref, trf_bound_jac_ref, trf_bound_grad_ref, &
        trf_bound_optimality_ref, trf_bound_active_ref, trf_bound_nfev_ref, &
        trf_bound_njev_ref, trf_bound_status_ref, 2.0e-8_dp, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        jac=rosen_jac, method='trf', lower=lb, upper=ub, tr_solver='lsmr')
    call check_result('TRF bounded LSMR', res, trf_lsmr_bound_x_ref, &
        trf_lsmr_bound_cost_ref, trf_lsmr_bound_fun_ref, trf_lsmr_bound_jac_ref, &
        trf_lsmr_bound_grad_ref, trf_lsmr_bound_optimality_ref, &
        trf_lsmr_bound_active_ref, trf_lsmr_bound_nfev_ref, &
        trf_lsmr_bound_njev_ref, trf_lsmr_bound_status_ref, 2.0e-7_dp, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        jac_matvec=rosen_jac_matvec, jac_rmatvec=rosen_jac_rmatvec, method='trf', &
        lower=lb, upper=ub)
    call check_operator_result('TRF bounded matrix-free LSMR', res, &
        trf_operator_bound_x_ref, trf_operator_bound_cost_ref, trf_operator_bound_fun_ref, &
        trf_operator_bound_grad_ref, trf_operator_bound_optimality_ref, &
        trf_operator_bound_active_ref, trf_operator_bound_nfev_ref, &
        trf_operator_bound_njev_ref, trf_operator_bound_status_ref, 2.0e-7_dp, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        method='trf', lower=lb, upper=ub, diff_scheme='3-point')
    call check_result('TRF bounded three-point', res, trf_bound_3_x_ref, &
        trf_bound_3_cost_ref, trf_bound_3_fun_ref, trf_bound_3_jac_ref, &
        trf_bound_3_grad_ref, trf_bound_3_optimality_ref, trf_bound_3_active_ref, &
        trf_bound_3_nfev_ref, trf_bound_3_njev_ref, trf_bound_3_status_ref, &
        2.0e-7_dp, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        method='trf', lower=lb, upper=ub, diff_scheme='cs', &
        complex_fun=rosen_residual_complex)
    call check_result('TRF bounded complex-step', res, trf_bound_cs_x_ref, &
        trf_bound_cs_cost_ref, trf_bound_cs_fun_ref, trf_bound_cs_jac_ref, &
        trf_bound_cs_grad_ref, trf_bound_cs_optimality_ref, trf_bound_cs_active_ref, &
        trf_bound_cs_nfev_ref, trf_bound_cs_njev_ref, trf_bound_cs_status_ref, &
        2.0e-8_dp, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        jac=rosen_jac, method='trf', x_scale=[2.0_dp, 0.5_dp])
    call check_result('TRF x_scale', res, trf_scale_x_ref, trf_scale_cost_ref, &
        trf_scale_fun_ref, trf_scale_jac_ref, trf_scale_grad_ref, &
        trf_scale_optimality_ref, trf_scale_active_ref, trf_scale_nfev_ref, &
        trf_scale_njev_ref, trf_scale_status_ref, 5.0e-10_dp, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        jac=rosen_jac, method='trf', x_scale='jac')
    call check_result('TRF Jacobian x_scale', res, trf_jac_scale_x_ref, &
        trf_jac_scale_cost_ref, trf_jac_scale_fun_ref, trf_jac_scale_jac_ref, &
        trf_jac_scale_grad_ref, trf_jac_scale_optimality_ref, &
        trf_jac_scale_active_ref, trf_jac_scale_nfev_ref, &
        trf_jac_scale_njev_ref, trf_jac_scale_status_ref, 5.0e-10_dp, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        jac=rosen_jac, method='trf', x_scale='jac', lower=lb, upper=ub)
    call check_result('TRF bounded Jacobian x_scale', res, trf_bound_jac_scale_x_ref, &
        trf_bound_jac_scale_cost_ref, trf_bound_jac_scale_fun_ref, &
        trf_bound_jac_scale_jac_ref, trf_bound_jac_scale_grad_ref, &
        trf_bound_jac_scale_optimality_ref, trf_bound_jac_scale_active_ref, &
        trf_bound_jac_scale_nfev_ref, trf_bound_jac_scale_njev_ref, &
        trf_bound_jac_scale_status_ref, 2.0e-8_dp, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        jac=rosen_jac, method='trf', max_nfev=2)
    call check_int('TRF max_nfev nfev', res%nfev, trf_limit_nfev_ref, failures)
    call check_int('TRF max_nfev njev', res%njev, trf_limit_njev_ref, failures)
    call check_int('TRF max_nfev status', res%status, trf_limit_status_ref, failures)
    call check_real('TRF max_nfev cost', res%cost, trf_limit_cost_ref, 5.0e-2_dp, failures)

    call reset_callback_state()
    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        jac=rosen_jac, method='trf', callback=stop_after_two)
    call check_result('TRF callback stop', res, trf_callback_x_ref, &
        trf_callback_cost_ref, trf_callback_fun_ref, trf_callback_jac_ref, &
        trf_callback_grad_ref, trf_callback_optimality_ref, trf_callback_active_ref, &
        trf_callback_nfev_ref, trf_callback_njev_ref, trf_callback_status_ref, &
        2.0e-1_dp, failures)
    call check_callback_state('TRF callback', trf_callback_x_ref, trf_callback_fun_ref, &
        trf_callback_cost_ref, trf_callback_nfev_ref, 2.0e-1_dp, failures)

    sparse_fun_calls = 0
    res = least_squares(sparse_residual, [0.5_dp, 1.5_dp, 2.5_dp, 3.5_dp], 4, &
        method='dogbox', jac_sparsity=sparse_pattern)
    call check_result('dogbox jac_sparsity two-point', res, dog_sparse_2_x_ref, &
        dog_sparse_2_cost_ref, dog_sparse_2_fun_ref, dog_sparse_2_jac_ref, &
        dog_sparse_2_grad_ref, dog_sparse_2_optimality_ref, dog_sparse_2_active_ref, &
        dog_sparse_2_nfev_ref, dog_sparse_2_njev_ref, dog_sparse_2_status_ref, &
        5.0e-8_dp, failures)
    call check_int('dogbox jac_sparsity two-point raw calls', sparse_fun_calls, &
        dog_sparse_2_calls_ref, failures)

    sparse_fun_calls = 0
    res = least_squares(sparse_residual, [0.5_dp, 1.5_dp, 2.5_dp, 3.5_dp], 4, &
        method='dogbox', jac_sparsity=sparse_pattern, diff_scheme='3-point')
    call check_result('dogbox jac_sparsity three-point', res, dog_sparse_3_x_ref, &
        dog_sparse_3_cost_ref, dog_sparse_3_fun_ref, dog_sparse_3_jac_ref, &
        dog_sparse_3_grad_ref, dog_sparse_3_optimality_ref, dog_sparse_3_active_ref, &
        dog_sparse_3_nfev_ref, dog_sparse_3_njev_ref, dog_sparse_3_status_ref, &
        5.0e-8_dp, failures)
    call check_int('dogbox jac_sparsity three-point raw calls', sparse_fun_calls, &
        dog_sparse_3_calls_ref, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, method='dogbox')
    call check_result('dogbox finite difference', res, dog_fd_x_ref, dog_fd_cost_ref, &
        dog_fd_fun_ref, dog_fd_jac_ref, dog_fd_grad_ref, dog_fd_optimality_ref, &
        dog_fd_active_ref, dog_fd_nfev_ref, dog_fd_njev_ref, dog_fd_status_ref, &
        5.0e-7_dp, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        method='dogbox', tr_solver='lsmr')
    call check_result('dogbox LSMR finite difference', res, dog_lsmr_fd_x_ref, &
        dog_lsmr_fd_cost_ref, dog_lsmr_fd_fun_ref, dog_lsmr_fd_jac_ref, &
        dog_lsmr_fd_grad_ref, dog_lsmr_fd_optimality_ref, dog_lsmr_fd_active_ref, &
        dog_lsmr_fd_nfev_ref, dog_lsmr_fd_njev_ref, dog_lsmr_fd_status_ref, &
        5.0e-7_dp, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        method='dogbox', diff_scheme='3-point')
    call check_result('dogbox three-point', res, dog_3_x_ref, dog_3_cost_ref, &
        dog_3_fun_ref, dog_3_jac_ref, dog_3_grad_ref, dog_3_optimality_ref, &
        dog_3_active_ref, dog_3_nfev_ref, dog_3_njev_ref, dog_3_status_ref, &
        5.0e-9_dp, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        method='dogbox', diff_scheme='cs', complex_fun=rosen_residual_complex)
    call check_result('dogbox complex-step', res, dog_cs_x_ref, dog_cs_cost_ref, &
        dog_cs_fun_ref, dog_cs_jac_ref, dog_cs_grad_ref, dog_cs_optimality_ref, &
        dog_cs_active_ref, dog_cs_nfev_ref, dog_cs_njev_ref, dog_cs_status_ref, &
        5.0e-12_dp, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        jac=rosen_jac, method='dogbox')
    call check_result('dogbox analytic', res, dog_jac_x_ref, dog_jac_cost_ref, &
        dog_jac_fun_ref, dog_jac_jac_ref, dog_jac_grad_ref, dog_jac_optimality_ref, &
        dog_jac_active_ref, dog_jac_nfev_ref, dog_jac_njev_ref, dog_jac_status_ref, &
        5.0e-10_dp, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        jac=rosen_jac, method='dogbox', tr_solver='lsmr')
    call check_result('dogbox LSMR analytic', res, dog_lsmr_jac_x_ref, &
        dog_lsmr_jac_cost_ref, dog_lsmr_jac_fun_ref, dog_lsmr_jac_jac_ref, &
        dog_lsmr_jac_grad_ref, dog_lsmr_jac_optimality_ref, dog_lsmr_jac_active_ref, &
        dog_lsmr_jac_nfev_ref, dog_lsmr_jac_njev_ref, dog_lsmr_jac_status_ref, &
        5.0e-10_dp, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        jac=rosen_jac, method='dogbox', lower=lb, upper=ub)
    call check_result('dogbox bounded', res, dog_bound_x_ref, dog_bound_cost_ref, &
        dog_bound_fun_ref, dog_bound_jac_ref, dog_bound_grad_ref, &
        dog_bound_optimality_ref, dog_bound_active_ref, dog_bound_nfev_ref, &
        dog_bound_njev_ref, dog_bound_status_ref, 5.0e-10_dp, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        jac=rosen_jac, method='dogbox', lower=lb, upper=ub, tr_solver='lsmr')
    call check_result('dogbox bounded LSMR', res, dog_lsmr_bound_x_ref, &
        dog_lsmr_bound_cost_ref, dog_lsmr_bound_fun_ref, dog_lsmr_bound_jac_ref, &
        dog_lsmr_bound_grad_ref, dog_lsmr_bound_optimality_ref, &
        dog_lsmr_bound_active_ref, dog_lsmr_bound_nfev_ref, dog_lsmr_bound_njev_ref, &
        dog_lsmr_bound_status_ref, 5.0e-10_dp, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        method='dogbox', lower=lb, upper=ub, diff_scheme='3-point')
    call check_result('dogbox bounded three-point', res, dog_bound_3_x_ref, &
        dog_bound_3_cost_ref, dog_bound_3_fun_ref, dog_bound_3_jac_ref, &
        dog_bound_3_grad_ref, dog_bound_3_optimality_ref, dog_bound_3_active_ref, &
        dog_bound_3_nfev_ref, dog_bound_3_njev_ref, dog_bound_3_status_ref, &
        5.0e-9_dp, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        method='dogbox', lower=lb, upper=ub, diff_scheme='cs', &
        complex_fun=rosen_residual_complex)
    call check_result('dogbox bounded complex-step', res, dog_bound_cs_x_ref, &
        dog_bound_cs_cost_ref, dog_bound_cs_fun_ref, dog_bound_cs_jac_ref, &
        dog_bound_cs_grad_ref, dog_bound_cs_optimality_ref, dog_bound_cs_active_ref, &
        dog_bound_cs_nfev_ref, dog_bound_cs_njev_ref, dog_bound_cs_status_ref, &
        5.0e-12_dp, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        jac=rosen_jac, method='dogbox', x_scale=[2.0_dp, 0.5_dp])
    call check_result('dogbox x_scale', res, dog_scale_x_ref, dog_scale_cost_ref, &
        dog_scale_fun_ref, dog_scale_jac_ref, dog_scale_grad_ref, &
        dog_scale_optimality_ref, dog_scale_active_ref, dog_scale_nfev_ref, &
        dog_scale_njev_ref, dog_scale_status_ref, 5.0e-10_dp, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        jac=rosen_jac, method='dogbox', x_scale='jac')
    call check_result('dogbox Jacobian x_scale', res, dog_jac_scale_x_ref, &
        dog_jac_scale_cost_ref, dog_jac_scale_fun_ref, dog_jac_scale_jac_ref, &
        dog_jac_scale_grad_ref, dog_jac_scale_optimality_ref, &
        dog_jac_scale_active_ref, dog_jac_scale_nfev_ref, &
        dog_jac_scale_njev_ref, dog_jac_scale_status_ref, 5.0e-10_dp, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        jac=rosen_jac, method='dogbox', max_nfev=2)
    call check_result('dogbox max_nfev', res, dog_limit_x_ref, dog_limit_cost_ref, &
        dog_limit_fun_ref, dog_limit_jac_ref, dog_limit_grad_ref, &
        dog_limit_optimality_ref, dog_limit_active_ref, dog_limit_nfev_ref, &
        dog_limit_njev_ref, dog_limit_status_ref, 5.0e-10_dp, failures)

    call reset_callback_state()
    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        jac=rosen_jac, method='dogbox', callback=stop_after_two)
    call check_result('dogbox callback stop', res, dog_callback_x_ref, &
        dog_callback_cost_ref, dog_callback_fun_ref, dog_callback_jac_ref, &
        dog_callback_grad_ref, dog_callback_optimality_ref, dog_callback_active_ref, &
        dog_callback_nfev_ref, dog_callback_njev_ref, dog_callback_status_ref, &
        5.0e-10_dp, failures)
    call check_callback_state('dogbox callback', dog_callback_x_ref, dog_callback_fun_ref, &
        dog_callback_cost_ref, dog_callback_nfev_ref, 5.0e-10_dp, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        method='trf', diff_scheme='cs')
    call check_int('complex-step missing callback status', res%status, -1, failures)

    res = least_squares(robust_residual, [0.0_dp], 4, jac=robust_jac, &
        method='trf', loss='soft_l1', f_scale=0.5_dp)
    call check_result('TRF soft_l1', res, robust_trf_soft_x_ref, &
        robust_trf_soft_cost_ref, robust_trf_soft_fun_ref, robust_trf_soft_jac_ref, &
        robust_trf_soft_grad_ref, robust_trf_soft_optimality_ref, &
        robust_trf_soft_active_ref, robust_trf_soft_nfev_ref, robust_trf_soft_njev_ref, &
        robust_trf_soft_status_ref, 5.0e-9_dp, failures)

    res = least_squares(robust_residual, [0.0_dp], 4, jac=robust_jac, &
        method='trf', loss_fun=custom_soft_l1, f_scale=0.5_dp)
    call check_result('TRF custom soft_l1', res, robust_trf_soft_x_ref, &
        robust_trf_soft_cost_ref, robust_trf_soft_fun_ref, robust_trf_soft_jac_ref, &
        robust_trf_soft_grad_ref, robust_trf_soft_optimality_ref, &
        robust_trf_soft_active_ref, robust_trf_soft_nfev_ref, robust_trf_soft_njev_ref, &
        robust_trf_soft_status_ref, 5.0e-9_dp, failures)

    res = least_squares(robust_residual, [0.0_dp], 4, jac=robust_jac, &
        method='dogbox', loss_fun=custom_soft_l1, f_scale=0.5_dp)
    call check_result('dogbox custom soft_l1', res, robust_trf_soft_x_ref, &
        robust_trf_soft_cost_ref, robust_trf_soft_fun_ref, robust_trf_soft_jac_ref, &
        robust_trf_soft_grad_ref, robust_trf_soft_optimality_ref, &
        robust_trf_soft_active_ref, robust_trf_soft_nfev_ref, robust_trf_soft_njev_ref, &
        robust_trf_soft_status_ref, 5.0e-9_dp, failures)

    res = least_squares(robust_residual, [0.0_dp], 4, jac=robust_jac, &
        method='trf', loss='cauchy', f_scale=0.5_dp)
    call check_result('TRF cauchy', res, robust_trf_cauchy_x_ref, &
        robust_trf_cauchy_cost_ref, robust_trf_cauchy_fun_ref, robust_trf_cauchy_jac_ref, &
        robust_trf_cauchy_grad_ref, robust_trf_cauchy_optimality_ref, &
        robust_trf_cauchy_active_ref, robust_trf_cauchy_nfev_ref, &
        robust_trf_cauchy_njev_ref, robust_trf_cauchy_status_ref, 1.0e-7_dp, failures)

    res = least_squares(robust_residual, [0.0_dp], 4, jac=robust_jac, &
        method='dogbox', loss='huber', f_scale=0.5_dp)
    call check_result('dogbox huber', res, robust_dog_huber_x_ref, &
        robust_dog_huber_cost_ref, robust_dog_huber_fun_ref, robust_dog_huber_jac_ref, &
        robust_dog_huber_grad_ref, robust_dog_huber_optimality_ref, &
        robust_dog_huber_active_ref, robust_dog_huber_nfev_ref, robust_dog_huber_njev_ref, &
        robust_dog_huber_status_ref, 5.0e-10_dp, failures)

    res = least_squares(robust_residual, [0.0_dp], 4, jac=robust_jac, &
        method='dogbox', loss='arctan', f_scale=0.5_dp)
    call check_result('dogbox arctan', res, robust_dog_arctan_x_ref, &
        robust_dog_arctan_cost_ref, robust_dog_arctan_fun_ref, robust_dog_arctan_jac_ref, &
        robust_dog_arctan_grad_ref, robust_dog_arctan_optimality_ref, &
        robust_dog_arctan_active_ref, robust_dog_arctan_nfev_ref, &
        robust_dog_arctan_njev_ref, robust_dog_arctan_status_ref, 5.0e-10_dp, failures)

    res = least_squares(robust_residual, [0.2_dp], 4, jac=robust_jac, &
        method='trf', loss='soft_l1', f_scale=0.5_dp, lower=robust_lb, upper=robust_ub)
    call check_result('bounded robust TRF', res, robust_trf_bound_x_ref, &
        robust_trf_bound_cost_ref, robust_trf_bound_fun_ref, robust_trf_bound_jac_ref, &
        robust_trf_bound_grad_ref, robust_trf_bound_optimality_ref, &
        robust_trf_bound_active_ref, robust_trf_bound_nfev_ref, robust_trf_bound_njev_ref, &
        robust_trf_bound_status_ref, 5.0e-9_dp, failures)

    bounded_model_calls = 0
    cf = curve_fit(counted_linear_model, xdata, ydata, [1.0_dp, 1.0_dp], &
        lower=fit_lb, upper=fit_ub)
    call check_curve('bounded curve_fit default TRF', cf, cf_bound_popt_ref, &
        cf_bound_pcov_ref, cf_bound_fvec_ref, cf_bound_nfev_ref, &
        cf_bound_status_ref, 5.0e-7_dp, failures, nfev_slack=1)
    call check_int('bounded TRF raw model calls', bounded_model_calls, &
        cf%nfev + size(cf%popt) * cf%njev, failures)

    cf = curve_fit(linear_model, xdata, ydata, [1.0_dp, 1.0_dp], &
        jac=linear_model_jac, method='trf', lower=fit_lb, upper=fit_ub)
    call check_curve('bounded curve_fit TRF analytic', cf, cf_bound_jac_popt_ref, &
        cf_bound_jac_pcov_ref, cf_bound_jac_fvec_ref, cf_bound_jac_nfev_ref, &
        cf_bound_jac_status_ref, 5.0e-9_dp, failures)

    cf = curve_fit(linear_model, xdata, ydata, [1.0_dp, 1.0_dp], &
        jac=linear_model_jac, method='dogbox', lower=fit_lb, upper=fit_ub)
    call check_curve('bounded curve_fit dogbox', cf, cf_bound_dog_popt_ref, &
        cf_bound_dog_pcov_ref, cf_bound_dog_fvec_ref, cf_bound_dog_nfev_ref, &
        cf_bound_dog_status_ref, 5.0e-10_dp, failures)

    cf = curve_fit(linear_model, xdata, ydata, [1.0_dp, 1.0_dp], &
        jac=linear_model_jac, method='trf')
    call check_curve('unbounded curve_fit TRF', cf, cf_trf_popt_ref, &
        cf_trf_pcov_ref, cf_trf_fvec_ref, cf_trf_nfev_ref, &
        cf_trf_status_ref, 5.0e-10_dp, failures)

    cf = curve_fit(linear_model, xdata, ydata, [1.0_dp, 1.0_dp], method='trf', &
        diff_scheme='3-point')
    call check_curve('curve_fit TRF three-point', cf, cf_trf_3_popt_ref, &
        cf_trf_3_pcov_ref, cf_trf_3_fvec_ref, cf_trf_3_nfev_ref, &
        cf_trf_3_status_ref, 5.0e-12_dp, failures)

    cf = curve_fit(linear_model, xdata, ydata, [1.0_dp, 1.0_dp], method='trf', &
        diff_scheme='cs', complex_model=linear_model_complex)
    call check_curve('curve_fit TRF complex-step', cf, cf_trf_cs_popt_ref, &
        cf_trf_cs_pcov_ref, cf_trf_cs_fvec_ref, cf_trf_cs_nfev_ref, &
        cf_trf_cs_status_ref, 5.0e-12_dp, failures)

    cf = curve_fit(linear_model, xdata, ydata, [1.0_dp, 1.0_dp], method='dogbox', &
        lower=fit_lb, upper=fit_ub, diff_scheme='3-point')
    call check_curve('bounded curve_fit dogbox three-point', cf, &
        cf_bound_dog_3_popt_ref, cf_bound_dog_3_pcov_ref, cf_bound_dog_3_fvec_ref, &
        cf_bound_dog_3_nfev_ref, cf_bound_dog_3_status_ref, 5.0e-12_dp, failures)

    cf = curve_fit(linear_model, xdata, ydata, [1.0_dp, 1.0_dp], method='dogbox', &
        lower=fit_lb, upper=fit_ub, diff_scheme='cs', complex_model=linear_model_complex)
    call check_curve('bounded curve_fit dogbox complex-step', cf, &
        cf_bound_dog_cs_popt_ref, cf_bound_dog_cs_pcov_ref, cf_bound_dog_cs_fvec_ref, &
        cf_bound_dog_cs_nfev_ref, cf_bound_dog_cs_status_ref, 5.0e-12_dp, failures)

    cf = curve_fit(linear_model, xdata, ydata_outlier, [2.0_dp, 1.0_dp], &
        jac=linear_model_jac, method='trf', loss='soft_l1', f_scale=0.5_dp)
    call check_curve('robust curve_fit TRF', cf, cf_robust_trf_popt_ref, &
        cf_robust_trf_pcov_ref, cf_robust_trf_fvec_ref, cf_robust_trf_nfev_ref, &
        cf_robust_trf_status_ref, 5.0e-8_dp, failures)

    cf = curve_fit(linear_model, xdata, ydata_outlier, [2.0_dp, 1.0_dp], &
        jac=linear_model_jac, method='dogbox', loss='soft_l1', f_scale=0.5_dp)
    call check_curve('robust curve_fit dogbox', cf, cf_robust_dog_popt_ref, &
        cf_robust_dog_pcov_ref, cf_robust_dog_fvec_ref, cf_robust_dog_nfev_ref, &
        cf_robust_dog_status_ref, 5.0e-8_dp, failures)

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        jac_matvec=rosen_jac_matvec, method='trf')
    if (res%status /= -1 .or. index(res%message, 'both jac_matvec and jac_rmatvec') == 0) then
        write(*, '(a)') 'matrix-free callback-pair validation failed'
        failures = failures + 1
    end if

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        jac_matvec=rosen_jac_matvec, jac_rmatvec=rosen_jac_rmatvec, method='trf', &
        tr_solver='exact')
    if (res%status /= -1 .or. index(res%message, 'incompatible') == 0) then
        write(*, '(a)') 'matrix-free exact-solver validation failed'
        failures = failures + 1
    end if

    res = least_squares(sparse_residual, [0.5_dp, 1.5_dp, 2.5_dp, 3.5_dp], 4, &
        method='trf', jac_sparsity=sparse_pattern, tr_solver='exact')
    if (res%status /= -1 .or. index(res%message, 'incompatible') == 0) then
        write(*, '(a)') 'jac_sparsity exact-solver validation failed'
        failures = failures + 1
    end if

    res = least_squares(sparse_residual, [0.5_dp, 1.5_dp, 2.5_dp, 3.5_dp], 4, &
        method='dogbox', jac_sparsity=sparse_pattern, tr_solver='exact')
    if (res%status /= -1 .or. index(res%message, 'incompatible') == 0) then
        write(*, '(a)') 'dogbox jac_sparsity exact-solver validation failed'
        failures = failures + 1
    end if

    res = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        jac_matvec=rosen_jac_matvec, jac_rmatvec=rosen_jac_rmatvec, method='trf', &
        x_scale='jac')
    if (res%status /= -1 .or. index(res%message, "x_scale='jac'") == 0) then
        write(*, '(a)') 'matrix-free Jacobian scaling validation failed'
        failures = failures + 1
    end if

    if (failures /= 0) then
        write(*, '(a,i0)') 'test_optimize_least_squares_trust: FAIL ', failures
        stop 1
    end if
    write(*, '(a)') 'test_optimize_least_squares_trust: PASS'

contains

    subroutine reset_callback_state()
        callback_calls = 0
        callback_last_nit = 0
        callback_last_nfev = 0
        callback_last_cost = 0.0_dp
        callback_last_x = 0.0_dp
        callback_last_fun = 0.0_dp
    end subroutine reset_callback_state

    logical function stop_after_two(x, fvec, cost, nit, nfev) result(stop)
        real(dp), intent(in) :: x(:) !! current parameter vector
        real(dp), intent(in) :: fvec(:) !! current true residual vector
        real(dp), intent(in) :: cost !! current cost
        integer, intent(in) :: nit !! completed outer iterations
        integer, intent(in) :: nfev !! residual evaluation count

        callback_calls = callback_calls + 1
        callback_last_nit = nit
        callback_last_nfev = nfev
        callback_last_cost = cost
        callback_last_x = x
        callback_last_fun = fvec
        stop = nit >= 2
    end function stop_after_two

    subroutine check_callback_state(name, expected_x, expected_fun, expected_cost, &
            expected_nfev, tol, nfail)
        character(len=*), intent(in) :: name !! check label
        real(dp), intent(in) :: expected_x(:), expected_fun(:), expected_cost
        integer, intent(in) :: expected_nfev
        real(dp), intent(in) :: tol !! absolute callback-state tolerance
        integer, intent(inout) :: nfail !! accumulated failure count

        call check_int(trim(name)//' calls', callback_calls, 2, nfail)
        call check_int(trim(name)//' nit', callback_last_nit, 2, nfail)
        call check_int(trim(name)//' nfev', callback_last_nfev, expected_nfev, nfail)
        call check_real(trim(name)//' cost', callback_last_cost, expected_cost, tol, nfail)
        call check_vector(trim(name)//' x', callback_last_x, expected_x, tol, nfail)
        call check_vector(trim(name)//' fun', callback_last_fun, expected_fun, tol, nfail)
    end subroutine check_callback_state

    subroutine rosen_residual(x, fvec)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: fvec(:) !! Rosenbrock residuals

        fvec = [10.0_dp * (x(2) - x(1)**2), 1.0_dp - x(1)]
    end subroutine rosen_residual

    subroutine rosen_residual_complex(x, fvec)
        complex(dp), intent(in) :: x(:) !! complex-step parameter vector
        complex(dp), intent(out) :: fvec(:) !! analytic Rosenbrock residual continuation

        fvec = [10.0_dp * (x(2) - x(1)**2), 1.0_dp - x(1)]
    end subroutine rosen_residual_complex

    subroutine rosen_jac(x, fjac)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: fjac(:, :) !! Rosenbrock residual Jacobian

        fjac(1, :) = [-20.0_dp * x(1), 10.0_dp]
        fjac(2, :) = [-1.0_dp, 0.0_dp]
    end subroutine rosen_jac

    subroutine rosen_jac_matvec(x, v, jv)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(in) :: v(:) !! parameter-space vector
        real(dp), intent(out) :: jv(:) !! Jacobian-vector product

        jv(1) = -20.0_dp * x(1) * v(1) + 10.0_dp * v(2)
        jv(2) = -v(1)
    end subroutine rosen_jac_matvec

    subroutine rosen_jac_rmatvec(x, w, jtw)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(in) :: w(:) !! residual-space vector
        real(dp), intent(out) :: jtw(:) !! transposed-Jacobian-vector product

        jtw(1) = -20.0_dp * x(1) * w(1) - w(2)
        jtw(2) = 10.0_dp * w(1)
    end subroutine rosen_jac_rmatvec

    subroutine sparse_residual(x, fvec)
        real(dp), intent(in) :: x(:) !! four-parameter block-diagonal test point
        real(dp), intent(out) :: fvec(:) !! residuals with two independent parameter pairs

        sparse_fun_calls = sparse_fun_calls + 1
        fvec = [x(1) + x(2) - 3.0_dp, x(1) - x(2) + 1.0_dp, &
            x(3) + x(4) - 7.0_dp, x(3) - x(4) + 1.0_dp]
    end subroutine sparse_residual

    subroutine sparse_residual_complex(x, fvec)
        complex(dp), intent(in) :: x(:) !! complex-step block-diagonal test point
        complex(dp), intent(out) :: fvec(:) !! analytic continuation of sparse residuals

        sparse_fun_calls = sparse_fun_calls + 1
        fvec = [x(1) + x(2) - 3.0_dp, x(1) - x(2) + 1.0_dp, &
            x(3) + x(4) - 7.0_dp, x(3) - x(4) + 1.0_dp]
    end subroutine sparse_residual_complex

    subroutine robust_residual(x, fvec)
        real(dp), intent(in) :: x(:) !! scalar location parameter
        real(dp), intent(out) :: fvec(:) !! residuals with one large outlier

        fvec = x(1) - [1.0_dp, 1.2_dp, 0.8_dp, 10.0_dp]
    end subroutine robust_residual

    subroutine robust_jac(x, fjac)
        real(dp), intent(in) :: x(:) !! scalar location parameter
        real(dp), intent(out) :: fjac(:, :) !! constant residual Jacobian

        fjac = 1.0_dp
        if (size(x) < 0) fjac = 0.0_dp
    end subroutine robust_jac

    subroutine custom_soft_l1(z, rho0, rho1, rho2)
        real(dp), intent(in) :: z(:) !! squared scaled residuals
        real(dp), intent(out) :: rho0(:) !! loss values
        real(dp), intent(out) :: rho1(:) !! first derivatives
        real(dp), intent(out) :: rho2(:) !! second derivatives

        real(dp), allocatable :: t(:)

        allocate(t(size(z)))
        t = 1.0_dp + z
        rho0 = 2.0_dp * (sqrt(t) - 1.0_dp)
        rho1 = 1.0_dp / sqrt(t)
        rho2 = -0.5_dp / (t * sqrt(t))
    end subroutine custom_soft_l1

    subroutine counted_linear_model(x, params, ymodel)
        implicit none
        real(dp), intent(in) :: x(:) !! independent data
        real(dp), intent(in) :: params(:) !! slope and intercept
        real(dp), intent(out) :: ymodel(:) !! model values

        bounded_model_calls = bounded_model_calls + 1
        call linear_model(x, params, ymodel)
    end subroutine counted_linear_model

    subroutine linear_model(x, params, ymodel)
        real(dp), intent(in) :: x(:) !! independent data
        real(dp), intent(in) :: params(:) !! slope and intercept
        real(dp), intent(out) :: ymodel(:) !! model values

        ymodel = params(1) * x + params(2)
    end subroutine linear_model

    subroutine linear_model_complex(x, params, ymodel)
        real(dp), intent(in) :: x(:) !! independent data
        complex(dp), intent(in) :: params(:) !! complex slope and intercept
        complex(dp), intent(out) :: ymodel(:) !! analytic model continuation

        ymodel = params(1) * x + params(2)
    end subroutine linear_model_complex

    subroutine linear_model_jac(x, params, fjac)
        real(dp), intent(in) :: x(:) !! independent data
        real(dp), intent(in) :: params(:) !! slope and intercept
        real(dp), intent(out) :: fjac(:, :) !! model Jacobian

        if (size(params) < 0) fjac = 0.0_dp
        fjac(:, 1) = x
        fjac(:, 2) = 1.0_dp
    end subroutine linear_model_jac

    subroutine check_operator_result(name, actual, expected_x, expected_cost, expected_fun, &
            expected_grad, expected_optimality, expected_active, expected_nfev, &
            expected_njev, expected_status, tol, nfail)
        character(len=*), intent(in) :: name !! check label
        type(least_squares_result), intent(in) :: actual !! actual result
        real(dp), intent(in) :: expected_x(:), expected_cost, expected_fun(:)
        real(dp), intent(in) :: expected_grad(:), expected_optimality
        integer, intent(in) :: expected_active(:), expected_nfev, expected_njev
        integer, intent(in) :: expected_status
        real(dp), intent(in) :: tol !! absolute tolerance
        integer, intent(inout) :: nfail !! accumulated failure count

        call check_vector(trim(name)//' x', actual%x, expected_x, tol, nfail)
        call check_real(trim(name)//' cost', actual%cost, expected_cost, tol, nfail)
        call check_vector(trim(name)//' fun', actual%fun, expected_fun, tol, nfail)
        call check_vector(trim(name)//' grad', actual%grad, expected_grad, tol, nfail)
        call check_real(trim(name)//' optimality', actual%optimality, expected_optimality, tol, nfail)
        call check_int_vector(trim(name)//' active', actual%active_mask, expected_active, nfail)
        call check_int(trim(name)//' nfev', actual%nfev, expected_nfev, nfail)
        call check_int(trim(name)//' njev', actual%njev, expected_njev, nfail)
        call check_int(trim(name)//' status', actual%status, expected_status, nfail)
        if (.not. actual%jac_is_operator .or. allocated(actual%jac)) then
            write(*, '(a)') trim(name)//' Jacobian representation mismatch'
            nfail = nfail + 1
        end if
    end subroutine check_operator_result

    subroutine check_result(name, actual, expected_x, expected_cost, expected_fun, &
            expected_jac, expected_grad, expected_optimality, expected_active, &
            expected_nfev, expected_njev, expected_status, tol, nfail)
        character(len=*), intent(in) :: name !! check label
        type(least_squares_result), intent(in) :: actual !! actual result
        real(dp), intent(in) :: expected_x(:), expected_cost, expected_fun(:)
        real(dp), intent(in) :: expected_jac(:, :), expected_grad(:), expected_optimality
        integer, intent(in) :: expected_active(:), expected_nfev, expected_njev
        integer, intent(in) :: expected_status
        real(dp), intent(in) :: tol !! absolute tolerance
        integer, intent(inout) :: nfail !! accumulated failure count

        call check_vector(trim(name)//' x', actual%x, expected_x, tol, nfail)
        call check_real(trim(name)//' cost', actual%cost, expected_cost, tol, nfail)
        call check_vector(trim(name)//' fun', actual%fun, expected_fun, tol, nfail)
        call check_matrix(trim(name)//' jac', actual%jac, expected_jac, tol, nfail)
        call check_vector(trim(name)//' grad', actual%grad, expected_grad, tol, nfail)
        call check_real(trim(name)//' optimality', actual%optimality, &
            expected_optimality, tol, nfail)
        call check_int_vector(trim(name)//' active', actual%active_mask, &
            expected_active, nfail)
        call check_int(trim(name)//' nfev', actual%nfev, expected_nfev, nfail)
        call check_int(trim(name)//' njev', actual%njev, expected_njev, nfail)
        call check_int(trim(name)//' status', actual%status, expected_status, nfail)
    end subroutine check_result

    subroutine check_curve(name, actual, expected_popt, expected_pcov, &
            expected_fvec, expected_nfev, expected_status, tol, nfail, nfev_slack)
        character(len=*), intent(in) :: name !! check label
        type(curve_fit_result), intent(in) :: actual !! actual curve-fit result
        real(dp), intent(in) :: expected_popt(:), expected_pcov(:, :), expected_fvec(:)
        integer, intent(in) :: expected_nfev, expected_status
        real(dp), intent(in) :: tol !! absolute tolerance
        integer, intent(inout) :: nfail !! accumulated failure count
        integer, intent(in), optional :: nfev_slack !! allowed platform work-count variation

        call check_vector(trim(name)//' popt', actual%popt, expected_popt, tol, nfail)
        call check_matrix(trim(name)//' pcov', actual%pcov, expected_pcov, tol, nfail)
        call check_vector(trim(name)//' fvec', actual%fvec, expected_fvec, tol, nfail)
        if (present(nfev_slack)) then
            ! Finite-difference rounding can move gradient termination by one iteration.
            if (abs(actual%nfev - expected_nfev) > nfev_slack) then
                call check_int(trim(name)//' nfev', actual%nfev, expected_nfev, nfail)
            end if
        else
            call check_int(trim(name)//' nfev', actual%nfev, expected_nfev, nfail)
        end if
        call check_int(trim(name)//' status', actual%status, expected_status, nfail)
    end subroutine check_curve

    subroutine check_vector(name, actual, expected, tol, nfail)
        character(len=*), intent(in) :: name !! check label
        real(dp), intent(in) :: actual(:), expected(:), tol
        integer, intent(inout) :: nfail !! accumulated failure count

        if (size(actual) /= size(expected) .or. any(abs(actual - expected) > tol)) then
            write(*, '(a,*(1x,es24.16))') trim(name)//' actual:', actual
            write(*, '(a,*(1x,es24.16))') trim(name)//' expected:', expected
            nfail = nfail + 1
        end if
    end subroutine check_vector

    subroutine check_matrix(name, actual, expected, tol, nfail)
        character(len=*), intent(in) :: name !! check label
        real(dp), intent(in) :: actual(:, :), expected(:, :), tol
        integer, intent(inout) :: nfail !! accumulated failure count

        if (any(shape(actual) /= shape(expected)) .or. any(abs(actual - expected) > tol)) then
            write(*, '(a)') trim(name)//' matrix mismatch'
            nfail = nfail + 1
        end if
    end subroutine check_matrix

    subroutine check_real(name, actual, expected, tol, nfail)
        character(len=*), intent(in) :: name !! check label
        real(dp), intent(in) :: actual, expected, tol
        integer, intent(inout) :: nfail !! accumulated failure count

        if (abs(actual - expected) > tol) then
            write(*, '(a,2(1x,es24.16))') trim(name)//' mismatch:', actual, expected
            nfail = nfail + 1
        end if
    end subroutine check_real

    subroutine check_int(name, actual, expected, nfail)
        character(len=*), intent(in) :: name !! check label
        integer, intent(in) :: actual, expected
        integer, intent(inout) :: nfail !! accumulated failure count

        if (actual /= expected) then
            write(*, '(a,2(1x,i0))') trim(name)//' mismatch:', actual, expected
            nfail = nfail + 1
        end if
    end subroutine check_int

    subroutine check_int_vector(name, actual, expected, nfail)
        character(len=*), intent(in) :: name !! check label
        integer, intent(in) :: actual(:), expected(:)
        integer, intent(inout) :: nfail !! accumulated failure count

        if (size(actual) /= size(expected) .or. any(actual /= expected)) then
            write(*, '(a,*(1x,i0))') trim(name)//' actual:', actual
            write(*, '(a,*(1x,i0))') trim(name)//' expected:', expected
            nfail = nfail + 1
        end if
    end subroutine check_int_vector

end program test_optimize_least_squares_trust
