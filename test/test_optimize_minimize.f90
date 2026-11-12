program test_optimize_minimize
    use scifort_kinds, only : dp
    use scifort_optimize, only : bfgs, cg, dogleg, minimize, &
        multivariate_minimize_result, nelder_mead, newton_cg, powell, trust_exact, &
        trust_krylov, trust_ncg
    use test_optimize_minimize_reference
    implicit none

    type(multivariate_minimize_result) :: r
    integer :: failures

    failures = 0

    r = minimize(rosen, [-1.2_dp, 1.0_dp], method='nelder-mead')
    call check_result('nelder-mead rosen', r, nm_rosen_x_ref, nm_rosen_fun_ref, &
        nm_rosen_nit_ref, nm_rosen_nfev_ref, nm_rosen_njev_ref, nm_rosen_status_ref, &
        5.0e-13_dp, failures)

    r = nelder_mead(rosen, [-1.2_dp, 1.0_dp], maxiter=2)
    call check_result('nelder-mead maxiter', r, nm_limit_x_ref, nm_limit_fun_ref, &
        nm_limit_nit_ref, nm_limit_nfev_ref, nm_limit_njev_ref, nm_limit_status_ref, &
        5.0e-13_dp, failures)

    r = minimize(rosen, [-1.2_dp, 1.0_dp], method='powell')
    call check_result('powell rosen', r, powell_rosen_x_ref, powell_rosen_fun_ref, &
        powell_rosen_nit_ref, powell_rosen_nfev_ref, powell_rosen_njev_ref, &
        powell_rosen_status_ref, 5.0e-12_dp, failures)

    r = powell(rosen, [-1.2_dp, 1.0_dp], maxiter=1)
    call check_result('powell maxiter', r, powell_limit_x_ref, powell_limit_fun_ref, &
        powell_limit_nit_ref, powell_limit_nfev_ref, powell_limit_njev_ref, &
        powell_limit_status_ref, 5.0e-12_dp, failures)

    r = minimize(box_quadratic, [0.0_dp, 0.0_dp], method='nelder-mead', &
        lower=[-1.0_dp, -2.0_dp], upper=[1.0_dp, 0.0_dp])
    call check_result('nelder-mead box bounds', r, nm_box_x_ref, nm_box_fun_ref, &
        nm_box_nit_ref, nm_box_nfev_ref, nm_box_njev_ref, nm_box_status_ref, &
        5.0e-12_dp, failures)

    r = nelder_mead(box_quadratic, [2.0_dp, -3.0_dp], &
        lower=[-1.0_dp, -2.0_dp], upper=[1.0_dp, 0.0_dp])
    call check_result('nelder-mead clipped start', r, nm_box_outside_x_ref, &
        nm_box_outside_fun_ref, nm_box_outside_nit_ref, nm_box_outside_nfev_ref, &
        nm_box_outside_njev_ref, nm_box_outside_status_ref, 5.0e-12_dp, failures)

    r = nelder_mead(box_quadratic, [0.5_dp, -0.5_dp], lower=[0.0_dp, -10.0_dp])
    call check_result('nelder-mead lower only', r, nm_lower_only_x_ref, &
        nm_lower_only_fun_ref, nm_lower_only_nit_ref, nm_lower_only_nfev_ref, &
        nm_lower_only_njev_ref, nm_lower_only_status_ref, 5.0e-12_dp, failures)

    r = minimize(box_quadratic, [0.0_dp, 0.0_dp], method='powell', &
        lower=[-1.0_dp, -2.0_dp], upper=[1.0_dp, 0.0_dp])
    call check_result('powell box bounds', r, powell_box_x_ref, powell_box_fun_ref, &
        powell_box_nit_ref, powell_box_nfev_ref, powell_box_njev_ref, &
        powell_box_status_ref, 5.0e-10_dp, failures)

    r = powell(box_quadratic, [2.0_dp, -3.0_dp], lower=[-1.0_dp, -2.0_dp], &
        upper=[1.0_dp, 0.0_dp])
    call check_result('powell outside start', r, powell_box_outside_x_ref, &
        powell_box_outside_fun_ref, powell_box_outside_nit_ref, &
        powell_box_outside_nfev_ref, powell_box_outside_njev_ref, &
        powell_box_outside_status_ref, 5.0e-8_dp, failures)

    r = powell(box_quadratic, [0.5_dp, -0.5_dp], lower=[0.0_dp, -10.0_dp])
    call check_result('powell lower only', r, powell_lower_only_x_ref, &
        powell_lower_only_fun_ref, powell_lower_only_nit_ref, &
        powell_lower_only_nfev_ref, powell_lower_only_njev_ref, &
        powell_lower_only_status_ref, 5.0e-8_dp, failures)

    r = powell(box_quadratic, [2.0_dp, -1.0_dp], &
        lower=[-1.0_dp, -2.0_dp], upper=[1.0_dp, 0.0_dp], &
        direc=reshape([0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2], order=[2, 1]))
    call check_result('powell out of bounds status', r, powell_out_of_bounds_x_ref, &
        powell_out_of_bounds_fun_ref, powell_out_of_bounds_nit_ref, &
        powell_out_of_bounds_nfev_ref, powell_out_of_bounds_njev_ref, &
        powell_out_of_bounds_status_ref, 5.0e-12_dp, failures)

    r = minimize(rosen, [-1.2_dp, 1.0_dp], method='bfgs', jac=rosen_gradient)
    call check_result('bfgs analytic rosen', r, bfgs_rosen_jac_x_ref, &
        bfgs_rosen_jac_fun_ref, bfgs_rosen_jac_nit_ref, bfgs_rosen_jac_nfev_ref, &
        bfgs_rosen_jac_njev_ref, bfgs_rosen_jac_status_ref, 5.0e-11_dp, failures)

    r = bfgs(rosen, [-1.2_dp, 1.0_dp], jac=rosen_gradient, maxiter=2)
    call check_result('bfgs maxiter', r, bfgs_limit_x_ref, bfgs_limit_fun_ref, &
        bfgs_limit_nit_ref, bfgs_limit_nfev_ref, bfgs_limit_njev_ref, &
        bfgs_limit_status_ref, 5.0e-11_dp, failures)

    r = minimize(rosen, [-1.2_dp, 1.0_dp], method='cg', jac=rosen_gradient)
    call check_result('cg analytic rosen', r, cg_rosen_jac_x_ref, cg_rosen_jac_fun_ref, &
        cg_rosen_jac_nit_ref, cg_rosen_jac_nfev_ref, cg_rosen_jac_njev_ref, &
        cg_rosen_jac_status_ref, 5.0e-11_dp, failures)

    r = cg(rosen, [-1.2_dp, 1.0_dp], jac=rosen_gradient, maxiter=2)
    call check_result('cg maxiter', r, cg_limit_x_ref, cg_limit_fun_ref, &
        cg_limit_nit_ref, cg_limit_nfev_ref, cg_limit_njev_ref, cg_limit_status_ref, &
        5.0e-11_dp, failures)

    r = nelder_mead(quadratic, [2.0_dp, -1.0_dp])
    call check_result('nelder-mead quadratic', r, nm_quad_x_ref, nm_quad_fun_ref, &
        nm_quad_nit_ref, nm_quad_nfev_ref, nm_quad_njev_ref, nm_quad_status_ref, &
        5.0e-12_dp, failures)

    r = powell(quadratic, [2.0_dp, -1.0_dp])
    call check_result('powell quadratic', r, powell_quad_x_ref, powell_quad_fun_ref, &
        powell_quad_nit_ref, -1, powell_quad_njev_ref, &
        powell_quad_status_ref, 5.0e-8_dp, failures)

    r = bfgs(quadratic, [2.0_dp, -1.0_dp])
    call check_result('bfgs finite-difference quadratic', r, bfgs_quad_fd_x_ref, &
        bfgs_quad_fd_fun_ref, bfgs_quad_fd_nit_ref, bfgs_quad_fd_nfev_ref, &
        bfgs_quad_fd_njev_ref, bfgs_quad_fd_status_ref, 5.0e-8_dp, failures)

    r = cg(quadratic, [2.0_dp, -1.0_dp])
    call check_result('cg finite-difference quadratic', r, cg_quad_fd_x_ref, &
        cg_quad_fd_fun_ref, cg_quad_fd_nit_ref, cg_quad_fd_nfev_ref, &
        cg_quad_fd_njev_ref, cg_quad_fd_status_ref, 5.0e-8_dp, failures)

    r = minimize(rosen, [-1.2_dp, 1.0_dp], method='newton-cg', &
        jac=rosen_gradient, hess=rosen_hessian)
    call check_result('newton-cg hess', r, newton_rosen_hess_x_ref, &
        newton_rosen_hess_fun_ref, newton_rosen_hess_nit_ref, &
        newton_rosen_hess_nfev_ref, newton_rosen_hess_njev_ref, &
        newton_rosen_hess_status_ref, 5.0e-11_dp, failures)
    call check_int('newton-cg hess nhev', r%nhev, newton_rosen_hess_nhev_ref, failures)

    r = newton_cg(rosen, [-1.2_dp, 1.0_dp], rosen_gradient, hessp=rosen_hess_product)
    call check_result('newton-cg hessp', r, newton_rosen_hessp_x_ref, &
        newton_rosen_hessp_fun_ref, newton_rosen_hessp_nit_ref, &
        newton_rosen_hessp_nfev_ref, newton_rosen_hessp_njev_ref, &
        newton_rosen_hessp_status_ref, 5.0e-11_dp, failures)
    call check_int('newton-cg hessp nhev', r%nhev, newton_rosen_hessp_nhev_ref, failures)

    r = newton_cg(rosen, [-1.2_dp, 1.0_dp], rosen_gradient, hess=rosen_hessian, &
        maxiter=1)
    call check_result('newton-cg maxiter', r, newton_limit_x_ref, newton_limit_fun_ref, &
        newton_limit_nit_ref, newton_limit_nfev_ref, newton_limit_njev_ref, &
        newton_limit_status_ref, 5.0e-12_dp, failures)
    call check_int('newton-cg maxiter nhev', r%nhev, newton_limit_nhev_ref, failures)

    r = newton_cg(quadratic, [2.0_dp, -1.0_dp], quadratic_gradient)
    call check_result('newton-cg fd hessp', r, newton_quad_fd_hessp_x_ref, &
        newton_quad_fd_hessp_fun_ref, newton_quad_fd_hessp_nit_ref, &
        -1, -1, newton_quad_fd_hessp_status_ref, 1.0e-8_dp, failures)
    call check_int('newton-cg fd hessp nhev', r%nhev, &
        newton_quad_fd_hessp_nhev_ref, failures)

    r = minimize(rosen, [-1.2_dp, 1.0_dp], method='dogleg', &
        jac=rosen_gradient, hess=rosen_hessian)
    call check_result('dogleg hess', r, dogleg_rosen_hess_x_ref, &
        dogleg_rosen_hess_fun_ref, dogleg_rosen_hess_nit_ref, &
        dogleg_rosen_hess_nfev_ref, dogleg_rosen_hess_njev_ref, &
        dogleg_rosen_hess_status_ref, 5.0e-11_dp, failures)
    call check_int('dogleg hess nhev', r%nhev, dogleg_rosen_hess_nhev_ref, failures)

    r = dogleg(rosen, [-1.2_dp, 1.0_dp], rosen_gradient, rosen_hessian, maxiter=1)
    call check_result('dogleg maxiter', r, dogleg_limit_x_ref, dogleg_limit_fun_ref, &
        dogleg_limit_nit_ref, dogleg_limit_nfev_ref, dogleg_limit_njev_ref, &
        dogleg_limit_status_ref, 5.0e-12_dp, failures)
    call check_int('dogleg maxiter nhev', r%nhev, dogleg_limit_nhev_ref, failures)

    r = dogleg(rosen, [-1.2_dp, 1.0_dp], rosen_gradient, rosen_hessian, &
        initial_trust_radius=0.25_dp, max_trust_radius=10.0_dp, eta=0.2_dp)
    call check_result('dogleg custom radii', r, dogleg_custom_x_ref, &
        dogleg_custom_fun_ref, dogleg_custom_nit_ref, dogleg_custom_nfev_ref, &
        dogleg_custom_njev_ref, dogleg_custom_status_ref, 5.0e-11_dp, failures)
    call check_int('dogleg custom nhev', r%nhev, dogleg_custom_nhev_ref, failures)

    r = dogleg(quadratic, [2.0_dp, -1.0_dp], quadratic_gradient, quadratic_hessian)
    call check_result('dogleg quadratic', r, dogleg_quad_hess_x_ref, &
        dogleg_quad_hess_fun_ref, dogleg_quad_hess_nit_ref, dogleg_quad_hess_nfev_ref, &
        dogleg_quad_hess_njev_ref, dogleg_quad_hess_status_ref, 5.0e-12_dp, failures)
    call check_int('dogleg quadratic nhev', r%nhev, dogleg_quad_hess_nhev_ref, failures)

    r = dogleg(indefinite_quadratic, [1.0_dp, 1.0_dp], &
        indefinite_quadratic_gradient, indefinite_quadratic_hessian, maxiter=3)
    call check_result('dogleg indefinite', r, dogleg_indefinite_x_ref, &
        dogleg_indefinite_fun_ref, dogleg_indefinite_nit_ref, &
        dogleg_indefinite_nfev_ref, dogleg_indefinite_njev_ref, &
        dogleg_indefinite_status_ref, 5.0e-12_dp, failures)
    call check_int('dogleg indefinite nhev', r%nhev, dogleg_indefinite_nhev_ref, failures)

    r = minimize(rosen, [-1.2_dp, 1.0_dp], method='trust-ncg', &
        jac=rosen_gradient, hess=rosen_hessian)
    call check_result('trust-ncg hess', r, trust_ncg_rosen_hess_x_ref, &
        trust_ncg_rosen_hess_fun_ref, trust_ncg_rosen_hess_nit_ref, &
        trust_ncg_rosen_hess_nfev_ref, trust_ncg_rosen_hess_njev_ref, &
        trust_ncg_rosen_hess_status_ref, 5.0e-11_dp, failures)
    call check_int('trust-ncg hess nhev', r%nhev, trust_ncg_rosen_hess_nhev_ref, failures)

    r = trust_ncg(rosen, [-1.2_dp, 1.0_dp], rosen_gradient, &
        hessp=rosen_hess_product)
    call check_result('trust-ncg hessp', r, trust_ncg_rosen_hessp_x_ref, &
        trust_ncg_rosen_hessp_fun_ref, trust_ncg_rosen_hessp_nit_ref, &
        trust_ncg_rosen_hessp_nfev_ref, trust_ncg_rosen_hessp_njev_ref, &
        trust_ncg_rosen_hessp_status_ref, 5.0e-11_dp, failures)
    call check_int('trust-ncg hessp nhev', r%nhev, trust_ncg_rosen_hessp_nhev_ref, &
        failures)

    r = trust_ncg(rosen, [-1.2_dp, 1.0_dp], rosen_gradient, hess=rosen_hessian, &
        maxiter=1)
    call check_result('trust-ncg maxiter', r, trust_ncg_limit_x_ref, &
        trust_ncg_limit_fun_ref, trust_ncg_limit_nit_ref, trust_ncg_limit_nfev_ref, &
        trust_ncg_limit_njev_ref, trust_ncg_limit_status_ref, 5.0e-12_dp, failures)
    call check_int('trust-ncg maxiter nhev', r%nhev, trust_ncg_limit_nhev_ref, failures)

    r = minimize(rosen, [-1.2_dp, 1.0_dp], method='trust_ncg', &
        jac=rosen_gradient, hessp=rosen_hess_product, initial_trust_radius=0.25_dp, &
        max_trust_radius=10.0_dp, eta=0.2_dp)
    call check_result('trust-ncg custom radii', r, trust_ncg_custom_x_ref, &
        trust_ncg_custom_fun_ref, trust_ncg_custom_nit_ref, trust_ncg_custom_nfev_ref, &
        trust_ncg_custom_njev_ref, trust_ncg_custom_status_ref, 5.0e-11_dp, failures)
    call check_int('trust-ncg custom nhev', r%nhev, trust_ncg_custom_nhev_ref, failures)

    r = trust_ncg(quadratic, [2.0_dp, -1.0_dp], quadratic_gradient, &
        hess=quadratic_hessian)
    call check_result('trust-ncg quadratic', r, trust_ncg_quad_hess_x_ref, &
        trust_ncg_quad_hess_fun_ref, trust_ncg_quad_hess_nit_ref, &
        trust_ncg_quad_hess_nfev_ref, trust_ncg_quad_hess_njev_ref, &
        trust_ncg_quad_hess_status_ref, 5.0e-12_dp, failures)
    call check_int('trust-ncg quadratic nhev', r%nhev, trust_ncg_quad_hess_nhev_ref, &
        failures)

    r = trust_ncg(indefinite_quadratic, [1.0_dp, 1.0_dp], &
        indefinite_quadratic_gradient, hess=indefinite_quadratic_hessian, maxiter=3)
    call check_result('trust-ncg indefinite', r, trust_ncg_indefinite_x_ref, &
        trust_ncg_indefinite_fun_ref, trust_ncg_indefinite_nit_ref, &
        trust_ncg_indefinite_nfev_ref, trust_ncg_indefinite_njev_ref, &
        trust_ncg_indefinite_status_ref, 5.0e-12_dp, failures)
    call check_int('trust-ncg indefinite nhev', r%nhev, trust_ncg_indefinite_nhev_ref, &
        failures)

    r = minimize(rosen, [-1.2_dp, 1.0_dp], method='trust-krylov', &
        jac=rosen_gradient, hess=rosen_hessian)
    call check_result('trust-krylov hess', r, trust_krylov_rosen_hess_x_ref, &
        trust_krylov_rosen_hess_fun_ref, trust_krylov_rosen_hess_nit_ref, &
        trust_krylov_rosen_hess_nfev_ref, trust_krylov_rosen_hess_njev_ref, &
        trust_krylov_rosen_hess_status_ref, 5.0e-10_dp, failures)
    call check_int('trust-krylov hess nhev', r%nhev, &
        trust_krylov_rosen_hess_nhev_ref, failures)

    r = trust_krylov(rosen, [-1.2_dp, 1.0_dp], rosen_gradient, &
        hessp=rosen_hess_product)
    call check_result('trust-krylov hessp', r, trust_krylov_rosen_hessp_x_ref, &
        trust_krylov_rosen_hessp_fun_ref, trust_krylov_rosen_hessp_nit_ref, &
        trust_krylov_rosen_hessp_nfev_ref, trust_krylov_rosen_hessp_njev_ref, &
        trust_krylov_rosen_hessp_status_ref, 5.0e-10_dp, failures)
    call check_int('trust-krylov hessp nhev', r%nhev, &
        trust_krylov_rosen_hessp_nhev_ref, failures)

    r = minimize(rosen, [-1.2_dp, 1.0_dp], method='trust-krylov', &
        jac=rosen_gradient, hessp=rosen_hess_product, inexact=.false.)
    call check_result('trust-krylov exact subproblem', r, trust_krylov_rosen_exact_x_ref, &
        trust_krylov_rosen_exact_fun_ref, trust_krylov_rosen_exact_nit_ref, &
        trust_krylov_rosen_exact_nfev_ref, trust_krylov_rosen_exact_njev_ref, &
        trust_krylov_rosen_exact_status_ref, 5.0e-10_dp, failures)
    call check_int('trust-krylov exact nhev', r%nhev, &
        trust_krylov_rosen_exact_nhev_ref, failures)

    r = trust_krylov(rosen, [-1.2_dp, 1.0_dp], rosen_gradient, &
        hess=rosen_hessian, maxiter=1)
    call check_result('trust-krylov maxiter', r, trust_krylov_limit_x_ref, &
        trust_krylov_limit_fun_ref, trust_krylov_limit_nit_ref, &
        trust_krylov_limit_nfev_ref, trust_krylov_limit_njev_ref, &
        trust_krylov_limit_status_ref, 5.0e-10_dp, failures)
    call check_int('trust-krylov maxiter nhev', r%nhev, &
        trust_krylov_limit_nhev_ref, failures)

    r = minimize(rosen, [-1.2_dp, 1.0_dp], method='trust_krylov', &
        jac=rosen_gradient, hessp=rosen_hess_product, initial_trust_radius=0.25_dp, &
        max_trust_radius=10.0_dp, eta=0.2_dp)
    call check_result('trust-krylov custom radii', r, trust_krylov_custom_x_ref, &
        trust_krylov_custom_fun_ref, trust_krylov_custom_nit_ref, &
        trust_krylov_custom_nfev_ref, trust_krylov_custom_njev_ref, &
        trust_krylov_custom_status_ref, 5.0e-10_dp, failures)
    call check_int('trust-krylov custom nhev', r%nhev, &
        trust_krylov_custom_nhev_ref, failures)

    r = trust_krylov(quadratic, [2.0_dp, -1.0_dp], quadratic_gradient, &
        hess=quadratic_hessian)
    call check_result('trust-krylov quadratic hess', r, trust_krylov_quad_hess_x_ref, &
        trust_krylov_quad_hess_fun_ref, trust_krylov_quad_hess_nit_ref, &
        trust_krylov_quad_hess_nfev_ref, trust_krylov_quad_hess_njev_ref, &
        trust_krylov_quad_hess_status_ref, 5.0e-11_dp, failures)
    call check_int('trust-krylov quadratic hess nhev', r%nhev, &
        trust_krylov_quad_hess_nhev_ref, failures)

    r = trust_krylov(quadratic, [2.0_dp, -1.0_dp], quadratic_gradient, &
        hessp=quadratic_hess_product)
    call check_result('trust-krylov quadratic hessp', r, trust_krylov_quad_hessp_x_ref, &
        trust_krylov_quad_hessp_fun_ref, trust_krylov_quad_hessp_nit_ref, &
        trust_krylov_quad_hessp_nfev_ref, trust_krylov_quad_hessp_njev_ref, &
        trust_krylov_quad_hessp_status_ref, 5.0e-11_dp, failures)
    call check_int('trust-krylov quadratic hessp nhev', r%nhev, &
        trust_krylov_quad_hessp_nhev_ref, failures)

    r = trust_krylov(rosen, [-1.2_dp, 1.0_dp, 0.5_dp, -0.5_dp, 1.2_dp], &
        rosen_gradient, hessp=rosen_hess_product)
    call check_result('trust-krylov rosen5', r, trust_krylov_rosen5_x_ref, &
        trust_krylov_rosen5_fun_ref, trust_krylov_rosen5_nit_ref, &
        trust_krylov_rosen5_nfev_ref, trust_krylov_rosen5_njev_ref, &
        trust_krylov_rosen5_status_ref, 2.0e-8_dp, failures)
    call check_int('trust-krylov rosen5 nhev', r%nhev, &
        trust_krylov_rosen5_nhev_ref, failures)

    r = minimize(rosen, [-1.2_dp, 1.0_dp], method='trust-exact', &
        jac=rosen_gradient, hess=rosen_hessian)
    call check_result('trust-exact hess', r, trust_exact_rosen_hess_x_ref, &
        trust_exact_rosen_hess_fun_ref, trust_exact_rosen_hess_nit_ref, &
        trust_exact_rosen_hess_nfev_ref, trust_exact_rosen_hess_njev_ref, &
        trust_exact_rosen_hess_status_ref, 5.0e-11_dp, failures)
    call check_int('trust-exact hess nhev', r%nhev, trust_exact_rosen_hess_nhev_ref, failures)

    r = trust_exact(rosen, [-1.2_dp, 1.0_dp], rosen_gradient, rosen_hessian, maxiter=1)
    call check_result('trust-exact maxiter', r, trust_exact_limit_x_ref, &
        trust_exact_limit_fun_ref, trust_exact_limit_nit_ref, trust_exact_limit_nfev_ref, &
        trust_exact_limit_njev_ref, trust_exact_limit_status_ref, 5.0e-12_dp, failures)
    call check_int('trust-exact maxiter nhev', r%nhev, trust_exact_limit_nhev_ref, failures)

    r = minimize(rosen, [-1.2_dp, 1.0_dp], method='trust_exact', &
        jac=rosen_gradient, hess=rosen_hessian, initial_trust_radius=0.25_dp, &
        max_trust_radius=10.0_dp, eta=0.2_dp, subproblem_maxiter=50)
    call check_result('trust-exact custom radii', r, trust_exact_custom_x_ref, &
        trust_exact_custom_fun_ref, trust_exact_custom_nit_ref, &
        trust_exact_custom_nfev_ref, trust_exact_custom_njev_ref, &
        trust_exact_custom_status_ref, 5.0e-11_dp, failures)
    call check_int('trust-exact custom nhev', r%nhev, trust_exact_custom_nhev_ref, failures)

    r = trust_exact(quadratic, [2.0_dp, -1.0_dp], quadratic_gradient, quadratic_hessian)
    call check_result('trust-exact quadratic', r, trust_exact_quad_hess_x_ref, &
        trust_exact_quad_hess_fun_ref, trust_exact_quad_hess_nit_ref, &
        trust_exact_quad_hess_nfev_ref, trust_exact_quad_hess_njev_ref, &
        trust_exact_quad_hess_status_ref, 5.0e-12_dp, failures)
    call check_int('trust-exact quadratic nhev', r%nhev, trust_exact_quad_hess_nhev_ref, failures)

    r = trust_exact(indefinite_quadratic, [1.0_dp, 1.0_dp], &
        indefinite_quadratic_gradient, indefinite_quadratic_hessian, maxiter=3)
    call check_result('trust-exact indefinite', r, trust_exact_indefinite_x_ref, &
        trust_exact_indefinite_fun_ref, trust_exact_indefinite_nit_ref, &
        trust_exact_indefinite_nfev_ref, trust_exact_indefinite_njev_ref, &
        trust_exact_indefinite_status_ref, 5.0e-12_dp, failures)
    call check_int('trust-exact indefinite nhev', r%nhev, &
        trust_exact_indefinite_nhev_ref, failures)

    r = trust_exact(trust_exact_hard_quadratic, [0.0_dp, 1.0_dp], &
        trust_exact_hard_quadratic_gradient, trust_exact_hard_quadratic_hessian, maxiter=1)
    call check_result('trust-exact hard case', r, trust_exact_hard_case_x_ref, &
        trust_exact_hard_case_fun_ref, trust_exact_hard_case_nit_ref, &
        trust_exact_hard_case_nfev_ref, trust_exact_hard_case_njev_ref, &
        trust_exact_hard_case_status_ref, 5.0e-12_dp, failures)
    call check_int('trust-exact hard case nhev', r%nhev, &
        trust_exact_hard_case_nhev_ref, failures)

    r = minimize(quadratic, [2.0_dp, -1.0_dp], jac=quadratic_gradient)
    call check_result('default minimize bfgs', r, default_quad_jac_x_ref, &
        default_quad_jac_fun_ref, default_quad_jac_nit_ref, default_quad_jac_nfev_ref, &
        default_quad_jac_njev_ref, default_quad_jac_status_ref, 5.0e-8_dp, failures)

    r = nelder_mead(rosen, [-1.2_dp, 1.0_dp, 0.5_dp], adaptive=.true.)
    call check_result('nelder-mead adaptive', r, nm_adaptive3_x_ref, nm_adaptive3_fun_ref, &
        nm_adaptive3_nit_ref, nm_adaptive3_nfev_ref, nm_adaptive3_njev_ref, &
        nm_adaptive3_status_ref, 5.0e-12_dp, failures)

    if (failures /= 0) then
        write(*, '(a,i0)') 'test_optimize_minimize: FAIL ', failures
        stop 1
    end if
    write(*, '(a)') 'test_optimize_minimize: PASS'

contains

    function box_quadratic(x) result(fx)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp) :: fx

        fx = (x(1) - 2.0_dp)**2 + (x(2) + 1.0_dp)**2
    end function box_quadratic

    function rosen(x) result(fx)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp) :: fx
        integer :: i

        fx = 0.0_dp
        do i = 1, size(x) - 1
            fx = fx + 100.0_dp * (x(i + 1) - x(i)**2)**2 + (1.0_dp - x(i))**2
        end do
    end function rosen

    subroutine rosen_gradient(x, grad)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: grad(:) !! Rosenbrock gradient
        integer :: i

        grad = 0.0_dp
        grad(1) = -400.0_dp * x(1) * (x(2) - x(1)**2) - 2.0_dp * (1.0_dp - x(1))
        do i = 2, size(x) - 1
            grad(i) = 200.0_dp * (x(i) - x(i - 1)**2) &
                - 400.0_dp * (x(i + 1) - x(i)**2) * x(i) &
                - 2.0_dp * (1.0_dp - x(i))
        end do
        grad(size(x)) = 200.0_dp * (x(size(x)) - x(size(x) - 1)**2)
    end subroutine rosen_gradient

    subroutine rosen_hessian(x, hess)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: hess(:, :) !! Rosenbrock Hessian
        integer :: i

        hess = 0.0_dp
        do i = 1, size(x)
            if (i == 1) then
                hess(i, i) = 1200.0_dp * x(i)**2 - 400.0_dp * x(i + 1) + 2.0_dp
            else if (i == size(x)) then
                hess(i, i) = 200.0_dp
            else
                hess(i, i) = 1200.0_dp * x(i)**2 - 400.0_dp * x(i + 1) + 202.0_dp
            end if
        end do
        do i = 1, size(x) - 1
            hess(i, i + 1) = -400.0_dp * x(i)
            hess(i + 1, i) = hess(i, i + 1)
        end do
    end subroutine rosen_hessian

    subroutine rosen_hess_product(x, p, hp)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(in) :: p(:) !! Hessian-product direction
        real(dp), intent(out) :: hp(:) !! Rosenbrock Hessian-vector product
        real(dp) :: hess(size(x), size(x))

        call rosen_hessian(x, hess)
        hp = matmul(hess, p)
    end subroutine rosen_hess_product

    function quadratic(x) result(fx)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp) :: fx

        fx = 0.5_dp * (3.0_dp * x(1)**2 + x(1) * x(2) + x(2)**2) &
            - x(1) + 2.0_dp * x(2)
    end function quadratic

    subroutine quadratic_gradient(x, grad)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: grad(:) !! quadratic gradient

        grad(1) = 3.0_dp * x(1) + 0.5_dp * x(2) - 1.0_dp
        grad(2) = 0.5_dp * x(1) + x(2) + 2.0_dp
    end subroutine quadratic_gradient

    subroutine quadratic_hessian(x, hess)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: hess(:, :) !! constant quadratic Hessian

        if (size(x) /= 2) error stop 'quadratic_hessian expects two parameters'
        hess = reshape([3.0_dp, 0.5_dp, 0.5_dp, 1.0_dp], [2, 2])
    end subroutine quadratic_hessian

    subroutine quadratic_hess_product(x, pvec, hp)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(in) :: pvec(:) !! Hessian-product direction
        real(dp), intent(out) :: hp(:) !! Hessian-vector product

        if (size(x) /= 2) error stop 'quadratic_hess_product expects two parameters'
        hp(1) = 3.0_dp * pvec(1) + 0.5_dp * pvec(2)
        hp(2) = 0.5_dp * pvec(1) + pvec(2)
    end subroutine quadratic_hess_product

    function indefinite_quadratic(x) result(fx)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp) :: fx

        fx = 0.5_dp * (x(1)**2 - x(2)**2)
    end function indefinite_quadratic

    subroutine indefinite_quadratic_gradient(x, grad)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: grad(:) !! indefinite quadratic gradient

        grad = [x(1), -x(2)]
    end subroutine indefinite_quadratic_gradient

    subroutine indefinite_quadratic_hessian(x, hess)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: hess(:, :) !! indefinite quadratic Hessian

        if (size(x) /= 2) error stop 'indefinite_quadratic_hessian expects two parameters'
        hess = 0.0_dp
        hess(1, 1) = 1.0_dp
        hess(2, 2) = -1.0_dp
    end subroutine indefinite_quadratic_hessian

    function trust_exact_hard_quadratic(x) result(fx)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp) :: fx

        fx = -0.5_dp * x(1)**2 + x(2)**2
    end function trust_exact_hard_quadratic

    subroutine trust_exact_hard_quadratic_gradient(x, grad)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: grad(:) !! hard-case quadratic gradient

        grad = [-x(1), 2.0_dp * x(2)]
    end subroutine trust_exact_hard_quadratic_gradient

    subroutine trust_exact_hard_quadratic_hessian(x, hess)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: hess(:, :) !! hard-case indefinite Hessian

        if (size(x) /= 2) error stop 'trust_exact_hard_quadratic_hessian expects two parameters'
        hess = 0.0_dp
        hess(1, 1) = -1.0_dp
        hess(2, 2) = 2.0_dp
    end subroutine trust_exact_hard_quadratic_hessian

    subroutine check_result(name, actual, expected_x, expected_fun, expected_nit, &
            expected_nfev, expected_njev, expected_status, tol, nfail)
        character(len=*), intent(in) :: name !! check label
        type(multivariate_minimize_result), intent(in) :: actual !! actual result
        real(dp), intent(in) :: expected_x(:) !! SciPy reference minimizer
        real(dp), intent(in) :: expected_fun !! SciPy reference objective
        integer, intent(in) :: expected_nit !! SciPy reference iterations
        integer, intent(in) :: expected_nfev !! SciPy reference function evaluations
        integer, intent(in) :: expected_njev !! SciPy reference gradient evaluations
        integer, intent(in) :: expected_status !! SciPy reference status
        real(dp), intent(in) :: tol !! absolute numerical tolerance
        integer, intent(inout) :: nfail !! accumulated failure count

        call check_vector(trim(name)//' x', actual%x, expected_x, tol, nfail)
        call check_real(trim(name)//' fun', actual%fun, expected_fun, tol, nfail)
        call check_int(trim(name)//' nit', actual%nit, expected_nit, nfail)
        if (expected_nfev >= 0) then
            call check_int(trim(name)//' nfev', actual%nfev, expected_nfev, nfail)
        end if
        if (expected_njev >= 0) then
            call check_int(trim(name)//' njev', actual%njev, expected_njev, nfail)
        end if
        call check_int(trim(name)//' status', actual%status, expected_status, nfail)
        call check_true(trim(name)//' success', actual%success .eqv. (expected_status == 0), &
            nfail)
    end subroutine check_result

    subroutine check_vector(name, actual, expected, tol, nfail)
        character(len=*), intent(in) :: name !! check label
        real(dp), intent(in) :: actual(:) !! actual vector
        real(dp), intent(in) :: expected(:) !! reference vector
        real(dp), intent(in) :: tol !! absolute tolerance
        integer, intent(inout) :: nfail !! accumulated failure count

        if (size(actual) /= size(expected)) then
            write(*, '(a)') trim(name)//' size mismatch'
            nfail = nfail + 1
        else if (maxval(abs(actual - expected)) > tol) then
            write(*, '(a,*(1x,es24.16))') trim(name)//' actual:', actual
            write(*, '(a,*(1x,es24.16))') trim(name)//' expected:', expected
            nfail = nfail + 1
        end if
    end subroutine check_vector

    subroutine check_real(name, actual, expected, tol, nfail)
        character(len=*), intent(in) :: name !! check label
        real(dp), intent(in) :: actual !! actual value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tol !! absolute tolerance
        integer, intent(inout) :: nfail !! accumulated failure count

        if (abs(actual - expected) > tol) then
            write(*, '(a,2(1x,es24.16))') trim(name)//' mismatch:', actual, expected
            nfail = nfail + 1
        end if
    end subroutine check_real

    subroutine check_int(name, actual, expected, nfail)
        character(len=*), intent(in) :: name !! check label
        integer, intent(in) :: actual !! actual integer
        integer, intent(in) :: expected !! reference integer
        integer, intent(inout) :: nfail !! accumulated failure count

        if (actual /= expected) then
            write(*, '(a,2(1x,i0))') trim(name)//' mismatch:', actual, expected
            nfail = nfail + 1
        end if
    end subroutine check_int

    subroutine check_true(name, condition, nfail)
        character(len=*), intent(in) :: name !! check label
        logical, intent(in) :: condition !! condition that must hold
        integer, intent(inout) :: nfail !! accumulated failure count

        if (.not. condition) then
            write(*, '(a)') trim(name)//' failed'
            nfail = nfail + 1
        end if
    end subroutine check_true

end program test_optimize_minimize
