! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_integrate
    use scifort_integrate_bvp, only : solve_bvp, solve_bvp_result, solve_bvp_complex_result, &
        bvp_solution, bvp_complex_solution, bvp_status_success, bvp_status_max_nodes, &
        bvp_status_singular, bvp_status_bc_tolerance, bvp_status_invalid_argument, &
        bvp_status_numerical_failure, bvp_status_setup_failure, bvp_fun_real, bvp_bc_real, &
        bvp_fun_complex, bvp_bc_complex, bvp_fun_jac_real, bvp_bc_jac_real, &
        bvp_fun_jac_complex, bvp_bc_jac_complex
    use scifort_kinds, only : dp
    use scifort_integrate_quadrature, only : cumulative_simpson, cumulative_trapezoid, &
        fixed_quad, integrate_status_invalid_argument, integrate_status_invalid_shape, &
        integrate_status_numerical_failure, integrate_status_success, newton_cotes, &
        romb, simpson, trapezoid
    use scifort_integrate_cubature, only : cubature, cubature_status_max_subdivisions
    use scifort_integrate_ivp, only : ivp_dense_solution, ivp_event_set, ivp_status_event, &
        ivp_status_failure, ivp_status_invalid_argument, ivp_status_max_steps, &
        ivp_status_success, solve_ivp_real => solve_ivp, solve_ivp_result
    use scifort_integrate_ivp_complex, only : ivp_complex_dense_solution, &
        ivp_complex_event_set, solve_ivp_complex, solve_ivp_complex_result
    use scifort_integrate_nested, only : dblquad, nquad, tplquad
    use scifort_integrate_odeint, only : odeint, odeint_result, odeint_tfirst
    use scifort_integrate_sparse_jacobian, only : ivp_complex_sparse_matrix, &
        ivp_csc_matrix, ivp_sparse_matrix
    use scifort_integrate_qmc_quad, only : qmc_quad
    use scifort_integrate_quad_vec, only : quad_vec
    use scifort_integrate_tanhsinh, only : tanhsinh, tanhsinh_status_maxlevel
    use scifort_integrate_quadpack, only : quad, quad_status_bad_integrand, &
        quad_status_cycle_failure, quad_status_divergent, quad_status_invalid_argument, &
        quad_status_max_subdivisions, &
        quad_status_nonconvergence, quad_status_roundoff, quad_status_success
    implicit none
    private

    public :: solve_bvp, solve_bvp_result, solve_bvp_complex_result
    public :: bvp_solution, bvp_complex_solution, bvp_status_success, bvp_status_max_nodes
    public :: bvp_status_singular, bvp_status_bc_tolerance, bvp_status_invalid_argument
    public :: bvp_status_numerical_failure, bvp_status_setup_failure
    public :: bvp_fun_real, bvp_bc_real, bvp_fun_complex, bvp_bc_complex
    public :: bvp_fun_jac_real, bvp_bc_jac_real, bvp_fun_jac_complex, bvp_bc_jac_complex

    interface solve_ivp
        procedure :: solve_ivp_real
        procedure :: solve_ivp_complex
    end interface solve_ivp

    public :: cumulative_simpson
    public :: cumulative_trapezoid
    public :: cubature
    public :: cubature_status_max_subdivisions
    public :: dblquad
    public :: dp
    public :: fixed_quad
    public :: integrate_status_invalid_argument
    public :: integrate_status_invalid_shape
    public :: integrate_status_numerical_failure
    public :: integrate_status_success
    public :: ivp_dense_solution
    public :: ivp_complex_dense_solution
    public :: ivp_complex_event_set
    public :: ivp_event_set
    public :: ivp_sparse_matrix
    public :: ivp_complex_sparse_matrix
    public :: ivp_csc_matrix
    public :: ivp_status_event
    public :: ivp_status_failure
    public :: ivp_status_invalid_argument
    public :: ivp_status_max_steps
    public :: ivp_status_success
    public :: newton_cotes
    public :: odeint
    public :: odeint_result
    public :: odeint_tfirst
    public :: nquad
    public :: qmc_quad
    public :: quad
    public :: quad_vec
    public :: quad_status_bad_integrand
    public :: quad_status_cycle_failure
    public :: quad_status_divergent
    public :: quad_status_invalid_argument
    public :: quad_status_max_subdivisions
    public :: quad_status_nonconvergence
    public :: quad_status_roundoff
    public :: quad_status_success
    public :: romb
    public :: simpson
    public :: solve_ivp
    public :: solve_ivp_result
    public :: solve_ivp_complex_result
    public :: tanhsinh
    public :: tanhsinh_status_maxlevel
    public :: tplquad
    public :: trapezoid

end module scifort_integrate
