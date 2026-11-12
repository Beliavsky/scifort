! SPDX-License-Identifier: MIT
! Work-aware refresh on a nonlinear system with an evolving rank-two Jacobian.
program test_optimize_root_krylov_work_refresh_nonlinear
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, vector_root_result
    implicit none

    integer, parameter :: n = 32
    real(dp), parameter :: rate = 2.0_dp
    real(dp) :: a(n), u(n,2), v(n,2), target(n), x0(n), target_value(n)
    real(dp) :: theta, error
    integer :: i, base_calls, jv_calls, jtv_calls
    integer :: base_stall, base_cheap, base_cost4, base_expensive
    integer :: jv_stall, jv_cheap
    type(vector_root_result) :: fit

    do i = 1, n
        theta = real(i, dp)
        a(i) = 10.0_dp**(-2.0_dp*real(i-1, dp)/real(n-1, dp))
        u(i,1) = 2.0_dp*a(i)*(0.8_dp + 0.22_dp*sin(0.53_dp*theta))
        u(i,2) = 2.0_dp*a(i)*(-0.6_dp + 0.17_dp*cos(0.61_dp*theta))
        v(i,1) = (0.9_dp + 0.30_dp*cos(0.29_dp*theta))/real(n, dp)
        v(i,2) = (0.7_dp + 0.42_dp*sin(0.37_dp*theta))/real(n, dp)
        target(i) = 0.08_dp*cos(0.31_dp*theta)
        x0(i) = target(i) + 0.95_dp + 0.25_dp*sin(0.23_dp*theta)
    end do
    call nonlinear_action(target, target_value)

    call run_case(1)
    base_stall = base_calls
    jv_stall = jv_calls
    call run_case(2)
    base_cheap = base_calls
    jv_cheap = jv_calls
    call run_case(3)
    base_cost4 = base_calls
    call run_case(4)
    base_expensive = base_calls

    if (base_stall /= 1) error stop 'unexpected stall-only refresh'
    if (base_cheap <= base_stall) error stop 'inner work did not trigger early refresh'
    if (base_cost4 >= base_cheap) error stop 'rebuild cost did not suppress a refresh'
    if (base_expensive /= base_stall) error stop 'expensive rebuild triggered unexpectedly'
    if (jv_cheap >= jv_stall) error stop 'work-aware refresh did not reduce forward products'
    print '(a)', 'nonlinear work-aware Krylov refresh tests passed'
contains
    subroutine run_case(mode)
        integer, intent(in) :: mode !! 1 stall, 2 cheap, 3 moderate, 4 expensive

        base_calls = 0
        jv_calls = 0
        jtv_calls = 0
        select case (mode)
        case (1)
            fit = root(residual, x0, method='krylov', &
                inner_method='gmres', inner_maxiter=60, maxiter=50, &
                fatol=1.0e-10_dp, tol_norm=scaled_norm, &
                inner_m_diag_auto=.true., inner_m_diag_update_interval=12, &
                inner_m_diag_update_adaptive=.true., &
                inner_m_diag_update_stall_ratio=0.8_dp, &
                inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
                inner_m_diag_rank_two_nonsymmetric_base_fun=base_diagonal, &
                inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=16.0_dp, &
                jacobian_product=jvp, jacobian_transpose_product=jtv)
        case default
            fit = root(residual, x0, method='krylov', &
                inner_method='gmres', inner_maxiter=60, maxiter=50, &
                fatol=1.0e-10_dp, tol_norm=scaled_norm, &
                inner_m_diag_auto=.true., inner_m_diag_update_interval=12, &
                inner_m_diag_update_adaptive=.true., &
                inner_m_diag_update_stall_ratio=0.8_dp, &
                inner_m_diag_update_work_growth=1.2_dp, &
                inner_m_diag_update_rebuild_cost=cost(mode), &
                inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
                inner_m_diag_rank_two_nonsymmetric_base_fun=base_diagonal, &
                inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=16.0_dp, &
                jacobian_product=jvp, jacobian_transpose_product=jtv)
        end select
        error = maxval(abs(fit%x-target))
        if (.not. fit%success) error stop 'nonlinear root did not converge'
        if (scaled_norm(fit%fun) > 1.0e-9_dp) error stop 'inaccurate residual'
        if (error > 1.0e-7_dp) error stop 'incorrect nonlinear solution'
        if (base_calls < 1 .or. jv_calls < 1 .or. jtv_calls < 1) then
            error stop 'missing Jacobian or refresh callbacks'
        end if
    end subroutine run_case

    integer function cost(mode) result(value)
        integer, intent(in) :: mode !! benchmark mode number
        value = 0
        if (mode == 3) value = 4
        if (mode == 4) value = 128
    end function cost

    subroutine nonlinear_action(x, value)
        real(dp), intent(in) :: x(:) !! current state
        real(dp), intent(out) :: value(:) !! model residual before target subtraction
        real(dp) :: s1, s2
        s1 = dot_product(v(:,1), x)
        s2 = dot_product(v(:,2), x)
        value = a*(x + 0.5_dp*rate*x**3) + &
            u(:,1)*sin(rate*s1)/rate + &
            u(:,2)*(s2 + 0.25_dp*rate*s2**3)
    end subroutine nonlinear_action

    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! current nonlinear state
        real(dp), intent(out) :: fvec(:) !! residual minus known target
        call nonlinear_action(x, fvec)
        fvec = fvec - target_value
    end subroutine residual

    subroutine base_diagonal(x, diagonal, ok)
        real(dp), intent(in) :: x(:) !! nonlinear state at preconditioner rebuild
        real(dp), intent(out) :: diagonal(:) !! known nonlinear diagonal D(x)
        logical, intent(out) :: ok !! valid diagonal flag
        base_calls = base_calls + 1
        diagonal = a*(1.0_dp + 1.5_dp*rate*x*x)
        ok = all(ieee_is_finite(diagonal))
    end subroutine base_diagonal

    subroutine jvp(x, direction, product, ok)
        real(dp), intent(in) :: x(:) !! Jacobian evaluation point
        real(dp), intent(in) :: direction(:) !! forward direction
        real(dp), intent(out) :: product(:) !! J(x) times direction
        logical, intent(out) :: ok !! valid forward product flag
        real(dp) :: s1, s2
        jv_calls = jv_calls + 1
        s1 = dot_product(v(:,1), x)
        s2 = dot_product(v(:,2), x)
        product = a*(1.0_dp + 1.5_dp*rate*x*x)*direction + &
            u(:,1)*cos(rate*s1)*dot_product(v(:,1), direction) + &
            u(:,2)*(1.0_dp + 0.75_dp*rate*s2*s2)*dot_product(v(:,2), direction)
        ok = all(ieee_is_finite(product))
    end subroutine jvp

    subroutine jtv(x, direction, product, ok)
        real(dp), intent(in) :: x(:) !! Jacobian evaluation point
        real(dp), intent(in) :: direction(:) !! transpose direction
        real(dp), intent(out) :: product(:) !! J(x)^T times direction
        logical, intent(out) :: ok !! valid transpose product flag
        real(dp) :: s1, s2
        jtv_calls = jtv_calls + 1
        s1 = dot_product(v(:,1), x)
        s2 = dot_product(v(:,2), x)
        product = a*(1.0_dp + 1.5_dp*rate*x*x)*direction + &
            v(:,1)*cos(rate*s1)*dot_product(u(:,1), direction) + &
            v(:,2)*(1.0_dp + 0.75_dp*rate*s2*s2)*dot_product(u(:,2), direction)
        ok = all(ieee_is_finite(product))
    end subroutine jtv

    function scaled_norm(fvec) result(value)
        real(dp), intent(in) :: fvec(:) !! residual to scale by nominal diagonal
        real(dp) :: value !! maximum scaled residual
        value = maxval(abs(fvec)/a)
    end function scaled_norm
end program test_optimize_root_krylov_work_refresh_nonlinear
