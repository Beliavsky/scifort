program test_integrate_lsoda_stiff_switch
    use scifort_integrate_ivp, only : solve_ivp, solve_ivp_result
    use scifort_integrate_lsoda, only : lsoda_advance, lsoda_initialize, lsoda_state, &
        lsoda_status_success
    use scifort_kinds, only : dp
    implicit none

    call test_deferred_stiff_handoff()
    call test_stiff_work_count_parity()
    call test_numerical_jacobian_accounting()

    print *, 'test_integrate_lsoda_stiff_switch: PASS'

contains

    subroutine test_deferred_stiff_handoff()
        type(lsoda_state) :: state
        real(dp) :: t
        real(dp) :: t_new
        real(dp) :: y(1)
        real(dp) :: y_new(1)
        real(dp) :: f(1)
        real(dp) :: f_new(1)
        real(dp) :: dense_d(1, 7)
        real(dp) :: dense_h
        real(dp) :: dense_anchor
        integer :: dense_order
        integer :: nfev
        integer :: njev
        integer :: nlu
        integer :: attempts
        integer :: rejected
        integer :: status
        integer :: k
        logical :: dense_valid

        t = 0.0_dp
        y = 1.0_dp
        call stiff_cosine_rhs(t, y, f)
        nfev = 1
        njev = 0
        nlu = 0
        call lsoda_initialize(state, y, f, 1.0e-3_dp)

        do k = 1, 21
            call lsoda_advance(state, stiff_cosine_rhs, t, y, f, 1.0_dp, 1.0_dp, &
                1.0e-6_dp, 1.0e-9_dp, huge(1.0_dp), y_new, f_new, t_new, dense_d, &
                dense_h, dense_order, dense_anchor, dense_valid, nfev, njev, nlu, &
                100, attempts, rejected, status, stiff_cosine_jacobian)
            call require(status == lsoda_status_success, 'LSODA stiff handoff status')
            t = t_new
            y = y_new
            f = f_new
        end do

        call check_abs(t, 1.1199980939832314e-2_dp, 5.0e-8_dp, &
            'LSODA stiff handoff accepted time')
        call require(state%bdf%order == 2, 'LSODA stiff handoff preserves BDF order')
        call require(njev == 0, 'LSODA stiff handoff defers Jacobian evaluation')
        call require(nlu == 0, 'LSODA stiff handoff defers LU factorization')
        call check_abs(state%bdf%h_abs, 5.7499188792070820e-3_dp, 2.0e-7_dp, &
            'LSODA stiff handoff selected BDF step')

        call lsoda_advance(state, stiff_cosine_rhs, t, y, f, 1.0_dp, 1.0_dp, &
            1.0e-6_dp, 1.0e-9_dp, huge(1.0_dp), y_new, f_new, t_new, dense_d, &
            dense_h, dense_order, dense_anchor, dense_valid, nfev, njev, nlu, &
            100, attempts, rejected, status, stiff_cosine_jacobian)
        call require(status == lsoda_status_success, 'LSODA first BDF trial status')
        call check_abs(t_new, 1.6949899819039396e-2_dp, 2.0e-7_dp, &
            'LSODA first BDF accepted time')
        call require(njev == 1, 'LSODA first BDF trial evaluates Jacobian')
        call require(nlu == 1, 'LSODA first BDF trial factors Newton matrix')
    end subroutine test_deferred_stiff_handoff

    subroutine test_stiff_work_count_parity()
        type(solve_ivp_result) :: sol

        sol = solve_ivp(stiff_cosine_rhs, [0.0_dp, 1.0_dp], [1.0_dp], &
            method='LSODA', jac=stiff_cosine_jacobian, rtol=1.0e-6_dp, atol=1.0e-9_dp)

        call require(sol%success, 'LSODA stiff reference trajectory success')
        call require(sol%nfev == 116, 'LSODA stiff reference RHS work count')
        call require(sol%njev == 10, 'LSODA stiff reference Jacobian work count')
        call require(sol%nlu == 10, 'LSODA stiff reference LU work count')
        call require(size(sol%t) == 59, 'LSODA stiff reference accepted-point count')
        call check_abs(sol%t(33), 2.4963393555355357e-1_dp, 5.0e-7_dp, &
            'LSODA stiff order-three transition time')
        call check_abs(sol%t(41), 5.2121319954882398e-1_dp, 5.0e-7_dp, &
            'LSODA stiff order-four transition time')
        call check_abs(sol%t(58), 9.8592731526334264e-1_dp, 5.0e-7_dp, &
            'LSODA stiff order-five transition time')
    end subroutine test_stiff_work_count_parity

    subroutine test_numerical_jacobian_accounting()
        type(solve_ivp_result) :: analytic
        type(solve_ivp_result) :: numerical

        analytic = solve_ivp(stiff_cosine_rhs, [0.0_dp, 1.0_dp], [1.0_dp], &
            method='LSODA', jac=stiff_cosine_jacobian, rtol=1.0e-6_dp, atol=1.0e-9_dp)
        numerical = solve_ivp(stiff_cosine_rhs, [0.0_dp, 1.0_dp], [1.0_dp], &
            method='LSODA', rtol=1.0e-6_dp, atol=1.0e-9_dp)

        call require(analytic%success, 'LSODA stiff analytic Jacobian success')
        call require(numerical%success, 'LSODA stiff numerical Jacobian success')
        call require(analytic%njev == numerical%njev, &
            'LSODA stiff analytic/numerical Jacobian evaluation parity')
        call require(analytic%nlu == numerical%nlu, &
            'LSODA stiff analytic/numerical LU parity')
        call require(analytic%nfev == 116 .and. numerical%nfev == 126, &
            'LSODA stiff analytic/numerical RHS reference counts')
        call require(analytic%njev == 10 .and. analytic%nlu == 10, &
            'LSODA stiff analytic Jacobian/LU reference counts')
        call require(numerical%nfev - analytic%nfev == numerical%njev, &
            'LSODA scalar finite-difference Jacobian RHS accounting')
        call check_abs(numerical%y(1, size(numerical%t)), &
            analytic%y(1, size(analytic%t)), 2.0e-12_dp, &
            'LSODA analytic/numerical stiff endpoint parity')
    end subroutine test_numerical_jacobian_accounting

    subroutine stiff_cosine_rhs(t, y, dydt)
        real(dp), intent(in) :: t !! independent variable
        real(dp), intent(in) :: y(:) !! state at t
        real(dp), intent(out) :: dydt(:) !! derivative at t

        dydt(1) = -1000.0_dp * (y(1) - cos(t)) - sin(t)
    end subroutine stiff_cosine_rhs

    subroutine stiff_cosine_jacobian(t, y, jac)
        real(dp), intent(in) :: t !! independent variable
        real(dp), intent(in) :: y(:) !! state at t
        real(dp), intent(out) :: jac(:, :) !! dense Jacobian at t and y

        if (t /= t .or. size(y) < 1) error stop 'invalid Jacobian input'
        jac = 0.0_dp
        jac(1, 1) = -1000.0_dp
    end subroutine stiff_cosine_jacobian

    subroutine require(condition, label)
        logical, intent(in) :: condition !! assertion condition
        character(len=*), intent(in) :: label !! diagnostic label

        if (.not. condition) then
            print *, 'FAIL: ', trim(label)
            error stop 1
        end if
    end subroutine require

    subroutine check_abs(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tolerance !! allowed absolute error
        character(len=*), intent(in) :: label !! diagnostic label

        if (abs(actual - expected) > tolerance) then
            print *, 'FAIL: ', trim(label)
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance= ', tolerance
            error stop 1
        end if
    end subroutine check_abs

end program test_integrate_lsoda_stiff_switch
