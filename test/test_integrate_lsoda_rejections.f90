program test_integrate_lsoda_rejections
    use scifort_kinds, only : dp
    use scifort_integrate_lsoda, only : lsoda_advance, lsoda_initialize, lsoda_state, &
        lsoda_status_max_attempts, lsoda_status_success
    implicit none

    call test_rejected_step_order_choice()
    call test_second_failure_step_cap()

    print *, 'test_integrate_lsoda_rejections: PASS'

contains

    subroutine test_rejected_step_order_choice()
        type(lsoda_state) :: state
        real(dp) :: t
        real(dp) :: t_new
        real(dp) :: y(1)
        real(dp) :: y_new(1)
        real(dp) :: f(1)
        real(dp) :: f_new(1)
        real(dp) :: dense_d(1, 13)
        real(dp) :: dense_h
        real(dp) :: dense_anchor
        integer :: dense_order
        integer :: nfev
        integer :: njev
        integer :: nlu
        integer :: attempts
        integer :: rejected
        integer :: status
        integer :: old_order
        logical :: dense_valid

        call build_high_order_state(state, t, y, f, nfev, njev, nlu)
        old_order = state%adams_method_order
        call require(old_order >= 4, 'LSODA rejection test reaches high Adams order')

        state%adams_yh(:, old_order + 1) = 0.0_dp
        state%h_abs = 0.20_dp
        call lsoda_advance(state, exponential_rhs, t, y, f, 10.0_dp, 1.0_dp, &
            3.0e-13_dp, 1.0e-14_dp, 0.20_dp, y_new, f_new, t_new, dense_d, &
            dense_h, dense_order, dense_anchor, dense_valid, nfev, njev, nlu, &
            1, attempts, rejected, status)

        call require(status == lsoda_status_max_attempts, &
            'LSODA rejected-order probe exhausts one attempt')
        call require(attempts == 1 .and. rejected == 1, &
            'LSODA rejected-order probe records one error failure')
        call require(state%adams_method_order == old_order - 1, &
            'LSODA rejected step chooses the lower Adams order')
        call require(state%h_abs > 0.19_dp .and. state%h_abs <= 0.20_dp, &
            'LSODA first rejection may retain nearly the same step at lower order')
        call require(state%adams_ialth == state%adams_method_order + 1, &
            'LSODA rejected-step order change resets IALTH')
    end subroutine test_rejected_step_order_choice

    subroutine test_second_failure_step_cap()
        type(lsoda_state) :: state
        real(dp) :: t
        real(dp) :: t_new
        real(dp) :: y(1)
        real(dp) :: y_new(1)
        real(dp) :: f(1)
        real(dp) :: f_new(1)
        real(dp) :: dense_d(1, 13)
        real(dp) :: dense_h
        real(dp) :: dense_anchor
        integer :: dense_order
        integer :: nfev
        integer :: njev
        integer :: nlu
        integer :: attempts
        integer :: rejected
        integer :: status
        integer :: old_order
        logical :: dense_valid

        call build_high_order_state(state, t, y, f, nfev, njev, nlu)
        old_order = state%adams_method_order
        state%adams_yh(:, old_order + 1) = 0.0_dp
        state%h_abs = 0.20_dp
        call lsoda_advance(state, exponential_rhs, t, y, f, 10.0_dp, 1.0_dp, &
            3.0e-13_dp, 1.0e-14_dp, 0.20_dp, y_new, f_new, t_new, dense_d, &
            dense_h, dense_order, dense_anchor, dense_valid, nfev, njev, nlu, &
            2, attempts, rejected, status)

        call require(status == lsoda_status_max_attempts, &
            'LSODA repeated-rejection probe exhausts two attempts')
        call require(attempts == 2 .and. rejected == 2, &
            'LSODA repeated-rejection probe records two error failures')
        call require(state%h_abs > 0.0_dp .and. state%h_abs <= 0.0400001_dp, &
            'LSODA second error failure caps retry factor at 0.2')
    end subroutine test_second_failure_step_cap

    subroutine build_high_order_state(state, t, y, f, nfev, njev, nlu)
        type(lsoda_state), intent(out) :: state !! built nonstiff LSODA state
        real(dp), intent(out) :: t !! endpoint after high-order Adams buildup
        real(dp), intent(out) :: y(1) !! solution value at the returned endpoint
        real(dp), intent(out) :: f(1) !! derivative value at the returned endpoint
        integer, intent(out) :: nfev !! accumulated right-hand-side evaluations
        integer, intent(out) :: njev !! accumulated Jacobian evaluations
        integer, intent(out) :: nlu !! accumulated LU factorizations

        real(dp) :: t_new
        real(dp) :: y_new(1)
        real(dp) :: f_new(1)
        real(dp) :: dense_d(1, 13)
        real(dp) :: dense_h
        real(dp) :: dense_anchor
        integer :: dense_order
        integer :: attempts
        integer :: rejected
        integer :: status
        integer :: k
        logical :: dense_valid

        t = 0.0_dp
        y = 1.0_dp
        f = 1.0_dp
        nfev = 1
        njev = 0
        nlu = 0
        call lsoda_initialize(state, y, f, 0.001_dp)
        do k = 1, 60
            call lsoda_advance(state, exponential_rhs, t, y, f, 10.0_dp, 1.0_dp, &
                1.0e-9_dp, 1.0e-12_dp, 0.02_dp, y_new, f_new, t_new, dense_d, &
                dense_h, dense_order, dense_anchor, dense_valid, nfev, njev, nlu, &
                100, attempts, rejected, status)
            call require(status == lsoda_status_success, &
                'LSODA rejected-step test buildup succeeds')
            t = t_new
            y = y_new
            f = f_new
        end do
    end subroutine build_high_order_state

    subroutine exponential_rhs(t, y, dydt)
        real(dp), intent(in) :: t !! independent variable, unused for autonomous system
        real(dp), intent(in) :: y(:) !! current solution state
        real(dp), intent(out) :: dydt(:) !! derivative dy/dt

        dydt = y
        if (t < -huge(t)) dydt = 0.0_dp
    end subroutine exponential_rhs

    subroutine require(condition, label)
        logical, intent(in) :: condition !! assertion condition
        character(len=*), intent(in) :: label !! assertion label printed on failure

        if (.not. condition) then
            print *, 'FAIL ', trim(label)
            error stop 1
        end if
    end subroutine require

end program test_integrate_lsoda_rejections
