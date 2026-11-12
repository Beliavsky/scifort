program test_integrate_lsoda_transfer
    use scifort_kinds, only : dp
    use scifort_integrate_bdf, only : bdf_advance, bdf_dense_state, bdf_initialize, &
        bdf_state, bdf_status_success
    use scifort_integrate_lsoda, only : lsoda_adams_dense_state, lsoda_advance, &
        lsoda_initialize, lsoda_state, lsoda_status_success
    implicit none

    call test_bdf_to_adams_state_transfer()
    call test_odepack_jacobian_refresh()

    print *, 'test_integrate_lsoda_transfer: PASS'

contains

    subroutine test_bdf_to_adams_state_transfer()
        type(lsoda_state) :: state
        real(dp) :: t
        real(dp) :: t_new
        real(dp) :: y(1)
        real(dp) :: y_new(1)
        real(dp) :: f(1)
        real(dp) :: f_new(1)
        real(dp) :: y_bdf_poly(1)
        real(dp) :: y_adams_poly(1)
        real(dp) :: dense_d(1, 13)
        real(dp) :: dense_h
        real(dp) :: dense_anchor
        real(dp) :: theta
        real(dp) :: tq
        real(dp) :: bdf_linear
        integer :: dense_order
        integer :: nfev
        integer :: njev
        integer :: nlu
        integer :: njev_before
        integer :: nlu_before
        integer :: attempts
        integer :: rejected
        integer :: status
        integer :: k
        integer :: switch_order
        integer :: j
        logical :: dense_valid

        t = 0.0_dp
        y = 1.0_dp
        f = 1.0_dp
        nfev = 1
        njev = 0
        nlu = 0
        call lsoda_initialize(state, y, f, 1.0e-3_dp)
        call bdf_initialize(state%bdf, exponential_rhs, t, y, f, 1.0_dp, &
            1.0e-3_dp, 1.0e-7_dp, 1.0e-10_dp, njev, exponential_jacobian)
        state%bdf%odepack_jacobian_refresh = .true.

        do k = 1, 20
            call bdf_advance(state%bdf, exponential_rhs, t, 1.0_dp, 1.0_dp, &
                1.0e-7_dp, 1.0e-10_dp, 1.0e-2_dp, y_new, t_new, dense_d, dense_h, &
                dense_order, dense_anchor, nfev, njev, nlu, 100, attempts, rejected, &
                status, exponential_jacobian, f_new_out=f_new)
            call require(status == bdf_status_success, &
                'LSODA transfer BDF buildup succeeds')
            t = t_new
            y = y_new
            f = f_new
        end do
        call require(state%bdf%last_step_order >= 3, &
            'LSODA transfer BDF buildup reaches multistep order')

        ! The public state component uses 1 for Adams and 2 for BDF internally.
        state%mode = 2
        state%switch_countdown = -1
        state%h_abs = state%bdf%h_abs
        call lsoda_advance(state, exponential_rhs, t, y, f, 1.0_dp, 1.0_dp, &
            1.0e-7_dp, 1.0e-10_dp, 1.0e-2_dp, y_new, f_new, t_new, dense_d, dense_h, &
            dense_order, dense_anchor, dense_valid, nfev, njev, nlu, 100, attempts, &
            rejected, status, exponential_jacobian)
        call require(status == lsoda_status_success, 'LSODA BDF-to-Adams switch succeeds')
        switch_order = min(state%bdf%last_step_order, 5)
        call require(state%mode == 1, 'LSODA BDF-to-Adams switch enters Adams mode')
        call require(state%adams_yh_valid, 'LSODA BDF-to-Adams transfer keeps Nordsieck state')
        call require(state%adams_method_order == switch_order, &
            'LSODA BDF-to-Adams transfer preserves switch order')
        call require(.not. state%force_history_build, &
            'LSODA BDF-to-Adams transfer avoids low-order history rebuild')
        call check_rel(state%adams_yh(1, 1), y_new(1), 4.0e-15_dp, &
            'LSODA BDF-to-Adams transferred endpoint')
        call check_rel(state%adams_yh(1, 2), state%adams_yh_h * f_new(1), 4.0e-15_dp, &
            'LSODA BDF-to-Adams transferred derivative')

        ! The converted higher coefficients reproduce the same truncated BDF
        ! polynomial, apart from replacing its linear term by the exact h*f.
        theta = -0.37_dp
        tq = t_new + theta * state%adams_yh_h
        call bdf_dense_state(dense_anchor, dense_h, switch_order, dense_d, tq, y_bdf_poly)
        call lsoda_adams_dense_state(t_new, state%adams_yh_h, switch_order, &
            state%adams_yh, tq, y_adams_poly)
        bdf_linear = 0.0_dp
        do j = 1, switch_order
            bdf_linear = bdf_linear + dense_d(1, j + 1) / real(j, dp)
        end do
        call check_rel(y_adams_poly(1), y_bdf_poly(1) + theta * &
            (state%adams_yh_h * f_new(1) - bdf_linear), 2.0e-14_dp, &
            'LSODA BDF-to-Adams higher polynomial transfer')

        t = t_new
        y = y_new
        f = f_new
        njev_before = njev
        nlu_before = nlu
        call lsoda_advance(state, exponential_rhs, t, y, f, 1.0_dp, 1.0_dp, &
            1.0e-7_dp, 1.0e-10_dp, 1.0e-2_dp, y_new, f_new, t_new, dense_d, dense_h, &
            dense_order, dense_anchor, dense_valid, nfev, njev, nlu, 100, attempts, &
            rejected, status, exponential_jacobian)
        call require(status == lsoda_status_success, &
            'LSODA transferred Adams state advances successfully')
        call require(njev == njev_before .and. nlu == nlu_before, &
            'LSODA transferred Adams step avoids stiff work')
        call require(state%adams_method_order == switch_order, &
            'LSODA transferred Adams step retains delayed switch order')
    end subroutine test_bdf_to_adams_state_transfer

    subroutine test_odepack_jacobian_refresh()
        type(bdf_state) :: state
        real(dp) :: t
        real(dp) :: t_new
        real(dp) :: y(1)
        real(dp) :: y_new(1)
        real(dp) :: f(1)
        real(dp) :: f_new(1)
        real(dp) :: dense_d(1, 8)
        real(dp) :: dense_h
        real(dp) :: dense_anchor
        integer :: dense_order
        integer :: nfev
        integer :: njev
        integer :: nlu
        integer :: njev_before
        integer :: attempts
        integer :: rejected
        integer :: status
        integer :: k

        t = 0.0_dp
        y = 1.0_dp
        f = -1.0_dp
        nfev = 1
        njev = 0
        nlu = 0
        call bdf_initialize(state, decay_rhs, t, y, f, 1.0_dp, 1.0e-4_dp, &
            1.0e-6_dp, 1.0e-9_dp, njev, decay_jacobian)
        state%odepack_jacobian_refresh = .true.

        do k = 1, 20
            call bdf_advance(state, decay_rhs, t, 1.0_dp, 1.0_dp, 1.0e-6_dp, &
                1.0e-9_dp, 1.0e-4_dp, y_new, t_new, dense_d, dense_h, dense_order, &
                dense_anchor, nfev, njev, nlu, 100, attempts, rejected, status, &
                decay_jacobian, f_new_out=f_new)
            call require(status == bdf_status_success, &
                'LSODA Jacobian-refresh interval buildup succeeds')
            t = t_new
            y = y_new
            f = f_new
            state%n_equal_steps = 0
        end do
        call require(njev == 1, 'LSODA BDF reuses Jacobian for first 20 steps')

        call bdf_advance(state, decay_rhs, t, 1.0_dp, 1.0_dp, 1.0e-6_dp, &
            1.0e-9_dp, 1.0e-4_dp, y_new, t_new, dense_d, dense_h, dense_order, &
            dense_anchor, nfev, njev, nlu, 100, attempts, rejected, status, &
            decay_jacobian, f_new_out=f_new)
        call require(status == bdf_status_success, &
            'LSODA periodic Jacobian-refresh step succeeds')
        call require(njev == 2, 'LSODA BDF refreshes Jacobian after 20 steps')

        t = t_new
        y = y_new
        f = f_new
        state%n_equal_steps = 0
        njev_before = njev
        call bdf_advance(state, decay_rhs, t, 1.0_dp, 1.0_dp, 1.0e-6_dp, &
            1.0e-9_dp, 4.0e-5_dp, y_new, t_new, dense_d, dense_h, dense_order, &
            dense_anchor, nfev, njev, nlu, 100, attempts, rejected, status, &
            decay_jacobian, f_new_out=f_new)
        call require(status == bdf_status_success, &
            'LSODA coefficient-change Jacobian-refresh step succeeds')
        call require(njev == njev_before + 1, &
            'LSODA BDF refreshes Jacobian after large coefficient change')
    end subroutine test_odepack_jacobian_refresh

    subroutine exponential_rhs(t, y, dydt)
        real(dp), intent(in) :: t !! independent variable, unused for autonomous system
        real(dp), intent(in) :: y(:) !! current exponential state
        real(dp), intent(out) :: dydt(:) !! derivative dy/dt

        dydt = y
        if (t < -huge(t)) dydt = 0.0_dp
    end subroutine exponential_rhs

    subroutine exponential_jacobian(t, y, jac)
        real(dp), intent(in) :: t !! independent variable, unused for autonomous system
        real(dp), intent(in) :: y(:) !! current exponential state, used for shape validation
        real(dp), intent(out) :: jac(:, :) !! scalar exponential Jacobian

        jac = 1.0_dp
        if (t < -huge(t) .or. size(y) < 1) jac = 0.0_dp
    end subroutine exponential_jacobian

    subroutine decay_rhs(t, y, dydt)
        real(dp), intent(in) :: t !! independent variable, unused for autonomous system
        real(dp), intent(in) :: y(:) !! current decay state
        real(dp), intent(out) :: dydt(:) !! derivative dy/dt

        dydt = -y
        if (t < -huge(t)) dydt = 0.0_dp
    end subroutine decay_rhs

    subroutine decay_jacobian(t, y, jac)
        real(dp), intent(in) :: t !! independent variable, unused for autonomous system
        real(dp), intent(in) :: y(:) !! current decay state, used for shape validation
        real(dp), intent(out) :: jac(:, :) !! scalar decay Jacobian

        jac = -1.0_dp
        if (t < -huge(t) .or. size(y) < 1) jac = 0.0_dp
    end subroutine decay_jacobian

    subroutine check_rel(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tolerance !! permitted relative error
        character(len=*), intent(in) :: label !! assertion label printed on failure

        real(dp) :: scale

        scale = max(abs(expected), tiny(1.0_dp))
        call require(abs(actual - expected) <= tolerance * scale, label)
    end subroutine check_rel

    subroutine require(condition, label)
        logical, intent(in) :: condition !! assertion condition
        character(len=*), intent(in) :: label !! assertion label printed on failure

        if (.not. condition) then
            print *, 'FAIL ', trim(label)
            error stop 1
        end if
    end subroutine require

end program test_integrate_lsoda_transfer
