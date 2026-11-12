program test_integrate_odeint
    use scifort_integrate, only : dp, ivp_status_invalid_argument, ivp_status_max_steps, &
        odeint, odeint_result, odeint_tfirst
    implicit none

    logical :: critical_seen = .false.
    logical :: crossed_before_critical = .false.

    call test_scipy_pendulum_reference()
    call test_vector_rtol()
    call test_backward_repeated_times()
    call test_degenerate_time_grids()
    call test_banded_jacobian()
    call test_col_deriv_banded_jacobian()
    call test_tfirst_callback()
    call test_controls()
    call test_mxstep_per_output()
    call test_full_output_diagnostics()
    call test_full_output_stiff_history()
    call test_full_output_repeated_backward()
    call test_tcrit()
    call test_order_controls()
    call test_printing_controls()
    call test_invalid_arguments()

    print *, 'test_integrate_odeint: PASS'

contains

    subroutine test_scipy_pendulum_reference()
        type(odeint_result) :: sol
        real(dp), parameter :: times(7) = [0.0_dp, 1.0_dp, 2.0_dp, 4.0_dp, &
            6.0_dp, 8.0_dp, 10.0_dp]
        real(dp), parameter :: reference(7, 2) = reshape([ &
            3.041592653589793_dp, 2.702109010032864_dp, 0.246765393631014_dp, &
            0.421440796203147_dp, -1.051941719663416_dp, 0.836462554680834_dp, &
            0.020011693078482_dp, 0.0_dp, -0.901061118862392_dp, &
            -4.003435766589569_dp, 3.026210612083851_dp, -1.077461734208435_dp, &
            -0.967176556573099_dp, 1.567818217897237_dp], [7, 2])
        integer :: i
        integer :: j

        sol = odeint(pendulum_rhs, [acos(-1.0_dp) - 0.1_dp, 0.0_dp], times)

        call require(sol%success, 'odeint pendulum success')
        call require(size(sol%y, 1) == size(times), 'odeint pendulum output rows')
        call require(size(sol%y, 2) == 2, 'odeint pendulum output columns')
        do j = 1, 2
            do i = 1, size(times)
                call check_abs(sol%y(i, j), reference(i, j), 3.0e-6_dp, &
                    'odeint SciPy pendulum reference')
            end do
        end do
        call require(sol%nfev > 0, 'odeint pendulum RHS work reported')
    end subroutine test_scipy_pendulum_reference

    subroutine test_vector_rtol()
        type(odeint_result) :: scalar_tol
        type(odeint_result) :: vector_tol
        real(dp), parameter :: times(2) = [0.0_dp, 0.2_dp]
        real(dp) :: expected_fast

        expected_fast = exp(-2.0_dp)
        scalar_tol = odeint(multiscale_decay_rhs_yt, [1.0_dp, 1.0_dp], times, &
            rtol=1.0e-3_dp, atol=1.0e-12_dp)
        vector_tol = odeint(multiscale_decay_rhs_yt, [1.0_dp, 1.0_dp], times, &
            rtol=[1.0e-3_dp, 1.0e-8_dp], atol=1.0e-12_dp)

        call require(scalar_tol%success, 'odeint scalar rtol success')
        call require(vector_tol%success, 'odeint vector rtol success')
        call check_abs(vector_tol%y(2, 2), expected_fast, 2.0e-8_dp, &
            'odeint vector rtol fast component')
        call require(vector_tol%nfev > scalar_tol%nfev, &
            'odeint tighter vector rtol increases RHS work')
    end subroutine test_vector_rtol

    subroutine test_backward_repeated_times()
        type(odeint_result) :: sol
        real(dp), parameter :: times(5) = [1.0_dp, 0.75_dp, 0.75_dp, 0.25_dp, 0.0_dp]
        real(dp), parameter :: reference(5) = [ &
            0.367879441171442_dp, 0.472366553674602_dp, 0.472366553674602_dp, &
            0.778800802387897_dp, 1.000000024272113_dp]
        integer :: i

        sol = odeint(decay_rhs_yt, [exp(-1.0_dp)], times)

        call require(sol%success, 'odeint backward repeated-time success')
        call require(size(sol%y, 1) == size(times), 'odeint repeated-time row count')
        call check_abs(sol%y(2, 1), sol%y(3, 1), 0.0_dp, &
            'odeint repeated time duplicates state')
        do i = 1, size(times)
            call check_abs(sol%y(i, 1), reference(i), 2.0e-7_dp, &
                'odeint backward SciPy reference')
        end do
    end subroutine test_backward_repeated_times

    subroutine test_degenerate_time_grids()
        type(odeint_result) :: sol

        sol = odeint(decay_rhs_yt, [2.0_dp], [1.0_dp])
        call require(sol%success, 'odeint single-time grid success')
        call require(all(shape(sol%y) == [1, 1]), 'odeint single-time output shape')
        call check_abs(sol%y(1, 1), 2.0_dp, 0.0_dp, 'odeint single-time state')

        sol = odeint(decay_rhs_yt, [2.0_dp], [1.0_dp, 1.0_dp, 1.0_dp])
        call require(sol%success, 'odeint all-equal grid success')
        call require(all(shape(sol%y) == [3, 1]), 'odeint all-equal output shape')
        call require(all(sol%y(:, 1) == 2.0_dp), 'odeint all-equal states unchanged')
    end subroutine test_degenerate_time_grids

    subroutine test_banded_jacobian()
        type(odeint_result) :: sol
        real(dp), parameter :: times(4) = [0.0_dp, 0.1_dp, 0.5_dp, 1.0_dp]
        integer :: i

        sol = odeint(stiff_cosine_rhs_yt, [1.0_dp], times, rtol=1.0e-8_dp, &
            atol=1.0e-10_dp, dfun=stiff_cosine_jacobian_yt, ml=0, mu=0)

        call require(sol%success, 'odeint banded Jacobian success')
        call require(sol%njev > 0, 'odeint banded Jacobian used')
        call require(sol%nlu > 0, 'odeint banded Jacobian LU work')
        do i = 1, size(times)
            call check_abs(sol%y(i, 1), cos(times(i)), 1.0e-7_dp, &
                'odeint banded Jacobian value')
        end do
    end subroutine test_banded_jacobian

    subroutine test_col_deriv_banded_jacobian()
        type(odeint_result) :: sol
        real(dp), parameter :: times(3) = [0.0_dp, 0.25_dp, 1.0_dp]
        real(dp) :: y0(2)
        real(dp) :: expected(2)
        integer :: i

        y0 = [100.0_dp / 999.0_dp, 1.0_dp]
        sol = odeint(upper_banded_rhs_yt, y0, times, rtol=1.0e-8_dp, &
            atol=1.0e-10_dp, dfun=upper_banded_jacobian_transposed_yt, &
            col_deriv=.true., ml=0, mu=1)

        call require(sol%success, 'odeint col_deriv banded Jacobian success')
        call require(sol%njev > 0, 'odeint col_deriv banded Jacobian used')
        do i = 1, size(times)
            expected = [100.0_dp / 999.0_dp, 1.0_dp] * exp(-times(i))
            call check_abs(sol%y(i, 1), expected(1), 2.0e-7_dp, &
                'odeint col_deriv first component')
            call check_abs(sol%y(i, 2), expected(2), 2.0e-7_dp, &
                'odeint col_deriv second component')
        end do
    end subroutine test_col_deriv_banded_jacobian

    subroutine test_tfirst_callback()
        type(odeint_result) :: sol

        sol = odeint_tfirst(decay_rhs_ty, [1.0_dp], [0.0_dp, 1.0_dp])

        call require(sol%success, 'odeint_tfirst success')
        call check_abs(sol%y(2, 1), exp(-1.0_dp), 2.0e-7_dp, &
            'odeint_tfirst value')
    end subroutine test_tfirst_callback

    subroutine test_controls()
        type(odeint_result) :: unrestricted
        type(odeint_result) :: limited
        real(dp), parameter :: times(2) = [0.0_dp, 1.0_dp]

        unrestricted = odeint(decay_rhs_yt, [1.0_dp], times)
        limited = odeint(decay_rhs_yt, [1.0_dp], times, h0=1.0e-3_dp, &
            hmax=0.05_dp, hmin=0.0_dp, mxstep=10000)

        call require(unrestricted%success, 'odeint unrestricted controls success')
        call require(limited%success, 'odeint explicit controls success')
        call check_abs(limited%y(2, 1), exp(-1.0_dp), 2.0e-7_dp, &
            'odeint explicit controls value')
        call require(limited%naccepted >= unrestricted%naccepted, &
            'odeint hmax limits accepted step size')
    end subroutine test_controls

    subroutine test_mxstep_per_output()
        type(odeint_result) :: too_small
        type(odeint_result) :: enough
        type(odeint_result) :: default_cap
        type(odeint_result) :: explicit_default
        real(dp), parameter :: times(3) = [0.0_dp, 0.1_dp, 1.0_dp]

        too_small = odeint(decay_rhs_yt, [1.0_dp], times, rtol=1.0e-4_dp, &
            atol=1.0e-8_dp, hmax=0.15_dp, mxstep=6)
        call require(.not. too_small%success, 'odeint mxstep enforced per output call')
        call require(too_small%status == ivp_status_max_steps, &
            'odeint per-output mxstep status')
        call require(too_small%naccepted == 10, &
            'odeint mxstep=6 SciPy accepted-step count')
        call require(too_small%nfev == 21, 'odeint mxstep=6 SciPy RHS work count')

        enough = odeint(decay_rhs_yt, [1.0_dp], times, rtol=1.0e-4_dp, &
            atol=1.0e-8_dp, hmax=0.15_dp, mxstep=7)
        call require(enough%success, 'odeint sufficient per-output mxstep succeeds')
        call require(enough%naccepted == 11, 'odeint mxstep=7 SciPy accepted-step count')
        call require(enough%nfev == 23, 'odeint mxstep=7 SciPy RHS work count')
        call check_abs(enough%y(3, 1), exp(-1.0_dp), 2.0e-4_dp, &
            'odeint sufficient mxstep value')

        default_cap = odeint(decay_rhs_yt, [1.0_dp], [0.0_dp, 1.0_dp], mxstep=0)
        explicit_default = odeint(decay_rhs_yt, [1.0_dp], [0.0_dp, 1.0_dp], mxstep=500)
        call require(default_cap%success .and. explicit_default%success, &
            'odeint default mxstep runs succeed')
        call require(default_cap%naccepted == explicit_default%naccepted, &
            'odeint mxstep=0 selects 500-step default')
        call require(default_cap%nfev == explicit_default%nfev, &
            'odeint mxstep default work matches explicit 500')
    end subroutine test_mxstep_per_output

    subroutine test_full_output_diagnostics()
        type(odeint_result) :: basic
        type(odeint_result) :: disabled
        type(odeint_result) :: tfirst
        real(dp), parameter :: times(3) = [0.0_dp, 0.1_dp, 0.2_dp]

        basic = odeint(decay_rhs_yt, [1.0_dp], times, full_output=.true.)
        call require(basic%success, 'odeint full_output success')
        call require(size(basic%hu) == 2, 'odeint full_output hu size')
        call require(size(basic%tcur) == 2, 'odeint full_output tcur size')
        call require(size(basic%tolsf) == 2, 'odeint full_output tolsf size')
        call require(size(basic%tsw) == 2, 'odeint full_output tsw size')
        call require(size(basic%nst) == 2, 'odeint full_output nst size')
        call require(size(basic%nfe) == 2, 'odeint full_output nfe size')
        call require(size(basic%nje) == 2, 'odeint full_output nje size')
        call require(size(basic%nqu) == 2, 'odeint full_output nqu size')
        call require(size(basic%mused) == 2, 'odeint full_output mused size')
        call require(all(basic%nst == [10, 13]), 'odeint full_output SciPy nst counts')
        call require(all(basic%nfe == [23, 29]), 'odeint full_output SciPy nfe counts')
        call require(all(basic%nje == [0, 0]), 'odeint full_output SciPy nje counts')
        call require(all(basic%nqu == [4, 4]), 'odeint full_output SciPy orders')
        call require(all(basic%mused == [1, 1]), 'odeint full_output Adams indicators')
        call require(all(basic%hu > 0.0_dp), 'odeint full_output positive forward steps')
        call require(all(basic%tcur >= times(2:)), 'odeint full_output internal times reach outputs')
        call require(all(basic%tolsf == 0.0_dp), 'odeint full_output normal tolerance scales')
        call require(all(basic%tsw == 0.0_dp), 'odeint full_output no-switch time')
        call require(basic%imxer == -1, 'odeint full_output imxer success value')
        call require(basic%lenrw == 36, 'odeint full_output ODEPACK lenrw')
        call require(basic%leniw == 21, 'odeint full_output ODEPACK leniw')

        tfirst = odeint_tfirst(decay_rhs_ty, [1.0_dp], times, full_output=.true.)
        call require(tfirst%success, 'odeint_tfirst full_output success')
        call require(all(tfirst%nst == basic%nst) .and. all(tfirst%nfe == basic%nfe), &
            'odeint_tfirst full_output work matches legacy callback')
        call require(all(tfirst%nqu == basic%nqu) .and. all(tfirst%mused == basic%mused), &
            'odeint_tfirst full_output method history matches legacy callback')

        disabled = odeint(decay_rhs_yt, [1.0_dp], times)
        call require(size(disabled%hu) == 0, 'odeint diagnostics disabled by default')
        call require(disabled%lenrw == 0 .and. disabled%leniw == 0, &
            'odeint work lengths disabled by default')
    end subroutine test_full_output_diagnostics

    subroutine test_full_output_stiff_history()
        type(odeint_result) :: sol
        real(dp), parameter :: times(4) = [0.0_dp, 0.1_dp, 0.5_dp, 1.0_dp]

        sol = odeint(stiff_cosine_rhs_yt, [1.0_dp], times, rtol=1.0e-8_dp, &
            atol=1.0e-10_dp, full_output=.true.)

        call require(sol%success, 'odeint stiff full_output success')
        call require(all(sol%mused == 2), 'odeint stiff full_output BDF indicators')
        call require(all(sol%nqu >= 1 .and. sol%nqu <= 5), &
            'odeint stiff full_output BDF order range')
        call require(all(sol%tsw > 0.0_dp .and. sol%tsw < times(2)), &
            'odeint stiff full_output switch time')
        call require(all(sol%tsw == sol%tsw(1)), &
            'odeint stiff full_output last switch retained')
        call require(all(sol%nst(2:) > sol%nst(:2)), &
            'odeint stiff full_output cumulative steps increase')
        call require(all(sol%nfe(2:) > sol%nfe(:2)), &
            'odeint stiff full_output cumulative RHS work increases')
        call require(all(sol%nje(2:) >= sol%nje(:2)), &
            'odeint stiff full_output cumulative Jacobian work')
    end subroutine test_full_output_stiff_history

    subroutine test_full_output_repeated_backward()
        type(odeint_result) :: sol
        real(dp), parameter :: times(4) = [1.0_dp, 1.0_dp, 0.5_dp, 0.0_dp]

        sol = odeint(decay_rhs_yt, [exp(-1.0_dp)], times, full_output=.true.)

        call require(sol%success, 'odeint backward full_output success')
        call require(size(sol%hu) == 3, 'odeint backward full_output size')
        call require(sol%hu(1) == 0.0_dp, &
            'odeint repeated initial output has deterministic zero step')
        call require(sol%tcur(1) == 1.0_dp, &
            'odeint repeated initial output reports initial internal time')
        call require(sol%nst(1) == 0 .and. sol%nfe(1) == 1, &
            'odeint repeated initial output reports initial work')
        call require(all(sol%hu(2:) < 0.0_dp), &
            'odeint backward full_output signed steps')
        call require(all(sol%tcur(2:) <= times(3:)), &
            'odeint backward full_output internal times reach outputs')
        call require(all(sol%mused == 1), 'odeint backward full_output Adams indicators')
        call require(all(sol%tsw == 1.0_dp), &
            'odeint backward full_output initial switch time')
    end subroutine test_full_output_repeated_backward

    subroutine test_tcrit()
        type(odeint_result) :: sol
        real(dp), parameter :: times(3) = [0.0_dp, 1.0_dp, 2.0_dp]

        critical_seen = .false.
        crossed_before_critical = .false.
        sol = odeint(tracked_decay_rhs_yt, [1.0_dp], times, tcrit=[1.0_dp])

        call require(sol%success, 'odeint tcrit success')
        call require(critical_seen, 'odeint tcrit reaches critical time exactly')
        call require(.not. crossed_before_critical, &
            'odeint tcrit prevents RHS evaluation past active critical time')
        call check_abs(sol%y(2, 1), exp(-1.0_dp), 2.0e-7_dp, &
            'odeint tcrit intermediate value')
        call check_abs(sol%y(3, 1), exp(-2.0_dp), 2.0e-7_dp, &
            'odeint tcrit final value')
    end subroutine test_tcrit

    subroutine test_order_controls()
        type(odeint_result) :: default_adams
        type(odeint_result) :: first_order_adams
        type(odeint_result) :: default_bdf
        type(odeint_result) :: first_order_bdf
        real(dp), parameter :: times(2) = [0.0_dp, 1.0_dp]

        default_adams = odeint(decay_rhs_yt, [1.0_dp], times, rtol=1.0e-4_dp, &
            atol=1.0e-8_dp)
        first_order_adams = odeint(decay_rhs_yt, [1.0_dp], times, rtol=1.0e-4_dp, &
            atol=1.0e-8_dp, mxordn=1)
        call require(default_adams%success .and. first_order_adams%success, &
            'odeint mxordn runs succeed')
        call require(first_order_adams%naccepted > default_adams%naccepted, &
            'odeint mxordn=1 increases nonstiff step count')
        call check_abs(first_order_adams%y(2, 1), exp(-1.0_dp), 3.0e-3_dp, &
            'odeint mxordn=1 value')

        default_bdf = odeint(stiff_cosine_rhs_yt, [1.0_dp], times, rtol=1.0e-4_dp, &
            atol=1.0e-8_dp)
        first_order_bdf = odeint(stiff_cosine_rhs_yt, [1.0_dp], times, rtol=1.0e-4_dp, &
            atol=1.0e-8_dp, mxords=1)
        call require(default_bdf%success .and. first_order_bdf%success, &
            'odeint mxords runs succeed')
        call require(first_order_bdf%naccepted > default_bdf%naccepted, &
            'odeint mxords=1 increases stiff step count')
        call check_abs(first_order_bdf%y(2, 1), cos(1.0_dp), 5.0e-5_dp, &
            'odeint mxords=1 value')
    end subroutine test_order_controls


    subroutine test_printing_controls()
        type(odeint_result) :: quiet
        type(odeint_result) :: verbose
        real(dp), parameter :: times(4) = [0.0_dp, 0.1_dp, 0.5_dp, 1.0_dp]

        quiet = odeint(stiff_cosine_rhs_yt, [1.0_dp], times, rtol=1.0e-8_dp, &
            atol=1.0e-10_dp, ixpr=.false., mxhnil=3, printmessg=.false.)
        verbose = odeint(stiff_cosine_rhs_yt, [1.0_dp], times, rtol=1.0e-8_dp, &
            atol=1.0e-10_dp, ixpr=.true., mxhnil=3, printmessg=.true.)

        call require(quiet%success .and. verbose%success, &
            'odeint printing-control runs succeed')
        call require(quiet%naccepted == verbose%naccepted, &
            'odeint printing controls do not change accepted steps')
        call require(quiet%nfev == verbose%nfev .and. quiet%njev == verbose%njev, &
            'odeint printing controls do not change work counts')
        call check_abs(verbose%y(size(times), 1), quiet%y(size(times), 1), 0.0_dp, &
            'odeint printing controls do not change solution')
    end subroutine test_printing_controls

    subroutine test_invalid_arguments()
        type(odeint_result) :: sol

        sol = odeint(decay_rhs_yt, [1.0_dp], [0.0_dp, 1.0_dp, 0.5_dp])
        call require(.not. sol%success, 'odeint nonmonotone times rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'odeint nonmonotone times status')

        sol = odeint(decay_rhs_yt, [1.0_dp], [0.0_dp, 1.0_dp], h0=-1.0_dp)
        call require(.not. sol%success, 'odeint negative h0 rejected')
        call require(sol%status == ivp_status_invalid_argument, 'odeint negative h0 status')

        sol = odeint(decay_rhs_yt, [1.0_dp], [0.0_dp, 1.0_dp], ml=1)
        call require(.not. sol%success, 'odeint invalid band rejected')
        call require(sol%status == ivp_status_invalid_argument, 'odeint invalid band status')

        sol = odeint(decay_rhs_yt, [1.0_dp], [0.0_dp, 1.0_dp], mxstep=-1)
        call require(.not. sol%success, 'odeint negative mxstep rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'odeint negative mxstep status')

        sol = odeint(decay_rhs_yt, [1.0_dp], [0.0_dp, 1.0_dp], mxhnil=-1)
        call require(.not. sol%success, 'odeint negative mxhnil rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'odeint negative mxhnil status')

        sol = odeint(multiscale_decay_rhs_yt, [1.0_dp, 1.0_dp], [0.0_dp, 0.2_dp], &
            rtol=[1.0e-3_dp])
        call require(.not. sol%success, 'odeint short vector rtol rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'odeint short vector rtol status')

        sol = odeint(decay_rhs_yt, [1.0_dp], [0.0_dp, 1.0_dp], mxordn=-1)
        call require(.not. sol%success, 'odeint negative mxordn rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'odeint negative mxordn status')

        sol = odeint(decay_rhs_yt, [1.0_dp], [0.0_dp, 1.0_dp], mxords=-1)
        call require(.not. sol%success, 'odeint negative mxords rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'odeint negative mxords status')

        sol = odeint(decay_rhs_yt, [1.0_dp], [0.0_dp, 1.0_dp], &
            tcrit=[1.0_dp, 0.5_dp])
        call require(.not. sol%success, 'odeint unordered tcrit rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'odeint unordered tcrit status')
    end subroutine test_invalid_arguments

    subroutine pendulum_rhs(y, t, dydt)
        real(dp), intent(in) :: y(:) !! pendulum angle and angular velocity
        real(dp), intent(in) :: t !! independent variable, unused by autonomous system
        real(dp), intent(out) :: dydt(:) !! pendulum state derivative

        if (t /= t .or. size(y) /= 2) error stop 'unexpected pendulum input'
        dydt(1) = y(2)
        dydt(2) = -0.25_dp * y(2) - 5.0_dp * sin(y(1))
    end subroutine pendulum_rhs

    subroutine multiscale_decay_rhs_yt(y, t, dydt)
        real(dp), intent(in) :: y(:) !! two-component state
        real(dp), intent(in) :: t !! independent variable, unused by autonomous system
        real(dp), intent(out) :: dydt(:) !! derivative with separated decay rates

        if (t /= t .or. size(y) /= 2) error stop 'unexpected multiscale odeint input'
        dydt(1) = -y(1)
        dydt(2) = -10.0_dp * y(2)
    end subroutine multiscale_decay_rhs_yt

    subroutine tracked_decay_rhs_yt(y, t, dydt)
        real(dp), intent(in) :: y(:) !! current exponential-decay state
        real(dp), intent(in) :: t !! independent variable
        real(dp), intent(out) :: dydt(:) !! derivative -y

        if (.not. critical_seen .and. t > 1.0_dp) crossed_before_critical = .true.
        if (abs(t - 1.0_dp) <= 32.0_dp * epsilon(1.0_dp)) critical_seen = .true.
        dydt = -y
    end subroutine tracked_decay_rhs_yt

    subroutine decay_rhs_yt(y, t, dydt)
        real(dp), intent(in) :: y(:) !! current exponential-decay state
        real(dp), intent(in) :: t !! independent variable, unused by autonomous system
        real(dp), intent(out) :: dydt(:) !! derivative -y

        if (t /= t) error stop 'unexpected NaN time'
        dydt = -y
    end subroutine decay_rhs_yt

    subroutine decay_rhs_ty(t, y, dydt)
        real(dp), intent(in) :: t !! independent variable, unused by autonomous system
        real(dp), intent(in) :: y(:) !! current exponential-decay state
        real(dp), intent(out) :: dydt(:) !! derivative -y

        if (t /= t) error stop 'unexpected NaN time'
        dydt = -y
    end subroutine decay_rhs_ty

    subroutine stiff_cosine_rhs_yt(y, t, dydt)
        real(dp), intent(in) :: y(:) !! current stiff scalar state
        real(dp), intent(in) :: t !! independent variable
        real(dp), intent(out) :: dydt(:) !! derivative with exact solution cos(t)

        dydt(1) = -1000.0_dp * (y(1) - cos(t)) - sin(t)
    end subroutine stiff_cosine_rhs_yt

    subroutine stiff_cosine_jacobian_yt(y, t, jac)
        real(dp), intent(in) :: y(:) !! current state, unused by constant Jacobian
        real(dp), intent(in) :: t !! independent variable, unused by constant Jacobian
        real(dp), intent(out) :: jac(:, :) !! packed diagonal Jacobian

        if (t /= t .or. size(y) /= 1) error stop 'unexpected stiff Jacobian input'
        jac = 0.0_dp
        jac(1, 1) = -1000.0_dp
    end subroutine stiff_cosine_jacobian_yt

    subroutine upper_banded_rhs_yt(y, t, dydt)
        real(dp), intent(in) :: y(:) !! two-component upper-triangular stiff state
        real(dp), intent(in) :: t !! independent variable, unused by autonomous system
        real(dp), intent(out) :: dydt(:) !! upper-triangular system derivative

        if (t /= t .or. size(y) /= 2) error stop 'unexpected upper-banded input'
        dydt(1) = -1000.0_dp * y(1) + 100.0_dp * y(2)
        dydt(2) = -y(2)
    end subroutine upper_banded_rhs_yt

    subroutine upper_banded_jacobian_transposed_yt(y, t, jac)
        real(dp), intent(in) :: y(:) !! current state, unused by constant Jacobian
        real(dp), intent(in) :: t !! independent variable, unused by constant Jacobian
        real(dp), intent(out) :: jac(:, :) !! transpose of packed upper-banded Jacobian

        if (t /= t .or. size(y) /= 2 .or. size(jac, 1) /= 2 .or. size(jac, 2) /= 2) then
            error stop 'unexpected transposed banded Jacobian shape'
        end if
        jac = 0.0_dp
        jac(1, 2) = -1000.0_dp
        jac(2, 1) = 100.0_dp
        jac(2, 2) = -1.0_dp
    end subroutine upper_banded_jacobian_transposed_yt

    subroutine require(condition, label)
        logical, intent(in) :: condition !! condition that must be true
        character(len=*), intent(in) :: label !! failure label

        if (.not. condition) then
            print *, 'FAIL: ', trim(label)
            error stop 1
        end if
    end subroutine require

    subroutine check_abs(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tolerance !! maximum absolute error
        character(len=*), intent(in) :: label !! failure label

        if (abs(actual - expected) > tolerance) then
            print *, 'FAIL: ', trim(label)
            print *, 'actual  = ', actual
            print *, 'expected= ', expected
            error stop 1
        end if
    end subroutine check_abs

end program test_integrate_odeint
