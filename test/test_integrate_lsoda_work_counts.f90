program test_integrate_lsoda_work_counts
    use scifort_integrate_ivp, only : solve_ivp, solve_ivp_result
    use scifort_kinds, only : dp
    implicit none

    call test_scalar_default_tolerance()
    call test_scalar_tight_tolerance()
    call test_vector_weighted_max_norm()

    print *, 'test_integrate_lsoda_work_counts: PASS'

contains

    subroutine test_scalar_default_tolerance()
        type(solve_ivp_result) :: sol

        sol = solve_ivp(exponential_rhs, [0.0_dp, 1.0_dp], [1.0_dp], method='LSODA', &
            rtol=1.0e-6_dp, atol=1.0e-9_dp)
        call require(sol%success, 'LSODA default-tolerance exponential success')
        call require(sol%nfev == 35, 'LSODA default-tolerance exponential nfev')
        call require(sol%njev == 0, 'LSODA default-tolerance exponential njev')
        call require(sol%nlu == 0, 'LSODA default-tolerance exponential nlu')
        call require(sol%naccepted == 17, &
            'LSODA default-tolerance exponential accepted-step count')
        call require(sol%nrejected == 0, &
            'LSODA default-tolerance exponential rejection count')
        call check_abs(sol%t(2), 7.0746006954465273e-4_dp, 2.0e-18_dp, &
            'LSODA ODEPACK initial step')
        call check_abs(sol%y(1, size(sol%t)), 2.7182852838645504_dp, 2.0e-14_dp, &
            'LSODA default-tolerance exponential endpoint')
    end subroutine test_scalar_default_tolerance

    subroutine test_scalar_tight_tolerance()
        type(solve_ivp_result) :: sol

        sol = solve_ivp(exponential_rhs, [0.0_dp, 1.0_dp], [1.0_dp], method='LSODA', &
            rtol=1.0e-8_dp, atol=1.0e-10_dp)
        call require(sol%success, 'LSODA tight-tolerance exponential success')
        call require(sol%nfev == 59, 'LSODA tight-tolerance exponential nfev')
        call require(sol%njev == 0, 'LSODA tight-tolerance exponential njev')
        call require(sol%nlu == 0, 'LSODA tight-tolerance exponential nlu')
        call require(sol%naccepted == 26, &
            'LSODA tight-tolerance exponential accepted-step count')
        call require(sol%nrejected == 3, &
            'LSODA tight-tolerance exponential rejection count')
        call check_abs(sol%t(2), 7.1061593089204135e-5_dp, 2.0e-19_dp, &
            'LSODA tight ODEPACK initial step')
        call check_abs(sol%y(1, size(sol%t)), 2.7182818625584528_dp, 2.0e-14_dp, &
            'LSODA tight-tolerance exponential endpoint')
    end subroutine test_scalar_tight_tolerance

    subroutine test_vector_weighted_max_norm()
        type(solve_ivp_result) :: sol
        real(dp) :: y0(8)

        y0 = 1.0_dp
        sol = solve_ivp(vector_rhs, [0.0_dp, 1.0_dp], y0, method='LSODA', &
            rtol=1.0e-6_dp, atol=1.0e-9_dp)
        call require(sol%success, 'LSODA vector max-norm success')
        call require(sol%nfev == 101, 'LSODA vector max-norm nfev')
        call require(sol%njev == 0, 'LSODA vector max-norm njev')
        call require(sol%nlu == 0, 'LSODA vector max-norm nlu')
        call require(sol%naccepted == 50, 'LSODA vector max-norm accepted-step count')
        call check_abs(sol%y(1, size(sol%t)), 2.718281867538598_dp, 2.0e-13_dp, &
            'LSODA vector max-norm endpoint')
    end subroutine test_vector_weighted_max_norm

    subroutine exponential_rhs(t, y, dydt)
        real(dp), intent(in) :: t !! independent variable, unused for autonomous system
        real(dp), intent(in) :: y(:) !! current solution state
        real(dp), intent(out) :: dydt(:) !! derivative dy/dt

        dydt = y
        if (t < -huge(t)) dydt = 0.0_dp
    end subroutine exponential_rhs

    subroutine vector_rhs(t, y, dydt)
        real(dp), intent(in) :: t !! independent variable, unused for autonomous system
        real(dp), intent(in) :: y(:) !! current solution state
        real(dp), intent(out) :: dydt(:) !! derivative dy/dt

        integer :: i

        do i = 1, size(y)
            dydt(i) = real(i, dp) * y(i)
        end do
        if (t < -huge(t)) dydt = 0.0_dp
    end subroutine vector_rhs

    subroutine check_abs(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tolerance !! permitted absolute error
        character(len=*), intent(in) :: label !! assertion label printed on failure

        if (abs(actual - expected) > tolerance) then
            print *, 'FAIL ', trim(label)
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance= ', tolerance
            error stop 1
        end if
    end subroutine check_abs

    subroutine require(condition, label)
        logical, intent(in) :: condition !! assertion condition
        character(len=*), intent(in) :: label !! assertion label printed on failure

        if (.not. condition) then
            print *, 'FAIL ', trim(label)
            error stop 1
        end if
    end subroutine require

end program test_integrate_lsoda_work_counts
