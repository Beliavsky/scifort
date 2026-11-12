! SPDX-License-Identifier: MIT
! Compare packed vector minimum bracketing with SciPy and the scalar kernel.
program test_optimize_elementwise_vectorized_bracket_minimum
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_minimize_scalar, only : bracket_minimum, bracket_minimum_vectorized, &
        bracket_minimum_result, minimize_status_invalid_input, minimize_status_value_error, &
        minimize_status_convergence_error, minimize_status_success
    use optimize_elementwise_vectorized_brackets_reference, only : vector_bracket_min_status, &
        vector_bracket_min_nit, vector_bracket_min_nfev
    implicit none

    real(dp) :: middle(7), lows(7), highs(7), factors(7), nan
    type(bracket_minimum_result) :: r(7), baseline, stopped(2), nonfinite(2)
    type(bracket_minimum_result), allocatable :: empty(:)
    integer :: failures, calls, points, initial_batch, smallest_batch, i

    failures = 0
    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    middle = [0.0_dp, 1.5_dp, 5.0_dp, -4.0_dp, 9.5_dp, 7.0_dp, 1.0_dp]
    lows = -10.0_dp
    highs = [10.0_dp, 10.0_dp, 10.0_dp, 10.0_dp, 10.0_dp, 10.0_dp, 1.1_dp]
    factors = 2.0_dp

    call reset_counts()
    r = bracket_minimum_vectorized(vector_function, middle, xmin=lows, xmax=highs, maxiter=25)
    call check(all(r%status == vector_bracket_min_status), 'SciPy 1.17 status codes')
    call check(all(r%nit == vector_bracket_min_nit), 'SciPy 1.17 iterations')
    call check(all(r%nfev == vector_bracket_min_nfev), 'SciPy 1.17 evaluations')
    call check(initial_batch == 21, 'initial triples evaluated together')
    call check(smallest_batch < initial_batch, 'evaluation batches shrink')
    call check(points == sum(r%nfev), 'function-value count across lanes')
    call check(calls < points, 'fewer callback invocations than values')
    call check(r(2)%status == minimize_status_success .and. r(2)%nit == 0, &
        'initial valid minimum bracket')
    call check(r(7)%status == minimize_status_convergence_error, &
        'one bounded lane reaches iteration cap')
    do i = 1, 7
        baseline = bracket_minimum(scalar_function, middle(i), &
            xmin=lows(i), xmax=highs(i), maxiter=25)
        call same_result(r(i), baseline, 'scalar bracketing parity')
    end do

    ! Scalar and length-one options are broadcast to every lane.
    call reset_counts()
    r = bracket_minimum_vectorized(vector_function, middle, factor=[2.0_dp], &
        xmin=[-10.0_dp], xmax=[10.0_dp], maxiter=12)
    do i = 1, 7
        baseline = bracket_minimum(scalar_function, middle(i), factor=2.0_dp, &
            xmin=-10.0_dp, xmax=10.0_dp, maxiter=12)
        call same_result(r(i), baseline, 'singleton bounds and factor')
    end do
    call check(points == sum(r%nfev), 'singleton broadcast accounting')

    factors(3) = 1.0_dp
    call reset_counts()
    r = bracket_minimum_vectorized(vector_function, middle, xmin=lows, &
        xmax=highs, factor=factors, maxiter=25)
    call check(r(3)%status == minimize_status_invalid_input .and. r(3)%nfev == 0, &
        'bad growth factor excludes one lane')
    call check(initial_batch == 18 .and. r(1)%success, 'other lanes unaffected')
    factors(3) = nan
    call reset_counts()
    r = bracket_minimum_vectorized(vector_function, middle, factor=factors, maxiter=8)
    call check(r(3)%status == minimize_status_invalid_input .and. initial_batch == 18, &
        'nonfinite factor rejected safely')
    factors(3) = 2.0_dp
    middle(3) = nan
    call reset_counts()
    r = bracket_minimum_vectorized(vector_function, middle, maxiter=8)
    call check(r(3)%status == minimize_status_invalid_input .and. initial_batch == 18, &
        'nonfinite middle point excludes lane')
    middle(3) = 5.0_dp

    call reset_counts()
    r = bracket_minimum_vectorized(vector_function, middle, factor=[2.0_dp, 3.0_dp])
    call check(all(r%status == minimize_status_invalid_input) .and. calls == 0, &
        'malformed broadcast rejects whole batch')
    call reset_counts()
    r = bracket_minimum_vectorized(vector_function, middle, maxiter=-1)
    call check(all(r%status == minimize_status_invalid_input) .and. calls == 0, &
        'negative iteration limit without function calls')
    call reset_counts()
    stopped = bracket_minimum_vectorized(vector_function, [1.5_dp, 0.0_dp], maxiter=0)
    call check(stopped(1)%status == minimize_status_success .and. stopped(1)%nfev == 3, &
        'initial minimum bracket with zero iterations')
    call check(stopped(2)%status == minimize_status_convergence_error .and. &
        stopped(2)%nit == 0 .and. stopped(2)%nfev == 3, &
        'zero-iteration termination')
    call check(points == 6, 'zero-iteration triple counts')

    call reset_counts()
    nonfinite = bracket_minimum_vectorized(nonfinite_vector, [0.0_dp, 7.0_dp], &
        xmin=-10.0_dp, xmax=10.0_dp, maxiter=10)
    call check(nonfinite(1)%success, 'finite lane unaffected by initial NaN')
    call check(nonfinite(2)%status == minimize_status_value_error .and. &
        nonfinite(2)%nfev == 3, 'nonfinite initial value isolated')
    call check(points == sum(nonfinite%nfev), 'nonfinite evaluation accounting')
    call reset_counts()
    nonfinite = bracket_minimum_vectorized(late_nonfinite_vector, [0.0_dp, 1.5_dp], &
        xmin=-10.0_dp, xmax=10.0_dp, maxiter=8)
    call check(nonfinite(1)%status == minimize_status_value_error .and. &
        nonfinite(1)%nit >= 1, 'late nonfinite expansion isolated')
    call check(nonfinite(2)%success, 'other minimum bracket not affected')

    call reset_counts()
    empty = bracket_minimum_vectorized(vector_function, [real(dp) :: ])
    call check(size(empty) == 0 .and. calls == 0, 'empty batch')

    if (failures /= 0) then
        print '(a,i0)', 'vectorized bracket_minimum regression FAIL: ', failures
        stop 1
    end if
    print '(a)', 'vectorized bracket_minimum regression PASS'
contains
    subroutine reset_counts()
        calls = 0
        points = 0
        initial_batch = 0
        smallest_batch = huge(0)
    end subroutine reset_counts

    subroutine check(ok, description)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: description
        if (.not. ok) then
            print '(a)', 'FAIL: '//description
            failures = failures + 1
        end if
    end subroutine check

    subroutine same_result(a, b, description)
        type(bracket_minimum_result), intent(in) :: a, b
        character(len=*), intent(in) :: description
        call check(a%status == b%status .and. a%nfev == b%nfev .and. &
            a%nit == b%nit .and. a%success .eqv. b%success, description//' codes')
        if (a%nfev > 0 .and. b%nfev > 0) then
            call check(abs(a%xl-b%xl) < 1.0e-13_dp .and. &
                abs(a%xm-b%xm) < 1.0e-13_dp .and. &
                abs(a%xr-b%xr) < 1.0e-13_dp, description//' points')
            call check(abs(a%fl-b%fl) < 1.0e-11_dp .and. &
                abs(a%fm-b%fm) < 1.0e-11_dp .and. &
                abs(a%fr-b%fr) < 1.0e-11_dp, description//' values')
        end if
    end subroutine same_result

    real(dp) function scalar_function(x) result(fx)
        real(dp), intent(in) :: x
        fx = (x-1.5_dp)**2 + 0.02_dp*(x-1.5_dp)**4
    end function scalar_function

    subroutine vector_function(x, fx)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: fx(:)
        call count_points(size(x))
        fx = (x-1.5_dp)**2 + 0.02_dp*(x-1.5_dp)**4
    end subroutine vector_function

    subroutine nonfinite_vector(x, fx)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: fx(:)
        integer :: k
        call count_points(size(x))
        fx = (x-1.5_dp)**2 + 0.02_dp*(x-1.5_dp)**4
        do k = 1, size(x)
            if (x(k) > 7.0_dp) fx(k) = nan
        end do
    end subroutine nonfinite_vector

    subroutine late_nonfinite_vector(x, fx)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: fx(:)
        integer :: k
        call count_points(size(x))
        fx = (x-1.5_dp)**2 + 0.02_dp*(x-1.5_dp)**4
        do k = 1, size(x)
            if (x(k) > 3.0_dp) fx(k) = nan
        end do
    end subroutine late_nonfinite_vector

    subroutine count_points(n)
        integer, intent(in) :: n
        calls = calls + 1
        points = points + n
        if (calls == 1) initial_batch = n
        smallest_batch = min(smallest_batch, n)
    end subroutine count_points
end program test_optimize_elementwise_vectorized_bracket_minimum
