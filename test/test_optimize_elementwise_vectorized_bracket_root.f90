! SPDX-License-Identifier: MIT
! Compare packed vector root bracketing with SciPy and the scalar kernel.
program test_optimize_elementwise_vectorized_bracket_root
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_root_scalar, only : bracket_root, bracket_root_vectorized, &
        bracket_root_result, root_status_invalid_input, root_status_value_error, &
        root_status_convergence_error, root_status_success
    use optimize_elementwise_vectorized_brackets_reference, only : vector_bracket_root_status, &
        vector_bracket_root_nit, vector_bracket_root_nfev
    implicit none

    real(dp) :: left(7), right(7), lows(7), highs(7), factors(7), nan
    type(bracket_root_result) :: r(7), baseline, stopped(2), nonfinite(2)
    type(bracket_root_result), allocatable :: empty(:)
    integer :: failures, calls, points, initial_batch, smallest_batch, i

    failures = 0
    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    left = [0.0_dp, 0.0_dp, 2.0_dp, -3.0_dp, 6.0_dp, 0.0_dp, 1.0_dp]
    right = left + 0.5_dp
    lows = -10.0_dp
    highs = [10.0_dp, 10.0_dp, 10.0_dp, 10.0_dp, 10.0_dp, 0.75_dp, 10.0_dp]
    factors = [2.0_dp, 1.5_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp]

    call reset_counts()
    r = bracket_root_vectorized(vector_function, left, xr0=right, xmin=lows, &
        xmax=highs, factor=factors, maxiter=30)
    call check(all(r%status == vector_bracket_root_status), 'SciPy 1.17 status codes')
    call check(all(r%nit == vector_bracket_root_nit), 'SciPy 1.17 iterations')
    call check(all(r%nfev == vector_bracket_root_nfev), 'SciPy 1.17 evaluations')
    call check(initial_batch == 14, 'initial endpoints evaluated together')
    call check(smallest_batch < initial_batch, 'evaluation batches shrink')
    call check(points == sum(r%nfev), 'function-value count across lanes')
    call check(calls < points, 'fewer callback invocations than values')
    call check(r(7)%status == root_status_success .and. r(7)%nit == 0, &
        'initial endpoint root stops immediately')
    call check(r(6)%status == root_status_convergence_error, 'one lane reaches iteration cap')
    do i = 1, 7
        baseline = bracket_root(scalar_function, left(i), xr0=right(i), &
            xmin=lows(i), xmax=highs(i), factor=factors(i), maxiter=30)
        call same_result(r(i), baseline, 'scalar bracketing parity')
    end do

    ! Singleton vector arguments broadcast to all lanes.
    call reset_counts()
    r = bracket_root_vectorized(vector_function, left, factor=[2.0_dp], &
        xmin=[-10.0_dp], xmax=[10.0_dp], maxiter=12)
    do i = 1, 7
        baseline = bracket_root(scalar_function, left(i), factor=2.0_dp, &
            xmin=-10.0_dp, xmax=10.0_dp, maxiter=12)
        call same_result(r(i), baseline, 'singleton bounds and factor')
    end do
    call check(points == sum(r%nfev), 'singleton broadcast point accounting')

    ! A bad parameter excludes only that lane from the packed callback.
    factors(3) = 1.0_dp
    call reset_counts()
    r = bracket_root_vectorized(vector_function, left, xr0=right, xmin=lows, &
        xmax=highs, factor=factors, maxiter=30)
    call check(r(3)%status == root_status_invalid_input .and. r(3)%nfev == 0, &
        'invalid growth factor excluded')
    call check(initial_batch == 12 .and. r(1)%success, 'other lanes still evaluated')
    factors(3) = nan
    call reset_counts()
    r = bracket_root_vectorized(vector_function, left, xr0=right, factor=factors, maxiter=8)
    call check(r(3)%status == root_status_invalid_input .and. initial_batch == 12, &
        'nonfinite factor rejected without trap')
    factors(3) = 2.0_dp
    left(3) = nan
    call reset_counts()
    r = bracket_root_vectorized(vector_function, left, xr0=right, factor=factors, maxiter=8)
    call check(r(3)%status == root_status_invalid_input .and. initial_batch == 12, &
        'nonfinite start never reaches callback')
    left(3) = 2.0_dp

    ! Invalid broadcast dimensions reject the entire call without evaluations.
    call reset_counts()
    r = bracket_root_vectorized(vector_function, left, factor=[2.0_dp, 3.0_dp])
    call check(all(r%status == root_status_invalid_input) .and. calls == 0, &
        'invalid broadcast shape')
    call reset_counts()
    r = bracket_root_vectorized(vector_function, left, maxiter=-1)
    call check(all(r%status == root_status_invalid_input) .and. calls == 0, &
        'negative iteration count')
    call reset_counts()
    stopped = bracket_root_vectorized(vector_function, [1.0_dp, 0.0_dp], &
        xr0=[1.5_dp, 0.5_dp], maxiter=0)
    call check(stopped(1)%status == root_status_success .and. stopped(1)%nfev == 2, &
        'endpoint root with zero-iteration limit')
    call check(stopped(2)%status == root_status_convergence_error .and. &
        stopped(2)%nit == 0 .and. stopped(2)%nfev == 2, 'zero-iteration termination')
    call check(points == 4, 'zero-iteration initial endpoint accounting')

    ! Nonfinite output is isolated; another lane can still bracket its root.
    call reset_counts()
    nonfinite = bracket_root_vectorized(nonfinite_vector, [0.0_dp, 7.0_dp], &
        xr0=[0.5_dp, 7.5_dp], xmin=-10.0_dp, xmax=10.0_dp)
    call check(nonfinite(1)%success, 'unaffected lane converges')
    call check(nonfinite(2)%status == root_status_value_error .and. &
        nonfinite(2)%nfev == 2, 'initial nonfinite objective isolated')
    call check(points == sum(nonfinite%nfev), 'nonfinite work accounting')
    call reset_counts()
    nonfinite = bracket_root_vectorized(late_nonfinite_vector, [0.0_dp, 1.0_dp], &
        xr0=[0.5_dp, 1.5_dp], xmin=-10.0_dp, xmax=10.0_dp, maxiter=4)
    call check(nonfinite(1)%status == root_status_value_error .and. &
        nonfinite(1)%nit >= 1, 'late nonfinite expansion isolated')
    call check(nonfinite(2)%success, 'endpoint-root lane unaffected by late failure')

    call reset_counts()
    empty = bracket_root_vectorized(vector_function, [real(dp) :: ])
    call check(size(empty) == 0 .and. calls == 0, 'empty batch')

    if (failures /= 0) then
        print '(a,i0)', 'vectorized bracket_root regression FAIL: ', failures
        stop 1
    end if
    print '(a)', 'vectorized bracket_root regression PASS'
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
        type(bracket_root_result), intent(in) :: a, b
        character(len=*), intent(in) :: description
        call check(a%status == b%status .and. a%nfev == b%nfev .and. &
            a%nit == b%nit .and. a%success .eqv. b%success, description//' codes')
        if (a%nfev > 0 .and. b%nfev > 0) then
            call check(abs(a%xl-b%xl) < 1.0e-13_dp .and. abs(a%xr-b%xr) < 1.0e-13_dp, &
                description//' points')
            call check(abs(a%fl-b%fl) < 1.0e-13_dp .and. abs(a%fr-b%fr) < 1.0e-13_dp, &
                description//' values')
        end if
    end subroutine same_result

    real(dp) function scalar_function(x) result(fx)
        real(dp), intent(in) :: x
        fx = x - 1.0_dp
    end function scalar_function

    subroutine vector_function(x, fx)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: fx(:)
        call count_points(size(x))
        fx = x - 1.0_dp
    end subroutine vector_function

    subroutine nonfinite_vector(x, fx)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: fx(:)
        integer :: k
        call count_points(size(x))
        fx = x - 1.0_dp
        do k = 1, size(x)
            if (x(k) > 6.8_dp) fx(k) = nan
        end do
    end subroutine nonfinite_vector

    subroutine late_nonfinite_vector(x, fx)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: fx(:)
        integer :: k
        call count_points(size(x))
        fx = x - 1.0_dp
        do k = 1, size(x)
            if (x(k) > 2.0_dp) fx(k) = nan
        end do
    end subroutine late_nonfinite_vector

    subroutine count_points(n)
        integer, intent(in) :: n
        calls = calls + 1
        points = points + n
        if (calls == 1) initial_batch = n
        smallest_batch = min(smallest_batch, n)
    end subroutine count_points
end program test_optimize_elementwise_vectorized_bracket_root
