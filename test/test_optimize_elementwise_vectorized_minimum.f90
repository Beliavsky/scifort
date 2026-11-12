! SPDX-License-Identifier: MIT
! Batched Chandrupatla minimizer: independent lanes and SciPy references.
program test_optimize_elementwise_vectorized_minimum
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_minimize_scalar, only : find_minimum, find_minimum_vectorized, &
        find_minimum_result, find_minimum_status_callback, minimize_status_invalid_input, &
        minimize_status_value_error, minimize_status_bracket_error, &
        minimize_status_convergence_error
    use optimize_elementwise_vectorized_minimum_reference, only : vector_min_ref_status, &
        vector_min_ref_nit, vector_min_ref_nfev
    implicit none

    integer :: failures, vector_calls, vector_points, initial_batch
    integer :: largest_active, smallest_active, callback_calls, i
    real(dp) :: brackets(3,5), nan, tolerance(5)
    type(find_minimum_result) :: r(5), ref, stopped(5), bad(2), zero_iter(2)
    type(find_minimum_result), allocatable :: empty(:)

    failures = 0
    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    brackets(:,1) = [ 0.0_dp, 1.0_dp, 3.0_dp]
    brackets(:,2) = [-1.0_dp, 2.0_dp, 6.0_dp]
    brackets(:,3) = [-4.0_dp, 1.0_dp, 4.0_dp]
    brackets(:,4) = [ 6.0_dp, 4.0_dp, 3.0_dp] ! not a minimum bracket
    brackets(:,5) = [ 3.0_dp, 1.0_dp, 0.0_dp] ! reversed

    call reset_counts()
    r = find_minimum_vectorized(vector_function, brackets, xatol=1.0e-12_dp)
    call check(all(r%status == vector_min_ref_status), 'SciPy 1.17 statuses')
    call check(all(r%nit == vector_min_ref_nit), 'SciPy 1.17 iteration counts')
    call check(all(r%nfev == vector_min_ref_nfev), 'SciPy 1.17 evaluation counts')
    call check(initial_batch == 15, 'all fifteen initial endpoint points batched')
    call check(vector_calls > 1 .and. vector_calls < vector_points, 'fewer callback calls than points')
    call check(smallest_active < largest_active, 'active-point batch shrinks')
    call check(vector_points == sum(r%nfev), 'vector evaluation point accounting')
    call check(r(4)%status == minimize_status_bracket_error .and. .not. r(4)%success, &
        'invalid bracket isolated')
    call check(all(r([1,2,3,5])%success), 'other problems converge')
    do i = 1, 5
        ref = find_minimum(scalar_function, brackets(:,i), xatol=1.0e-12_dp)
        call compare(r(i), ref, 'scalar parity for every lane')
        if (r(i)%success) then
            call check(abs(r(i)%x - 1.5_dp) < 1.0e-7_dp, 'SciPy reference minimizer')
        end if
    end do

    ! Broadcast independent tolerances; isolate a negative or NaN tolerance.
    tolerance = [1.0e-10_dp, 1.0e-8_dp, 1.0e-5_dp, 1.0e-12_dp, 1.0e-12_dp]
    call reset_counts()
    r = find_minimum_vectorized(vector_function, brackets, xatol=tolerance)
    do i = 1, 5
        ref = find_minimum(scalar_function, brackets(:,i), xatol=tolerance(i))
        call compare(r(i), ref, 'per-lane tolerance parity')
    end do
    tolerance(3) = -1.0_dp
    call reset_counts()
    r = find_minimum_vectorized(vector_function, brackets, xatol=tolerance)
    call check(r(3)%status == minimize_status_invalid_input .and. r(3)%nfev == 0, &
        'invalid tolerance excluded without evaluation')
    call check(r(2)%success .and. initial_batch == 12, 'other lanes still run')
    tolerance(3) = nan
    call reset_counts()
    r = find_minimum_vectorized(vector_function, brackets, xatol=tolerance)
    call check(r(3)%status == minimize_status_invalid_input .and. initial_batch == 12, &
        'NaN tolerance excluded without trap')

    call reset_counts()
    r = find_minimum_vectorized(vector_function, brackets, xatol=[1.0_dp, 2.0_dp])
    call check(all(r%status == minimize_status_invalid_input) .and. vector_calls == 0, &
        'malformed broadcast rejects entire batch')
    call reset_counts()
    zero_iter = find_minimum_vectorized(vector_function, brackets(:,1:2), maxiter=0)
    call check(all(zero_iter%status == minimize_status_convergence_error), 'zero iterations status')
    call check(all(zero_iter%nfev == 3) .and. vector_points == 6, 'zero iterations only evaluate bracket')
    call reset_counts()
    r = find_minimum_vectorized(vector_function, brackets, maxiter=-1)
    call check(all(r%status == minimize_status_invalid_input) .and. vector_calls == 0, &
        'negative maximum iterations no evaluations')

    ! Stop a single lane from the initial callback without stopping the rest.
    call reset_counts()
    callback_calls = 0
    stopped = find_minimum_vectorized(vector_function, brackets, &
        xatol=1.0e-12_dp, callback=stop_first)
    call check(stopped(1)%status == find_minimum_status_callback .and. &
        stopped(1)%nfev == 3, 'per-lane initial callback termination')
    call check(stopped(2)%success .and. callback_calls > 0, 'other lanes proceed')
    call check(vector_points == sum(stopped%nfev), 'callback accounting')

    ! An invalid starting coordinate never reaches the objective callback.
    brackets(1,3) = nan
    call reset_counts()
    r = find_minimum_vectorized(vector_function, brackets)
    call check(r(3)%status == minimize_status_invalid_input .and. initial_batch == 12, &
        'NaN bracket isolated')

    ! Bad function values at a starting endpoint and at a later probe.
    block
        real(dp) :: cases(3,2)
        cases(:,1) = [0.0_dp, 1.0_dp, 3.0_dp]
        cases(:,2) = [-1.0_dp, 2.0_dp, 7.0_dp]
        call reset_counts()
        bad = find_minimum_vectorized(nonfinite_vector, cases)
        call check(bad(1)%success, 'first lane unaffected by initial NaN')
        call check(bad(2)%status == minimize_status_value_error .and. bad(2)%nfev == 3, &
            'invalid initial objective isolated')
        call check(ieee_is_nan(bad(2)%x), 'invalid initial objective reports NaN result')
        call check(vector_points == sum(bad%nfev), 'nonfinite accounting')
        call reset_counts()
        bad = find_minimum_vectorized(late_nonfinite_vector, cases)
        call check(bad(1)%status == minimize_status_value_error .and. bad(1)%nit == 1, &
            'late nonfinite trial isolated')
        call check(ieee_is_nan(bad(1)%x), 'late nonfinite trial reports NaN minimizer')
    end block

    call reset_counts()
    empty = find_minimum_vectorized(vector_function, reshape([real(dp)::], [3,0]))
    call check(size(empty) == 0 .and. vector_calls == 0, 'empty grid')
    call reset_counts()
    r = find_minimum_vectorized(vector_function, reshape([real(dp)::], [0,5]))
    call check(all(r%status == minimize_status_invalid_input) .and. vector_calls == 0, &
        'invalid leading endpoint dimension')

    if (failures /= 0) then
        print '(a,i0)', 'vectorized minimum regression FAIL: ', failures
        stop 1
    end if
    print '(a)', 'vectorized elementwise minimum regression PASS'
contains
    subroutine reset_counts()
        vector_calls = 0
        vector_points = 0
        initial_batch = 0
        largest_active = 0
        smallest_active = huge(0)
    end subroutine reset_counts

    subroutine check(ok, label)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: label
        if (.not. ok) then
            print '(a)', 'FAIL: '//label
            failures = failures + 1
        end if
    end subroutine check

    subroutine compare(a, b, label)
        type(find_minimum_result), intent(in) :: a
        type(find_minimum_result), intent(in) :: b
        character(len=*), intent(in) :: label
        call check(a%status == b%status .and. a%nit == b%nit .and. a%nfev == b%nfev, label)
        if (a%success .and. b%success) then
            call check(abs(a%x - b%x) < 1.0e-13_dp .and. &
                abs(a%f_x - b%f_x) < 1.0e-13_dp, label//' minimizer')
            call check(abs(a%xl - b%xl) < 1.0e-13_dp .and. &
                abs(a%xr - b%xr) < 1.0e-13_dp, label//' interval')
        end if
    end subroutine compare

    real(dp) function scalar_function(x) result(fx)
        real(dp), intent(in) :: x
        fx = (x - 1.5_dp)**2 + 0.02_dp*(x - 1.5_dp)**4
    end function scalar_function

    subroutine vector_function(x, fx)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: fx(:)
        vector_calls = vector_calls + 1
        vector_points = vector_points + size(x)
        if (vector_calls == 1) initial_batch = size(x)
        if (vector_calls > 1) then
            largest_active = max(largest_active, size(x))
            smallest_active = min(smallest_active, size(x))
        end if
        fx = (x - 1.5_dp)**2 + 0.02_dp*(x - 1.5_dp)**4
    end subroutine vector_function

    subroutine nonfinite_vector(x, fx)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: fx(:)
        integer :: j
        vector_calls = vector_calls + 1
        vector_points = vector_points + size(x)
        do j = 1, size(x)
            if (x(j) > 6.0_dp) then
                fx(j) = nan
            else
                fx(j) = (x(j) - 1.5_dp)**2 + 0.02_dp*(x(j) - 1.5_dp)**4
            end if
        end do
    end subroutine nonfinite_vector

    subroutine late_nonfinite_vector(x, fx)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: fx(:)
        integer :: j
        vector_calls = vector_calls + 1
        vector_points = vector_points + size(x)
        do j = 1, size(x)
            if (x(j) > 1.5_dp .and. x(j) < 2.0_dp) then
                fx(j) = nan
            else
                fx(j) = (x(j) - 1.5_dp)**2 + 0.02_dp*(x(j) - 1.5_dp)**4
            end if
        end do
    end subroutine late_nonfinite_vector

    subroutine stop_first(state, stop)
        type(find_minimum_result), intent(in) :: state
        logical, intent(out) :: stop
        callback_calls = callback_calls + 1
        stop = state%nit == 0 .and. state%xl == 0.0_dp .and. state%xr == 3.0_dp
    end subroutine stop_first
end program test_optimize_elementwise_vectorized_minimum
