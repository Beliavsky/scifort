! SPDX-License-Identifier: MIT
! Compare active-set vectorized Chandrupatla refinement with scalar solves.
program test_optimize_elementwise_vectorized_root
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_root_scalar, only : find_root, find_root_vectorized, &
        find_root_result, root_status_invalid_input, root_status_value_error, &
        root_status_convergence_error, root_status_sign_error, find_root_status_callback
    use optimize_elementwise_vectorized_root_reference, only : vector_root_ref_status, &
        vector_root_ref_nit, vector_root_ref_nfev
    implicit none

    integer :: failures, vector_calls, scalar_calls, vector_points, callback_calls
    integer :: initial_batch, max_active, minimum_batch, i
    real(dp) :: brackets(2,5), nan, tolerance(5)
    type(find_root_result) :: r(5), ref, stopped(5), nonfinite(2), zero_iter(2)
    type(find_root_result), allocatable :: empty(:)

    failures = 0
    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    brackets(:,1) = [0.0_dp, 1.0_dp]  ! endpoint root
    brackets(:,2) = [0.0_dp, 2.0_dp]  ! normal refinement
    brackets(:,3) = [-4.0_dp, 4.0_dp] ! wider interval
    brackets(:,4) = [3.0_dp, 4.0_dp]  ! no sign change
    brackets(:,5) = [2.0_dp, -1.0_dp]  ! reversed interval

    call reset_counts()
    r = find_root_vectorized(vector_function, brackets, xatol=1.0e-12_dp)
    call check(all(r%status == vector_root_ref_status), 'SciPy 1.17 batch statuses')
    call check(all(r%nit == vector_root_ref_nit), 'SciPy 1.17 batch iterations')
    call check(all(r%nfev == vector_root_ref_nfev), 'SciPy 1.17 batch evaluation counts')
    call check(vector_calls > 1, 'more than one vector callback')
    call check(initial_batch == 10, 'initial call packs ten endpoint values')
    call check(max_active >= 3, 'multiple independent active lanes')
    call check(minimum_batch < max_active, 'active callback batch shrinks')
    call check(vector_points == sum(r%nfev), 'vector point accounting')
    call check(.not. r(4)%success .and. r(4)%status == root_status_sign_error, &
        'invalid bracket does not stop other lanes')
    call check(all(r([1,2,3,5])%success), 'successful root solves')
    do i = 1, 5
        ref = find_root(scalar_function, brackets(:,i), xatol=1.0e-12_dp)
        call check(r(i)%status == ref%status, 'scalar status for each lane')
        call check(r(i)%nit == ref%nit, 'scalar iteration counts')
        call check(r(i)%nfev == ref%nfev, 'scalar point counts')
        if (r(i)%success) then
            call check(abs(r(i)%x-ref%x) < 1.0e-14_dp, 'scalar root positions')
            call check(abs(r(i)%f_x-ref%f_x) < 1.0e-14_dp, 'scalar root residuals')
        end if
    end do
    call check(vector_calls < vector_points, 'vectorization reduces callback invocations')

    ! Per-lane valid tolerance vectors and invalid tolerance isolation.
    tolerance = [1.0e-12_dp, 1.0e-10_dp, 1.0e-8_dp, 1.0e-12_dp, 1.0e-12_dp]
    call reset_counts()
    r = find_root_vectorized(vector_function, brackets, xatol=tolerance)
    do i = 1, 5
        ref = find_root(scalar_function, brackets(:,i), xatol=tolerance(i))
        call check(r(i)%status == ref%status, 'array tolerance status')
        call check(r(i)%nfev == ref%nfev, 'array tolerance work count')
    end do
    tolerance(3) = -1.0_dp
    call reset_counts()
    r = find_root_vectorized(vector_function, brackets, xatol=tolerance)
    call check(r(3)%status == root_status_invalid_input .and. r(3)%nfev == 0, &
        'negative per-lane tolerance is excluded')
    call check(r(2)%success, 'valid lanes still converge with one invalid tolerance')
    call check(initial_batch == 8, 'invalid tolerance not sent to callback')
    tolerance(3) = nan
    call reset_counts()
    r = find_root_vectorized(vector_function, brackets, xatol=tolerance)
    call check(r(3)%status == root_status_invalid_input, 'NaN tolerance rejected without trap')
    call check(initial_batch == 8, 'NaN tolerance not evaluated')

    ! A malformed broadcast shape rejects the batch before any evaluation.
    call reset_counts()
    r = find_root_vectorized(vector_function, brackets, xatol=[1.0_dp, 2.0_dp])
    call check(all(r%status == root_status_invalid_input), 'invalid broadcast')
    call check(vector_calls == 0, 'no callback for invalid broadcast')
    call reset_counts()
    zero_iter = find_root_vectorized(vector_function, brackets(:,1:2), maxiter=0)
    call check(zero_iter(1)%success, 'zero-iteration endpoint root')
    call check(zero_iter(2)%status == root_status_convergence_error, 'zero iteration limit')
    call check(vector_points == 4, 'zero-iteration endpoint evaluations')

    ! Terminate only one lane on the initial callback, keep other lanes active.
    call reset_counts()
    callback_calls = 0
    stopped = find_root_vectorized(vector_function, brackets, xatol=1.0e-12_dp, callback=stop_lane)
    call check(stopped(5)%status == find_root_status_callback, 'one callback-stopped lane')
    call check(stopped(2)%success, 'other lanes unaffected by stopped lane')
    call check(callback_calls > 0, 'scalar per-lane callback invoked')
    call check(vector_points == sum(stopped%nfev), 'callback termination work accounting')

    ! Input NaN and nonfinite function value are handled per lane.
    call reset_counts()
    block
        real(dp) :: nf_brackets(2,2)
        nf_brackets(:,1) = [0.0_dp, 2.0_dp]
        nf_brackets(:,2) = [0.0_dp, 7.0_dp]
        nonfinite = find_root_vectorized(nonfinite_vector, nf_brackets)
    end block
    call check(nonfinite(1)%success, 'unaffected lane with nonfinite neighbor')
    call check(nonfinite(2)%status == root_status_value_error, 'nonfinite callback output')
    call check(vector_points == sum(nonfinite%nfev), 'nonfinite value work accounting')
    brackets(1,3) = nan
    call reset_counts()
    r = find_root_vectorized(vector_function, brackets)
    call check(r(3)%status == root_status_invalid_input, 'NaN bracket rejected')
    call check(initial_batch == 8, 'NaN bracket not sent to vector callback')

    ! Zero-size batches and malformed leading endpoint dimension are safe.
    call reset_counts()
    empty = find_root_vectorized(vector_function, reshape([real(dp)::], [2,0]))
    call check(size(empty) == 0 .and. vector_calls == 0, 'empty batch')
    call reset_counts()
    r = find_root_vectorized(vector_function, reshape([real(dp)::], [0,5]))
    call check(all(r%status == root_status_invalid_input) .and. vector_calls == 0, &
        'invalid endpoint shape')

    if (failures /= 0) then
        print '(a,i0)', 'vectorized root regression FAIL: ', failures
        stop 1
    end if
    print '(a)', 'vectorized elementwise root regression PASS'
contains
    subroutine reset_counts()
        vector_calls = 0
        scalar_calls = 0
        vector_points = 0
        initial_batch = 0
        max_active = 0
        minimum_batch = huge(0)
    end subroutine reset_counts

    subroutine check(ok, label)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: label
        if (.not. ok) then
            print '(a)', 'FAIL: '//label
            failures = failures + 1
        end if
    end subroutine check

    real(dp) function scalar_function(x) result(fx)
        real(dp), intent(in) :: x
        scalar_calls = scalar_calls + 1
        fx = x**3 - 1.0_dp
    end function scalar_function

    subroutine vector_function(x, fx)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: fx(:)
        vector_calls = vector_calls + 1
        vector_points = vector_points + size(x)
        if (vector_calls == 1) initial_batch = size(x)
        if (vector_calls > 1) then
            max_active = max(max_active, size(x))
            minimum_batch = min(minimum_batch, size(x))
        end if
        fx = x**3 - 1.0_dp
    end subroutine vector_function

    subroutine nonfinite_vector(x, fx)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: fx(:)
        integer :: j
        vector_calls = vector_calls + 1
        vector_points = vector_points + size(x)
        do j = 1, size(x)
            if (x(j) > 5.0_dp) then
                fx(j) = nan
            else
                fx(j) = x(j)**3 - 1.0_dp
            end if
        end do
    end subroutine nonfinite_vector

    subroutine stop_lane(state, stop)
        type(find_root_result), intent(in) :: state
        logical, intent(out) :: stop
        callback_calls = callback_calls + 1
        stop = (state%nit == 0 .and. state%xl == -1.0_dp .and. state%xr == 2.0_dp)
    end subroutine stop_lane
end program test_optimize_elementwise_vectorized_root
