! SPDX-License-Identifier: MIT
! Rank-two packed-callback bracketing and refinement, with trailing broadcasting.
program test_optimize_elementwise_vectorized_rank2
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_root_scalar, only : find_root_vectorized, bracket_root_vectorized, &
        find_root_result, bracket_root_result, find_root_status_callback, root_status_invalid_input
    use scifort_optimize_minimize_scalar, only : find_minimum_vectorized, &
        bracket_minimum_vectorized, find_minimum_result, bracket_minimum_result, &
        find_minimum_status_callback, minimize_status_invalid_input
    implicit none

    integer, parameter :: m = 2, n = 3, lanes = m*n
    real(dp) :: root_brackets(2,m,n), min_brackets(3,m,n)
    real(dp) :: root_left(m,n), min_middle(m,n)
    real(dp) :: atols(n), scalar_atols(m,n), limits(m,1), factors(1,n), bad_factors(m,n)
    real(dp) :: nan
    type(find_root_result) :: roots(m,n), roots_stop(m,n), roots_bad(m,n)
    type(find_root_result) :: roots_one(lanes)
    type(find_minimum_result) :: mins(m,n), mins_stop(m,n), mins_bad(m,n)
    type(find_minimum_result) :: mins_one(lanes)
    type(bracket_root_result) :: root_bs(m,n), root_bs_one(lanes)
    type(bracket_minimum_result) :: min_bs(m,n), min_bs_one(lanes)
    type(find_root_result), allocatable :: empty_roots(:,:)
    type(find_minimum_result), allocatable :: empty_mins(:,:)
    type(bracket_root_result), allocatable :: empty_root_bs(:,:)
    type(bracket_minimum_result), allocatable :: empty_min_bs(:,:)
    integer :: failures, calls, points, first_batch, smallest_batch, stop_calls
    integer :: i, j

    failures = 0
    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    root_brackets(:,1,1) = [ 0.0_dp, 2.0_dp]
    root_brackets(:,2,1) = [-1.0_dp, 1.2_dp]
    root_brackets(:,1,2) = [ 1.0_dp, 2.0_dp]
    root_brackets(:,2,2) = [ 3.0_dp, 4.0_dp]
    root_brackets(:,1,3) = [-2.0_dp, 0.0_dp]
    root_brackets(:,2,3) = [ 0.5_dp, 2.0_dp]
    min_brackets(:,1,1) = [ 0.0_dp, 1.0_dp, 3.0_dp]
    min_brackets(:,2,1) = [-1.0_dp, 2.0_dp, 6.0_dp]
    min_brackets(:,1,2) = [-4.0_dp, 1.0_dp, 4.0_dp]
    min_brackets(:,2,2) = [ 6.0_dp, 4.0_dp, 3.0_dp]
    min_brackets(:,1,3) = [ 3.0_dp, 1.0_dp, 0.0_dp]
    min_brackets(:,2,3) = [ 0.0_dp, 1.5_dp, 3.0_dp]
    atols = [1.0e-12_dp, 1.0e-9_dp, 1.0e-11_dp]
    scalar_atols = spread(atols, dim=1, ncopies=m)
    limits(:,1) = [-6.0_dp, -8.0_dp]
    factors(1,:) = [2.0_dp, 3.0_dp, 2.0_dp]

    ! Rank-one vector callbacks are the exact numerical oracle: the rank-two
    ! wrappers must preserve every coordinate, status, and evaluation count.
    call reset()
    roots = find_root_vectorized(root_function, root_brackets, xatol=atols, &
        xrtol=1.0e-12_dp, fatol=0.0_dp)
    call check(first_batch == 2*lanes, 'root initial endpoint packing')
    call check(points == sum(roots%nfev), 'root per-lane accounting')
    call check(smallest_batch < first_batch .and. calls < points, 'root active batch shrinks')
    call check(roots(2,2)%status /= 0 .and. roots(1,1)%success, 'root lane isolation')
    call reset()
    roots_one = find_root_vectorized(root_function, reshape(root_brackets, [2,lanes]), &
        xatol=reshape(scalar_atols, [lanes]), xrtol=1.0e-12_dp, fatol=0.0_dp)
    call same_roots(roots, reshape(roots_one, [m,n]), 'root rank-one equivalence')

    call reset()
    mins = find_minimum_vectorized(minimum_function, min_brackets, xatol=atols, &
        fatol=reshape([0.0_dp, 1.0e-15_dp], [m,1]))
    call check(first_batch == 3*lanes, 'minimum initial triple packing')
    call check(points == sum(mins%nfev), 'minimum per-lane accounting')
    call check(mins(2,2)%status /= 0 .and. mins(1,1)%success, 'minimum lane isolation')
    call reset()
    mins_one = find_minimum_vectorized(minimum_function, reshape(min_brackets, [3,lanes]), &
        xatol=reshape(scalar_atols, [lanes]), &
        fatol=reshape(spread([0.0_dp,1.0e-15_dp],2,n), [lanes]))
    call same_mins(mins, reshape(mins_one, [m,n]), 'minimum rank-one equivalence')

    root_left = reshape([-3.0_dp, 0.0_dp, 5.0_dp, 1.0_dp, 6.0_dp, 0.5_dp], [m,n])
    min_middle = reshape([0.0_dp, 1.5_dp, 5.0_dp, -4.0_dp, 9.5_dp, 7.0_dp], [m,n])
    call reset()
    root_bs = bracket_root_vectorized(root_function, root_left, xmin=limits, &
        xmax=10.0_dp, factor=factors)
    call check(first_batch == 2*lanes, 'bracket root initial endpoints packed')
    call check(points == sum(root_bs%nfev), 'bracket root counts')
    call reset()
    root_bs_one = bracket_root_vectorized(root_function, reshape(root_left, [lanes]), &
        xmin=reshape(spread(limits,2,n), [lanes]), xmax=10.0_dp, &
        factor=reshape(spread(factors,1,m), [lanes]))
    call same_root_brackets(root_bs, reshape(root_bs_one,[m,n]), 'bracket root parity')

    call reset()
    min_bs = bracket_minimum_vectorized(minimum_function, min_middle, &
        xmin=-10.0_dp, xmax=10.0_dp, factor=factors)
    call check(first_batch == 3*lanes, 'bracket minimum initial triples packed')
    call check(points == sum(min_bs%nfev), 'bracket minimum counts')
    call reset()
    min_bs_one = bracket_minimum_vectorized(minimum_function, reshape(min_middle,[lanes]), &
        xmin=-10.0_dp, xmax=10.0_dp, factor=reshape(spread(factors,1,m),[lanes]))
    call same_min_brackets(min_bs, reshape(min_bs_one,[m,n]), 'bracket minimum parity')

    ! Present optional endpoints with a singleton axis must be forwarded as
    ! present. Absent endpoints must retain the rank-one defaults.
    call reset()
    root_bs = bracket_root_vectorized(root_function, root_left, &
        xr0=reshape([0.1_dp, 1.2_dp, 1.5_dp], [1,n]), factor=2.0_dp)
    call reset()
    root_bs_one = bracket_root_vectorized(root_function, reshape(root_left,[lanes]), &
        xr0=reshape(spread([0.1_dp,1.2_dp,1.5_dp],1,m),[lanes]), factor=2.0_dp)
    call same_root_brackets(root_bs, reshape(root_bs_one,[m,n]), 'root explicit endpoints')

    call reset()
    min_bs = bracket_minimum_vectorized(minimum_function, min_middle, &
        xl0=min_middle-0.25_dp, xr0=min_middle+0.3_dp, factor=2.0_dp)
    call reset()
    min_bs_one = bracket_minimum_vectorized(minimum_function, reshape(min_middle,[lanes]), &
        xl0=reshape(min_middle-0.25_dp,[lanes]), &
        xr0=reshape(min_middle+0.3_dp,[lanes]), factor=2.0_dp)
    call same_min_brackets(min_bs, reshape(min_bs_one,[m,n]), 'minimum explicit endpoints')

    ! A callback stops a single root or minimum, without aborting other lanes.
    call reset()
    stop_calls = 0
    roots_stop = find_root_vectorized(root_function, root_brackets, &
        xatol=1.0e-12_dp, callback=root_stop)
    call check(stop_calls > 0 .and. roots_stop(2,1)%status == find_root_status_callback, &
        'root scalar-state callback forwarded')
    call check(roots_stop(1,1)%success, 'other root survives callback')
    call check(points == sum(roots_stop%nfev), 'stopped root accounting')
    call reset()
    stop_calls = 0
    mins_stop = find_minimum_vectorized(minimum_function, min_brackets, &
        xatol=1.0e-12_dp, callback=minimum_stop)
    call check(stop_calls > 0 .and. mins_stop(1,2)%status == find_minimum_status_callback, &
        'minimum scalar-state callback forwarded')
    call check(mins_stop(1,1)%success, 'other minimum survives callback')
    call check(points == sum(mins_stop%nfev), 'stopped minimum accounting')

    ! Invalid broadcasting, input extents, or tolerances must not cause
    ! objective evaluations for rejected lanes.
    call reset()
    roots_bad = find_root_vectorized(root_function, root_brackets, xatol=[1.0_dp,2.0_dp])
    call check(all(roots_bad%status == root_status_invalid_input) .and. calls == 0, &
        'root incompatible trailing axis')
    call reset()
    mins_bad = find_minimum_vectorized(minimum_function, min_brackets, xatol=[1.0_dp,2.0_dp])
    call check(all(mins_bad%status == minimize_status_invalid_input) .and. calls == 0, &
        'minimum incompatible trailing axis')
    call reset()
    root_bs = bracket_root_vectorized(root_function, root_left, factor=[1.5_dp,2.0_dp])
    call check(all(root_bs%status == root_status_invalid_input) .and. calls == 0, &
        'root bracket incompatible factor')
    call reset()
    min_bs = bracket_minimum_vectorized(minimum_function, min_middle, xl0=[0.0_dp,1.0_dp])
    call check(all(min_bs%status == minimize_status_invalid_input) .and. calls == 0, &
        'minimum bracket incompatible endpoint')
    ! One invalid growth factor excludes only its corresponding lane, with
    ! other lanes still packed into the same initial vector evaluation.
    bad_factors = 2.0_dp
    bad_factors(1,3) = nan
    call reset()
    root_bs = bracket_root_vectorized(root_function, root_left, factor=bad_factors)
    call check(root_bs(1,3)%status == root_status_invalid_input .and. &
        root_bs(1,3)%nfev == 0 .and. first_batch == 2*(lanes-1), &
        'root bracket NaN factor isolated without trap')
    call reset()
    min_bs = bracket_minimum_vectorized(minimum_function, min_middle, factor=bad_factors)
    call check(min_bs(1,3)%status == minimize_status_invalid_input .and. &
        min_bs(1,3)%nfev == 0 .and. first_batch == 3*(lanes-1), &
        'minimum bracket NaN factor isolated without trap')
    call reset()
    scalar_atols = 1.0e-12_dp
    scalar_atols(1,2) = nan
    roots_bad = find_root_vectorized(root_function, root_brackets, xatol=scalar_atols)
    call check(roots_bad(1,2)%status == root_status_invalid_input .and. &
        roots_bad(1,2)%nfev == 0 .and. first_batch == 2*(lanes-1), &
        'root NaN tolerance excluded without trap')
    call reset()
    scalar_atols(1,2) = -1.0_dp
    mins_bad = find_minimum_vectorized(minimum_function, min_brackets, xatol=scalar_atols)
    call check(mins_bad(1,2)%status == minimize_status_invalid_input .and. &
        mins_bad(1,2)%nfev == 0 .and. first_batch == 3*(lanes-1), &
        'minimum negative tolerance excludes one lane')
    call reset()
    root_bs = bracket_root_vectorized(root_function, root_left, maxiter=-1)
    call check(all(root_bs%status == root_status_invalid_input) .and. calls == 0, &
        'root bracket negative maxiter')
    call reset()
    min_bs = bracket_minimum_vectorized(minimum_function, min_middle, maxiter=-1)
    call check(all(min_bs%status == minimize_status_invalid_input) .and. calls == 0, &
        'minimum bracket negative maxiter')

    call reset()
    roots_bad = find_root_vectorized(root_function, reshape([1.0_dp], [1,m,n], pad=[1.0_dp]))
    call check(all(roots_bad%status == root_status_invalid_input) .and. calls == 0, &
        'invalid root endpoint extent')
    call reset()
    mins_bad = find_minimum_vectorized(minimum_function, reshape([1.0_dp], [2,m,n], pad=[1.0_dp]))
    call check(all(mins_bad%status == minimize_status_invalid_input) .and. calls == 0, &
        'invalid minimum endpoint extent')
    call reset()
    empty_roots = find_root_vectorized(root_function, reshape([real(dp)::], [2,0,3]))
    empty_mins = find_minimum_vectorized(minimum_function, reshape([real(dp)::], [3,0,3]))
    empty_root_bs = bracket_root_vectorized(root_function, reshape([real(dp)::], [0,3]))
    empty_min_bs = bracket_minimum_vectorized(minimum_function, reshape([real(dp)::], [0,3]))
    call check(all(shape(empty_roots) == [0,3]) .and. all(shape(empty_mins) == [0,3]), &
        'empty refinement result shapes')
    call check(all(shape(empty_root_bs) == [0,3]) .and. &
        all(shape(empty_min_bs) == [0,3]) .and. calls == 0, &
        'empty bracketing result shapes without callback')

    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_elementwise_vectorized_rank2: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_elementwise_vectorized_rank2: PASS'
contains
    subroutine reset()
        calls = 0
        points = 0
        first_batch = 0
        smallest_batch = huge(0)
    end subroutine reset

    subroutine count_evaluations(npoints)
        integer, intent(in) :: npoints !! packed function evaluations in this call
        calls = calls + 1
        points = points + npoints
        if (calls == 1) first_batch = npoints
        smallest_batch = min(smallest_batch, npoints)
    end subroutine count_evaluations

    subroutine check(ok, description)
        logical, intent(in) :: ok !! regression predicate
        character(len=*), intent(in) :: description !! failure description
        if (.not. ok) then
            print '(a)', 'FAIL: '//description
            failures = failures + 1
        end if
    end subroutine check

    subroutine root_function(x, fx)
        real(dp), intent(in) :: x(:) !! packed abscissae
        real(dp), intent(out) :: fx(:) !! corresponding residual values
        call count_evaluations(size(x))
        fx = x**3 - 1.0_dp
    end subroutine root_function

    subroutine minimum_function(x, fx)
        real(dp), intent(in) :: x(:) !! packed abscissae
        real(dp), intent(out) :: fx(:) !! corresponding objective values
        call count_evaluations(size(x))
        fx = (x-1.5_dp)**2 + 0.02_dp*(x-1.5_dp)**4
    end subroutine minimum_function

    subroutine root_stop(state, stop)
        type(find_root_result), intent(in) :: state !! active root state
        logical, intent(out) :: stop !! whether to stop this lane
        stop_calls = stop_calls + 1
        stop = state%xl < -0.5_dp .and. state%xr > 1.0_dp
    end subroutine root_stop

    subroutine minimum_stop(state, stop)
        type(find_minimum_result), intent(in) :: state !! active minimizer state
        logical, intent(out) :: stop !! whether to stop this lane
        stop_calls = stop_calls + 1
        stop = state%xl < -2.0_dp
    end subroutine minimum_stop

    subroutine same_roots(a, b, description)
        type(find_root_result), intent(in) :: a(:,:) !! results to compare
        type(find_root_result), intent(in) :: b(:,:) !! comparison results
        character(len=*), intent(in) :: description !! test description
        call check(all(a%status == b%status) .and. all(a%nit == b%nit) .and. &
            all(a%nfev == b%nfev) .and. all(a%success .eqv. b%success), description//' codes')
        do j = 1, n
            do i = 1, m
                if (.not. a(i,j)%success) cycle
                call check(a(i,j)%x == b(i,j)%x .and. a(i,j)%xl == b(i,j)%xl .and. &
                    a(i,j)%xr == b(i,j)%xr, description//' coordinates')
            end do
        end do
    end subroutine same_roots

    subroutine same_mins(a, b, description)
        type(find_minimum_result), intent(in) :: a(:,:) !! results to compare
        type(find_minimum_result), intent(in) :: b(:,:) !! comparison results
        character(len=*), intent(in) :: description !! test description
        call check(all(a%status == b%status) .and. all(a%nit == b%nit) .and. &
            all(a%nfev == b%nfev) .and. all(a%success .eqv. b%success), description//' codes')
        do j = 1, n
            do i = 1, m
                if (.not. a(i,j)%success) cycle
                call check(a(i,j)%x == b(i,j)%x .and. a(i,j)%xl == b(i,j)%xl .and. &
                    a(i,j)%xr == b(i,j)%xr, description//' coordinates')
            end do
        end do
    end subroutine same_mins

    subroutine same_root_brackets(a, b, description)
        type(bracket_root_result), intent(in) :: a(:,:) !! bracket outputs
        type(bracket_root_result), intent(in) :: b(:,:) !! comparison results
        character(len=*), intent(in) :: description !! test description
        call check(all(a%status == b%status) .and. all(a%nit == b%nit) .and. &
            all(a%nfev == b%nfev) .and. all(a%success .eqv. b%success), description//' codes')
        do j = 1, n
            do i = 1, m
                if (.not. a(i,j)%success) cycle
                call check(a(i,j)%xl == b(i,j)%xl .and. a(i,j)%xr == b(i,j)%xr, &
                    description//' coordinates')
            end do
        end do
    end subroutine same_root_brackets

    subroutine same_min_brackets(a, b, description)
        type(bracket_minimum_result), intent(in) :: a(:,:) !! bracket outputs
        type(bracket_minimum_result), intent(in) :: b(:,:) !! comparison results
        character(len=*), intent(in) :: description !! test description
        call check(all(a%status == b%status) .and. all(a%nit == b%nit) .and. &
            all(a%nfev == b%nfev) .and. all(a%success .eqv. b%success), description//' codes')
        do j = 1, n
            do i = 1, m
                if (.not. a(i,j)%success) cycle
                call check(a(i,j)%xl == b(i,j)%xl .and. a(i,j)%xm == b(i,j)%xm .and. &
                    a(i,j)%xr == b(i,j)%xr, description//' coordinates')
            end do
        end do
    end subroutine same_min_brackets
end program test_optimize_elementwise_vectorized_rank2
