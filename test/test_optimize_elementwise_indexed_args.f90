! SPDX-License-Identifier: MIT
! Zero-copy indexed mixed-type callbacks: broadcasting, lane packing and parity.
program test_optimize_elementwise_indexed_args
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_elementwise_batch, only : elementwise_typed_argument, &
        elementwise_typed_argument_from, elementwise_gather_typed_indices
    use scifort_optimize_root_scalar, only : find_root_vectorized_indexed_args_into, &
        bracket_root_vectorized_indexed_args_into, find_root_vectorized_args_into, &
        bracket_root_vectorized_args_into, find_root_vectorized_into, &
        find_root_result, bracket_root_result, root_status_invalid_input
    use scifort_optimize_minimize_scalar, only : find_minimum_vectorized_indexed_args_into, &
        bracket_minimum_vectorized_indexed_args_into, find_minimum_vectorized_args_into, &
        bracket_minimum_vectorized_args_into, find_minimum_result, bracket_minimum_result, &
        minimize_status_invalid_input
    implicit none

    type :: metadata
        real(dp) :: bias = 0.0_dp
        real(dp) :: scale = 1.0_dp
    end type metadata

    real(dp) :: roots(2,2,3), minima(3,2,3), xl(2,3), xm(2,3)
    real(dp) :: centers(3), expected(2,3), nan, one_bracket(2)
    real(dp) :: root14_bracket(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: start15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    integer :: shifts(2,1), i, j, errors, points, callbacks
    integer :: source_indices(3,4), selected_indices(3,6), empty_indices(3,0)
    integer :: selected_lanes(6)
    logical :: enabled(1,3)
    type(metadata) :: extras(2,3)
    type(elementwise_typed_argument) :: args(4), bad(1), single(1), empty(0)
    type(find_root_result) :: ri(2,3), rt(2,3), scalar_root
    type(find_minimum_result) :: mi(2,3), mt(2,3)
    type(bracket_root_result) :: bri(2,3), brt(2,3)
    type(bracket_minimum_result) :: bmi(2,3), bmt(2,3)
    type(find_root_result) :: high_root(1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    type(bracket_root_result) :: high_bracket(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    type(bracket_root_result) :: sliced(2,3)
    type(find_root_result), allocatable :: no_root(:,:)

    errors = 0
    ! Exercise the same gather helper used for initial repeated endpoints,
    ! later active-lane compaction and empty callback batches.
    source_indices = reshape([1,2,3, 4,5,6, 7,8,9, 10,11,12], [3,4])
    selected_lanes = [4,4,2,1,3,2]
    call elementwise_gather_typed_indices(source_indices, selected_lanes, selected_indices)
    do i = 1, size(selected_lanes)
        call check(all(selected_indices(:,i) == source_indices(:,selected_lanes(i))), &
            'repeated and shuffled typed index gather')
    end do
    call elementwise_gather_typed_indices(source_indices, [integer ::], empty_indices)
    call check(size(empty_indices,2) == 0, 'empty typed index gather')

    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    centers = [-0.5_dp, 0.2_dp, 0.6_dp]
    shifts(:,1) = [-1,1]
    enabled(1,:) = [.false.,.true.,.true.]
    do j = 1, 3
        do i = 1, 2
            extras(i,j) = metadata(0.04_dp*real(i+j,dp), 1.0_dp+real(i+j,dp))
            expected(i,j) = centers(j) + real(shifts(i,1),dp) + extras(i,j)%bias
            roots(:,i,j) = [-4.0_dp,4.0_dp]
            minima(:,i,j) = [-5.0_dp,0.0_dp,5.0_dp]
        end do
    end do
    xl = -3.0_dp
    xm = 0.0_dp
    args = [elementwise_typed_argument_from(centers), &
        elementwise_typed_argument_from(shifts), &
        elementwise_typed_argument_from(enabled), &
        elementwise_typed_argument_from(extras)]

    ! An invalid lane must not enter either callback and must not move the
    ! original indices of subsequent lanes in the active batch.
    roots(1,2,1) = nan
    minima(1,2,1) = nan
    xl(2,1) = nan
    xm(2,1) = nan

    call reset()
    call find_root_vectorized_indexed_args_into(root_indexed, roots, args, ri, maxiter=150)
    call check(ri(2,1)%status == root_status_invalid_input, 'invalid root lane')
    call check(callbacks > 1 .and. points == sum(ri%nfev), 'root packed accounting')
    call find_root_vectorized_args_into(root_typed, roots, args, rt, maxiter=150)
    do j = 1, 3
        do i = 1, 2
            call same_root(ri(i,j), rt(i,j), 'root typed/indexed parity')
            if (i == 2 .and. j == 1) cycle
            call check(ri(i,j)%success, 'root converged')
            call check(abs(ri(i,j)%x-expected(i,j)) < 1.0e-8_dp, 'root target')
        end do
    end do

    call reset()
    call find_minimum_vectorized_indexed_args_into(minimum_indexed, minima, args, mi, maxiter=150)
    call check(mi(2,1)%status == minimize_status_invalid_input, 'invalid minimum lane')
    call check(callbacks > 1 .and. points == sum(mi%nfev), 'minimum packed accounting')
    call find_minimum_vectorized_args_into(minimum_typed, minima, args, mt, maxiter=150)
    do j = 1, 3
        do i = 1, 2
            call same_minimum(mi(i,j), mt(i,j), 'minimum typed/indexed parity')
            if (i == 2 .and. j == 1) cycle
            call check(mi(i,j)%success, 'minimum converged')
            call check(abs(mi(i,j)%x-expected(i,j)) < 2.0e-5_dp, 'minimum target')
        end do
    end do

    call reset()
    call bracket_root_vectorized_indexed_args_into(root_indexed, xl, args, bri, xr0=-2.5_dp)
    call check(bri(2,1)%status == root_status_invalid_input, 'invalid root bracket lane')
    call check(callbacks > 1 .and. points == sum(bri%nfev), 'root bracket accounting')
    call bracket_root_vectorized_args_into(root_typed, xl, args, brt, xr0=-2.5_dp)
    do j = 1, 3
        do i = 1, 2
            call same_root_bracket(bri(i,j), brt(i,j), 'root bracket parity')
            if (i == 2 .and. j == 1) cycle
            call check(bri(i,j)%success, 'bracketed root')
            call check(bri(i,j)%xl <= expected(i,j) .and. &
                bri(i,j)%xr >= expected(i,j), 'root bracket contains target')
        end do
    end do

    call reset()
    call bracket_minimum_vectorized_indexed_args_into(minimum_indexed, xm, args, bmi)
    call check(bmi(2,1)%status == minimize_status_invalid_input, 'invalid minimum bracket lane')
    call check(callbacks > 1 .and. points == sum(bmi%nfev), 'minimum bracket accounting')
    call bracket_minimum_vectorized_args_into(minimum_typed, xm, args, bmt)
    do j = 1, 3
        do i = 1, 2
            call same_minimum_bracket(bmi(i,j), bmt(i,j), 'minimum bracket parity')
            if (i == 2 .and. j == 1) cycle
            call check(bmi(i,j)%success, 'bracketed minimum')
            call check(bmi(i,j)%xl < expected(i,j) .and. &
                bmi(i,j)%xr > expected(i,j), 'minimum bracket contains target')
        end do
    end do

    ! A section on output must not overwrite its neighbors. The input
    ! arguments below only vary along the remaining trailing-axis grid.
    sliced%status = 123
    xl(2,1) = -3.0_dp
    call bracket_root_vectorized_indexed_args_into(root_indexed, xl(1:1,:), &
        [elementwise_typed_argument_from(centers), &
         elementwise_typed_argument_from(-1), &
         elementwise_typed_argument_from(.true.), &
         elementwise_typed_argument_from(extras(1:1,:))], sliced(1:1,:))
    call check(all(sliced(1,:)%success), 'noncontiguous indexed output')
    call check(all(sliced(2,:)%status == 123), 'noncontiguous neighbor preserved')

    ! Rank-14 refinement and rank-15 bracketing with no polymorphic reshaping.
    root14_bracket(1,:,:,:,:,:,:,:,:,:,:,:,:,:,:) = -3.0_dp
    root14_bracket(2,:,:,:,:,:,:,:,:,:,:,:,:,:,:) = 3.0_dp
    single(1) = elementwise_typed_argument_from(0.5_dp)
    call find_root_vectorized_indexed_args_into(simple_indexed, root14_bracket, single, high_root)
    call check(all(high_root%success), 'rank14 indexed root')
    start15 = -3.0_dp
    call bracket_root_vectorized_indexed_args_into(simple_indexed, start15, single, high_bracket, &
        xr0=-2.5_dp)
    call check(all(high_bracket%success), 'rank15 indexed bracket')

    ! Invalid lists, ambiguous callback combinations and empty result grids.
    one_bracket = [-2.0_dp,2.0_dp]
    bad(1) = elementwise_typed_argument_from([1,2,3,4])
    call reset()
    call find_root_vectorized_indexed_args_into(simple_indexed, roots, bad, ri)
    call check(all(ri%status == root_status_invalid_input), 'malformed indexed broadcast')
    call find_minimum_vectorized_indexed_args_into(minimum_indexed, minima, bad, mi)
    call check(all(mi%status == minimize_status_invalid_input), 'malformed indexed minimum')
    call bracket_root_vectorized_indexed_args_into(simple_indexed, xl, bad, bri)
    call check(all(bri%status == root_status_invalid_input), 'malformed indexed root bracket')
    call bracket_minimum_vectorized_indexed_args_into(minimum_indexed, xm, bad, bmi)
    call check(all(bmi%status == minimize_status_invalid_input), 'malformed indexed minimum bracket')
    call find_root_vectorized_indexed_args_into(simple_indexed, one_bracket, empty, scalar_root)
    call check(scalar_root%status == root_status_invalid_input, 'empty typed argument list')
    call find_root_vectorized_into(bracket=one_bracket, res=scalar_root, typed_args=single, &
        f_typed=typed_simple, f_indexed=simple_indexed)
    call check(scalar_root%status == root_status_invalid_input, 'ambiguous typed callbacks rejected')
    call find_root_vectorized_into(bracket=one_bracket, res=scalar_root, typed_args=single, &
        f_indexed=simple_indexed)
    call check(scalar_root%success, 'indexed low-level dispatcher')
    allocate(no_root(0,3))
    call find_root_vectorized_indexed_args_into(root_indexed, roots(:,1:0,:), args, no_root)
    call check(size(no_root) == 0, 'zero-sized indexed result')
    call check(points == 0 .and. callbacks == 0, 'no callback for malformed/empty grids')

    if (errors /= 0) then
        print '(a,i0)', 'test_optimize_elementwise_indexed_args: failures=', errors
        error stop 1
    end if
    print '(a)', 'test_optimize_elementwise_indexed_args: PASS'
contains
    subroutine check(ok, label)
        logical, intent(in) :: ok !! assertion condition
        character(len=*), intent(in) :: label !! assertion description
        if (.not. ok) then
            errors = errors + 1
            print '(2a)', 'FAIL: ', label
        end if
    end subroutine check

    subroutine reset()
        points = 0
        callbacks = 0
    end subroutine reset

    function target(args, ids, k) result(x0)
        type(elementwise_typed_argument), intent(in) :: args(:) !! source argument list
        integer, intent(in) :: ids(:,:) !! source indices for each packed abscissa
        integer, intent(in) :: k !! packed position
        real(dp) :: x0, center, bias
        integer :: shift
        logical :: active
        center = 0.0_dp
        shift = 0
        bias = 0.0_dp
        active = .false.
        select type(v => args(1)%values(ids(1,k))%value)
        type is(real(dp))
            center = v
        class default
            call check(.false., 'center dynamic type')
        end select
        if (size(args) > 1) then
            select type(v => args(2)%values(ids(2,k))%value)
            type is(integer)
                shift = v
            class default
                call check(.false., 'offset dynamic type')
            end select
            select type(v => args(3)%values(ids(3,k))%value)
            type is(logical)
                active = v
            class default
                call check(.false., 'logical dynamic type')
            end select
            select type(v => args(4)%values(ids(4,k))%value)
            type is(metadata)
                bias = v%bias
                call check(v%scale > 0.0_dp, 'derived context preserved')
            class default
                call check(.false., 'derived context dynamic type')
            end select
        end if
        x0 = center + real(shift,dp) + bias
        if (active) x0 = x0 + 0.0_dp
    end function target

    subroutine root_indexed(x, data, ids, fx)
        real(dp), intent(in) :: x(:) !! packed evaluation points
        type(elementwise_typed_argument), intent(in) :: data(:) !! original typed arguments
        integer, intent(in) :: ids(:,:) !! source lookup indices
        real(dp), intent(out) :: fx(:) !! root objective values
        integer :: k
        callbacks = callbacks + 1
        points = points + size(x)
        call check(size(ids,1) == 4 .and. size(ids,2) == size(x), 'indexed root offsets')
        call check(size(data(1)%values) == 3 .and. &
            (size(data(2)%values) == 2 .or. size(data(2)%values) == 1) .and. &
            (size(data(4)%values) == 6 .or. size(data(4)%values) == 3), 'original boxes reused')
        do k = 1, size(x)
            fx(k) = x(k) - target(data,ids,k)
        end do
    end subroutine root_indexed

    subroutine minimum_indexed(x, data, ids, fx)
        real(dp), intent(in) :: x(:) !! packed evaluation points
        type(elementwise_typed_argument), intent(in) :: data(:) !! original typed arguments
        integer, intent(in) :: ids(:,:) !! source lookup indices
        real(dp), intent(out) :: fx(:) !! minimum objective values
        integer :: k
        real(dp) :: diff
        callbacks = callbacks + 1
        points = points + size(x)
        call check(size(ids,1) == size(data) .and. size(ids,2) == size(x), 'indexed minimum offsets')
        do k = 1, size(x)
            diff = x(k) - target(data,ids,k)
            fx(k) = diff*diff
        end do
    end subroutine minimum_indexed

    subroutine root_typed(x, data, fx)
        real(dp), intent(in) :: x(:) !! packed evaluation points
        type(elementwise_typed_argument), intent(in) :: data(:) !! copied, packed typed args
        real(dp), intent(out) :: fx(:) !! root objective values
        integer :: ids(size(data),size(x)), k
        do k = 1, size(x)
            ids(:,k) = k
            fx(k) = x(k) - target(data,ids,k)
        end do
    end subroutine root_typed

    subroutine minimum_typed(x, data, fx)
        real(dp), intent(in) :: x(:) !! packed evaluation points
        type(elementwise_typed_argument), intent(in) :: data(:) !! copied, packed typed args
        real(dp), intent(out) :: fx(:) !! minimum objective values
        integer :: ids(size(data),size(x)), k
        real(dp) :: diff
        do k = 1, size(x)
            ids(:,k) = k
            diff = x(k) - target(data,ids,k)
            fx(k) = diff*diff
        end do
    end subroutine minimum_typed

    subroutine simple_indexed(x, data, ids, fx)
        real(dp), intent(in) :: x(:) !! packed evaluation points
        type(elementwise_typed_argument), intent(in) :: data(:) !! original typed scalar
        integer, intent(in) :: ids(:,:) !! source index per abscissa
        real(dp), intent(out) :: fx(:) !! linear objective values
        integer :: k
        do k = 1, size(x)
            fx(k) = x(k)-target(data,ids,k)
        end do
    end subroutine simple_indexed

    subroutine typed_simple(x, data, fx)
        real(dp), intent(in) :: x(:) !! packed evaluation points
        type(elementwise_typed_argument), intent(in) :: data(:) !! copied typed scalars
        real(dp), intent(out) :: fx(:) !! objective values
        fx = x-0.5_dp
    end subroutine typed_simple

    subroutine same_root(a, b, label)
        type(find_root_result), intent(in) :: a, b !! indexed and copied typed states
        character(len=*), intent(in) :: label !! assertion description
        call check(a%status == b%status .and. a%nfev == b%nfev .and. a%nit == b%nit, label)
        if (a%success) call check(abs(a%x-b%x) < 1.0e-12_dp, label)
    end subroutine same_root

    subroutine same_minimum(a, b, label)
        type(find_minimum_result), intent(in) :: a, b !! indexed and copied typed states
        character(len=*), intent(in) :: label !! assertion description
        call check(a%status == b%status .and. a%nfev == b%nfev .and. a%nit == b%nit, label)
        if (a%success) call check(abs(a%x-b%x) < 1.0e-12_dp, label)
    end subroutine same_minimum

    subroutine same_root_bracket(a, b, label)
        type(bracket_root_result), intent(in) :: a, b !! indexed and copied typed brackets
        character(len=*), intent(in) :: label !! assertion description
        call check(a%status == b%status .and. a%nfev == b%nfev .and. a%nit == b%nit, label)
        if (a%success) call check(abs(a%xl-b%xl) < 1.0e-12_dp .and. &
            abs(a%xr-b%xr) < 1.0e-12_dp, label)
    end subroutine same_root_bracket

    subroutine same_minimum_bracket(a, b, label)
        type(bracket_minimum_result), intent(in) :: a, b !! indexed and copied typed brackets
        character(len=*), intent(in) :: label !! assertion description
        call check(a%status == b%status .and. a%nfev == b%nfev .and. a%nit == b%nit, label)
        if (a%success) call check(abs(a%xl-b%xl) < 1.0e-12_dp .and. &
            abs(a%xr-b%xr) < 1.0e-12_dp, label)
    end subroutine same_minimum_bracket
end program test_optimize_elementwise_indexed_args
