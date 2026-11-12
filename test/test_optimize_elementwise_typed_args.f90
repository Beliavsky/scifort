! SPDX-License-Identifier: MIT
! Typed objective arguments: exact dynamic types, independent broadcast axes,
! lane compaction, invalid input, scalar/high-rank and noncontiguous grids.
program test_optimize_elementwise_typed_args
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_optimize_elementwise_batch, only : elementwise_typed_argument, &
        elementwise_typed_argument_from
    use scifort_optimize_root_scalar, only : find_root_vectorized_args_into, &
        find_root_vectorized_into, &
        bracket_root_vectorized_args_into, find_root_result, bracket_root_result, &
        root_status_invalid_input
    use scifort_optimize_minimize_scalar, only : find_minimum_vectorized_args_into, &
        bracket_minimum_vectorized_args_into, find_minimum_result, &
        bracket_minimum_result, minimize_status_invalid_input
    implicit none

    type :: lane_metadata
        integer :: label = -1
        real(dp) :: bias = 0.0_dp
        real(dp) :: gain = 1.0_dp
    end type lane_metadata

    real(dp) :: rb(2,2,3), mb(3,2,3), xl(2,3), xm(2,3)
    real(dp) :: center(3), expected(2,3), rb_scalar(2)
    real(dp) :: rb14(2,2,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: target14(2,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: start15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    integer :: shift(2,1)
    logical :: negative(1,3)
    type(lane_metadata) :: metadata(2,3)
    type(lane_metadata) :: metadata15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    type(elementwise_typed_argument) :: params(4), sliced(4), broken(1), empty(0)
    type(find_root_result) :: rr(2,3), rs, root14(2,1,1,1,1,1,1,1,1,1,1,1,1,1)
    type(bracket_root_result) :: br(2,3), bracket15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    type(find_minimum_result) :: mr(2,3)
    type(bracket_minimum_result) :: bm(2,3)
    type(find_root_result), allocatable :: none(:,:)
    type(find_root_result) :: labelled(2)
    real(dp) :: labels_bracket(2,2)
    integer :: i, j, errors, evaluations, callbacks

    errors = 0
    evaluations = 0
    callbacks = 0
    center = [-0.85_dp, 0.15_dp, 0.65_dp]
    shift(:,1) = [-1,1]
    negative(1,:) = [.false.,.true.,.false.]
    do j = 1, 3
        do i = 1, 2
            metadata(i,j) = lane_metadata(10*j+i, 0.11_dp*i+0.03_dp*j, &
                1.0_dp+0.2_dp*i+0.1_dp*j)
            expected(i,j) = center(j)+real(shift(i,1),dp)+metadata(i,j)%bias
            rb(:,i,j) = [-3.0_dp, 3.0_dp]
            mb(:,i,j) = [-6.0_dp, 0.0_dp, 6.0_dp]
        end do
    end do
    xl = -3.0_dp
    xm = 0.0_dp
    params = [elementwise_typed_argument_from(center), &
        elementwise_typed_argument_from(shift), &
        elementwise_typed_argument_from(negative), &
        elementwise_typed_argument_from(metadata)]

    ! Invalid input lane forces packing of only the surviving problems.
    rb(1,2,1) = ieee_value(0.0_dp, ieee_quiet_nan)
    mb(1,2,1) = ieee_value(0.0_dp, ieee_quiet_nan)
    xl(2,1) = ieee_value(0.0_dp, ieee_quiet_nan)
    xm(2,1) = ieee_value(0.0_dp, ieee_quiet_nan)

    evaluations = 0
    callbacks = 0
    call find_root_vectorized_args_into(root_typed, rb, params, rr, maxiter=150)
    call check(rr(2,1)%status == root_status_invalid_input, 'root invalid lane')
    call check(rr(2,1)%nfev == 0, 'root invalid lane no evaluations')
    call check(callbacks > 1, 'root callback lane compaction')
    call check(evaluations == sum(rr%nfev), 'root function counts')
    do j = 1, 3
        do i = 1, 2
            if (i == 2 .and. j == 1) cycle
            call check(rr(i,j)%success, 'typed root success')
            call check(abs(rr(i,j)%x-expected(i,j)) < 1.0e-8_dp, 'typed root target')
        end do
    end do

    evaluations = 0
    callbacks = 0
    call find_minimum_vectorized_args_into(minimum_typed, mb, params, mr, maxiter=150)
    call check(mr(2,1)%status == minimize_status_invalid_input, 'minimum invalid lane')
    call check(callbacks > 1, 'minimum callback lane compaction')
    call check(evaluations == sum(mr%nfev), 'minimum function counts')
    do j = 1, 3
        do i = 1, 2
            if (i == 2 .and. j == 1) cycle
            call check(mr(i,j)%success, 'typed minimum success')
            call check(abs(mr(i,j)%x-expected(i,j)) < 2.0e-5_dp, 'typed minimum target')
        end do
    end do

    evaluations = 0
    callbacks = 0
    call bracket_root_vectorized_args_into(root_typed, xl, params, br, xr0=-2.5_dp)
    call check(br(2,1)%status == root_status_invalid_input, 'root bracket invalid lane')
    call check(callbacks > 1, 'root bracket expansion')
    call check(evaluations == sum(br%nfev), 'root bracket counts')
    do j = 1, 3
        do i = 1, 2
            if (i == 2 .and. j == 1) cycle
            call check(br(i,j)%success, 'typed root bracket success')
            call check(br(i,j)%xl <= expected(i,j) .and. &
                br(i,j)%xr >= expected(i,j), 'typed root bracket contains root')
        end do
    end do

    evaluations = 0
    callbacks = 0
    call bracket_minimum_vectorized_args_into(minimum_typed, xm, params, bm)
    call check(bm(2,1)%status == minimize_status_invalid_input, 'minimum bracket invalid lane')
    call check(callbacks > 1, 'minimum bracket expansion')
    call check(evaluations == sum(bm%nfev), 'minimum bracket counts')
    do j = 1, 3
        do i = 1, 2
            if (i == 2 .and. j == 1) cycle
            call check(bm(i,j)%success, 'typed minimum bracket success')
            call check(bm(i,j)%xl < expected(i,j) .and. &
                bm(i,j)%xr > expected(i,j), 'typed minimum bracket contains minimizer')
        end do
    end do

    ! Noncontiguous output sections must not overwrite neighboring results.
    rb(1,2,1) = -3.0_dp
    sliced = [elementwise_typed_argument_from(center), &
        elementwise_typed_argument_from(shift(1,1)), &
        elementwise_typed_argument_from(negative), &
        elementwise_typed_argument_from(metadata(1:1,:))]
    rr%status = 197
    call find_root_vectorized_args_into(root_typed, rb(:,1:2:2,:), sliced, rr(1:2:2,:))
    call check(all(rr(1,:)%success), 'noncontiguous typed output success')
    call check(all(rr(2,:)%status == 197), 'noncontiguous neighbor untouched')

    ! Scalar grids with independently typed scalar arguments.
    rb_scalar = [-3.0_dp, 3.0_dp]
    call find_root_vectorized_args_into(root_typed, rb_scalar, &
        [elementwise_typed_argument_from(-0.25_dp), &
         elementwise_typed_argument_from(1), &
         elementwise_typed_argument_from(.true.), &
         elementwise_typed_argument_from(lane_metadata(42,0.0_dp,1.0_dp))], rs)
    call check(rs%success, 'scalar typed grid success')
    call check(abs(rs%x-0.75_dp) < 1.0e-8_dp, 'scalar typed grid target')

    ! Maximum ranks, including a derived-type input of rank 15.
    rb14(1,:,:,:,:,:,:,:,:,:,:,:,:,:,:) = -3.0_dp
    rb14(2,:,:,:,:,:,:,:,:,:,:,:,:,:,:) = 3.0_dp
    target14 = 0.9_dp
    start15 = -3.0_dp
    metadata15 = lane_metadata(1,0.9_dp,1.0_dp)
    call find_root_vectorized_args_into(simple_typed, rb14, &
        [elementwise_typed_argument_from(target14)], root14)
    call check(all(root14%success), 'rank14 typed root')
    call bracket_root_vectorized_args_into(simple_typed, start15, &
        [elementwise_typed_argument_from(metadata15)], bracket15, xr0=3.0_dp)
    call check(all(bracket15%success), 'rank15 derived-type bracket')

    ! Invalid shapes, empty arrays and unfinished argument objects must fail
    ! before evaluating the objective, with initialized invalid statuses.
    evaluations = 0
    broken(1) = elementwise_typed_argument_from([1,2,3,4])
    rr%status = 197
    call find_root_vectorized_args_into(simple_typed, rb, broken, rr)
    call check(all(rr%status == root_status_invalid_input), 'bad typed broadcast shape')
    mr%status = 197
    call find_minimum_vectorized_args_into(minimum_typed, mb, broken, mr)
    call check(all(mr%status == minimize_status_invalid_input), 'bad min broadcast shape')
    br%status = 197
    call bracket_root_vectorized_args_into(simple_typed, xl, broken, br)
    call check(all(br%status == root_status_invalid_input), 'bad root bracket shape')
    bm%status = 197
    call bracket_minimum_vectorized_args_into(minimum_typed, xm, broken, bm)
    call check(all(bm%status == minimize_status_invalid_input), 'bad min bracket shape')
    broken(1) = elementwise_typed_argument()
    call find_root_vectorized_args_into(simple_typed, rb, broken, rr)
    call check(all(rr%status == root_status_invalid_input), 'default typed object invalid')
    broken(1) = elementwise_typed_argument_from(2)
    deallocate(broken(1)%values)
    call find_root_vectorized_args_into(simple_typed, rb, broken, rr)
    call check(all(rr%status == root_status_invalid_input), 'malformed typed object invalid')
    broken(1) = elementwise_typed_argument_from(2)
    deallocate(broken(1)%values(1)%value)
    call find_root_vectorized_args_into(simple_typed, rb, broken, rr)
    call check(all(rr%status == root_status_invalid_input), 'unallocated payload invalid')
    call find_root_vectorized_args_into(simple_typed, rb, empty, rr)
    call check(all(rr%status == root_status_invalid_input), 'empty typed list invalid')
    call check(evaluations == 0, 'invalid typed arguments skip callback')
    allocate(none(0,3))
    call find_root_vectorized_args_into(root_typed, rb(:,1:0,:), params, none)
    call check(size(none) == 0, 'zero-sized typed grid accepted')
    call check(evaluations == 0, 'zero-sized typed grid no callbacks')

    ! Character data also retain their exact type and per-lane values.
    labels_bracket(:,1) = [-3.0_dp,3.0_dp]
    labels_bracket(:,2) = [-3.0_dp,3.0_dp]
    call find_root_vectorized_args_into(labelled_root, labels_bracket, &
        [elementwise_typed_argument_from(['low ', 'high'])], labelled)
    call check(all(labelled%success), 'character typed roots succeed')
    call check(abs(labelled(1)%x-0.5_dp) < 1.0e-8_dp, 'character low root')
    call check(abs(labelled(2)%x-1.5_dp) < 1.0e-8_dp, 'character high root')

    ! Heterogeneous values are not forced through real(dp): an int64 that
    ! cannot be represented exactly as double precision remains exact.
    call find_root_vectorized_args_into(large_integer_complex, rb_scalar, &
        [elementwise_typed_argument_from(9007199254740993_int64), &
         elementwise_typed_argument_from(cmplx(1.0_dp,-2.0_dp,kind=dp))], rs)
    call check(rs%success, 'integer64 and complex scalar success')
    call check(abs(rs%x-1.0_dp) < 1.0e-8_dp, 'integer64 and complex root')

    ! The lower-level dispatcher accepts f_typed plus typed_args directly.
    call find_root_vectorized_into(bracket=rb_scalar, res=rs, &
        f_typed=simple_typed, &
        typed_args=[elementwise_typed_argument_from(0.4_dp)])
    call check(rs%success, 'direct typed dispatcher success')
    call check(abs(rs%x-0.4_dp) < 1.0e-8_dp, 'direct typed dispatcher target')
    call find_root_vectorized_into(bracket=rb_scalar, res=rs, &
        typed_args=[elementwise_typed_argument_from(0.4_dp)])
    call check(rs%status == root_status_invalid_input, 'typed dispatcher requires callback')


    if (errors /= 0) then
        print '(a,i0)', 'typed argument regression failures: ', errors
        error stop 1
    end if
    print '(a)', 'test_optimize_elementwise_typed_args: PASS'
contains
    subroutine check(ok, label)
        logical, intent(in) :: ok !! assertion condition
        character(len=*), intent(in) :: label !! assertion message
        if (.not. ok) then
            print '(2a)', 'FAIL: ', label
            errors = errors + 1
        end if
    end subroutine check

    subroutine extract_factors(args, target, slope)
        type(elementwise_typed_argument), intent(in) :: args(:) !! active typed parameters
        real(dp), intent(out) :: target(:) !! root/minimum location per lane
        real(dp), intent(out) :: slope(:) !! signed objective slope per lane
        integer :: k

        call check(size(args) == 4, 'typed callback four arguments')
        do k = 1, size(target)
            select type(v => args(1)%values(k)%value)
            type is (real(dp))
                target(k) = v
            class default
                error stop 'expected real(dp) argument'
            end select
            select type(v => args(2)%values(k)%value)
            type is (integer)
                target(k) = target(k) + real(v, dp)
            class default
                error stop 'expected integer argument'
            end select
            select type(v => args(3)%values(k)%value)
            type is (logical)
                slope(k) = 1.0_dp
                if (v) slope(k) = -1.0_dp
            class default
                error stop 'expected logical argument'
            end select
            select type(v => args(4)%values(k)%value)
            type is (lane_metadata)
                target(k) = target(k) + v%bias
                slope(k) = slope(k) * v%gain
                call check(v%label > 0, 'derived type label')
            class default
                error stop 'expected derived-type argument'
            end select
        end do
    end subroutine extract_factors

    subroutine root_typed(x, args, fx)
        real(dp), intent(in) :: x(:) !! packed abscissae
        type(elementwise_typed_argument), intent(in) :: args(:) !! packed heterogeneous arguments
        real(dp), intent(out) :: fx(:) !! evaluated roots
        real(dp) :: target(size(x)), slope(size(x))
        call extract_factors(args, target, slope)
        evaluations = evaluations + size(x)
        callbacks = callbacks + 1
        fx = slope * (x-target)
    end subroutine root_typed

    subroutine minimum_typed(x, args, fx)
        real(dp), intent(in) :: x(:) !! packed abscissae
        type(elementwise_typed_argument), intent(in) :: args(:) !! packed heterogeneous arguments
        real(dp), intent(out) :: fx(:) !! evaluated minima
        real(dp) :: target(size(x)), slope(size(x))
        call extract_factors(args, target, slope)
        evaluations = evaluations + size(x)
        callbacks = callbacks + 1
        fx = (slope * (x-target))**2
    end subroutine minimum_typed

    subroutine simple_typed(x, args, fx)
        real(dp), intent(in) :: x(:) !! packed abscissae
        type(elementwise_typed_argument), intent(in) :: args(:) !! typed parameters
        real(dp), intent(out) :: fx(:) !! returned values
        integer :: k

        call check(size(args) == 1, 'simple typed argument count')
        evaluations = evaluations + size(x)
        do k = 1, size(x)
            select type(v => args(1)%values(k)%value)
            type is (real(dp))
                fx(k) = x(k) - v
            type is (lane_metadata)
                fx(k) = x(k) - v%bias
            class default
                error stop 'unsupported simple typed parameter'
            end select
        end do
    end subroutine simple_typed

    subroutine labelled_root(x, args, fx)
        real(dp), intent(in) :: x(:) !! packed abscissae
        type(elementwise_typed_argument), intent(in) :: args(:) !! character parameters
        real(dp), intent(out) :: fx(:) !! returned values
        integer :: k

        do k = 1, size(x)
            select type(v => args(1)%values(k)%value)
            type is (character(len=*))
                if (v == 'low ') then
                    fx(k) = x(k)-0.5_dp
                else if (v == 'high') then
                    fx(k) = x(k)-1.5_dp
                else
                    error stop 'unexpected label'
                end if
            class default
                error stop 'expected character argument'
            end select
        end do
    end subroutine labelled_root

    subroutine large_integer_complex(x, args, fx)
        real(dp), intent(in) :: x(:) !! packed abscissae
        type(elementwise_typed_argument), intent(in) :: args(:) !! typed inputs
        real(dp), intent(out) :: fx(:) !! results
        integer :: k

        do k = 1, size(x)
            select type(v => args(1)%values(k)%value)
            type is (integer(int64))
                call check(v == 9007199254740993_int64, 'int64 not rounded to real(dp)')
            class default
                error stop 'expected int64'
            end select
            select type(v => args(2)%values(k)%value)
            type is (complex(dp))
                call check(abs(v-cmplx(1.0_dp,-2.0_dp,kind=dp)) < &
                    1.0e-14_dp, 'complex dynamic type retained')
            class default
                error stop 'expected complex(dp)'
            end select
            fx(k) = x(k) - 1.0_dp
        end do
    end subroutine large_integer_complex
end program test_optimize_elementwise_typed_args
