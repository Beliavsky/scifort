! SPDX-License-Identifier: MIT
! Independently broadcast objective arguments, generic compatibility, ranks,
! lane compaction, invalid parameters, and noncontiguous sections.
program test_optimize_elementwise_vectorized_arg_list
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_elementwise_batch, only : elementwise_argument, &
        elementwise_argument_from
    use scifort_optimize_root_scalar, only : find_root_vectorized_args_into, &
        bracket_root_vectorized_args_into, find_root_vectorized_into, &
        find_root_result, bracket_root_result, root_status_invalid_input
    use scifort_optimize_minimize_scalar, only : find_minimum_vectorized_args_into, &
        bracket_minimum_vectorized_args_into, find_minimum_result, &
        bracket_minimum_result, minimize_status_invalid_input
    implicit none

    real(dp) :: rb(2,2,3), mb(3,2,3), xl(2,3), xm(2,3)
    real(dp) :: center(3), offset(2,1), scale(1,3), stacked(3,2,3)
    real(dp) :: expected(2,3), nan, bad_axis(4)
    real(dp) :: rb14(2,2,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: cent15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: cent14(2,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: start15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    type(elementwise_argument) :: args(3), invalid(1), empty(0)
    type(find_root_result) :: rr(2,3), rr_stacked(2,3), rscalar
    type(find_minimum_result) :: mr(2,3), mr_stacked(2,3)
    type(bracket_root_result) :: br(2,3), br_stacked(2,3)
    type(bracket_minimum_result) :: bm(2,3), bm_stacked(2,3)
    type(find_root_result) :: root14(2,1,1,1,1,1,1,1,1,1,1,1,1,1)
    type(bracket_root_result) :: bracket15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    type(find_root_result), allocatable :: zero(:,:)
    integer :: i, j, errors, evaluations, calls

    errors = 0
    evaluations = 0
    calls = 0
    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    center = [0.20_dp, 0.65_dp, 1.15_dp]
    offset(:,1) = [0.15_dp, 0.3_dp]
    scale(1,:) = [1.0_dp, 2.5_dp, 0.8_dp]
    args = [elementwise_argument_from(center), &
            elementwise_argument_from(offset), elementwise_argument_from(scale)]
    do j = 1, 3
        do i = 1, 2
            expected(i,j) = center(j) + offset(i,1)
            stacked(:,i,j) = [center(j), offset(i,1), scale(1,j)]
            rb(:,i,j) = [-1.0_dp,2.0_dp]
            mb(:,i,j) = [-1.0_dp,1.0_dp,2.0_dp]
        end do
    end do
    xl = -1.0_dp
    xm = 1.0_dp

    ! One disabled lane forces active-lane compaction.  The three distinct
    ! argument shapes must still follow the same lane over later iterations.
    rb(1,2,1) = nan
    mb(1,2,1) = nan
    xl(2,1) = nan
    xm(2,1) = nan
    evaluations = 0
    calls = 0
    call find_root_vectorized_args_into(root_args, rb, args, rr, maxiter=150)
    call check(rr(2,1)%status == root_status_invalid_input, 'root invalid lane')
    call check(rr(2,1)%nfev == 0, 'invalid root lane no work')
    call check(evaluations == sum(rr%nfev), 'root packed nfev accounting')
    call check(calls > 1, 'root active-lane iterations')
    do j = 1, 3
        do i = 1, 2
            if (i == 2 .and. j == 1) cycle
            call check(rr(i,j)%success, 'root list success')
            call check(abs(rr(i,j)%x - expected(i,j)) < 1.0e-8_dp, 'root list targets')
        end do
    end do
    call find_root_vectorized_args_into(root_args, rb, stacked, rr_stacked, maxiter=150)
    call check(all(rr%status == rr_stacked%status), 'root stack/list status parity')
    call check(all(rr%nfev == rr_stacked%nfev), 'root stack/list counts parity')
    call check(all(rr%x == rr_stacked%x), 'root stack/list solution parity')

    evaluations = 0
    calls = 0
    call find_minimum_vectorized_args_into(minimum_args, mb, args, mr, maxiter=150)
    call check(mr(2,1)%status == minimize_status_invalid_input, 'minimum invalid lane')
    call check(evaluations == sum(mr%nfev), 'minimum packed nfev accounting')
    call check(calls > 1, 'minimum active-lane iterations')
    do j = 1, 3
        do i = 1, 2
            if (i == 2 .and. j == 1) cycle
            call check(mr(i,j)%success, 'minimum list success')
            call check(abs(mr(i,j)%x - expected(i,j)) < 5.0e-6_dp, 'minimum list targets')
        end do
    end do
    call find_minimum_vectorized_args_into(minimum_args, mb, stacked, mr_stacked, maxiter=150)
    call check(all(mr%status == mr_stacked%status), 'minimum stack/list status parity')
    call check(all(mr%nfev == mr_stacked%nfev), 'minimum stack/list counts parity')
    call check(all(mr%x == mr_stacked%x), 'minimum stack/list solution parity')

    evaluations = 0
    calls = 0
    call bracket_root_vectorized_args_into(root_args, xl, args, br, xr0=0.0_dp)
    call check(br(2,1)%status == root_status_invalid_input, 'bracket root invalid lane')
    call check(evaluations == sum(br%nfev), 'bracket root packed accounting')
    call check(calls > 1, 'bracket root active-lane iterations')
    do j = 1, 3
        do i = 1, 2
            if (i == 2 .and. j == 1) cycle
            call check(br(i,j)%success, 'bracket root list success')
            call check(br(i,j)%xl <= expected(i,j) .and. &
                br(i,j)%xr >= expected(i,j), 'bracket root contains target')
        end do
    end do
    call bracket_root_vectorized_args_into(root_args, xl, stacked, br_stacked, xr0=0.0_dp)
    call check(all(br%status == br_stacked%status), 'bracket root status parity')
    call check(all(br%nfev == br_stacked%nfev), 'bracket root counts parity')

    evaluations = 0
    calls = 0
    call bracket_minimum_vectorized_args_into(minimum_args, xm, args, bm)
    call check(bm(2,1)%status == minimize_status_invalid_input, 'bracket minimum invalid lane')
    call check(evaluations == sum(bm%nfev), 'bracket minimum packed accounting')
    call check(calls > 1, 'bracket minimum active-lane iterations')
    do j = 1, 3
        do i = 1, 2
            if (i == 2 .and. j == 1) cycle
            call check(bm(i,j)%success, 'bracket minimum list success')
            call check(bm(i,j)%xl < expected(i,j) .and. &
                bm(i,j)%xr > expected(i,j), 'bracket minimum contains target')
        end do
    end do
    call bracket_minimum_vectorized_args_into(minimum_args, xm, stacked, bm_stacked)
    call check(all(bm%status == bm_stacked%status), 'bracket min status parity')
    call check(all(bm%nfev == bm_stacked%nfev), 'bracket min counts parity')

    ! A scalar argument, a trailing vector and a rank-two grid section can
    ! coexist.  A noncontiguous result section must preserve neighboring cells.
    args(1) = elementwise_argument_from(center)
    args(2) = elementwise_argument_from(offset(1,1))
    args(3) = elementwise_argument_from(1.25_dp)
    rr%status = 197
    rb(1,2,1) = -1.0_dp
    call find_root_vectorized_args_into(root_args, rb(:,1:2:2,:), args, rr(1:2:2,:))
    call check(all(rr(1,:)%success), 'noncontiguous result computed')
    call check(all(rr(2,:)%status == 197), 'noncontiguous neighbor untouched')

    call find_root_vectorized_args_into(root_args, [-1.0_dp, 2.0_dp], &
        [elementwise_argument_from(0.65_dp), &
         elementwise_argument_from(0.10_dp), &
         elementwise_argument_from(1.0_dp)], rscalar)
    call check(rscalar%success, 'scalar grid with three scalar args')
    call check(abs(rscalar%x - 0.75_dp) < 1.0e-8_dp, 'scalar root value')

    ! Independent high-rank shapes are legal up to Fortran's rank limit.
    rb14(1,:,:,:,:,:,:,:,:,:,:,:,:,:,:) = -1.0_dp
    rb14(2,:,:,:,:,:,:,:,:,:,:,:,:,:,:) = 2.0_dp
    cent15 = 0.75_dp
    cent14 = 0.75_dp
    start15 = -1.0_dp
    call find_root_vectorized_args_into(root_args, rb14, &
        [elementwise_argument_from(cent14), &
         elementwise_argument_from(0.0_dp), elementwise_argument_from(1.0_dp)], root14)
    call check(all(root14%success), 'rank14 args list refinement')
    call bracket_root_vectorized_args_into(root_args, start15, &
        [elementwise_argument_from(cent15), elementwise_argument_from(0.0_dp), &
         elementwise_argument_from(1.0_dp)], bracket15)
    call check(all(bracket15%success), 'rank15 args list bracketing')

    ! Invalid shapes and uninitialized argument objects reject the entire
    ! call without executing the objective or leaving previous results intact.
    bad_axis = 0.0_dp
    args(1) = elementwise_argument_from(bad_axis)
    evaluations = 0
    rr%status = 197
    call find_root_vectorized_args_into(root_args, rb, args, rr)
    call check(all(rr%status == root_status_invalid_input), 'bad arg shape')
    call check(evaluations == 0, 'bad arg shape no objective calls')
    mr%status = 197
    call find_minimum_vectorized_args_into(minimum_args, mb, args, mr)
    call check(all(mr%status == minimize_status_invalid_input), 'bad min arg shape')
    br%status = 197
    call bracket_root_vectorized_args_into(root_args, xl, args, br)
    call check(all(br%status == root_status_invalid_input), 'bad root bracket arg shape')
    bm%status = 197
    call bracket_minimum_vectorized_args_into(minimum_args, xm, args, bm)
    call check(all(bm%status == minimize_status_invalid_input), 'bad min bracket arg shape')
    call check(evaluations == 0, 'all shape mismatches skip objective')

    invalid(1) = elementwise_argument()
    rr%status = 197
    call find_root_vectorized_args_into(root_args, rb, invalid, rr)
    call check(all(rr%status == root_status_invalid_input), 'default arg invalid')
    rr%status = 197
    call find_root_vectorized_args_into(root_args, rb, empty, rr)
    call check(all(rr%status == root_status_invalid_input), 'empty args list invalid')
    call check(evaluations == 0, 'uninitialized and empty args skip objective')
    rr%status = 197
    call find_root_vectorized_into(bracket=rb, res=rr, f_args=root_args, &
        args=stacked, args_list=invalid)
    call check(all(rr%status == root_status_invalid_input), 'cannot mix stack and list')
    call check(evaluations == 0, 'mixed representations skip objective')

    allocate(zero(0,3))
    call find_root_vectorized_args_into(root_args, rb(:,1:0,:), args, zero)
    call check(size(zero) == 0, 'empty grid accepted')
    call check(evaluations == 0, 'empty grid no objective calls')

    if (errors /= 0) then
        print '(a,i0)', 'arg-list regression failures: ', errors
        error stop 1
    end if
    print '(a)', 'test_optimize_elementwise_vectorized_arg_list: PASS'
contains
    subroutine check(ok, label)
        logical, intent(in) :: ok !! assertion condition
        character(len=*), intent(in) :: label !! assertion name
        if (.not. ok) then
            print '(2a)', 'FAIL: ', label
            errors = errors + 1
        end if
    end subroutine check

    subroutine root_args(x, a, fx)
        real(dp), intent(in) :: x(:) !! packed abscissae
        real(dp), intent(in) :: a(:,:) !! independently broadcast parameter rows
        real(dp), intent(out) :: fx(:) !! function values
        call check(size(a,1) == 3, 'root parameter count')
        call check(size(a,2) == size(x), 'root parameter packing')
        evaluations = evaluations + size(x)
        calls = calls + 1
        fx = a(3,:) * (x - a(1,:) - a(2,:))
    end subroutine root_args

    subroutine minimum_args(x, a, fx)
        real(dp), intent(in) :: x(:) !! packed abscissae
        real(dp), intent(in) :: a(:,:) !! independently broadcast parameter rows
        real(dp), intent(out) :: fx(:) !! function values
        call check(size(a,1) == 3, 'minimum parameter count')
        call check(size(a,2) == size(x), 'minimum parameter packing')
        evaluations = evaluations + size(x)
        calls = calls + 1
        fx = a(3,:) * (x - a(1,:) - a(2,:))**2
    end subroutine minimum_args
end program test_optimize_elementwise_vectorized_arg_list
