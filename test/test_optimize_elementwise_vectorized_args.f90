! SPDX-License-Identifier: MIT
! Objective-parameter broadcasting and packed active-lane mapping regressions.
program test_optimize_elementwise_vectorized_args
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_optimize_root_scalar, only : find_root_vectorized_args_into, &
        bracket_root_vectorized_args_into, find_root_vectorized, &
        find_root_result, bracket_root_result, root_status_invalid_input
    use scifort_optimize_minimize_scalar, only : find_minimum_vectorized_args_into, &
        bracket_minimum_vectorized_args_into, find_minimum_vectorized, &
        find_minimum_result, bracket_minimum_result, minimize_status_invalid_input
    implicit none

    real(dp) :: rb(2,2,3), mb(3,2,3), starts(2,3), mids(2,3)
    real(dp) :: args(2,1,3), full_args(2,2,3), bad_shape(2,3,3), no_args(0,1,3)
    real(dp) :: trailing_args(2,3)
    real(dp) :: rb5(2,2,1,2,1,2), args5(2,1,2,1,2)
    real(dp) :: rb14(2,2,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: start15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    type(find_root_result) :: rr(2,3), bad_res(2,2), one_root
    type(find_root_result) :: r5(2,1,2,1,2), r14(2,1,1,1,1,1,1,1,1,1,1,1,1,1)
    type(find_minimum_result) :: mr(2,3), one_min, old_min(2)
    type(bracket_root_result) :: br(2,3), b15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    type(bracket_minimum_result) :: bm(2,3)
    type(find_root_result) :: old_root(2)
    integer :: i, j, k, errors, evals, batch_calls
    real(dp) :: nan

    errors = 0
    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    args(1,1,:) = [0.2_dp, 0.8_dp, 1.4_dp]
    args(2,1,:) = [2.0_dp, 1.0_dp, 3.0_dp]
    do j = 1, 3
        do i = 1, 2
            rb(:,i,j) = [-1.0_dp, 2.0_dp]
            mb(:,i,j) = [-1.0_dp, 1.0_dp, 2.0_dp]
        end do
    end do
    starts = -1.0_dp
    mids = 1.0_dp

    ! One invalid lane forces active-point compaction. Parameter column j
    ! must continue to follow the same lane after other lanes stop.
    rb(1,2,1) = nan
    evals = 0
    batch_calls = 0
    call find_root_vectorized_args_into(root_args, rb, args, rr, maxiter=150)
    call check(rr(2,1)%status == root_status_invalid_input, 'root invalid lane status')
    call check(rr(2,1)%nfev == 0, 'root invalid lane no work')
    do j = 1, 3
        do i = 1, 2
            if (i == 2 .and. j == 1) cycle
            call check(rr(i,j)%success, 'root parameterized convergence')
            call check(abs(rr(i,j)%x - args(1,1,j)) < 1.0e-8_dp, 'root parameter lane mapping')
        end do
    end do
    call check(evals == sum(rr%nfev), 'root packed evaluation accounting')
    call check(batch_calls > 1, 'root callback invoked over iterations')

    mb(1,2,1) = nan
    evals = 0
    batch_calls = 0
    call find_minimum_vectorized_args_into(minimum_args, mb, args, mr, maxiter=150)
    call check(mr(2,1)%status == minimize_status_invalid_input, 'minimum invalid lane')
    call check(mr(2,1)%nfev == 0, 'minimum invalid lane no work')
    do j = 1, 3
        do i = 1, 2
            if (i == 2 .and. j == 1) cycle
            call check(mr(i,j)%success, 'minimum parameterized convergence')
            call check(abs(mr(i,j)%x - args(1,1,j)) < 5.0e-6_dp, 'minimum parameter lane mapping')
        end do
    end do
    call check(evals == sum(mr%nfev), 'minimum packed evaluation accounting')
    call check(batch_calls > 1, 'minimum callback invoked over iterations')

    starts(2,1) = nan
    evals = 0
    batch_calls = 0
    call bracket_root_vectorized_args_into(root_args, starts, args, br, &
        xr0=0.0_dp, factor=2.5_dp)
    call check(br(2,1)%status == root_status_invalid_input, 'bracket root invalid lane')
    do j = 1, 3
        do i = 1, 2
            if (i == 2 .and. j == 1) cycle
            call check(br(i,j)%success, 'bracket root parameterized success')
            call check(br(i,j)%xl <= args(1,1,j) .and. &
                br(i,j)%xr >= args(1,1,j), 'root bracket encloses parameter target')
        end do
    end do
    call check(evals == sum(br%nfev), 'root bracket packed evaluation accounting')
    call check(batch_calls > 1, 'root bracket active expansion')

    mids(2,1) = nan
    evals = 0
    batch_calls = 0
    call bracket_minimum_vectorized_args_into(minimum_args, mids, args, bm)
    call check(bm(2,1)%status == minimize_status_invalid_input, 'bracket min invalid lane')
    do j = 1, 3
        do i = 1, 2
            if (i == 2 .and. j == 1) cycle
            call check(bm(i,j)%success, 'bracket minimum parameterized success')
            call check(bm(i,j)%xl < args(1,1,j) .and. &
                bm(i,j)%xr > args(1,1,j), 'minimum bracket encloses parameter target')
        end do
    end do
    call check(evals == sum(bm%nfev), 'minimum bracket packed evaluation accounting')
    call check(batch_calls > 1, 'minimum bracket active expansion')

    ! Full lane-specific arguments and the shorter trailing-vector form.
    do j = 1, 3
        do i = 1, 2
            full_args(1,i,j) = 0.2_dp * real(i,dp) + 0.3_dp * real(j,dp)
            full_args(2,i,j) = 1.0_dp + real(i+j,dp)
        end do
    end do
    rb(1,:,:) = -1.0_dp
    call find_root_vectorized_args_into(root_args, rb, full_args, rr)
    call check(all(rr%success), 'full-grid parameter broadcasting')
    do j = 1, 3
        do i = 1, 2
            call check(abs(rr(i,j)%x - full_args(1,i,j)) < 1.0e-8_dp, &
                'full-grid per-lane target')
        end do
    end do
    trailing_args = args(:,1,:)
    call find_root_vectorized_args_into(root_args, rb, trailing_args, rr)
    call check(all(rr%success), 'trailing-vector parameter broadcasting')
    do j = 1, 3
        call check(all(abs(rr(:,j)%x - args(1,1,j)) < 1.0e-8_dp), &
            'trailing-vector parameter targets')
    end do

    ! A rank-one parameter stack broadcasts over a scalar result.
    call find_root_vectorized_args_into(root_args, [-1.0_dp,2.0_dp], &
        [0.75_dp,2.0_dp], one_root)
    call check(one_root%success, 'scalar parameterized root success')
    call check(abs(one_root%x - 0.75_dp) < 1.0e-8_dp, 'scalar parameterized root value')
    call find_minimum_vectorized_args_into(minimum_args, [-1.0_dp,1.0_dp,2.0_dp], &
        [0.75_dp,2.0_dp], one_min)
    call check(one_min%success, 'scalar parameterized minimum success')
    call check(abs(one_min%x - 0.75_dp) < 1.0e-6_dp, 'scalar parameterized min value')

    ! Rank-five result grid, with an implicit (trailing-axis) broadcast.
    do k = 1, 2
        do j = 1, 2
            args5(1,1,j,1,k) = 0.25_dp + real(j+k,dp) / 4.0_dp
            args5(2,1,j,1,k) = real(j+k,dp)
        end do
    end do
    rb5(1,:,:,:,:,:) = -1.0_dp
    rb5(2,:,:,:,:,:) = 2.0_dp
    call find_root_vectorized_args_into(root_args, rb5, args5, r5)
    call check(all(r5%success), 'rank-five args root')
    do k = 1, 2
        do j = 1, 2
            call check(all(abs(r5(:,1,j,1,k)%x-args5(1,1,j,1,k))<1.0e-8_dp), &
                'rank-five parameter indices')
        end do
    end do

    ! Maximum available ranks: args with only the first parameter axis.
    rb14(1,:,:,:,:,:,:,:,:,:,:,:,:,:,:) = -1.0_dp
    rb14(2,:,:,:,:,:,:,:,:,:,:,:,:,:,:) = 2.0_dp
    call find_root_vectorized_args_into(root_args, rb14, [0.75_dp,1.0_dp], r14)
    call check(all(r14%success), 'rank14 parameterized root')
    start15 = -1.0_dp
    call bracket_root_vectorized_args_into(root_args, start15, &
        [0.75_dp,1.0_dp], b15)
    call check(all(b15%success), 'rank15 parameterized bracketing')

    ! Bad parameter shapes do not evaluate the function and reset results.
    bad_shape = 0.0_dp
    evals = 0
    rr%status = 77
    call find_root_vectorized_args_into(root_args, rb, bad_shape, rr)
    call check(all(rr%status == root_status_invalid_input), 'bad args shape status')
    call check(evals == 0, 'bad args shape no callbacks')
    rr%status = 77
    call find_root_vectorized_args_into(root_args, rb, no_args, rr)
    call check(all(rr%status == root_status_invalid_input), 'zero parameter count')
    rr%status = 77
    call find_root_vectorized_args_into(root_args, rb, 1.0_dp, rr)
    call check(all(rr%status == root_status_invalid_input), 'rank-zero args invalid')
    bad_res%status = 77
    call find_root_vectorized_args_into(root_args, rb, args, bad_res)
    call check(all(bad_res%status == root_status_invalid_input), 'bad result shape')
    call check(evals == 0, 'invalid args/results never call callback')
    evals = 0
    mr%status = 77
    call find_minimum_vectorized_args_into(minimum_args, mb, bad_shape, mr)
    call check(all(mr%status == minimize_status_invalid_input), 'minimum bad arguments')
    br%status = 77
    call bracket_root_vectorized_args_into(root_args, starts, bad_shape, br)
    call check(all(br%status == root_status_invalid_input), 'bracket root bad arguments')
    bm%status = 77
    call bracket_minimum_vectorized_args_into(minimum_args, mids, bad_shape, bm)
    call check(all(bm%status == minimize_status_invalid_input), 'bracket minimum bad arguments')
    call check(evals == 0, 'invalid argument shape skips all callbacks')

    ! Noncontiguous result slices retain untouched cells.
    mr%status = 77
    call find_minimum_vectorized_args_into(minimum_args, mb(:,1:2:2,:), args, &
        mr(1:2:2,:))
    call check(all(mr(2,:)%status == 77), 'noncontiguous neighbors preserved')
    call check(all(mr(1,:)%success), 'noncontiguous computed neighbors')

    ! Rank-one APIs with their original two-argument callbacks are unchanged.
    old_root = find_root_vectorized(plain_root, reshape([-1.0_dp,2.0_dp, &
        -1.0_dp,2.0_dp],[2,2]))
    old_min = find_minimum_vectorized(plain_minimum, reshape([-1.0_dp,1.0_dp,2.0_dp, &
        -1.0_dp,1.0_dp,2.0_dp],[3,2]))
    call check(all(old_root%success), 'legacy root callback')
    call check(all(old_min%success), 'legacy minimum callback')

    if (errors /= 0) then
        print '(a,i0)', 'parameterized elementwise failures: ', errors
        error stop 1
    end if
    print '(a)', 'test_optimize_elementwise_vectorized_args: PASS'
contains
    subroutine check(ok, label)
        logical, intent(in) :: ok !! test condition
        character(len=*), intent(in) :: label !! test description
        if (.not. ok) then
            print '(2a)', 'FAIL: ', label
            errors = errors + 1
        end if
    end subroutine check

    subroutine root_args(x, a, fx)
        real(dp), intent(in) :: x(:) !! packed trial values
        real(dp), intent(in) :: a(:,:) !! first row center, second row scale
        real(dp), intent(out) :: fx(:) !! evaluated roots
        call check(size(a,1) == 2 .and. size(a,2) == size(x), 'root packed args dimensions')
        evals = evals + size(x)
        batch_calls = batch_calls + 1
        fx = a(2,:) * (x - a(1,:))
    end subroutine root_args

    subroutine minimum_args(x, a, fx)
        real(dp), intent(in) :: x(:) !! packed trial values
        real(dp), intent(in) :: a(:,:) !! first row center, second row scale
        real(dp), intent(out) :: fx(:) !! evaluated minimum
        call check(size(a,1) == 2 .and. size(a,2) == size(x), 'minimum packed args dimensions')
        evals = evals + size(x)
        batch_calls = batch_calls + 1
        fx = a(2,:) * (x - a(1,:))**2
    end subroutine minimum_args

    subroutine plain_root(x, fx)
        real(dp), intent(in) :: x(:) !! packed trial values
        real(dp), intent(out) :: fx(:) !! function values
        fx = x - 0.75_dp
    end subroutine plain_root

    subroutine plain_minimum(x, fx)
        real(dp), intent(in) :: x(:) !! packed trial values
        real(dp), intent(out) :: fx(:) !! function values
        fx = (x - 0.75_dp)**2
    end subroutine plain_minimum
end program test_optimize_elementwise_vectorized_args
