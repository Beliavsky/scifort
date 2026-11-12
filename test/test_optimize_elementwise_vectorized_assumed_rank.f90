! SPDX-License-Identifier: MIT
! Assumed-rank packed callbacks, including rank-0, rank-5, rank-14/15,
! empty/noncontiguous grids, broadcasting, and input validation.
program test_optimize_elementwise_vectorized_assumed_rank
    use scifort_optimize_root_scalar, only : find_root_vectorized, &
        find_root_vectorized_into, bracket_root_vectorized, bracket_root_vectorized_into, &
        find_root_result, bracket_root_result, root_status_invalid_input
    use scifort_optimize_minimize_scalar, only : find_minimum_vectorized, &
        find_minimum_vectorized_into, bracket_minimum_vectorized, bracket_minimum_vectorized_into, &
        find_minimum_result, bracket_minimum_result, minimize_status_invalid_input
    use scifort_kinds, only : dp
    implicit none

    integer, parameter :: n = 8
    real(dp) :: rb5(2,2,1,2,1,2), mb5(3,2,1,2,1,2)
    real(dp) :: xl5(2,1,2,1,2), xm5(2,1,2,1,2)
    real(dp) :: rb_bad(3,2,1,2,1,2)
    real(dp) :: br14(2,2,1,1,1,1,1,1,1,1,1,1,1,1,1), bm14(3,2,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: starts15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: bad_tol(3), tol(2)
    real(dp) :: factor5(1,1,2,1,2), factor_full(2,1,2,1,2)
    type(find_root_result) :: rr5(2,1,2,1,2), rr_ref(n), rscalar
    type(find_root_result) :: rr14(2,1,1,1,1,1,1,1,1,1,1,1,1,1), rr_bad(1,1,2,1,2)
    type(find_minimum_result) :: mr5(2,1,2,1,2), mr_ref(n), mscalar
    type(find_minimum_result) :: mr14(2,1,1,1,1,1,1,1,1,1,1,1,1,1)
    type(bracket_root_result) :: br5(2,1,2,1,2), br_ref(n), brscalar
    type(bracket_root_result) :: br15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    type(bracket_minimum_result) :: bm5(2,1,2,1,2), bm_ref(n), bmscalar
    type(bracket_minimum_result) :: bm15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    type(find_root_result), allocatable :: empty(:,:,:,:,:)
    integer :: i, j, k, calls, fails

    fails = 0
    calls = 0
    tol = [1.0e-12_dp, 1.0e-11_dp]
    bad_tol = 1.0e-12_dp
    do k = 1, 2
        do j = 1, 2
            factor5(1,1,j,1,k) = 2.0_dp + real(j+k,dp)
            do i = 1, 2
                factor_full(i,1,j,1,k) = factor5(1,1,j,1,k)
                rb5(:,i,1,j,1,k) = [-1.0_dp,2.0_dp]
                mb5(:,i,1,j,1,k) = [-1.0_dp,0.75_dp,2.0_dp]
                xl5(i,1,j,1,k) = -1.0_dp
                xm5(i,1,j,1,k) = 0.75_dp
            end do
        end do
    end do

    ! Rank-five grids: all four must be identical to flattened rank-one APIs.
    calls = 0
    call find_root_vectorized_into(root_batch, rb5, rr5, xatol=tol)
    call check(calls == sum(rr5%nfev), 'rank5 packed root work accounting')
    rr_ref = find_root_vectorized(root_batch, reshape(rb5,[2,n]), &
        xatol=reshape(spread(tol,1,4),[n]))
    call check(all(rr5%success), 'rank5 root convergence')
    call check(all(reshape(rr5%nfev,[n]) == rr_ref%nfev), 'rank5 root counts')
    call check(maxval(abs(reshape(rr5%x,[n])-rr_ref%x)) < 1.0e-12_dp, 'rank5 root parity')

    call find_minimum_vectorized_into(min_batch, mb5, mr5, xatol=tol)
    mr_ref = find_minimum_vectorized(min_batch, reshape(mb5,[3,n]), &
        xatol=reshape(spread(tol,1,4),[n]))
    call check(all(mr5%success), 'rank5 minimum convergence')
    call check(all(reshape(mr5%nfev,[n]) == mr_ref%nfev), 'rank5 minimum counts')
    call check(maxval(abs(reshape(mr5%x,[n])-mr_ref%x)) < 1.0e-12_dp, 'rank5 minimum parity')

    call bracket_root_vectorized_into(root_batch, xl5, br5, xmin=-5.0_dp, factor=factor5)
    br_ref = bracket_root_vectorized(root_batch, reshape(xl5,[n]), &
        xmin=-5.0_dp, factor=reshape(factor_full,[n]))
    call check(all(br5%success), 'rank5 bracket root success')
    call check(all(reshape(br5%nfev,[n]) == br_ref%nfev), 'rank5 bracket root counts')
    call check(all(reshape(br5%xl,[n]) == br_ref%xl), 'rank5 bracket root parity')

    call bracket_minimum_vectorized_into(min_batch, xm5, bm5, factor=factor5)
    bm_ref = bracket_minimum_vectorized(min_batch, reshape(xm5,[n]), &
        factor=reshape(factor_full,[n]))
    call check(all(bm5%success), 'rank5 bracket minimum success')
    call check(all(reshape(bm5%nfev,[n]) == bm_ref%nfev), 'rank5 bracket minimum counts')
    call check(all(reshape(bm5%xm,[n]) == bm_ref%xm), 'rank5 bracket minimum parity')

    ! Scalar (rank zero) output works without a new scalar wrapper.
    call find_root_vectorized_into(root_batch, [-1.0_dp,2.0_dp], rscalar)
    call check(rscalar%success .and. abs(rscalar%x-0.75_dp)<1.0e-8_dp, 'rank0 root')
    call find_minimum_vectorized_into(min_batch, [-1.0_dp,0.75_dp,2.0_dp], mscalar)
    call check(mscalar%success .and. abs(mscalar%x-0.75_dp)<1.0e-6_dp, 'rank0 minimum')
    call bracket_root_vectorized_into(root_batch, -1.0_dp, brscalar)
    call check(brscalar%success, 'rank0 bracket root')
    call bracket_minimum_vectorized_into(min_batch, 0.75_dp, bmscalar)
    call check(bmscalar%success, 'rank0 bracket minimum')

    ! Maximum Fortran 2018 ranks: 15 for starting points, 14 for bracketed
    ! refinement because the initial bracket occupies one leading dimension.
    br14 = 0.0_dp
    br14(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1) = -1.0_dp
    br14(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1) = 2.0_dp
    br14(1,2,1,1,1,1,1,1,1,1,1,1,1,1,1) = -1.0_dp
    br14(2,2,1,1,1,1,1,1,1,1,1,1,1,1,1) = 2.0_dp
    bm14 = 0.0_dp
    bm14(1,:,:,:,:,:,:,:,:,:,:,:,:,:,:) = -1.0_dp
    bm14(2,:,:,:,:,:,:,:,:,:,:,:,:,:,:) = 0.75_dp
    bm14(3,:,:,:,:,:,:,:,:,:,:,:,:,:,:) = 2.0_dp
    starts15 = -1.0_dp
    call find_root_vectorized_into(root_batch, br14, rr14)
    call check(all(rr14%success), 'rank14 root refinement')
    call find_minimum_vectorized_into(min_batch, bm14, mr14)
    call check(all(mr14%success), 'rank14 minimum refinement')
    call bracket_root_vectorized_into(root_batch, starts15, br15)
    call check(all(br15%success), 'rank15 root bracketing')
    starts15 = 0.75_dp
    call bracket_minimum_vectorized_into(min_batch, starts15, bm15)
    call check(all(bm15%success), 'rank15 minimum bracketing')

    ! Invalid shapes and parameters are diagnosed without invoking callback,
    ! and previous results are cleared to invalid (INOUT workaround).
    rr_bad%status = 101
    calls = 0
    call find_root_vectorized_into(root_batch, rb5, rr_bad)
    call check(calls == 0, 'mismatched output no calls')
    call check(all(rr_bad%status == root_status_invalid_input), 'mismatched output status')
    rr5%status = 101
    calls = 0
    call find_root_vectorized_into(root_batch, rb_bad, rr5)
    call check(calls == 0, 'wrong bracket leading axis no calls')
    call check(all(rr5%status == root_status_invalid_input), 'wrong bracket axis status')
    rr5%status = 101
    calls = 0
    call find_root_vectorized_into(root_batch, rb5, rr5, xatol=bad_tol)
    call check(calls == 0, 'incompatible trailing broadcast no calls')
    call check(all(rr5%status == root_status_invalid_input), 'incompatible broadcast status')

    ! A single bad high-rank growth factor invalidates only its own lane.
    factor_full = 2.0_dp
    factor_full(1,1,1,1,1) = 1.0_dp
    call bracket_minimum_vectorized_into(min_batch, xm5, bm5, factor=factor_full)
    call check(bm5(1,1,1,1,1)%status == minimize_status_invalid_input, &
        'invalid factor isolated to lane')
    call check(all(bm5(2,:,:,:,:)%success), 'valid factor lanes unaffected')

    ! Noncontiguous input and output sections preserve unaffected lanes.
    rr5%status = 101
    calls = 0
    call find_root_vectorized_into(root_batch, rb5(:,1:2:2,:,:,:,:), &
        rr5(1:2:2,:,:,:,:))
    call check(all(rr5(1:2:2,:,:,:,:)%success), 'noncontiguous output computed')
    call check(all(rr5(2:2:2,:,:,:,:)%status == 101), 'noncontiguous neighbors unaffected')

    ! Empty grids preserve shape and perform zero objective evaluations.
    allocate(empty(0,1,2,1,2))
    calls = 0
    call find_root_vectorized_into(root_batch, rb5(:,1:0,:,:,:,:), empty)
    call check(calls == 0 .and. size(empty)==0, 'empty grid')

    if (fails > 0) then
        print '(a,i0)', 'assumed-rank vectorized failures: ', fails
        error stop 1
    end if
    print '(a)', 'test_optimize_elementwise_vectorized_assumed_rank: PASS'
contains
    subroutine check(ok, message)
        logical, intent(in) :: ok !! true if test passed
        character(len=*), intent(in) :: message !! diagnostic
        if (.not. ok) then
            print '(2a)', 'FAIL: ', message
            fails = fails + 1
        end if
    end subroutine check

    subroutine root_batch(x, fx)
        real(dp), intent(in) :: x(:) !! active abscissae
        real(dp), intent(out) :: fx(:) !! function values
        calls = calls + size(x)
        fx = x - 0.75_dp
    end subroutine root_batch

    subroutine min_batch(x, fx)
        real(dp), intent(in) :: x(:) !! active abscissae
        real(dp), intent(out) :: fx(:) !! function values
        calls = calls + size(x)
        fx = (x - 0.75_dp)**2
    end subroutine min_batch
end program test_optimize_elementwise_vectorized_assumed_rank
