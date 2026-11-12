! SPDX-License-Identifier: MIT
! High-rank packed vector callbacks: broadcasting and rank-one equivalence.
program test_optimize_elementwise_vectorized_rank34
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_root_scalar, only : find_root_vectorized, bracket_root_vectorized, &
        find_root_result, bracket_root_result, root_status_invalid_input, &
        find_root_status_callback
    use scifort_optimize_minimize_scalar, only : find_minimum_vectorized, &
        bracket_minimum_vectorized, find_minimum_result, bracket_minimum_result, &
        minimize_status_invalid_input, find_minimum_status_callback
    implicit none

    integer, parameter :: m = 2, n = 2, k = 3, p = 2
    integer, parameter :: count3 = m*n*k, count4 = m*n*k*p
    real(dp), parameter :: root_cases(2,6) = reshape([ &
        0.0_dp, 2.0_dp, -1.0_dp, 1.2_dp, 1.0_dp, 2.0_dp, &
        3.0_dp, 4.0_dp, -2.0_dp, 0.0_dp, 0.5_dp, 2.0_dp], [2,6])
    real(dp), parameter :: min_cases(3,6) = reshape([ &
        0.0_dp, 1.0_dp, 3.0_dp, -1.0_dp, 2.0_dp, 6.0_dp, &
        -4.0_dp, 1.0_dp, 4.0_dp, 6.0_dp, 4.0_dp, 3.0_dp, &
        3.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.5_dp, 3.0_dp], [3,6])
    real(dp), parameter :: root_starts(6) = [-3.0_dp,0.0_dp,5.0_dp,1.0_dp,6.0_dp,0.5_dp]
    real(dp), parameter :: min_starts(6) = [0.0_dp,1.5_dp,5.0_dp,-4.0_dp,9.5_dp,7.0_dp]
    real(dp) :: rb3(2,m,n,k), mb3(3,m,n,k), xl3(m,n,k), xm3(m,n,k)
    real(dp) :: rb4(2,m,n,k,p), mb4(3,m,n,k,p), xl4(m,n,k,p), xm4(m,n,k,p)
    real(dp) :: tol_last(k), tol3(m,n,k), tol4(1,n,1,p), tol_full4(m,n,k,p)
    real(dp) :: fatol3(m,1,1), fatol_full3(m,n,k)
    real(dp) :: factor3(1,n,k), factor_full3(m,n,k)
    real(dp) :: factor4(1,n,1,p), factor_full4(m,n,k,p)
    real(dp) :: lower3(m,1,1), lower_full3(m,n,k)
    real(dp) :: lower4(m,1,1,1), lower_full4(m,n,k,p)
    real(dp) :: right3(1,n,k), right_full3(m,n,k)
    real(dp) :: bad_factor4(m,n,k,p), bad_tol3(m,n,k), nan
    type(find_root_result) :: roots3(m,n,k), roots4(m,n,k,p), roots_one3(count3)
    type(find_root_result) :: roots_one4(count4)
    type(find_minimum_result) :: mins3(m,n,k), mins4(m,n,k,p), mins_one3(count3)
    type(find_minimum_result) :: mins_one4(count4)
    type(bracket_root_result) :: br3(m,n,k), br4(m,n,k,p), br_one3(count3)
    type(bracket_root_result) :: br_one4(count4)
    type(bracket_minimum_result) :: bm3(m,n,k), bm4(m,n,k,p), bm_one3(count3)
    type(bracket_minimum_result) :: bm_one4(count4)
    type(find_root_result), allocatable :: empty_r3(:,:,:), empty_r4(:,:,:,:)
    type(find_minimum_result), allocatable :: empty_m3(:,:,:), empty_m4(:,:,:,:)
    type(bracket_root_result), allocatable :: empty_br3(:,:,:), empty_br4(:,:,:,:)
    type(bracket_minimum_result), allocatable :: empty_bm3(:,:,:), empty_bm4(:,:,:,:)
    integer :: i, j, h, q, t, failures, calls, points, first_batch, stop_calls

    failures = 0
    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    tol_last = [1.0e-12_dp, 1.0e-9_dp, 1.0e-11_dp]
    fatol3(:,1,1) = [0.0_dp,1.0e-15_dp]
    lower3(:,1,1) = [-6.0_dp,-8.0_dp]
    lower4(:,1,1,1) = lower3(:,1,1)
    do h = 1, k
        do j = 1, n
            factor3(1,j,h) = 2.0_dp + real(mod(j+h, 2),dp)
            right3(1,j,h) = 0.5_dp + real(j+h,dp) / 4.0_dp
            do i = 1, m
                t = 1 + mod(i + 2*j + 3*h, 6)
                rb3(:,i,j,h) = root_cases(:,t)
                mb3(:,i,j,h) = min_cases(:,t)
                xl3(i,j,h) = root_starts(t)
                xm3(i,j,h) = min_starts(t)
                tol3(i,j,h) = tol_last(h)
                fatol_full3(i,j,h) = fatol3(i,1,1)
                factor_full3(i,j,h) = factor3(1,j,h)
                right_full3(i,j,h) = right3(1,j,h)
                lower_full3(i,j,h) = lower3(i,1,1)
            end do
        end do
    end do
    do q = 1, p
        do j = 1, n
            tol4(1,j,1,q) = 1.0e-12_dp * real(2**(j+q),dp)
            factor4(1,j,1,q) = 2.0_dp + real(mod(j+q,2),dp)
        end do
        do h = 1, k
            do j = 1, n
                do i = 1, m
                    t = i
                    if (q == 2) t = m - i + 1
                    rb4(:,i,j,h,q) = rb3(:,t,j,h)
                    mb4(:,i,j,h,q) = mb3(:,t,j,h)
                    xl4(i,j,h,q) = xl3(t,j,h)
                    xm4(i,j,h,q) = xm3(t,j,h)
                    tol_full4(i,j,h,q) = tol4(1,j,1,q)
                    factor_full4(i,j,h,q) = factor4(1,j,1,q)
                    lower_full4(i,j,h,q) = lower4(i,1,1,1)
                end do
            end do
        end do
    end do

    ! Rank-three: all four algorithms must match a flattened packed batch.
    call reset()
    roots3 = find_root_vectorized(root_function, rb3, xatol=tol_last, &
        xrtol=1.0e-12_dp, fatol=0.0_dp)
    call check(first_batch == 2*count3, 'rank3 root endpoints packed')
    call check(points == sum(roots3%nfev), 'rank3 root work accounting')
    call reset()
    roots_one3 = find_root_vectorized(root_function, reshape(rb3,[2,count3]), &
        xatol=reshape(tol3,[count3]), xrtol=1.0e-12_dp, fatol=0.0_dp)
    call same_roots(reshape(roots3,[count3]), roots_one3, 'rank3 root parity')

    call reset()
    mins3 = find_minimum_vectorized(minimum_function, mb3, xatol=tol_last, fatol=fatol3)
    call check(first_batch == 3*count3, 'rank3 minimum triples packed')
    call check(points == sum(mins3%nfev), 'rank3 minimum work accounting')
    call reset()
    mins_one3 = find_minimum_vectorized(minimum_function, reshape(mb3,[3,count3]), &
        xatol=reshape(tol3,[count3]), fatol=reshape(fatol_full3,[count3]))
    call same_mins(reshape(mins3,[count3]), mins_one3, 'rank3 minimum parity')

    call reset()
    br3 = bracket_root_vectorized(root_function, xl3, xmin=lower3, &
        xmax=10.0_dp, factor=factor3)
    call check(first_batch == 2*count3, 'rank3 root brackets packed')
    call check(points == sum(br3%nfev), 'rank3 root bracket accounting')
    call reset()
    br_one3 = bracket_root_vectorized(root_function, reshape(xl3,[count3]), &
        xmin=reshape(lower_full3,[count3]), xmax=10.0_dp, &
        factor=reshape(factor_full3,[count3]))
    call same_br_roots(reshape(br3,[count3]), br_one3, 'rank3 root bracket parity')

    call reset()
    bm3 = bracket_minimum_vectorized(minimum_function, xm3, xmin=-10.0_dp, &
        xmax=10.0_dp, factor=factor3)
    call check(first_batch == 3*count3, 'rank3 minimum brackets packed')
    call check(points == sum(bm3%nfev), 'rank3 minimum bracket accounting')
    call reset()
    bm_one3 = bracket_minimum_vectorized(minimum_function, reshape(xm3,[count3]), &
        xmin=-10.0_dp, xmax=10.0_dp, factor=reshape(factor_full3,[count3]))
    call same_br_mins(reshape(bm3,[count3]), bm_one3, 'rank3 minimum bracket parity')

    ! Rank-four: singleton axes not adjacent to the trailing axis broadcast.
    call reset()
    roots4 = find_root_vectorized(root_function, rb4, xatol=tol4, fatol=0.0_dp)
    call check(first_batch == 2*count4, 'rank4 root endpoints packed')
    call check(points == sum(roots4%nfev), 'rank4 root work accounting')
    call reset()
    roots_one4 = find_root_vectorized(root_function, reshape(rb4,[2,count4]), &
        xatol=reshape(tol_full4,[count4]), fatol=0.0_dp)
    call same_roots(reshape(roots4,[count4]), roots_one4, 'rank4 root parity')

    call reset()
    mins4 = find_minimum_vectorized(minimum_function, mb4, xatol=tol4)
    call check(first_batch == 3*count4, 'rank4 minimum triples packed')
    call check(points == sum(mins4%nfev), 'rank4 minimum work accounting')
    call reset()
    mins_one4 = find_minimum_vectorized(minimum_function, reshape(mb4,[3,count4]), &
        xatol=reshape(tol_full4,[count4]))
    call same_mins(reshape(mins4,[count4]), mins_one4, 'rank4 minimum parity')

    call reset()
    br4 = bracket_root_vectorized(root_function, xl4, xmin=lower4, &
        xmax=10.0_dp, factor=factor4)
    call check(first_batch == 2*count4, 'rank4 root brackets packed')
    call check(points == sum(br4%nfev), 'rank4 root bracket accounting')
    call reset()
    br_one4 = bracket_root_vectorized(root_function, reshape(xl4,[count4]), &
        xmin=reshape(lower_full4,[count4]), xmax=10.0_dp, &
        factor=reshape(factor_full4,[count4]))
    call same_br_roots(reshape(br4,[count4]), br_one4, 'rank4 root bracket parity')

    call reset()
    bm4 = bracket_minimum_vectorized(minimum_function, xm4, xmin=-10.0_dp, &
        xmax=10.0_dp, factor=factor4)
    call check(first_batch == 3*count4, 'rank4 minimum brackets packed')
    call check(points == sum(bm4%nfev), 'rank4 minimum bracket accounting')
    call reset()
    bm_one4 = bracket_minimum_vectorized(minimum_function, reshape(xm4,[count4]), &
        xmin=-10.0_dp, xmax=10.0_dp, factor=reshape(factor_full4,[count4]))
    call same_br_mins(reshape(bm4,[count4]), bm_one4, 'rank4 minimum bracket parity')

    ! Present optional starting endpoints, in addition to absent defaults.
    call reset()
    br3 = bracket_root_vectorized(root_function, xl3, xr0=right3, factor=2.0_dp)
    call reset()
    br_one3 = bracket_root_vectorized(root_function, reshape(xl3,[count3]), &
        xr0=reshape(right_full3,[count3]), factor=2.0_dp)
    call same_br_roots(reshape(br3,[count3]), br_one3, 'rank3 explicit right endpoints')
    call reset()
    bm4 = bracket_minimum_vectorized(minimum_function, xm4, &
        xl0=xm4-0.25_dp, xr0=xm4+0.3_dp, factor=2.0_dp)
    call reset()
    bm_one4 = bracket_minimum_vectorized(minimum_function, reshape(xm4,[count4]), &
        xl0=reshape(xm4-0.25_dp,[count4]), xr0=reshape(xm4+0.3_dp,[count4]), &
        factor=2.0_dp)
    call same_br_mins(reshape(bm4,[count4]), bm_one4, 'rank4 explicit endpoints')

    ! Optional per-lane callbacks must not lose their grid location or counts.
    call reset()
    stop_calls = 0
    roots3 = find_root_vectorized(root_function, rb3, callback=root_stop)
    call check(stop_calls > 0 .and. any(roots3%status == find_root_status_callback), &
        'rank3 root callback stop')
    call check(any(roots3%success) .and. points == sum(roots3%nfev), &
        'rank3 other root lanes complete')
    call reset()
    stop_calls = 0
    mins4 = find_minimum_vectorized(minimum_function, mb4, callback=minimum_stop)
    call check(stop_calls > 0 .and. any(mins4%status == find_minimum_status_callback), &
        'rank4 minimum callback stop')
    call check(any(mins4%success) .and. points == sum(mins4%nfev), &
        'rank4 other minimum lanes complete')

    ! Invalid shapes reject the complete grid *before* any callback invocation.
    call reset()
    roots3 = find_root_vectorized(root_function, rb3, xatol=[1.0_dp,2.0_dp])
    mins3 = find_minimum_vectorized(minimum_function, mb3, fatol=[1.0_dp,2.0_dp])
    br3 = bracket_root_vectorized(root_function, xl3, factor=[2.0_dp,3.0_dp])
    bm3 = bracket_minimum_vectorized(minimum_function, xm3, xr0=[2.0_dp,3.0_dp])
    call check(all(roots3%status == root_status_invalid_input) .and. &
        all(mins3%status == minimize_status_invalid_input) .and. &
        all(br3%status == root_status_invalid_input) .and. &
        all(bm3%status == minimize_status_invalid_input) .and. calls == 0, &
        'rank3 incompatible broadcast shapes')
    roots4 = find_root_vectorized(root_function, rb4, xatol=[1.0_dp,2.0_dp,3.0_dp])
    mins4 = find_minimum_vectorized(minimum_function, mb4, frtol=[1.0_dp,2.0_dp,3.0_dp])
    br4 = bracket_root_vectorized(root_function, xl4, factor=[2.0_dp,3.0_dp,4.0_dp])
    bm4 = bracket_minimum_vectorized(minimum_function, xm4, xmin=[0.0_dp,1.0_dp,2.0_dp])
    call check(all(roots4%status == root_status_invalid_input) .and. &
        all(mins4%status == minimize_status_invalid_input) .and. &
        all(br4%status == root_status_invalid_input) .and. &
        all(bm4%status == minimize_status_invalid_input) .and. calls == 0, &
        'rank4 incompatible broadcast shapes')

    ! A nonfinite parameter invalidates exactly one lane under IEEE traps.
    bad_tol3 = tol3
    bad_tol3(2,1,2) = nan
    call reset()
    roots3 = find_root_vectorized(root_function, rb3, xatol=bad_tol3)
    call check(roots3(2,1,2)%status == root_status_invalid_input .and. &
        roots3(2,1,2)%nfev == 0 .and. first_batch == 2*(count3-1), &
        'rank3 nonfinite tolerance stays in its lane')
    bad_factor4 = factor_full4
    bad_factor4(2,1,1,2) = nan
    call reset()
    bm4 = bracket_minimum_vectorized(minimum_function, xm4, factor=bad_factor4)
    call check(bm4(2,1,1,2)%status == minimize_status_invalid_input .and. &
        bm4(2,1,1,2)%nfev == 0 .and. first_batch == 3*(count4-1), &
        'rank4 nonfinite factor stays in its lane')

    call reset()
    roots3 = find_root_vectorized(root_function, &
        reshape([1.0_dp],[1,m,n,k],pad=[1.0_dp]))
    mins4 = find_minimum_vectorized(minimum_function, &
        reshape([1.0_dp],[2,m,n,k,p],pad=[1.0_dp]))
    call check(all(roots3%status == root_status_invalid_input) .and. &
        all(mins4%status == minimize_status_invalid_input) .and. calls == 0, &
        'invalid leading bracket extents')

    empty_r3 = find_root_vectorized(root_function, reshape([real(dp)::],[2,0,n,k]))
    empty_m3 = find_minimum_vectorized(minimum_function, reshape([real(dp)::],[3,m,0,k]))
    empty_br3 = bracket_root_vectorized(root_function, reshape([real(dp)::],[m,0,k]))
    empty_bm3 = bracket_minimum_vectorized(minimum_function, reshape([real(dp)::],[0,n,k]))
    empty_r4 = find_root_vectorized(root_function, reshape([real(dp)::],[2,m,n,0,p]))
    empty_m4 = find_minimum_vectorized(minimum_function, reshape([real(dp)::],[3,0,n,k,p]))
    empty_br4 = bracket_root_vectorized(root_function, reshape([real(dp)::],[m,n,k,0]))
    empty_bm4 = bracket_minimum_vectorized(minimum_function, reshape([real(dp)::],[m,0,k,p]))
    call check(all(shape(empty_r3) == [0,n,k]) .and. &
        all(shape(empty_m3) == [m,0,k]) .and. &
        all(shape(empty_br3) == [m,0,k]) .and. &
        all(shape(empty_bm3) == [0,n,k]) .and. &
        all(shape(empty_r4) == [m,n,0,p]) .and. &
        all(shape(empty_m4) == [0,n,k,p]) .and. &
        all(shape(empty_br4) == [m,n,k,0]) .and. &
        all(shape(empty_bm4) == [m,0,k,p]) .and. calls == 0, &
        'zero-size grids keep all original dimensions')

    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_elementwise_vectorized_rank34: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_elementwise_vectorized_rank34: PASS'
contains
    subroutine reset()
        calls = 0
        points = 0
        first_batch = 0
    end subroutine reset

    subroutine check(ok, description)
        logical, intent(in) :: ok !! regression predicate
        character(len=*), intent(in) :: description !! failure description
        if (.not. ok) then
            print '(a)', 'FAIL: '//description
            failures = failures + 1
        end if
    end subroutine check

    subroutine root_function(x, fx)
        real(dp), intent(in) :: x(:) !! packed points
        real(dp), intent(out) :: fx(:) !! residuals at packed points
        calls = calls + 1
        points = points + size(x)
        if (calls == 1) first_batch = size(x)
        fx = x**3 - 1.0_dp
    end subroutine root_function

    subroutine minimum_function(x, fx)
        real(dp), intent(in) :: x(:) !! packed points
        real(dp), intent(out) :: fx(:) !! objective at packed points
        calls = calls + 1
        points = points + size(x)
        if (calls == 1) first_batch = size(x)
        fx = (x-1.5_dp)**2 + 0.02_dp*(x-1.5_dp)**4
    end subroutine minimum_function

    subroutine root_stop(state, stop)
        type(find_root_result), intent(in) :: state !! current root bracket
        logical, intent(out) :: stop !! stop this individual problem
        stop_calls = stop_calls + 1
        stop = state%xl < -0.5_dp .and. state%xr > 1.0_dp
    end subroutine root_stop

    subroutine minimum_stop(state, stop)
        type(find_minimum_result), intent(in) :: state !! current minimum bracket
        logical, intent(out) :: stop !! stop this individual problem
        stop_calls = stop_calls + 1
        stop = state%xl < -2.0_dp
    end subroutine minimum_stop

    subroutine same_roots(a, b, description)
        type(find_root_result), intent(in) :: a(:) !! rank-three or -four flattened results
        type(find_root_result), intent(in) :: b(:) !! direct rank-one results
        character(len=*), intent(in) :: description !! comparison label
        integer :: ix
        call check(all(a%status == b%status) .and. all(a%nit == b%nit) .and. &
            all(a%nfev == b%nfev) .and. all(a%success .eqv. b%success), description//' status')
        do ix = 1, size(a)
            if (.not. a(ix)%success) cycle
            call check(a(ix)%x == b(ix)%x .and. a(ix)%f_x == b(ix)%f_x .and. &
                a(ix)%xl == b(ix)%xl .and. a(ix)%xr == b(ix)%xr, description//' values')
        end do
    end subroutine same_roots

    subroutine same_mins(a, b, description)
        type(find_minimum_result), intent(in) :: a(:) !! rank-three or -four flattened results
        type(find_minimum_result), intent(in) :: b(:) !! direct rank-one results
        character(len=*), intent(in) :: description !! comparison label
        integer :: ix
        call check(all(a%status == b%status) .and. all(a%nit == b%nit) .and. &
            all(a%nfev == b%nfev) .and. all(a%success .eqv. b%success), description//' status')
        do ix = 1, size(a)
            if (.not. a(ix)%success) cycle
            call check(a(ix)%x == b(ix)%x .and. a(ix)%f_x == b(ix)%f_x .and. &
                a(ix)%xl == b(ix)%xl .and. a(ix)%xm == b(ix)%xm .and. &
                a(ix)%xr == b(ix)%xr, description//' values')
        end do
    end subroutine same_mins

    subroutine same_br_roots(a, b, description)
        type(bracket_root_result), intent(in) :: a(:) !! flattened grid results
        type(bracket_root_result), intent(in) :: b(:) !! direct rank-one results
        character(len=*), intent(in) :: description !! comparison label
        integer :: ix
        call check(all(a%status == b%status) .and. all(a%nit == b%nit) .and. &
            all(a%nfev == b%nfev) .and. all(a%success .eqv. b%success), description//' status')
        do ix = 1, size(a)
            if (.not. a(ix)%success) cycle
            call check(a(ix)%xl == b(ix)%xl .and. a(ix)%xr == b(ix)%xr .and. &
                a(ix)%fl == b(ix)%fl .and. a(ix)%fr == b(ix)%fr, description//' values')
        end do
    end subroutine same_br_roots

    subroutine same_br_mins(a, b, description)
        type(bracket_minimum_result), intent(in) :: a(:) !! flattened grid results
        type(bracket_minimum_result), intent(in) :: b(:) !! direct rank-one results
        character(len=*), intent(in) :: description !! comparison label
        integer :: ix
        call check(all(a%status == b%status) .and. all(a%nit == b%nit) .and. &
            all(a%nfev == b%nfev) .and. all(a%success .eqv. b%success), description//' status')
        do ix = 1, size(a)
            if (.not. a(ix)%success) cycle
            call check(a(ix)%xl == b(ix)%xl .and. a(ix)%xm == b(ix)%xm .and. &
                a(ix)%xr == b(ix)%xr .and. a(ix)%fl == b(ix)%fl .and. &
                a(ix)%fm == b(ix)%fm .and. a(ix)%fr == b(ix)%fr, description//' values')
        end do
    end subroutine same_br_mins
end program test_optimize_elementwise_vectorized_rank34
