! SPDX-License-Identifier: MIT
! Rank-two independent elementwise optimization and NumPy-style broadcasting.
program test_optimize_elementwise_rank2
    use scifort_kinds, only : dp
    use scifort_optimize, only : bracket_root, bracket_root_result, &
        find_root, find_root_result, root_status_invalid_input, find_root_status_callback, &
        bracket_minimum, bracket_minimum_result, find_minimum, find_minimum_result, &
        minimize_status_invalid_input, find_minimum_status_callback
    use optimize_elementwise_rank2_reference
    implicit none

    type(bracket_root_result), allocatable :: br(:, :)
    type(find_root_result), allocatable :: rr(:, :)
    type(bracket_minimum_result), allocatable :: bm(:, :)
    type(find_minimum_result), allocatable :: rm(:, :)
    type(bracket_root_result) :: one_br
    type(find_root_result) :: one_rr
    type(bracket_minimum_result) :: one_bm
    type(find_minimum_result) :: one_rm
    real(dp) :: root_left(2,3), root_right(2,3), root_lower(2,3), root_upper(2,3)
    real(dp) :: root_pairs(2,2,3), minimum_middle(2,3), minimum_triples(3,2,3)
    real(dp) :: valid_root_pairs(2,2,3), valid_minimum_triples(3,2,3)
    real(dp) :: xt(3)
    integer :: calls, callback_calls, failures, i, j, batch_calls

    failures = 0
    calls = 0
    callback_calls = 0
    root_left = reshape([-4.0_dp, -4.0_dp, 0.0_dp, 0.0_dp, &
        3.0_dp, 3.0_dp], [2,3])
    root_right = reshape([-3.0_dp, -3.0_dp, 1.0_dp, 1.0_dp, &
        3.5_dp, 3.5_dp], [2,3])
    root_lower = reshape([-10.0_dp, -10.0_dp, 0.0_dp, 0.0_dp, &
        1.0_dp, 1.0_dp], [2,3])
    root_upper = reshape([-1.0_dp, -1.0_dp, 10.0_dp, 10.0_dp, &
        4.0_dp, 4.0_dp], [2,3])

    br = bracket_root(root_fun, root_left, xr0=root_right, &
        xmin=root_lower, xmax=root_upper, factor=2.0_dp, maxiter=15)
    batch_calls = calls
    call check(all(shape(br) == [2,3]), 'bracket_root rank-two output shape')
    call check(all(br%status == br_root_status_ref), 'SciPy root bracket status')
    call check(all(br%nit == br_root_nit_ref), 'SciPy root bracket iterations')
    call check(all(br%nfev == br_root_nfev_ref), 'SciPy root bracket evaluations')
    call check(all(abs(br%xl - br_root_xl_ref) < 1.0e-12_dp), 'SciPy root left endpoints')
    call check(all(abs(br%xr - br_root_xr_ref) < 1.0e-12_dp), 'SciPy root right endpoints')
    call check(batch_calls == sum(br%nfev), 'root bracket total function calls')
    do j = 1, 3
        do i = 1, 2
            one_br = bracket_root(root_fun, root_left(i,j), xr0=root_right(i,j), &
                xmin=root_lower(i,j), xmax=root_upper(i,j), factor=2.0_dp, maxiter=15)
            call compare_root_bracket(br(i,j), one_br)
        end do
    end do

    root_pairs(:,1,1) = [-4.0_dp, -1.0_dp]
    root_pairs(:,2,1) = [4.0_dp, 1.0_dp]
    root_pairs(:,1,2) = [1.0_dp, 4.0_dp]
    root_pairs(:,2,2) = [-4.0_dp, -1.0_dp]
    root_pairs(:,1,3) = [3.0_dp, 4.0_dp]
    root_pairs(:,2,3) = [1.0_dp, 4.0_dp]
    calls = 0
    rr = find_root(root_fun, root_pairs, xatol=1.0e-12_dp)
    batch_calls = calls
    call check(all(shape(rr) == [2,3]), 'find_root rank-two output shape')
    call check(all(rr%status == find_root_status_ref), 'SciPy root refinement status')
    call check(all(rr%nit == find_root_nit_ref), 'SciPy root refinement iterations')
    call check(all(rr%nfev == find_root_nfev_ref), 'SciPy root refinement evaluations')
    call check(batch_calls == sum(rr%nfev), 'root refinement total function calls')
    do j = 1, 3
        do i = 1, 2
            one_rr = find_root(root_fun, root_pairs(:,i,j), xatol=1.0e-12_dp)
            call compare_root_refinement(rr(i,j), one_rr)
            if (rr(i,j)%success) then
                call check(abs(rr(i,j)%x - find_root_x_ref(i,j)) < 1.0e-10_dp, &
                    'SciPy rank-two root estimate')
            end if
        end do
    end do
    call check(.not. rr(1,3)%success .and. rr(2,3)%success, &
        'invalid root bracket isolated from other elements')

    minimum_middle = reshape([-3.0_dp, -3.0_dp, 0.5_dp, 0.5_dp, &
        1.5_dp, 1.5_dp], [2,3])
    calls = 0
    bm = bracket_minimum(min_fun, minimum_middle, xl0=minimum_middle-0.3_dp, &
        xr0=minimum_middle+0.3_dp, xmin=-10.0_dp, xmax=10.0_dp)
    batch_calls = calls
    call check(all(shape(bm) == [2,3]), 'bracket_minimum rank-two output shape')
    call check(all(bm%status == br_min_status_ref), 'SciPy minimum bracket status')
    call check(all(bm%nit == br_min_nit_ref), 'SciPy minimum bracket iterations')
    call check(all(bm%nfev == br_min_nfev_ref), 'SciPy minimum bracket evaluations')
    call check(all(abs(bm%xm-br_min_xm_ref) < 1.0e-12_dp), 'SciPy minimum middle points')
    call check(batch_calls == sum(bm%nfev), 'minimum bracket total function calls')
    do j = 1, 3
        do i = 1, 2
            one_bm = bracket_minimum(min_fun, minimum_middle(i,j), &
                xl0=minimum_middle(i,j)-0.3_dp, xr0=minimum_middle(i,j)+0.3_dp, &
                xmin=-10.0_dp, xmax=10.0_dp)
            call compare_minimum_bracket(bm(i,j), one_bm)
        end do
    end do

    minimum_triples(:,1,1) = [0.0_dp, 1.0_dp, 3.0_dp]
    minimum_triples(:,2,1) = [3.0_dp, 1.0_dp, 0.0_dp]
    minimum_triples(:,1,2) = [-3.0_dp, 0.0_dp, 4.0_dp]
    minimum_triples(:,2,2) = [0.0_dp, 1.0_dp, 3.0_dp]
    minimum_triples(:,1,3) = [3.0_dp, 4.0_dp, 5.0_dp]
    minimum_triples(:,2,3) = [-3.0_dp, 0.0_dp, 4.0_dp]
    calls = 0
    rm = find_minimum(min_fun, minimum_triples, xatol=1.0e-9_dp)
    batch_calls = calls
    call check(all(shape(rm) == [2,3]), 'find_minimum rank-two output shape')
    call check(all(rm%status == find_min_status_ref), 'SciPy minimum refinement status')
    call check(all(rm%nit == find_min_nit_ref), 'SciPy minimum refinement iterations')
    call check(all(rm%nfev == find_min_nfev_ref), 'SciPy minimum refinement evaluations')
    call check(batch_calls == sum(rm%nfev), 'minimum refinement total function calls')
    do j = 1, 3
        do i = 1, 2
            one_rm = find_minimum(min_fun, minimum_triples(:,i,j), xatol=1.0e-9_dp)
            call compare_minimum_refinement(rm(i,j), one_rm)
            if (rm(i,j)%success) then
                call check(abs(rm(i,j)%x-find_min_x_ref(i,j)) < 1.0e-10_dp, &
                    'SciPy rank-two minimum estimate')
            end if
        end do
    end do
    call check(.not. rm(1,3)%success .and. rm(2,3)%success, &
        'invalid minimum bracket isolated from other elements')

    ! A rank-one argument broadcasts along the trailing (column) axis.
    ! Singleton rank-two axes expand independently like NumPy arrays.
    calls = 0
    br = bracket_root(root_fun, root_left, xr0=[-3.0_dp, 1.0_dp, 3.5_dp], &
        xmin=reshape([-10.0_dp, -10.0_dp], [2,1]), &
        xmax=reshape([-1.0_dp, 10.0_dp, 4.0_dp], [1,3]))
    batch_calls = calls
    call check(batch_calls == sum(br%nfev), 'broadcast root bracket evaluation count')
    do j = 1, 3
        do i = 1, 2
            one_br = bracket_root(root_fun, root_left(i,j), xr0=root_right(i,j), &
                xmin=-10.0_dp, xmax=root_upper(i,j))
            call compare_root_bracket(br(i,j), one_br)
        end do
    end do

    calls = 0
    bm = bracket_minimum(min_fun, minimum_middle, xl0=minimum_middle-0.3_dp, &
        xr0=[-2.7_dp, 0.8_dp, 1.8_dp], &
        xmin=reshape([-10.0_dp, -5.0_dp], [2,1]), xmax=10.0_dp)
    batch_calls = calls
    call check(batch_calls == sum(bm%nfev), 'broadcast minimum bracket evaluation count')
    do j = 1, 3
        do i = 1, 2
            one_bm = bracket_minimum(min_fun, minimum_middle(i,j), &
                xl0=minimum_middle(i,j)-0.3_dp, &
                xr0=minimum_middle(i,j)+0.3_dp, &
                xmin=merge(-10.0_dp, -5.0_dp, i == 1), xmax=10.0_dp)
            call compare_minimum_bracket(bm(i,j), one_bm)
        end do
    end do

    xt = [1.0e-12_dp, 1.0e-8_dp, 1.0e-10_dp]
    calls = 0
    rr = find_root(root_fun, root_pairs, xatol=xt, &
        xrtol=reshape([0.0_dp, 0.0_dp], [2,1]), &
        fatol=0.0_dp, frtol=reshape([0.0_dp, 0.0_dp, 0.0_dp], [1,3]))
    batch_calls = calls
    call check(batch_calls == sum(rr%nfev), 'broadcast root tolerances evaluation count')
    do j = 1, 3
        do i = 1, 2
            one_rr = find_root(root_fun, root_pairs(:,i,j), &
                xatol=xt(j), xrtol=0.0_dp, fatol=0.0_dp, frtol=0.0_dp)
            call compare_root_refinement(rr(i,j), one_rr)
        end do
    end do

    calls = 0
    rm = find_minimum(min_fun, minimum_triples, &
        xatol=reshape([1.0e-9_dp, 1.0e-6_dp, 1.0e-8_dp], [1,3]), &
        xrtol=reshape([0.0_dp, 0.0_dp], [2,1]), &
        fatol=0.0_dp, frtol=0.0_dp)
    batch_calls = calls
    call check(batch_calls == sum(rm%nfev), 'broadcast minimum tolerances evaluation count')
    do j = 1, 3
        do i = 1, 2
            xt = [1.0e-9_dp, 1.0e-6_dp, 1.0e-8_dp]
            one_rm = find_minimum(min_fun, minimum_triples(:,i,j), &
                xatol=xt(j), xrtol=0.0_dp, fatol=0.0_dp, frtol=0.0_dp)
            call compare_minimum_refinement(rm(i,j), one_rm)
        end do
    end do

    ! Length-one vectors and (1,1) matrices also broadcast to the grid.
    calls = 0
    br = bracket_root(root_fun, root_left, xr0=[1.0_dp], &
        xmin=reshape([-10.0_dp], [1,1]), xmax=20.0_dp)
    batch_calls = calls
    call check(batch_calls == sum(br%nfev), 'singleton root broadcast evaluation count')
    do j = 1, 3
        do i = 1, 2
            one_br = bracket_root(root_fun, root_left(i,j), xr0=1.0_dp, &
                xmin=-10.0_dp, xmax=20.0_dp)
            call compare_root_bracket(br(i,j), one_br)
        end do
    end do

    calls = 0
    rr = find_root(root_fun, root_pairs, xatol=[1.0e-12_dp], &
        xrtol=reshape([0.0_dp], [1,1]))
    batch_calls = calls
    call check(batch_calls == sum(rr%nfev), 'singleton root tolerance evaluation count')
    do j = 1, 3
        do i = 1, 2
            one_rr = find_root(root_fun, root_pairs(:,i,j), xatol=1.0e-12_dp, xrtol=0.0_dp)
            call compare_root_refinement(rr(i,j), one_rr)
        end do
    end do

    ! Each callback receives a scalar state, including an initial state.
    valid_root_pairs(1,:,:) = -4.0_dp
    valid_root_pairs(2,:,:) = -1.0_dp
    callback_calls = 0
    rr = find_root(root_fun, valid_root_pairs, callback=root_callback)
    call check(all(rr%status == find_root_status_callback), 'rank-two root callback stops each')
    call check(all(rr%nit == 1), 'rank-two root callback iteration')
    call check(callback_calls == 12, 'rank-two root callback count')
    valid_minimum_triples(1,:,:) = 0.0_dp
    valid_minimum_triples(2,:,:) = 1.0_dp
    valid_minimum_triples(3,:,:) = 3.0_dp
    callback_calls = 0
    rm = find_minimum(min_fun, valid_minimum_triples, callback=minimum_callback)
    call check(all(rm%status == find_minimum_status_callback), 'rank-two minimum callback stops each')
    call check(all(rm%nit == 1), 'rank-two minimum callback iteration')
    call check(callback_calls == 12, 'rank-two minimum callback count')

    calls = 0
    br = bracket_root(root_fun, root_left, xr0=[-3.0_dp, -2.0_dp])
    call check(all(br%status == root_status_invalid_input), 'root incorrect trailing-axis size')
    bm = bracket_minimum(min_fun, minimum_middle, xmin=reshape([0.0_dp, 1.0_dp, 2.0_dp], [3,1]))
    call check(all(bm%status == minimize_status_invalid_input), 'minimum incorrect matrix size')
    rr = find_root(root_fun, root_pairs, xatol=reshape([0.0_dp, 1.0_dp], [1,1,2]))
    call check(all(rr%status == root_status_invalid_input), 'root rank-three tolerance invalid')
    rm = find_minimum(min_fun, minimum_triples, xatol=[1.0_dp, 2.0_dp])
    call check(all(rm%status == minimize_status_invalid_input), 'minimum wrong vector tolerance')
    rr = find_root(root_fun, reshape([(real(i,dp), i=1,18)], [3,2,3]))
    call check(all(rr%status == root_status_invalid_input), 'root incorrect bracket leading size')
    rm = find_minimum(min_fun, reshape([(real(i,dp), i=1,12)], [2,2,3]))
    call check(all(rm%status == minimize_status_invalid_input), 'minimum incorrect bracket leading size')
    call check(calls == 0, 'no objective evaluations for invalid shapes')

    br = bracket_root(root_fun, root_left(:,1:0))
    bm = bracket_minimum(min_fun, minimum_middle(1:0,:))
    rr = find_root(root_fun, root_pairs(:,:,1:0))
    rm = find_minimum(min_fun, minimum_triples(:,1:0,:))
    call check(all(shape(br) == [2,0]) .and. all(shape(bm) == [0,3]), &
        'empty bracketing grids preserve shape')
    call check(all(shape(rr) == [2,0]) .and. all(shape(rm) == [0,3]), &
        'empty refinement grids preserve shape')
    call check(calls == 0, 'no objective evaluations for empty grids')

    if (failures > 0) then
        print '(a,i0)', 'test_optimize_elementwise_rank2: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_elementwise_rank2: PASS'

contains

    function root_fun(x) result(fx)
        real(dp), intent(in) :: x !! independent root-function argument
        real(dp) :: fx !! corresponding function value
        calls = calls + 1
        fx = x*x - 4.0_dp
    end function root_fun

    function min_fun(x) result(fx)
        real(dp), intent(in) :: x !! independent minimization argument
        real(dp) :: fx !! objective value
        calls = calls + 1
        fx = (x - 1.5_dp)**2
    end function min_fun

    subroutine root_callback(state, stop)
        type(find_root_result), intent(in) :: state !! current scalar state
        logical, intent(out) :: stop !! true after the first iteration
        callback_calls = callback_calls + 1
        stop = state%nit >= 1
    end subroutine root_callback

    subroutine minimum_callback(state, stop)
        type(find_minimum_result), intent(in) :: state !! current scalar state
        logical, intent(out) :: stop !! true after the first iteration
        callback_calls = callback_calls + 1
        stop = state%nit >= 1
    end subroutine minimum_callback

    subroutine check(condition, message)
        logical, intent(in) :: condition !! whether a regression condition held
        character(len=*), intent(in) :: message !! regression case label
        if (.not. condition) then
            failures = failures + 1
            print '(a)', 'FAIL: '//message
        end if
    end subroutine check

    subroutine compare_root_bracket(a, b)
        type(bracket_root_result), intent(in) :: a !! rank-two grid element
        type(bracket_root_result), intent(in) :: b !! scalar reference result
        call check(a%status == b%status .and. (a%success .eqv. b%success), 'root bracket status')
        call check(a%nit == b%nit .and. a%nfev == b%nfev, 'root bracket counts')
        call check(all(abs([a%xl-b%xl, a%xr-b%xr, a%fl-b%fl, a%fr-b%fr]) < 1.0e-12_dp), &
            'root bracket data')
    end subroutine compare_root_bracket

    subroutine compare_root_refinement(a, b)
        type(find_root_result), intent(in) :: a !! rank-two grid element
        type(find_root_result), intent(in) :: b !! scalar reference result
        call check(a%status == b%status .and. (a%success .eqv. b%success), 'root refinement status')
        call check(a%nit == b%nit .and. a%nfev == b%nfev, 'root refinement counts')
        if (a%success) then
            call check(abs(a%x-b%x) < 1.0e-12_dp, 'root refinement estimate')
        end if
    end subroutine compare_root_refinement

    subroutine compare_minimum_bracket(a, b)
        type(bracket_minimum_result), intent(in) :: a !! rank-two grid element
        type(bracket_minimum_result), intent(in) :: b !! scalar reference result
        call check(a%status == b%status .and. (a%success .eqv. b%success), 'minimum bracket status')
        call check(a%nit == b%nit .and. a%nfev == b%nfev, 'minimum bracket counts')
        call check(all(abs([a%xl-b%xl, a%xm-b%xm, a%xr-b%xr]) < 1.0e-12_dp), &
            'minimum bracket positions')
    end subroutine compare_minimum_bracket

    subroutine compare_minimum_refinement(a, b)
        type(find_minimum_result), intent(in) :: a !! rank-two grid element
        type(find_minimum_result), intent(in) :: b !! scalar reference result
        call check(a%status == b%status .and. (a%success .eqv. b%success), 'minimum refinement status')
        call check(a%nit == b%nit .and. a%nfev == b%nfev, 'minimum refinement counts')
        if (a%success) then
            call check(abs(a%x-b%x) < 1.0e-12_dp, 'minimum refinement estimate')
        end if
    end subroutine compare_minimum_refinement
end program test_optimize_elementwise_rank2
