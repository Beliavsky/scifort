! SPDX-License-Identifier: MIT
! Independent rank-three problem grids with NumPy-style trailing broadcasting.
program test_optimize_elementwise_rank3
    use scifort_kinds, only : dp
    use scifort_optimize, only : bracket_root, bracket_root_result, &
        find_root, find_root_result, root_status_invalid_input, find_root_status_callback, &
        bracket_minimum, bracket_minimum_result, find_minimum, find_minimum_result, &
        minimize_status_invalid_input, find_minimum_status_callback
    use optimize_elementwise_rank3_reference
    implicit none

    type(bracket_root_result), allocatable :: br(:, :, :)
    type(find_root_result), allocatable :: rr(:, :, :)
    type(bracket_minimum_result), allocatable :: bm(:, :, :)
    type(find_minimum_result), allocatable :: rm(:, :, :)
    type(bracket_root_result) :: one_br
    type(find_root_result) :: one_rr
    type(bracket_minimum_result) :: one_bm
    type(find_minimum_result) :: one_rm
    real(dp) :: roots_valid(2,2,3,2), minimum_valid(3,2,3,2)
    real(dp) :: tol_matrix(3,2), upper_tensor(2,1,2)
    real(dp) :: zero_bound(1,3,1)
    integer :: calls, callback_calls, failures, i, j, k, batch_calls

    failures = 0
    calls = 0
    callback_calls = 0

    br = bracket_root(root_fun, root_left_ref, xr0=root_right_ref, &
        xmin=root_lower_ref, xmax=root_upper_ref, factor=2.0_dp, maxiter=15)
    batch_calls = calls
    call check(all(shape(br) == [2,3,2]), 'rank-three bracket_root shape')
    call check(all(br%status == br_root_status_ref), 'SciPy rank-three root bracket status')
    call check(all(br%nit == br_root_nit_ref), 'SciPy rank-three root bracket iterations')
    call check(all(br%nfev == br_root_nfev_ref), 'SciPy rank-three root bracket evaluations')
    call check(all(abs(br%xl-br_root_xl_ref) < 1.0e-12_dp), 'SciPy rank-three left root brackets')
    call check(all(abs(br%xr-br_root_xr_ref) < 1.0e-12_dp), 'SciPy rank-three right root brackets')
    call check(batch_calls == sum(br%nfev), 'rank-three root bracket call accounting')
    do k = 1, 2
        do j = 1, 3
            do i = 1, 2
                one_br = bracket_root(root_fun, root_left_ref(i,j,k), &
                    xr0=root_right_ref(i,j,k), xmin=root_lower_ref(i,j,k), &
                    xmax=root_upper_ref(i,j,k), factor=2.0_dp, maxiter=15)
                call compare_root_bracket(br(i,j,k), one_br)
            end do
        end do
    end do

    calls = 0
    rr = find_root(root_fun, root_pairs_ref, xatol=1.0e-12_dp)
    batch_calls = calls
    call check(all(shape(rr) == [2,3,2]), 'rank-three find_root shape')
    call check(all(rr%status == find_root_status_ref), 'SciPy rank-three find_root status')
    call check(all(rr%nit == find_root_nit_ref), 'SciPy rank-three find_root iterations')
    call check(all(rr%nfev == find_root_nfev_ref), 'SciPy rank-three find_root evaluations')
    call check(batch_calls == sum(rr%nfev), 'rank-three find_root call accounting')
    do k = 1, 2
        do j = 1, 3
            do i = 1, 2
                one_rr = find_root(root_fun, root_pairs_ref(:,i,j,k), xatol=1.0e-12_dp)
                call compare_root_refinement(rr(i,j,k), one_rr)
                if (rr(i,j,k)%success) then
                    call check(abs(rr(i,j,k)%x-find_root_x_ref(i,j,k)) < 1.0e-10_dp, &
                        'SciPy rank-three refined root value')
                end if
            end do
        end do
    end do
    call check(any(rr%status /= 0) .and. any(rr%status == 0), &
        'rank-three root failures do not affect successful elements')

    calls = 0
    bm = bracket_minimum(min_fun, minimum_middle_ref, xl0=minimum_middle_ref-0.3_dp, &
        xr0=minimum_middle_ref+0.3_dp, xmin=-10.0_dp, xmax=10.0_dp)
    batch_calls = calls
    call check(all(shape(bm) == [2,3,2]), 'rank-three bracket_minimum shape')
    call check(all(bm%status == br_min_status_ref), 'SciPy rank-three minimum bracket status')
    call check(all(bm%nit == br_min_nit_ref), 'SciPy rank-three minimum bracket iterations')
    call check(all(bm%nfev == br_min_nfev_ref), 'SciPy rank-three minimum bracket evaluations')
    call check(all(abs(bm%xm-br_min_xm_ref) < 1.0e-12_dp), 'SciPy rank-three middle minima')
    call check(batch_calls == sum(bm%nfev), 'rank-three minimum bracket call accounting')
    do k = 1, 2
        do j = 1, 3
            do i = 1, 2
                one_bm = bracket_minimum(min_fun, minimum_middle_ref(i,j,k), &
                    xl0=minimum_middle_ref(i,j,k)-0.3_dp, &
                    xr0=minimum_middle_ref(i,j,k)+0.3_dp, xmin=-10.0_dp, xmax=10.0_dp)
                call compare_minimum_bracket(bm(i,j,k), one_bm)
            end do
        end do
    end do

    calls = 0
    rm = find_minimum(min_fun, minimum_triples_ref, xatol=1.0e-9_dp)
    batch_calls = calls
    call check(all(shape(rm) == [2,3,2]), 'rank-three find_minimum shape')
    call check(all(rm%status == find_min_status_ref), 'SciPy rank-three minimum status')
    call check(all(rm%nit == find_min_nit_ref), 'SciPy rank-three minimum iterations')
    call check(all(rm%nfev == find_min_nfev_ref), 'SciPy rank-three minimum evaluations')
    call check(batch_calls == sum(rm%nfev), 'rank-three minimum refinement call accounting')
    do k = 1, 2
        do j = 1, 3
            do i = 1, 2
                one_rm = find_minimum(min_fun, minimum_triples_ref(:,i,j,k), xatol=1.0e-9_dp)
                call compare_minimum_refinement(rm(i,j,k), one_rm)
                if (rm(i,j,k)%success) then
                    call check(abs(rm(i,j,k)%x-find_min_x_ref(i,j,k)) < 1.0e-10_dp, &
                        'SciPy rank-three refined minimum value')
                end if
            end do
        end do
    end do
    call check(any(rm%status /= 0) .and. any(rm%status == 0), &
        'rank-three minimum failures do not affect successful elements')

    ! Optional arrays of rank 1, 2, and 3 align on the trailing axes.
    ! The rank-1 argument indexes k, a rank-2 argument indexes (j,k).
    upper_tensor = reshape([-1.0_dp, -1.0_dp, 4.0_dp, 4.0_dp], [2,1,2])
    calls = 0
    br = bracket_root(root_fun, root_left_ref, xr0=[-3.0_dp, -2.5_dp], &
        xmin=reshape([-10.0_dp, -1.0_dp, 1.0_dp], [3,1]), &
        xmax=upper_tensor)
    batch_calls = calls
    call check(batch_calls == sum(br%nfev), 'rank-three broadcast bracket count')
    do k = 1, 2
        do j = 1, 3
            do i = 1, 2
                one_br = bracket_root(root_fun, root_left_ref(i,j,k), &
                    xr0=merge(-3.0_dp, -2.5_dp, k == 1), &
                    xmin=merge(-10.0_dp, merge(-1.0_dp, 1.0_dp, j == 2), j == 1), &
                    xmax=merge(-1.0_dp, 4.0_dp, k == 1))
                call compare_root_bracket(br(i,j,k), one_br)
            end do
        end do
    end do

    tol_matrix = reshape([1.0e-12_dp, 1.0e-9_dp, 1.0e-10_dp, &
        1.0e-8_dp, 1.0e-11_dp, 1.0e-7_dp], [3,2])
    zero_bound = 0.0_dp
    calls = 0
    rr = find_root(root_fun, root_pairs_ref, xatol=tol_matrix, &
        xrtol=zero_bound, fatol=[0.0_dp], frtol=0.0_dp)
    batch_calls = calls
    call check(batch_calls == sum(rr%nfev), 'rank-three broadcast root tolerance count')
    do k = 1, 2
        do j = 1, 3
            do i = 1, 2
                one_rr = find_root(root_fun, root_pairs_ref(:,i,j,k), &
                    xatol=tol_matrix(j,k), xrtol=0.0_dp, fatol=0.0_dp, frtol=0.0_dp)
                call compare_root_refinement(rr(i,j,k), one_rr)
            end do
        end do
    end do

    calls = 0
    bm = bracket_minimum(min_fun, minimum_middle_ref, &
        xl0=minimum_middle_ref-0.3_dp, &
        xr0=reshape([-2.7_dp, 0.8_dp, 1.8_dp], [3,1]), &
        xmin=reshape([-10.0_dp, -5.0_dp], [2,1,1]), xmax=10.0_dp)
    batch_calls = calls
    call check(batch_calls == sum(bm%nfev), 'rank-three broadcast minimum bracket count')
    do k = 1, 2
        do j = 1, 3
            do i = 1, 2
                one_bm = bracket_minimum(min_fun, minimum_middle_ref(i,j,k), &
                    xl0=minimum_middle_ref(i,j,k)-0.3_dp, &
                    xr0=merge(-2.7_dp, merge(0.8_dp, 1.8_dp, j == 2), j == 1), &
                    xmin=merge(-10.0_dp, -5.0_dp, i == 1), xmax=10.0_dp)
                call compare_minimum_bracket(bm(i,j,k), one_bm)
            end do
        end do
    end do

    calls = 0
    rm = find_minimum(min_fun, minimum_triples_ref, &
        xatol=reshape([1.0e-9_dp, 1.0e-7_dp], [1,1,2]), &
        xrtol=reshape([0.0_dp, 0.0_dp, 0.0_dp], [3,1]), &
        fatol=[0.0_dp], frtol=0.0_dp)
    batch_calls = calls
    call check(batch_calls == sum(rm%nfev), 'rank-three broadcast minimum tolerance count')
    do k = 1, 2
        do j = 1, 3
            do i = 1, 2
                one_rm = find_minimum(min_fun, minimum_triples_ref(:,i,j,k), &
                    xatol=merge(1.0e-9_dp, 1.0e-7_dp, k == 1), &
                    xrtol=0.0_dp, fatol=0.0_dp, frtol=0.0_dp)
                call compare_minimum_refinement(rm(i,j,k), one_rm)
            end do
        end do
    end do

    ! No supplied endpoints must be passed through as absent (not zeros).
    calls = 0
    br = bracket_root(root_fun, root_left_ref)
    bm = bracket_minimum(min_fun, minimum_middle_ref)
    batch_calls = calls
    call check(batch_calls == sum(br%nfev) + sum(bm%nfev), &
        'rank-three default endpoint calls counted')
    do k = 1, 2
        do j = 1, 3
            do i = 1, 2
                one_br = bracket_root(root_fun, root_left_ref(i,j,k))
                one_bm = bracket_minimum(min_fun, minimum_middle_ref(i,j,k))
                call compare_root_bracket(br(i,j,k), one_br)
                call compare_minimum_bracket(bm(i,j,k), one_bm)
            end do
        end do
    end do

    roots_valid(1,:,:,:) = -4.0_dp
    roots_valid(2,:,:,:) = -1.0_dp
    callback_calls = 0
    rr = find_root(root_fun, roots_valid, callback=root_callback)
    call check(all(rr%status == find_root_status_callback), 'rank-three root callback statuses')
    call check(all(rr%nit == 1), 'rank-three root callback iteration count')
    call check(callback_calls == 24, 'rank-three root callbacks are independent')
    minimum_valid(1,:,:,:) = 0.0_dp
    minimum_valid(2,:,:,:) = 1.0_dp
    minimum_valid(3,:,:,:) = 3.0_dp
    callback_calls = 0
    rm = find_minimum(min_fun, minimum_valid, callback=minimum_callback)
    call check(all(rm%status == find_minimum_status_callback), 'rank-three minimum callback statuses')
    call check(all(rm%nit == 1), 'rank-three minimum callback iterations')
    call check(callback_calls == 24, 'rank-three minimum callbacks are independent')

    calls = 0
    br = bracket_root(root_fun, root_left_ref, xr0=[-3.0_dp, -2.5_dp, 1.0_dp])
    call check(all(br%status == root_status_invalid_input), 'root rank-1 size mismatch')
    bm = bracket_minimum(min_fun, minimum_middle_ref, &
        xmin=reshape([-10.0_dp, 0.0_dp, 10.0_dp, 20.0_dp], [2,2]))
    call check(all(bm%status == minimize_status_invalid_input), 'minimum rank-2 shape mismatch')
    rr = find_root(root_fun, root_pairs_ref, xatol=reshape([0.0_dp, 0.0_dp, 0.0_dp], [3,1,1]))
    call check(all(rr%status == root_status_invalid_input), 'root rank-3 shape mismatch')
    rm = find_minimum(min_fun, minimum_triples_ref, &
        xatol=reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [1,1,1,4]))
    call check(all(rm%status == minimize_status_invalid_input), 'minimum rank-4 tolerance mismatch')
    rr = find_root(root_fun, reshape([(real(i,dp), i=1,36)], [3,2,3,2]))
    call check(all(rr%status == root_status_invalid_input), 'root leading bracket dimension invalid')
    rm = find_minimum(min_fun, reshape([(real(i,dp), i=1,24)], [2,2,3,2]))
    call check(all(rm%status == minimize_status_invalid_input), &
        'minimum leading bracket dimension invalid')
    call check(calls == 0, 'rank-three invalid shapes never evaluate objectives')

    br = bracket_root(root_fun, root_left_ref(:,1:0,:))
    bm = bracket_minimum(min_fun, minimum_middle_ref(1:0,:,:))
    rr = find_root(root_fun, root_pairs_ref(:,:,:,1:0))
    rm = find_minimum(min_fun, minimum_triples_ref(:,1:0,:,:))
    call check(all(shape(br) == [2,0,2]) .and. all(shape(bm) == [0,3,2]), &
        'empty rank-three bracket grids preserve shapes')
    call check(all(shape(rr) == [2,3,0]) .and. all(shape(rm) == [0,3,2]), &
        'empty rank-three refinement grids preserve shapes')
    call check(calls == 0, 'empty rank-three grids never evaluate objectives')

    if (failures > 0) then
        print '(a,i0)', 'test_optimize_elementwise_rank3: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_elementwise_rank3: PASS'

contains

    function root_fun(x) result(fx)
        real(dp), intent(in) :: x !! scalar root-function argument
        real(dp) :: fx !! scalar function value
        calls = calls + 1
        fx = x*x-4.0_dp
    end function root_fun

    function min_fun(x) result(fx)
        real(dp), intent(in) :: x !! scalar minimization argument
        real(dp) :: fx !! scalar objective value
        calls = calls + 1
        fx = (x-1.5_dp)**2
    end function min_fun

    subroutine root_callback(state, stop)
        type(find_root_result), intent(in) :: state !! current scalar root state
        logical, intent(out) :: stop !! request early termination
        callback_calls = callback_calls + 1
        stop = state%nit >= 1
    end subroutine root_callback

    subroutine minimum_callback(state, stop)
        type(find_minimum_result), intent(in) :: state !! current scalar minimum state
        logical, intent(out) :: stop !! request early termination
        callback_calls = callback_calls + 1
        stop = state%nit >= 1
    end subroutine minimum_callback

    subroutine check(condition, message)
        logical, intent(in) :: condition !! expected condition
        character(len=*), intent(in) :: message !! regression diagnostic
        if (.not. condition) then
            failures = failures+1
            print '(a)', 'FAIL: '//message
        end if
    end subroutine check

    subroutine compare_root_bracket(a, b)
        type(bracket_root_result), intent(in) :: a !! rank-three value
        type(bracket_root_result), intent(in) :: b !! scalar reference
        call check(a%status == b%status .and. (a%success .eqv. b%success), 'root bracket status')
        call check(a%nit == b%nit .and. a%nfev == b%nfev, 'root bracket counts')
        if (a%success) call check(all(abs([a%xl-b%xl, a%xr-b%xr, a%fl-b%fl, a%fr-b%fr]) < 1.0e-12_dp), &
            'root bracket values')
    end subroutine compare_root_bracket

    subroutine compare_root_refinement(a, b)
        type(find_root_result), intent(in) :: a !! rank-three value
        type(find_root_result), intent(in) :: b !! scalar reference
        call check(a%status == b%status .and. (a%success .eqv. b%success), 'root refinement status')
        call check(a%nit == b%nit .and. a%nfev == b%nfev, 'root refinement counts')
        if (a%success) call check(abs(a%x-b%x) < 1.0e-12_dp, 'root refinement value')
    end subroutine compare_root_refinement

    subroutine compare_minimum_bracket(a, b)
        type(bracket_minimum_result), intent(in) :: a !! rank-three value
        type(bracket_minimum_result), intent(in) :: b !! scalar reference
        call check(a%status == b%status .and. (a%success .eqv. b%success), 'minimum bracket status')
        call check(a%nit == b%nit .and. a%nfev == b%nfev, 'minimum bracket counts')
        if (a%success) call check(all(abs([a%xl-b%xl, a%xm-b%xm, a%xr-b%xr]) < 1.0e-12_dp), &
            'minimum bracket positions')
    end subroutine compare_minimum_bracket

    subroutine compare_minimum_refinement(a, b)
        type(find_minimum_result), intent(in) :: a !! rank-three value
        type(find_minimum_result), intent(in) :: b !! scalar reference
        call check(a%status == b%status .and. (a%success .eqv. b%success), 'minimum refinement status')
        call check(a%nit == b%nit .and. a%nfev == b%nfev, 'minimum refinement counts')
        if (a%success) call check(abs(a%x-b%x) < 1.0e-12_dp, 'minimum refinement value')
    end subroutine compare_minimum_refinement
end program test_optimize_elementwise_rank3
