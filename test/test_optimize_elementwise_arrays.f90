! SPDX-License-Identifier: MIT
! Regression for independent rank-one SciFort elementwise root/minimum solvers.
program test_optimize_elementwise_arrays
    use scifort_kinds, only : dp
    use scifort_optimize, only : bracket_root, bracket_root_result, find_root, find_root_result, &
        bracket_minimum, bracket_minimum_result, find_minimum, find_minimum_result, &
        root_status_invalid_input, minimize_status_invalid_input, &
        find_root_status_callback, find_minimum_status_callback
    use test_optimize_elementwise_arrays_reference
    implicit none

    type(bracket_root_result), allocatable :: br(:)
    type(bracket_minimum_result), allocatable :: bm(:)
    type(find_root_result), allocatable :: rr(:)
    type(find_minimum_result), allocatable :: rm(:)
    type(bracket_root_result) :: br_single
    type(bracket_minimum_result) :: bm_single
    type(find_root_result) :: rr_single
    type(find_minimum_result) :: rm_single
    real(dp) :: root_pairs(2,4), minimum_triples(3,4)
    real(dp) :: root_start(4), root_right(4), root_lower(4), root_upper(4)
    real(dp) :: minimum_start(4), minimum_left(4), minimum_right(4)
    real(dp) :: minimum_lower(4), minimum_upper(4)
    integer :: failures, calls, batch_calls, callback_calls, i, work

    failures = 0
    calls = 0
    callback_calls = 0

    root_start = [-4.0_dp, 0.0_dp, 3.0_dp, 2.0_dp]
    root_right = [-3.0_dp, 1.0_dp, 3.5_dp, 3.0_dp]
    root_lower = [-10.0_dp, 0.0_dp, 1.0_dp, 3.0_dp]
    root_upper = [-1.0_dp, 10.0_dp, 4.0_dp, 4.0_dp]

    br = bracket_root(root_fun, root_start, xr0=root_right, xmin=root_lower, &
        xmax=root_upper, factor=2.0_dp, maxiter=15)
    batch_calls = calls
    work = 0
    call check(size(br) == 4, 'bracket_root array length')
    do i = 1, 4
        br_single = bracket_root(root_fun, root_start(i), xr0=root_right(i), &
            xmin=root_lower(i), xmax=root_upper(i), factor=2.0_dp, maxiter=15)
        call compare_root_brackets(br(i), br_single, i)
        work = work + br(i)%nfev
    end do
    call check(batch_calls == work, 'root bracket independent evaluation count')
    call check(br(1)%success .and. br(2)%success .and. br(3)%success, 'root three brackets found')
    call check(all(br(1:3)%status == br_root_status_ref), 'SciPy root bracket statuses')
    call check(all(br(1:3)%nit == br_root_nit_ref), 'SciPy root bracket iterations')
    call check(all(br(1:3)%nfev == br_root_nfev_ref), 'SciPy root bracket evaluations')
    call check(all(abs(br(1:3)%xl-br_root_xl_ref) < 1.e-12_dp), 'SciPy root left brackets')
    call check(all(abs(br(1:3)%xr-br_root_xr_ref) < 1.e-12_dp), 'SciPy root right brackets')
    call check(br(4)%status == root_status_invalid_input, 'root invalid element isolated')

    calls = 0
    br = bracket_root(root_fun, root_start(1:3), xmin=root_lower(1:3), xmax=root_upper(1:3))
    batch_calls = calls
    work = 0
    do i = 1, 3
        br_single = bracket_root(root_fun, root_start(i), xmin=root_lower(i), xmax=root_upper(i))
        call compare_root_brackets(br(i), br_single, i)
        work = work + br(i)%nfev
    end do
    call check(batch_calls == work, 'root default right point matches scalar')

    root_pairs(:, 1) = [-4.0_dp, -1.0_dp]
    root_pairs(:, 2) = [1.0_dp, 4.0_dp]
    root_pairs(:, 3) = [3.0_dp, 4.0_dp]
    root_pairs(:, 4) = [4.0_dp, 1.0_dp]
    calls = 0
    rr = find_root(root_fun, root_pairs, xatol=1.0e-12_dp)
    batch_calls = calls
    work = 0
    do i = 1, 4
        rr_single = find_root(root_fun, root_pairs(:,i), xatol=1.0e-12_dp)
        call compare_root_results(rr(i), rr_single, i)
        work = work + rr(i)%nfev
    end do
    call check(batch_calls == work, 'root refinement independent evaluation count')
    call check(rr(1)%success .and. rr(2)%success .and. rr(4)%success, 'root refinements converge')
    call check(.not. rr(3)%success, 'invalid root bracket isolated')
    call check(all(rr%status == root_status_ref), 'SciPy root refinement statuses')
    call check(all(rr%nit == root_nit_ref), 'SciPy root refinement iterations')
    call check(all(rr%nfev == root_nfev_ref), 'SciPy root refinement evaluations')
    call check(all(abs(rr([1,2,4])%x-root_x_valid_ref) < 1.e-10_dp), 'SciPy root estimates')
    call check(abs(rr(1)%x+2.0_dp) < 1.0e-10_dp, 'negative root')
    call check(abs(rr(2)%x-2.0_dp) < 1.0e-10_dp, 'positive root')

    callback_calls = 0
    rr = find_root(root_fun, root_pairs(:,1:2), callback=root_callback)
    call check(all(rr%status == find_root_status_callback), 'root callback stops independently')
    call check(all(rr%nit == 1), 'root callback iteration preserved')
    call check(callback_calls == 4, 'root callback called once initially and after iteration')

    minimum_start = [-3.0_dp, 0.5_dp, 1.5_dp, 5.0_dp]
    minimum_left = minimum_start - 0.3_dp
    minimum_right = minimum_start + 0.3_dp
    minimum_lower = [-10.0_dp, -5.0_dp, -5.0_dp, 4.9_dp]
    minimum_upper = [10.0_dp, 10.0_dp, 10.0_dp, 6.0_dp]
    calls = 0
    bm = bracket_minimum(min_fun, minimum_start, xl0=minimum_left, xr0=minimum_right, &
        xmin=minimum_lower, xmax=minimum_upper)
    batch_calls = calls
    work = 0
    do i = 1, 4
        bm_single = bracket_minimum(min_fun, minimum_start(i), xl0=minimum_left(i), &
            xr0=minimum_right(i), xmin=minimum_lower(i), xmax=minimum_upper(i))
        call compare_minimum_brackets(bm(i), bm_single, i)
        work = work + bm(i)%nfev
    end do
    call check(batch_calls == work, 'minimum bracket independent evaluation count')
    call check(bm(1)%success .and. bm(2)%success .and. bm(3)%success, 'minimum bracket success')
    call check(bm(4)%status == minimize_status_invalid_input, 'minimum invalid element isolated')
    call check(all(bm(1:3)%status == br_min_status_ref), 'SciPy minimum bracket statuses')
    call check(all(bm(1:3)%nit == br_min_nit_ref), 'SciPy minimum bracket iterations')
    call check(all(bm(1:3)%nfev == br_min_nfev_ref), 'SciPy minimum bracket evaluations')
    call check(all(abs(bm(1:3)%xm-br_min_xm_ref) < 1.e-12_dp), 'SciPy minimum bracket positions')

    calls = 0
    bm = bracket_minimum(min_fun, minimum_start(1:3), &
        xmin=minimum_lower(1:3), xmax=minimum_upper(1:3))
    batch_calls = calls
    work = 0
    do i = 1, 3
        bm_single = bracket_minimum(min_fun, minimum_start(i), &
            xmin=minimum_lower(i), xmax=minimum_upper(i))
        call compare_minimum_brackets(bm(i), bm_single, i)
        work = work + bm(i)%nfev
    end do
    call check(batch_calls == work, 'minimum default triple matches scalar')

    minimum_triples(:,1) = [0.0_dp, 1.0_dp, 3.0_dp]
    minimum_triples(:,2) = [-3.0_dp, 0.0_dp, 4.0_dp]
    minimum_triples(:,3) = [3.0_dp, 4.0_dp, 5.0_dp]
    minimum_triples(:,4) = [3.0_dp, 1.0_dp, 0.0_dp]
    calls = 0
    rm = find_minimum(min_fun, minimum_triples, xatol=1.0e-9_dp)
    batch_calls = calls
    work = 0
    do i = 1, 4
        rm_single = find_minimum(min_fun, minimum_triples(:,i), xatol=1.0e-9_dp)
        call compare_minimum_results(rm(i), rm_single, i)
        work = work + rm(i)%nfev
    end do
    call check(batch_calls == work, 'minimum refinement independent evaluation count')
    call check(rm(1)%success .and. rm(2)%success .and. rm(4)%success, 'minimum refinements converge')
    call check(.not. rm(3)%success, 'invalid minimum bracket isolated')
    call check(all(rm%status == min_status_ref), 'SciPy minimum refinement statuses')
    call check(all(rm%nit == min_nit_ref), 'SciPy minimum refinement iterations')
    call check(all(rm%nfev == min_nfev_ref), 'SciPy minimum refinement evaluations')
    call check(all(abs(rm([1,2,4])%x-min_x_valid_ref) < 1.e-10_dp), 'SciPy minimum estimates')
    call check(abs(rm(1)%x-1.5_dp) < 1.0e-6_dp, 'minimum location')

    callback_calls = 0
    rm = find_minimum(min_fun, minimum_triples(:,1:2), callback=minimum_callback)
    call check(all(rm%status == find_minimum_status_callback), 'minimum callback stops independently')
    call check(all(rm%nit == 1), 'minimum callback iteration preserved')
    call check(callback_calls == 4, 'minimum callback called once initially and after iteration')

    calls = 0
    br = bracket_root(root_fun, root_start(1:2), xr0=[1.0_dp, 2.0_dp, 3.0_dp])
    call check(all(br%status == root_status_invalid_input), 'root mismatched optional size')
    call check(calls == 0, 'root mismatched size no objective evaluations')
    bm = bracket_minimum(min_fun, minimum_start(1:2), xmin=[-2.0_dp, -3.0_dp, -4.0_dp])
    call check(all(bm%status == minimize_status_invalid_input), 'minimum mismatched optional size')
    call check(calls == 0, 'minimum mismatched size no objective evaluations')
    rr = find_root(root_fun, reshape([0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], [3,2]))
    call check(all(rr%status == root_status_invalid_input), 'root wrong bracket leading dimension')
    rm = find_minimum(min_fun, reshape([0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp], [2,2]))
    call check(all(rm%status == minimize_status_invalid_input), 'minimum wrong bracket leading dimension')
    call check(calls == 0, 'invalid shapes do not call objective')

    ! A scalar bound or endpoint broadcasts to each independent problem.
    calls = 0
    br = bracket_root(root_fun, root_start(1:3), xmin=-10.0_dp, xmax=10.0_dp)
    batch_calls = calls
    work = 0
    do i = 1, 3
        br_single = bracket_root(root_fun, root_start(i), xmin=-10.0_dp, xmax=10.0_dp)
        call compare_root_brackets(br(i), br_single, i)
        work = work + br(i)%nfev
    end do
    call check(batch_calls == work, 'broadcast root scalar bounds')

    br = bracket_root(root_fun, [0.0_dp, 0.0_dp], xr0=1.0_dp, &
        xmin=0.0_dp, xmax=[5.0_dp, 10.0_dp])
    call check(all(br%success), 'broadcast root scalar right endpoint')

    calls = 0
    bm = bracket_minimum(min_fun, minimum_start(1:3), xmin=-10.0_dp, xmax=10.0_dp)
    batch_calls = calls
    work = 0
    do i = 1, 3
        bm_single = bracket_minimum(min_fun, minimum_start(i), xmin=-10.0_dp, xmax=10.0_dp)
        call compare_minimum_brackets(bm(i), bm_single, i)
        work = work + bm(i)%nfev
    end do
    call check(batch_calls == work, 'broadcast minimum scalar bounds')

    calls = 0
    rr = find_root(root_fun, root_pairs, xatol=[1.0e-12_dp, 1.0e-7_dp, &
        1.0e-12_dp, 1.0e-6_dp], xrtol=0.0_dp)
    batch_calls = calls
    work = 0
    do i = 1, 4
        select case(i)
        case(1,3)
            rr_single = find_root(root_fun, root_pairs(:,i), xatol=1.0e-12_dp, xrtol=0.0_dp)
        case(2)
            rr_single = find_root(root_fun, root_pairs(:,i), xatol=1.0e-7_dp, xrtol=0.0_dp)
        case default
            rr_single = find_root(root_fun, root_pairs(:,i), xatol=1.0e-6_dp, xrtol=0.0_dp)
        end select
        call compare_root_results(rr(i), rr_single, i)
        work = work + rr(i)%nfev
    end do
    call check(batch_calls == work, 'root vector tolerances and scalar broadcasting')

    calls = 0
    rm = find_minimum(min_fun, minimum_triples, xatol=[1.0e-9_dp, 1.0e-6_dp, &
        1.0e-9_dp, 1.0e-8_dp], fatol=0.0_dp)
    batch_calls = calls
    work = 0
    do i = 1, 4
        select case(i)
        case(1,3)
            rm_single = find_minimum(min_fun, minimum_triples(:,i), xatol=1.0e-9_dp, fatol=0.0_dp)
        case(2)
            rm_single = find_minimum(min_fun, minimum_triples(:,i), xatol=1.0e-6_dp, fatol=0.0_dp)
        case default
            rm_single = find_minimum(min_fun, minimum_triples(:,i), xatol=1.0e-8_dp, fatol=0.0_dp)
        end select
        call compare_minimum_results(rm(i), rm_single, i)
        work = work + rm(i)%nfev
    end do
    call check(batch_calls == work, 'minimum vector tolerances and scalar broadcasting')

    calls = 0
    rr = find_root(root_fun, root_pairs(:,1:2), &
        xatol=reshape([1.0_dp, 2.0_dp], [1,2]))
    call check(all(rr%status == root_status_invalid_input), 'root unsupported tolerance rank')
    rm = find_minimum(min_fun, minimum_triples(:,1:2), xatol=[1.0_dp, 2.0_dp, 3.0_dp])
    call check(all(rm%status == minimize_status_invalid_input), 'minimum tolerance wrong size')
    br = bracket_root(root_fun, root_start(1:2), xmax=reshape([1.0_dp, 2.0_dp], [2,1]))
    call check(all(br%status == root_status_invalid_input), 'root bound unsupported rank')
    bm = bracket_minimum(min_fun, minimum_start(1:2), xmax=reshape([1.0_dp, 2.0_dp], [2,1]))
    call check(all(bm%status == minimize_status_invalid_input), 'minimum bound unsupported rank')
    call check(calls == 0, 'invalid broadcast shapes do not evaluate objective')

    br = bracket_root(root_fun, root_start(1:0))
    bm = bracket_minimum(min_fun, minimum_start(1:0))
    rr = find_root(root_fun, root_pairs(:,1:0))
    rm = find_minimum(min_fun, minimum_triples(:,1:0))
    call check(size(br) == 0 .and. size(bm) == 0 .and. &
        size(rr) == 0 .and. size(rm) == 0, 'empty problems give empty results')

    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_elementwise_arrays: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_elementwise_arrays: PASS'

contains
    function root_fun(x) result(fx)
        real(dp), intent(in) :: x !! candidate abscissa
        real(dp) :: fx !! root objective value
        calls = calls + 1
        fx = x*x-4.0_dp
    end function root_fun

    function min_fun(x) result(fx)
        real(dp), intent(in) :: x !! candidate abscissa
        real(dp) :: fx !! minimization objective value
        calls = calls + 1
        fx = (x-1.5_dp)**2
    end function min_fun

    subroutine root_callback(state, stop)
        type(find_root_result), intent(in) :: state !! current root state
        logical, intent(out) :: stop !! stop each element after first iteration
        callback_calls = callback_calls + 1
        stop = state%nit >= 1
    end subroutine root_callback

    subroutine minimum_callback(state, stop)
        type(find_minimum_result), intent(in) :: state !! current minimum state
        logical, intent(out) :: stop !! stop each element after first iteration
        callback_calls = callback_calls + 1
        stop = state%nit >= 1
    end subroutine minimum_callback

    subroutine check(ok, name)
        logical, intent(in) :: ok !! truth of the condition being tested
        character(len=*), intent(in) :: name !! case identifier
        if (.not. ok) then
            print '(a)', 'FAIL: '//name
            failures = failures + 1
        end if
    end subroutine check

    subroutine compare_root_brackets(a, b, n)
        type(bracket_root_result), intent(in) :: a !! vector solution element
        type(bracket_root_result), intent(in) :: b !! scalar reference solution
        integer, intent(in) :: n !! element index
        call check(a%status == b%status .and. (a%success .eqv. b%success), 'root bracket status')
        call check(a%nfev == b%nfev .and. a%nit == b%nit, 'root bracket counts')
        call check(all(abs([a%xl-b%xl, a%xr-b%xr, a%fl-b%fl, a%fr-b%fr]) < 1.e-12_dp), &
            'root bracket data')
    end subroutine compare_root_brackets

    subroutine compare_root_results(a, b, n)
        type(find_root_result), intent(in) :: a !! vector solution element
        type(find_root_result), intent(in) :: b !! scalar reference solution
        integer, intent(in) :: n !! element index
        call check(a%status == b%status .and. (a%success .eqv. b%success), 'root result status')
        call check(a%nfev == b%nfev .and. a%nit == b%nit, 'root result counts')
        if (a%success) then
            call check(abs(a%x-b%x) < 1.e-12_dp .and. abs(a%f_x-b%f_x) < 1.e-12_dp, &
                'root result location')
        end if
    end subroutine compare_root_results

    subroutine compare_minimum_brackets(a, b, n)
        type(bracket_minimum_result), intent(in) :: a !! vector solution element
        type(bracket_minimum_result), intent(in) :: b !! scalar reference solution
        integer, intent(in) :: n !! element index
        call check(a%status == b%status .and. (a%success .eqv. b%success), 'minimum bracket status')
        call check(a%nfev == b%nfev .and. a%nit == b%nit, 'minimum bracket counts')
        call check(all(abs([a%xl-b%xl, a%xm-b%xm, a%xr-b%xr]) < 1.e-12_dp), 'minimum bracket data')
    end subroutine compare_minimum_brackets

    subroutine compare_minimum_results(a, b, n)
        type(find_minimum_result), intent(in) :: a !! vector solution element
        type(find_minimum_result), intent(in) :: b !! scalar reference solution
        integer, intent(in) :: n !! element index
        call check(a%status == b%status .and. (a%success .eqv. b%success), 'minimum result status')
        call check(a%nfev == b%nfev .and. a%nit == b%nit, 'minimum result counts')
        if (a%success) then
            call check(abs(a%x-b%x) < 1.e-12_dp .and. abs(a%f_x-b%f_x) < 1.e-12_dp, &
                'minimum result location')
        end if
    end subroutine compare_minimum_results
end program test_optimize_elementwise_arrays
