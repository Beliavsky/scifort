! SPDX-License-Identifier: MIT
program test_optimize_basinhopping_reference
    use scifort_kinds, only : dp
    use scifort_optimize_basinhopping, only : basinhopping_core, basinhopping_result, &
        bh_completed, bh_no_improvement, bh_callback_stop, bh_invalid, bh_local_failure
    implicit none
    type(basinhopping_result) :: a, b
    integer :: checks, callback_calls

    checks = 0
    callback_calls = 0
    a = basinhopping_core(doublewell, [2.4_dp], local_search, niter=110, seed=321, &
        temperature=0.3_dp, stepsize=1.5_dp)
    call must(a%success .and. a%status == bh_completed, 'double-well completes')
    call must(a%fun < -0.20_dp, 'double-well reaches deep global basin')
    call must(a%nlocal == 111 .and. a%naccepted + a%nrejected == 110, 'move accounting')
    call must(a%nfev > a%nlocal, 'objective evaluations counted')
    call must(abs(a%fun - doublewell(a%x)) <= 1.0e-12_dp, 'reported objective verified')
    call must(a%fun <= a%current_fun + 1.0e-12_dp, 'best is no worse than chain state')
    b = basinhopping_core(doublewell, [2.4_dp], local_search, niter=110, seed=321, &
        temperature=0.3_dp, stepsize=1.5_dp)
    call must(all(a%x == b%x) .and. a%fun == b%fun, 'seed reproducibility: solution')
    call must(a%naccepted == b%naccepted .and. a%nfev == b%nfev, 'seed reproducibility: counters')

    a = basinhopping_core(sphere, [0.0_dp], identity_local, niter=6, &
        temperature=0.0_dp, stepsize=1.0_dp, take_step=up_one)
    call must(a%status == bh_completed .and. a%nrejected == 6, 'T=0 rejects uphill moves')
    call must(a%naccepted == 0 .and. a%fun == 0.0_dp, 'T=0 global best')
    a = basinhopping_core(sphere, [0.0_dp], identity_local, niter=4, &
        temperature=0.0_dp, take_step=up_one, accept_test=force_yes)
    call must(a%naccepted == 4 .and. a%nrejected == 0, 'custom force accept')
    call must(a%current_x(1) == 4.0_dp .and. a%fun == 0.0_dp, 'best survives uphill walks')
    a = basinhopping_core(sphere, [0.0_dp], identity_local, niter=5, &
        temperature=1.0e100_dp, take_step=up_one, accept_test=force_no)
    call must(a%naccepted == 0 .and. a%nrejected == 5, 'custom force reject')

    a = basinhopping_core(sphere, [0.0_dp], identity_local, niter=10, &
        niter_success=2, stepsize=0.0_dp, seed=18)
    call must(a%status == bh_no_improvement .and. a%nit == 3, 'no improvement stop')
    a = basinhopping_core(sphere, [0.0_dp], identity_local, niter=5, &
        callback=stop_at_three, seed=1)
    call must(a%status == bh_callback_stop .and. a%nit == 3, 'callback stop at iteration 3')
    call must(callback_calls == 4, 'callback includes initial minimum')
    a = basinhopping_core(sphere, [0.0_dp], identity_local, niter=4, &
        callback=stop_at_zero)
    call must(a%status == bh_callback_stop .and. a%nit == 0, 'initial callback stop')
    a = basinhopping_core(sphere, [0.0_dp], identity_local, niter=1, &
        interval=2, stepsize=1.0_dp, stepwise_factor=0.5_dp, &
        take_step=up_one)
    call must(a%final_stepsize == 1.0_dp, 'custom step retains chosen step size')
    a = basinhopping_core(sphere, [0.0_dp], identity_local, niter=4, &
        interval=2, stepsize=1.0_dp, stepwise_factor=0.5_dp, &
        accept_test=force_yes)
    call must(a%final_stepsize == 1.0_dp, 'adaptive step both shrinks and grows')
    a = basinhopping_core(sphere, [0.0_dp], identity_local, niter=0, seed=1)
    call must(a%nlocal == 1 .and. a%nfev == 1 .and. a%nit == 0, 'niter=0 still minimizes')
    a = basinhopping_core(sphere, [0.0_dp], fail_positive, niter=10, &
        temperature=0.0_dp, niter_success=1, take_step=up_one)
    call must(a%status == bh_no_improvement .and. a%nfailed_local == 2, &
        'failed trials count toward no-improvement stop')
    a = basinhopping_core(sphere, [0.0_dp], failure_local, niter=10)
    call must(a%status == bh_local_failure .and. .not. a%success, 'initial local failure')
    a = basinhopping_core(sphere, [0.0_dp], identity_local, niter=-1)
    call must(a%status == bh_invalid .and. .not. allocated(a%x), 'invalid negative niter')
    a = basinhopping_core(sphere, [0.0_dp], identity_local, interval=0)
    call must(a%status == bh_invalid, 'invalid adaptation interval')
    a = basinhopping_core(sphere, [0.0_dp], identity_local, temperature=-0.1_dp)
    call must(a%status == bh_invalid, 'negative temperature rejected')
    a = basinhopping_core(sphere, [0.0_dp], identity_local, target_accept_rate=1.0_dp)
    call must(a%status == bh_invalid, 'invalid target acceptance rate')
    a = basinhopping_core(sphere, [0.0_dp], identity_local, stepwise_factor=1.0_dp)
    call must(a%status == bh_invalid, 'invalid adaptation factor')
    a = basinhopping_core(sphere, [0.0_dp], identity_local, stepsize=-1.0_dp)
    call must(a%status == bh_invalid, 'negative stepsize rejected')
    a = basinhopping_core(sphere, [0.0_dp], identity_local, niter_success=-1)
    call must(a%status == bh_invalid, 'negative niter_success rejected')
    a = basinhopping_core(sphere, [real(dp) ::], identity_local)
    call must(a%status == bh_invalid, 'empty x rejected')
    a = basinhopping_core(sphere, [0.0_dp], identity_local, niter=3, &
        take_step=up_one, accept_test=invalid_accept)
    call must(a%status == bh_invalid .and. .not.a%success, 'invalid accept decision rejected')
    write(*,'(a,1x,i0)') 'BASINHOPPING_NATIVE_CHECKS', checks
contains
    subroutine must(ok, label)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: label
        if (.not.ok) then
            write(*,'(a)') 'FAIL: '//label
            error stop 1
        end if
        checks = checks + 1
    end subroutine must

    function sphere(x) result(fx)
        real(dp), intent(in) :: x(:)
        real(dp) :: fx
        fx = sum(x*x)
    end function sphere

    function doublewell(x) result(fx)
        real(dp), intent(in) :: x(:)
        real(dp) :: fx
        fx = (x(1)*x(1) - 1.0_dp)**2 + 0.25_dp*x(1)
    end function doublewell

    subroutine identity_local(fun, x0, x, fx, nfev, success)
        interface
            function fun(x) result(fx)
                import dp
                real(dp), intent(in) :: x(:)
                real(dp) :: fx
            end function fun
        end interface
        real(dp), intent(in) :: x0(:)
        real(dp), intent(out) :: x(:)
        real(dp), intent(out) :: fx
        integer, intent(out) :: nfev
        logical, intent(out) :: success
        x = x0
        fx = fun(x)
        nfev = 1
        success = .true.
    end subroutine identity_local

    subroutine fail_positive(fun, x0, x, fx, nfev, success)
        interface
            function fun(x) result(fx)
                import dp
                real(dp), intent(in) :: x(:)
                real(dp) :: fx
            end function fun
        end interface
        real(dp), intent(in) :: x0(:)
        real(dp), intent(out) :: x(:)
        real(dp), intent(out) :: fx
        integer, intent(out) :: nfev
        logical, intent(out) :: success
        x = x0
        fx = fun(x)
        nfev = 1
        success = x0(1) <= 0.0_dp
    end subroutine fail_positive

    subroutine failure_local(fun, x0, x, fx, nfev, success)
        interface
            function fun(x) result(fx)
                import dp
                real(dp), intent(in) :: x(:)
                real(dp) :: fx
            end function fun
        end interface
        real(dp), intent(in) :: x0(:)
        real(dp), intent(out) :: x(:)
        real(dp), intent(out) :: fx
        integer, intent(out) :: nfev
        logical, intent(out) :: success
        x = x0
        fx = huge(1.0_dp)
        nfev = 0
        success = .false.
    end subroutine failure_local

    subroutine local_search(fun, x0, x, fx, nfev, success)
        interface
            function fun(x) result(fx)
                import dp
                real(dp), intent(in) :: x(:)
                real(dp) :: fx
            end function fun
        end interface
        real(dp), intent(in) :: x0(:)
        real(dp), intent(out) :: x(:)
        real(dp), intent(out) :: fx
        integer, intent(out) :: nfev
        logical, intent(out) :: success
        real(dp) :: delta, trial, y(size(x))
        integer :: i, j, s
        logical :: improved
        x = x0
        fx = fun(x)
        nfev = 1
        delta = 0.5_dp
        do i = 1, 600
            improved = .false.
            do j = 1, size(x)
                do s = -1, 1, 2
                    y = x
                    y(j) = x(j) + real(s, dp)*delta
                    trial = fun(y)
                    nfev = nfev + 1
                    if (trial < fx) then
                        x = y
                        fx = trial
                        improved = .true.
                    end if
                end do
            end do
            if (.not.improved) delta = delta*0.5_dp
            if (delta < 1.0e-8_dp) exit
        end do
        success = .true.
    end subroutine local_search

    subroutine up_one(x, stepsize)
        real(dp), intent(inout) :: x(:)
        real(dp), intent(in) :: stepsize
        x = x + 1.0_dp + 0.0_dp*stepsize
    end subroutine up_one

    subroutine force_yes(fnew, xnew, fold, xold, decision)
        real(dp), intent(in) :: fnew, fold, xnew(:), xold(:)
        integer, intent(out) :: decision
        decision = 1
    end subroutine force_yes

    subroutine force_no(fnew, xnew, fold, xold, decision)
        real(dp), intent(in) :: fnew, fold, xnew(:), xold(:)
        integer, intent(out) :: decision
        decision = -1
    end subroutine force_no

    subroutine invalid_accept(fnew, xnew, fold, xold, decision)
        real(dp), intent(in) :: fnew, fold, xnew(:), xold(:)
        integer, intent(out) :: decision
        decision = 2
    end subroutine invalid_accept

    subroutine stop_at_three(x, value, accepted, iteration, stop)
        real(dp), intent(in) :: x(:), value
        logical, intent(in) :: accepted
        integer, intent(in) :: iteration
        logical, intent(out) :: stop
        callback_calls = callback_calls + 1
        stop = iteration == 3
    end subroutine stop_at_three

    subroutine stop_at_zero(x, value, accepted, iteration, stop)
        real(dp), intent(in) :: x(:), value
        logical, intent(in) :: accepted
        integer, intent(in) :: iteration
        logical, intent(out) :: stop
        stop = iteration == 0
    end subroutine stop_at_zero
end program test_optimize_basinhopping_reference
