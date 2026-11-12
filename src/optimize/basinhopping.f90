! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Basin hopping: stochastic displacement, local minimization, and
! Metropolis acceptance. The algorithm follows Wales and Doye (1997)
! and the public behavior documented for scipy.optimize.basinhopping.
! This code is an independent implementation, not copied from SciPy.
!
module scifort_optimize_basinhopping
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: bh_completed = 0
    integer, parameter, public :: bh_no_improvement = 1
    integer, parameter, public :: bh_callback_stop = 2
    integer, parameter, public :: bh_invalid = -1
    integer, parameter, public :: bh_local_failure = -2
    integer(int64), parameter :: rng_modulus = 2147483647_int64
    integer(int64), parameter :: rng_multiplier = 48271_int64

    type, public :: basinhopping_result
        real(dp), allocatable :: x(:) !! best minimum found; shape (n)
        real(dp) :: fun = huge(1.0_dp) !! objective at the best minimum
        real(dp), allocatable :: current_x(:) !! final accepted chain state; shape (n)
        real(dp) :: current_fun = huge(1.0_dp) !! value at final chain state
        real(dp) :: final_stepsize = 0.0_dp !! adapted random-displacement amplitude
        integer :: nit = 0 !! number of attempted basin-hopping steps
        integer :: nfev = 0 !! objective calls, as reported by local minimizer
        integer :: nlocal = 0 !! local-minimizer invocations including initial minimization
        integer :: naccepted = 0 !! accepted moves (excluding initial minimization)
        integer :: nrejected = 0 !! rejected moves
        integer :: nfailed_local = 0 !! failed trial minimizations
        integer :: nbest_updates = 0 !! strict improvements after initial minimization
        integer :: status = bh_invalid !! termination code
        logical :: success = .false. !! true if at least one valid local minimum was obtained
        character(len=120) :: message = 'invalid input' !! termination description
    end type basinhopping_result

    abstract interface
        function bh_objective(x) result(value)
            import dp
            real(dp), intent(in) :: x(:) !! objective location, shape (n)
            real(dp) :: value !! scalar objective value
        end function bh_objective

        subroutine bh_local_minimizer(fun, x0, x, value, nfev, success)
            import dp, bh_objective
            procedure(bh_objective) :: fun !! objective to minimize
            real(dp), intent(in) :: x0(:) !! local search start, shape (n)
            real(dp), intent(out) :: x(:) !! local minimum, shape (n)
            real(dp), intent(out) :: value !! objective at x
            integer, intent(out) :: nfev !! objective calls spent by the local solver
            logical, intent(out) :: success !! local convergence flag
        end subroutine bh_local_minimizer

        subroutine bh_take_step(x, stepsize)
            import dp
            real(dp), intent(inout) :: x(:) !! input state / proposed step, shape (n)
            real(dp), intent(in) :: stepsize !! current step amplitude
        end subroutine bh_take_step

        subroutine bh_accept_test(fnew, xnew, fold, xold, decision)
            import dp
            real(dp), intent(in) :: fnew !! trial local objective
            real(dp), intent(in) :: xnew(:) !! trial local minimum
            real(dp), intent(in) :: fold !! current local objective
            real(dp), intent(in) :: xold(:) !! current accepted local minimum
            integer, intent(out) :: decision !! -1 reject, 0 Metropolis, +1 force accept
        end subroutine bh_accept_test

        subroutine bh_callback(x, value, accepted, iteration, stop)
            import dp
            real(dp), intent(in) :: x(:) !! latest local minimum (initial when iteration = 0)
            real(dp), intent(in) :: value !! objective at that point
            logical, intent(in) :: accepted !! whether the latest trial was accepted
            integer, intent(in) :: iteration !! attempt index, 0 for the initial minimum
            logical, intent(out) :: stop !! stop the optimization if true
        end subroutine bh_callback
    end interface

    public :: basinhopping_core, bh_objective, bh_local_minimizer
    public :: bh_take_step, bh_accept_test, bh_callback

contains

    function basinhopping_core(fun, x0, minimizer, niter, temperature, stepsize, &
            interval, niter_success, seed, target_accept_rate, stepwise_factor, &
            take_step, accept_test, callback) result(res)
        procedure(bh_objective) :: fun !! objective function to minimize
        real(dp), intent(in) :: x0(:) !! initial state, shape (n), must be finite
        procedure(bh_local_minimizer) :: minimizer !! caller-provided local solver
        integer, intent(in), optional :: niter !! number of attempted moves, default 100
        real(dp), intent(in), optional :: temperature !! Metropolis temperature, default 1
        real(dp), intent(in), optional :: stepsize !! initial displacement radius, default 0.5
        integer, intent(in), optional :: interval !! step adaptation interval, default 50
        integer, intent(in), optional :: niter_success !! early stop after this many stale moves
        integer, intent(in), optional :: seed !! seed for local independent RNG
        real(dp), intent(in), optional :: target_accept_rate !! in (0,1), default 0.5
        real(dp), intent(in), optional :: stepwise_factor !! adaptation multiplier in (0,1), default 0.9
        procedure(bh_take_step), optional :: take_step !! custom trial step (must manage its own RNG)
        procedure(bh_accept_test), optional :: accept_test !! -1 reject, 0 default, 1 force accept
        procedure(bh_callback), optional :: callback !! initial and per-step callback
        type(basinhopping_result) :: res !! best solution and diagnostics

        integer :: iterations, adapt_interval, no_improve_limit, i, j, neval, decision, stale
        integer(int64) :: state, clock_count
        real(dp) :: temp, step, target, factor, trial_value, accepted_before
        real(dp), allocatable :: xtrial(:), xmin(:)
        logical :: local_ok, accepted, stop, improves

        iterations = 100
        if (present(niter)) iterations = niter
        adapt_interval = 50
        if (present(interval)) adapt_interval = interval
        temp = 1.0_dp
        if (present(temperature)) temp = temperature
        step = 0.5_dp
        if (present(stepsize)) step = stepsize
        target = 0.5_dp
        if (present(target_accept_rate)) target = target_accept_rate
        factor = 0.9_dp
        if (present(stepwise_factor)) factor = stepwise_factor
        no_improve_limit = huge(0)
        if (present(niter_success)) no_improve_limit = niter_success

        if (size(x0) < 1) then
            res%message = 'x0 must contain at least one coordinate'
            return
        end if
        if (.not. all(ieee_is_finite(x0))) then
            res%message = 'x0 must be finite'
            return
        end if
        if (iterations < 0 .or. adapt_interval < 1 .or. no_improve_limit < 0) then
            res%message = 'invalid iteration or adaptation count'
            return
        end if
        if (.not. ieee_is_finite(temp) .or. .not. ieee_is_finite(step) .or. &
            .not. ieee_is_finite(target) .or. .not. ieee_is_finite(factor)) then
            res%message = 'algorithm parameters must be finite'
            return
        end if
        if (temp < 0.0_dp .or. step < 0.0_dp .or. &
            target <= 0.0_dp .or. target >= 1.0_dp .or. &
            factor <= 0.0_dp .or. factor >= 1.0_dp) then
            res%message = 'invalid temperature, step, or adaptation parameter'
            return
        end if

        if (present(seed)) then
            state = 1_int64 + modulo(int(seed, int64), rng_modulus - 1_int64)
        else
            call system_clock(count=clock_count)
            state = 1_int64 + modulo(clock_count, rng_modulus - 1_int64)
        end if

        allocate(res%x(size(x0)), res%current_x(size(x0)))
        allocate(xtrial(size(x0)), xmin(size(x0)))
        call minimizer(fun, x0, xmin, trial_value, neval, local_ok)
        res%nlocal = 1
        res%nfev = max(0, neval)
        if (neval < 0 .or. .not. local_ok .or. .not. valid_minimum(xmin, trial_value)) then
            res%message = 'initial local minimization failed'
            res%status = bh_local_failure
            res%x = x0
            res%current_x = x0
            return
        end if
        res%x = xmin
        res%fun = trial_value
        res%current_x = xmin
        res%current_fun = trial_value
        res%final_stepsize = step
        res%status = bh_completed
        res%success = .true.
        res%message = 'completed requested basin-hopping iterations'

        if (present(callback)) then
            call callback(xmin, trial_value, .true., 0, stop)
            if (stop) then
                res%status = bh_callback_stop
                res%message = 'stopped by callback at initial minimum'
                return
            end if
        end if

        stale = 0
        do i = 1, iterations
            res%nit = i
            if (mod(i, adapt_interval) == 0) then
                accepted_before = real(res%naccepted, dp)/real(i, dp)
                if (accepted_before > target) then
                    step = min(huge(1.0_dp)*factor, step)/factor
                else
                    step = step*factor
                end if
            end if
            res%final_stepsize = step
            xtrial = res%current_x
            if (present(take_step)) then
                call take_step(xtrial, step)
            else
                do j = 1, size(xtrial)
                    xtrial(j) = xtrial(j) + step*(2.0_dp*uniform(state) - 1.0_dp)
                end do
            end if
            if (.not. all(ieee_is_finite(xtrial))) then
                res%nrejected = res%nrejected + 1
                res%nfailed_local = res%nfailed_local + 1
                stale = stale + 1
                if (stale > no_improve_limit) then
                    res%status = bh_no_improvement
                    res%message = 'no improved minimum within requested iteration limit'
                    return
                end if
                cycle
            end if

            call minimizer(fun, xtrial, xmin, trial_value, neval, local_ok)
            res%nlocal = res%nlocal + 1
            if (neval > 0) res%nfev = res%nfev + neval
            if (neval < 0 .or. .not. local_ok .or. .not. valid_minimum(xmin, trial_value)) then
                res%nrejected = res%nrejected + 1
                res%nfailed_local = res%nfailed_local + 1
                stale = stale + 1
                if (stale > no_improve_limit) then
                    res%status = bh_no_improvement
                    res%message = 'no improved minimum within requested iteration limit'
                    return
                end if
                cycle
            end if

            decision = 0
            if (present(accept_test)) then
                call accept_test(trial_value, xmin, res%current_fun, res%current_x, decision)
                if (decision < -1 .or. decision > 1) then
                    res%status = bh_invalid
                    res%success = .false.
                    res%message = 'accept_test returned an invalid decision'
                    return
                end if
            end if
            accepted = .false.
            if (decision == 1) then
                accepted = .true.
            else if (decision == 0) then
                if (trial_value <= res%current_fun) then
                    accepted = .true.
                else if (temp > 0.0_dp) then
                    accepted = log(uniform(state)) < -(trial_value - res%current_fun)/temp
                end if
            end if

            if (accepted) then
                res%naccepted = res%naccepted + 1
                res%current_x = xmin
                res%current_fun = trial_value
            else
                res%nrejected = res%nrejected + 1
            end if
            improves = accepted .and. trial_value < res%fun
            if (improves) then
                res%x = xmin
                res%fun = trial_value
                res%nbest_updates = res%nbest_updates + 1
                stale = 0
            else
                stale = stale + 1
            end if
            if (present(callback)) then
                call callback(xmin, trial_value, accepted, i, stop)
                if (stop) then
                    res%status = bh_callback_stop
                    res%message = 'stopped by callback'
                    return
                end if
            end if
            if (stale > no_improve_limit) then
                res%status = bh_no_improvement
                res%message = 'no improved minimum within requested iteration limit'
                return
            end if
        end do
    end function basinhopping_core

    pure logical function valid_minimum(x, value)
        real(dp), intent(in) :: x(:) !! candidate minimum
        real(dp), intent(in) :: value !! candidate objective
        valid_minimum = all(ieee_is_finite(x)) .and. ieee_is_finite(value)
    end function valid_minimum

    function uniform(state) result(value)
        integer(int64), intent(inout) :: state !! Park-Miller state, in [1, modulus - 1]
        real(dp) :: value !! pseudorandom number in (0,1)
        state = modulo(rng_multiplier*state, rng_modulus)
        value = real(state, dp)/real(rng_modulus, dp)
    end function uniform

end module scifort_optimize_basinhopping
