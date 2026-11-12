! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Generalized simulated annealing with Tsallis visiting and acceptance laws,
! geometric cooling in the generalized-temperature coordinate, and reannealing.
! Inspired by the published Xiang et al. method and the documented semantics
! of scipy.optimize.dual_annealing. This implementation is independent;
! no SciPy or other third-party implementation code has been copied.
!
module scifort_optimize_dual_annealing
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: da_maxiter = 0
    integer, parameter, public :: da_maxfun = 1
    integer, parameter, public :: da_callback_stop = 2
    integer, parameter, public :: da_invalid = -1
    integer, parameter, public :: da_no_initial = -2
    integer(int64), parameter :: rng_modulus = 2147483647_int64
    integer(int64), parameter :: rng_multiplier = 48271_int64
    real(dp), parameter :: tail_limit = 1.0e8_dp
    real(dp), parameter :: pi = acos(-1.0_dp)

    type, public :: dual_annealing_result
        real(dp), allocatable :: x(:) !! best variables, shape (n)
        real(dp) :: fun = huge(1.0_dp) !! best finite objective
        real(dp), allocatable :: current_x(:) !! final accepted chain state
        real(dp) :: current_fun = huge(1.0_dp) !! final chain energy
        integer :: nit = 0 !! completed annealing temperature steps
        integer :: nfev = 0 !! counted objective calls, including local search
        integer :: nlocal = 0 !! local minimizer invocations
        integer :: naccepted = 0 !! accepted trial moves
        integer :: nrestarts = 0 !! temperature restarts
        integer :: nbest_updates = 0 !! strict new global minima after initialization
        integer :: status = da_invalid !! termination code
        logical :: success = .false. !! true on normal maxiter termination
        character(len=120) :: message = 'invalid input' !! termination reason
    end type dual_annealing_result

    abstract interface
        function da_objective(x) result(f)
            import dp
            real(dp), intent(in) :: x(:) !! point in the original domain, shape (n)
            real(dp) :: f !! objective value
        end function da_objective

        subroutine da_local_minimizer(fun, lower, upper, start, x, value, nfev, success)
            import dp, da_objective
            procedure(da_objective) :: fun !! objective to minimize
            real(dp), intent(in) :: lower(:) !! finite lower bounds, shape (n)
            real(dp), intent(in) :: upper(:) !! finite upper bounds, shape (n)
            real(dp), intent(in) :: start(:) !! starting point, shape (n)
            real(dp), intent(out) :: x(:) !! returned minimum, shape (n)
            real(dp), intent(out) :: value !! objective at returned minimum
            integer, intent(out) :: nfev !! objective calls made by local minimizer
            logical, intent(out) :: success !! whether local solver returned usable point
        end subroutine da_local_minimizer

        subroutine da_callback(x, value, context, stop)
            import dp
            real(dp), intent(in) :: x(:) !! new best point, shape (n)
            real(dp), intent(in) :: value !! new best objective
            integer, intent(in) :: context !! 0: annealing, 1: local minimization
            logical, intent(out) :: stop !! return true to stop search
        end subroutine da_callback
    end interface

    public :: dual_annealing_core, da_objective, da_local_minimizer, da_callback

contains

    function dual_annealing_core(fun, bounds, maxiter, initial_temp, restart_temp_ratio, &
            visit, accept, maxfun, seed, x0, minimizer, no_local_search, &
            local_interval, callback) result(res)
        procedure(da_objective) :: fun !! scalar objective callback
        real(dp), intent(in) :: bounds(:,:) !! finite bounds, shape (n,2)
        integer, intent(in), optional :: maxiter !! outer iteration cap, default 1000
        real(dp), intent(in), optional :: initial_temp !! initial temperature, default 5230
        real(dp), intent(in), optional :: restart_temp_ratio !! cooling restart ratio, default 2e-5
        real(dp), intent(in), optional :: visit !! Tsallis visiting parameter, default 2.62
        real(dp), intent(in), optional :: accept !! Tsallis acceptance parameter, default -5
        integer, intent(in), optional :: maxfun !! objective evaluation cap, default 10000000
        integer, intent(in), optional :: seed !! private reproducible RNG seed
        real(dp), intent(in), optional :: x0(:) !! starting point, default random within bounds
        procedure(da_local_minimizer), optional :: minimizer !! optional bounded local optimizer
        logical, intent(in), optional :: no_local_search !! disable local search, default false
        integer, intent(in), optional :: local_interval !! local search after this many steps, default 20
        procedure(da_callback), optional :: callback !! notified on strictly improved minima
        type(dual_annealing_result) :: res !! best point and diagnostic counters

        integer :: n, iter_limit, eval_limit, ls_interval, outer, age, j, k, tries
        integer :: local_nfev
        integer(int64) :: state, clock_count
        real(dp) :: temp0, ratio, qv, qa, temp, accept_temp, temperature_numer
        real(dp) :: energy, trial_energy, local_energy, delta, probability
        real(dp) :: log_factor, step, u, width
        real(dp), allocatable :: lower(:), upper(:), current(:), trial(:), local_x(:)
        logical :: use_local, valid, improved, local_ok, stop

        n = size(bounds, 1)
        if (n < 1 .or. size(bounds, 2) /= 2) then
            res%message = 'bounds must have shape (n,2), n >= 1'
            return
        end if
        if (.not. all(ieee_is_finite(bounds))) then
            res%message = 'bounds must be finite'
            return
        end if
        if (any(bounds(:,2) <= bounds(:,1))) then
            res%message = 'upper bounds must be strictly greater than lower bounds'
            return
        end if
        do j = 1, n
            if (bounds(j,1) < 0.0_dp .and. bounds(j,2) > 0.0_dp) then
                if (bounds(j,2) > huge(1.0_dp) + bounds(j,1)) then
                    res%message = 'bound widths overflow double precision'
                    return
                end if
            end if
        end do
        if (present(x0)) then
            if (size(x0) /= n) then
                res%message = 'x0 has wrong size'
                return
            end if
            if (.not. all(ieee_is_finite(x0))) then
                res%message = 'x0 must be finite'
                return
            end if
            if (any(x0 < bounds(:,1)) .or. any(x0 > bounds(:,2))) then
                res%message = 'x0 must satisfy bounds'
                return
            end if
        end if
        iter_limit = 1000
        if (present(maxiter)) iter_limit = maxiter
        eval_limit = 10000000
        if (present(maxfun)) eval_limit = maxfun
        ls_interval = 20
        if (present(local_interval)) ls_interval = local_interval
        if (iter_limit < 0 .or. eval_limit < 1 .or. ls_interval < 1) then
            res%message = 'invalid iteration, evaluation or local search limit'
            return
        end if
        temp0 = 5230.0_dp
        if (present(initial_temp)) temp0 = initial_temp
        ratio = 2.0e-5_dp
        if (present(restart_temp_ratio)) ratio = restart_temp_ratio
        qv = 2.62_dp
        if (present(visit)) qv = visit
        qa = -5.0_dp
        if (present(accept)) qa = accept
        if (.not. ieee_is_finite(temp0) .or. .not. ieee_is_finite(ratio) .or. &
            .not. ieee_is_finite(qv) .or. .not. ieee_is_finite(qa)) then
            res%message = 'temperature and distribution controls must be finite'
            return
        end if
        if (temp0 <= 0.01_dp .or. temp0 > 50000.0_dp .or. ratio <= 0.0_dp .or. &
            ratio >= 1.0_dp .or. qv <= 1.0_dp .or. qv >= 3.0_dp .or. &
            qa < -10000.0_dp .or. qa > -5.0_dp) then
            res%message = 'invalid temperature, restart, visit or acceptance parameter'
            return
        end if
        use_local = present(minimizer)
        if (present(no_local_search)) then
            if (no_local_search) use_local = .false.
            if (.not. no_local_search .and. .not. present(minimizer)) then
                res%message = 'local search requested without minimizer callback'
                return
            end if
        end if
        if (present(seed)) then
            state = 1_int64 + modulo(int(seed, int64), rng_modulus - 1_int64)
        else
            call system_clock(count=clock_count)
            state = 1_int64 + modulo(clock_count, rng_modulus - 1_int64)
        end if
        allocate(lower(n), upper(n), current(n), trial(n), local_x(n), res%x(n), res%current_x(n))
        lower = bounds(:,1)
        upper = bounds(:,2)
        res%x = lower
        res%current_x = lower

        if (present(x0)) then
            current = x0
        else
            call random_position(current)
        end if
        valid = .false.
        do tries = 1, 100
            call evaluate(current, energy, valid)
            if (valid) exit
            if (res%nfev >= eval_limit) exit
            call random_position(current)
        end do
        if (.not. valid) then
            res%status = da_no_initial
            res%message = 'could not find a finite initial objective'
            return
        end if
        res%x = current
        res%fun = energy
        res%current_x = current
        res%current_fun = energy
        res%status = da_maxiter
        res%message = 'maximum number of iterations reached'

        temperature_numer = exp((qv - 1.0_dp)*log(2.0_dp)) - 1.0_dp
        log_factor = visiting_log_constant(qv)
        age = 0
        do outer = 1, iter_limit
            if (res%nfev >= eval_limit) then
                res%status = da_maxfun
                res%message = 'maximum number of objective calls reached'
                res%current_x = current
                res%current_fun = energy
                return
            end if
            temp = temp0*temperature_numer / &
                (exp((qv - 1.0_dp)*log(real(age + 2, dp))) - 1.0_dp)
            if (temp < temp0*ratio) then
                valid = .false.
                do tries = 1, 100
                    if (res%nfev >= eval_limit) exit
                    call random_position(current)
                    call evaluate(current, energy, valid)
                    if (valid) exit
                end do
                res%nrestarts = res%nrestarts + 1
                age = 0
                if (.not. valid) then
                    if (res%nfev >= eval_limit) then
                        res%status = da_maxfun
                        res%message = 'maximum number of objective calls reached'
                        res%current_x = res%x
                        res%current_fun = res%fun
                        return
                    end if
                    current = res%x
                    energy = res%fun
                end if
                if (energy < res%fun) then
                    call best_update(current, energy, 0, stop, callback)
                    if (stop) return
                end if
                temp = temp0
            end if
            age = age + 1
            res%nit = outer
            accept_temp = temp/real(age,dp)
            improved = .false.
            do k = 1, 2*n
                if (res%nfev >= eval_limit) then
                    res%status = da_maxfun
                    res%message = 'maximum number of objective calls reached'
                    res%current_x = current
                    res%current_fun = energy
                    return
                end if
                trial = current
                if (k <= n) then
                    do j = 1, n
                        call visiting_step(temp, qv, log_factor, state, step)
                        width = upper(j) - lower(j)
                        trial(j) = lower(j) + modulo(modulo(current(j) - lower(j), width) + &
                            modulo(step, width), width)
                    end do
                else
                    j = k - n
                    call visiting_step(temp, qv, log_factor, state, step)
                    width = upper(j) - lower(j)
                    trial(j) = lower(j) + modulo(modulo(current(j) - lower(j), width) + &
                        modulo(step, width), width)
                end if
                call evaluate(trial, trial_energy, valid)
                if (.not. valid) cycle
                if (trial_energy < energy) then
                    current = trial
                    energy = trial_energy
                    res%naccepted = res%naccepted + 1
                    if (energy < res%fun) then
                        improved = .true.
                        call best_update(current, energy, 0, stop, callback)
                        if (stop) return
                    end if
                else
                    ! Avoid overflow if an extreme positive energy is compared
                    ! with an extreme negative one: this jump cannot be accepted.
                    if (trial_energy > 0.0_dp .and. energy < 0.0_dp) then
                        if (trial_energy > huge(1.0_dp) + energy) cycle
                    end if
                    delta = trial_energy - energy
                    ! Generalized Tsallis acceptance: take the logarithm to avoid overflow.
                    if (delta <= 0.0_dp) then
                        probability = 1.0_dp
                    else if (delta >= accept_temp/(1.0_dp - qa)) then
                        probability = 0.0_dp
                    else
                        u = (1.0_dp - qa)*(delta/accept_temp)
                        probability = exp(log(1.0_dp - u)/(1.0_dp - qa))
                    end if
                    if (uniform(state) < probability) then
                        current = trial
                        energy = trial_energy
                        res%naccepted = res%naccepted + 1
                    end if
                end if
            end do
            if (use_local .and. (improved .or. mod(outer, ls_interval) == 0)) then
                if (res%nfev >= eval_limit) cycle
                ! Search around the best known point, not an arbitrary uphill state.
                call minimizer(fun, lower, upper, res%x, local_x, local_energy, local_nfev, local_ok)
                res%nlocal = res%nlocal + 1
                res%nfev = res%nfev + max(0, local_nfev)
                if (local_nfev >= 0 .and. local_ok) then
                    valid = all(ieee_is_finite(local_x)) .and. ieee_is_finite(local_energy)
                    if (valid) valid = all(local_x >= lower) .and. all(local_x <= upper)
                    if (valid .and. local_energy < res%fun) then
                        current = local_x
                        energy = local_energy
                        call best_update(local_x, local_energy, 1, stop, callback)
                        if (stop) return
                    end if
                end if
            end if
            res%current_x = current
            res%current_fun = energy
        end do
        if (res%nfev >= eval_limit) then
            res%status = da_maxfun
            res%message = 'maximum number of objective calls reached'
        else
            res%success = .true.
        end if

    contains

        subroutine random_position(x)
            real(dp), intent(out) :: x(:) !! random feasible point
            integer :: m
            do m = 1, size(x)
                x(m) = lower(m) + (upper(m) - lower(m))*uniform(state)
            end do
        end subroutine random_position

        subroutine evaluate(x, val, finite)
            real(dp), intent(in) :: x(:) !! argument to objective
            real(dp), intent(out) :: val !! resulting function value
            logical, intent(out) :: finite !! whether objective result is finite
            val = fun(x)
            res%nfev = res%nfev + 1
            finite = ieee_is_finite(val)
        end subroutine evaluate

        subroutine best_update(x, val, context, request_stop, callback)
            procedure(da_callback), optional :: callback !! notified on improved minima
            real(dp), intent(in) :: x(:) !! improving vector
            real(dp), intent(in) :: val !! improved objective
            integer, intent(in) :: context !! 0 = annealing, 1 = local
            logical, intent(out) :: request_stop !! callback result
            res%x = x
            res%fun = val
            res%nbest_updates = res%nbest_updates + 1
            res%current_x = current
            res%current_fun = energy
            request_stop = .false.
            if (present(callback)) call callback(x, val, context, request_stop)
            if (request_stop) then
                res%status = da_callback_stop
                res%message = 'stopped by callback'
            end if
        end subroutine best_update

    end function dual_annealing_core

    function visiting_log_constant(q) result(log_constant)
        real(dp), intent(in) :: q !! Tsallis visiting parameter, 1 < q < 3
        real(dp) :: log_constant !! log of the visiting scale independent of T
        real(dp) :: factor5, d1, z, log_factor4, log_factor6, log_ratio

        factor5 = 1.0_dp/(q - 1.0_dp) - 0.5_dp
        d1 = 2.0_dp - factor5
        z = 1.0_dp - factor5
        log_factor4 = log(sqrt(pi)) + (4.0_dp - q)*log(q - 1.0_dp) - &
            ((2.0_dp - q)/(q - 1.0_dp))*log(2.0_dp) - log(3.0_dp - q)
        if (abs(z) < 1.0e-7_dp) then
            log_ratio = (pi*z)**2/6.0_dp
        else
            log_ratio = log(pi*abs(z)) - log(max(abs(sin(pi*z)),tiny(1.0_dp)))
        end if
        if (d1 > 0.0_dp) then
            log_factor6 = log_ratio - log_gamma(d1)
        else
            ! Reflection computes log(abs(Gamma(d1))) for negative d1.
            log_factor6 = log_ratio - log(pi) + &
                log(max(abs(sin(pi*d1)),tiny(1.0_dp))) + log_gamma(1.0_dp - d1)
        end if
        log_constant = -(q - 1.0_dp)/(3.0_dp - q)*(log_factor6 - log_factor4)
    end function visiting_log_constant

    subroutine visiting_step(temp, q, log_constant, state, step)
        real(dp), intent(in) :: temp !! current positive annealing temperature
        real(dp), intent(in) :: q !! visiting distribution shape
        real(dp), intent(in) :: log_constant !! precomputed scale constant
        integer(int64), intent(inout) :: state !! private RNG state
        real(dp), intent(out) :: step !! bounded generalized visiting displacement
        real(dp) :: a, b, exponent, log_step

        a = normal(state)
        b = normal(state)
        if (abs(a) <= tiny(1.0_dp)) then
            step = 0.0_dp
            return
        end if
        exponent = (q - 1.0_dp)/(3.0_dp - q)
        log_step = log(abs(a)) + log_constant + log(temp)/(3.0_dp - q) - &
            exponent*log(max(abs(b),tiny(1.0_dp)))
        if (log_step >= log(tail_limit)) then
            step = sign(tail_limit*uniform(state), a)
        else if (log_step <= log(tiny(1.0_dp))) then
            step = 0.0_dp
        else
            step = sign(exp(log_step),a)
        end if
    end subroutine visiting_step

    function uniform(state) result(u)
        integer(int64), intent(inout) :: state !! Park-Miller generator state
        real(dp) :: u !! uniform variate strictly between zero and one
        state = modulo(rng_multiplier*state, rng_modulus)
        u = real(state,dp)/real(rng_modulus,dp)
    end function uniform

    function normal(state) result(x)
        integer(int64), intent(inout) :: state !! private RNG state
        real(dp) :: x !! standard normal variate
        x = sqrt(-2.0_dp*log(uniform(state)))*cos(2.0_dp*pi*uniform(state))
    end function normal

end module scifort_optimize_dual_annealing
