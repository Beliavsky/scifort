! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Bounded local-search adapter for the stand-alone generalized annealing core.
! This module uses the already-declared SciFort optimizer dependencies.
!
module scifort_optimize_dual_annealing_minimize
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_optimize_dual_annealing, only : dual_annealing_core, &
        dual_annealing_result, da_objective, da_local_minimizer, da_callback
    use scifort_optimize_minimize, only : lbfgsb, powell, multivariate_minimize_result
    implicit none
    private
    public :: dual_annealing

contains

    function dual_annealing(fun, bounds, maxiter, initial_temp, restart_temp_ratio, &
            visit, accept, maxfun, seed, x0, no_local_search, method, &
            local_maxiter, local_interval, minimizer, callback) result(res)
        procedure(da_objective) :: fun !! scalar objective callback
        real(dp), intent(in) :: bounds(:,:) !! bounds, shape (n,2)
        integer, intent(in), optional :: maxiter !! outer iteration budget, default 1000
        real(dp), intent(in), optional :: initial_temp !! initial temperature, default 5230
        real(dp), intent(in), optional :: restart_temp_ratio !! temperature restart ratio
        real(dp), intent(in), optional :: visit !! Tsallis visiting parameter
        real(dp), intent(in), optional :: accept !! acceptance parameter
        integer, intent(in), optional :: maxfun !! soft objective evaluation budget
        integer, intent(in), optional :: seed !! reproducible seed
        real(dp), intent(in), optional :: x0(:) !! optional feasible initial location
        logical, intent(in), optional :: no_local_search !! skip polishing, default false
        character(len=*), intent(in), optional :: method !! 'l-bfgs-b' (default), 'powell'
        integer, intent(in), optional :: local_maxiter !! local iteration cap, default 200
        integer, intent(in), optional :: local_interval !! local search frequency, default 20
        procedure(da_local_minimizer), optional :: minimizer !! custom bounded local minimizer
        procedure(da_callback), optional :: callback !! triggered for improved best points
        type(dual_annealing_result) :: res !! best point and diagnostic counts

        character(len=24) :: local_method
        integer :: local_limit
        logical :: do_local

        local_method = 'l-bfgs-b'
        if (present(method)) local_method = trim(method)
        local_limit = 200
        if (present(local_maxiter)) local_limit = local_maxiter
        do_local = .true.
        if (present(no_local_search)) do_local = .not. no_local_search

        if (local_limit < 0) then
            res%message = 'local_maxiter must be nonnegative'
            return
        end if
        if (do_local .and. .not. present(minimizer)) then
            if (trim(local_method) /= 'l-bfgs-b' .and. trim(local_method) /= 'powell') then
                res%message = 'unsupported local search method'
                return
            end if
        end if

        if (.not. do_local) then
            res = dual_annealing_core(fun, bounds, maxiter=maxiter, initial_temp=initial_temp, &
                restart_temp_ratio=restart_temp_ratio, visit=visit, accept=accept, &
                maxfun=maxfun, seed=seed, x0=x0, no_local_search=.true., &
                local_interval=local_interval, callback=callback)
        else if (present(minimizer)) then
            res = dual_annealing_core(fun, bounds, maxiter=maxiter, initial_temp=initial_temp, &
                restart_temp_ratio=restart_temp_ratio, visit=visit, accept=accept, &
                maxfun=maxfun, seed=seed, x0=x0, minimizer=minimizer, &
                local_interval=local_interval, callback=callback)
        else
            res = dual_annealing_core(fun, bounds, maxiter=maxiter, initial_temp=initial_temp, &
                restart_temp_ratio=restart_temp_ratio, visit=visit, accept=accept, &
                maxfun=maxfun, seed=seed, x0=x0, minimizer=local_search, &
                local_interval=local_interval, callback=callback)
        end if

    contains

        subroutine local_search(objective, lower, upper, start, x, value, nfev, success)
            procedure(da_objective) :: objective !! scalar objective callback
            real(dp), intent(in) :: lower(:) !! lower bounds, shape (n)
            real(dp), intent(in) :: upper(:) !! upper bounds, shape (n)
            real(dp), intent(in) :: start(:) !! initial feasible point, shape (n)
            real(dp), intent(out) :: x(:) !! final local point, shape (n)
            real(dp), intent(out) :: value !! objective at final local point
            integer, intent(out) :: nfev !! objective evaluation count
            logical, intent(out) :: success !! valid converged local result
            type(multivariate_minimize_result) :: local

            select case (trim(local_method))
            case ('l-bfgs-b')
                local = lbfgsb(objective, start, lower=lower, upper=upper, maxiter=local_limit)
            case ('powell')
                local = powell(objective, start, lower=lower, upper=upper, maxiter=local_limit)
            end select
            x = start
            value = huge(1.0_dp)
            nfev = max(0, local%nfev)
            success = .false.
            if (.not. allocated(local%x)) return
            if (size(local%x) /= size(start)) return
            if (.not. ieee_is_finite(local%fun)) return
            if (.not. all(ieee_is_finite(local%x))) return
            if (any(local%x < lower) .or. any(local%x > upper)) return
            x = local%x
            value = local%fun
            success = local%success
        end subroutine local_search

    end function dual_annealing

end module scifort_optimize_dual_annealing_minimize
