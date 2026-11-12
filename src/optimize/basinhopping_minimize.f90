! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Connect the independent basin-hopping engine to existing SciFort
! Nelder-Mead, Powell and BFGS local minimizers.
!
module scifort_optimize_basinhopping_minimize
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_optimize_basinhopping, only : basinhopping_core, basinhopping_result, &
        bh_objective, bh_local_minimizer, bh_take_step, bh_accept_test, bh_callback
    use scifort_optimize_minimize, only : nelder_mead, powell, bfgs, multivariate_minimize_result
    implicit none
    private
    public :: basinhopping

contains

    function basinhopping(fun, x0, niter, temperature, stepsize, interval, &
            niter_success, seed, target_accept_rate, stepwise_factor, method, &
            local_maxiter, minimizer, take_step, accept_test, callback) result(res)
        procedure(bh_objective) :: fun !! scalar objective to minimize
        real(dp), intent(in) :: x0(:) !! initial variables, size n
        integer, intent(in), optional :: niter !! outer moves (default 100)
        real(dp), intent(in), optional :: temperature !! Metropolis temperature (default 1)
        real(dp), intent(in), optional :: stepsize !! maximum coordinate step (default 0.5)
        integer, intent(in), optional :: interval !! adaptation interval (default 50)
        integer, intent(in), optional :: niter_success !! optional no-improvement stopping limit
        integer, intent(in), optional :: seed !! reproducible seed, independent from global RNG
        real(dp), intent(in), optional :: target_accept_rate !! desired acceptance fraction (default 0.5)
        real(dp), intent(in), optional :: stepwise_factor !! adaptation fraction (default 0.9)
        character(len=*), intent(in), optional :: method !! 'nelder-mead', 'powell', 'bfgs'
        integer, intent(in), optional :: local_maxiter !! optional maximum local iterations
        procedure(bh_local_minimizer), optional :: minimizer !! custom local minimizer
        procedure(bh_take_step), optional :: take_step !! custom perturbation
        procedure(bh_accept_test), optional :: accept_test !! optional acceptance override
        procedure(bh_callback), optional :: callback !! early stopping / inspection
        type(basinhopping_result) :: res !! best minimum and diagnostic counts

        character(len=32) :: local_method

        local_method = 'nelder-mead'
        if (present(method)) local_method = method
        if (.not. present(minimizer)) then
            select case (trim(local_method))
            case ('nelder-mead', 'powell', 'bfgs')
                continue
            case default
                res%message = 'unsupported local method'
                return
            end select
        end if
        if (present(local_maxiter)) then
            if (local_maxiter < 0) then
                res%message = 'local_maxiter must be nonnegative'
                return
            end if
        end if

        if (present(minimizer)) then
            res = basinhopping_core(fun, x0, minimizer, niter=niter, temperature=temperature, &
                stepsize=stepsize, interval=interval, niter_success=niter_success, seed=seed, &
                target_accept_rate=target_accept_rate, stepwise_factor=stepwise_factor, &
                take_step=take_step, accept_test=accept_test, callback=callback)
        else
            res = basinhopping_core(fun, x0, local_min, niter=niter, temperature=temperature, &
                stepsize=stepsize, interval=interval, niter_success=niter_success, seed=seed, &
                target_accept_rate=target_accept_rate, stepwise_factor=stepwise_factor, &
                take_step=take_step, accept_test=accept_test, callback=callback)
        end if

    contains

        subroutine local_min(objective, start, x, value, nfev, success)
            procedure(bh_objective) :: objective !! local objective
            real(dp), intent(in) :: start(:) !! local search starting point
            real(dp), intent(out) :: x(:) !! final local point
            real(dp), intent(out) :: value !! objective at local minimum
            integer, intent(out) :: nfev !! local objective evaluation count
            logical, intent(out) :: success !! local solver converged
            type(multivariate_minimize_result) :: local

            select case (trim(local_method))
            case ('nelder-mead')
                local = nelder_mead(objective, start, maxiter=local_maxiter)
            case ('powell')
                local = powell(objective, start, maxiter=local_maxiter)
            case ('bfgs')
                local = bfgs(objective, start, maxiter=local_maxiter)
            end select
            nfev = max(0, local%nfev)
            success = local%success
            x = start
            value = huge(1.0_dp)
            if (.not. allocated(local%x)) then
                success = .false.
                return
            end if
            if (size(local%x) /= size(start)) then
                success = .false.
                return
            end if
            if (.not. ieee_is_finite(local%fun)) then
                success = .false.
                return
            end if
            x = local%x
            value = local%fun
        end subroutine local_min

    end function basinhopping

end module scifort_optimize_basinhopping_minimize
