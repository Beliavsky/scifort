! SPDX-License-Identifier: MIT
! TEST-ONLY replacement for the large optimizer dependency in isolated builds.
! Never compile or link this file in the SciFort library/FPM build.
module scifort_optimize_minimize
    use scifort_kinds, only : dp
    implicit none
    private
    type, public :: multivariate_minimize_result
        real(dp), allocatable :: x(:)
        real(dp) :: fun = 0.0_dp
        integer :: nfev = 0
        logical :: success = .false.
    end type multivariate_minimize_result
    abstract interface
        function objective(x) result(fx)
            import dp
            real(dp), intent(in) :: x(:)
            real(dp) :: fx
        end function objective
    end interface
    public :: nelder_mead, powell, bfgs
contains
    function nelder_mead(fun, x0, maxiter) result(res)
        procedure(objective) :: fun
        real(dp), intent(in) :: x0(:)
        integer, intent(in), optional :: maxiter
        type(multivariate_minimize_result) :: res
        call search(fun,x0,res)
    end function nelder_mead

    function powell(fun, x0, maxiter) result(res)
        procedure(objective) :: fun
        real(dp), intent(in) :: x0(:)
        integer, intent(in), optional :: maxiter
        type(multivariate_minimize_result) :: res
        call search(fun,x0,res)
    end function powell

    function bfgs(fun, x0, maxiter) result(res)
        procedure(objective) :: fun
        real(dp), intent(in) :: x0(:)
        integer, intent(in), optional :: maxiter
        type(multivariate_minimize_result) :: res
        call search(fun,x0,res)
    end function bfgs

    subroutine search(fun, x0, res)
        procedure(objective) :: fun
        real(dp), intent(in) :: x0(:)
        type(multivariate_minimize_result), intent(out) :: res
        real(dp) :: x(size(x0)), trial(size(x0)), f, value, scale
        integer :: i, j, s
        logical :: improved
        x = x0
        f = fun(x)
        res%nfev = 1
        scale = 0.5_dp
        do i = 1, 500
            improved = .false.
            do j = 1, size(x)
                do s = -1, 1, 2
                    trial = x
                    trial(j) = x(j) + real(s, dp)*scale
                    value = fun(trial)
                    res%nfev = res%nfev + 1
                    if (value < f) then
                        x = trial
                        f = value
                        improved = .true.
                    end if
                end do
            end do
            if (.not.improved) scale = scale*0.5_dp
            if (scale < 1.0e-8_dp) exit
        end do
        res%x = x
        res%fun = f
        res%success = .true.
    end subroutine search
end module scifort_optimize_minimize
