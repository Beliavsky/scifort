! SPDX-License-Identifier: MIT
! TEST-ONLY optimizer substitute, not a library source or numerical validation.
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
        function objective(x) result(f)
            import dp
            real(dp), intent(in) :: x(:)
            real(dp) :: f
        end function objective
    end interface

    public :: lbfgsb, powell
contains
    function lbfgsb(fun, x0, lower, upper, maxiter) result(res)
        procedure(objective) :: fun !! objective function
        real(dp), intent(in) :: x0(:) !! starting position
        real(dp), intent(in), optional :: lower(:) !! lower bounds
        real(dp), intent(in), optional :: upper(:) !! upper bounds
        integer, intent(in), optional :: maxiter !! local iteration cap
        type(multivariate_minimize_result) :: res !! local result
        call coordinate_search(fun, x0, lower, upper, res)
    end function lbfgsb

    function powell(fun, x0, lower, upper, maxiter) result(res)
        procedure(objective) :: fun !! objective function
        real(dp), intent(in) :: x0(:) !! starting position
        real(dp), intent(in), optional :: lower(:) !! lower bounds
        real(dp), intent(in), optional :: upper(:) !! upper bounds
        integer, intent(in), optional :: maxiter !! local iteration cap
        type(multivariate_minimize_result) :: res !! local result
        call coordinate_search(fun, x0, lower, upper, res)
    end function powell

    subroutine coordinate_search(fun, x0, lower, upper, res)
        procedure(objective) :: fun !! objective function
        real(dp), intent(in) :: x0(:) !! starting position
        real(dp), intent(in), optional :: lower(:) !! lower bounds
        real(dp), intent(in), optional :: upper(:) !! upper bounds
        type(multivariate_minimize_result), intent(out) :: res !! result
        real(dp) :: x(size(x0)), candidate(size(x0)), best, value, step
        integer :: i,j,k
        logical :: better

        x = x0
        best = fun(x)
        res%nfev = 1
        step = 0.5_dp
        do i = 1, 150
            better = .false.
            do j = 1, size(x)
                do k = -1,1,2
                    candidate = x
                    candidate(j) = x(j) + real(k,dp)*step
                    if (present(lower)) candidate = max(lower,candidate)
                    if (present(upper)) candidate = min(upper,candidate)
                    value = fun(candidate)
                    res%nfev = res%nfev + 1
                    if (value < best) then
                        x = candidate
                        best = value
                        better = .true.
                    end if
                end do
            end do
            if (.not. better) step = step*0.5_dp
            if (step < 1.0e-8_dp) exit
        end do
        res%x = x
        res%fun = best
        res%success = .true.
    end subroutine coordinate_search
end module scifort_optimize_minimize
