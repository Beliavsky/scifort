! SPDX-License-Identifier: MIT
! TEST ONLY. Not a numerical substitute for SciFort's actual local optimizers.
! Compile in a separate temporary directory to check SHGO adapter signatures.
module scifort_optimize_minimize
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private
    public :: nelder_mead, powell, slsqp, multivariate_minimize_result

    type :: multivariate_minimize_result
        real(dp), allocatable :: x(:)
        real(dp) :: fun = huge(1.0_dp)
        integer :: nfev = 0
    end type multivariate_minimize_result

    abstract interface
        function scalar_objective(x) result(f)
            import dp
            real(dp), intent(in) :: x(:) !! objective coordinate vector
            real(dp) :: f !! scalar function value
        end function scalar_objective
    end interface
contains
    function nelder_mead(fun, x0, lower, upper, maxiter, maxfev) result(res)
        procedure(scalar_objective) :: fun !! objective to evaluate
        real(dp), intent(in) :: x0(:) !! local search start
        real(dp), intent(in), optional :: lower(:) !! bound floor
        real(dp), intent(in), optional :: upper(:) !! bound ceiling
        integer, intent(in), optional :: maxiter !! unused iteration limit
        integer, intent(in), optional :: maxfev !! evaluation limit
        type(multivariate_minimize_result) :: res !! dummy initial evaluation
        res%x = x0
        res%fun = fun(x0)
        res%nfev = 1
        if (present(lower)) res%x = max(res%x,lower)
        if (present(upper)) res%x = min(res%x,upper)
        if (present(maxiter)) then
            if (maxiter < 0) error stop 'invalid test maxiter'
        end if
        if (present(maxfev)) then
            if (maxfev < 1) error stop 'invalid test maxfev'
        end if
    end function nelder_mead

    function powell(fun, x0, lower, upper, maxiter, maxfev) result(res)
        procedure(scalar_objective) :: fun !! objective to evaluate
        real(dp), intent(in) :: x0(:) !! local search start
        real(dp), intent(in), optional :: lower(:) !! bound floor
        real(dp), intent(in), optional :: upper(:) !! bound ceiling
        integer, intent(in), optional :: maxiter !! unused iteration limit
        integer, intent(in), optional :: maxfev !! evaluation limit
        type(multivariate_minimize_result) :: res !! dummy initial evaluation
        res = nelder_mead(fun,x0,lower,upper,maxiter,maxfev)
    end function powell
    function slsqp(fun, x0, lower, upper, constraints, constraint_jac, &
        meq, mineq, maxiter, ftol) result(res)
        procedure(scalar_objective) :: fun
        real(dp), intent(in) :: x0(:)
        real(dp), intent(in), optional :: lower(:)
        real(dp), intent(in), optional :: upper(:)
        integer, intent(in), optional :: meq
        integer, intent(in), optional :: mineq
        integer, intent(in), optional :: maxiter
        real(dp), intent(in), optional :: ftol
        interface
            subroutine test_constraints(x,c)
                import dp
                real(dp), intent(in) :: x(:)
                real(dp), intent(out) :: c(:)
            end subroutine test_constraints
            subroutine test_constraint_jac(x,a)
                import dp
                real(dp), intent(in) :: x(:) !! bounded coordinates
                real(dp), intent(out) :: a(:,:) !! one row per canonical constraint
            end subroutine test_constraint_jac
        end interface
        procedure(test_constraints), optional :: constraints
        procedure(test_constraint_jac), optional :: constraint_jac
        type(multivariate_minimize_result) :: res
        real(dp) :: cval(64), cnext(64), ajac(64,size(x0)), xtest(size(x0))
        integer :: m, j, k
        m = 0
        if (present(meq)) m = m+meq
        if (present(mineq)) m = m+mineq
        if (present(constraints) .and. m>0) call constraints(x0,cval(:m))
        if (present(constraint_jac) .and. m>0) then
            call constraint_jac(x0,ajac(:m,:))
            if (.not. all(ieee_is_finite(ajac(:m,:)))) error stop 'nonfinite SLSQP Jacobian'
            do j = 1, size(x0)
                xtest = x0
                xtest(j) = x0(j)+1.0e-6_dp
                call constraints(xtest,cnext(:m))
                do k = 1, m
                    ! Finite differences cannot resolve a finite bound near HUGE.
                    if (abs(cval(k)) > 1.0e12_dp) cycle
                    if (abs((cnext(k)-cval(k))/1.0e-6_dp-ajac(k,j)) > 1.0e-4_dp) &
                        error stop 'incorrect canonical SLSQP Jacobian'
                end do
            end do
        end if
        if (present(maxiter)) then
            if (maxiter<1) error stop 'invalid SLSQP stub maxiter'
        end if
        if (present(ftol)) then
            if (ftol<0.0_dp) error stop 'invalid SLSQP stub ftol'
        end if
        res = nelder_mead(fun,x0,lower,upper,maxiter)
    end function slsqp
end module scifort_optimize_minimize
