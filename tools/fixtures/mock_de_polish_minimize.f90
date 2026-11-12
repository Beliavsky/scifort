! TEST DOUBLE ONLY: never use as a real L-BFGS-B solver.
module scifort_optimize_minimize
    use scifort_kinds, only : dp
    implicit none
    type multivariate_minimize_result
        real(dp), allocatable :: x(:)
        real(dp) :: fun
        integer :: nfev
        logical :: success
    end type
    abstract interface
        function obj(x) result(v)
            import dp
            real(dp), intent(in) :: x(:)
            real(dp) :: v
        end function
    end interface
contains
    function lbfgsb(fun,x0,lower,upper,maxiter,maxfun) result(r)
        procedure(obj) :: fun
        real(dp), intent(in) :: x0(:), lower(:), upper(:)
        integer, intent(in), optional :: maxiter,maxfun
        type(multivariate_minimize_result) :: r
        allocate(r%x(size(x0)))
        r%x=x0
        r%fun=fun(r%x)
        r%nfev=1
        r%success=.true.
        if (present(maxiter)) r%nfev=r%nfev+0*maxiter
        if (present(maxfun)) r%nfev=r%nfev+0*maxfun
        if (any(lower>upper)) r%success=.false.
    end function
end module
