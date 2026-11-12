! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Dense nonnegative and bounded-variable linear least squares.
! Independent active-set implementation based on the Lawson-Hanson/BVLS
! strategy, using SciFort's Jacobi-SVD least-squares kernel for free columns.
! No SLSQP source code is reproduced in this file.
module scifort_optimize_linear_least_squares
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_value, ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_svd, only : svd_lstsq
    implicit none
    private

    integer, parameter, public :: linear_lsq_success = 0
    integer, parameter, public :: linear_lsq_invalid = -1
    integer, parameter, public :: linear_lsq_maxiter = -2
    integer, parameter, public :: linear_lsq_numerical = -3

    type, public :: linear_lsq_result
        real(dp), allocatable :: x(:) !! optimal coefficients
        real(dp), allocatable :: fun(:) !! residual A*x - b
        real(dp), allocatable :: grad(:) !! A^T * fun
        integer, allocatable :: active_mask(:) !! -1 lower, 0 free, +1 upper
        real(dp) :: cost = 0.0_dp !! 0.5 * sum(fun**2)
        real(dp) :: optimality = 0.0_dp !! maximum KKT violation in gradient units
        integer :: nit = 0 !! active-set iterations (not inner linear algebra steps)
        integer :: status = linear_lsq_invalid !! 1 optimal, 3 unconstrained feasible, negative failure
        logical :: success = .false. !! true for valid KKT point
        character(len=160) :: message = 'invalid input' !! solver termination message
    end type linear_lsq_result

    public :: nnls, lsq_linear

contains

    subroutine nnls(a, b, x, rnorm, status, maxiter, tol)
        real(dp), intent(in) :: a(:, :) !! m-by-n coefficient matrix
        real(dp), intent(in) :: b(:) !! length-m right-hand side
        real(dp), intent(out) :: x(:) !! length-n nonnegative solution
        real(dp), intent(out) :: rnorm !! Euclidean residual norm
        integer, intent(out) :: status !! 0 success; negative for invalid/nonconvergence
        integer, intent(in), optional :: maxiter !! active-set iteration limit
        real(dp), intent(in), optional :: tol !! KKT stopping tolerance
        type(linear_lsq_result) :: res
        real(dp), allocatable :: lower(:), upper(:)
        real(dp) :: inf
        integer :: n

        x = 0.0_dp
        rnorm = 0.0_dp
        status = linear_lsq_invalid
        n = size(a,2)
        if (size(x) /= n) return
        inf = ieee_value(1.0_dp, ieee_positive_inf)
        allocate(lower(n), upper(n))
        lower = 0.0_dp
        upper = inf
        res = lsq_linear(a, b, lower, upper, tol=tol, maxiter=maxiter)
        if (allocated(res%x)) then
            x = res%x
            rnorm = norm2(res%fun)
        end if
        if (res%success) then
            status = linear_lsq_success
        else
            status = res%status
        end if
    end subroutine nnls

    function lsq_linear(a, b, lower, upper, tol, maxiter, method) result(res)
        real(dp), intent(in) :: a(:, :) !! m-by-n coefficient matrix
        real(dp), intent(in) :: b(:) !! length-m right-hand side
        real(dp), intent(in), optional :: lower(:) !! n lower bounds, default -infinity
        real(dp), intent(in), optional :: upper(:) !! n upper bounds, default +infinity
        real(dp), intent(in), optional :: tol !! relative KKT tolerance, default 1e-10
        integer, intent(in), optional :: maxiter !! outer active-set limit, default max(100,20*n)
        character(len=*), intent(in), optional :: method !! 'bvls' only (TRF not yet implemented)
        type(linear_lsq_result) :: res

        real(dp), allocatable :: lo(:), hi(:), xfree(:,:), xfull(:,:), &
            rhs(:,:), a_free(:,:), candidate(:), direction(:)
        integer, allocatable :: free_ids(:)
        integer :: m, n, k, j, nfree, info, its, limit, entered
        integer :: hit
        real(dp) :: inf, threshold, scale, worst, violation, alpha, alpha_j, &
            optimality, atolfac, rtol, bound_eps
        logical :: needs_refit, feasible

        res%status = linear_lsq_invalid
        res%message = 'invalid input'
        m = size(a,1)
        n = size(a,2)
        if (size(b) /= m) return
        if (present(lower)) then
            if (size(lower) /= n) return
        end if
        if (present(upper)) then
            if (size(upper) /= n) return
        end if
        if (present(method)) then
            if (trim(method) /= 'bvls') then
                res%message = 'only method=bvls is implemented'
                return
            end if
        end if
        rtol = 1.0e-10_dp
        if (present(tol)) rtol = tol
        if (.not. ieee_is_finite(rtol)) return
        if (rtol < 0.0_dp) return
        limit = max(100, 20*n)
        if (present(maxiter)) limit = maxiter
        if (limit < 0) return
        if (.not. all(ieee_is_finite(a)) .or. .not. all(ieee_is_finite(b))) return

        inf = ieee_value(1.0_dp, ieee_positive_inf)
        allocate(lo(n), hi(n))
        lo = -inf
        hi = inf
        if (present(lower)) lo = lower
        if (present(upper)) hi = upper
        if (any(ieee_is_nan(lo)) .or. any(ieee_is_nan(hi))) return
        if (any(lo > hi) .or. any(lo >= inf) .or. any(hi <= -inf)) return

        allocate(res%x(n), res%fun(m), res%grad(n), res%active_mask(n))
        res%x = 0.0_dp
        res%active_mask = 0
        res%fun = -b
        res%grad = 0.0_dp

        if (n == 0 .or. m == 0) then
            res%x = max(lo,min(hi,res%x))
            call finish_metrics()
            res%status = 3
            res%success = .true.
            res%message = 'empty system, feasible minimum'
            return
        end if

        ! Unconstrained minimum-norm solution. SciPy reports status 3 if it
        ! already satisfies the constraints. Using SVD handles wide/rank-
        ! deficient systems without forming normal equations.
        allocate(xfull(n,1), rhs(m,1))
        rhs(:,1) = b
        call svd_lstsq(a, rhs, xfull, info)
        if (info /= linalg_status_success) then
            res%status = linear_lsq_numerical
            res%message = 'unconstrained SVD did not converge'
            return
        end if
        scale = max(1.0_dp, maxval(abs(xfull(:,1))))
        bound_eps = 32.0_dp*epsilon(1.0_dp)*scale
        feasible = all(xfull(:,1) >= lo - bound_eps) .and. &
            all(xfull(:,1) <= hi + bound_eps)
        if (feasible) then
            res%x = max(lo,min(hi,xfull(:,1)))
            call finish_metrics()
            if (res%optimality <= rtol*max(1.0_dp,maxval(abs(matmul(transpose(a),b)))) &
                    + bound_eps) then
                res%status = 3
                res%success = .true.
                res%message = 'unconstrained solution is optimal'
                return
            end if
        end if

        res%x = max(lo,min(hi,xfull(:,1)))
        do j = 1,n
            if (res%x(j) <= lo(j)) res%active_mask(j) = -1
            if (res%x(j) >= hi(j)) res%active_mask(j) = 1
            if (lo(j) >= hi(j)) res%active_mask(j) = -1
        end do
        atolfac = max(1.0_dp,maxval(abs(matmul(transpose(a),b))))
        threshold = max(16.0_dp*epsilon(1.0_dp)*atolfac, rtol*atolfac)
        needs_refit = .true.
        res%status = linear_lsq_maxiter
        res%message = 'active-set iteration limit reached'
        do its = 0,limit
            if (needs_refit) then
                ! Solve the unconstrained least-squares problem over the free
                ! variables, keeping active variables fixed.
                do
                    nfree = count(res%active_mask == 0)
                    if (nfree == 0) exit
                    free_ids = pack([(j,j=1,n)],res%active_mask == 0)
                    allocate(a_free(m,nfree),xfree(nfree,1),candidate(nfree),direction(nfree))
                    a_free = a(:,free_ids)
                    rhs(:,1) = b
                    do j=1,n
                        if (res%active_mask(j) /= 0) rhs(:,1) = rhs(:,1)-a(:,j)*res%x(j)
                    end do
                    call svd_lstsq(a_free,rhs,xfree,info)
                    if (info /= linalg_status_success) then
                        res%status = linear_lsq_numerical
                        res%message = 'free-variable least-squares solve failed'
                        return
                    end if
                    candidate = xfree(:,1)
                    direction = candidate-res%x(free_ids)
                    alpha = 1.0_dp
                    hit = 0
                    do k=1,nfree
                        j = free_ids(k)
                        if (direction(k) > 0.0_dp .and. candidate(k) > hi(j)) then
                            alpha_j = (hi(j)-res%x(j))/direction(k)
                            if (alpha_j < alpha) then
                                alpha = max(0.0_dp,alpha_j)
                                hit = j
                            end if
                        else if (direction(k) < 0.0_dp .and. candidate(k) < lo(j)) then
                            alpha_j = (lo(j)-res%x(j))/direction(k)
                            if (alpha_j < alpha) then
                                alpha = max(0.0_dp,alpha_j)
                                hit = j
                            end if
                        end if
                    end do
                    res%x(free_ids) = res%x(free_ids)+alpha*direction
                    res%x(free_ids) = max(lo(free_ids),min(hi(free_ids),res%x(free_ids)))
                    if (hit /= 0) then
                        if (direction(findloc(free_ids,hit,dim=1)) > 0.0_dp) then
                            res%x(hit)=hi(hit)
                            res%active_mask(hit)=1
                        else
                            res%x(hit)=lo(hit)
                            res%active_mask(hit)=-1
                        end if
                    end if
                    deallocate(a_free,xfree,candidate,direction)
                    if (hit == 0) exit
                    if (its == limit) exit
                end do
            end if
            call finish_metrics()
            if (.not. ieee_is_finite(res%cost)) then
                res%status=linear_lsq_numerical
                res%message='nonfinite cost'
                return
            end if
            worst = threshold
            entered = 0
            do j=1,n
                if (lo(j) >= hi(j)) cycle
                violation = 0.0_dp
                if (res%active_mask(j) == -1) violation = -res%grad(j)
                if (res%active_mask(j) == 1) violation = res%grad(j)
                if (violation > worst) then
                    worst = violation
                    entered = j
                end if
            end do
            if (entered == 0 .and. res%optimality <= threshold) then
                res%status=1
                res%message='first-order KKT conditions satisfied'
                res%success=.true.
                return
            end if
            if (entered == 0) then
                res%status=linear_lsq_numerical
                res%message='free-variable KKT residual too large'
                return
            end if
            if (its == limit) exit
            res%active_mask(entered) = 0
            res%nit = res%nit+1
            needs_refit=.true.
        end do
        call finish_metrics()

    contains
        subroutine finish_metrics()
            integer :: p
            real(dp) :: v
            res%fun = matmul(a,res%x)-b
            res%grad = matmul(transpose(a),res%fun)
            res%cost = 0.5_dp*sum(res%fun**2)
            optimality = 0.0_dp
            do p=1,n
                v = abs(res%grad(p))
                if (res%active_mask(p) == -1) v = max(0.0_dp,-res%grad(p))
                if (res%active_mask(p) == 1) v = max(0.0_dp,res%grad(p))
                if (lo(p) >= hi(p)) v=0.0_dp
                optimality = max(optimality,v)
            end do
            res%optimality = optimality
        end subroutine finish_metrics
    end function lsq_linear

end module scifort_optimize_linear_least_squares
