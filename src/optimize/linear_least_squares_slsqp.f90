! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Thin adapter to the BSD-style Lawson-Hanson BVLS implementation exposed
! by SciFort's already-pinned jacobwilliams/slsqp 1.6.0 FPM dependency.
! No third-party source is copied into SciFort. The dependency retains its
! original license and notices. Use a default real64 SLSQP build.
module scifort_optimize_linear_least_squares_slsqp
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_value, ieee_positive_inf
    use scifort_kinds, only : dp
    use slsqp_kinds, only : wp
    use bvls_module, only : bvls
    use scifort_optimize_linear_least_squares, only : &
        linear_lsq_result, linear_lsq_success, linear_lsq_invalid, &
        linear_lsq_maxiter, linear_lsq_numerical
    implicit none
    private
    public :: nnls_slsqp, lsq_linear_slsqp

contains

    subroutine nnls_slsqp(a,b,x,rnorm,status,maxiter)
        real(dp), intent(in) :: a(:, :) !! m-by-n real coefficient matrix
        real(dp), intent(in) :: b(:) !! length-m right-hand side
        real(dp), intent(out) :: x(:) !! length-n solution constrained x>=0
        real(dp), intent(out) :: rnorm !! residual norm
        integer, intent(out) :: status !! zero success, negative error
        integer, intent(in), optional :: maxiter !! BVLS iteration limit
        real(dp), allocatable :: lo(:),hi(:)
        type(linear_lsq_result) :: result
        integer :: n
        x=0.0_dp
        rnorm=0.0_dp
        status=linear_lsq_invalid
        n=size(a,2)
        if (size(x) /= n) return
        allocate(lo(n),hi(n))
        lo=0.0_dp
        hi=ieee_value(1.0_dp,ieee_positive_inf)
        result=lsq_linear_slsqp(a,b,lo,hi,maxiter)
        if (allocated(result%x)) then
            x=result%x
            rnorm=norm2(result%fun)
        end if
        if (result%success) then
            status=linear_lsq_success
        else
            status=result%status
        end if
    end subroutine nnls_slsqp

    function lsq_linear_slsqp(a,b,lower,upper,maxiter) result(res)
        real(dp), intent(in) :: a(:, :) !! m-by-n input matrix
        real(dp), intent(in) :: b(:) !! length-m target
        real(dp), intent(in) :: lower(:) !! n lower bounds (may contain -infinity)
        real(dp), intent(in) :: upper(:) !! n upper bounds (may contain +infinity)
        integer, intent(in), optional :: maxiter !! BVLS iteration cap
        type(linear_lsq_result) :: res
        real(wp), allocatable :: awork(:,:),bwork(:),bnd(:,:),xwork(:),w(:)
        integer, allocatable :: index(:)
        integer :: m,n,ierr,nsetp,limit,j
        real(wp) :: rnorm

        m=size(a,1)
        n=size(a,2)
        if (size(b) /= m .or. size(lower) /= n .or. size(upper) /= n) return
        if (m<1 .or. n<1) return ! For empty systems use the native implementation.
        if (.not.all(ieee_is_finite(a)) .or. .not.all(ieee_is_finite(b))) return
        if (any(ieee_is_nan(lower)) .or. any(ieee_is_nan(upper))) return
        if (any(lower>upper)) return
        if (any(lower >= ieee_value(1.0_dp,ieee_positive_inf))) return
        if (any(upper <= -ieee_value(1.0_dp,ieee_positive_inf))) return
        limit=3*n
        if (present(maxiter)) limit=maxiter
        if (limit<0) return

        allocate(awork(m,n),bwork(m),bnd(2,n),xwork(n),w(n),index(n))
        awork=real(a,wp)
        bwork=real(b,wp)
        bnd(1,:)=real(lower,wp)
        bnd(2,:)=real(upper,wp)
        where (.not.ieee_is_finite(lower)) bnd(1,:)=-huge(1.0_wp)
        where (.not.ieee_is_finite(upper)) bnd(2,:)=huge(1.0_wp)
        call bvls(awork,bwork,bnd,xwork,rnorm,nsetp,w,index,ierr,limit)

        allocate(res%x(n),res%fun(m),res%grad(n),res%active_mask(n))
        res%x=0.0_dp
        if (ierr == 0 .or. ierr == 4) then
            if (.not.all(ieee_is_finite(xwork))) then
                res%status=linear_lsq_numerical
                res%message='SLSQP BVLS returned nonfinite coefficients'
                return
            end if
            res%x=real(xwork,dp)
        end if
        res%fun=matmul(a,res%x)-b
        res%grad=matmul(transpose(a),res%fun)
        res%cost=0.5_dp*sum(res%fun**2)
        res%active_mask=0
        res%optimality=0.0_dp
        do j=1,n
            if (res%x(j)<=lower(j)) res%active_mask(j)=-1
            if (res%x(j)>=upper(j)) res%active_mask(j)=1
            if (lower(j)>=upper(j)) res%active_mask(j)=-1
        end do
        do j=1,n
            if (res%active_mask(j)==-1) then
                res%optimality=max(res%optimality,max(0.0_dp,-res%grad(j)))
            else if (res%active_mask(j)==1) then
                res%optimality=max(res%optimality,max(0.0_dp,res%grad(j)))
            else
                res%optimality=max(res%optimality,abs(res%grad(j)))
            end if
        end do
        res%nit=-1 ! The dependency does not expose its iteration count.
        select case(ierr)
        case(0)
            res%status=1
            res%success=.true.
            res%message='SLSQP BVLS completed'
        case(4)
            res%status=linear_lsq_maxiter
            res%message='SLSQP BVLS maximum iterations'
        case default
            res%status=linear_lsq_numerical
            res%message='SLSQP BVLS failed'
        end select
    end function lsq_linear_slsqp

end module scifort_optimize_linear_least_squares_slsqp
