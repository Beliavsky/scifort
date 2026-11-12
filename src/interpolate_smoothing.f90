! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Natural cubic smoothing splines, including generalized cross validation.
! Independent O(n) Reinsch / roughness-penalty formulation.  No FITPACK,
! BLAS, or LAPACK dependency.  Smoothing weights are NOT squared.
module scifort_interpolate_smoothing
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_interpolate_piecewise, only : ppoly_type, ppoly
    implicit none
    private

    integer, parameter, public :: smoothing_status_success = 0
    integer, parameter, public :: smoothing_status_invalid_input = 1
    integer, parameter, public :: smoothing_status_singular = 2

    public :: make_smoothing_spline
    interface make_smoothing_spline
        module procedure smoothing_spline_scalar
        module procedure smoothing_spline_channels
    end interface
contains

    function smoothing_spline_scalar(x, y, w, lam, status, selected_lam) result(spline)
        real(dp), intent(in) :: x(:) !! increasing x coordinates, length n >= 5
        real(dp), intent(in) :: y(:) !! scalar observations at x, length n
        real(dp), intent(in), optional :: w(:) !! positive weights, length n; default one
        real(dp), intent(in), optional :: lam !! smoothing penalty, >=0; absent selects GCV optimum
        integer, intent(out), optional :: status !! 0 success, 1 invalid inputs, 2 singular solve
        real(dp), intent(out), optional :: selected_lam !! penalty chosen by GCV, or supplied lam
        type(ppoly_type) :: spline
        real(dp), allocatable :: coeff(:,:,:), y_matrix(:,:), lambdas(:)
        integer :: code

        allocate(y_matrix(size(y),1))
        y_matrix(:,1) = y
        call make_coefficients(x,y_matrix,w,lam,coeff,lambdas,code)
        if (code == smoothing_status_success) spline = ppoly(coeff,x,status=code)
        if (present(status)) status = code
        if (present(selected_lam)) then
            selected_lam = 0.0_dp
            if (code == smoothing_status_success) selected_lam = lambdas(1)
        end if
    end function smoothing_spline_scalar

    function smoothing_spline_channels(x, y, w, lam, status, selected_lam) result(spline)
        real(dp), intent(in) :: x(:) !! increasing x coordinates, length n >= 5
        real(dp), intent(in) :: y(:,:) !! observations (n, nchannels), matching SciPy axis=0
        real(dp), intent(in), optional :: w(:) !! shared positive point weights, length n
        real(dp), intent(in), optional :: lam !! common penalty, >=0; absent chooses GCV per channel
        integer, intent(out), optional :: status !! 0 success, 1 invalid inputs, 2 singular solve
        real(dp), intent(out), optional :: selected_lam(:) !! chosen penalty for each channel
        type(ppoly_type) :: spline
        real(dp), allocatable :: coeff(:,:,:), lambdas(:)
        integer :: code

        code = smoothing_status_invalid_input
        if (present(selected_lam)) then
            if (size(selected_lam) /= size(y,2)) goto 900
        end if
        call make_coefficients(x,y,w,lam,coeff,lambdas,code)
        if (code == smoothing_status_success) spline = ppoly(coeff,x,status=code)
900     if (present(status)) status = code
        if (present(selected_lam)) then
            selected_lam = 0.0_dp
            if (code == smoothing_status_success) selected_lam = lambdas
        end if
    end function smoothing_spline_channels

    subroutine make_coefficients(x,y,w_in,lam_in,coeff,lambdas,status)
        real(dp), intent(in) :: x(:) !! strictly increasing finite observation sites
        real(dp), intent(in) :: y(:,:) !! finite observations (n, channels)
        real(dp), intent(in), optional :: w_in(:) !! positive observation weights
        real(dp), intent(in), optional :: lam_in !! nonnegative penalty; absent uses GCV
        real(dp), allocatable, intent(out) :: coeff(:,:,:) !! descending cubic coefficients (4,n-1,channels)
        real(dp), allocatable, intent(out) :: lambdas(:) !! chosen penalties per channel
        integer, intent(out) :: status !! success, invalid, or singular
        real(dp), allocatable :: h(:), w(:), q(:,:), b0(:), b1(:), b2(:), r0(:), r1(:)
        real(dp), allocatable :: rhs(:), m(:), fit(:), curv(:), residual(:)
        real(dp) :: t, lambda
        integer :: n, nc, msize, i, j

        status = smoothing_status_invalid_input
        n = size(x)
        nc = size(y,2)
        if (n < 5 .or. nc < 1 .or. size(y,1) /= n) return
        if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(y))) return
        if (.not. all(x(2:) > x(:n-1))) return
        if (present(w_in)) then
            if (size(w_in) /= n) return
            if (.not. all(ieee_is_finite(w_in))) return
            if (.not. all(w_in > 0.0_dp)) return
            w = w_in
        else
            allocate(w(n))
            w = 1.0_dp
        end if
        if (present(lam_in)) then
            if (.not. ieee_is_finite(lam_in)) return
            if (lam_in < 0.0_dp) return
        end if
        h = x(2:) - x(:n-1)
        msize = n-2
        allocate(q(3,msize), b0(msize), b1(max(0,msize-1)), b2(max(0,msize-2)), &
            r0(msize), r1(max(0,msize-1)),rhs(msize),m(msize),fit(n),curv(n),residual(n))
        q(1,:) = 1.0_dp/h(:n-2)
        q(2,:) = -1.0_dp/h(:n-2) - 1.0_dp/h(2:)
        q(3,:) = 1.0_dp/h(2:)
        r0 = (h(:n-2)+h(2:))/3.0_dp
        r1 = h(2:n-2)/6.0_dp
        do i=1,msize
            b0(i) = q(1,i)**2/w(i) + q(2,i)**2/w(i+1) + q(3,i)**2/w(i+2)
        end do
        do i=1,msize-1
            b1(i) = q(2,i)*q(1,i+1)/w(i+1) + q(3,i)*q(2,i+1)/w(i+2)
        end do
        do i=1,msize-2
            b2(i) = q(3,i)*q(1,i+2)/w(i+2)
        end do
        if (.not. all(ieee_is_finite(b0)) .or. .not. all(ieee_is_finite(b1)) &
            .or. .not. all(ieee_is_finite(b2))) return
        allocate(coeff(4,n-1,nc),lambdas(nc))
        do j=1,nc
            do i=1,msize
                rhs(i) = sum(q(:,i)*y(i:i+2,j))
            end do
            lambda = 0.0_dp
            if (present(lam_in)) then
                lambda = lam_in
            else
                call minimize_gcv(real(n,dp),y(:,j),w,q,r0,r1,b0,b1,b2,rhs,lambda,status)
                if (status /= smoothing_status_success) return
            end if
            lambdas(j) = lambda
            call solve_penalty(lambda,r0,r1,b0,b1,b2,rhs,m,status)
            if (status /= smoothing_status_success) return
            residual=0.0_dp
            do i=1,msize
                residual(i:i+2) = residual(i:i+2) + q(:,i)*m(i)
            end do
            fit = y(:,j) - lambda*residual/w
            if (.not. all(ieee_is_finite(fit))) then
                status=smoothing_status_singular
                return
            end if
            curv=0.0_dp
            curv(2:n-1)=m
            do i=1,n-1
                t=h(i)
                coeff(1,i,j)=(curv(i+1)-curv(i))/(6.0_dp*t)
                coeff(2,i,j)=curv(i)/2.0_dp
                coeff(3,i,j)=(fit(i+1)-fit(i))/t - t*(2.0_dp*curv(i)+curv(i+1))/6.0_dp
                coeff(4,i,j)=fit(i)
            end do
        end do
        status=smoothing_status_success
    end subroutine make_coefficients

    subroutine factor_penalty(lambda,r0,r1,b0,b1,b2,d,l1,l2,status)
        real(dp), intent(in) :: lambda !! smoothing penalty
        real(dp), intent(in) :: r0(:) !! R main diagonal
        real(dp), intent(in) :: r1(:) !! R first upper diagonal
        real(dp), intent(in) :: b0(:) !! Q^T W^-1 Q main diagonal
        real(dp), intent(in) :: b1(:) !! Q^T W^-1 Q first upper diagonal
        real(dp), intent(in) :: b2(:) !! Q^T W^-1 Q second upper diagonal
        real(dp), intent(out) :: d(:) !! Cholesky diagonal (positive)
        real(dp), intent(out) :: l1(:) !! Cholesky first lower diagonal, l1(1)=0
        real(dp), intent(out) :: l2(:) !! Cholesky second lower diagonal, l2(1:2)=0
        integer, intent(out) :: status !! 0 for positive definite, 2 otherwise
        real(dp) :: pivot
        integer :: i, n

        n=size(r0)
        d=0.0_dp
        l1=0.0_dp
        l2=0.0_dp
        status=smoothing_status_singular
        do i=1,n
            if (i>2) l2(i)=lambda*b2(i-2)/d(i-2)
            if (i>1) then
                l1(i)=r1(i-1)+lambda*b1(i-1)
                if (i>2) l1(i)=l1(i)-l2(i)*l1(i-1)
                l1(i)=l1(i)/d(i-1)
            end if
            pivot=r0(i)+lambda*b0(i)-l1(i)*l1(i)-l2(i)*l2(i)
            if (.not. ieee_is_finite(pivot)) return
            if (pivot<=0.0_dp) return
            d(i)=sqrt(pivot)
        end do
        status=smoothing_status_success
    end subroutine factor_penalty

    subroutine solve_factored(d,l1,l2,rhs,x)
        real(dp), intent(in) :: d(:) !! lower Cholesky diagonal
        real(dp), intent(in) :: l1(:) !! first subdiagonal
        real(dp), intent(in) :: l2(:) !! second subdiagonal
        real(dp), intent(in) :: rhs(:) !! system right-hand side
        real(dp), intent(out) :: x(:) !! solution to L L^T x = rhs
        real(dp) :: y(size(rhs)), z
        integer :: i, n

        n=size(rhs)
        do i=1,n
            z=rhs(i)
            if (i>1) z=z-l1(i)*y(i-1)
            if (i>2) z=z-l2(i)*y(i-2)
            y(i)=z/d(i)
        end do
        do i=n,1,-1
            z=y(i)
            if (i<n) z=z-l1(i+1)*x(i+1)
            if (i<n-1) z=z-l2(i+2)*x(i+2)
            x(i)=z/d(i)
        end do
    end subroutine solve_factored

    subroutine solve_penalty(lambda,r0,r1,b0,b1,b2,rhs,x,status)
        real(dp), intent(in) :: lambda !! smoothing penalty
        real(dp), intent(in) :: r0(:) !! R main diagonal
        real(dp), intent(in) :: r1(:) !! R offdiagonal
        real(dp), intent(in) :: b0(:) !! quadratic penalty main diagonal
        real(dp), intent(in) :: b1(:) !! quadratic penalty first offdiagonal
        real(dp), intent(in) :: b2(:) !! quadratic penalty second offdiagonal
        real(dp), intent(in) :: rhs(:) !! system rhs
        real(dp), intent(out) :: x(:) !! fitted inner second derivatives
        integer, intent(out) :: status !! solve status
        real(dp) :: d(size(rhs)),l1(size(rhs)),l2(size(rhs))

        call factor_penalty(lambda,r0,r1,b0,b1,b2,d,l1,l2,status)
        if (status/=smoothing_status_success) return
        call solve_factored(d,l1,l2,rhs,x)
        if (.not. all(ieee_is_finite(x))) status=smoothing_status_singular
    end subroutine solve_penalty

    function gcv_score(lambda,y,w,q,r0,r1,b0,b1,b2,rhs,status) result(score)
        real(dp), intent(in) :: lambda !! candidate penalty
        real(dp), intent(in) :: y(:) !! observed data (used for size)
        real(dp), intent(in) :: w(:) !! positive observation weights
        real(dp), intent(in) :: q(:,:) !! three-point second-difference rows
        real(dp), intent(in) :: r0(:) !! R diagonal
        real(dp), intent(in) :: r1(:) !! R first upper diagonal
        real(dp), intent(in) :: b0(:) !! penalty main diagonal
        real(dp), intent(in) :: b1(:) !! penalty first upper diagonal
        real(dp), intent(in) :: b2(:) !! penalty second upper diagonal
        real(dp), intent(in) :: rhs(:) !! second-difference RHS
        integer, intent(out) :: status !! 0 for successful evaluation
        real(dp) :: score
        real(dp) :: d(size(rhs)),l1(size(rhs)),l2(size(rhs)),m(size(rhs)),r(size(y))
        real(dp) :: zd(size(rhs)), z1(size(rhs)), z2(size(rhs)),u1,u2,trace_ab
        integer :: i,n

        score=huge(1.0_dp)
        call factor_penalty(lambda,r0,r1,b0,b1,b2,d,l1,l2,status)
        if (status/=smoothing_status_success) return
        call solve_factored(d,l1,l2,rhs,m)
        n=size(rhs)
        zd=0.0_dp
        z1=0.0_dp
        z2=0.0_dp
        do i=n,1,-1
            u1=0.0_dp
            u2=0.0_dp
            if (i<n) u1=l1(i+1)/d(i)
            if (i<n-1) u2=l2(i+2)/d(i)
            if (i<n-1) z2(i)=-u1*z1(i+1)-u2*zd(i+2)
            if (i<n) z1(i)=-u1*zd(i+1)-u2*z1(i+1)
            zd(i)=1.0_dp/d(i)**2-u1*z1(i)-u2*z2(i)
        end do
        trace_ab = sum(zd*b0) + 2.0_dp*sum(z1(:n-1)*b1) + &
            2.0_dp*sum(z2(:n-2)*b2)
        if (trace_ab<=0.0_dp .or. .not. ieee_is_finite(trace_ab)) then
            status=smoothing_status_singular
            return
        end if
        r=0.0_dp
        do i=1,n
            r(i:i+2)=r(i:i+2)+q(:,i)*m(i)
        end do
        score=real(size(y),dp)*sum((r/w)**2)/(trace_ab**2)
        if (.not. ieee_is_finite(score)) status=smoothing_status_singular
    end function gcv_score

    subroutine minimize_gcv(upper,y,w,q,r0,r1,b0,b1,b2,rhs,lambda,status)
        real(dp), intent(in) :: upper !! upper bound on lambda, SciPy uses n
        real(dp), intent(in) :: y(:) !! data
        real(dp), intent(in) :: w(:) !! weights
        real(dp), intent(in) :: q(:,:) !! second difference coefficients
        real(dp), intent(in) :: r0(:) !! R diagonal
        real(dp), intent(in) :: r1(:) !! R first offdiagonal
        real(dp), intent(in) :: b0(:) !! Q^T W^-1 Q diagonal
        real(dp), intent(in) :: b1(:) !! Q^T W^-1 Q first offdiagonal
        real(dp), intent(in) :: b2(:) !! Q^T W^-1 Q second offdiagonal
        real(dp), intent(in) :: rhs(:) !! second difference right-hand side
        real(dp), intent(out) :: lambda !! minimizing lambda in [0,upper]
        integer, intent(out) :: status !! zero if optimization successful
        real(dp), parameter :: invphi = 0.6180339887498948482_dp
        real(dp) :: a,b,c,d,fc,fd,ft
        integer :: iter

        a=0.0_dp
        b=upper
        c=b-invphi*(b-a)
        d=a+invphi*(b-a)
        fc=gcv_score(c,y,w,q,r0,r1,b0,b1,b2,rhs,status)
        if (status/=smoothing_status_success) return
        fd=gcv_score(d,y,w,q,r0,r1,b0,b1,b2,rhs,status)
        if (status/=smoothing_status_success) return
        do iter=1,70
            if (fc < fd) then
                b=d
                d=c
                fd=fc
                c=b-invphi*(b-a)
                fc=gcv_score(c,y,w,q,r0,r1,b0,b1,b2,rhs,status)
                if (status/=smoothing_status_success) return
            else
                a=c
                c=d
                fc=fd
                d=a+invphi*(b-a)
                fd=gcv_score(d,y,w,q,r0,r1,b0,b1,b2,rhs,status)
                if (status/=smoothing_status_success) return
            end if
            if ((b-a)<2.0e-9_dp*max(1.0_dp,b)) exit
        end do
        lambda=0.5_dp*(a+b)
        ft=gcv_score(lambda,y,w,q,r0,r1,b0,b1,b2,rhs,status)
        if (status/=smoothing_status_success) return
        if (.not. ieee_is_finite(ft)) status=smoothing_status_singular
    end subroutine minimize_gcv
end module scifort_interpolate_smoothing
