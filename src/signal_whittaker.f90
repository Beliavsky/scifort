! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Whittaker-Henderson smoothing for a fixed penalty or automatic REML penalty.
! Independent banded Cholesky implementation of the difference-penalty model.
! The profiled REML criterion follows Biessy (2023) and the algorithmic
! description in scipy.signal.whittaker_henderson (SciPy 1.18, BSD-3-Clause).
! See CODE_PROVENANCE.md for attribution.
module scifort_signal_whittaker
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: whittaker_status_success = 0
    integer, parameter, public :: whittaker_status_invalid_input = 1
    integer, parameter, public :: whittaker_status_singular = 2

    type, public :: whittaker_result
        real(dp), allocatable :: x(:) !! smoothed samples, empty on failure
        real(dp) :: lamb = 0.0_dp !! supplied or REML-optimized penalty
        integer :: status = whittaker_status_invalid_input !! zero on success
    end type whittaker_result

    public :: whittaker_henderson
contains

    function whittaker_henderson(signal, lamb, order, weights) result(out)
        real(dp), intent(in) :: signal(:) !! observations, equally spaced
        real(dp), intent(in), optional :: lamb !! nonnegative penalty, absent selects REML
        integer, intent(in), optional :: order !! difference order, default 2
        real(dp), intent(in), optional :: weights(:) !! nonnegative weights; zero permits NaN
        type(whittaker_result) :: out
        real(dp), allocatable :: w(:)
        real(dp) :: penalty, logdet
        integer :: n, p

        allocate(out%x(0))
        n = size(signal)
        p = 2
        if (present(order)) p = order
        ! High-order differences are ill-conditioned in binary64 arithmetic.
        if (p < 1 .or. p > 12 .or. n <= p) return
        allocate(w(n))
        w = 1.0_dp
        if (present(weights)) then
            if (size(weights) /= n) return
            if (any(.not. ieee_is_finite(weights))) return
            if (any(weights < 0.0_dp)) return
            w = weights
        end if
        ! Missing observations are allowed only with exactly zero weight.
        if (any(.not. ieee_is_finite(signal) .and. w > 0.0_dp)) return
        if (present(lamb)) then
            penalty = lamb
        else
            call choose_penalty_reml(signal,p,w,penalty,out%status)
            if (out%status /= whittaker_status_success) return
        end if
        out%lamb = penalty
        if (.not. ieee_is_finite(penalty) .or. penalty < 0.0_dp) then
            out%status = whittaker_status_invalid_input
            return
        end if
        call solve_whittaker(signal,penalty,p,w,out%x,logdet,out%status)
    end function whittaker_henderson

    subroutine solve_whittaker(signal,lamb,p,w,x,logdet,status)
        real(dp), intent(in) :: signal(:) !! samples; nonfinite permitted at zero weight
        real(dp), intent(in) :: lamb !! penalty, nonnegative
        integer, intent(in) :: p !! difference order
        real(dp), intent(in) :: w(:) !! observation weights
        real(dp), allocatable, intent(out) :: x(:) !! solution
        real(dp), intent(out) :: logdet !! log determinant of W + lamb * D^T D
        integer, intent(out) :: status !! zero on success
        real(dp), allocatable :: band(:,:), rhs(:), coeff(:)
        real(dp) :: pivot, sum_term, coeff_j
        integer :: n, j, k, i, l, lower

        allocate(x(0))
        logdet = 0.0_dp
        status = whittaker_status_singular
        n = size(signal)
        allocate(band(0:p,n),rhs(n),coeff(0:p))
        band = 0.0_dp
        rhs = 0.0_dp
        do i=1,n
            band(0,i) = w(i)
            if (w(i) > 0.0_dp) rhs(i) = w(i)*signal(i)
        end do
        ! Forward finite-difference coefficients with alternating signs.
        coeff(0) = real(merge(1,-1,mod(p,2)==0),dp)
        do j=1,p
            coeff(j) = -coeff(j-1)*real(p-j+1,dp)/real(j,dp)
        end do
        if (lamb > 0.0_dp) then
            do k=1,n-p
                do j=0,p
                    coeff_j = lamb*coeff(j)
                    do l=0,j
                        band(j-l,k+j) = band(j-l,k+j)+coeff_j*coeff(l)
                    end do
                end do
            end do
        end if
        ! Factor symmetric positive definite band matrix in place.
        do i=1,n
            lower = max(1,i-p)
            do j=lower,i-1
                sum_term = band(i-j,i)
                do k=max(lower,j-p),j-1
                    sum_term = sum_term-band(i-k,i)*band(j-k,j)
                end do
                band(i-j,i) = sum_term/band(0,j)
            end do
            pivot = band(0,i)
            do k=lower,i-1
                pivot = pivot-band(i-k,i)**2
            end do
            if (.not. ieee_is_finite(pivot) .or. pivot <= 0.0_dp) return
            band(0,i) = sqrt(pivot)
            logdet = logdet+2.0_dp*log(band(0,i))
        end do
        ! Forward and backward triangular solves.
        do i=1,n
            do k=max(1,i-p),i-1
                rhs(i) = rhs(i)-band(i-k,i)*rhs(k)
            end do
            rhs(i) = rhs(i)/band(0,i)
        end do
        do i=n,1,-1
            do k=i+1,min(n,i+p)
                rhs(i) = rhs(i)-band(k-i,k)*rhs(k)
            end do
            rhs(i) = rhs(i)/band(0,i)
        end do
        if (any(.not. ieee_is_finite(rhs))) return
        x = rhs
        status = whittaker_status_success
    end subroutine solve_whittaker

    real(dp) function reml_cost(log_lamb,signal,p,w) result(cost)
        real(dp), intent(in) :: log_lamb !! natural logarithm of the penalty
        real(dp), intent(in) :: signal(:) !! observed samples
        integer, intent(in) :: p !! difference penalty order
        real(dp), intent(in) :: w(:) !! observation weights
        real(dp), allocatable :: x(:), coeff(:)
        real(dp) :: lamb, logdet, residual, difference, r2
        integer :: i, j, k, n, status

        cost = huge(1.0_dp)
        n = size(signal)
        lamb = exp(log_lamb)
        call solve_whittaker(signal,lamb,p,w,x,logdet,status)
        if (status /= whittaker_status_success) return
        r2 = 0.0_dp
        do i=1,n
            if (w(i) > 0.0_dp) then
                residual = signal(i)-x(i)
                r2 = r2+w(i)*residual**2
            end if
        end do
        allocate(coeff(0:p))
        coeff(0) = real(merge(1,-1,mod(p,2)==0),dp)
        do j=1,p
            coeff(j) = -coeff(j-1)*real(p-j+1,dp)/real(j,dp)
        end do
        do k=1,n-p
            difference = 0.0_dp
            do j=0,p
                difference = difference+coeff(j)*x(k+j)
            end do
            r2 = r2+lamb*difference**2
        end do
        if (.not. ieee_is_finite(r2) .or. r2 <= 0.0_dp) return
        ! Profiled restricted log likelihood, omitting constants independent
        ! of lambda: (n-p)*log(r2/(n-p)) + log|W+lambda D'D| - (n-p)*log(lambda).
        cost = real(n-p,dp)*(log(r2/real(n-p,dp))-log_lamb)+logdet
        if (.not. ieee_is_finite(cost)) cost = huge(1.0_dp)
    end function reml_cost

    subroutine choose_penalty_reml(signal,p,w,lamb,status)
        real(dp), intent(in) :: signal(:) !! observations
        integer, intent(in) :: p !! finite difference order
        real(dp), intent(in) :: w(:) !! case weights
        real(dp), intent(out) :: lamb !! optimized smoothing parameter
        integer, intent(out) :: status !! zero on success
        integer, parameter :: ngrid = 33
        real(dp), parameter :: left_bound = -24.0_dp, right_bound = 24.0_dp
        real(dp), parameter :: golden = 0.6180339887498948482_dp
        real(dp) :: grid(ngrid), cost(ngrid), left, right, a, b, fa, fb, optimum
        integer :: j, imin, iter

        lamb = 0.0_dp
        status = whittaker_status_singular
        do j=1,ngrid
            grid(j) = left_bound+real(j-1,dp)*(right_bound-left_bound)/real(ngrid-1,dp)
            cost(j) = reml_cost(grid(j),signal,p,w)
        end do
        imin = minloc(cost,dim=1)
        if (cost(imin) == huge(1.0_dp)) return
        ! Bracket a minimum in log(lambda), then use golden-section refinement.
        left = grid(max(1,imin-1))
        right = grid(min(ngrid,imin+1))
        a = right-golden*(right-left)
        b = left+golden*(right-left)
        fa = reml_cost(a,signal,p,w)
        fb = reml_cost(b,signal,p,w)
        do iter=1,70
            if (right-left <= 2.0e-9_dp) exit
            if (fa < fb) then
                right = b
                b = a
                fb = fa
                a = right-golden*(right-left)
                fa = reml_cost(a,signal,p,w)
            else
                left = a
                a = b
                fa = fb
                b = left+golden*(right-left)
                fb = reml_cost(b,signal,p,w)
            end if
        end do
        optimum = 0.5_dp*(left+right)
        if (reml_cost(optimum,signal,p,w) == huge(1.0_dp)) return
        lamb = exp(optimum)
        status = whittaker_status_success
    end subroutine choose_penalty_reml
end module scifort_signal_whittaker
