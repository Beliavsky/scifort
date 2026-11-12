! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Separable cubic/quadratic spline coefficients in two dimensions.
! Independent Fortran implementation of the symmetric IIR recurrences
! used by scipy.signal.cspline2d / qspline2d / spline_filter.
! See CODE_PROVENANCE.md for SciPy BSD-3-Clause attribution.
module scifort_signal_spline2d
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_spatial, only : sepfir2d
    use scifort_signal_spline, only : spline_status_success, spline_status_invalid_input
    implicit none
    private

    public :: cspline2d, qspline2d, spline_filter, symmetric_iir1, symmetric_iir2

    interface cspline2d
        module procedure cspline2d_real
        module procedure cspline2d_complex
    end interface
    interface qspline2d
        module procedure qspline2d_real
        module procedure qspline2d_complex
    end interface
    interface spline_filter
        module procedure spline_filter_real
        module procedure spline_filter_complex
    end interface

contains

    function cspline2d_real(signal, lamb, precision, status) result(coeff)
        real(dp), intent(in) :: signal(:,:) !! real sampled image
        real(dp), optional, intent(in) :: lamb !! smoothing parameter, default zero
        real(dp), optional, intent(in) :: precision !! mirror-boundary truncation tolerance
        integer, optional, intent(out) :: status !! zero for success
        real(dp), allocatable :: coeff(:,:)
        complex(dp), allocatable :: work(:,:)
        integer :: code

        work = spline2d_core(cmplx(signal, 0.0_dp, dp), 3, lamb, precision, code)
        coeff = real(work, dp)
        if (present(status)) status = code
    end function cspline2d_real

    function cspline2d_complex(signal, lamb, precision, status) result(coeff)
        complex(dp), intent(in) :: signal(:,:) !! complex sampled image
        real(dp), optional, intent(in) :: lamb !! smoothing parameter, default zero
        real(dp), optional, intent(in) :: precision !! mirror-boundary truncation tolerance
        integer, optional, intent(out) :: status !! zero for success
        complex(dp), allocatable :: coeff(:,:)
        integer :: code

        coeff = spline2d_core(signal, 3, lamb, precision, code)
        if (present(status)) status = code
    end function cspline2d_complex

    function qspline2d_real(signal, lamb, precision, status) result(coeff)
        real(dp), intent(in) :: signal(:,:) !! real sampled image
        real(dp), optional, intent(in) :: lamb !! quadratic smoothing must be zero
        real(dp), optional, intent(in) :: precision !! mirror-boundary truncation tolerance
        integer, optional, intent(out) :: status !! zero for success
        real(dp), allocatable :: coeff(:,:)
        complex(dp), allocatable :: work(:,:)
        integer :: code

        work = spline2d_core(cmplx(signal, 0.0_dp, dp), 2, lamb, precision, code)
        coeff = real(work, dp)
        if (present(status)) status = code
    end function qspline2d_real

    function qspline2d_complex(signal, lamb, precision, status) result(coeff)
        complex(dp), intent(in) :: signal(:,:) !! complex sampled image
        real(dp), optional, intent(in) :: lamb !! quadratic smoothing must be zero
        real(dp), optional, intent(in) :: precision !! mirror-boundary truncation tolerance
        integer, optional, intent(out) :: status !! zero for success
        complex(dp), allocatable :: coeff(:,:)
        integer :: code

        coeff = spline2d_core(signal, 2, lamb, precision, code)
        if (present(status)) status = code
    end function qspline2d_complex

    function spline_filter_real(signal, lmbda, status) result(filtered)
        real(dp), intent(in) :: signal(:,:) !! image to smooth
        real(dp), optional, intent(in) :: lmbda !! cubic smoothing parameter, default five
        integer, optional, intent(out) :: status !! zero for success
        real(dp), allocatable :: filtered(:,:), coeff(:,:)
        real(dp), parameter :: taps(3) = [1.0_dp, 4.0_dp, 1.0_dp] / 6.0_dp
        real(dp) :: lam
        integer :: code

        lam = 5.0_dp
        if (present(lmbda)) lam = lmbda
        coeff = cspline2d_real(signal, lamb=lam, status=code)
        if (code == spline_status_success) then
            filtered = sepfir2d(coeff, taps, taps, status=code)
        else
            allocate(filtered(0,0))
        end if
        if (present(status)) status = code
    end function spline_filter_real

    function spline_filter_complex(signal, lmbda, status) result(filtered)
        complex(dp), intent(in) :: signal(:,:) !! complex image to smooth
        real(dp), optional, intent(in) :: lmbda !! cubic smoothing parameter, default five
        integer, optional, intent(out) :: status !! zero for success
        complex(dp), allocatable :: filtered(:,:), coeff(:,:)
        real(dp), parameter :: taps(3) = [1.0_dp, 4.0_dp, 1.0_dp] / 6.0_dp
        real(dp) :: lam
        integer :: code

        lam = 5.0_dp
        if (present(lmbda)) lam = lmbda
        coeff = cspline2d_complex(signal, lamb=lam, status=code)
        if (code == spline_status_success) then
            filtered = sepfir2d(coeff, taps, taps, status=code)
        else
            allocate(filtered(0,0))
        end if
        if (present(status)) status = code
    end function spline_filter_complex

    function spline2d_core(signal, degree, lamb, precision, status) result(coeff)
        complex(dp), intent(in) :: signal(:,:) !! 2-D input image
        integer, intent(in) :: degree !! degree two or three
        real(dp), optional, intent(in) :: lamb !! optional smoothing strength
        real(dp), optional, intent(in) :: precision !! boundary-series tolerance
        integer, intent(out) :: status !! return status
        complex(dp), allocatable :: coeff(:,:), temp(:,:), line(:)
        real(dp) :: lam, eps
        integer :: i, j, code, nrow, ncol

        status = spline_status_invalid_input
        allocate(coeff(0,0))
        nrow = size(signal,1)
        ncol = size(signal,2)
        if (nrow == 0 .or. ncol == 0) return
        if (any(.not. ieee_is_finite(real(signal,dp))) .or. &
            any(.not. ieee_is_finite(aimag(signal)))) return
        if (degree /= 2 .and. degree /= 3) return
        lam = 0.0_dp
        if (present(lamb)) lam = lamb
        if (.not. ieee_is_finite(lam) .or. lam < 0.0_dp) return
        if (degree == 2 .and. lam /= 0.0_dp) return

        eps = 1.0e-6_dp
        if (present(precision)) then
            if (.not. ieee_is_finite(precision)) return
            if (precision > 0.0_dp .and. precision < 1.0_dp) eps = precision
        end if
        if (degree == 3 .and. lam > (1.0_dp/144.0_dp)) then
            if (nrow < 2 .or. ncol < 2) return
        end if
        allocate(temp(nrow,ncol))
        do i = 1, nrow
            line = spline_line(signal(i,:), degree, lam, eps, code)
            if (code /= spline_status_success) return
            temp(i,:) = line
        end do
        deallocate(coeff)
        allocate(coeff(nrow,ncol))
        do j = 1, ncol
            line = spline_line(temp(:,j), degree, lam, eps, code)
            if (code /= spline_status_success) then
                deallocate(coeff)
                allocate(coeff(0,0))
                return
            end if
            coeff(:,j) = line
        end do
        if (any(.not. ieee_is_finite(real(coeff,dp))) .or. &
            any(.not. ieee_is_finite(aimag(coeff)))) then
            deallocate(coeff)
            allocate(coeff(0,0))
            return
        end if
        status = spline_status_success
    end function spline2d_core

    function spline_line(signal, degree, lam, eps, status) result(coeff)
        complex(dp), intent(in) :: signal(:) !! one row or column
        integer, intent(in) :: degree !! two or three
        real(dp), intent(in) :: lam !! cubic smoothing strength
        real(dp), intent(in) :: eps !! truncation tolerance for the IIR tail
        integer, intent(out) :: status !! return status
        complex(dp), allocatable :: coeff(:)
        real(dp) :: r, scale, xi, omega

        status = spline_status_success
        if (size(signal) == 1) then
            ! A constant one-pixel axis is invariant under interpolation.
            coeff = signal
        else if (degree == 2) then
            r = -3.0_dp + 2.0_dp * sqrt(2.0_dp)
            scale = -8.0_dp*r
            coeff = symmetric_iir1(signal, r, scale, eps)
        else if (lam <= 1.0_dp / 144.0_dp) then
            r = -2.0_dp + sqrt(3.0_dp)
            scale = -6.0_dp*r
            coeff = symmetric_iir1(signal, r, scale, eps)
        else
            xi = 1.0_dp - 96.0_dp * lam + &
                24.0_dp * lam * sqrt(3.0_dp + 144.0_dp * lam)
            if (xi <= 0.0_dp) then
                status = spline_status_invalid_input
                allocate(coeff(0))
                return
            end if
            omega = atan(sqrt((144.0_dp*lam - 1.0_dp)/xi))
            r = ((24.0_dp*lam - 1.0_dp - sqrt(xi)) / (24.0_dp*lam)) * &
                sqrt(48.0_dp*lam + 24.0_dp*lam*sqrt(3.0_dp+144.0_dp*lam)) / sqrt(xi)
            if (.not. ieee_is_finite(r) .or. abs(r) >= 1.0_dp .or. &
                abs(sin(omega)) < tiny(omega)) then
                status = spline_status_invalid_input
                allocate(coeff(0))
                return
            end if
            coeff = symmetric_iir2(signal, r, omega, eps)
        end if
    end function spline_line

    function symmetric_iir1(signal, r, gain, eps) result(y)
        complex(dp), intent(in) :: signal(:) !! signal with symmetric endpoints
        real(dp), intent(in) :: r !! stable real IIR pole
        real(dp), intent(in) :: gain !! numerator of the symmetric filter
        real(dp), intent(in) :: eps !! boundary truncation tolerance
        complex(dp), allocatable :: y(:), causal(:)
        complex(dp) :: init
        real(dp) :: p
        integer :: n, i

        n = size(signal)
        allocate(causal(n), y(n))
        p = 1.0_dp
        init = signal(1)
        do i = 1, n
            init = init + r*p*signal(i)
            p = p*r
            if (abs(p) < eps) exit
        end do
        causal(1) = init
        do i = 2, n
            causal(i) = signal(i) + r*causal(i-1)
        end do
        y(n) = gain/(1.0_dp-r) * causal(n)
        do i = n-1, 1, -1
            y(i) = gain * causal(i) + r * y(i+1)
        end do
    end function symmetric_iir1

    function symmetric_iir2(signal, r, omega, eps) result(y)
        complex(dp), intent(in) :: signal(:) !! input signal
        real(dp), intent(in) :: r !! stable complex-pole radius
        real(dp), intent(in) :: omega !! pole angle in radians
        real(dp), intent(in) :: eps !! boundary truncation tolerance
        complex(dp), allocatable :: y(:), causal(:)
        complex(dp) :: initial0, initial1, back0, back1
        logical :: done0, done1
        real(dp) :: sn, cs, a2, a3, p, hscale, gamma, r2, pwr
        real(dp) :: h0, h1, h2, hm1
        integer :: n, i, k

        n = size(signal)
        allocate(causal(n),y(n))
        sn = sin(omega)
        r2 = r*r
        cs = 1.0_dp - 2.0_dp*r*cos(omega) + r2
        a2 = 2.0_dp*r*cos(omega)
        a3 = -r2

        initial0 = cs*signal(1)
        initial1 = cs*signal(2) + a2*cs*signal(1)
        p = r
        done0 = .false.
        done1 = .false.
        do k = 0, n-1
            if (.not. done0) then
                initial0 = initial0 + cs/sn*p*sin(real(k+2,dp)*omega)*signal(k+1)
                if (abs(p*sin(real(k+2,dp)*omega)) < eps) done0 = .true.
            end if
            if (.not. done1) then
                initial1 = initial1 + cs/sn*p*r*sin(real(k+3,dp)*omega)*signal(k+1)
                if (abs(p*r*sin(real(k+3,dp)*omega)) < eps) done1 = .true.
            end if
            p = p*r
            if (done0 .and. done1) exit
        end do
        causal(1) = initial0
        causal(2) = initial1
        do i = 3, n
            causal(i) = cs*signal(i) + a2*causal(i-1) + a3*causal(i-2)
        end do

        hscale = cs*cs * (1.0_dp+r2)/(1.0_dp-r2) / &
            (1.0_dp - 2.0_dp*r2*cos(2.0_dp*omega) + r2*r2)
        gamma = (1.0_dp-r2)/(1.0_dp+r2)/tan(omega)
        back0 = cmplx(0.0_dp, 0.0_dp, dp)
        back1 = cmplx(0.0_dp, 0.0_dp, dp)
        pwr = 1.0_dp
        done0 = .false.
        done1 = .false.
        do k = 0, n-1
            h0 = hscale*pwr*(cos(real(k,dp)*omega)+gamma*sin(real(k,dp)*omega))
            h1 = hscale*pwr*r*(cos(real(k+1,dp)*omega)+gamma*sin(real(k+1,dp)*omega))
            h2 = hscale*pwr*r*r*(cos(real(k+2,dp)*omega)+gamma*sin(real(k+2,dp)*omega))
            if (k == 0) then
                hm1 = h1
            else
                hm1 = hscale*pwr/r*(cos(real(k-1,dp)*omega)+gamma*sin(real(k-1,dp)*omega))
            end if
            if (.not. done0) then
                back0 = back0 + (h0+h1)*signal(n-k)
                if ((h0+h1)**2 < eps) done0 = .true.
            end if
            if (.not. done1) then
                back1 = back1 + (hm1+h2)*signal(n-k)
                if ((hm1+h2)**2 < eps) done1 = .true.
            end if
            pwr = pwr*r
            if (done0 .and. done1) exit
        end do
        y(n) = back0
        y(n-1) = back1
        do i = n-2, 1, -1
            y(i) = cs*causal(i) + a2*y(i+1) + a3*y(i+2)
        end do
    end function symmetric_iir2

end module scifort_signal_spline2d
