! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Analog Butterworth, Chebyshev I, and Chebyshev II filter design.
! Independent prototype and low/high/band transformations, not a
! translation of external implementation code. See CODE_PROVENANCE.md.
module scifort_signal_analog_design
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter, public :: analog_design_success = 0
    integer, parameter, public :: analog_design_invalid_input = 1

    type, public :: analog_tf_result
        real(dp), allocatable :: b(:)
        real(dp), allocatable :: a(:)
        complex(dp), allocatable :: zeros(:) !! finite zeros, excluding roots at infinity
        complex(dp), allocatable :: poles(:) !! finite analog poles
        real(dp) :: gain = 0.0_dp !! ZPK gain, for descending powers of s
        integer :: status = analog_design_invalid_input
    end type analog_tf_result

    public :: analog_iir

contains

    ! Construct an analog transfer function H(s) from normalized prototypes.
    function analog_iir(order, cutoff, btype, family, attenuation, roots_only) result(out)
        integer, intent(in) :: order !! analog lowpass prototype order (1 through 20)
        real(dp), intent(in) :: cutoff(:) !! rad/s critical frequencies
        character(*), intent(in), optional :: btype !! lowpass, highpass, bandpass, bandstop
        character(*), intent(in) :: family !! butter, cheby1, or cheby2
        real(dp), intent(in), optional :: attenuation !! Chebyshev pass/stop dB
        logical, intent(in), optional :: roots_only !! skip potentially overflowing BA expansion
        type(analog_tf_result) :: out
        complex(dp), allocatable :: p(:), z(:), transformed_p(:), transformed_z(:)
        complex(dp) :: pole, v, disc, gain
        real(dp) :: epsilon, mu, theta, m, omega, bw, center
        real(dp), allocatable :: a(:), b(:)
        character(:), allocatable :: mode
        integer :: i, k, nz, degree, nfinal

        if (order < 1 .or. order > 20) return
        if (size(cutoff) < 1 .or. size(cutoff) > 2) return
        if (.not. all(ieee_is_finite(cutoff))) return
        if (any(cutoff <= 0.0_dp)) return
        mode = 'lowpass'
        if (present(btype)) mode = trim(btype)
        select case (mode)
        case ('low', 'lowpass', 'high', 'highpass')
            if (size(cutoff) /= 1) return
        case ('band', 'bandpass', 'stop', 'bandstop')
            if (size(cutoff) /= 2) return
            if (cutoff(1) >= cutoff(2)) return
        case default
            return
        end select
        omega = cutoff(1)
        bw = 0.0_dp
        center = 0.0_dp
        if (size(cutoff) == 2) then
            bw = cutoff(2) - cutoff(1)
            center = sqrt(cutoff(1)) * sqrt(cutoff(2))
            if (.not. ieee_is_finite(bw) .or. .not. ieee_is_finite(center)) return
        end if
        select case (family)
        case ('butter', 'butterworth')
            if (present(attenuation)) return
        case ('cheby1', 'chebyshev1', 'cheby2', 'chebyshev2')
            if (.not. present(attenuation)) return
            if (.not. ieee_is_finite(attenuation)) return
            if (attenuation <= 0.0_dp .or. attenuation > 100.0_dp) return
        case default
            return
        end select

        allocate(p(order), z(order))
        nz = 0
        gain = cmplx(1.0_dp, 0.0_dp, dp)
        select case (family)
        case ('butter', 'butterworth')
            do k = 0, order-1
                theta = pi * real(2*k+order+1, dp) / real(2*order, dp)
                p(k+1) = cmplx(cos(theta), sin(theta), dp)
            end do
        case ('cheby1', 'chebyshev1')
            epsilon = sqrt(10.0_dp**(attenuation/10.0_dp) - 1.0_dp)
            mu = asinh(1.0_dp/epsilon) / real(order, dp)
            do k = 0, order-1
                theta = pi * real(2*k+1, dp) / real(2*order, dp)
                p(k+1) = cmplx(-sinh(mu)*sin(theta), cosh(mu)*cos(theta), dp)
            end do
        case ('cheby2', 'chebyshev2')
            epsilon = 1.0_dp / sqrt(10.0_dp**(attenuation/10.0_dp) - 1.0_dp)
            mu = asinh(1.0_dp/epsilon) / real(order, dp)
            do k = 0, order-1
                m = real(2*k-order+1, dp)
                theta = pi*m / real(2*order, dp)
                pole = cmplx(-sinh(mu)*cos(theta), -cosh(mu)*sin(theta), dp)
                p(k+1) = 1.0_dp / pole
                if (abs(m) > 0.5_dp) then
                    nz = nz + 1
                    z(nz) = cmplx(0.0_dp, 1.0_dp/sin(theta), dp)
                end if
            end do
        end select

        ! Normalize the lowpass prototype to unit DC gain. Even Chebyshev I
        ! prototypes instead have the prescribed -rp dB gain at DC.
        do i = 1, order
            gain = gain * (-p(i))
        end do
        do i = 1, nz
            gain = gain / (-z(i))
        end do
        if (family == 'cheby1' .or. family == 'chebyshev1') then
            if (mod(order, 2) == 0) gain = gain / sqrt(1.0_dp + epsilon**2)
        end if
        degree = order - nz
        nfinal = order
        if (size(cutoff) == 2) nfinal = 2*order
        allocate(transformed_p(nfinal))
        select case (mode)
        case ('low', 'lowpass')
            transformed_p = omega*p
            allocate(transformed_z(nz))
            transformed_z = omega*z(:nz)
            gain = gain * omega**degree
        case ('high', 'highpass')
            transformed_p = omega/p
            allocate(transformed_z(order))
            transformed_z(:nz) = omega/z(:nz)
            transformed_z(nz+1:) = cmplx(0.0_dp, 0.0_dp, dp)
            gain = gain * gain_at_infinity(z(:nz), p)
        case ('band', 'bandpass')
            do i = 1, order
                v = (bw/2.0_dp)*p(i)
                disc = sqrt(v*v - cmplx(center**2, 0.0_dp, dp))
                transformed_p(2*i-1) = v + disc
                transformed_p(2*i) = v - disc
            end do
            allocate(transformed_z(2*nz+degree))
            do i = 1, nz
                v = (bw/2.0_dp)*z(i)
                disc = sqrt(v*v - cmplx(center**2, 0.0_dp, dp))
                transformed_z(2*i-1) = v + disc
                transformed_z(2*i) = v - disc
            end do
            transformed_z(2*nz+1:) = cmplx(0.0_dp, 0.0_dp, dp)
            gain = gain * bw**degree
        case ('stop', 'bandstop')
            do i = 1, order
                v = (bw/2.0_dp)/p(i)
                disc = sqrt(v*v - cmplx(center**2, 0.0_dp, dp))
                transformed_p(2*i-1) = v + disc
                transformed_p(2*i) = v - disc
            end do
            allocate(transformed_z(2*nz+2*degree))
            do i = 1, nz
                v = (bw/2.0_dp)/z(i)
                disc = sqrt(v*v - cmplx(center**2, 0.0_dp, dp))
                transformed_z(2*i-1) = v + disc
                transformed_z(2*i) = v - disc
            end do
            do i = 1, degree
                transformed_z(2*nz+2*i-1) = cmplx(0.0_dp, center, dp)
                transformed_z(2*nz+2*i) = cmplx(0.0_dp, -center, dp)
            end do
            gain = gain * gain_at_infinity(z(:nz), p)
        end select

        if (.not. valid_complex(gain)) return
        if (.not. all(ieee_is_finite(real(transformed_p, dp)))) return
        if (.not. all(ieee_is_finite(aimag(transformed_p)))) return
        if (.not. all(ieee_is_finite(real(transformed_z, dp)))) return
        if (.not. all(ieee_is_finite(aimag(transformed_z)))) return
        out%zeros = transformed_z
        out%poles = transformed_p
        out%gain = real(gain, dp)
        if (present(roots_only)) then
            if (roots_only) then
                out%status = analog_design_success
                return
            end if
        end if
        call polynomial(transformed_p, a)
        call polynomial(transformed_z, b)
        b = b * real(gain, dp)
        if (.not. all(ieee_is_finite(a))) return
        if (.not. all(ieee_is_finite(b))) return
        out%a = a
        out%b = b
        out%status = analog_design_success
    end function analog_iir

    pure function gain_at_infinity(z, p) result(r)
        complex(dp), intent(in) :: z(:) !! finite lowpass prototype zeros
        complex(dp), intent(in) :: p(:) !! lowpass prototype poles
        complex(dp) :: r
        integer :: i
        r = cmplx(1.0_dp, 0.0_dp, dp)
        do i = 1, size(z)
            r = r*(-z(i))
        end do
        do i = 1, size(p)
            r = r/(-p(i))
        end do
    end function gain_at_infinity

    pure subroutine polynomial(roots, coeff)
        complex(dp), intent(in) :: roots(:) !! complex roots of real polynomial
        real(dp), allocatable, intent(out) :: coeff(:) !! descending real coefficients
        complex(dp) :: c(size(roots)+1)
        integer :: i, j
        c = cmplx(0.0_dp, 0.0_dp, dp)
        c(1) = cmplx(1.0_dp, 0.0_dp, dp)
        do i = 1, size(roots)
            do j = i+1, 2, -1
                c(j) = c(j) - roots(i)*c(j-1)
            end do
        end do
        coeff = real(c, dp)
    end subroutine polynomial

    pure logical function valid_complex(z) result(ok)
        complex(dp), intent(in) :: z !! candidate complex gain
        ok = ieee_is_finite(real(z, dp)) .and. ieee_is_finite(aimag(z))
    end function valid_complex

end module scifort_signal_analog_design
