! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Complex digamma/psi support.
!
! This is an independent implementation based on the NIST DLMF recurrence,
! reflection, and asymptotic expansion for the psi function. Exact real-axis
! nonpole inputs reuse SciFort's real digamma kernel so existing real behavior
! and accuracy are preserved.

module scifort_complex_digamma
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : positive_infinity, quiet_nan
    use scifort_special_elementary, only : digamma_real => digamma
    implicit none
    private

    public :: digamma_complex
    public :: psi_complex

    real(dp), parameter :: asymptotic_threshold = 8.0_dp

contains

    pure elemental function digamma_complex(z) result(value)
        complex(dp), intent(in) :: z !! complex argument of the digamma function
        complex(dp) :: value

        complex(dp) :: reflected
        real(dp) :: x
        real(dp) :: y

        x = real(z, dp)
        y = aimag(z)

        if (ieee_is_nan(x) .or. ieee_is_nan(y)) then
            value = complex_nan()
            return
        end if

        if (.not. ieee_is_finite(x) .or. .not. ieee_is_finite(y)) then
            value = digamma_nonfinite(x, y)
            return
        end if

        if (y == 0.0_dp) then
            if (x <= 0.0_dp .and. x == anint(x)) then
                value = complex_nan()
            else
                value = cmplx(digamma_real(x), sign(0.0_dp, y), kind=dp)
            end if
            return
        end if

        if (x < 0.5_dp) then
            reflected = digamma_right_half(cmplx(1.0_dp - x, -y, kind=dp))
            value = reflected - scifort_pi * cotpi_complex(z)
        else
            value = digamma_right_half(z)
        end if
    end function digamma_complex

    pure elemental function psi_complex(z) result(value)
        complex(dp), intent(in) :: z !! complex argument; alias of digamma
        complex(dp) :: value

        value = digamma_complex(z)
    end function psi_complex

    pure function digamma_right_half(z) result(value)
        complex(dp), intent(in) :: z !! finite complex argument with real part at least 0.5
        complex(dp) :: value

        complex(dp) :: inverse
        complex(dp) :: inverse2
        complex(dp) :: shifted

        shifted = z
        value = cmplx(0.0_dp, 0.0_dp, kind=dp)

        do while (abs(shifted) < asymptotic_threshold)
            value = value - 1.0_dp / shifted
            shifted = shifted + 1.0_dp
        end do

        inverse = 1.0_dp / shifted
        inverse2 = inverse * inverse
        value = value + log(shifted) - 0.5_dp * inverse + inverse2 * ( &
            -1.0_dp / 12.0_dp + inverse2 * ( &
            1.0_dp / 120.0_dp + inverse2 * ( &
            -1.0_dp / 252.0_dp + inverse2 * ( &
            1.0_dp / 240.0_dp + inverse2 * ( &
            -1.0_dp / 132.0_dp + inverse2 * ( &
            691.0_dp / 32760.0_dp - inverse2 / 12.0_dp))))))
    end function digamma_right_half

    pure elemental function cotpi_complex(z) result(value)
        complex(dp), intent(in) :: z !! finite nonreal argument of cot(pi*z)
        complex(dp) :: value

        real(dp) :: angle
        real(dp) :: c
        real(dp) :: ch
        real(dp) :: denominator
        real(dp) :: q
        real(dp) :: r
        real(dp) :: s
        real(dp) :: scale
        real(dp) :: sh
        real(dp) :: x
        real(dp) :: y

        x = real(z, dp)
        y = aimag(z)
        r = modulo(x, 1.0_dp)
        if (r > 0.5_dp) r = r - 1.0_dp
        angle = scifort_pi * y

        if (abs(angle) <= 20.0_dp) then
            s = sin(scifort_pi * r)
            c = cos(scifort_pi * r)
            sh = sinh(angle)
            ch = cosh(angle)
            scale = max(abs(s), abs(sh))
            denominator = (s / scale)**2 + (sh / scale)**2
            denominator = scale * denominator
            value = cmplx( &
                safe_real_ratio((s / scale) * c, denominator), &
                safe_real_ratio(-(sh / scale) * ch, denominator), kind=dp)
        else
            q = exp(-2.0_dp * abs(angle))
            denominator = 1.0_dp + q * q - 2.0_dp * q * cos(2.0_dp * scifort_pi * r)
            value = cmplx( &
                2.0_dp * q * sin(2.0_dp * scifort_pi * r) / denominator, &
                -sign(1.0_dp, y) * (1.0_dp - q * q) / denominator, kind=dp)
        end if
    end function cotpi_complex

    pure elemental function digamma_nonfinite(x, y) result(value)
        real(dp), intent(in) :: x !! real component, possibly infinite
        real(dp), intent(in) :: y !! imaginary component, possibly infinite
        complex(dp) :: value

        real(dp) :: angle
        real(dp) :: infinity

        infinity = positive_infinity(0.0_dp)

        if (.not. ieee_is_finite(y)) then
            if (.not. ieee_is_finite(x)) then
                if (x > 0.0_dp) then
                    angle = sign(0.25_dp * scifort_pi, y)
                else
                    angle = sign(0.75_dp * scifort_pi, y)
                end if
            else
                angle = sign(0.5_dp * scifort_pi, y)
            end if
            value = cmplx(infinity, angle, kind=dp)
        else if (x > 0.0_dp) then
            value = cmplx(infinity, sign(0.0_dp, y), kind=dp)
        else
            value = complex_nan()
        end if
    end function digamma_nonfinite

    pure function safe_real_ratio(numerator, denominator) result(value)
        real(dp), intent(in) :: numerator !! finite numerator
        real(dp), intent(in) :: denominator !! positive finite denominator
        real(dp) :: value

        real(dp) :: infinity

        infinity = positive_infinity(0.0_dp)
        if (denominator == 0.0_dp) then
            value = sign(infinity, numerator)
        else if (denominator < 1.0_dp) then
            if (abs(numerator) > huge(1.0_dp) * denominator) then
                value = sign(infinity, numerator)
            else
                value = numerator / denominator
            end if
        else
            value = numerator / denominator
        end if
    end function safe_real_ratio

    pure function complex_nan() result(value)
        complex(dp) :: value

        real(dp) :: nan_value

        nan_value = quiet_nan(0.0_dp)
        value = cmplx(nan_value, nan_value, kind=dp)
    end function complex_nan

end module scifort_complex_digamma
