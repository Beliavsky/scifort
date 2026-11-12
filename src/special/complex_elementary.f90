! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Complex extensions of SciPy-style elementary special functions.  The
! implementation is independent and uses the defining analytic identities
! together with SciFort's Faddeeva/error-function kernels.

module scifort_complex_elementary
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi, scifort_sqrt_two
    use scifort_faddeeva, only : erfc_complex, erfcx_complex
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, &
        positive_infinity, quiet_nan
    use scifort_normal, only : normal_cdf, normal_logcdf
    implicit none
    private

    real(dp), parameter :: log_two = log(2.0_dp)

    public :: expm1_complex
    public :: log1p_complex
    public :: log_ndtr_complex
    public :: ndtr_complex
    public :: xlog1py_complex
    public :: xlogy_complex

contains

    pure elemental function ndtr_complex(z) result(value)
        complex(dp), intent(in) :: z !! complex standard-normal variate
        complex(dp) :: value

        real(dp) :: x
        real(dp) :: y

        x = real(z, dp)
        y = aimag(z)

        if (ieee_is_nan(x) .or. ieee_is_nan(y)) then
            value = complex_nan()
        else if (.not. ieee_is_finite(x)) then
            if (x > 0.0_dp .and. ieee_is_finite(y)) then
                value = cmplx(1.0_dp, 0.0_dp, kind=dp)
            else if (x < 0.0_dp .and. ieee_is_finite(y)) then
                value = cmplx(0.0_dp, 0.0_dp, kind=dp)
            else
                value = complex_nan()
            end if
        else if (.not. ieee_is_finite(y)) then
            if (x == 0.0_dp) then
                value = cmplx(0.5_dp, sign(positive_infinity(0.0_dp), y), kind=dp)
            else
                value = complex_nan()
            end if
        else if (y == 0.0_dp) then
            value = cmplx(normal_cdf(x), y, kind=dp)
        else
            value = 0.5_dp * erfc_complex(-z / scifort_sqrt_two)
        end if
    end function ndtr_complex

    pure elemental function log_ndtr_complex(z) result(value)
        complex(dp), intent(in) :: z !! complex standard-normal variate
        complex(dp) :: value

        complex(dp) :: complement
        complex(dp) :: raw
        complex(dp) :: scaled
        complex(dp) :: w
        real(dp) :: x
        real(dp) :: y

        x = real(z, dp)
        y = aimag(z)

        if (ieee_is_nan(x) .or. ieee_is_nan(y)) then
            value = complex_nan()
            return
        end if
        if (.not. ieee_is_finite(x)) then
            if (x > 0.0_dp .and. ieee_is_finite(y)) then
                value = cmplx(-0.0_dp, -0.0_dp, kind=dp)
            else if (x < 0.0_dp .and. ieee_is_finite(y)) then
                value = cmplx(negative_infinity(0.0_dp), quiet_nan(0.0_dp), kind=dp)
            else
                value = complex_nan()
            end if
            return
        end if
        if (.not. ieee_is_finite(y)) then
            value = complex_nan()
            return
        end if
        if (y == 0.0_dp) then
            value = cmplx(normal_logcdf(x), y, kind=dp)
            return
        end if

        ! In the right half-plane, reflection is accurate whenever the
        ! complementary CDF is small.  This avoids losing the tiny logarithm
        ! when Phi(z) is close to one.
        if (x >= 0.0_dp) then
            complement = ndtr_complex(-z)
            if (complex_is_finite(complement) .and. abs(complement) < 0.5_dp) then
                value = log1p_complex(-complement)
                return
            end if
        end if

        ! In the left half-plane use erfcx to avoid underflow in Phi(z):
        !   Phi(z) = exp(-z**2/2) * erfcx(-z/sqrt(2)) / 2.
        ! The imaginary part is reduced to the principal logarithm branch.
        if (x < -1.0_dp) then
            w = -z / scifort_sqrt_two
            scaled = erfcx_complex(w)
            raw = -0.5_dp * z * z + log(scaled) - log_two
            value = cmplx(real(raw, dp), principal_angle(aimag(raw)), kind=dp)
        else
            value = log(ndtr_complex(z))
        end if
    end function log_ndtr_complex

    pure elemental function expm1_complex(z) result(value)
        complex(dp), intent(in) :: z !! complex exponent in exp(z)-1
        complex(dp) :: value

        real(dp) :: cos_y
        real(dp) :: exp_x
        real(dp) :: expm1_x
        real(dp) :: real_part
        real(dp) :: sin_y
        real(dp) :: x
        real(dp) :: y

        x = real(z, dp)
        y = aimag(z)

        if (ieee_is_nan(x) .or. ieee_is_nan(y)) then
            value = complex_nan()
            return
        end if
        if (.not. ieee_is_finite(y)) then
            value = complex_nan()
            return
        end if
        if (y == 0.0_dp) then
            if (x > log(huge(1.0_dp))) then
                value = cmplx(positive_infinity(0.0_dp), y, kind=dp)
            else if (.not. ieee_is_finite(x)) then
                if (x > 0.0_dp) then
                    value = cmplx(positive_infinity(0.0_dp), y, kind=dp)
                else
                    value = cmplx(-1.0_dp, y, kind=dp)
                end if
            else
                value = cmplx(expm1_safe(x), y, kind=dp)
            end if
            return
        end if
        if (.not. ieee_is_finite(x)) then
            if (x < 0.0_dp) then
                value = cmplx(-1.0_dp, 0.0_dp, kind=dp)
            else
                cos_y = cos(y)
                sin_y = sin(y)
                value = cmplx(infinite_component(cos_y, -1.0_dp), &
                    infinite_component(sin_y, 0.0_dp), kind=dp)
            end if
            return
        end if

        cos_y = cos(y)
        sin_y = sin(y)
        if (x > log(huge(1.0_dp))) then
            value = cmplx(infinite_component(cos_y, -1.0_dp), &
                infinite_component(sin_y, 0.0_dp), kind=dp)
            return
        end if

        expm1_x = expm1_safe(x)
        real_part = expm1_x * cos_y - 2.0_dp * sin(0.5_dp * y) ** 2
        exp_x = exp(x)
        value = cmplx(real_part, exp_x * sin_y, kind=dp)
    end function expm1_complex

    pure elemental function log1p_complex(z) result(value)
        complex(dp), intent(in) :: z !! complex argument in the principal log(1+z)
        complex(dp) :: value

        real(dp) :: angle
        real(dp) :: magnitude
        real(dp) :: q
        real(dp) :: x
        real(dp) :: y

        x = real(z, dp)
        y = aimag(z)

        if (ieee_is_nan(x) .or. ieee_is_nan(y)) then
            value = complex_nan()
            return
        end if
        if (.not. ieee_is_finite(x) .or. .not. ieee_is_finite(y)) then
            if (.not. ieee_is_finite(x) .and. .not. ieee_is_finite(y)) then
                angle = atan2(sign(1.0_dp, y), sign(1.0_dp, x))
            else if (.not. ieee_is_finite(y)) then
                angle = sign(0.5_dp * scifort_pi, y)
            else if (x > 0.0_dp) then
                angle = sign(0.0_dp, y)
            else
                angle = sign(scifort_pi, y)
            end if
            value = cmplx(positive_infinity(0.0_dp), angle, kind=dp)
            return
        end if
        if (x == -1.0_dp .and. y == 0.0_dp) then
            value = cmplx(negative_infinity(0.0_dp), y, kind=dp)
            return
        end if

        angle = atan2(y, 1.0_dp + x)
        if (abs(z) <= 0.5_dp) then
            q = 2.0_dp * x + x * x + y * y
            value = cmplx(0.5_dp * log1p_safe(q), angle, kind=dp)
        else
            magnitude = hypot(1.0_dp + x, y)
            value = cmplx(log(magnitude), angle, kind=dp)
        end if
    end function log1p_complex

    pure elemental function xlogy_complex(x, y_arg) result(value)
        complex(dp), intent(in) :: x !! complex multiplier
        complex(dp), intent(in) :: y_arg !! complex logarithm argument
        complex(dp) :: value

        if (complex_has_nan(x) .or. complex_has_nan(y_arg)) then
            value = complex_nan()
        else if (x == cmplx(0.0_dp, 0.0_dp, kind=dp)) then
            value = cmplx(0.0_dp, 0.0_dp, kind=dp)
        else
            value = x * log(y_arg)
        end if
    end function xlogy_complex

    pure elemental function xlog1py_complex(x, y_arg) result(value)
        complex(dp), intent(in) :: x !! complex multiplier
        complex(dp), intent(in) :: y_arg !! complex argument in log(1+y_arg)
        complex(dp) :: value

        if (complex_has_nan(x) .or. complex_has_nan(y_arg)) then
            value = complex_nan()
        else if (x == cmplx(0.0_dp, 0.0_dp, kind=dp)) then
            value = cmplx(0.0_dp, 0.0_dp, kind=dp)
        else
            value = x * log1p_complex(y_arg)
        end if
    end function xlog1py_complex

    pure elemental logical function complex_has_nan(z) result(has_nan)
        complex(dp), intent(in) :: z !! complex value to classify

        has_nan = ieee_is_nan(real(z, dp)) .or. ieee_is_nan(aimag(z))
    end function complex_has_nan

    pure elemental logical function complex_is_finite(z) result(is_finite)
        complex(dp), intent(in) :: z !! complex value to classify

        is_finite = ieee_is_finite(real(z, dp)) .and. ieee_is_finite(aimag(z))
    end function complex_is_finite

    pure elemental function complex_nan() result(value)
        complex(dp) :: value

        value = cmplx(quiet_nan(0.0_dp), quiet_nan(0.0_dp), kind=dp)
    end function complex_nan

    pure elemental function infinite_component(factor, offset) result(value)
        real(dp), intent(in) :: factor !! finite sign factor multiplying positive infinity
        real(dp), intent(in) :: offset !! finite result used when factor is exactly zero
        real(dp) :: value

        if (factor > 0.0_dp) then
            value = positive_infinity(0.0_dp)
        else if (factor < 0.0_dp) then
            value = negative_infinity(0.0_dp)
        else
            value = offset
        end if
    end function infinite_component

    pure elemental function principal_angle(theta) result(angle)
        real(dp), intent(in) :: theta !! finite angle in radians
        real(dp) :: angle

        angle = atan2(sin(theta), cos(theta))
        if (angle <= -scifort_pi) angle = scifort_pi
    end function principal_angle

end module scifort_complex_elementary
