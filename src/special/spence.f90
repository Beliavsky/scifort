! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Real and complex Spence function scipy.special.spence(z) = Li_2(1-z).
! The implementation is independent of SciPy source code and uses the
! dilogarithm series and connection formulas DLMF 25.12.1, 25.12.4,
! 25.12.6, and the logarithmic expansion obtained from DLMF 25.12.12.

module scifort_spence
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, quiet_nan
    implicit none
    private

    public :: spence

    real(dp), parameter :: pi_squared_over_six = scifort_pi * scifort_pi / 6.0_dp
    real(dp), parameter :: series_radius = 0.85_dp
    integer, parameter :: max_series_iterations = 2048

    ! Coefficients -B_(2k)/[(2k)(2k+1)(2k)!], k=1,...,24, for the
    ! s -> 2 logarithmic expansion of Li_s(exp(mu)); see DLMF 25.12.12.
    real(dp), parameter :: log_series_coefficients(24) = [ &
        -1.388888888888888812e-02_dp, 6.944444444444444438e-05_dp, &
        -7.873519778281682970e-07_dp, 1.148221634332745510e-08_dp, &
        -1.897886998897099899e-10_dp, 3.387301370953521202e-12_dp, &
        -6.372636443183180759e-14_dp, 1.246205991295067152e-15_dp, &
        -2.510544460899954549e-17_dp, 5.178258806090623204e-19_dp, &
        -1.088735736830084924e-20_dp, 2.325744114302087076e-22_dp, &
        -5.035195213147389651e-24_dp, 1.102649929438121498e-25_dp, &
        -2.438658550900734404e-27_dp, 5.440142678856252739e-29_dp, &
        -1.222834013121735176e-30_dp, 2.767263468967950833e-32_dp, &
        -6.300090591832013554e-34_dp, 1.442086838841847638e-35_dp, &
        -3.317093999159542757e-37_dp, 7.663913557920658378e-39_dp, &
        -1.777871473383065865e-40_dp, 4.139605898234137515e-42_dp]

    interface spence
        module procedure spence_real
        module procedure spence_complex
    end interface spence

contains

    pure elemental function spence_real(x) result(value)
        real(dp), intent(in) :: x !! real argument; supported on the nonnegative real axis
        real(dp) :: value
        real(dp) :: z

        if (ieee_is_nan(x)) then
            value = x
            return
        end if
        if (.not. ieee_is_finite(x)) then
            value = quiet_nan(x)
            return
        end if
        if (x < 0.0_dp) then
            value = quiet_nan(x)
            return
        end if
        if (x == 0.0_dp) then
            value = pi_squared_over_six
            return
        end if

        z = 1.0_dp - x
        value = dilogarithm_real(z)
    end function spence_real

    pure elemental function spence_complex(z) result(value)
        complex(dp), intent(in) :: z !! complex argument on the principal Spence branch
        complex(dp) :: value
        complex(dp) :: dilog_argument
        real(dp) :: x
        real(dp) :: y

        x = real(z, dp)
        y = aimag(z)

        if (ieee_is_nan(x) .or. ieee_is_nan(y) .or. &
                .not. ieee_is_finite(x) .or. .not. ieee_is_finite(y)) then
            value = complex_nan()
            return
        end if

        ! Preserve the existing real branch exactly for nonnegative real
        ! complex inputs. Negative real complex inputs instead take the
        ! principal cut limit, with signed zero selecting the side of the cut.
        if (y == 0.0_dp .and. x >= 0.0_dp) then
            if (x > 1.0_dp) then
                value = cmplx(spence_real(x), -0.0_dp, kind=dp)
            else
                value = cmplx(spence_real(x), 0.0_dp, kind=dp)
            end if
            return
        end if

        ! Li_2(1-z). Construct the imaginary part explicitly so signed zero is
        ! retained on the negative-real Spence branch cut.
        dilog_argument = cmplx(1.0_dp - x, -y, kind=dp)
        value = dilogarithm_complex(dilog_argument)
    end function spence_complex

    pure elemental function dilogarithm_real(z) result(value)
        real(dp), intent(in) :: z !! real dilogarithm argument, <= 1
        real(dp) :: value
        real(dp) :: transformed
        real(dp) :: log_term

        if (z == 1.0_dp) then
            value = pi_squared_over_six
        else if (z == 0.0_dp) then
            value = 0.0_dp
        else if (z < -1.0_dp) then
            transformed = 1.0_dp / z
            log_term = log(-z)
            value = -pi_squared_over_six - 0.5_dp * log_term * log_term - &
                dilogarithm_reduced(transformed)
        else
            value = dilogarithm_reduced(z)
        end if
    end function dilogarithm_real

    pure elemental function dilogarithm_reduced(z) result(value)
        real(dp), intent(in) :: z !! real dilogarithm argument in [-1, 1)
        real(dp) :: value
        real(dp) :: transformed
        real(dp) :: log_term

        if (z < -0.5_dp) then
            transformed = z / (z - 1.0_dp)
            log_term = log1p_safe(-z)
            value = -0.5_dp * log_term * log_term - dilogarithm_series_real(transformed)
        else if (z <= 0.5_dp) then
            value = dilogarithm_series_real(z)
        else
            transformed = 1.0_dp - z
            value = pi_squared_over_six - log(z) * log1p_safe(-z) - &
                dilogarithm_series_real(transformed)
        end if
    end function dilogarithm_reduced

    pure elemental function dilogarithm_series_real(z) result(value)
        real(dp), intent(in) :: z !! series argument with abs(z) <= 0.5
        real(dp) :: value
        real(dp) :: power
        real(dp) :: addend
        integer :: k

        if (z == 0.0_dp) then
            value = 0.0_dp
            return
        end if

        power = z
        value = z
        do k = 2, 256
            power = power * z
            addend = power / real(k * k, dp)
            value = value + addend
            if (abs(addend) <= epsilon(1.0_dp) * max(1.0_dp, abs(value))) exit
        end do
    end function dilogarithm_series_real

    pure elemental function dilogarithm_complex(z) result(value)
        complex(dp), intent(in) :: z !! finite principal-branch dilogarithm argument
        complex(dp) :: value
        complex(dp) :: inverse_z
        real(dp) :: radius

        if (real(z, dp) == 0.0_dp .and. aimag(z) == 0.0_dp) then
            value = cmplx(0.0_dp, aimag(z), kind=dp)
            return
        end if
        if (real(z, dp) == 1.0_dp .and. aimag(z) == 0.0_dp) then
            value = cmplx(pi_squared_over_six, aimag(z), kind=dp)
            return
        end if

        radius = bounded_complex_magnitude(z)
        if (radius <= series_radius) then
            value = dilogarithm_series_complex(z)
        else if (radius >= 1.0_dp / series_radius) then
            inverse_z = reciprocal_complex(z)
            value = cmplx(-pi_squared_over_six, 0.0_dp, kind=dp) - &
                0.5_dp * log(-z) * log(-z) - dilogarithm_series_complex(inverse_z)
        else
            value = dilogarithm_log_series(z)
        end if
    end function dilogarithm_complex

    pure elemental function dilogarithm_series_complex(z) result(value)
        complex(dp), intent(in) :: z !! complex series argument with abs(z) <= 0.85
        complex(dp) :: value
        complex(dp) :: power
        complex(dp) :: addend
        integer :: k

        if (real(z, dp) == 0.0_dp .and. aimag(z) == 0.0_dp) then
            value = z
            return
        end if

        power = z
        value = z
        do k = 2, max_series_iterations
            power = power * z
            addend = power / real(k * k, dp)
            value = value + addend
            if (abs(addend) <= epsilon(1.0_dp) * max(1.0_dp, abs(value))) exit
        end do
    end function dilogarithm_series_complex

    pure elemental function dilogarithm_log_series(z) result(value)
        complex(dp), intent(in) :: z !! complex argument in the annulus around abs(z)=1
        complex(dp) :: value
        complex(dp) :: mu
        complex(dp) :: mu_squared
        complex(dp) :: power
        integer :: k

        mu = log(z)
        value = cmplx(pi_squared_over_six, 0.0_dp, kind=dp) + &
            mu * (1.0_dp - log(-mu)) - 0.25_dp * mu * mu
        mu_squared = mu * mu
        power = mu * mu_squared

        do k = 1, size(log_series_coefficients)
            value = value + log_series_coefficients(k) * power
            power = power * mu_squared
        end do
    end function dilogarithm_log_series

    pure elemental function bounded_complex_magnitude(z) result(radius)
        complex(dp), intent(in) :: z !! finite complex value used for numerical region selection
        real(dp) :: radius
        real(dp) :: scale
        real(dp) :: x
        real(dp) :: y

        x = abs(real(z, dp))
        y = abs(aimag(z))
        scale = max(x, y)
        if (scale == 0.0_dp) then
            radius = 0.0_dp
        else if (scale > 2.0_dp) then
            radius = scale
        else
            radius = scale * sqrt((x / scale) ** 2 + (y / scale) ** 2)
        end if
    end function bounded_complex_magnitude

    pure elemental function reciprocal_complex(z) result(value)
        complex(dp), intent(in) :: z !! nonzero finite complex argument
        complex(dp) :: value
        real(dp) :: ax
        real(dp) :: ay
        real(dp) :: ratio
        real(dp) :: denominator
        real(dp) :: x
        real(dp) :: y

        x = real(z, dp)
        y = aimag(z)
        ax = abs(x)
        ay = abs(y)

        if (ax >= ay) then
            ratio = y / x
            denominator = 1.0_dp + ratio * ratio
            value = cmplx((1.0_dp / x) / denominator, &
                -(ratio / x) / denominator, kind=dp)
        else
            ratio = x / y
            denominator = 1.0_dp + ratio * ratio
            value = cmplx((ratio / y) / denominator, &
                -(1.0_dp / y) / denominator, kind=dp)
        end if
    end function reciprocal_complex

    pure function complex_nan() result(value)
        complex(dp) :: value

        value = cmplx(quiet_nan(0.0_dp), quiet_nan(0.0_dp), kind=dp)
    end function complex_nan

end module scifort_spence
