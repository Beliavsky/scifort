! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Complex E1 and Ei exponential integrals.
!
! This is an independent implementation based on NIST DLMF 6.2.1--6.2.7
! and the continued fraction in DLMF 6.9.1. Small and near-cut arguments
! use the convergent entire-function series for Ein; the complementary
! region uses the continued fraction. Exact real-axis inputs delegate to
! SciFort's existing real kernels so their established behavior is preserved.

module scifort_complex_exponential_integrals
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_exponential_integrals, only : exp1_real => exp1, expi_real => expi
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan
    implicit none
    private

    public :: exp1_complex
    public :: expi_complex

    real(dp), parameter :: euler_gamma = &
        0.577215664901532860606512090082402431_dp
    integer, parameter :: max_iterations = 10000

contains

    pure elemental function exp1_complex(z) result(value)
        complex(dp), intent(in) :: z !! complex argument of the principal E1 function
        complex(dp) :: value

        real(dp) :: x
        real(dp) :: y

        x = real(z, dp)
        y = aimag(z)

        if (ieee_is_nan(x) .or. ieee_is_nan(y) .or. &
                .not. ieee_is_finite(x) .or. .not. ieee_is_finite(y)) then
            value = complex_nan()
            return
        end if

        if (y == 0.0_dp) then
            if (x == 0.0_dp) then
                value = cmplx(positive_infinity(0.0_dp), 0.0_dp, kind=dp)
            else if (x > 0.0_dp) then
                value = cmplx(exp1_real(x), 0.0_dp, kind=dp)
            else
                value = cmplx(-expi_real(-x), -sign(scifort_pi, y), kind=dp)
            end if
            return
        end if

        if (use_e1_series(x, y)) then
            value = exp1_series(z)
        else
            value = exp1_continued_fraction(z)
        end if
    end function exp1_complex

    pure elemental function expi_complex(z) result(value)
        complex(dp), intent(in) :: z !! complex argument of the principal Ei function
        complex(dp) :: value

        real(dp) :: x
        real(dp) :: y

        x = real(z, dp)
        y = aimag(z)

        if (ieee_is_nan(x) .or. ieee_is_nan(y) .or. &
                .not. ieee_is_finite(x) .or. .not. ieee_is_finite(y)) then
            value = complex_nan()
            return
        end if

        if (y == 0.0_dp) then
            if (x > 0.0_dp) then
                value = cmplx(expi_real(x), 0.0_dp, kind=dp)
            else
                value = cmplx(expi_real(x), -0.0_dp, kind=dp)
            end if
            return
        end if

        if (use_ei_series(x, y)) then
            value = expi_series(z)
        else
            value = -exp1_complex(-z) + &
                cmplx(0.0_dp, sign(scifort_pi, y), kind=dp)
        end if
    end function expi_complex

    pure elemental logical function use_e1_series(x, y) result(use_series)
        real(dp), intent(in) :: x !! finite real component of the argument
        real(dp), intent(in) :: y !! finite nonzero imaginary component of the argument

        use_series = small_complex_argument(x, y)
        if (.not. use_series .and. x < 0.0_dp .and. -x < 40.0_dp) then
            use_series = abs(y) <= 0.35_dp * (-x)
        end if
    end function use_e1_series

    pure elemental logical function use_ei_series(x, y) result(use_series)
        real(dp), intent(in) :: x !! finite real component of the argument
        real(dp), intent(in) :: y !! finite nonzero imaginary component of the argument

        use_series = small_complex_argument(x, y)
        if (.not. use_series .and. x > 0.0_dp .and. x <= 700.0_dp) then
            use_series = abs(y) <= 0.10_dp * x
        end if
    end function use_ei_series

    pure elemental logical function small_complex_argument(x, y) result(is_small)
        real(dp), intent(in) :: x !! finite real component of the argument
        real(dp), intent(in) :: y !! finite imaginary component of the argument

        if (abs(x) > 8.0_dp .or. abs(y) > 8.0_dp) then
            is_small = .false.
        else
            is_small = x * x + y * y <= 64.0_dp
        end if
    end function small_complex_argument

    pure elemental function exp1_series(z) result(value)
        complex(dp), intent(in) :: z !! finite nonzero argument in the series region
        complex(dp) :: value

        complex(dp) :: sum_value
        complex(dp) :: term
        integer :: k

        term = z
        sum_value = term
        do k = 2, max_iterations
            term = term * (-z) * real(k - 1, dp) / real(k * k, dp)
            sum_value = sum_value + term
            if (abs(term) <= 2.0_dp * epsilon(1.0_dp) * &
                    max(1.0_dp, abs(sum_value))) exit
        end do

        value = sum_value - log(z) - euler_gamma
    end function exp1_series

    pure elemental function expi_series(z) result(value)
        complex(dp), intent(in) :: z !! finite nonzero argument in the Ei series region
        complex(dp) :: value

        complex(dp) :: sum_value
        complex(dp) :: term
        integer :: k

        term = z
        sum_value = term
        do k = 2, max_iterations
            term = term * z * real(k - 1, dp) / real(k * k, dp)
            sum_value = sum_value + term
            if (abs(term) <= 2.0_dp * epsilon(1.0_dp) * &
                    max(1.0_dp, abs(sum_value))) exit
        end do

        value = euler_gamma + log(z) + sum_value
    end function expi_series

    pure elemental function exp1_continued_fraction(z) result(value)
        complex(dp), intent(in) :: z !! finite nonzero argument away from the slow near-cut region
        complex(dp) :: value

        complex(dp) :: c
        complex(dp) :: d
        complex(dp) :: delta
        complex(dp) :: denominator
        complex(dp) :: factor
        complex(dp) :: b
        real(dp) :: fpmin
        integer :: k
        integer :: n

        fpmin = tiny(1.0_dp) / epsilon(1.0_dp)
        denominator = z
        if (abs(denominator) < fpmin) denominator = cmplx(fpmin, 0.0_dp, kind=dp)
        c = denominator
        d = cmplx(0.0_dp, 0.0_dp, kind=dp)

        do n = 1, max_iterations
            k = (n + 1) / 2
            if (mod(n, 2) == 1) then
                b = cmplx(1.0_dp, 0.0_dp, kind=dp)
            else
                b = z
            end if

            d = b + real(k, dp) * d
            if (abs(d) < fpmin) d = cmplx(fpmin, 0.0_dp, kind=dp)
            c = b + real(k, dp) / c
            if (abs(c) < fpmin) c = cmplx(fpmin, 0.0_dp, kind=dp)
            d = 1.0_dp / d
            delta = c * d
            denominator = denominator * delta
            if (abs(delta - 1.0_dp) <= 8.0_dp * epsilon(1.0_dp)) exit
        end do

        factor = 1.0_dp / denominator
        value = exp_minus_z_times(z, factor)
    end function exp1_continued_fraction

    pure elemental function exp_minus_z_times(z, factor) result(value)
        complex(dp), intent(in) :: z !! finite exponent argument
        complex(dp), intent(in) :: factor !! finite multiplicative factor
        complex(dp) :: value

        real(dp), parameter :: log_huge = log(huge(1.0_dp))
        real(dp) :: angle
        real(dp) :: cosine
        real(dp) :: magnitude
        real(dp) :: magnitude_log
        real(dp) :: sine

        if (real(factor, dp) == 0.0_dp .and. aimag(factor) == 0.0_dp) then
            value = cmplx(0.0_dp, 0.0_dp, kind=dp)
            return
        end if

        magnitude_log = -real(z, dp) + log(abs(factor))
        angle = -aimag(z) + atan2(aimag(factor), real(factor, dp))
        cosine = cos(angle)
        sine = sin(angle)

        if (magnitude_log > log_huge) then
            value = cmplx(infinite_component(cosine), &
                infinite_component(sine), kind=dp)
        else
            magnitude = exp(magnitude_log)
            value = cmplx(magnitude * cosine, magnitude * sine, kind=dp)
        end if
    end function exp_minus_z_times

    pure elemental function infinite_component(factor) result(value)
        real(dp), intent(in) :: factor !! finite trigonometric factor multiplying infinity
        real(dp) :: value

        if (factor == 0.0_dp) then
            value = factor
        else
            value = sign(positive_infinity(0.0_dp), factor)
        end if
    end function infinite_component

    pure function complex_nan() result(value)
        complex(dp) :: value

        real(dp) :: nan_value

        nan_value = quiet_nan(0.0_dp)
        value = cmplx(nan_value, nan_value, kind=dp)
    end function complex_nan

end module scifort_complex_exponential_integrals
