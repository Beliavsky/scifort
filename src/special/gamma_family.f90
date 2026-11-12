! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Real and complex gamma-family helpers matching SciPy behavior where practical.
!
! The implementation is independent of SciPy source code. Real gamma magnitudes
! reuse SciFort's gammaln kernel; complex log-gamma uses DLMF recurrence,
! reflection, and Stirling expansion with explicit principal-branch tracking.

module scifort_gamma_family
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_log_sqrt_two_pi, scifort_log_two, scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, positive_infinity, quiet_nan
    use scifort_special_elementary, only : gammaln
    implicit none
    private

    public :: gamma
    public :: gammasgn
    public :: loggamma
    public :: multigammaln
    public :: rgamma

    interface gamma
        module procedure gamma_real
        module procedure gamma_complex
    end interface gamma

    interface loggamma
        module procedure loggamma_real
        module procedure loggamma_complex
    end interface loggamma

    interface rgamma
        module procedure rgamma_real
        module procedure rgamma_complex
    end interface rgamma

    ! B(2k) / (2k (2k - 1)) for k = 1, ..., 7 (DLMF 5.11.1).
    real(dp), parameter :: complex_stirling_coefficients(7) = [ &
        1.0_dp / 12.0_dp, &
        -1.0_dp / 360.0_dp, &
        1.0_dp / 1260.0_dp, &
        -1.0_dp / 1680.0_dp, &
        1.0_dp / 1188.0_dp, &
        -691.0_dp / 360360.0_dp, &
        1.0_dp / 156.0_dp]
    real(dp), parameter :: complex_stirling_threshold = 10.0_dp

contains

    pure elemental function gammasgn(x) result(y)
        real(dp), intent(in) :: x !! real gamma argument
        real(dp) :: y

        real(dp) :: interval_index

        if (ieee_is_nan(x)) then
            y = x
        else if (x > 0.0_dp) then
            y = 1.0_dp
        else if (x == 0.0_dp) then
            y = sign(1.0_dp, x)
        else if (.not. ieee_is_finite(x)) then
            y = quiet_nan(x)
        else if (x == anint(x)) then
            y = quiet_nan(x)
        else
            ! For x in (-(k+1),-k), sign(Gamma(x)) = (-1)**(k+1).
            ! Every finite binary64 negative noninteger has magnitude below
            ! 2**52, so the integer interval index is represented exactly.
            interval_index = aint(-x)
            if (modulo(interval_index, 2.0_dp) == 0.0_dp) then
                y = -1.0_dp
            else
                y = 1.0_dp
            end if
        end if
    end function gammasgn

    pure elemental function gamma_real(x) result(y)
        real(dp), intent(in) :: x !! real gamma argument
        real(dp) :: y

        real(dp) :: log_abs
        real(dp) :: sign_value

        if (ieee_is_nan(x)) then
            y = x
            return
        end if
        if (x == 0.0_dp) then
            y = sign(positive_infinity(x), x)
            return
        end if
        if (.not. ieee_is_finite(x)) then
            if (x > 0.0_dp) then
                y = positive_infinity(x)
            else
                y = quiet_nan(x)
            end if
            return
        end if
        if (x < 0.0_dp .and. x == anint(x)) then
            y = quiet_nan(x)
            return
        end if

        sign_value = gammasgn(x)
        log_abs = gammaln(x)
        y = signed_exp(log_abs, sign_value, x)
    end function gamma_real

    pure elemental function rgamma_real(x) result(y)
        real(dp), intent(in) :: x !! real argument of the reciprocal gamma function
        real(dp) :: y

        real(dp) :: log_abs
        real(dp) :: sign_value

        if (ieee_is_nan(x)) then
            y = x
            return
        end if
        if (x == 0.0_dp) then
            y = sign(0.0_dp, x)
            return
        end if
        if (.not. ieee_is_finite(x)) then
            y = 0.0_dp
            return
        end if
        if (x < 0.0_dp .and. x == anint(x)) then
            y = 0.0_dp
            return
        end if

        sign_value = gammasgn(x)
        log_abs = -gammaln(x)
        y = signed_exp(log_abs, sign_value, x)
    end function rgamma_real

    pure elemental function loggamma_real(x) result(y)
        real(dp), intent(in) :: x !! real argument; negative values use SciPy's NaN real branch
        real(dp) :: y

        if (ieee_is_nan(x)) then
            y = x
        else if (x == 0.0_dp) then
            y = positive_infinity(x)
        else if (x > 0.0_dp) then
            y = gammaln(x)
        else
            y = quiet_nan(x)
        end if
    end function loggamma_real

    pure elemental function loggamma_complex(z) result(value)
        complex(dp), intent(in) :: z !! complex argument on the principal log-gamma branch
        complex(dp) :: value

        complex(dp) :: reflected
        real(dp) :: branch_count
        real(dp) :: t
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
                value = cmplx(gammaln(x), y, kind=dp)
            else if (x == anint(x)) then
                value = complex_nan()
            else
                branch_count = aint(-x) + 1.0_dp
                value = cmplx(gammaln(x), &
                    -sign(1.0_dp, y) * scifort_pi * branch_count, kind=dp)
            end if
            return
        end if

        if (x >= 0.5_dp) then
            value = loggamma_right_half_plane(z)
            return
        end if

        reflected = cmplx(log(scifort_pi), 0.0_dp, kind=dp) - &
            log_sin_pi(z) - loggamma_right_half_plane(1.0_dp - z)

        ! The principal log used by log(sin(pi*z)) changes sheet in vertical
        ! strips. Correct it so continuation from the positive real axis has
        ! the single negative-real branch cut documented for log Gamma.
        t = max(0.0_dp, -0.5_dp * (x + 0.5_dp))
        branch_count = aint(t)
        if (branch_count < t) branch_count = branch_count + 1.0_dp
        if (y > 0.0_dp) then
            value = reflected - cmplx(0.0_dp, &
                2.0_dp * scifort_pi * branch_count, kind=dp)
        else
            value = reflected + cmplx(0.0_dp, &
                2.0_dp * scifort_pi * branch_count, kind=dp)
        end if
    end function loggamma_complex

    pure elemental function gamma_complex(z) result(value)
        complex(dp), intent(in) :: z !! complex argument of the gamma function
        complex(dp) :: value

        real(dp) :: x
        real(dp) :: y

        x = real(z, dp)
        y = aimag(z)
        if (ieee_is_nan(x) .or. ieee_is_nan(y) .or. &
                .not. ieee_is_finite(x) .or. .not. ieee_is_finite(y)) then
            value = complex_nan()
        else if (y == 0.0_dp .and. x <= 0.0_dp .and. x == anint(x)) then
            value = complex_nan()
        else if (y == 0.0_dp) then
            value = cmplx(gamma_real(x), y, kind=dp)
        else
            value = safe_complex_exp(loggamma_complex(z))
        end if
    end function gamma_complex

    pure elemental function rgamma_complex(z) result(value)
        complex(dp), intent(in) :: z !! complex argument of the reciprocal gamma function
        complex(dp) :: value

        real(dp) :: x
        real(dp) :: y

        x = real(z, dp)
        y = aimag(z)
        if (ieee_is_nan(x) .or. ieee_is_nan(y)) then
            value = complex_nan()
        else if (.not. ieee_is_finite(x) .or. .not. ieee_is_finite(y)) then
            if (x < 0.0_dp .and. ieee_is_finite(y) .and. y == 0.0_dp) then
                value = cmplx(0.0_dp, y, kind=dp)
            else
                value = complex_nan()
            end if
        else if (y == 0.0_dp .and. x <= 0.0_dp .and. x == anint(x)) then
            value = cmplx(0.0_dp, y, kind=dp)
        else if (y == 0.0_dp) then
            value = cmplx(rgamma_real(x), y, kind=dp)
        else
            value = safe_complex_exp(-loggamma_complex(z))
        end if
    end function rgamma_complex

    pure elemental function loggamma_right_half_plane(z) result(value)
        complex(dp), intent(in) :: z !! finite argument with real part at least 0.5
        complex(dp) :: value

        complex(dp) :: compensation
        complex(dp) :: correction
        complex(dp) :: inv_w
        complex(dp) :: inv_w2
        complex(dp) :: series
        complex(dp) :: term
        complex(dp) :: updated
        complex(dp) :: w
        integer :: j
        integer :: shift

        shift = max(0, ceiling(complex_stirling_threshold - real(z, dp)))
        w = z + real(shift, dp)
        inv_w = 1.0_dp / w
        inv_w2 = inv_w * inv_w
        series = cmplx(complex_stirling_coefficients(size(complex_stirling_coefficients)), &
            0.0_dp, kind=dp)
        do j = size(complex_stirling_coefficients) - 1, 1, -1
            series = cmplx(complex_stirling_coefficients(j), 0.0_dp, kind=dp) + &
                inv_w2 * series
        end do
        series = inv_w * series
        value = (w - 0.5_dp) * log(w) - w + &
            cmplx(scifort_log_sqrt_two_pi, 0.0_dp, kind=dp) + series

        correction = cmplx(0.0_dp, 0.0_dp, kind=dp)
        compensation = cmplx(0.0_dp, 0.0_dp, kind=dp)
        do j = 0, shift - 1
            term = log(z + real(j, dp)) - compensation
            updated = correction + term
            compensation = (updated - correction) - term
            correction = updated
        end do
        value = value - correction
    end function loggamma_right_half_plane

    pure elemental function log_sin_pi(z) result(value)
        complex(dp), intent(in) :: z !! finite nonreal complex argument of log(sin(pi*z))
        complex(dp) :: value

        complex(dp) :: scaled
        real(dp) :: a
        real(dp) :: c
        real(dp) :: one_minus_q
        real(dp) :: q
        real(dp) :: s
        real(dp) :: x
        real(dp) :: y

        x = real(z, dp)
        y = aimag(z)
        s = sinpi_reduced(x)
        c = cospi_reduced(x)
        a = scifort_pi * abs(y)
        if (a > 0.25_dp * log(huge(1.0_dp))) then
            q = 0.0_dp
            one_minus_q = 1.0_dp
        else
            q = exp(-2.0_dp * a)
            one_minus_q = -expm1_safe(-2.0_dp * a)
        end if
        scaled = cmplx(s * (1.0_dp + q), &
            sign(1.0_dp, y) * c * one_minus_q, kind=dp)
        value = cmplx(a - scifort_log_two, 0.0_dp, kind=dp) + log(scaled)
    end function log_sin_pi

    pure elemental function sinpi_reduced(x) result(value)
        real(dp), intent(in) :: x !! finite real argument of sin(pi*x)
        real(dp) :: value

        real(dp) :: r

        r = modulo(x, 2.0_dp)
        if (r > 1.0_dp) r = r - 2.0_dp
        if (r == 0.0_dp .or. abs(r) == 1.0_dp) then
            value = sign(0.0_dp, r)
        else if (abs(r) == 0.5_dp) then
            value = sign(1.0_dp, r)
        else
            value = sin(scifort_pi * r)
        end if
    end function sinpi_reduced

    pure elemental function cospi_reduced(x) result(value)
        real(dp), intent(in) :: x !! finite real argument of cos(pi*x)
        real(dp) :: value

        real(dp) :: r

        r = modulo(x, 2.0_dp)
        if (r > 1.0_dp) r = r - 2.0_dp
        if (r == 0.0_dp) then
            value = 1.0_dp
        else if (abs(r) == 1.0_dp) then
            value = -1.0_dp
        else if (abs(r) == 0.5_dp) then
            value = 0.0_dp
        else
            value = cos(scifort_pi * r)
        end if
    end function cospi_reduced

    pure elemental function safe_complex_exp(log_value) result(value)
        complex(dp), intent(in) :: log_value !! complex logarithm to exponentiate safely
        complex(dp) :: value

        real(dp), parameter :: log_huge = log(huge(1.0_dp))
        real(dp) :: angle
        real(dp) :: cosine
        real(dp) :: magnitude_log
        real(dp) :: sine

        magnitude_log = real(log_value, dp)
        angle = aimag(log_value)
        if (ieee_is_nan(magnitude_log) .or. ieee_is_nan(angle)) then
            value = complex_nan()
        else if (magnitude_log <= log_huge) then
            value = exp(log_value)
        else
            cosine = cos(angle)
            sine = sin(angle)
            value = cmplx(infinite_component(cosine), &
                infinite_component(sine), kind=dp)
        end if
    end function safe_complex_exp

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

        value = cmplx(quiet_nan(0.0_dp), quiet_nan(0.0_dp), kind=dp)
    end function complex_nan

    pure elemental function multigammaln(a, d) result(y)
        real(dp), intent(in) :: a !! multivariate-gamma argument, greater than (d-1)/2
        integer, intent(in) :: d !! positive matrix dimension
        real(dp) :: y

        integer :: j
        real(dp) :: d_real
        real(dp) :: term

        if (ieee_is_nan(a) .or. d < 1) then
            y = quiet_nan(a)
            return
        end if

        d_real = real(d, dp)
        if (.not. ieee_is_finite(a)) then
            if (a > 0.0_dp) then
                y = positive_infinity(a)
            else
                y = quiet_nan(a)
            end if
            return
        end if
        if (a <= 0.5_dp * real(d - 1, dp)) then
            y = quiet_nan(a)
            return
        end if

        y = 0.25_dp * d_real * real(d - 1, dp) * log(scifort_pi)
        do j = 1, d
            term = a + 0.5_dp * real(1 - j, dp)
            y = y + gammaln(term)
        end do
    end function multigammaln

    pure elemental function signed_exp(log_abs, sign_value, kind_arg) result(y)
        real(dp), intent(in) :: log_abs !! logarithm of the requested magnitude
        real(dp), intent(in) :: sign_value !! requested sign, normally -1 or +1
        real(dp), intent(in) :: kind_arg !! kind selector for IEEE infinity construction
        real(dp) :: y

        real(dp) :: log_huge

        log_huge = log(huge(1.0_dp))
        if (log_abs > log_huge) then
            y = sign(positive_infinity(kind_arg), sign_value)
        else
            y = sign(exp(log_abs), sign_value)
        end if
    end function signed_exp

end module scifort_gamma_family
