! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Complex Faddeeva/error/Fresnel functions. This independent implementation
! follows DLMF definitions, interrelations, series, integral representations,
! continued fractions, and symmetry identities from Chapter 7.

module scifort_faddeeva
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi, scifort_sqrt_two
    use scifort_error_fresnel, only : dawsn_real => dawsn, erfi_real => erfi, &
        fresnel_real => fresnel
    use scifort_kinds, only : dp
    use scifort_math, only : positive_infinity, quiet_nan
    implicit none
    private

    real(dp), parameter :: gl16_nodes(8) = [ &
        0.09501250983763744051_dp, &
        0.2816035507792589154_dp, &
        0.4580167776572273697_dp, &
        0.6178762444026437706_dp, &
        0.7554044083550029987_dp, &
        0.8656312023878317552_dp, &
        0.9445750230732326003_dp, &
        0.9894009349916499385_dp]
    real(dp), parameter :: gl16_weights(8) = [ &
        0.1894506104550686409_dp, &
        0.1826034150449236393_dp, &
        0.1691565193950026469_dp, &
        0.1495959888165767082_dp, &
        0.1246289712555340712_dp, &
        0.09515851168249260528_dp, &
        0.06225352393864745648_dp, &
        0.02715245941175417621_dp]

    complex(dp), parameter :: imaginary_unit = cmplx(0.0_dp, 1.0_dp, kind=dp)
    real(dp), parameter :: inv_sqrt_pi = 1.0_dp / sqrt(scifort_pi)
    real(dp), parameter :: two_inv_sqrt_pi = 2.0_dp / sqrt(scifort_pi)

    public :: dawsn_complex
    public :: erf_complex
    public :: erfc_complex
    public :: erfcx_complex
    public :: erfi_complex
    public :: fresnel_complex
    public :: modfresnelm
    public :: modfresnelp
    public :: voigt_profile
    public :: wofz

contains

    pure elemental function voigt_profile(x, sigma, gamma) result(value)
        real(dp), intent(in) :: x !! real profile coordinate
        real(dp), intent(in) :: sigma !! Gaussian standard-deviation parameter; SciPy permits either sign
        real(dp), intent(in) :: gamma !! Cauchy half-width parameter; SciPy permits either sign
        real(dp) :: value

        real(dp), parameter :: sqrt_two_pi = sqrt(2.0_dp * scifort_pi)
        real(dp) :: scale
        real(dp) :: sx
        real(dp) :: sg
        real(dp) :: sumsq
        complex(dp) :: z

        if (ieee_is_nan(x) .or. ieee_is_nan(sigma) .or. ieee_is_nan(gamma)) then
            value = quiet_nan(x + sigma + gamma)
            return
        end if

        if (sigma == 0.0_dp) then
            if (gamma == 0.0_dp) then
                if (x == 0.0_dp) then
                    value = positive_infinity(x)
                else
                    value = 0.0_dp
                end if
                return
            end if
            if (.not. ieee_is_finite(gamma)) then
                value = quiet_nan(gamma)
                return
            end if
            if (.not. ieee_is_finite(x)) then
                value = 0.0_dp
                return
            end if
            scale = max(abs(x), abs(gamma))
            sx = x / scale
            sg = gamma / scale
            sumsq = sx * sx + sg * sg
            value = ((gamma / scale) / sumsq / scifort_pi) / scale
            return
        end if

        z = cmplx(x, gamma, kind=dp) / (scifort_sqrt_two * sigma)
        value = real(wofz(z), dp) / sigma / sqrt_two_pi
    end function voigt_profile

    pure elemental subroutine modfresnelp(x, fp, kp)
        real(dp), intent(in) :: x !! real lower limit of the modified Fresnel tail integral
        complex(dp), intent(out) :: fp !! integral of exp(i*t**2) from x to positive infinity
        complex(dp), intent(out) :: kp !! scaled modified positive Fresnel integral

        real(dp), parameter :: inv_sqrt_two = 1.0_dp / sqrt(2.0_dp)
        real(dp), parameter :: sqrt_huge = sqrt(huge(1.0_dp))
        real(dp), parameter :: sqrt_pi = sqrt(scifort_pi)
        real(dp), parameter :: sqrt_pi_over_two = sqrt(0.5_dp * scifort_pi)
        complex(dp) :: k_positive
        complex(dp) :: phase_positive
        complex(dp) :: rotated
        complex(dp) :: total_integral
        real(dp) :: ax
        real(dp) :: cos_square
        real(dp) :: sin_square

        if (ieee_is_nan(x) .or. .not. ieee_is_finite(x)) then
            fp = cmplx(quiet_nan(0.0_dp), quiet_nan(0.0_dp), kind=dp)
            kp = fp
            return
        end if

        ax = abs(x)
        if (ax > sqrt_huge) then
            fp = cmplx(quiet_nan(0.0_dp), quiet_nan(0.0_dp), kind=dp)
            kp = fp
            return
        end if

        rotated = cmplx(ax * inv_sqrt_two, ax * inv_sqrt_two, kind=dp)
        k_positive = 0.5_dp * wofz(rotated)
        call square_phase(ax, sin_square, cos_square)
        phase_positive = cmplx((cos_square - sin_square) * inv_sqrt_two, &
            (sin_square + cos_square) * inv_sqrt_two, kind=dp)

        if (x >= 0.0_dp) then
            kp = k_positive
            fp = sqrt_pi * phase_positive * kp
        else
            total_integral = sqrt_pi_over_two * cmplx(1.0_dp, 1.0_dp, kind=dp)
            fp = total_integral - sqrt_pi * phase_positive * k_positive
            kp = cmplx(cos_square, -sin_square, kind=dp) - k_positive
        end if
    end subroutine modfresnelp

    pure elemental subroutine modfresnelm(x, fm, km)
        real(dp), intent(in) :: x !! real lower limit of the modified Fresnel tail integral
        complex(dp), intent(out) :: fm !! integral of exp(-i*t**2) from x to positive infinity
        complex(dp), intent(out) :: km !! scaled modified negative Fresnel integral

        complex(dp) :: fp
        complex(dp) :: kp

        call modfresnelp(x, fp, kp)
        fm = conjg(fp)
        km = conjg(kp)
    end subroutine modfresnelm

    pure subroutine square_phase(x, sin_square, cos_square)
        real(dp), intent(in) :: x !! finite nonnegative argument with x**2 representable
        real(dp), intent(out) :: sin_square !! sine of the exact binary64 square x**2
        real(dp), intent(out) :: cos_square !! cosine of the exact binary64 square x**2

        real(dp), parameter :: splitter = 134217729.0_dp
        real(dp) :: cos_error
        real(dp) :: cos_product
        real(dp) :: error_term
        real(dp) :: high
        real(dp) :: low
        real(dp) :: product
        real(dp) :: sin_error
        real(dp) :: sin_product
        real(dp) :: split

        split = splitter * x
        high = split - (split - x)
        low = x - high
        product = x * x
        error_term = ((high * high - product) + 2.0_dp * high * low) + low * low

        sin_product = sin(product)
        cos_product = cos(product)
        sin_error = sin(error_term)
        cos_error = cos(error_term)
        sin_square = sin_product * cos_error + cos_product * sin_error
        cos_square = cos_product * cos_error - sin_product * sin_error
    end subroutine square_phase

    pure elemental function wofz(z) result(w)
        complex(dp), intent(in) :: z !! complex argument of the Faddeeva function
        complex(dp) :: w

        real(dp), parameter :: log_huge = log(huge(1.0_dp))
        real(dp), parameter :: log_tiny = log(tiny(1.0_dp))
        real(dp), parameter :: sqrt_huge = sqrt(huge(1.0_dp))
        complex(dp) :: reflected
        real(dp) :: ax
        real(dp) :: ay
        real(dp) :: exponent
        real(dp) :: phase
        real(dp) :: rexp
        real(dp) :: x
        real(dp) :: y

        x = real(z, dp)
        y = aimag(z)

        if (ieee_is_nan(x) .or. ieee_is_nan(y)) then
            w = cmplx(quiet_nan(0.0_dp), quiet_nan(0.0_dp), kind=dp)
            return
        end if

        if (.not. ieee_is_finite(x) .or. .not. ieee_is_finite(y)) then
            call wofz_infinite(x, y, w)
            return
        end if

        if (y >= 0.0_dp) then
            w = wofz_upper(z)
            return
        end if

        if (abs(z) <= 1.5_dp) then
            w = wofz_series(z)
            return
        end if

        reflected = wofz_upper(-z)
        ax = abs(x)
        ay = abs(y)
        if (max(ax, ay) > 0.5_dp * sqrt_huge) then
            if (ay < ax) then
                w = -reflected
            else if (ay > ax) then
                rexp = positive_infinity(0.0_dp)
                w = cmplx(rexp, rexp, kind=dp)
            else
                w = cmplx(quiet_nan(0.0_dp), quiet_nan(0.0_dp), kind=dp)
            end if
            return
        end if

        exponent = (ay - ax) * (ay + ax)
        if (exponent < log_tiny) then
            w = -reflected
            return
        end if

        phase = -2.0_dp * x * y
        if (exponent <= log_huge - log(2.0_dp)) then
            w = 2.0_dp * exp(cmplx(exponent, phase, kind=dp)) - reflected
        else
            rexp = positive_infinity(0.0_dp)
            w = cmplx(infinite_component(cos(phase), rexp), &
                infinite_component(sin(phase), rexp), kind=dp) - reflected
        end if
    end function wofz


    pure elemental function erf_complex(z) result(value)
        complex(dp), intent(in) :: z !! complex argument of the error function
        complex(dp) :: value

        real(dp) :: x
        real(dp) :: y

        x = real(z, dp)
        y = aimag(z)
        if (ieee_is_nan(x) .or. ieee_is_nan(y)) then
            value = cmplx(quiet_nan(0.0_dp), quiet_nan(0.0_dp), kind=dp)
        else if (.not. ieee_is_finite(x) .or. .not. ieee_is_finite(y)) then
            value = cmplx(quiet_nan(0.0_dp), quiet_nan(0.0_dp), kind=dp)
        else if (y == 0.0_dp) then
            value = cmplx(erf(x), y, kind=dp)
        else if (x == 0.0_dp) then
            value = cmplx(x, erfi_real(y), kind=dp)
        else if (abs(z) <= 1.5_dp) then
            value = erf_series(z)
        else if (x < 0.0_dp) then
            value = -erf_positive(-z)
        else
            value = erf_positive(z)
        end if
    end function erf_complex

    pure elemental function erfc_complex(z) result(value)
        complex(dp), intent(in) :: z !! complex argument of the complementary error function
        complex(dp) :: value

        real(dp) :: x
        real(dp) :: y

        x = real(z, dp)
        y = aimag(z)
        if (ieee_is_nan(x) .or. ieee_is_nan(y)) then
            value = cmplx(quiet_nan(0.0_dp), quiet_nan(0.0_dp), kind=dp)
        else if (.not. ieee_is_finite(x) .or. .not. ieee_is_finite(y)) then
            value = cmplx(quiet_nan(0.0_dp), quiet_nan(0.0_dp), kind=dp)
        else if (y == 0.0_dp) then
            value = cmplx(erfc(x), -y, kind=dp)
        else if (x < 0.0_dp) then
            value = cmplx(2.0_dp, 0.0_dp, kind=dp) - erfc_positive(-z)
        else
            value = erfc_positive(z)
        end if
    end function erfc_complex

    pure elemental function erfi_complex(z) result(value)
        complex(dp), intent(in) :: z !! complex argument of the imaginary error function
        complex(dp) :: value

        real(dp) :: x
        real(dp) :: y

        x = real(z, dp)
        y = aimag(z)
        if (ieee_is_nan(x) .or. ieee_is_nan(y)) then
            value = cmplx(quiet_nan(0.0_dp), quiet_nan(0.0_dp), kind=dp)
        else if (.not. ieee_is_finite(x) .or. .not. ieee_is_finite(y)) then
            value = cmplx(quiet_nan(0.0_dp), quiet_nan(0.0_dp), kind=dp)
        else if (y == 0.0_dp) then
            value = cmplx(erfi_real(x), y, kind=dp)
        else
            value = -imaginary_unit * erf_complex(imaginary_unit * z)
        end if
    end function erfi_complex

    pure elemental subroutine fresnel_complex(z, s, c)
        complex(dp), intent(in) :: z !! complex argument of the Fresnel integrals
        complex(dp), intent(out) :: s !! Fresnel S integral from zero to z
        complex(dp), intent(out) :: c !! Fresnel C integral from zero to z

        real(dp), parameter :: half_sqrt_pi = 0.5_dp * sqrt(scifort_pi)
        complex(dp) :: minus_combination
        complex(dp) :: plus_combination
        complex(dp) :: zeta_minus
        complex(dp) :: zeta_plus
        real(dp) :: c_axis
        real(dp) :: s_axis
        real(dp) :: x
        real(dp) :: y

        x = real(z, dp)
        y = aimag(z)
        if (ieee_is_nan(x) .or. ieee_is_nan(y) .or. &
                .not. ieee_is_finite(x) .or. .not. ieee_is_finite(y)) then
            s = cmplx(quiet_nan(0.0_dp), quiet_nan(0.0_dp), kind=dp)
            c = s
            return
        end if

        if (y == 0.0_dp) then
            call fresnel_real(x, s_axis, c_axis)
            s = cmplx(s_axis, y, kind=dp)
            c = cmplx(c_axis, y, kind=dp)
            return
        end if

        if (x == 0.0_dp) then
            call fresnel_real(y, s_axis, c_axis)
            s = cmplx(x, -s_axis, kind=dp)
            c = cmplx(x, c_axis, kind=dp)
            return
        end if

        if (abs(z) <= 20.0_dp .and. abs(abs(x) - abs(y)) <= 1.5_dp) then
            call fresnel_series_complex(z, s, c)
            return
        end if

        zeta_plus = half_sqrt_pi * cmplx(1.0_dp, -1.0_dp, kind=dp) * z
        zeta_minus = half_sqrt_pi * cmplx(1.0_dp, 1.0_dp, kind=dp) * z
        plus_combination = 0.5_dp * cmplx(1.0_dp, 1.0_dp, kind=dp) * &
            erf_complex(zeta_plus)
        minus_combination = 0.5_dp * cmplx(1.0_dp, -1.0_dp, kind=dp) * &
            erf_complex(zeta_minus)
        c = 0.5_dp * (plus_combination + minus_combination)
        s = -0.5_dp * imaginary_unit * (plus_combination - minus_combination)
    end subroutine fresnel_complex

    pure subroutine fresnel_series_complex(z, s, c)
        complex(dp), intent(in) :: z !! finite complex argument in the well-conditioned series sector
        complex(dp), intent(out) :: s !! Fresnel S integral from zero to z
        complex(dp), intent(out) :: c !! Fresnel C integral from zero to z

        integer :: n
        complex(dp) :: tc
        complex(dp) :: ts
        complex(dp) :: z4
        real(dp), parameter :: half_pi = 0.5_dp * scifort_pi

        z4 = z**4
        tc = z
        c = tc
        ts = half_pi * z**3 / 3.0_dp
        s = ts
        do n = 0, 500
            tc = tc * (-half_pi**2 * z4) * real(4 * n + 1, dp) / &
                real((2 * n + 1) * (2 * n + 2) * (4 * n + 5), dp)
            c = c + tc
            ts = ts * (-half_pi**2 * z4) * real(4 * n + 3, dp) / &
                real((2 * n + 2) * (2 * n + 3) * (4 * n + 7), dp)
            s = s + ts
            if (max(abs(tc), abs(ts)) <= &
                    epsilon(1.0_dp) * max(1.0_dp, abs(s), abs(c))) exit
        end do
    end subroutine fresnel_series_complex

    pure elemental function dawsn_complex(z) result(value)
        complex(dp), intent(in) :: z !! complex argument of Dawson's integral
        complex(dp) :: value

        real(dp), parameter :: half_sqrt_pi = 0.5_dp * sqrt(scifort_pi)
        real(dp) :: x
        real(dp) :: y

        x = real(z, dp)
        y = aimag(z)
        if (ieee_is_nan(x) .or. ieee_is_nan(y)) then
            value = cmplx(quiet_nan(0.0_dp), quiet_nan(0.0_dp), kind=dp)
        else if (.not. ieee_is_finite(x) .or. .not. ieee_is_finite(y)) then
            value = cmplx(quiet_nan(0.0_dp), quiet_nan(0.0_dp), kind=dp)
        else if (y == 0.0_dp) then
            value = cmplx(dawsn_real(x), y, kind=dp)
        else if (abs(z) <= 1.5_dp) then
            value = dawsn_series(z)
        else
            value = imaginary_unit * half_sqrt_pi * (exp(-z * z) - wofz(z))
        end if
    end function dawsn_complex

    pure function erf_positive(z) result(value)
        complex(dp), intent(in) :: z !! finite argument with nonnegative real part
        complex(dp) :: value

        value = cmplx(1.0_dp, 0.0_dp, kind=dp) - erfc_positive(z)
    end function erf_positive

    pure function erfc_positive(z) result(value)
        complex(dp), intent(in) :: z !! finite argument with nonnegative real part
        complex(dp) :: value

        value = exp(-z * z) * wofz(imaginary_unit * z)
    end function erfc_positive

    pure function erf_series(z) result(value)
        complex(dp), intent(in) :: z !! finite argument with magnitude at most 1.5
        complex(dp) :: value

        integer :: n
        complex(dp) :: sumv
        complex(dp) :: term
        complex(dp) :: z2

        term = z
        sumv = term
        z2 = z * z
        do n = 1, 200
            term = term * (-z2) / real(n, dp)
            sumv = sumv + term / real(2 * n + 1, dp)
            if (abs(term) / real(2 * n + 1, dp) <= &
                    epsilon(1.0_dp) * max(1.0_dp, abs(sumv))) exit
        end do
        value = two_inv_sqrt_pi * sumv
    end function erf_series

    pure function dawsn_series(z) result(value)
        complex(dp), intent(in) :: z !! finite argument with magnitude at most 1.5
        complex(dp) :: value

        integer :: n
        complex(dp) :: sumv
        complex(dp) :: term
        complex(dp) :: z2

        term = z
        sumv = term
        z2 = z * z
        do n = 1, 200
            term = term * (-2.0_dp * z2) / real(2 * n + 1, dp)
            sumv = sumv + term
            if (abs(term) <= epsilon(1.0_dp) * max(1.0_dp, abs(sumv))) exit
        end do
        value = sumv
    end function dawsn_series

    pure elemental function erfcx_complex(z) result(value)
        complex(dp), intent(in) :: z !! complex argument of exp(z**2)*erfc(z)
        complex(dp) :: value

        value = wofz(imaginary_unit * z)
    end function erfcx_complex

    pure function wofz_upper(z) result(w)
        complex(dp), intent(in) :: z !! finite argument in the closed upper half-plane
        complex(dp) :: w

        real(dp) :: radius
        real(dp) :: x
        real(dp) :: y

        x = real(z, dp)
        y = aimag(z)
        radius = abs(z)

        if (y == 0.0_dp .and. abs(x) <= 26.0_dp) then
            w = cmplx(exp(-x * x), two_inv_sqrt_pi * dawsn_real(x), kind=dp)
        else if (radius <= 1.5_dp) then
            w = wofz_series(z)
        else if (radius >= 6.0_dp) then
            w = wofz_continued_fraction(z)
        else if (y <= 0.5_dp) then
            w = wofz_real_axis_taylor(x, y)
        else
            w = wofz_integral(z)
        end if
    end function wofz_upper

    pure function wofz_series(z) result(w)
        complex(dp), intent(in) :: z !! finite complex argument with modest magnitude
        complex(dp) :: w

        integer :: n
        complex(dp) :: current
        complex(dp) :: previous
        complex(dp) :: sumv
        complex(dp) :: term
        complex(dp) :: z2

        previous = cmplx(1.0_dp, 0.0_dp, kind=dp)
        current = imaginary_unit * two_inv_sqrt_pi * z
        sumv = previous + current
        z2 = z * z

        do n = 0, 198
            if (mod(n, 2) == 0) then
                term = -z2 * previous / (0.5_dp * real(n, dp) + 1.0_dp)
                previous = term
            else
                term = -z2 * current / (0.5_dp * real(n, dp) + 1.0_dp)
                current = term
            end if
            sumv = sumv + term
            if (n >= 8 .and. abs(term) <= epsilon(1.0_dp) * max(1.0_dp, abs(sumv))) exit
        end do
        w = sumv
    end function wofz_series

    pure function wofz_real_axis_taylor(x, y) result(w)
        real(dp), intent(in) :: x !! finite real-axis anchor
        real(dp), intent(in) :: y !! positive imaginary displacement, at most 0.5
        complex(dp) :: w

        integer :: n
        complex(dp) :: h
        complex(dp) :: next_term
        complex(dp) :: previous_term
        complex(dp) :: current_term
        complex(dp) :: sumv
        complex(dp) :: w0

        h = cmplx(0.0_dp, y, kind=dp)
        w0 = cmplx(exp(-x * x), two_inv_sqrt_pi * dawsn_real(x), kind=dp)
        previous_term = w0
        current_term = (-2.0_dp * x * w0 + imaginary_unit * two_inv_sqrt_pi) * h
        sumv = previous_term + current_term

        do n = 1, 100
            next_term = (-2.0_dp * x * h / real(n + 1, dp)) * current_term - &
                (2.0_dp * h * h / real(n + 1, dp)) * previous_term
            sumv = sumv + next_term
            if (abs(next_term) <= epsilon(1.0_dp) * max(1.0_dp, abs(sumv))) exit
            previous_term = current_term
            current_term = next_term
        end do
        w = sumv
    end function wofz_real_axis_taylor

    pure function wofz_integral(z) result(w)
        complex(dp), intent(in) :: z !! finite upper-half-plane argument in the transition region
        complex(dp) :: w

        integer, parameter :: nseg = 72
        real(dp), parameter :: bound = 9.0_dp
        integer :: j
        integer :: k
        real(dp) :: a
        real(dp) :: b
        real(dp) :: half
        real(dp) :: mid
        real(dp) :: t1
        real(dp) :: t2
        complex(dp) :: sumv

        sumv = cmplx(0.0_dp, 0.0_dp, kind=dp)
        do j = 0, nseg - 1
            a = -bound + 2.0_dp * bound * real(j, dp) / real(nseg, dp)
            b = -bound + 2.0_dp * bound * real(j + 1, dp) / real(nseg, dp)
            mid = 0.5_dp * (a + b)
            half = 0.5_dp * (b - a)
            do k = 1, size(gl16_nodes)
                t1 = mid - half * gl16_nodes(k)
                t2 = mid + half * gl16_nodes(k)
                sumv = sumv + half * gl16_weights(k) * ( &
                    exp(-t1 * t1) / (z - t1) + exp(-t2 * t2) / (z - t2))
            end do
        end do
        w = imaginary_unit * sumv / scifort_pi
    end function wofz_integral

    pure function wofz_continued_fraction(z) result(w)
        complex(dp), intent(in) :: z !! finite upper-half-plane argument with magnitude at least 6
        complex(dp) :: w

        integer, parameter :: n_terms = 80
        integer :: k
        complex(dp) :: tail

        tail = cmplx(0.0_dp, 0.0_dp, kind=dp)
        do k = n_terms, 1, -1
            tail = (0.5_dp * real(k, dp)) / (z - tail)
        end do
        w = imaginary_unit * inv_sqrt_pi / (z - tail)
    end function wofz_continued_fraction

    pure subroutine wofz_infinite(x, y, w)
        real(dp), intent(in) :: x !! real component, possibly infinite
        real(dp), intent(in) :: y !! imaginary component, possibly infinite
        complex(dp), intent(out) :: w !! limiting Faddeeva value

        real(dp) :: inf

        inf = positive_infinity(0.0_dp)
        if (ieee_is_nan(x) .or. ieee_is_nan(y)) then
            w = cmplx(quiet_nan(0.0_dp), quiet_nan(0.0_dp), kind=dp)
        else if (y > 0.0_dp) then
            w = cmplx(0.0_dp, 0.0_dp, kind=dp)
        else if (y < 0.0_dp .and. ieee_is_finite(x)) then
            w = cmplx(inf, 0.0_dp, kind=dp)
        else if (.not. ieee_is_finite(x) .and. ieee_is_finite(y)) then
            w = cmplx(0.0_dp, sign(0.0_dp, x), kind=dp)
        else
            w = cmplx(quiet_nan(0.0_dp), quiet_nan(0.0_dp), kind=dp)
        end if
    end subroutine wofz_infinite

    pure function infinite_component(direction, inf) result(value)
        real(dp), intent(in) :: direction !! finite trigonometric direction cosine or sine
        real(dp), intent(in) :: inf !! positive infinity of the working kind
        real(dp) :: value

        if (direction == 0.0_dp) then
            value = direction
        else
            value = sign(inf, direction)
        end if
    end function infinite_component

end module scifort_faddeeva
