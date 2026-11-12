! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Real Riemann/Hurwitz zeta functions plus the exponent derivative of the
! positive-q Hurwitz zeta.  Positive-q Euler-Maclaurin evaluation follows
! DLMF 25.11.43.  Riemann-zeta continuation on the negative axis uses the
! reflection formula DLMF 25.4.2.  The implementation is independent of
! SciPy source code.

module scifort_zeta
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : positive_infinity, quiet_nan
    implicit none
    private

    public :: hurwitz_zeta
    public :: hurwitz_zeta_derivative
    public :: zeta
    public :: zetac

    real(dp), parameter :: em_coeff(6) = [ &
        1.0_dp / 12.0_dp, &
        -1.0_dp / 720.0_dp, &
        1.0_dp / 30240.0_dp, &
        -1.0_dp / 1209600.0_dp, &
        1.0_dp / 47900160.0_dp, &
        -691.0_dp / 1307674368000.0_dp]
    real(dp), parameter :: log_two = log(2.0_dp)
    real(dp), parameter :: log_two_pi = log(2.0_dp * scifort_pi)
    real(dp), parameter :: log_huge_dp = log(huge(1.0_dp))
    real(dp), parameter :: log_smallest_dp = log(tiny(1.0_dp)) + log(epsilon(1.0_dp))

contains

    pure elemental function zeta(s, q) result(value)
        real(dp), intent(in) :: s !! real zeta exponent
        real(dp), intent(in), optional :: q !! real Hurwitz shift; omitted for Riemann zeta
        real(dp) :: value

        if (present(q)) then
            value = scipy_hurwitz_zeta(s, q)
        else
            value = riemann_zeta(s)
        end if
    end function zeta

    pure elemental function zetac(s) result(value)
        real(dp), intent(in) :: s !! real Riemann-zeta exponent
        real(dp) :: value

        if (ieee_is_nan(s)) then
            value = s
            return
        end if
        if (.not. ieee_is_finite(s)) then
            if (s > 0.0_dp) then
                value = 0.0_dp
            else
                value = quiet_nan(s)
            end if
            return
        end if
        if (s == 1.0_dp) then
            value = positive_infinity(s)
        else if (s >= 48.0_dp) then
            value = zetac_large_positive(s)
        else if (s > 1.0_dp) then
            value = hurwitz_zeta(s, 2.0_dp)
        else
            value = riemann_zeta(s) - 1.0_dp
        end if
    end function zetac


    pure elemental function zetac_large_positive(s) result(value)
        real(dp), intent(in) :: s !! finite exponent, at least 48
        real(dp) :: value
        real(dp) :: fractional
        real(dp) :: term
        integer :: exponent_two
        integer :: k

        if (s > 1075.0_dp) then
            value = 0.0_dp
            return
        end if

        exponent_two = int(floor(s))
        fractional = s - real(exponent_two, dp)
        value = scale(exp(-fractional * log_two), -exponent_two)

        do k = 3, 32
            term = exp(-s * log(real(k, dp)))
            if (term == 0.0_dp) exit
            value = value + term
            if (term <= 0.25_dp * epsilon(1.0_dp) * value) exit
        end do
    end function zetac_large_positive

    pure elemental function hurwitz_zeta(s, q) result(value)
        real(dp), intent(in) :: s !! exponent, strictly greater than one
        real(dp), intent(in) :: q !! positive Hurwitz shift
        real(dp) :: value
        real(dp) :: derivative

        call zeta_pair(s, q, value, derivative)
    end function hurwitz_zeta

    pure elemental function hurwitz_zeta_derivative(s, q) result(value)
        real(dp), intent(in) :: s !! exponent, strictly greater than one
        real(dp), intent(in) :: q !! positive Hurwitz shift
        real(dp) :: value
        real(dp) :: zeta_value

        call zeta_pair(s, q, zeta_value, value)
    end function hurwitz_zeta_derivative

    pure elemental function riemann_zeta(s) result(value)
        real(dp), intent(in) :: s !! real Riemann-zeta exponent
        real(dp) :: value
        real(dp) :: phase
        real(dp) :: reflected
        real(dp) :: t
        real(dp) :: log_abs

        if (ieee_is_nan(s)) then
            value = s
            return
        end if
        if (.not. ieee_is_finite(s)) then
            if (s > 0.0_dp) then
                value = 1.0_dp
            else
                value = quiet_nan(s)
            end if
            return
        end if
        if (s == 1.0_dp) then
            value = positive_infinity(s)
            return
        end if
        if (s > 1.0_dp) then
            value = hurwitz_zeta(s, 1.0_dp)
            return
        end if
        if (s == 0.0_dp) then
            value = -0.5_dp
            return
        end if
        if (s > -0.01_dp) then
            value = zeta_em_value(s, 1.0_dp)
            return
        end if

        phase = sin_half_pi_reduced(s)
        if (phase == 0.0_dp) then
            value = 0.0_dp
            return
        end if

        t = 1.0_dp - s
        if (t > 64.0_dp) then
            reflected = 1.0_dp
        else
            reflected = hurwitz_zeta(t, 1.0_dp)
        end if

        if (t > 1.0e305_dp) then
            value = sign(positive_infinity(s), phase)
            return
        end if

        log_abs = log_two + (s - 1.0_dp) * log_two_pi + log(abs(phase)) + &
            log_gamma(t) + log(reflected)
        value = signed_exp(log_abs, phase, s)
    end function riemann_zeta

    pure elemental function scipy_hurwitz_zeta(s, q) result(value)
        real(dp), intent(in) :: s !! real exponent; SciPy-compatible real subset requires s >= 1
        real(dp), intent(in) :: q !! real Hurwitz shift
        real(dp) :: value
        real(dp) :: a
        real(dp) :: f
        real(dp) :: n_floor
        real(dp) :: negative_sum
        real(dp) :: parity

        if (ieee_is_nan(s)) then
            value = s
            return
        end if
        if (ieee_is_nan(q)) then
            value = q
            return
        end if
        if (s < 1.0_dp) then
            value = quiet_nan(s)
            return
        end if
        if (s == 1.0_dp) then
            value = positive_infinity(s)
            return
        end if
        if (.not. ieee_is_finite(s)) then
            if (.not. ieee_is_finite(q)) then
                if (q > 0.0_dp) then
                    value = 0.0_dp
                else
                    value = positive_infinity(s)
                end if
            else if (q == 1.0_dp) then
                value = 1.0_dp
            else if (q < 1.0_dp) then
                value = positive_infinity(s)
            else
                value = quiet_nan(s)
            end if
            return
        end if
        if (.not. ieee_is_finite(q)) then
            if (q > 0.0_dp) then
                value = 0.0_dp
            else
                value = positive_infinity(q)
            end if
            return
        end if
        if (q > 0.0_dp) then
            value = hurwitz_zeta(s, q)
            return
        end if
        if (q == anint(q)) then
            value = positive_infinity(q)
            return
        end if
        if (s /= anint(s)) then
            value = quiet_nan(s)
            return
        end if

        ! For q = -(N+f), N >= 0 and 0 < f < 1, separate the finite
        ! negative-base terms from the positive tail. Integer s keeps those
        ! terms real and their common sign is (-1)^s.
        a = -q
        n_floor = aint(a)
        f = a - n_floor
        if (modulo(s, 2.0_dp) == 0.0_dp) then
            parity = 1.0_dp
        else
            parity = -1.0_dp
        end if
        negative_sum = hurwitz_zeta(s, f) - hurwitz_zeta(s, a + 1.0_dp)
        value = hurwitz_zeta(s, 1.0_dp - f) + parity * negative_sum
    end function scipy_hurwitz_zeta

    pure elemental function zeta_em_value(s, q) result(value)
        real(dp), intent(in) :: s !! real exponent different from one
        real(dp), intent(in) :: q !! positive Hurwitz shift
        real(dp) :: value
        integer, parameter :: n_shift = 24
        integer :: i
        integer :: j
        integer :: r
        real(dp) :: correction
        real(dp) :: rising
        real(dp) :: term
        real(dp) :: x

        value = 0.0_dp
        do i = 0, n_shift - 1
            x = q + real(i, dp)
            term = exp(-s * log(x))
            value = value + term
        end do

        x = q + real(n_shift, dp)
        value = value + exp((1.0_dp - s) * log(x)) / (s - 1.0_dp)
        value = value + 0.5_dp * exp(-s * log(x))

        do i = 1, size(em_coeff)
            r = 2 * i - 1
            rising = 1.0_dp
            do j = 0, r - 1
                rising = rising * (s + real(j, dp))
            end do
            correction = em_coeff(i) * rising * exp(-(s + real(r, dp)) * log(x))
            value = value + correction
        end do
    end function zeta_em_value

    pure elemental function sin_half_pi_reduced(s) result(value)
        real(dp), intent(in) :: s !! finite real argument of sin(pi*s/2)
        real(dp) :: value
        real(dp) :: r

        r = modulo(s, 4.0_dp)
        if (r == 0.0_dp .or. r == 2.0_dp) then
            value = 0.0_dp
        else if (r == 1.0_dp) then
            value = 1.0_dp
        else if (r == 3.0_dp) then
            value = -1.0_dp
        else if (r < 1.0_dp) then
            value = sin(0.5_dp * scifort_pi * r)
        else if (r < 2.0_dp) then
            value = sin(0.5_dp * scifort_pi * (2.0_dp - r))
        else if (r < 3.0_dp) then
            value = -sin(0.5_dp * scifort_pi * (r - 2.0_dp))
        else
            value = -sin(0.5_dp * scifort_pi * (4.0_dp - r))
        end if
    end function sin_half_pi_reduced

    pure elemental function signed_exp(log_abs, sign_source, prototype) result(value)
        real(dp), intent(in) :: log_abs !! logarithm of the desired magnitude
        real(dp), intent(in) :: sign_source !! nonzero value supplying the result sign
        real(dp), intent(in) :: prototype !! value supplying result kind for IEEE helpers
        real(dp) :: value

        if (log_abs > log_huge_dp) then
            value = sign(positive_infinity(prototype), sign_source)
        else if (log_abs < log_smallest_dp) then
            value = sign(0.0_dp, sign_source)
        else
            value = sign(exp(log_abs), sign_source)
        end if
    end function signed_exp

    pure subroutine zeta_pair(s, q, value, derivative)
        real(dp), intent(in) :: s !! exponent, strictly greater than one
        real(dp), intent(in) :: q !! positive Hurwitz shift
        real(dp), intent(out) :: value !! Hurwitz zeta zeta(s,q)
        real(dp), intent(out) :: derivative !! derivative with respect to s
        integer, parameter :: n_shift = 24
        integer :: i
        integer :: j
        integer :: r
        real(dp) :: x
        real(dp) :: term
        real(dp) :: rising
        real(dp) :: harmonic
        real(dp) :: correction

        if (.not. ieee_is_finite(s) .or. .not. ieee_is_finite(q) .or. &
            s <= 1.0_dp .or. q <= 0.0_dp) then
            value = quiet_nan(s + q)
            derivative = value
            return
        end if

        value = 0.0_dp
        derivative = 0.0_dp
        do i = 0, n_shift - 1
            x = q + real(i, dp)
            term = exp(-s * log(x))
            value = value + term
            derivative = derivative - log(x) * term
        end do

        x = q + real(n_shift, dp)
        term = exp((1.0_dp - s) * log(x)) / (s - 1.0_dp)
        value = value + term
        derivative = derivative + term * (-log(x) - 1.0_dp / (s - 1.0_dp))

        term = 0.5_dp * exp(-s * log(x))
        value = value + term
        derivative = derivative - log(x) * term

        do i = 1, size(em_coeff)
            r = 2 * i - 1
            rising = 1.0_dp
            harmonic = 0.0_dp
            do j = 0, r - 1
                rising = rising * (s + real(j, dp))
                harmonic = harmonic + 1.0_dp / (s + real(j, dp))
            end do
            correction = em_coeff(i) * rising * exp(-(s + real(r, dp)) * log(x))
            value = value + correction
            derivative = derivative + correction * (harmonic - log(x))
        end do
    end subroutine zeta_pair

end module scifort_zeta
