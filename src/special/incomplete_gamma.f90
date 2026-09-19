! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Regularized incomplete gamma functions P(a, x) and Q(a, x) and their
! inverses with respect to x.
!
! Independent implementation from the following public specifications:
! - DLMF 8.7.1: series for P(a, x);
! - DLMF 8.7.3 with 8.7.1: series used for Q(a, x) when a < 1 and x is small;
! - DLMF 8.9.2: continued fraction for Q(a, x), evaluated with the modified
!   Lentz method (Lentz 1976; Thompson and Barnett 1986);
! - DLMF 5.7.3 and 5.11.1 through scifort_log_gamma;
! - Wilson and Hilferty (1931): starting values for the inverse.
!
! This module is an implementation layer. The supported entry points are
! re-exported by scifort_special.

module scifort_incomplete_gamma
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_log_sqrt_two_pi
    use scifort_kinds, only : dp
    use scifort_log_gamma, only : log1p_minus_x, log_gamma_one_plus, &
        stirling_remainder
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, &
        positive_infinity, quiet_nan
    use scifort_normal, only : normal_isf, normal_ppf
    implicit none
    private

    public :: gammainc
    public :: gammaincc
    public :: gammainccinv
    public :: gammaincinv
    public :: log_gamma_kernel
    public :: log_gammainc
    public :: log_gammaincc

    ! Shape parameter above which the Stirling form of the prefactor is used.
    real(dp), parameter :: large_shape = 10.0_dp

    integer, parameter :: max_inverse_iterations = 200

contains

    ! Regularized lower incomplete gamma function P(a, x).
    pure elemental function gammainc(a, x) result(p)
        real(dp), intent(in) :: a !! shape, >= 0
        real(dp), intent(in) :: x !! upper limit of integration, >= 0
        real(dp) :: p

        real(dp) :: logp
        real(dp) :: logq
        real(dp) :: q

        if (ieee_is_nan(a) .or. ieee_is_nan(x)) then
            p = quiet_nan(x)
        else if (a < 0.0_dp .or. x < 0.0_dp) then
            p = quiet_nan(x)
        else if (a == 0.0_dp) then
            if (x > 0.0_dp) then
                p = 1.0_dp
            else
                p = quiet_nan(x)
            end if
        else if (.not. ieee_is_finite(a)) then
            if (ieee_is_finite(x)) then
                p = 0.0_dp
            else
                p = quiet_nan(x)
            end if
        else if (x == 0.0_dp) then
            p = 0.0_dp
        else if (.not. ieee_is_finite(x)) then
            p = 1.0_dp
        else
            call incomplete_gamma_core(a, x, p, q, logp, logq)
        end if
    end function gammainc

    ! Regularized upper incomplete gamma function Q(a, x) = 1 - P(a, x),
    ! computed without forming the difference.
    pure elemental function gammaincc(a, x) result(q)
        real(dp), intent(in) :: a !! shape, >= 0
        real(dp), intent(in) :: x !! lower limit of integration, >= 0
        real(dp) :: q

        real(dp) :: logp
        real(dp) :: logq
        real(dp) :: p

        if (ieee_is_nan(a) .or. ieee_is_nan(x)) then
            q = quiet_nan(x)
        else if (a < 0.0_dp .or. x < 0.0_dp) then
            q = quiet_nan(x)
        else if (a == 0.0_dp) then
            if (x > 0.0_dp) then
                q = 0.0_dp
            else
                q = quiet_nan(x)
            end if
        else if (.not. ieee_is_finite(a)) then
            if (ieee_is_finite(x)) then
                q = 1.0_dp
            else
                q = quiet_nan(x)
            end if
        else if (x == 0.0_dp) then
            q = 1.0_dp
        else if (.not. ieee_is_finite(x)) then
            q = 0.0_dp
        else
            call incomplete_gamma_core(a, x, p, q, logp, logq)
        end if
    end function gammaincc

    ! log(P(a, x)), accurate when P(a, x) underflows or is close to one.
    pure elemental function log_gammainc(a, x) result(logp)
        real(dp), intent(in) :: a !! shape, >= 0
        real(dp), intent(in) :: x !! upper limit of integration, >= 0
        real(dp) :: logp

        real(dp) :: logq
        real(dp) :: p
        real(dp) :: q

        if (ieee_is_nan(a) .or. ieee_is_nan(x)) then
            logp = quiet_nan(x)
        else if (a <= 0.0_dp .or. x < 0.0_dp .or. .not. ieee_is_finite(a)) then
            logp = log_nonnegative(gammainc(a, x))
        else if (x == 0.0_dp) then
            logp = negative_infinity(x)
        else if (.not. ieee_is_finite(x)) then
            logp = 0.0_dp
        else
            call incomplete_gamma_core(a, x, p, q, logp, logq)
        end if
    end function log_gammainc

    ! log(Q(a, x)), accurate when Q(a, x) underflows or is close to one.
    pure elemental function log_gammaincc(a, x) result(logq)
        real(dp), intent(in) :: a !! shape, >= 0
        real(dp), intent(in) :: x !! lower limit of integration, >= 0
        real(dp) :: logq

        real(dp) :: logp
        real(dp) :: p
        real(dp) :: q

        if (ieee_is_nan(a) .or. ieee_is_nan(x)) then
            logq = quiet_nan(x)
        else if (a <= 0.0_dp .or. x < 0.0_dp .or. .not. ieee_is_finite(a)) then
            logq = log_nonnegative(gammaincc(a, x))
        else if (x == 0.0_dp) then
            logq = 0.0_dp
        else if (.not. ieee_is_finite(x)) then
            logq = negative_infinity(x)
        else
            call incomplete_gamma_core(a, x, p, q, logp, logq)
        end if
    end function log_gammaincc

    ! Inverse of P(a, x) with respect to x.
    pure elemental function gammaincinv(a, p) result(x)
        real(dp), intent(in) :: a !! shape, finite and > 0
        real(dp), intent(in) :: p !! value of P(a, x) in [0, 1]
        real(dp) :: x

        if (.not. valid_inverse_arguments(a, p)) then
            x = quiet_nan(p)
        else if (p <= 0.0_dp) then
            x = 0.0_dp
        else if (p >= 1.0_dp) then
            x = positive_infinity(p)
        else if (p <= 0.5_dp) then
            x = solve_incomplete_gamma(a, p, 1.0_dp - p, .true.)
        else
            x = solve_incomplete_gamma(a, p, 1.0_dp - p, .false.)
        end if
    end function gammaincinv

    ! Inverse of Q(a, x) with respect to x.
    pure elemental function gammainccinv(a, q) result(x)
        real(dp), intent(in) :: a !! shape, finite and > 0
        real(dp), intent(in) :: q !! value of Q(a, x) in [0, 1]
        real(dp) :: x

        if (.not. valid_inverse_arguments(a, q)) then
            x = quiet_nan(q)
        else if (q <= 0.0_dp) then
            x = positive_infinity(q)
        else if (q >= 1.0_dp) then
            x = 0.0_dp
        else if (q <= 0.5_dp) then
            x = solve_incomplete_gamma(a, 1.0_dp - q, q, .false.)
        else
            x = solve_incomplete_gamma(a, 1.0_dp - q, q, .true.)
        end if
    end function gammainccinv

    ! log(x**a * exp(-x) / Gamma(a + 1)) for finite a > 0 and x >= 0.
    ! For large a the Stirling series is used so that the result is formed
    ! from a * (log(1 + t) - t), t = (x - a) / a, without cancellation.
    pure elemental function log_gamma_kernel(a, x) result(y)
        real(dp), intent(in) :: a !! shape, finite and > 0
        real(dp), intent(in) :: x !! argument, >= 0
        real(dp) :: y

        real(dp) :: ratio
        real(dp) :: t

        if (x <= 0.0_dp) then
            y = negative_infinity(x)
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else if (a < large_shape) then
            y = a * log(x) - x - log_gamma_one_plus(a)
        else
            t = (x - a) / a
            if (abs(t) <= 0.25_dp) then
                y = a * log1p_minus_x(t)
            else
                ratio = x / a
                if (ratio > tiny(1.0_dp)) then
                    y = a * log(ratio) - (x - a)
                else
                    y = a * (log(x) - log(a)) - (x - a)
                end if
            end if
            y = y - scifort_log_sqrt_two_pi - 0.5_dp * log(a) - stirling_remainder(a)
        end if
    end function log_gamma_kernel

    ! x**a * exp(-x) / Gamma(a + 1) and its logarithm. When no intermediate
    ! can overflow or underflow and a is moderate, the kernel is formed as a
    ! product so that its relative error does not grow with |log(kernel)|.
    pure elemental subroutine gamma_kernel(a, x, d, logd)
        real(dp), intent(in) :: a !! shape, finite and > 0
        real(dp), intent(in) :: x !! argument, finite and > 0
        real(dp), intent(out) :: d !! x**a exp(-x) / Gamma(a + 1)
        real(dp), intent(out) :: logd !! log(d)

        real(dp), parameter :: safe_exponent = 690.0_dp

        logd = log_gamma_kernel(a, x)
        if (a < large_shape .and. x <= safe_exponent .and. &
            abs(a * log(x)) <= safe_exponent .and. logd >= -safe_exponent) then
            d = x**a * exp(-x) / gamma(1.0_dp + a)
            logd = log(d)
        else
            d = exp(logd)
        end if
    end subroutine gamma_kernel

    ! Compute P, Q, log(P), and log(Q) for finite a > 0 and finite x > 0.
    ! In each region the smaller tail is computed directly and the other is
    ! its complement, which is then at least about 0.37 and well conditioned.
    pure elemental subroutine incomplete_gamma_core(a, x, p, q, logp, logq)
        real(dp), intent(in) :: a !! shape, finite and > 0
        real(dp), intent(in) :: x !! argument, finite and > 0
        real(dp), intent(out) :: p !! P(a, x)
        real(dp), intent(out) :: q !! Q(a, x)
        real(dp), intent(out) :: logp !! log(P(a, x))
        real(dp), intent(out) :: logq !! log(Q(a, x))

        logical :: use_series
        real(dp) :: d
        real(dp) :: logd
        real(dp) :: s

        if (a < 1.0_dp) then
            use_series = x <= 1.5_dp
        else
            use_series = x < a
        end if

        call gamma_kernel(a, x, d, logd)
        if (use_series) then
            s = lower_series(a, x)
            p = d * s
            logp = logd + log(s)
            if (a < 1.0_dp) then
                q = small_shape_upper(a, x)
                if (q < 0.5_dp) then
                    p = 1.0_dp - q
                    logq = log(q)
                    logp = log1p_safe(-q)
                else
                    logq = log1p_safe(-p)
                end if
            else
                q = 1.0_dp - p
                logq = log1p_safe(-p)
            end if
        else
            s = upper_continued_fraction(a, x)
            q = a * d * s
            logq = logd + log(a) + log(s)
            p = 1.0_dp - q
            logp = log1p_safe(-q)
        end if
    end subroutine incomplete_gamma_core

    ! Sum over k >= 0 of x**k / ((a + 1) (a + 2) ... (a + k)), DLMF 8.7.1.
    pure function lower_series(a, x) result(s)
        real(dp), intent(in) :: a !! shape, > 0
        real(dp), intent(in) :: x !! argument, > 0
        real(dp) :: s

        integer :: k
        integer :: max_terms
        real(dp) :: term

        max_terms = iteration_limit(a)
        s = 1.0_dp
        term = 1.0_dp
        do k = 1, max_terms
            term = term * x / (a + real(k, dp))
            s = s + term
            if (term <= 0.5_dp * epsilon(1.0_dp) * s) exit
        end do
    end function lower_series

    ! Q(a, x) for 0 < a < 1 and 0 < x <= 1.5. From DLMF 8.7.1 and 8.7.3,
    ! P(a, x) = x**a / Gamma(a + 1) * (1 + a * T), where
    ! T = sum over k >= 1 of (-x)**k / (k! (a + k)). With
    ! L = a log(x) - log(Gamma(1 + a)), Q = -expm1(L) - exp(L) * a * T.
    ! This avoids the cancellation in 1 - P when a is small.
    pure function small_shape_upper(a, x) result(q)
        real(dp), intent(in) :: a !! shape, 0 < a < 1
        real(dp), intent(in) :: x !! argument, 0 < x <= 1.5
        real(dp) :: q

        integer :: k
        real(dp) :: big_l
        real(dp) :: power
        real(dp) :: t
        real(dp) :: term

        t = 0.0_dp
        power = 1.0_dp
        do k = 1, 60
            power = -power * x / real(k, dp)
            term = power / (a + real(k, dp))
            t = t + term
            if (abs(term) <= 0.25_dp * epsilon(1.0_dp) * abs(t)) exit
        end do

        big_l = a * log(x) - log_gamma_one_plus(a)
        q = -expm1_safe(big_l) - exp(big_l) * a * t
    end function small_shape_upper

    ! Continued fraction F with Q(a, x) = x**a exp(-x) F / Gamma(a), from
    ! DLMF 8.9.2: F = (1/x) / (1 + ((1 - a)/x) / (1 + (1/x) / (1 + ...))).
    ! Partial numerators are 1/x, then (k - a)/x and k/x for k = 1, 2, ...
    ! Every partial denominator is 1. Evaluated by the modified Lentz method.
    pure function upper_continued_fraction(a, x) result(f)
        real(dp), intent(in) :: a !! shape, > 0
        real(dp), intent(in) :: x !! argument, > 0, in the continued-fraction region
        real(dp) :: f

        real(dp), parameter :: tiny_value = 1.0e-300_dp
        integer :: k
        integer :: n
        integer :: max_terms
        real(dp) :: c
        real(dp) :: d
        real(dp) :: delta
        real(dp) :: numerator

        max_terms = 2 * iteration_limit(a)
        f = tiny_value
        c = f
        d = 0.0_dp
        do n = 1, max_terms
            if (n == 1) then
                numerator = 1.0_dp / x
            else if (mod(n, 2) == 0) then
                k = n / 2
                numerator = (real(k, dp) - a) / x
            else
                k = (n - 1) / 2
                numerator = real(k, dp) / x
            end if
            d = 1.0_dp + numerator * d
            if (d == 0.0_dp) d = tiny_value
            c = 1.0_dp + numerator / c
            if (c == 0.0_dp) c = tiny_value
            d = 1.0_dp / d
            delta = c * d
            f = f * delta
            if (abs(delta - 1.0_dp) <= 0.5_dp * epsilon(1.0_dp)) exit
        end do
    end function upper_continued_fraction

    ! Term limit for the series and continued fraction. Near x = a both
    ! need a number of terms proportional to sqrt(a).
    pure function iteration_limit(a) result(n)
        real(dp), intent(in) :: a !! shape, > 0
        integer :: n

        real(dp) :: estimate

        estimate = 500.0_dp + 20.0_dp * sqrt(a)
        if (estimate > real(huge(n) / 4, dp)) then
            n = huge(n) / 4
        else
            n = int(estimate)
        end if
    end function iteration_limit

    ! log(p) for p >= 0 without evaluating log(0).
    pure elemental function log_nonnegative(p) result(y)
        real(dp), intent(in) :: p !! value, >= 0, to take the log of
        real(dp) :: y

        if (p > 0.0_dp) then
            y = log(p)
        else if (p == 0.0_dp) then
            y = negative_infinity(p)
        else
            y = quiet_nan(p)
        end if
    end function log_nonnegative

    pure elemental logical function valid_inverse_arguments(a, p) result(valid)
        real(dp), intent(in) :: a !! shape to check
        real(dp), intent(in) :: p !! probability to check

        valid = ieee_is_finite(a) .and. a > 0.0_dp .and. &
            p >= 0.0_dp .and. p <= 1.0_dp
    end function valid_inverse_arguments

    ! Solve P(a, x) = p when lower is true, or Q(a, x) = q otherwise, for
    ! 0 < p < 1 and q = 1 - p. The residual is taken in logarithmic form,
    ! h = log(P) - log(p) or log(q) - log(Q), which is increasing in
    ! u = log(x). Newton steps in u are safeguarded by a bracket. Because
    ! the log-gamma density is log-concave, log(P) and log(Q) are concave in
    ! u, so the Newton iteration cannot oscillate once it is near the root.
    pure function solve_incomplete_gamma(a, p, q, lower) result(x)
        real(dp), intent(in) :: a !! shape, finite and > 0
        real(dp), intent(in) :: p !! target value of P(a, x), 0 < p < 1
        real(dp), intent(in) :: q !! target value of Q(a, x), 1 - p
        logical, intent(in) :: lower !! .true. to match P, .false. to match Q
        real(dp) :: x

        real(dp), parameter :: max_log_step = 50.0_dp
        integer :: iteration
        real(dp) :: hi
        real(dp) :: lo
        real(dp) :: logd
        real(dp) :: logp
        real(dp) :: logq
        real(dp) :: log_target
        real(dp) :: p_x
        real(dp) :: q_x
        real(dp) :: residual
        real(dp) :: slope
        real(dp) :: step
        real(dp) :: x_new

        if (lower) then
            log_target = log(p)
        else
            log_target = log(q)
        end if

        x = inverse_start(a, p, q, lower)
        lo = 0.0_dp
        hi = positive_infinity(x)

        do iteration = 1, max_inverse_iterations
            call incomplete_gamma_core(a, x, p_x, q_x, logp, logq)
            logd = log_gamma_kernel(a, x)
            if (lower) then
                residual = logp - log_target
                slope = exp(log(a) + logd - logp)
            else
                residual = log_target - logq
                slope = exp(log(a) + logd - logq)
            end if

            if (residual == 0.0_dp) exit
            if (residual < 0.0_dp) then
                lo = x
            else
                hi = x
            end if

            if (slope * max_log_step > abs(residual)) then
                step = -residual / slope
            else
                ! Far from the root the slope can underflow; take the largest
                ! permitted step toward the root without dividing by it.
                step = sign(max_log_step, -residual)
            end if
            x_new = x * exp(step)

            if (.not. (x_new > lo .and. x_new < hi)) then
                if (lo > 0.0_dp .and. ieee_is_finite(hi)) then
                    if (hi > 4.0_dp * lo) then
                        x_new = sqrt(lo) * sqrt(hi)
                    else
                        x_new = lo + 0.5_dp * (hi - lo)
                    end if
                else if (ieee_is_finite(hi)) then
                    x_new = hi * exp(-max_log_step)
                    if (x_new <= 0.0_dp) then
                        ! The root is below hi < 1e-300 and rounds to zero.
                        x = 0.0_dp
                        return
                    end if
                else
                    x_new = lo * exp(max_log_step)
                    if (.not. ieee_is_finite(x_new)) then
                        ! The root lies above the largest finite number.
                        x = positive_infinity(x)
                        return
                    end if
                end if
            end if

            if (abs(x_new - x) <= 2.0_dp * epsilon(1.0_dp) * x) then
                x = x_new
                exit
            end if
            if (ieee_is_finite(hi) .and. hi - lo <= 2.0_dp * epsilon(1.0_dp) * hi) then
                x = x_new
                exit
            end if
            x = x_new
        end do
    end function solve_incomplete_gamma

    ! Starting value for the inverse. For a >= 1 the Wilson-Hilferty
    ! approximation x = a (1 - 1/(9a) + z / (3 sqrt(a)))**3 is used. Otherwise,
    ! and when that approximation is not positive, the small-x behavior
    ! P(a, x) ~ x**a / Gamma(a + 1) or the large-x behavior Q(a, x) ~ exp(-x)
    ! provides the start.
    pure function inverse_start(a, p, q, lower) result(x)
        real(dp), intent(in) :: a !! shape, finite and > 0
        real(dp), intent(in) :: p !! target value of P(a, x)
        real(dp), intent(in) :: q !! target value of Q(a, x), 1 - p
        logical, intent(in) :: lower !! .true. when the root is sought through P
        real(dp) :: x

        real(dp) :: log_x
        real(dp) :: w
        real(dp) :: z

        x = 0.0_dp
        if (a >= 1.0_dp) then
            if (lower) then
                z = normal_ppf(p)
            else
                z = normal_isf(q)
            end if
            w = 1.0_dp - 1.0_dp / (9.0_dp * a) + z / (3.0_dp * sqrt(a))
            if (w > 0.0_dp) x = a * w * w * w
        end if

        if (.not. (x > tiny(1.0_dp) .and. ieee_is_finite(x))) then
            if (lower) then
                log_x = (log(p) + log_gamma_one_plus(a)) / a
                log_x = max(log(tiny(1.0_dp)), min(log(huge(1.0_dp)) - 1.0_dp, log_x))
                x = exp(log_x)
            else
                x = max(1.0_dp, a - log(q))
            end if
        end if
    end function inverse_start

end module scifort_incomplete_gamma
