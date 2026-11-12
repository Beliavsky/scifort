! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Real spherical Bessel functions following scipy.special naming and the
! definitions, series, recurrences, and explicit formulas in DLMF 10.47,
! 10.49, 10.51, and 10.53. Complex arguments are outside this binary64 API.
module scifort_spherical_bessel
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_bessel_i, only : log_besseli_scaled
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan
    implicit none
    private

    integer, parameter, public :: spherical_bessel_status_success = 0
    integer, parameter, public :: spherical_bessel_status_invalid_argument = 1

    public :: riccati_jn
    public :: riccati_yn
    public :: spherical_in
    public :: spherical_jn
    public :: spherical_kn
    public :: spherical_yn

contains

    pure subroutine riccati_jn(n, x, values, derivatives, status)
        integer, intent(in) :: n !! highest nonnegative Riccati-Bessel order
        real(dp), intent(in) :: x !! finite real evaluation argument
        real(dp), allocatable, intent(out) :: values(:) !! x*j_k(x), indexed from 0 through n
        real(dp), allocatable, intent(out) :: derivatives(:) !! d[x*j_k(x)]/dx, indexed 0 through n
        integer, intent(out), optional :: status !! zero on success; nonzero on invalid input

        integer :: k

        if (n < 0 .or. .not. ieee_is_finite(x)) then
            allocate(values(0:-1), derivatives(0:-1))
            if (present(status)) status = spherical_bessel_status_invalid_argument
            return
        end if

        allocate(values(0:n), derivatives(0:n))
        if (x == 0.0_dp) then
            values = 0.0_dp
            derivatives = 0.0_dp
            derivatives(0) = 1.0_dp
        else
            do k = 0, n
                values(k) = x * spherical_jn(k, x)
                derivatives(k) = spherical_jn(k, x) + x * spherical_jn(k, x, .true.)
            end do
        end if
        if (present(status)) status = spherical_bessel_status_success
    end subroutine riccati_jn

    pure subroutine riccati_yn(n, x, values, derivatives, status)
        integer, intent(in) :: n !! highest nonnegative Riccati-Bessel order
        real(dp), intent(in) :: x !! finite real evaluation argument
        real(dp), allocatable, intent(out) :: values(:) !! x*y_k(x), indexed from 0 through n
        real(dp), allocatable, intent(out) :: derivatives(:) !! d[x*y_k(x)]/dx, indexed 0 through n
        integer, intent(out), optional :: status !! zero on success; nonzero on invalid input

        integer :: k

        if (n < 0 .or. .not. ieee_is_finite(x)) then
            allocate(values(0:-1), derivatives(0:-1))
            if (present(status)) status = spherical_bessel_status_invalid_argument
            return
        end if

        allocate(values(0:n), derivatives(0:n))
        if (x == 0.0_dp) then
            values(0) = -1.0_dp
            derivatives(0) = 0.0_dp
            do k = 1, n
                values(k) = negative_infinity(x)
                derivatives(k) = positive_infinity(x)
            end do
        else
            do k = 0, n
                values(k) = x * spherical_yn(k, x)
                derivatives(k) = spherical_yn(k, x) + x * spherical_yn(k, x, .true.)
            end do
        end if
        if (present(status)) status = spherical_bessel_status_success
    end subroutine riccati_yn

    pure elemental function spherical_jn(n, x, derivative) result(y)
        integer, intent(in) :: n !! nonnegative spherical-Bessel order
        real(dp), intent(in) :: x !! real argument
        logical, intent(in), optional :: derivative !! return d/dx when true; default false
        real(dp) :: y

        logical :: want_derivative
        real(dp) :: ax
        real(dp) :: log_sum
        integer :: sign_value

        want_derivative = .false.
        if (present(derivative)) want_derivative = derivative

        if (n < 0 .or. ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if
        if (.not. ieee_is_finite(x)) then
            y = 0.0_dp
            return
        end if
        if (x == 0.0_dp) then
            if (want_derivative) then
                if (n == 1) then
                    y = 1.0_dp / 3.0_dp
                else
                    y = 0.0_dp
                end if
            else if (n == 0) then
                y = 1.0_dp
            else
                y = 0.0_dp
            end if
            return
        end if

        ax = abs(x)
        if (.not. want_derivative) then
            y = spherical_jn_positive(n, ax)
            if (x < 0.0_dp .and. modulo(n, 2) == 1) y = -y
            return
        end if

        if (n == 0) then
            y = -spherical_jn_positive(1, ax)
        else if (ax < 1.0_dp) then
            log_sum = log_spherical_jn_derivative_series(n, ax)
            sign_value = 1
            y = signed_exp_from_log(log_sum, sign_value, x)
        else
            y = spherical_jn_positive(n - 1, ax) - &
                (real(n + 1, dp) / ax) * spherical_jn_positive(n, ax)
        end if
        if (x < 0.0_dp .and. modulo(n, 2) == 0) y = -y
    end function spherical_jn

    pure elemental function spherical_yn(n, x, derivative) result(y)
        integer, intent(in) :: n !! nonnegative spherical-Bessel order
        real(dp), intent(in) :: x !! real argument
        logical, intent(in), optional :: derivative !! return d/dx when true; default false
        real(dp) :: y

        logical :: want_derivative
        real(dp) :: ax

        want_derivative = .false.
        if (present(derivative)) want_derivative = derivative

        if (n < 0 .or. ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if
        if (.not. ieee_is_finite(x)) then
            y = 0.0_dp
            return
        end if
        if (x == 0.0_dp) then
            if (want_derivative) then
                y = positive_infinity(x)
            else
                y = negative_infinity(x)
            end if
            return
        end if

        ax = abs(x)
        if (want_derivative) then
            y = spherical_yn_derivative_positive(n, ax)
            if (x < 0.0_dp .and. modulo(n, 2) == 1) y = -y
        else
            y = spherical_yn_positive(n, ax)
            if (x < 0.0_dp .and. modulo(n, 2) == 0) y = -y
        end if
    end function spherical_yn

    pure elemental function spherical_in(n, x, derivative) result(y)
        integer, intent(in) :: n !! nonnegative modified spherical-Bessel order
        real(dp), intent(in) :: x !! real argument
        logical, intent(in), optional :: derivative !! return d/dx when true; default false
        real(dp) :: y

        logical :: want_derivative
        real(dp) :: ax

        want_derivative = .false.
        if (present(derivative)) want_derivative = derivative

        if (n < 0 .or. ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if
        if (.not. ieee_is_finite(x)) then
            y = positive_infinity(x)
            if (x < 0.0_dp) then
                if (want_derivative) then
                    if (modulo(n, 2) == 0) y = -y
                else
                    if (modulo(n, 2) == 1) y = -y
                end if
            end if
            return
        end if
        if (x == 0.0_dp) then
            if (want_derivative) then
                if (n == 1) then
                    y = 1.0_dp / 3.0_dp
                else
                    y = 0.0_dp
                end if
            else if (n == 0) then
                y = 1.0_dp
            else
                y = 0.0_dp
            end if
            return
        end if

        ax = abs(x)
        if (.not. want_derivative) then
            y = positive_exp_from_log(log_spherical_in_positive(n, ax), x)
            if (x < 0.0_dp .and. modulo(n, 2) == 1) y = -y
            return
        end if

        y = positive_exp_from_log(log_spherical_in_derivative_positive(n, ax), x)
        if (x < 0.0_dp .and. modulo(n, 2) == 0) y = -y
    end function spherical_in

    pure elemental function spherical_kn(n, x, derivative) result(y)
        integer, intent(in) :: n !! nonnegative modified spherical-Bessel order
        real(dp), intent(in) :: x !! real argument; zero is singular
        logical, intent(in), optional :: derivative !! return d/dx when true; default false
        real(dp) :: y

        logical :: want_derivative
        real(dp) :: ax
        real(dp) :: log_a
        real(dp) :: log_b
        real(dp) :: log_sum
        integer :: sign_b
        integer :: sign_value

        want_derivative = .false.
        if (present(derivative)) want_derivative = derivative

        if (n < 0 .or. ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if
        if (.not. ieee_is_finite(x)) then
            if (x > 0.0_dp) then
                if (want_derivative) then
                    y = -0.0_dp
                else
                    y = 0.0_dp
                end if
            else if (want_derivative) then
                y = positive_infinity(x)
            else
                y = negative_infinity(x)
            end if
            return
        end if
        if (x == 0.0_dp) then
            if (want_derivative) then
                y = negative_infinity(x)
            else
                y = positive_infinity(x)
            end if
            return
        end if

        ax = abs(x)
        if (x > 0.0_dp) then
            if (want_derivative) then
                log_a = log_spherical_kn_derivative_positive(n, ax)
                y = -positive_exp_from_log(log_a, x)
            else
                y = positive_exp_from_log(log_spherical_kn_positive(n, ax), x)
            end if
            return
        end if

        ! DLMF 10.49.8 and 10.49.12 imply, for a > 0,
        ! k_n(-a) = -pi*i_n(a) + (-1)^(n+1)*k_n(a).
        ! Combine the two terms in logarithmic form because they can nearly cancel.
        if (want_derivative) then
            log_a = log(scifort_pi) + log_spherical_in_derivative_positive(n, ax)
            log_b = log_spherical_kn_derivative_positive(n, ax)
            sign_b = parity_sign(n + 1)
            call signed_log_add(log_a, 1, log_b, sign_b, log_sum, sign_value)
        else
            log_a = log(scifort_pi) + log_spherical_in_positive(n, ax)
            log_b = log_spherical_kn_positive(n, ax)
            sign_b = parity_sign(n + 1)
            call signed_log_add(log_a, -1, log_b, sign_b, log_sum, sign_value)
        end if
        y = signed_exp_from_log(log_sum, sign_value, x)
    end function spherical_kn

    pure function spherical_jn_positive(n, x) result(y)
        integer, intent(in) :: n !! nonnegative order
        real(dp), intent(in) :: x !! positive finite argument
        real(dp) :: y

        integer :: k
        real(dp) :: jm1
        real(dp) :: jnow
        real(dp) :: jnext

        if (n == 0) then
            y = spherical_j0_base(x)
        else if (n == 1) then
            y = spherical_j1_base(x)
        else if (x < 1.0_dp) then
            y = positive_exp_from_log(log_spherical_jn_series(n, x), x)
        else if (x > real(n, dp)) then
            jm1 = spherical_j0_base(x)
            jnow = spherical_j1_base(x)
            do k = 1, n - 1
                jnext = ((2.0_dp * real(k, dp) + 1.0_dp) / x) * jnow - jm1
                jm1 = jnow
                jnow = jnext
            end do
            y = jnow
        else
            y = spherical_jn_miller(n, x)
        end if
    end function spherical_jn_positive

    pure function spherical_jn_miller(n, x) result(y)
        integer, intent(in) :: n !! order >= 2
        real(dp), intent(in) :: x !! finite argument with 1 <= x <= n
        real(dp) :: y

        integer :: extra
        integer :: k
        integer :: m
        logical :: have_f1
        logical :: have_fn
        real(dp), parameter :: rescale = 1.0e-100_dp
        real(dp), parameter :: rescale_trigger = 1.0e100_dp
        real(dp) :: b0
        real(dp) :: b1
        real(dp) :: f0
        real(dp) :: f1
        real(dp) :: fk
        real(dp) :: fkm1
        real(dp) :: fkp1
        real(dp) :: fn

        extra = max(50, int(2.0_dp * sqrt(40.0_dp * (real(n, dp) + 1.0_dp))) + 10)
        if (n > huge(n) - extra) then
            y = 0.0_dp
            return
        end if
        m = n + extra
        fkp1 = 0.0_dp
        fk = 1.0_dp
        fn = 0.0_dp
        f1 = 0.0_dp
        have_fn = .false.
        have_f1 = .false.

        do k = m, 1, -1
            fkm1 = ((2.0_dp * real(k, dp) + 1.0_dp) / x) * fk - fkp1
            if (k - 1 == n) then
                fn = fkm1
                have_fn = .true.
            end if
            if (k - 1 == 1) then
                f1 = fkm1
                have_f1 = .true.
            end if
            if (abs(fkm1) > rescale_trigger) then
                fkm1 = fkm1 * rescale
                fk = fk * rescale
                if (have_fn) fn = fn * rescale
                if (have_f1) f1 = f1 * rescale
            end if
            fkp1 = fk
            fk = fkm1
        end do

        f0 = fk
        b0 = spherical_j0_base(x)
        b1 = spherical_j1_base(x)
        if (abs(b0) >= abs(b1)) then
            y = fn * (b0 / f0)
        else
            y = fn * (b1 / f1)
        end if
    end function spherical_jn_miller

    pure function spherical_yn_positive(n, x) result(y)
        integer, intent(in) :: n !! nonnegative order
        real(dp), intent(in) :: x !! positive finite argument
        real(dp) :: y

        integer :: k
        real(dp) :: factor
        real(dp) :: product
        real(dp) :: ym1
        real(dp) :: ynow
        real(dp) :: ynext

        ym1 = spherical_y0_base(x)
        if (n == 0) then
            y = ym1
            return
        end if
        ynow = spherical_y1_base(x)
        if (n == 1 .or. .not. ieee_is_finite(ynow)) then
            y = ynow
            return
        end if

        do k = 1, n - 1
            factor = (2.0_dp * real(k, dp) + 1.0_dp) / x
            if (abs(factor) > 1.0_dp) then
                if (abs(ynow) > huge(1.0_dp) / abs(factor)) then
                    y = signed_infinity(ynow, x)
                    return
                end if
            end if
            product = factor * ynow
            if (((product > 0.0_dp .and. ym1 < 0.0_dp) .or. &
                    (product < 0.0_dp .and. ym1 > 0.0_dp)) .and. &
                    abs(product) > huge(1.0_dp) - abs(ym1)) then
                y = signed_infinity(product, x)
                return
            end if
            ynext = product - ym1
            ym1 = ynow
            ynow = ynext
        end do
        y = ynow
    end function spherical_yn_positive

    pure function spherical_yn_derivative_positive(n, x) result(y)
        integer, intent(in) :: n !! nonnegative order
        real(dp), intent(in) :: x !! positive finite argument
        real(dp) :: y

        real(dp) :: a
        real(dp) :: b
        real(dp) :: factor
        real(dp) :: ynow

        if (n == 0) then
            y = -spherical_yn_positive(1, x)
            return
        end if
        a = spherical_yn_positive(n - 1, x)
        ynow = spherical_yn_positive(n, x)
        factor = real(n + 1, dp) / x
        if (.not. ieee_is_finite(ynow)) then
            y = -signed_infinity(ynow, x)
            return
        end if
        if (abs(factor) > 1.0_dp) then
            if (abs(ynow) > huge(1.0_dp) / abs(factor)) then
                y = -signed_infinity(ynow, x)
                return
            end if
        end if
        b = factor * ynow
        if (.not. ieee_is_finite(a)) then
            y = a
        else if (((a > 0.0_dp .and. b < 0.0_dp) .or. &
                (a < 0.0_dp .and. b > 0.0_dp)) .and. &
                abs(a) > huge(1.0_dp) - abs(b)) then
            y = signed_infinity(a, x)
        else
            y = a - b
        end if
    end function spherical_yn_derivative_positive

    pure function log_spherical_jn_series(n, x) result(log_value)
        integer, intent(in) :: n !! nonnegative order
        real(dp), intent(in) :: x !! argument with 0 < x < 1
        real(dp) :: log_value

        integer :: k
        integer, parameter :: max_terms = 200
        real(dp) :: sum_value
        real(dp) :: term

        term = 1.0_dp
        sum_value = 1.0_dp
        do k = 1, max_terms
            term = term * (-0.5_dp * x * x) / &
                (real(k, dp) * (2.0_dp * real(n + k, dp) + 1.0_dp))
            sum_value = sum_value + term
            if (abs(term) <= 0.25_dp * epsilon(1.0_dp) * abs(sum_value)) exit
        end do
        log_value = spherical_series_log_prefactor(n, x) + log(sum_value)
    end function log_spherical_jn_series

    pure function log_spherical_jn_derivative_series(n, x) result(log_value)
        integer, intent(in) :: n !! positive order
        real(dp), intent(in) :: x !! argument with 0 < x < 1
        real(dp) :: log_value

        integer :: k
        integer, parameter :: max_terms = 200
        real(dp) :: sum_value
        real(dp) :: term

        term = 1.0_dp
        sum_value = real(n, dp)
        do k = 1, max_terms
            term = term * (-0.5_dp * x * x) / &
                (real(k, dp) * (2.0_dp * real(n + k, dp) + 1.0_dp))
            sum_value = sum_value + real(n + 2 * k, dp) * term
            if (abs(real(n + 2 * k, dp) * term) <= &
                    0.25_dp * epsilon(1.0_dp) * abs(sum_value)) exit
        end do
        log_value = spherical_series_log_prefactor(n, x) - log(x) + log(sum_value)
    end function log_spherical_jn_derivative_series

    pure function log_spherical_in_positive(n, x) result(log_value)
        integer, intent(in) :: n !! nonnegative order
        real(dp), intent(in) :: x !! positive finite argument
        real(dp) :: log_value

        integer :: k
        integer, parameter :: max_terms = 200
        real(dp) :: sum_value
        real(dp) :: term

        if (x < 1.0_dp) then
            term = 1.0_dp
            sum_value = 1.0_dp
            do k = 1, max_terms
                term = term * (0.5_dp * x * x) / &
                    (real(k, dp) * (2.0_dp * real(n + k, dp) + 1.0_dp))
                sum_value = sum_value + term
                if (abs(term) <= 0.25_dp * epsilon(1.0_dp) * sum_value) exit
            end do
            log_value = spherical_series_log_prefactor(n, x) + log(sum_value)
        else
            log_value = log_besseli_scaled(real(n, dp) + 0.5_dp, x) + x + &
                0.5_dp * (log(scifort_pi) - log(2.0_dp) - log(x))
        end if
    end function log_spherical_in_positive

    pure function spherical_series_log_prefactor(n, x) result(log_value)
        integer, intent(in) :: n !! nonnegative series order
        real(dp), intent(in) :: x !! positive series argument
        real(dp) :: log_value
        real(dp) :: rn

        rn = real(n, dp)
        log_value = rn * log(x) + rn * log(2.0_dp) + log_gamma(rn + 1.0_dp) - &
            log_gamma(2.0_dp * rn + 2.0_dp)
    end function spherical_series_log_prefactor

    pure function log_spherical_in_derivative_series(n, x) result(log_value)
        integer, intent(in) :: n !! positive modified spherical-Bessel order
        real(dp), intent(in) :: x !! argument with 0 < x < 1
        real(dp) :: log_value

        integer :: k
        integer, parameter :: max_terms = 200
        real(dp) :: sum_value
        real(dp) :: term

        term = 1.0_dp
        sum_value = real(n, dp)
        do k = 1, max_terms
            term = term * (0.5_dp * x * x) / &
                (real(k, dp) * (2.0_dp * real(n + k, dp) + 1.0_dp))
            sum_value = sum_value + real(n + 2 * k, dp) * term
            if (abs(real(n + 2 * k, dp) * term) <= &
                    0.25_dp * epsilon(1.0_dp) * sum_value) exit
        end do
        log_value = spherical_series_log_prefactor(n, x) - log(x) + log(sum_value)
    end function log_spherical_in_derivative_series

    pure function log_spherical_in_derivative_positive(n, x) result(log_value)
        integer, intent(in) :: n !! nonnegative modified spherical-Bessel order
        real(dp), intent(in) :: x !! positive finite argument
        real(dp) :: log_value

        real(dp) :: log_a
        real(dp) :: log_b
        integer :: sign_value

        if (n == 0) then
            log_value = log_spherical_in_positive(1, x)
            return
        end if
        if (x < 1.0_dp) then
            log_value = log_spherical_in_derivative_series(n, x)
            return
        end if
        log_a = log_spherical_in_positive(n - 1, x)
        log_b = log(real(n + 1, dp)) - log(x) + log_spherical_in_positive(n, x)
        call signed_log_add(log_a, 1, log_b, -1, log_value, sign_value)
        if (sign_value <= 0) log_value = quiet_nan(x)
    end function log_spherical_in_derivative_positive

    pure function log_spherical_kn_positive(n, x) result(log_value)
        integer, intent(in) :: n !! nonnegative modified spherical-Bessel order
        real(dp), intent(in) :: x !! positive finite argument
        real(dp) :: log_value

        integer :: k
        integer :: sign_value
        real(dp) :: log_a
        real(dp) :: log_sum
        real(dp) :: log_term
        real(dp) :: new_log
        integer :: new_sign

        log_a = 0.0_dp
        log_sum = negative_infinity(x)
        sign_value = 0
        do k = 0, n
            log_term = log_a - real(k + 1, dp) * log(x)
            call signed_log_add(log_sum, sign_value, log_term, 1, new_log, new_sign)
            log_sum = new_log
            sign_value = new_sign
            if (k < n) then
                log_a = log_a + log(real(n + k + 1, dp)) + &
                    log(real(n - k, dp)) - log(2.0_dp * real(k + 1, dp))
            end if
        end do
        log_value = log(0.5_dp * scifort_pi) - x + log_sum
    end function log_spherical_kn_positive

    pure function log_spherical_kn_derivative_positive(n, x) result(log_value)
        integer, intent(in) :: n !! nonnegative modified spherical-Bessel order
        real(dp), intent(in) :: x !! positive finite argument
        real(dp) :: log_value

        real(dp) :: log_a
        real(dp) :: log_b
        integer :: sign_value

        if (n == 0) then
            log_value = log_spherical_kn_positive(1, x)
            return
        end if
        log_a = log_spherical_kn_positive(n - 1, x)
        log_b = log(real(n + 1, dp)) - log(x) + log_spherical_kn_positive(n, x)
        call signed_log_add(log_a, 1, log_b, 1, log_value, sign_value)
    end function log_spherical_kn_derivative_positive

    pure integer function parity_sign(n) result(value)
        integer, intent(in) :: n !! integer exponent of -1

        if (modulo(n, 2) == 0) then
            value = 1
        else
            value = -1
        end if
    end function parity_sign

    pure subroutine signed_log_add(log_a, sign_a, log_b, sign_b, log_value, sign_value)
        real(dp), intent(in) :: log_a !! logarithm of first absolute term
        integer, intent(in) :: sign_a !! first term sign: -1, 0, or 1
        real(dp), intent(in) :: log_b !! logarithm of second absolute term
        integer, intent(in) :: sign_b !! second term sign: -1, 0, or 1
        real(dp), intent(out) :: log_value !! logarithm of absolute sum
        integer, intent(out) :: sign_value !! sum sign: -1, 0, or 1

        real(dp) :: hi
        real(dp) :: lo
        integer :: hi_sign

        if (sign_a == 0) then
            log_value = log_b
            sign_value = sign_b
            return
        else if (sign_b == 0) then
            log_value = log_a
            sign_value = sign_a
            return
        end if

        if (sign_a == sign_b) then
            hi = max(log_a, log_b)
            lo = min(log_a, log_b)
            log_value = hi + log1p_safe(exp(lo - hi))
            sign_value = sign_a
            return
        end if

        if (log_a == log_b) then
            log_value = negative_infinity(log_a)
            sign_value = 0
            return
        end if
        if (log_a > log_b) then
            hi = log_a
            lo = log_b
            hi_sign = sign_a
        else
            hi = log_b
            lo = log_a
            hi_sign = sign_b
        end if
        log_value = hi + log1p_safe(-exp(lo - hi))
        sign_value = hi_sign
    end subroutine signed_log_add

    pure function spherical_j0_base(x) result(y)
        real(dp), intent(in) :: x !! positive finite argument
        real(dp) :: y
        real(dp) :: z

        if (x < 0.1_dp) then
            z = x * x
            y = 1.0_dp + z * (-1.0_dp / 6.0_dp + z * (1.0_dp / 120.0_dp + &
                z * (-1.0_dp / 5040.0_dp + z / 362880.0_dp)))
        else
            y = sin(x) / x
        end if
    end function spherical_j0_base

    pure function spherical_j1_base(x) result(y)
        real(dp), intent(in) :: x !! positive finite argument
        real(dp) :: y
        real(dp) :: z

        if (x < 0.1_dp) then
            z = x * x
            y = x * (1.0_dp / 3.0_dp + z * (-1.0_dp / 30.0_dp + &
                z * (1.0_dp / 840.0_dp - z / 45360.0_dp)))
        else
            y = (sin(x) / x - cos(x)) / x
        end if
    end function spherical_j1_base

    pure function spherical_y0_base(x) result(y)
        real(dp), intent(in) :: x !! positive finite argument
        real(dp) :: y

        if (x < 1.0_dp / huge(1.0_dp)) then
            y = negative_infinity(x)
        else
            y = -cos(x) / x
        end if
    end function spherical_y0_base

    pure function spherical_y1_base(x) result(y)
        real(dp), intent(in) :: x !! positive finite argument
        real(dp) :: y
        real(dp) :: inv_x
        real(dp) :: inner

        inv_x = 1.0_dp / x
        if (inv_x > sqrt(huge(1.0_dp))) then
            y = negative_infinity(x)
            return
        end if
        inner = cos(x) * inv_x + sin(x)
        if (inv_x > 1.0_dp) then
            if (abs(inner) > huge(1.0_dp) / inv_x) then
                y = -signed_infinity(inner, x)
                return
            end if
        end if
        y = -inv_x * inner
    end function spherical_y1_base

    pure function positive_exp_from_log(log_value, kind_arg) result(y)
        real(dp), intent(in) :: log_value !! logarithm of a nonnegative result
        real(dp), intent(in) :: kind_arg !! real-kind selector for IEEE special values
        real(dp) :: y

        if (ieee_is_nan(log_value)) then
            y = log_value
        else if (log_value > log(huge(1.0_dp))) then
            y = positive_infinity(kind_arg)
        else
            y = exp(log_value)
        end if
    end function positive_exp_from_log

    pure function signed_exp_from_log(log_value, sign_value, kind_arg) result(y)
        real(dp), intent(in) :: log_value !! logarithm of an absolute result
        integer, intent(in) :: sign_value !! result sign: -1, 0, or 1
        real(dp), intent(in) :: kind_arg !! real-kind selector for IEEE special values
        real(dp) :: y

        if (sign_value == 0) then
            y = 0.0_dp
        else
            y = positive_exp_from_log(log_value, kind_arg)
            if (sign_value < 0) y = -y
        end if
    end function signed_exp_from_log

    pure function signed_infinity(sign_source, kind_arg) result(y)
        real(dp), intent(in) :: sign_source !! nonzero value whose sign is copied
        real(dp), intent(in) :: kind_arg !! real-kind selector for IEEE infinity
        real(dp) :: y

        if (sign_source < 0.0_dp) then
            y = negative_infinity(kind_arg)
        else
            y = positive_infinity(kind_arg)
        end if
    end function signed_infinity

end module scifort_spherical_bessel
