! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Regularized incomplete beta function I_x(a, b), its complement, and their
! inverses with respect to x.
!
! Independent implementation from the following public specifications:
! - DLMF 8.17.4: symmetry I_x(a, b) = 1 - I_(1-x)(b, a);
! - DLMF 8.17.22 and 8.17.23: continued fraction for I_x(a, b), evaluated
!   with the modified Lentz method (Lentz 1976; Thompson and Barnett 1986);
! - DLMF 5.11.1 through scifort_log_gamma for the prefactor with large
!   parameters;
! - for one large and one small parameter, an expansion in incomplete gamma
!   functions derived here from the integral substitution t = exp(-v) and
!   the Bernoulli generating function (DLMF 24.2.1); see gamma_expansion.
!
! Internal procedures take both x and y = 1 - x so that callers that know
! y more accurately than 1 - x (for example the Student t and F
! distributions) do not lose accuracy near x = 1.
!
! This module is an implementation layer. The supported entry points are
! re-exported by scifort_special.

module scifort_incomplete_beta
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_constants, only : scifort_log_sqrt_two_pi
    use scifort_incomplete_gamma, only : log_gamma_kernel, log_gammainc, &
        log_gammaincc
    use scifort_kinds, only : dp
    use scifort_log_gamma, only : log1p_minus_x, log_beta, log_gamma_difference, &
        log_gamma_one_plus, log_gamma_ratio_scaled, stirling_remainder, &
        stirling_threshold
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, &
        positive_infinity, quiet_nan
    implicit none
    private

    public :: beta_inverse_xy
    public :: betainc
    public :: betaincc
    public :: betainccinv
    public :: betaincinv
    public :: incomplete_beta_xy
    public :: log_beta_kernel

    integer, parameter :: max_inverse_iterations = 200

    ! B(2n) / (2n (2n)!) for n = 1, ..., 30: the coefficients of v**(2n) in
    ! log(sinh(v/2) / (v/2)) = log((1 - exp(-v)) / v) + v/2
    !                        = sum B(2n) v**(2n) / (2n (2n)!),
    ! which follows from DLMF 24.2.1 because the derivative of
    ! log((1 - exp(-v)) / v) is 1 / (exp(v) - 1) - 1 / v. Generated with
    ! mpmath 1.3.0 at 40 digits and checked against the left side to 1e-40 at
    ! v = 1.
    real(dp), parameter :: log_ratio_coefficients(30) = [ &
        0.04166666666666666666666667_dp, &
        -0.0003472222222222222222222222_dp, &
        0.000005511463844797178130511464_dp, &
        -0.0000001033399470899470899470899_dp, &
        2.087675698786809897921009e-9_dp, &
        -4.403491782239577654039735e-11_dp, &
        9.558954664774770594876415e-13_dp, &
        -2.118550185201614291768872e-14_dp, &
        4.77003447570991364674217e-16_dp, &
        -1.087434349279030936520758e-17_dp, &
        2.504092194709195234183024e-19_dp, &
        -5.814360285755218058628203e-21_dp, &
        1.359502707549795181420473e-22_dp, &
        -3.197684795370552446572361e-24_dp, &
        7.559841507792276867703169e-26_dp, &
        -1.795247084022563264144963e-27_dp, &
        4.279919045926073240783132e-29_dp, &
        -1.023887483518141716161606e-30_dp, &
        2.457035330814485440008567e-32_dp, &
        -5.912556039251574836399088e-34_dp, &
        1.426350419638603405871408e-35_dp, &
        -3.448761101064296049342757e-37_dp, &
        8.355995924900409200498794e-39_dp, &
        -2.028406890134727299014088e-40_dp, &
        4.932494088401361914212801e-42_dp, &
        -1.201360911888604556374761e-43_dp, &
        2.930375971230539684724225e-45_dp, &
        -7.157631582051671372375022e-47_dp, &
        1.750530319925785570970813e-48_dp, &
        -4.286340263736452916541366e-50_dp]

    ! The gamma expansion is used when the large parameter is at least
    ! expansion_large and at least expansion_ratio times the small one, and
    ! w = -log of the coordinate of the large parameter satisfies w <= 1 and
    ! small * w**2 <= expansion_spread (the expansion coefficients behave like
    ! (small w**2 / 24)**n / n!).
    real(dp), parameter :: expansion_large = 20.0_dp
    real(dp), parameter :: expansion_ratio = 8.0_dp
    real(dp), parameter :: expansion_spread = 4.0_dp
    real(dp), parameter :: inverse_e = 0.367879441171442321595523770161460867_dp

contains

    ! Regularized incomplete beta function I_x(a, b).
    pure elemental function betainc(a, b, x) result(p)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        real(dp), intent(in) :: x
        real(dp) :: p

        real(dp) :: logp
        real(dp) :: logq
        real(dp) :: q

        if (.not. valid_arguments(a, b, x)) then
            p = quiet_nan(x)
        else if (x <= 0.0_dp) then
            p = 0.0_dp
        else if (x >= 1.0_dp) then
            p = 1.0_dp
        else
            call incomplete_beta_xy(a, b, x, 1.0_dp - x, p, q, logp, logq)
        end if
    end function betainc

    ! Complement 1 - I_x(a, b), computed without forming the difference.
    pure elemental function betaincc(a, b, x) result(q)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        real(dp), intent(in) :: x
        real(dp) :: q

        real(dp) :: logp
        real(dp) :: logq
        real(dp) :: p

        if (.not. valid_arguments(a, b, x)) then
            q = quiet_nan(x)
        else if (x <= 0.0_dp) then
            q = 1.0_dp
        else if (x >= 1.0_dp) then
            q = 0.0_dp
        else
            call incomplete_beta_xy(a, b, x, 1.0_dp - x, p, q, logp, logq)
        end if
    end function betaincc

    ! Inverse of I_x(a, b) with respect to x.
    pure elemental function betaincinv(a, b, p) result(x)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        real(dp), intent(in) :: p
        real(dp) :: x

        real(dp) :: y

        if (.not. valid_arguments(a, b, p)) then
            x = quiet_nan(p)
        else if (p <= 0.0_dp) then
            x = 0.0_dp
        else if (p >= 1.0_dp) then
            x = 1.0_dp
        else
            call beta_inverse_xy(a, b, p, 1.0_dp - p, x, y)
        end if
    end function betaincinv

    ! Inverse of 1 - I_x(a, b) with respect to x.
    pure elemental function betainccinv(a, b, q) result(x)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        real(dp), intent(in) :: q
        real(dp) :: x

        real(dp) :: y

        if (.not. valid_arguments(a, b, q)) then
            x = quiet_nan(q)
        else if (q <= 0.0_dp) then
            x = 1.0_dp
        else if (q >= 1.0_dp) then
            x = 0.0_dp
        else
            call beta_inverse_xy(a, b, 1.0_dp - q, q, x, y)
        end if
    end function betainccinv

    ! I_x(a, b), its complement, and their logarithms for finite a, b > 0 and
    ! 0 < x < 1 with y = 1 - x supplied by the caller. When one parameter is
    ! much larger than the other and the coordinate is in the transition
    ! region, gamma_expansion is used. Otherwise the tail on the side where
    ! the continued fraction converges quickly (DLMF 8.17.22) is computed
    ! first. That tail exceeds one half only when its first parameter is
    ! small; the complement is then computed directly by
    ! small_parameter_complement so that it is not formed by cancellation.
    pure elemental subroutine incomplete_beta_xy(a, b, x, y, p, q, logp, logq)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y
        real(dp), intent(out) :: p
        real(dp), intent(out) :: q
        real(dp), intent(out) :: logp
        real(dp), intent(out) :: logq

        if (use_gamma_expansion(a, b, x, y)) then
            call gamma_expansion(a, b, x, y, p, q, logp, logq)
        else if (use_gamma_expansion(b, a, y, x)) then
            call gamma_expansion(b, a, y, x, q, p, logq, logp)
        else if (x * (a + b + 2.0_dp) < a + 1.0_dp) then
            call lower_tail(a, b, x, y, p, logp)
            if (p > 0.5_dp .and. a < 1.0_dp) then
                q = small_parameter_complement(a, b, x)
                p = 1.0_dp - q
                logq = log(q)
                logp = log1p_safe(-q)
            else
                q = 1.0_dp - p
                logq = log1p_safe(-p)
            end if
        else
            call lower_tail(b, a, y, x, q, logq)
            if (q > 0.5_dp .and. b < 1.0_dp) then
                p = small_parameter_complement(b, a, y)
                q = 1.0_dp - p
                logp = log(p)
                logq = log1p_safe(-p)
            else
                p = 1.0_dp - q
                logp = log1p_safe(-q)
            end if
        end if
    end subroutine incomplete_beta_xy

    ! I_x(a, b) = x**a y**b / (a B(a, b)) / F with F the continued fraction
    ! 1 + d(1) / (1 + d(2) / (1 + ...)) of DLMF 8.17.22 and 8.17.23:
    ! d(2m) = m (b - m) x / ((a + 2m - 1) (a + 2m)),
    ! d(2m + 1) = -(a + m) (a + b + m) x / ((a + 2m) (a + 2m + 1)).
    pure elemental subroutine lower_tail(a, b, x, y, t, logt)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y
        real(dp), intent(out) :: t
        real(dp), intent(out) :: logt

        real(dp), parameter :: tiny_value = 1.0e-300_dp
        integer :: max_terms
        integer :: m
        integer :: n
        real(dp) :: c
        real(dp) :: d
        real(dp) :: delta
        real(dp) :: f
        real(dp) :: k
        real(dp) :: logk
        real(dp) :: numerator
        real(dp) :: rm

        max_terms = iteration_limit(a, b)
        f = 1.0_dp
        c = f
        d = 0.0_dp
        do n = 1, max_terms
            if (mod(n, 2) == 0) then
                m = n / 2
                rm = real(m, dp)
                numerator = rm * (b - rm) * x / ((a + 2.0_dp * rm - 1.0_dp) * &
                    (a + 2.0_dp * rm))
            else
                m = (n - 1) / 2
                rm = real(m, dp)
                numerator = -(a + rm) * (a + b + rm) * x / ((a + 2.0_dp * rm) * &
                    (a + 2.0_dp * rm + 1.0_dp))
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

        call beta_kernel(a, b, x, y, k, logk)
        t = k / (a * f)
        logt = logk - log(a) - log(f)
    end subroutine lower_tail

    ! 1 - I_z(s, t) for 0 < s < 1 and z below the continued-fraction switch
    ! point. From DLMF 8.17.7 with (s)_k / (s + 1)_k = s / (s + k),
    ! I_z(s, t) = exp(L) (1 + s T), where
    ! L = s log(z) + log(Gamma(s + t)) - log(Gamma(t)) - log(Gamma(1 + s)) and
    ! T = sum over k >= 1 of (1 - t)_k z**k / (k! (s + k)). The complement
    ! -expm1(L) - exp(L) s T keeps full relative accuracy as s -> 0.
    pure function small_parameter_complement(s, t, z) result(c)
        real(dp), intent(in) :: s
        real(dp), intent(in) :: t
        real(dp), intent(in) :: z
        real(dp) :: c

        integer :: k
        real(dp) :: big_l
        real(dp) :: power
        real(dp) :: series
        real(dp) :: term

        series = 0.0_dp
        power = 1.0_dp
        do k = 1, 1000
            power = power * (real(k, dp) - t) * z / real(k, dp)
            term = power / (s + real(k, dp))
            series = series + term
            if (abs(term) <= 0.25_dp * epsilon(1.0_dp) * abs(series)) exit
            if (power == 0.0_dp) exit
        end do

        big_l = s * log(z) + log_gamma_difference(t, s) - log_gamma_one_plus(s)
        c = -expm1_safe(big_l) - exp(big_l) * s * series
    end function small_parameter_complement

    pure elemental logical function use_gamma_expansion(big, small, x, y) result(use)
        real(dp), intent(in) :: big
        real(dp), intent(in) :: small
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y

        real(dp) :: w

        use = big >= expansion_large .and. big >= expansion_ratio * small .and. &
            x >= inverse_e
        if (use) then
            w = minus_log(x, y)
            use = small * w * w <= expansion_spread
        end if
    end function use_gamma_expansion

    ! -log(x) given y = 1 - x, using the complement when x > 1/2.
    pure elemental function minus_log(x, y) result(w)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y
        real(dp) :: w

        if (x > 0.5_dp) then
            w = -log1p_safe(-y)
        else
            w = -log(x)
        end if
    end function minus_log

    ! I_x(a, b) for large a, small b, and x >= exp(-1), with y = 1 - x.
    ! Substituting t = exp(-v) in the integral for B_x(a, b) gives
    ! B_x(a, b) = integral over v >= w of exp(-a v) (1 - exp(-v))**(b - 1) dv,
    ! w = -log(x). Writing 1 - exp(-v) = v exp(-v/2) sinh(v/2) / (v/2) and
    ! s = a + (b - 1)/2,
    ! B_x(a, b) = integral over v >= w of exp(-s v) v**(b - 1) g(v) dv,
    ! g(v) = (sinh(v/2) / (v/2))**(b - 1) = sum c_k v**k,
    ! an even power series with radius 2 pi whose coefficients are small.
    ! Integrating term by term,
    ! I_x(a, b) = R sum c_k Gamma(b + k, u) / (Gamma(b) s**k), u = s w,
    ! R = Gamma(a + b) / (Gamma(a) s**b). Over v >= w this is asymptotic in
    ! 1/s (the neglected part is of order exp(-2 pi s)); over 0 <= v <= w < 2 pi
    ! it converges and gives the complement with lower incomplete gamma
    ! functions. The smaller tail is summed and the other is its complement.
    pure elemental subroutine gamma_expansion(a, b, x, y, p, q, logp, logq)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y
        real(dp), intent(out) :: p
        real(dp), intent(out) :: q
        real(dp), intent(out) :: logp
        real(dp), intent(out) :: logq

        integer, parameter :: max_terms = 60
        integer :: j
        integer :: k
        real(dp) :: c(0:max_terms)
        real(dp) :: density_ratio
        real(dp) :: log_scale
        real(dp) :: logp0
        real(dp) :: logq0
        real(dp) :: logr
        real(dp) :: phi(max_terms)
        real(dp) :: ratio_k
        real(dp) :: series
        real(dp) :: shifted
        real(dp) :: term
        real(dp) :: u
        real(dp) :: w
        real(dp) :: w_power

        w = minus_log(x, y)
        shifted = a + 0.5_dp * (b - 1.0_dp)
        u = shifted * w

        ! Coefficients of g = exp(phi), phi(v) = (b - 1) log(sinh(v/2) / (v/2)),
        ! from k c_k = sum over j = 1..k of j phi_j c_(k-j).
        phi = 0.0_dp
        do j = 1, max_terms / 2
            phi(2 * j) = (b - 1.0_dp) * log_ratio_coefficients(j)
        end do
        c(0) = 1.0_dp
        do k = 1, max_terms
            c(k) = 0.0_dp
            do j = 1, k
                c(k) = c(k) + real(j, dp) * phi(j) * c(k - j)
            end do
            c(k) = c(k) / real(k, dp)
        end do

        ! log(Gamma(a + b) / (Gamma(a) a**b)) - b log(s / a).
        logr = log_gamma_ratio_scaled(a, b) - b * log1p_safe(0.5_dp * (b - 1.0_dp) / a)
        logq0 = log_gammaincc(b, u)
        logp0 = log_gammainc(b, u)

        if (logq0 <= logp0) then
            ! Upper series with J_k = Gamma(b + k, u) / (Gamma(b) s**k) scaled by
            ! J_0 = Q(b, u), from Gamma(s + 1, u) = s Gamma(s, u) + u**s exp(-u)
            ! (DLMF 8.8.2), which is stable upward.
            density_ratio = exp(log_gamma_kernel(b, u) + log(b) - logq0)
            ratio_k = 1.0_dp
            series = 1.0_dp
            w_power = 1.0_dp
            do k = 0, max_terms - 1
                ratio_k = ((b + real(k, dp)) / shifted) * ratio_k + &
                    density_ratio * w_power / shifted
                w_power = w_power * w
                ! The odd coefficients of the even series g vanish; only the
                ! nonzero terms are used to judge convergence.
                if (c(k + 1) == 0.0_dp) cycle
                term = c(k + 1) * ratio_k
                series = series + term
                if (abs(term) <= 0.25_dp * epsilon(1.0_dp) * abs(series)) exit
            end do
            logp = logr + logq0 + log(series)
            p = exp(logp)
            q = 1.0_dp - p
            logq = log1p_safe(-p)
        else
            ! Lower series with L_k = gamma(b + k, u) / (Gamma(b) s**k), each
            ! formed from P(b + k, u) because the upward recurrence for the
            ! lower function is unstable.
            series = 1.0_dp
            log_scale = 0.0_dp
            do k = 1, max_terms
                log_scale = log_scale + log(b + real(k - 1, dp)) - log(shifted)
                if (c(k) == 0.0_dp) cycle
                term = c(k) * exp(log_gammainc(b + real(k, dp), u) - logp0 + log_scale)
                series = series + term
                if (abs(term) <= 0.25_dp * epsilon(1.0_dp) * abs(series)) exit
            end do
            logq = logr + logp0 + log(series)
            q = exp(logq)
            p = 1.0_dp - q
            logp = log1p_safe(-q)
        end if
    end subroutine gamma_expansion

    ! Term limit for the continued fraction. Near the mean the number of
    ! terms needed grows like sqrt(a + b).
    pure function iteration_limit(a, b) result(n)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        integer :: n

        real(dp) :: estimate

        estimate = 1000.0_dp + 40.0_dp * sqrt(a + b)
        if (estimate > real(huge(n) / 4, dp)) then
            n = huge(n) / 4
        else
            n = int(estimate)
        end if
    end function iteration_limit

    ! log(x**a y**b / B(a, b)) for finite a, b > 0 and x, y in (0, 1] with
    ! x + y = 1. With n = a + b, p0 = a / n, and q0 = b / n, the Stirling
    ! series gives, when both parameters are large,
    ! a (log(x/p0) - (x/p0 - 1)) + b (log(y/q0) - (y/q0 - 1))
    !     + log(sqrt(a b / (2 pi n))) - S(a) - S(b) + S(n),
    ! in which no large terms cancel. When one parameter is large the
    ! n log(n) terms of log(B) cancel analytically (see log_beta).
    pure elemental function log_beta_kernel(a, b, x, y) result(logk)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y
        real(dp) :: logk

        real(dp) :: n

        if (x <= 0.0_dp .or. y <= 0.0_dp) then
            logk = negative_infinity(x)
            return
        end if

        n = a + b
        if (min(a, b) >= stirling_threshold) then
            logk = a * relative_log_term(x, y, a / n, b / n) + &
                b * relative_log_term(y, x, b / n, a / n) + &
                0.5_dp * (log(a) + log(b) - log(n)) - scifort_log_sqrt_two_pi - &
                stirling_remainder(a) - stirling_remainder(b) + stirling_remainder(n)
        else if (a >= stirling_threshold) then
            logk = b * log_pair(y, x) + a * log_pair(x, y) + b * log(n) - log_gamma(b) - &
                (a - 0.5_dp) * log1p_safe(-b / n) - b - stirling_remainder(a) + &
                stirling_remainder(n)
        else if (b >= stirling_threshold) then
            logk = a * log_pair(x, y) + b * log_pair(y, x) + a * log(n) - log_gamma(a) - &
                (b - 0.5_dp) * log1p_safe(-a / n) - a - stirling_remainder(b) + &
                stirling_remainder(n)
        else
            logk = a * log_pair(x, y) + b * log_pair(y, x) - log_beta(a, b)
        end if
    end function log_beta_kernel

    ! x**a y**b / B(a, b) and its logarithm. For moderate parameters and
    ! when no intermediate can overflow or underflow, the kernel is formed
    ! as a product so that its relative error does not grow with its
    ! logarithm.
    pure elemental subroutine beta_kernel(a, b, x, y, k, logk)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y
        real(dp), intent(out) :: k
        real(dp), intent(out) :: logk

        real(dp), parameter :: safe_exponent = 690.0_dp

        logk = log_beta_kernel(a, b, x, y)
        if (max(a, b) < stirling_threshold .and. min(a, b) >= 1.0e-100_dp .and. &
            logk >= -safe_exponent .and. abs(a * log(x)) <= safe_exponent .and. &
            abs(b * log(y)) <= safe_exponent) then
            k = x**a * y**b * (gamma(a + b) / gamma(a) / gamma(b))
            logk = log(k)
        else
            k = exp(logk)
        end if
    end subroutine beta_kernel

    ! log(z) for z in (0, 1] given its complement w = 1 - z. Near z = 1 the
    ! complement carries the information, so log1p(-w) is used.
    pure elemental function log_pair(z, w) result(y)
        real(dp), intent(in) :: z
        real(dp), intent(in) :: w
        real(dp) :: y

        if (z <= 0.5_dp) then
            y = log(z)
        else
            y = log1p_safe(-w)
        end if
    end function log_pair

    ! log(z / z0) - (z / z0 - 1) without cancellation near z = z0, where
    ! w = 1 - z and w0 = 1 - z0. When z > 1/2 the difference z - z0 is
    ! formed as w0 - w from the complements.
    pure elemental function relative_log_term(z, w, z0, w0) result(y)
        real(dp), intent(in) :: z
        real(dp), intent(in) :: w
        real(dp), intent(in) :: z0
        real(dp), intent(in) :: w0
        real(dp) :: y

        real(dp) :: t

        if (z <= 0.5_dp) then
            t = (z - z0) / z0
        else
            t = (w0 - w) / z0
        end if
        if (abs(t) <= 0.25_dp) then
            y = log1p_minus_x(t)
        else
            y = log(z / z0) - t
        end if
    end function relative_log_term

    pure elemental logical function valid_arguments(a, b, x) result(valid)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        real(dp), intent(in) :: x

        valid = ieee_is_finite(a) .and. a > 0.0_dp .and. &
            ieee_is_finite(b) .and. b > 0.0_dp .and. &
            x >= 0.0_dp .and. x <= 1.0_dp
    end function valid_arguments

    ! Solve I_x(a, b) = p, equivalently 1 - I_x(a, b) = q, for 0 < p < 1 and
    ! q = 1 - p, returning x and y = 1 - x. The root is located relative to
    ! x = 1/2 first; a root above one half is found as the solution y of
    ! I_y(b, a) = q (DLMF 8.17.4), so that both x and y keep full relative
    ! accuracy.
    pure elemental subroutine beta_inverse_xy(a, b, p, q, x, y)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        real(dp), intent(in) :: p
        real(dp), intent(in) :: q
        real(dp), intent(out) :: x
        real(dp), intent(out) :: y

        logical :: root_below_half
        real(dp) :: logp_half
        real(dp) :: logq_half
        real(dp) :: p_half
        real(dp) :: q_half

        call incomplete_beta_xy(a, b, 0.5_dp, 0.5_dp, p_half, q_half, logp_half, &
            logq_half)
        if (p <= q) then
            root_below_half = p <= p_half
        else
            root_below_half = q >= q_half
        end if

        if (root_below_half) then
            x = solve_below_half(a, b, p, q)
            y = 1.0_dp - x
        else
            y = solve_below_half(b, a, q, p)
            x = 1.0_dp - y
        end if
    end subroutine beta_inverse_xy

    ! Solve I_x(a, b) = p with q = 1 - p for a root known to lie in
    ! (0, 1/2]. The residual is logarithmic in the smaller of p and q and is
    ! increasing in x. Newton steps in log(x) are safeguarded by a bracket
    ! with geometric bisection as the fallback.
    pure function solve_below_half(a, b, p, q) result(x)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        real(dp), intent(in) :: p
        real(dp), intent(in) :: q
        real(dp) :: x

        real(dp), parameter :: max_log_step = 50.0_dp
        integer :: iteration
        logical :: use_lower
        real(dp) :: hi
        real(dp) :: lo
        real(dp) :: logk
        real(dp) :: logp_x
        real(dp) :: logq_x
        real(dp) :: log_target
        real(dp) :: p_x
        real(dp) :: q_x
        real(dp) :: residual
        real(dp) :: slope
        real(dp) :: step
        real(dp) :: x_new
        real(dp) :: y

        use_lower = p <= q
        if (use_lower) then
            log_target = log(p)
        else
            log_target = log(q)
        end if

        ! Start from the small-x behavior I_x(a, b) ~ x**a / (a B(a, b)).
        x = exp(max(log(tiny(1.0_dp)), (log(p) + log(a) + log_beta(a, b)) / a))
        x = min(x, 0.5_dp)
        lo = 0.0_dp
        hi = 0.5_dp

        do iteration = 1, max_inverse_iterations
            y = 1.0_dp - x
            call incomplete_beta_xy(a, b, x, y, p_x, q_x, logp_x, logq_x)
            logk = log_beta_kernel(a, b, x, y)
            ! d/d(log x) of log(I) is x f(x) / I with x f(x) = kernel / y.
            if (use_lower) then
                residual = logp_x - log_target
                slope = exp(logk - log(y) - logp_x)
            else
                residual = log_target - logq_x
                slope = exp(logk - log(y) - logq_x)
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
                step = sign(max_log_step, -residual)
            end if
            x_new = x * exp(step)

            if (.not. (x_new > lo .and. x_new < hi)) then
                if (lo > 0.0_dp) then
                    if (hi > 4.0_dp * lo) then
                        x_new = sqrt(lo) * sqrt(hi)
                    else
                        x_new = lo + 0.5_dp * (hi - lo)
                    end if
                else
                    x_new = hi * exp(-max_log_step)
                    if (x_new <= 0.0_dp) then
                        ! The root is below hi < 1e-300 and rounds to zero.
                        x = 0.0_dp
                        return
                    end if
                end if
            end if

            if (abs(x_new - x) <= 2.0_dp * epsilon(1.0_dp) * x) then
                x = x_new
                exit
            end if
            if (hi - lo <= 2.0_dp * epsilon(1.0_dp) * hi) then
                x = x_new
                exit
            end if
            x = x_new
        end do
    end function solve_below_half

end module scifort_incomplete_beta
