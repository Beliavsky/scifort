! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! SciPy-compatible legacy distribution helpers built from SciFort's regularized
! incomplete beta/gamma and normal-quantile kernels.  No probability mass or
! density is summed directly.

module scifort_special_distribution_helpers
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_incomplete_beta, only : beta_inverse_xy, betainc, betaincc, &
        betainccinv, betaincinv
    use scifort_incomplete_gamma, only : gammainc, gammaincc, gammainccinv, &
        gammaincinv
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan
    use scifort_special_elementary, only : ndtr, ndtri
    implicit none
    private

    public :: bdtr
    public :: bdtrc
    public :: bdtri
    public :: bdtrik
    public :: bdtrin
    public :: chdtr
    public :: chdtrc
    public :: btdtria
    public :: btdtrib
    public :: chdtri
    public :: chdtriv
    public :: fdtr
    public :: fdtrc
    public :: fdtri
    public :: fdtridfd
    public :: gdtr
    public :: gdtrc
    public :: gdtria
    public :: gdtrib
    public :: gdtrix
    public :: nbdtr
    public :: nbdtrc
    public :: nbdtri
    public :: nbdtrik
    public :: nbdtrin
    public :: pdtr
    public :: pdtrc
    public :: pdtri
    public :: pdtrik
    public :: nrdtrimn
    public :: nrdtrisd
    public :: stdtr
    public :: stdtridf
    public :: stdtrit

    integer, parameter :: root_bdtrin = 1
    integer, parameter :: root_btdtria = 2
    integer, parameter :: root_btdtrib = 3
    integer, parameter :: root_chdtriv = 4
    integer, parameter :: root_gdtrib = 5
    integer, parameter :: root_nbdtrik = 6
    integer, parameter :: root_nbdtrin = 7
    integer, parameter :: root_pdtrik = 8
    integer, parameter :: root_stdtridf = 9

contains

    ! Binomial CDF P(X <= floor(k)); n is the number of Bernoulli trials.
    pure elemental function bdtr(k, n, p) result(value)
        real(dp), intent(in) :: k !! success-count threshold, finite and >= 0
        integer, intent(in) :: n !! number of trials, >= 0
        real(dp), intent(in) :: p !! event probability in [0, 1]
        real(dp) :: value

        real(dp) :: kf

        if (.not. valid_probability(p) .or. n < 0 .or. ieee_is_nan(k) .or. &
                .not. ieee_is_finite(k) .or. k < 0.0_dp) then
            value = quiet_nan(p)
            return
        end if
        kf = floor(k)
        if (kf > real(n, dp)) then
            value = quiet_nan(p)
        else if (kf == real(n, dp)) then
            value = 1.0_dp
        else
            value = betaincc(kf + 1.0_dp, real(n, dp) - kf, p)
        end if
    end function bdtr

    ! Binomial survival probability P(X > floor(k)).
    pure elemental function bdtrc(k, n, p) result(value)
        real(dp), intent(in) :: k !! success-count threshold, finite and >= 0
        integer, intent(in) :: n !! number of trials, >= 0
        real(dp), intent(in) :: p !! event probability in [0, 1]
        real(dp) :: value

        real(dp) :: kf

        if (.not. valid_probability(p) .or. n < 0 .or. ieee_is_nan(k) .or. &
                .not. ieee_is_finite(k) .or. k < 0.0_dp) then
            value = quiet_nan(p)
            return
        end if
        kf = floor(k)
        if (kf > real(n, dp)) then
            value = quiet_nan(p)
        else if (kf == real(n, dp)) then
            value = 0.0_dp
        else
            value = betainc(kf + 1.0_dp, real(n, dp) - kf, p)
        end if
    end function bdtrc

    ! Inverse of bdtr with respect to the Bernoulli probability p.
    pure elemental function bdtri(k, n, y) result(p)
        real(dp), intent(in) :: k !! success-count threshold, finite and >= 0
        integer, intent(in) :: n !! number of trials, > 0
        real(dp), intent(in) :: y !! target binomial CDF in [0, 1]
        real(dp) :: p

        real(dp) :: kf

        if (.not. valid_probability(y) .or. n <= 0 .or. ieee_is_nan(k) .or. &
                .not. ieee_is_finite(k) .or. k < 0.0_dp) then
            p = quiet_nan(y)
            return
        end if
        kf = floor(k)
        if (kf >= real(n, dp)) then
            p = quiet_nan(y)
        else
            p = betainccinv(kf + 1.0_dp, real(n, dp) - kf, y)
        end if
    end function bdtri

    ! Negative-binomial CDF for failures k before the nth success.
    pure elemental function nbdtr(k, n, p) result(value)
        integer, intent(in) :: k !! maximum number of failures, >= 0
        integer, intent(in) :: n !! target number of successes, > 0
        real(dp), intent(in) :: p !! Bernoulli success probability in [0, 1]
        real(dp) :: value

        if (k < 0 .or. n <= 0 .or. .not. valid_probability(p)) then
            value = quiet_nan(p)
        else
            value = betainc(real(n, dp), real(k, dp) + 1.0_dp, p)
        end if
    end function nbdtr

    ! Negative-binomial survival probability beyond k failures.
    pure elemental function nbdtrc(k, n, p) result(value)
        integer, intent(in) :: k !! maximum number of failures, >= 0
        integer, intent(in) :: n !! target number of successes, > 0
        real(dp), intent(in) :: p !! Bernoulli success probability in [0, 1]
        real(dp) :: value

        if (k < 0 .or. n <= 0 .or. .not. valid_probability(p)) then
            value = quiet_nan(p)
        else
            value = betaincc(real(n, dp), real(k, dp) + 1.0_dp, p)
        end if
    end function nbdtrc

    ! Inverse of nbdtr with respect to the Bernoulli success probability p.
    pure elemental function nbdtri(k, n, y) result(p)
        integer, intent(in) :: k !! maximum number of failures, >= 0
        integer, intent(in) :: n !! target number of successes, > 0
        real(dp), intent(in) :: y !! target negative-binomial CDF in [0, 1]
        real(dp) :: p

        if (k < 0 .or. n <= 0 .or. .not. valid_probability(y)) then
            p = quiet_nan(y)
        else
            p = betaincinv(real(n, dp), real(k, dp) + 1.0_dp, y)
        end if
    end function nbdtri

    ! Poisson CDF P(X <= floor(k)) for mean m.
    pure elemental function pdtr(k, m) result(value)
        real(dp), intent(in) :: k !! occurrence-count threshold, >= 0; +infinity allowed
        real(dp), intent(in) :: m !! Poisson mean, >= 0; +infinity allowed
        real(dp) :: value

        real(dp) :: shape

        if (ieee_is_nan(k) .or. ieee_is_nan(m) .or. k < 0.0_dp .or. m < 0.0_dp) then
            value = quiet_nan(m)
        else if (.not. ieee_is_finite(k)) then
            if (.not. ieee_is_finite(m)) then
                value = quiet_nan(m)
            else
                value = 1.0_dp
            end if
        else
            shape = floor(k) + 1.0_dp
            value = gammaincc(shape, m)
        end if
    end function pdtr

    ! Poisson survival probability P(X > floor(k)).
    pure elemental function pdtrc(k, m) result(value)
        real(dp), intent(in) :: k !! occurrence-count threshold, >= 0; +infinity allowed
        real(dp), intent(in) :: m !! Poisson mean, >= 0; +infinity allowed
        real(dp) :: value

        real(dp) :: shape

        if (ieee_is_nan(k) .or. ieee_is_nan(m) .or. k < 0.0_dp .or. m < 0.0_dp) then
            value = quiet_nan(m)
        else if (.not. ieee_is_finite(k)) then
            if (.not. ieee_is_finite(m)) then
                value = quiet_nan(m)
            else
                value = 0.0_dp
            end if
        else
            shape = floor(k) + 1.0_dp
            value = gammainc(shape, m)
        end if
    end function pdtrc

    ! Inverse of pdtr with respect to the Poisson mean m.
    pure elemental function pdtri(k, y) result(m)
        real(dp), intent(in) :: k !! occurrence-count threshold, finite and >= 0
        real(dp), intent(in) :: y !! target Poisson CDF in [0, 1]
        real(dp) :: m

        real(dp) :: shape

        if (ieee_is_nan(k) .or. .not. ieee_is_finite(k) .or. k < 0.0_dp .or. &
                .not. valid_probability(y)) then
            m = quiet_nan(y)
        else
            shape = floor(k) + 1.0_dp
            m = gammainccinv(shape, y)
        end if
    end function pdtri

    ! Chi-square CDF for v degrees of freedom.
    pure elemental function chdtr(v, x) result(value)
        real(dp), intent(in) :: v !! degrees of freedom, > 0
        real(dp), intent(in) :: x !! variate, >= 0; +infinity allowed
        real(dp) :: value

        if (.not. valid_positive_parameter(v) .or. ieee_is_nan(x) .or. x < 0.0_dp) then
            value = quiet_nan(x)
        else
            value = gammainc(0.5_dp * v, 0.5_dp * x)
        end if
    end function chdtr

    ! Chi-square survival function for v degrees of freedom.
    pure elemental function chdtrc(v, x) result(value)
        real(dp), intent(in) :: v !! degrees of freedom, > 0
        real(dp), intent(in) :: x !! variate, >= 0; +infinity allowed
        real(dp) :: value

        if (.not. valid_positive_parameter(v) .or. ieee_is_nan(x) .or. x < 0.0_dp) then
            value = quiet_nan(x)
        else
            value = gammaincc(0.5_dp * v, 0.5_dp * x)
        end if
    end function chdtrc

    ! Inverse chi-square survival function: chdtrc(v, x) = p.
    pure elemental function chdtri(v, p) result(x)
        real(dp), intent(in) :: v !! degrees of freedom, finite and > 0
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp) :: x

        real(dp) :: half_x

        if (.not. valid_positive_finite_parameter(v) .or. .not. valid_probability(p)) then
            x = quiet_nan(p)
        else
            half_x = gammainccinv(0.5_dp * v, p)
            x = scale_nonnegative(half_x, 2.0_dp)
        end if
    end function chdtri

    ! F-distribution CDF.
    pure elemental function fdtr(dfn, dfd, x) result(value)
        real(dp), intent(in) :: dfn !! numerator degrees of freedom, finite and > 0
        real(dp), intent(in) :: dfd !! denominator degrees of freedom, finite and > 0
        real(dp), intent(in) :: x !! variate, >= 0; +infinity allowed
        real(dp) :: value

        real(dp) :: q
        real(dp) :: z

        if (.not. valid_positive_finite_parameter(dfn) .or. &
                .not. valid_positive_finite_parameter(dfd) .or. &
                ieee_is_nan(x) .or. x < 0.0_dp) then
            value = quiet_nan(x)
        else if (x == 0.0_dp) then
            value = 0.0_dp
        else if (.not. ieee_is_finite(x)) then
            value = 1.0_dp
        else
            call f_beta_coordinates(dfn, dfd, x, z, q)
            if (z <= 0.5_dp) then
                value = betainc(0.5_dp * dfn, 0.5_dp * dfd, z)
            else
                value = betaincc(0.5_dp * dfd, 0.5_dp * dfn, q)
            end if
        end if
    end function fdtr

    ! F-distribution survival function.
    pure elemental function fdtrc(dfn, dfd, x) result(value)
        real(dp), intent(in) :: dfn !! numerator degrees of freedom, finite and > 0
        real(dp), intent(in) :: dfd !! denominator degrees of freedom, finite and > 0
        real(dp), intent(in) :: x !! variate, >= 0; +infinity allowed
        real(dp) :: value

        real(dp) :: q
        real(dp) :: z

        if (.not. valid_positive_finite_parameter(dfn) .or. &
                .not. valid_positive_finite_parameter(dfd) .or. &
                ieee_is_nan(x) .or. x < 0.0_dp) then
            value = quiet_nan(x)
        else if (x == 0.0_dp) then
            value = 1.0_dp
        else if (.not. ieee_is_finite(x)) then
            value = 0.0_dp
        else
            call f_beta_coordinates(dfn, dfd, x, z, q)
            if (z <= 0.5_dp) then
                value = betaincc(0.5_dp * dfn, 0.5_dp * dfd, z)
            else
                value = betainc(0.5_dp * dfd, 0.5_dp * dfn, q)
            end if
        end if
    end function fdtrc

    ! F-distribution quantile, inverse of fdtr with respect to x.
    pure elemental function fdtri(dfn, dfd, p) result(x)
        real(dp), intent(in) :: dfn !! numerator degrees of freedom, finite and > 0
        real(dp), intent(in) :: dfd !! denominator degrees of freedom, finite and > 0
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp) :: x

        real(dp) :: log_x
        real(dp) :: q
        real(dp) :: z

        if (.not. valid_positive_finite_parameter(dfn) .or. &
                .not. valid_positive_finite_parameter(dfd) .or. &
                .not. valid_probability(p)) then
            x = quiet_nan(p)
        else if (p <= 0.0_dp) then
            x = 0.0_dp
        else if (p >= 1.0_dp) then
            x = positive_infinity(p)
        else
            call beta_inverse_xy(0.5_dp * dfn, 0.5_dp * dfd, p, 1.0_dp - p, z, q)
            log_x = log(dfd) - log(dfn) + log(z) - log(q)
            if (log_x >= log(huge(1.0_dp))) then
                x = positive_infinity(p)
            else
                x = exp(log_x)
            end if
        end if
    end function fdtri

    ! Gamma-distribution CDF with rate a and shape b.
    pure elemental function gdtr(a, b, x) result(value)
        real(dp), intent(in) :: a !! rate parameter, finite and > 0
        real(dp), intent(in) :: b !! shape parameter, finite and > 0
        real(dp), intent(in) :: x !! variate, >= 0; +infinity allowed
        real(dp) :: value

        real(dp) :: ax

        if (.not. valid_positive_finite_parameter(a) .or. &
                .not. valid_positive_finite_parameter(b) .or. &
                ieee_is_nan(x) .or. x < 0.0_dp) then
            value = quiet_nan(x)
        else
            ax = safe_product(a, x)
            value = gammainc(b, ax)
        end if
    end function gdtr

    ! Gamma-distribution survival function with rate a and shape b.
    pure elemental function gdtrc(a, b, x) result(value)
        real(dp), intent(in) :: a !! rate parameter, finite and > 0
        real(dp), intent(in) :: b !! shape parameter, finite and > 0
        real(dp), intent(in) :: x !! variate, >= 0; +infinity allowed
        real(dp) :: value

        real(dp) :: ax

        if (.not. valid_positive_finite_parameter(a) .or. &
                .not. valid_positive_finite_parameter(b) .or. &
                ieee_is_nan(x) .or. x < 0.0_dp) then
            value = quiet_nan(x)
        else
            ax = safe_product(a, x)
            value = gammaincc(b, ax)
        end if
    end function gdtrc

    ! Inverse of gdtr with respect to the rate parameter a.
    pure elemental function gdtria(p, b, x) result(a)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: b !! shape parameter, finite and > 0
        real(dp), intent(in) :: x !! variate, finite and > 0
        real(dp) :: a

        real(dp) :: ax

        if (.not. valid_probability(p) .or. .not. valid_positive_finite_parameter(b) .or. &
                .not. valid_positive_finite_parameter(x)) then
            a = quiet_nan(p)
        else
            ax = gammaincinv(b, p)
            a = positive_ratio(ax, x)
        end if
    end function gdtria

    ! Gamma-distribution quantile, inverse of gdtr with respect to x.
    pure elemental function gdtrix(a, b, p) result(x)
        real(dp), intent(in) :: a !! rate parameter, finite and > 0
        real(dp), intent(in) :: b !! shape parameter, finite and > 0
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp) :: x

        real(dp) :: ax

        if (.not. valid_positive_finite_parameter(a) .or. &
                .not. valid_positive_finite_parameter(b) .or. &
                .not. valid_probability(p)) then
            x = quiet_nan(p)
        else
            ax = gammaincinv(b, p)
            x = positive_ratio(ax, a)
        end if
    end function gdtrix

    ! Student-t cumulative distribution function.
    pure elemental function stdtr(df, t) result(value)
        real(dp), intent(in) :: df !! degrees of freedom, > 0; +infinity gives normal
        real(dp), intent(in) :: t !! Student-t variate
        real(dp) :: value

        real(dp) :: beta_x
        real(dp) :: inverse_ratio
        real(dp) :: log_ratio
        real(dp) :: ratio
        real(dp) :: tail

        if (ieee_is_nan(df) .or. df <= 0.0_dp .or. ieee_is_nan(t)) then
            value = quiet_nan(t)
        else if (.not. ieee_is_finite(df)) then
            value = ndtr(t)
        else if (.not. ieee_is_finite(t)) then
            if (t > 0.0_dp) then
                value = 1.0_dp
            else
                value = 0.0_dp
            end if
        else if (t == 0.0_dp) then
            value = 0.5_dp
        else
            log_ratio = log(abs(t)) - 0.5_dp * log(df)
            if (log_ratio <= 0.0_dp) then
                ratio = exp(log_ratio)
                beta_x = 1.0_dp / (1.0_dp + ratio * ratio)
            else
                inverse_ratio = exp(-log_ratio)
                beta_x = inverse_ratio * inverse_ratio / &
                    (1.0_dp + inverse_ratio * inverse_ratio)
            end if
            tail = 0.5_dp * betainc(0.5_dp * df, 0.5_dp, beta_x)
            if (t < 0.0_dp) then
                value = tail
            else
                value = 1.0_dp - tail
            end if
        end if
    end function stdtr

    ! Student-t quantile, inverse of stdtr with respect to t.
    pure elemental function stdtrit(df, p) result(t)
        real(dp), intent(in) :: df !! degrees of freedom, > 0; +infinity gives normal
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp) :: t

        real(dp) :: beta_x
        real(dp) :: beta_y
        real(dp) :: log_magnitude
        real(dp) :: lower_tail
        real(dp) :: magnitude

        if (ieee_is_nan(df) .or. df <= 0.0_dp .or. .not. valid_probability(p)) then
            t = quiet_nan(p)
        else if (p <= 0.0_dp) then
            t = negative_infinity(p)
        else if (p >= 1.0_dp) then
            t = positive_infinity(p)
        else if (.not. ieee_is_finite(df)) then
            t = ndtri(p)
        else if (p == 0.5_dp) then
            t = 0.0_dp
        else
            lower_tail = 2.0_dp * min(p, 1.0_dp - p)
            call beta_inverse_xy(0.5_dp * df, 0.5_dp, lower_tail, &
                1.0_dp - lower_tail, beta_x, beta_y)
            if (beta_x <= 0.0_dp) then
                magnitude = positive_infinity(p)
            else
                log_magnitude = 0.5_dp * (log(df) + log(beta_y) - log(beta_x))
                if (log_magnitude >= log(huge(1.0_dp))) then
                    magnitude = positive_infinity(p)
                else
                    magnitude = exp(log_magnitude)
                end if
            end if
            if (p < 0.5_dp) then
                t = -magnitude
            else
                t = magnitude
            end if
        end if
    end function stdtrit


    ! Inverse of the continuous binomial-CDF extension with respect to k.
    pure elemental function bdtrik(y, n, p) result(k)
        real(dp), intent(in) :: y !! target binomial CDF in [0, 1]
        real(dp), intent(in) :: n !! number of trials, finite and >= 0
        real(dp), intent(in) :: p !! Bernoulli success probability in [0, 1]
        real(dp) :: k

        integer :: iteration
        real(dp) :: cdf0
        real(dp) :: hi
        real(dp) :: lo
        real(dp) :: mid
        real(dp) :: value

        if (.not. valid_probability(y) .or. .not. ieee_is_finite(n) .or. n < 0.0_dp .or. &
                .not. valid_probability(p)) then
            k = quiet_nan(y)
            return
        end if
        if (n == 0.0_dp .or. p == 0.0_dp) then
            k = 0.0_dp
            return
        end if
        if (p == 1.0_dp) then
            if (y <= 0.0_dp) then
                k = 0.0_dp
            else
                k = n
            end if
            return
        end if

        cdf0 = exp(n * log(1.0_dp - p))
        if (y <= cdf0) then
            k = 0.0_dp
            return
        else if (y >= 1.0_dp) then
            k = n
            return
        end if

        lo = 0.0_dp
        hi = n
        do iteration = 1, 100
            mid = lo + 0.5_dp * (hi - lo)
            value = betaincc(mid + 1.0_dp, n - mid, p)
            if (value < y) then
                lo = mid
            else
                hi = mid
            end if
            if (hi - lo <= 16.0_dp * epsilon(1.0_dp) * max(1.0_dp, hi)) exit
        end do
        k = lo + 0.5_dp * (hi - lo)
    end function bdtrik

    ! Inverse of the continuous binomial-CDF extension with respect to n.
    pure elemental function bdtrin(k, y, p) result(n)
        real(dp), intent(in) :: k !! success threshold, finite and >= 0
        real(dp), intent(in) :: y !! target binomial CDF in [0, 1]
        real(dp), intent(in) :: p !! Bernoulli success probability in [0, 1]
        real(dp) :: n

        real(dp) :: q

        if (ieee_is_nan(k) .or. .not. ieee_is_finite(k) .or. k < 0.0_dp .or. &
                .not. valid_probability(y) .or. .not. valid_probability(p)) then
            n = quiet_nan(y)
            return
        end if
        if (p <= 0.0_dp) then
            if (y == 1.0_dp) then
                n = k
            else
                n = quiet_nan(y)
            end if
            return
        else if (p >= 1.0_dp) then
            if (y == 1.0_dp) then
                n = k
            else if (y == 0.0_dp) then
                n = positive_infinity(y)
            else
                n = quiet_nan(y)
            end if
            return
        end if
        if (y >= 1.0_dp) then
            n = k
        else if (y <= 0.0_dp) then
            n = positive_infinity(y)
        else
            q = solve_positive_probability_parameter(y, root_bdtrin, k, p, .false.)
            n = k + q
        end if
    end function bdtrin

    ! Inverse of betainc(a, b, x) with respect to a.
    pure elemental function btdtria(p, b, x) result(a)
        real(dp), intent(in) :: p !! target regularized beta value in [0, 1]
        real(dp), intent(in) :: b !! second beta shape, finite and > 0
        real(dp), intent(in) :: x !! beta argument, strictly between 0 and 1
        real(dp) :: a

        if (.not. valid_probability(p) .or. .not. valid_positive_finite_parameter(b) .or. &
                .not. ieee_is_finite(x) .or. x <= 0.0_dp .or. x >= 1.0_dp) then
            a = quiet_nan(p)
        else if (p <= 0.0_dp) then
            a = positive_infinity(p)
        else if (p >= 1.0_dp) then
            a = 0.0_dp
        else
            a = solve_positive_probability_parameter(p, root_btdtria, b, x, .false.)
        end if
    end function btdtria

    ! Inverse of betainc(a, b, x) with respect to b.
    pure elemental function btdtrib(a, p, x) result(b)
        real(dp), intent(in) :: a !! first beta shape, finite and > 0
        real(dp), intent(in) :: p !! target regularized beta value in [0, 1]
        real(dp), intent(in) :: x !! beta argument, strictly between 0 and 1
        real(dp) :: b

        if (.not. valid_positive_finite_parameter(a) .or. .not. valid_probability(p) .or. &
                .not. ieee_is_finite(x) .or. x <= 0.0_dp .or. x >= 1.0_dp) then
            b = quiet_nan(p)
        else if (p <= 0.0_dp) then
            b = 0.0_dp
        else if (p >= 1.0_dp) then
            b = positive_infinity(p)
        else
            b = solve_positive_probability_parameter(p, root_btdtrib, a, x, .true.)
        end if
    end function btdtrib

    ! Inverse of chdtr with respect to the degrees of freedom v.
    pure elemental function chdtriv(p, x) result(v)
        real(dp), intent(in) :: p !! target chi-square CDF in [0, 1]
        real(dp), intent(in) :: x !! positive chi-square variate
        real(dp) :: v

        if (.not. valid_probability(p) .or. .not. valid_positive_finite_parameter(x)) then
            v = quiet_nan(p)
        else if (p <= 0.0_dp) then
            v = positive_infinity(p)
        else if (p >= 1.0_dp) then
            v = 0.0_dp
        else
            v = solve_positive_probability_parameter(p, root_chdtriv, x, 0.0_dp, .false.)
        end if
    end function chdtriv

    ! Inverse of fdtr with respect to denominator degrees of freedom. When the
    ! equation has two positive roots, return the smaller one.
    pure elemental function fdtridfd(dfn, p, x) result(dfd)
        real(dp), intent(in) :: dfn !! numerator degrees of freedom, finite and > 0
        real(dp), intent(in) :: p !! target F CDF in [0, 1]
        real(dp), intent(in) :: x !! positive finite F variate
        real(dp) :: dfd

        integer :: iteration
        real(dp), parameter :: growth = 1.4142135623730950488_dp
        real(dp) :: hi
        real(dp) :: lo
        real(dp) :: mid
        real(dp) :: value

        if (.not. valid_positive_finite_parameter(dfn) .or. .not. valid_probability(p) .or. &
                .not. valid_positive_finite_parameter(x)) then
            dfd = quiet_nan(p)
            return
        end if
        if (p <= 0.0_dp) then
            dfd = 0.0_dp
            return
        else if (p >= 1.0_dp) then
            dfd = quiet_nan(p)
            return
        end if

        lo = 1.0e-12_dp
        value = fdtr(dfn, lo, x)
        do while (value >= p .and. lo > tiny(1.0_dp) * growth)
            hi = lo
            lo = lo / growth
            value = fdtr(dfn, lo, x)
        end do
        hi = lo * growth
        do iteration = 1, 4096
            value = fdtr(dfn, hi, x)
            if (value >= p) exit
            lo = hi
            if (hi > huge(1.0_dp) / growth) then
                dfd = quiet_nan(p)
                return
            end if
            hi = hi * growth
        end do
        if (value < p) then
            dfd = quiet_nan(p)
            return
        end if

        do iteration = 1, 120
            mid = sqrt(lo) * sqrt(hi)
            value = fdtr(dfn, mid, x)
            if (value < p) then
                lo = mid
            else
                hi = mid
            end if
            if (hi <= lo * (1.0_dp + 16.0_dp * epsilon(1.0_dp))) exit
        end do
        dfd = sqrt(lo) * sqrt(hi)
    end function fdtridfd

    ! Inverse of gdtr with respect to gamma shape b.
    pure elemental function gdtrib(a, p, x) result(b)
        real(dp), intent(in) :: a !! gamma rate, finite and > 0
        real(dp), intent(in) :: p !! target gamma CDF in [0, 1]
        real(dp), intent(in) :: x !! positive finite variate
        real(dp) :: b

        if (.not. valid_positive_finite_parameter(a) .or. .not. valid_probability(p) .or. &
                .not. valid_positive_finite_parameter(x)) then
            b = quiet_nan(p)
        else if (p <= 0.0_dp) then
            b = positive_infinity(p)
        else if (p >= 1.0_dp) then
            b = 0.0_dp
        else
            b = solve_positive_probability_parameter(p, root_gdtrib, a, x, .false.)
        end if
    end function gdtrib

    ! Inverse of the continuous negative-binomial CDF with respect to k.
    pure elemental function nbdtrik(y, n, p) result(k)
        real(dp), intent(in) :: y !! target negative-binomial CDF in [0, 1]
        real(dp), intent(in) :: n !! success-shape parameter, finite and > 0
        real(dp), intent(in) :: p !! Bernoulli success probability in [0, 1]
        real(dp) :: k

        real(dp) :: cdf0
        real(dp) :: shape

        if (.not. valid_probability(y) .or. .not. valid_positive_finite_parameter(n) .or. &
                .not. valid_probability(p)) then
            k = quiet_nan(y)
            return
        end if
        if (p <= 0.0_dp) then
            k = positive_infinity(y)
            return
        else if (p >= 1.0_dp) then
            k = 0.0_dp
            return
        end if
        cdf0 = exp(n * log(p))
        if (y <= cdf0) then
            k = 0.0_dp
        else if (y >= 1.0_dp) then
            k = positive_infinity(y)
        else
            shape = solve_positive_probability_parameter(y, root_nbdtrik, n, p, .true.)
            k = max(0.0_dp, shape - 1.0_dp)
        end if
    end function nbdtrik

    ! Inverse of the continuous negative-binomial CDF with respect to n.
    pure elemental function nbdtrin(k, y, p) result(n)
        real(dp), intent(in) :: k !! failure threshold, finite and >= 0
        real(dp), intent(in) :: y !! target negative-binomial CDF in [0, 1]
        real(dp), intent(in) :: p !! Bernoulli success probability in [0, 1]
        real(dp) :: n

        if (ieee_is_nan(k) .or. .not. ieee_is_finite(k) .or. k < 0.0_dp .or. &
                .not. valid_probability(y) .or. .not. valid_probability(p)) then
            n = quiet_nan(y)
            return
        end if
        if (p <= 0.0_dp) then
            n = positive_infinity(y)
        else if (p >= 1.0_dp) then
            n = 0.0_dp
        else if (y <= 0.0_dp) then
            n = positive_infinity(y)
        else if (y >= 1.0_dp) then
            n = 0.0_dp
        else
            n = solve_positive_probability_parameter(y, root_nbdtrin, k, p, .false.)
        end if
    end function nbdtrin

    ! Inverse of the continuous Poisson-CDF extension with respect to k.
    pure elemental function pdtrik(p, m) result(k)
        real(dp), intent(in) :: p !! target Poisson CDF in [0, 1]
        real(dp), intent(in) :: m !! Poisson mean, finite and >= 0
        real(dp) :: k

        real(dp) :: cdf0
        real(dp) :: shape

        if (.not. valid_probability(p) .or. ieee_is_nan(m) .or. &
                .not. ieee_is_finite(m) .or. m < 0.0_dp) then
            k = quiet_nan(p)
            return
        end if
        if (m == 0.0_dp) then
            if (p >= 1.0_dp) then
                k = 0.0_dp
            else
                k = quiet_nan(p)
            end if
            return
        end if
        cdf0 = exp(-m)
        if (p <= cdf0) then
            k = 0.0_dp
        else if (p >= 1.0_dp) then
            k = positive_infinity(p)
        else
            shape = solve_positive_probability_parameter(p, root_pdtrik, m, 0.0_dp, .true.)
            k = max(0.0_dp, shape - 1.0_dp)
        end if
    end function pdtrik

    ! Inverse of stdtr with respect to degrees of freedom.
    pure elemental function stdtridf(p, t) result(df)
        real(dp), intent(in) :: p !! target Student-t CDF in [0, 1]
        real(dp), intent(in) :: t !! finite Student-t variate
        real(dp) :: df

        real(dp) :: normal_limit
        real(dp) :: p_work
        real(dp) :: t_work

        if (.not. valid_probability(p) .or. ieee_is_nan(t) .or. .not. ieee_is_finite(t)) then
            df = quiet_nan(p)
            return
        end if
        if (t == 0.0_dp) then
            df = quiet_nan(p)
            return
        end if
        t_work = abs(t)
        if (t > 0.0_dp) then
            p_work = p
        else
            p_work = 1.0_dp - p
        end if
        normal_limit = ndtr(t_work)
        if (p_work < 0.5_dp .or. p_work > normal_limit) then
            df = quiet_nan(p)
        else if (p_work == 0.5_dp) then
            df = 0.0_dp
        else if (p_work == normal_limit) then
            df = positive_infinity(p)
        else
            df = solve_positive_probability_parameter(p_work, root_stdtridf, t_work, 0.0_dp, .true.)
        end if
    end function stdtridf

    ! Mean of a normal distribution from a CDF value, standard deviation, and x.
    pure elemental function nrdtrimn(p, std, x) result(mn)
        real(dp), intent(in) :: p !! normal CDF value in [0, 1]
        real(dp), intent(in) :: std !! standard deviation, finite and > 0
        real(dp), intent(in) :: x !! finite quantile
        real(dp) :: mn

        real(dp) :: z

        if (.not. valid_probability(p) .or. .not. valid_positive_finite_parameter(std) .or. &
                ieee_is_nan(x) .or. .not. ieee_is_finite(x)) then
            mn = quiet_nan(p)
            return
        end if
        z = ndtri(p)
        if (.not. ieee_is_finite(z)) then
            if (z < 0.0_dp) then
                mn = positive_infinity(p)
            else
                mn = negative_infinity(p)
            end if
        else if (abs(z) > 1.0_dp) then
            if (std > huge(1.0_dp) / abs(z)) then
                if (z > 0.0_dp) then
                    mn = negative_infinity(p)
                else
                    mn = positive_infinity(p)
                end if
            else
                mn = x - std * z
            end if
        else
            mn = x - std * z
        end if
    end function nrdtrimn

    ! Standard deviation of a normal distribution from mean, CDF value, and x.
    pure elemental function nrdtrisd(mn, p, x) result(std)
        real(dp), intent(in) :: mn !! finite normal mean
        real(dp), intent(in) :: p !! normal CDF value strictly between 0 and 1
        real(dp), intent(in) :: x !! finite quantile
        real(dp) :: std

        real(dp) :: delta
        real(dp) :: z

        if (ieee_is_nan(mn) .or. .not. ieee_is_finite(mn) .or. &
                .not. valid_probability(p) .or. p <= 0.0_dp .or. p >= 1.0_dp .or. &
                ieee_is_nan(x) .or. .not. ieee_is_finite(x)) then
            std = quiet_nan(p)
            return
        end if
        z = ndtri(p)
        delta = x - mn
        if (z == 0.0_dp .or. delta == 0.0_dp .or. sign(1.0_dp, delta) /= sign(1.0_dp, z)) then
            std = quiet_nan(p)
        else if (abs(z) < 1.0_dp) then
            if (abs(delta) > huge(1.0_dp) * abs(z)) then
                std = positive_infinity(p)
            else
                std = delta / z
            end if
        else
            std = delta / z
        end if
    end function nrdtrisd


    ! Solve a monotone probability equation over a positive parameter. The
    ! bracket is expanded multiplicatively and then bisected geometrically.
    pure function solve_positive_probability_parameter(target, family, arg1, arg2, increasing) result(root)
        real(dp), intent(in) :: target !! target probability strictly between 0 and 1
        integer, intent(in) :: family !! selector for the probability equation
        real(dp), intent(in) :: arg1 !! first fixed argument for the selected equation
        real(dp), intent(in) :: arg2 !! second fixed argument for the selected equation
        logical, intent(in) :: increasing !! true when the probability rises with the parameter
        real(dp) :: root

        integer :: iteration
        real(dp) :: hi
        real(dp) :: lo
        real(dp) :: mid
        real(dp) :: value

        root = 1.0_dp
        value = parameter_probability(root, family, arg1, arg2)
        if (parameter_too_low(value, target, increasing)) then
            lo = root
            hi = 2.0_dp
            do iteration = 1, 2048
                value = parameter_probability(hi, family, arg1, arg2)
                if (.not. parameter_too_low(value, target, increasing)) exit
                lo = hi
                if (hi > huge(1.0_dp) / 2.0_dp) then
                    root = positive_infinity(target)
                    return
                end if
                hi = 2.0_dp * hi
            end do
        else
            hi = root
            lo = 0.5_dp
            do iteration = 1, 2048
                value = parameter_probability(lo, family, arg1, arg2)
                if (parameter_too_low(value, target, increasing)) exit
                hi = lo
                if (lo <= 2.0_dp * tiny(1.0_dp)) then
                    root = 0.0_dp
                    return
                end if
                lo = 0.5_dp * lo
            end do
        end if

        do iteration = 1, 120
            mid = sqrt(lo) * sqrt(hi)
            value = parameter_probability(mid, family, arg1, arg2)
            if (parameter_too_low(value, target, increasing)) then
                lo = mid
            else
                hi = mid
            end if
            if (hi <= lo * (1.0_dp + 16.0_dp * epsilon(1.0_dp))) exit
        end do
        root = sqrt(lo) * sqrt(hi)
    end function solve_positive_probability_parameter

    pure function parameter_probability(parameter, family, arg1, arg2) result(value)
        real(dp), intent(in) :: parameter !! positive parameter being inverted
        integer, intent(in) :: family !! selector for the probability equation
        real(dp), intent(in) :: arg1 !! first fixed argument
        real(dp), intent(in) :: arg2 !! second fixed argument
        real(dp) :: value

        select case (family)
        case (root_bdtrin)
            value = betaincc(arg1 + 1.0_dp, parameter, arg2)
        case (root_btdtria)
            value = betainc(parameter, arg1, arg2)
        case (root_btdtrib)
            value = betainc(arg1, parameter, arg2)
        case (root_chdtriv)
            value = gammainc(0.5_dp * parameter, 0.5_dp * arg1)
        case (root_gdtrib)
            value = gammainc(parameter, safe_product(arg1, arg2))
        case (root_nbdtrik)
            value = betainc(arg1, parameter, arg2)
        case (root_nbdtrin)
            value = betainc(parameter, arg1 + 1.0_dp, arg2)
        case (root_pdtrik)
            value = gammaincc(parameter, arg1)
        case (root_stdtridf)
            value = stdtr(parameter, arg1)
        case default
            value = quiet_nan(parameter)
        end select
    end function parameter_probability

    pure logical function parameter_too_low(value, target, increasing) result(too_low)
        real(dp), intent(in) :: value !! probability at the current positive parameter
        real(dp), intent(in) :: target !! target probability
        logical, intent(in) :: increasing !! true when probability increases with parameter

        if (increasing) then
            too_low = value < target
        else
            too_low = value > target
        end if
    end function parameter_too_low

    pure elemental logical function valid_probability(p) result(valid)
        real(dp), intent(in) :: p !! probability candidate

        valid = .not. ieee_is_nan(p) .and. ieee_is_finite(p) .and. &
            p >= 0.0_dp .and. p <= 1.0_dp
    end function valid_probability

    pure elemental logical function valid_positive_parameter(x) result(valid)
        real(dp), intent(in) :: x !! positive parameter candidate; +infinity allowed

        valid = .not. ieee_is_nan(x) .and. x > 0.0_dp
    end function valid_positive_parameter

    pure elemental logical function valid_positive_finite_parameter(x) result(valid)
        real(dp), intent(in) :: x !! positive finite parameter candidate

        valid = ieee_is_finite(x) .and. x > 0.0_dp
    end function valid_positive_finite_parameter

    pure elemental subroutine f_beta_coordinates(dfn, dfd, x, z, q)
        real(dp), intent(in) :: dfn !! numerator degrees of freedom, finite and > 0
        real(dp), intent(in) :: dfd !! denominator degrees of freedom, finite and > 0
        real(dp), intent(in) :: x !! positive finite F variate
        real(dp), intent(out) :: z !! dfn*x/(dfd+dfn*x)
        real(dp), intent(out) :: q !! complementary coordinate 1-z

        real(dp) :: inverse_ratio
        real(dp) :: log_ratio
        real(dp) :: ratio

        log_ratio = log(dfn) + log(x) - log(dfd)
        if (log_ratio <= 0.0_dp) then
            ratio = exp(log_ratio)
            z = ratio / (1.0_dp + ratio)
            q = 1.0_dp / (1.0_dp + ratio)
        else
            inverse_ratio = exp(-log_ratio)
            z = 1.0_dp / (1.0_dp + inverse_ratio)
            q = inverse_ratio / (1.0_dp + inverse_ratio)
        end if
    end subroutine f_beta_coordinates

    pure elemental function positive_ratio(numerator, denominator) result(value)
        real(dp), intent(in) :: numerator !! nonnegative numerator; +infinity allowed
        real(dp), intent(in) :: denominator !! positive finite denominator
        real(dp) :: value

        if (numerator == 0.0_dp) then
            value = 0.0_dp
        else if (.not. ieee_is_finite(numerator)) then
            value = positive_infinity(denominator)
        else if (denominator >= 1.0_dp) then
            value = numerator / denominator
        else if (numerator > huge(1.0_dp) * denominator) then
            value = positive_infinity(denominator)
        else
            value = numerator / denominator
        end if
    end function positive_ratio

    pure elemental function safe_product(a, x) result(value)
        real(dp), intent(in) :: a !! positive finite multiplier
        real(dp), intent(in) :: x !! nonnegative factor; +infinity allowed
        real(dp) :: value

        if (x == 0.0_dp) then
            value = 0.0_dp
        else if (.not. ieee_is_finite(x)) then
            value = positive_infinity(a)
        else if (a <= 1.0_dp) then
            value = a * x
        else if (x > huge(1.0_dp) / a) then
            value = positive_infinity(a)
        else
            value = a * x
        end if
    end function safe_product

    pure elemental function scale_nonnegative(x, factor) result(value)
        real(dp), intent(in) :: x !! nonnegative value; +infinity allowed
        real(dp), intent(in) :: factor !! positive finite scale factor
        real(dp) :: value

        if (x == 0.0_dp) then
            value = 0.0_dp
        else if (.not. ieee_is_finite(x)) then
            value = positive_infinity(factor)
        else if (x > huge(1.0_dp) / factor) then
            value = positive_infinity(factor)
        else
            value = factor * x
        end if
    end function scale_nonnegative

end module scifort_special_distribution_helpers
