! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Beta-negative-binomial distribution matching scipy.stats.betanbinom.
module scifort_betanbinom
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_log_gamma, only : log_beta
    use scifort_math, only : expm1_safe, negative_infinity, positive_infinity, quiet_nan
    implicit none
    private

    public :: betanbinom_cdf, betanbinom_isf, betanbinom_logcdf, betanbinom_logpmf
    public :: betanbinom_logsf, betanbinom_pmf, betanbinom_ppf, betanbinom_sf

    interface betanbinom_pmf
        module procedure betanbinom_pmf_real
        module procedure betanbinom_pmf_int
    end interface
    interface betanbinom_logpmf
        module procedure betanbinom_logpmf_real
        module procedure betanbinom_logpmf_int
    end interface
    interface betanbinom_cdf
        module procedure betanbinom_cdf_real
        module procedure betanbinom_cdf_int
    end interface
    interface betanbinom_sf
        module procedure betanbinom_sf_real
        module procedure betanbinom_sf_int
    end interface
    interface betanbinom_logcdf
        module procedure betanbinom_logcdf_real
        module procedure betanbinom_logcdf_int
    end interface
    interface betanbinom_logsf
        module procedure betanbinom_logsf_real
        module procedure betanbinom_logsf_int
    end interface
contains
    pure elemental function betanbinom_logpmf_real(k, n, a, b, loc) result(y)
        real(dp), intent(in) :: k !! lattice observation
        real(dp), intent(in) :: n !! positive integer success count
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, count
        count = k - optional_loc(loc)
        if (.not. valid_parameters(n, a, b) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (.not. ieee_is_finite(count) .or. count < 0.0_dp .or. count /= aint(count)) then
            y = negative_infinity(k)
        else
            y = log_gamma(n + count) - log_gamma(n) - log_gamma(count + 1.0_dp) + &
                log_beta(a + n, b + count) - log_beta(a, b)
        end if
    end function
    pure elemental function betanbinom_logpmf_int(k, n, a, b, loc) result(y)
        integer, intent(in) :: k !! integer observation
        integer, intent(in) :: n !! positive integer success count
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = betanbinom_logpmf_real(real(k, dp), real(n, dp), a, b, loc)
    end function
    pure elemental function betanbinom_pmf_real(k, n, a, b, loc) result(y)
        real(dp), intent(in) :: k !! lattice observation
        real(dp), intent(in) :: n !! positive integer success count
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = exp(betanbinom_logpmf_real(k, n, a, b, loc))
    end function
    pure elemental function betanbinom_pmf_int(k, n, a, b, loc) result(y)
        integer, intent(in) :: k !! integer observation
        integer, intent(in) :: n !! positive integer success count
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = betanbinom_pmf_real(real(k, dp), real(n, dp), a, b, loc)
    end function
    pure elemental function betanbinom_logcdf_real(k, n, a, b, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of log(P(X <= k)); shifted floor is used
        real(dp), intent(in) :: n !! positive integer success count
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, count
        count = floor(k - optional_loc(loc))
        if (.not. valid_parameters(n, a, b) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (count < 0.0_dp) then
            y = negative_infinity(k)
        else if (.not. ieee_is_finite(count)) then
            y = 0.0_dp
        else
            y = log_forward_sum(int(count), n, a, b)
        end if
    end function
    pure elemental function betanbinom_logcdf_int(k, n, a, b, loc) result(y)
        integer, intent(in) :: k !! upper integer limit
        integer, intent(in) :: n !! positive integer success count
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = betanbinom_logcdf_real(real(k, dp), real(n, dp), a, b, loc)
    end function
    pure elemental function betanbinom_cdf_real(k, n, a, b, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of P(X <= k)
        real(dp), intent(in) :: n !! positive integer success count
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, ly
        ly = betanbinom_logcdf_real(k, n, a, b, loc)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = min(1.0_dp, exp(ly))
        end if
    end function
    pure elemental function betanbinom_cdf_int(k, n, a, b, loc) result(y)
        integer, intent(in) :: k !! upper integer limit
        integer, intent(in) :: n !! positive integer success count
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = betanbinom_cdf_real(real(k, dp), real(n, dp), a, b, loc)
    end function
    pure elemental function betanbinom_logsf_real(k, n, a, b, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of log(P(X > k)); shifted floor is used
        real(dp), intent(in) :: n !! positive integer success count
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, count, lcdf
        count = floor(k - optional_loc(loc))
        if (.not. valid_parameters(n, a, b) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (count < 0.0_dp) then
            y = 0.0_dp
        else if (.not. ieee_is_finite(count)) then
            y = negative_infinity(k)
        else
            lcdf = log_forward_sum(int(count), n, a, b)
            if (lcdf >= 0.0_dp) then
                y = negative_infinity(k)
            else
                y = log(-expm1_safe(lcdf))
            end if
        end if
    end function
    pure elemental function betanbinom_logsf_int(k, n, a, b, loc) result(y)
        integer, intent(in) :: k !! lower integer limit
        integer, intent(in) :: n !! positive integer success count
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = betanbinom_logsf_real(real(k, dp), real(n, dp), a, b, loc)
    end function
    pure elemental function betanbinom_sf_real(k, n, a, b, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of P(X > k)
        real(dp), intent(in) :: n !! positive integer success count
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, ly
        ly = betanbinom_logsf_real(k, n, a, b, loc)
        if (ieee_is_nan(ly)) then
            y = ly
        else
            y = exp(ly)
        end if
    end function
    pure elemental function betanbinom_sf_int(k, n, a, b, loc) result(y)
        integer, intent(in) :: k !! lower integer limit
        integer, intent(in) :: n !! positive integer success count
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = betanbinom_sf_real(real(k, dp), real(n, dp), a, b, loc)
    end function
    pure elemental function betanbinom_ppf(probability, n, a, b, loc) result(y)
        real(dp), intent(in) :: probability !! lower-tail probability in [0,1]
        real(dp), intent(in) :: n !! positive integer success count
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, shift
        integer :: lo, hi, mid
        shift = optional_loc(loc)
        if (.not. valid_parameters(n, a, b) .or. .not. ieee_is_finite(shift) .or. .not. valid_probability(probability)) then
            y = quiet_nan(probability)
        else if (probability <= 0.0_dp) then
            y = shift - 1.0_dp
        else if (probability >= 1.0_dp) then
            y = positive_infinity(probability)
        else
            lo = 0; hi = 1
            do while (betanbinom_cdf_real(shift + real(hi, dp), n, a, b, shift) < probability .and. hi < huge(hi)/2)
                hi = hi * 2
            end do
            do while (lo < hi)
                mid = lo + (hi-lo)/2
                if (betanbinom_cdf_real(shift + real(mid, dp), n, a, b, shift) >= probability) then
                    hi = mid
                else
                    lo = mid + 1
                end if
            end do
            y = shift + real(lo, dp)
        end if
    end function
    pure elemental function betanbinom_isf(probability, n, a, b, loc) result(y)
        real(dp), intent(in) :: probability !! upper-tail probability in [0,1]
        real(dp), intent(in) :: n !! positive integer success count
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, shift
        integer :: lo, hi, mid
        shift = optional_loc(loc)
        if (.not. valid_parameters(n, a, b) .or. .not. ieee_is_finite(shift) .or. .not. valid_probability(probability)) then
            y = quiet_nan(probability)
        else if (probability >= 1.0_dp) then
            y = shift - 1.0_dp
        else if (probability <= 0.0_dp) then
            y = positive_infinity(probability)
        else
            lo = 0; hi = 1
            do while (betanbinom_sf_real(shift + real(hi, dp), n, a, b, shift) > probability .and. hi < huge(hi)/2)
                hi = hi * 2
            end do
            do while (lo < hi)
                mid = lo + (hi-lo)/2
                if (betanbinom_sf_real(shift + real(mid, dp), n, a, b, shift) <= probability) then
                    hi = mid
                else
                    lo = mid + 1
                end if
            end do
            y = shift + real(lo, dp)
        end if
    end function
    pure function log_forward_sum(last, n, a, b) result(logsum)
        integer, intent(in) :: last !! last included count, nonnegative
        real(dp), intent(in) :: n !! positive integer success count
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp) :: logsum, logterm, logratio
        integer :: j
        logterm = log_beta(a+n, b) - log_beta(a, b)
        logsum = logterm
        do j = 0, last-1
            logratio = log(n + real(j,dp)) - log(real(j+1,dp)) + &
                log(b + real(j,dp)) - log(a+b+n+real(j,dp))
            logterm = logterm + logratio
            logsum = logadd(logsum, logterm)
        end do
        if (logsum > 0.0_dp .and. logsum < 64.0_dp*epsilon(1.0_dp)) logsum = 0.0_dp
    end function
    pure elemental function logadd(x, y) result(z)
        real(dp), intent(in) :: x !! first log-value
        real(dp), intent(in) :: y !! second log-value
        real(dp) :: z, m
        m = max(x,y)
        if (m == negative_infinity(m)) then
            z = m
        else
            z = m + log(exp(x-m)+exp(y-m))
        end if
    end function
    pure elemental function optional_loc(loc) result(value)
        real(dp), intent(in), optional :: loc !! support shift
        real(dp) :: value
        value = 0.0_dp
        if (present(loc)) value = loc
    end function
    pure elemental logical function valid_parameters(n, a, b)
        real(dp), intent(in) :: n !! success count to validate
        real(dp), intent(in) :: a !! first beta shape to validate
        real(dp), intent(in) :: b !! second beta shape to validate
        valid_parameters = ieee_is_finite(n) .and. ieee_is_finite(a) .and. ieee_is_finite(b) .and. &
            n >= 1.0_dp .and. n == aint(n) .and. a > 0.0_dp .and. b > 0.0_dp
    end function
    pure elemental logical function valid_probability(p)
        real(dp), intent(in) :: p !! probability to validate
        valid_probability = p >= 0.0_dp .and. p <= 1.0_dp
    end function
end module scifort_betanbinom
