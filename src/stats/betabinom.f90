! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Beta-binomial distribution matching scipy.stats.betabinom.
! The finite tails are accumulated in log space, preserving small probabilities
! without forming one tail by subtracting the other from one.

module scifort_betabinom
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_log_gamma, only : log_beta
    use scifort_math, only : log1p_safe, negative_infinity, quiet_nan
    implicit none
    private

    public :: betabinom_cdf, betabinom_isf, betabinom_logcdf, betabinom_logpmf
    public :: betabinom_logsf, betabinom_pmf, betabinom_ppf, betabinom_sf

    interface betabinom_pmf
        module procedure betabinom_pmf_real
        module procedure betabinom_pmf_int
    end interface betabinom_pmf
    interface betabinom_logpmf
        module procedure betabinom_logpmf_real
        module procedure betabinom_logpmf_int
    end interface betabinom_logpmf
    interface betabinom_cdf
        module procedure betabinom_cdf_real
        module procedure betabinom_cdf_int
    end interface betabinom_cdf
    interface betabinom_sf
        module procedure betabinom_sf_real
        module procedure betabinom_sf_int
    end interface betabinom_sf
    interface betabinom_logcdf
        module procedure betabinom_logcdf_real
        module procedure betabinom_logcdf_int
    end interface betabinom_logcdf
    interface betabinom_logsf
        module procedure betabinom_logsf_real
        module procedure betabinom_logsf_int
    end interface betabinom_logsf

contains

    pure elemental function betabinom_logpmf_real(k, n, a, b, loc) result(y)
        real(dp), intent(in) :: k !! lattice observation
        real(dp), intent(in) :: n !! nonnegative integer number of trials
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: count

        count = shifted_count(k, loc)
        if (.not. valid_parameters(n, a, b) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (.not. ieee_is_finite(count) .or. count /= aint(count) .or. &
                 count < 0.0_dp .or. count > n) then
            y = negative_infinity(k)
        else
            y = log_choose(n, count) + log_beta(count + a, n - count + b) - log_beta(a, b)
        end if
    end function betabinom_logpmf_real

    pure elemental function betabinom_logpmf_int(k, n, a, b, loc) result(y)
        integer, intent(in) :: k !! integer observation
        integer, intent(in) :: n !! nonnegative number of trials
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = betabinom_logpmf_real(real(k, dp), real(n, dp), a, b, loc)
    end function betabinom_logpmf_int

    pure elemental function betabinom_pmf_real(k, n, a, b, loc) result(y)
        real(dp), intent(in) :: k !! lattice observation
        real(dp), intent(in) :: n !! nonnegative integer number of trials
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = exp(betabinom_logpmf_real(k, n, a, b, loc))
    end function betabinom_pmf_real

    pure elemental function betabinom_pmf_int(k, n, a, b, loc) result(y)
        integer, intent(in) :: k !! integer observation
        integer, intent(in) :: n !! nonnegative number of trials
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = betabinom_pmf_real(real(k, dp), real(n, dp), a, b, loc)
    end function betabinom_pmf_int

    pure elemental function betabinom_cdf_real(k, n, a, b, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of P(X <= k); shifted floor is used
        real(dp), intent(in) :: n !! nonnegative integer number of trials
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: logcdf, logsf, sf
        call tails(floor_count(k, loc), n, a, b, y, sf, logcdf, logsf)
    end function betabinom_cdf_real

    pure elemental function betabinom_cdf_int(k, n, a, b, loc) result(y)
        integer, intent(in) :: k !! upper integer limit
        integer, intent(in) :: n !! nonnegative number of trials
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = betabinom_cdf_real(real(k, dp), real(n, dp), a, b, loc)
    end function betabinom_cdf_int

    pure elemental function betabinom_sf_real(k, n, a, b, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of P(X > k); shifted floor is used
        real(dp), intent(in) :: n !! nonnegative integer number of trials
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: cdf, logcdf, logsf
        call tails(floor_count(k, loc), n, a, b, cdf, y, logcdf, logsf)
    end function betabinom_sf_real

    pure elemental function betabinom_sf_int(k, n, a, b, loc) result(y)
        integer, intent(in) :: k !! lower integer limit
        integer, intent(in) :: n !! nonnegative number of trials
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = betabinom_sf_real(real(k, dp), real(n, dp), a, b, loc)
    end function betabinom_sf_int

    pure elemental function betabinom_logcdf_real(k, n, a, b, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of log(P(X <= k))
        real(dp), intent(in) :: n !! nonnegative integer number of trials
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: cdf, sf, logsf
        call tails(floor_count(k, loc), n, a, b, cdf, sf, y, logsf)
    end function betabinom_logcdf_real

    pure elemental function betabinom_logcdf_int(k, n, a, b, loc) result(y)
        integer, intent(in) :: k !! upper integer limit
        integer, intent(in) :: n !! nonnegative number of trials
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = betabinom_logcdf_real(real(k, dp), real(n, dp), a, b, loc)
    end function betabinom_logcdf_int

    pure elemental function betabinom_logsf_real(k, n, a, b, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of log(P(X > k))
        real(dp), intent(in) :: n !! nonnegative integer number of trials
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: cdf, sf, logcdf
        call tails(floor_count(k, loc), n, a, b, cdf, sf, logcdf, y)
    end function betabinom_logsf_real

    pure elemental function betabinom_logsf_int(k, n, a, b, loc) result(y)
        integer, intent(in) :: k !! lower integer limit
        integer, intent(in) :: n !! nonnegative number of trials
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = betabinom_logsf_real(real(k, dp), real(n, dp), a, b, loc)
    end function betabinom_logsf_int

    pure elemental function betabinom_ppf(probability, n, a, b, loc) result(y)
        real(dp), intent(in) :: probability !! lower-tail probability in [0,1]
        real(dp), intent(in) :: n !! nonnegative integer number of trials
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: shift
        integer :: lo, hi, mid

        shift = optional_loc(loc)
        if (.not. valid_parameters(n, a, b) .or. .not. ieee_is_finite(shift) .or. &
            .not. valid_probability(probability)) then
            y = quiet_nan(probability)
        else if (probability <= 0.0_dp) then
            y = shift - 1.0_dp
        else if (probability >= 1.0_dp) then
            y = shift + n
        else
            lo = 0
            hi = int(n)
            do while (lo < hi)
                mid = lo + (hi - lo) / 2
                if (betabinom_cdf_real(shift + real(mid, dp), n, a, b, shift) >= probability) then
                    hi = mid
                else
                    lo = mid + 1
                end if
            end do
            y = shift + real(lo, dp)
        end if
    end function betabinom_ppf

    pure elemental function betabinom_isf(probability, n, a, b, loc) result(y)
        real(dp), intent(in) :: probability !! upper-tail probability in [0,1]
        real(dp), intent(in) :: n !! nonnegative integer number of trials
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: shift
        integer :: lo, hi, mid

        shift = optional_loc(loc)
        if (.not. valid_parameters(n, a, b) .or. .not. ieee_is_finite(shift) .or. &
            .not. valid_probability(probability)) then
            y = quiet_nan(probability)
        else if (probability >= 1.0_dp) then
            y = shift - 1.0_dp
        else if (probability <= 0.0_dp) then
            y = shift + n
        else
            lo = 0
            hi = int(n)
            do while (lo < hi)
                mid = lo + (hi - lo) / 2
                if (betabinom_sf_real(shift + real(mid, dp), n, a, b, shift) <= probability) then
                    hi = mid
                else
                    lo = mid + 1
                end if
            end do
            y = shift + real(lo, dp)
        end if
    end function betabinom_isf

    pure subroutine tails(k, n, a, b, cdf, sf, logcdf, logsf)
        real(dp), intent(in) :: k !! shifted integer tail split, represented as real
        real(dp), intent(in) :: n !! nonnegative integer number of trials
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(out) :: cdf !! lower-tail probability
        real(dp), intent(out) :: sf !! upper-tail probability
        real(dp), intent(out) :: logcdf !! logarithm of lower tail
        real(dp), intent(out) :: logsf !! logarithm of upper tail
        integer :: ik, in

        if (.not. valid_parameters(n, a, b) .or. ieee_is_nan(k)) then
            cdf = quiet_nan(k); sf = cdf; logcdf = cdf; logsf = cdf
            return
        end if
        if (k < 0.0_dp) then
            cdf = 0.0_dp; sf = 1.0_dp; logcdf = negative_infinity(k); logsf = 0.0_dp
            return
        end if
        if (.not. ieee_is_finite(k) .or. k >= n) then
            cdf = 1.0_dp; sf = 0.0_dp; logcdf = 0.0_dp; logsf = negative_infinity(k)
            return
        end if
        ik = int(k)
        in = int(n)
        logcdf = logsum_range(0, ik, n, a, b)
        logsf = logsum_range(ik + 1, in, n, a, b)
        cdf = exp(logcdf)
        sf = exp(logsf)
    end subroutine tails

    pure function logsum_range(lo, hi, n, a, b) result(value)
        integer, intent(in) :: lo !! inclusive lower support index
        integer, intent(in) :: hi !! inclusive upper support index
        real(dp), intent(in) :: n !! number of trials
        real(dp), intent(in) :: a !! first beta shape
        real(dp), intent(in) :: b !! second beta shape
        real(dp) :: value
        real(dp) :: term
        integer :: j

        value = negative_infinity(1.0_dp)
        do j = lo, hi
            term = log_choose(n, real(j, dp)) + log_beta(real(j, dp) + a, n - real(j, dp) + b) - &
                log_beta(a, b)
            value = logadd(value, term)
        end do
    end function logsum_range

    pure elemental function log_choose(n, k) result(y)
        real(dp), intent(in) :: n !! nonnegative integer upper argument
        real(dp), intent(in) :: k !! integer lower argument in [0,n]
        real(dp) :: y
        y = log_gamma(n + 1.0_dp) - log_gamma(k + 1.0_dp) - log_gamma(n - k + 1.0_dp)
    end function log_choose

    pure elemental function logadd(a, b) result(y)
        real(dp), intent(in) :: a !! first logarithm
        real(dp), intent(in) :: b !! second logarithm
        real(dp) :: y
        real(dp) :: hi, lo
        hi = max(a, b)
        lo = min(a, b)
        if (.not. ieee_is_finite(hi)) then
            y = hi
        else
            y = hi + log1p_safe(exp(lo - hi))
        end if
    end function logadd

    pure elemental function shifted_count(k, loc) result(count)
        real(dp), intent(in) :: k !! observation
        real(dp), intent(in), optional :: loc !! support shift
        real(dp) :: count
        count = k - optional_loc(loc)
    end function shifted_count

    pure elemental function floor_count(k, loc) result(count)
        real(dp), intent(in) :: k !! CDF/SF argument
        real(dp), intent(in), optional :: loc !! support shift
        real(dp) :: count
        count = k - optional_loc(loc)
        if (ieee_is_finite(count)) count = floor(count)
    end function floor_count

    pure elemental function optional_loc(loc) result(shift)
        real(dp), intent(in), optional :: loc !! support shift
        real(dp) :: shift
        shift = 0.0_dp
        if (present(loc)) shift = loc
    end function optional_loc

    pure elemental logical function valid_parameters(n, a, b)
        real(dp), intent(in) :: n !! number of trials
        real(dp), intent(in) :: a !! first beta shape
        real(dp), intent(in) :: b !! second beta shape
        valid_parameters = ieee_is_finite(n) .and. n >= 0.0_dp .and. n == aint(n) .and. &
            n <= real(huge(0), dp) .and. ieee_is_finite(a) .and. a > 0.0_dp .and. &
            ieee_is_finite(b) .and. b > 0.0_dp
    end function valid_parameters

    pure elemental logical function valid_probability(p)
        real(dp), intent(in) :: p !! probability candidate
        valid_probability = ieee_is_finite(p) .and. p >= 0.0_dp .and. p <= 1.0_dp
    end function valid_probability

end module scifort_betabinom
