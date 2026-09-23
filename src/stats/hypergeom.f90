! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Hypergeometric distribution matching scipy.stats.hypergeom.
! Finite lower and upper tails are summed independently in log space.

module scifort_hypergeom
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : floor_real, log1p_safe, negative_infinity, quiet_nan
    implicit none
    private

    public :: hypergeom_cdf, hypergeom_isf, hypergeom_logcdf, hypergeom_logpmf
    public :: hypergeom_logsf, hypergeom_pmf, hypergeom_ppf, hypergeom_sf

    interface hypergeom_pmf
        module procedure hypergeom_pmf_real
        module procedure hypergeom_pmf_int
    end interface hypergeom_pmf
    interface hypergeom_logpmf
        module procedure hypergeom_logpmf_real
        module procedure hypergeom_logpmf_int
    end interface hypergeom_logpmf
    interface hypergeom_cdf
        module procedure hypergeom_cdf_real
        module procedure hypergeom_cdf_int
    end interface hypergeom_cdf
    interface hypergeom_sf
        module procedure hypergeom_sf_real
        module procedure hypergeom_sf_int
    end interface hypergeom_sf
    interface hypergeom_logcdf
        module procedure hypergeom_logcdf_real
        module procedure hypergeom_logcdf_int
    end interface hypergeom_logcdf
    interface hypergeom_logsf
        module procedure hypergeom_logsf_real
        module procedure hypergeom_logsf_int
    end interface hypergeom_logsf

contains

    pure elemental function hypergeom_logpmf_real(k, m, n, draws, loc) result(y)
        real(dp), intent(in) :: k !! observed Type-I count
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects in [0,m]
        real(dp), intent(in) :: draws !! integer sample size in [0,m]
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: count, lower, upper

        count = shifted_count(k, loc)
        if (.not. valid_parameters(m, n, draws) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
            return
        end if
        lower = support_lower(m, n, draws)
        upper = support_upper(n, draws)
        if (.not. ieee_is_finite(count) .or. count /= aint(count) .or. count < lower .or. count > upper) then
            y = negative_infinity(k)
        else
            y = log_choose(n, count) + log_choose(m - n, draws - count) - log_choose(m, draws)
        end if
    end function hypergeom_logpmf_real

    pure elemental function hypergeom_logpmf_int(k, m, n, draws, loc) result(y)
        integer, intent(in) :: k !! observed Type-I count
        integer, intent(in) :: m !! positive population size
        integer, intent(in) :: n !! number of Type-I objects
        integer, intent(in) :: draws !! sample size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = hypergeom_logpmf_real(real(k, dp), real(m, dp), real(n, dp), real(draws, dp), loc)
    end function hypergeom_logpmf_int

    pure elemental function hypergeom_pmf_real(k, m, n, draws, loc) result(y)
        real(dp), intent(in) :: k !! observed Type-I count
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects in [0,m]
        real(dp), intent(in) :: draws !! integer sample size in [0,m]
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = exp(hypergeom_logpmf_real(k, m, n, draws, loc))
    end function hypergeom_pmf_real

    pure elemental function hypergeom_pmf_int(k, m, n, draws, loc) result(y)
        integer, intent(in) :: k !! observed Type-I count
        integer, intent(in) :: m !! positive population size
        integer, intent(in) :: n !! number of Type-I objects
        integer, intent(in) :: draws !! sample size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = hypergeom_pmf_real(real(k, dp), real(m, dp), real(n, dp), real(draws, dp), loc)
    end function hypergeom_pmf_int

    pure elemental function hypergeom_cdf_real(k, m, n, draws, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of P(X <= k)
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: logcdf, logsf, sf
        call tails(floor_count(k, loc), m, n, draws, y, sf, logcdf, logsf)
    end function hypergeom_cdf_real

    pure elemental function hypergeom_cdf_int(k, m, n, draws, loc) result(y)
        integer, intent(in) :: k !! upper integer limit
        integer, intent(in) :: m !! positive population size
        integer, intent(in) :: n !! number of Type-I objects
        integer, intent(in) :: draws !! sample size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = hypergeom_cdf_real(real(k, dp), real(m, dp), real(n, dp), real(draws, dp), loc)
    end function hypergeom_cdf_int

    pure elemental function hypergeom_sf_real(k, m, n, draws, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of P(X > k)
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: cdf, logcdf, logsf
        call tails(floor_count(k, loc), m, n, draws, cdf, y, logcdf, logsf)
    end function hypergeom_sf_real

    pure elemental function hypergeom_sf_int(k, m, n, draws, loc) result(y)
        integer, intent(in) :: k !! lower integer limit
        integer, intent(in) :: m !! positive population size
        integer, intent(in) :: n !! number of Type-I objects
        integer, intent(in) :: draws !! sample size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = hypergeom_sf_real(real(k, dp), real(m, dp), real(n, dp), real(draws, dp), loc)
    end function hypergeom_sf_int

    pure elemental function hypergeom_logcdf_real(k, m, n, draws, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of log(P(X <= k))
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: cdf, sf, logsf
        call tails(floor_count(k, loc), m, n, draws, cdf, sf, y, logsf)
    end function hypergeom_logcdf_real

    pure elemental function hypergeom_logcdf_int(k, m, n, draws, loc) result(y)
        integer, intent(in) :: k !! upper integer limit
        integer, intent(in) :: m !! positive population size
        integer, intent(in) :: n !! number of Type-I objects
        integer, intent(in) :: draws !! sample size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = hypergeom_logcdf_real(real(k, dp), real(m, dp), real(n, dp), real(draws, dp), loc)
    end function hypergeom_logcdf_int

    pure elemental function hypergeom_logsf_real(k, m, n, draws, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of log(P(X > k))
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: cdf, sf, logcdf
        call tails(floor_count(k, loc), m, n, draws, cdf, sf, logcdf, y)
    end function hypergeom_logsf_real

    pure elemental function hypergeom_logsf_int(k, m, n, draws, loc) result(y)
        integer, intent(in) :: k !! lower integer limit
        integer, intent(in) :: m !! positive population size
        integer, intent(in) :: n !! number of Type-I objects
        integer, intent(in) :: draws !! sample size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = hypergeom_logsf_real(real(k, dp), real(m, dp), real(n, dp), real(draws, dp), loc)
    end function hypergeom_logsf_int

    pure elemental function hypergeom_ppf(probability, m, n, draws, loc) result(y)
        real(dp), intent(in) :: probability !! lower-tail probability in [0,1]
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: shift, lower, upper
        integer :: lo, hi, mid

        shift = optional_loc(loc)
        if (.not. valid_parameters(m, n, draws) .or. .not. ieee_is_finite(shift) .or. &
            .not. valid_probability(probability)) then
            y = quiet_nan(probability)
            return
        end if
        lower = support_lower(m, n, draws)
        upper = support_upper(n, draws)
        if (probability <= 0.0_dp) then
            y = shift + lower - 1.0_dp
        else if (probability >= 1.0_dp) then
            y = shift + upper
        else
            lo = int(lower); hi = int(upper)
            do while (lo < hi)
                mid = lo + (hi - lo) / 2
                if (hypergeom_cdf_real(shift + real(mid, dp), m, n, draws, shift) >= probability) then
                    hi = mid
                else
                    lo = mid + 1
                end if
            end do
            y = shift + real(lo, dp)
        end if
    end function hypergeom_ppf

    pure elemental function hypergeom_isf(probability, m, n, draws, loc) result(y)
        real(dp), intent(in) :: probability !! upper-tail probability in [0,1]
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: shift, lower, upper
        integer :: lo, hi, mid

        shift = optional_loc(loc)
        if (.not. valid_parameters(m, n, draws) .or. .not. ieee_is_finite(shift) .or. &
            .not. valid_probability(probability)) then
            y = quiet_nan(probability)
            return
        end if
        lower = support_lower(m, n, draws)
        upper = support_upper(n, draws)
        if (probability >= 1.0_dp) then
            y = shift + lower - 1.0_dp
        else if (probability <= 0.0_dp) then
            y = shift + upper
        else
            lo = int(lower); hi = int(upper)
            do while (lo < hi)
                mid = lo + (hi - lo) / 2
                if (hypergeom_sf_real(shift + real(mid, dp), m, n, draws, shift) <= probability) then
                    hi = mid
                else
                    lo = mid + 1
                end if
            end do
            y = shift + real(lo, dp)
        end if
    end function hypergeom_isf

    pure subroutine tails(k, m, n, draws, cdf, sf, logcdf, logsf)
        real(dp), intent(in) :: k !! shifted integer tail split
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(out) :: cdf !! lower-tail probability
        real(dp), intent(out) :: sf !! upper-tail probability
        real(dp), intent(out) :: logcdf !! logarithm of lower tail
        real(dp), intent(out) :: logsf !! logarithm of upper tail
        real(dp) :: lower, upper
        integer :: ik, ilo, ihi

        if (.not. valid_parameters(m, n, draws) .or. ieee_is_nan(k)) then
            cdf = quiet_nan(k); sf = cdf; logcdf = cdf; logsf = cdf
            return
        end if
        lower = support_lower(m, n, draws)
        upper = support_upper(n, draws)
        if (k < lower) then
            cdf = 0.0_dp; sf = 1.0_dp; logcdf = negative_infinity(k); logsf = 0.0_dp
            return
        end if
        if (.not. ieee_is_finite(k) .or. k >= upper) then
            cdf = 1.0_dp; sf = 0.0_dp; logcdf = 0.0_dp; logsf = negative_infinity(k)
            return
        end if
        ik = int(k); ilo = int(lower); ihi = int(upper)
        logcdf = logsum_range(ilo, ik, m, n, draws)
        logsf = logsum_range(ik + 1, ihi, m, n, draws)
        cdf = exp(logcdf)
        sf = exp(logsf)
    end subroutine tails

    pure function logsum_range(lo, hi, m, n, draws) result(value)
        integer, intent(in) :: lo !! inclusive lower support index
        integer, intent(in) :: hi !! inclusive upper support index
        real(dp), intent(in) :: m !! population size
        real(dp), intent(in) :: n !! Type-I count
        real(dp), intent(in) :: draws !! sample size
        real(dp) :: value, term
        integer :: j
        value = negative_infinity(1.0_dp)
        do j = lo, hi
            term = log_choose(n, real(j, dp)) + log_choose(m - n, draws - real(j, dp)) - log_choose(m, draws)
            value = logadd(value, term)
        end do
    end function logsum_range

    pure elemental function log_choose(total, chosen) result(y)
        real(dp), intent(in) :: total !! nonnegative integer upper argument
        real(dp), intent(in) :: chosen !! integer lower argument in [0,total]
        real(dp) :: y
        y = log_gamma(total + 1.0_dp) - log_gamma(chosen + 1.0_dp) - log_gamma(total - chosen + 1.0_dp)
    end function log_choose

    pure elemental function logadd(a, b) result(y)
        real(dp), intent(in) :: a !! first logarithm
        real(dp), intent(in) :: b !! second logarithm
        real(dp) :: y
        real(dp) :: hi, lo
        hi = max(a, b); lo = min(a, b)
        if (.not. ieee_is_finite(hi)) then
            y = hi
        else
            y = hi + log1p_safe(exp(lo - hi))
        end if
    end function logadd

    pure elemental function support_lower(m, n, draws) result(y)
        real(dp), intent(in) :: m !! population size
        real(dp), intent(in) :: n !! Type-I count
        real(dp), intent(in) :: draws !! sample size
        real(dp) :: y
        y = max(0.0_dp, draws - (m - n))
    end function support_lower

    pure elemental function support_upper(n, draws) result(y)
        real(dp), intent(in) :: n !! Type-I count
        real(dp), intent(in) :: draws !! sample size
        real(dp) :: y
        y = min(n, draws)
    end function support_upper

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
        count = floor_real(count)
    end function floor_count

    pure elemental function optional_loc(loc) result(shift)
        real(dp), intent(in), optional :: loc !! support shift
        real(dp) :: shift
        shift = 0.0_dp
        if (present(loc)) shift = loc
    end function optional_loc

    pure elemental logical function valid_parameters(m, n, draws)
        real(dp), intent(in) :: m !! population size
        real(dp), intent(in) :: n !! Type-I count
        real(dp), intent(in) :: draws !! sample size
        valid_parameters = ieee_is_finite(m) .and. ieee_is_finite(n) .and. ieee_is_finite(draws) .and. &
            m > 0.0_dp .and. n >= 0.0_dp .and. draws >= 0.0_dp .and. n <= m .and. draws <= m .and. &
            m == aint(m) .and. n == aint(n) .and. draws == aint(draws) .and. m <= real(huge(0), dp)
    end function valid_parameters

    pure elemental logical function valid_probability(p)
        real(dp), intent(in) :: p !! probability candidate
        valid_probability = ieee_is_finite(p) .and. p >= 0.0_dp .and. p <= 1.0_dp
    end function valid_probability

end module scifort_hypergeom
