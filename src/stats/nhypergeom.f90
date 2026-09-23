! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Negative hypergeometric distribution matching scipy.stats.nhypergeom.
! X counts successes drawn before the r-th failure when sampling without
! replacement from M objects containing n successes.

module scifort_nhypergeom
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, quiet_nan
    implicit none
    private

    public :: nhypergeom_cdf, nhypergeom_isf, nhypergeom_logcdf, nhypergeom_logpmf
    public :: nhypergeom_logsf, nhypergeom_pmf, nhypergeom_ppf, nhypergeom_sf

    interface nhypergeom_pmf
        module procedure nhypergeom_pmf_real
        module procedure nhypergeom_pmf_int
    end interface nhypergeom_pmf
    interface nhypergeom_logpmf
        module procedure nhypergeom_logpmf_real
        module procedure nhypergeom_logpmf_int
    end interface nhypergeom_logpmf
    interface nhypergeom_cdf
        module procedure nhypergeom_cdf_real
        module procedure nhypergeom_cdf_int
    end interface nhypergeom_cdf
    interface nhypergeom_sf
        module procedure nhypergeom_sf_real
        module procedure nhypergeom_sf_int
    end interface nhypergeom_sf
    interface nhypergeom_logcdf
        module procedure nhypergeom_logcdf_real
        module procedure nhypergeom_logcdf_int
    end interface nhypergeom_logcdf
    interface nhypergeom_logsf
        module procedure nhypergeom_logsf_real
        module procedure nhypergeom_logsf_int
    end interface nhypergeom_logsf

contains

    pure elemental function nhypergeom_logpmf_real(k, m, n, r, loc) result(y)
        real(dp), intent(in) :: k !! number of successes before the r-th failure
        real(dp), intent(in) :: m !! nonnegative integer population size
        real(dp), intent(in) :: n !! integer number of successes in [0,m]
        real(dp), intent(in) :: r !! integer failures required in [0,m-n]
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: count

        count = shifted_count(k, loc)
        if (.not. valid_parameters(m, n, r) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (.not. ieee_is_finite(count) .or. count /= aint(count) .or. &
                 count < 0.0_dp .or. count > n) then
            y = negative_infinity(k)
        else if (r == 0.0_dp) then
            if (count == 0.0_dp) then
                y = 0.0_dp
            else
                y = negative_infinity(k)
            end if
        else
            y = log_choose(count + r - 1.0_dp, count) + &
                log_choose(m - r - count, n - count) - log_choose(m, n)
        end if
    end function nhypergeom_logpmf_real

    pure elemental function nhypergeom_logpmf_int(k, m, n, r, loc) result(y)
        integer, intent(in) :: k !! number of successes before the r-th failure
        integer, intent(in) :: m !! population size
        integer, intent(in) :: n !! number of successes
        integer, intent(in) :: r !! failures required
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = nhypergeom_logpmf_real(real(k, dp), real(m, dp), real(n, dp), real(r, dp), loc)
    end function nhypergeom_logpmf_int

    pure elemental function nhypergeom_pmf_real(k, m, n, r, loc) result(y)
        real(dp), intent(in) :: k !! number of successes before the r-th failure
        real(dp), intent(in) :: m !! nonnegative integer population size
        real(dp), intent(in) :: n !! integer number of successes
        real(dp), intent(in) :: r !! integer failures required
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = exp(nhypergeom_logpmf_real(k, m, n, r, loc))
    end function nhypergeom_pmf_real

    pure elemental function nhypergeom_pmf_int(k, m, n, r, loc) result(y)
        integer, intent(in) :: k !! number of successes before the r-th failure
        integer, intent(in) :: m !! population size
        integer, intent(in) :: n !! number of successes
        integer, intent(in) :: r !! failures required
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = nhypergeom_pmf_real(real(k, dp), real(m, dp), real(n, dp), real(r, dp), loc)
    end function nhypergeom_pmf_int

    pure elemental function nhypergeom_cdf_real(k, m, n, r, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of P(X <= k)
        real(dp), intent(in) :: m !! nonnegative integer population size
        real(dp), intent(in) :: n !! integer number of successes
        real(dp), intent(in) :: r !! integer failures required
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: logcdf, logsf, sf
        call tails(floor_count(k, loc), m, n, r, y, sf, logcdf, logsf)
    end function nhypergeom_cdf_real

    pure elemental function nhypergeom_cdf_int(k, m, n, r, loc) result(y)
        integer, intent(in) :: k !! upper integer limit
        integer, intent(in) :: m !! population size
        integer, intent(in) :: n !! number of successes
        integer, intent(in) :: r !! failures required
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = nhypergeom_cdf_real(real(k, dp), real(m, dp), real(n, dp), real(r, dp), loc)
    end function nhypergeom_cdf_int

    pure elemental function nhypergeom_sf_real(k, m, n, r, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of P(X > k)
        real(dp), intent(in) :: m !! nonnegative integer population size
        real(dp), intent(in) :: n !! integer number of successes
        real(dp), intent(in) :: r !! integer failures required
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: cdf, logcdf, logsf
        call tails(floor_count(k, loc), m, n, r, cdf, y, logcdf, logsf)
    end function nhypergeom_sf_real

    pure elemental function nhypergeom_sf_int(k, m, n, r, loc) result(y)
        integer, intent(in) :: k !! lower integer limit
        integer, intent(in) :: m !! population size
        integer, intent(in) :: n !! number of successes
        integer, intent(in) :: r !! failures required
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = nhypergeom_sf_real(real(k, dp), real(m, dp), real(n, dp), real(r, dp), loc)
    end function nhypergeom_sf_int

    pure elemental function nhypergeom_logcdf_real(k, m, n, r, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of log(P(X <= k))
        real(dp), intent(in) :: m !! nonnegative integer population size
        real(dp), intent(in) :: n !! integer number of successes
        real(dp), intent(in) :: r !! integer failures required
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: cdf, sf, logsf
        call tails(floor_count(k, loc), m, n, r, cdf, sf, y, logsf)
    end function nhypergeom_logcdf_real

    pure elemental function nhypergeom_logcdf_int(k, m, n, r, loc) result(y)
        integer, intent(in) :: k !! upper integer limit
        integer, intent(in) :: m !! population size
        integer, intent(in) :: n !! number of successes
        integer, intent(in) :: r !! failures required
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = nhypergeom_logcdf_real(real(k, dp), real(m, dp), real(n, dp), real(r, dp), loc)
    end function nhypergeom_logcdf_int

    pure elemental function nhypergeom_logsf_real(k, m, n, r, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of log(P(X > k))
        real(dp), intent(in) :: m !! nonnegative integer population size
        real(dp), intent(in) :: n !! integer number of successes
        real(dp), intent(in) :: r !! integer failures required
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: cdf, sf, logcdf
        call tails(floor_count(k, loc), m, n, r, cdf, sf, logcdf, y)
    end function nhypergeom_logsf_real

    pure elemental function nhypergeom_logsf_int(k, m, n, r, loc) result(y)
        integer, intent(in) :: k !! lower integer limit
        integer, intent(in) :: m !! population size
        integer, intent(in) :: n !! number of successes
        integer, intent(in) :: r !! failures required
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = nhypergeom_logsf_real(real(k, dp), real(m, dp), real(n, dp), real(r, dp), loc)
    end function nhypergeom_logsf_int

    pure elemental function nhypergeom_ppf(probability, m, n, r, loc) result(y)
        real(dp), intent(in) :: probability !! lower-tail probability in [0,1]
        real(dp), intent(in) :: m !! nonnegative integer population size
        real(dp), intent(in) :: n !! integer number of successes
        real(dp), intent(in) :: r !! integer failures required
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: shift
        integer :: lo, hi, mid

        shift = optional_loc(loc)
        if (.not. valid_parameters(m, n, r) .or. .not. ieee_is_finite(shift) .or. &
            .not. valid_probability(probability)) then
            y = quiet_nan(probability)
        else if (probability <= 0.0_dp) then
            y = shift - 1.0_dp
        else if (probability >= 1.0_dp) then
            y = shift + n
        else
            lo = 0; hi = int(n)
            do while (lo < hi)
                mid = lo + (hi - lo) / 2
                if (nhypergeom_cdf_real(shift + real(mid, dp), m, n, r, shift) >= probability) then
                    hi = mid
                else
                    lo = mid + 1
                end if
            end do
            y = shift + real(lo, dp)
        end if
    end function nhypergeom_ppf

    pure elemental function nhypergeom_isf(probability, m, n, r, loc) result(y)
        real(dp), intent(in) :: probability !! upper-tail probability in [0,1]
        real(dp), intent(in) :: m !! nonnegative integer population size
        real(dp), intent(in) :: n !! integer number of successes
        real(dp), intent(in) :: r !! integer failures required
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: shift
        integer :: lo, hi, mid

        shift = optional_loc(loc)
        if (.not. valid_parameters(m, n, r) .or. .not. ieee_is_finite(shift) .or. &
            .not. valid_probability(probability)) then
            y = quiet_nan(probability)
        else if (probability >= 1.0_dp) then
            y = shift - 1.0_dp
        else if (probability <= 0.0_dp) then
            y = shift + n
        else
            lo = 0; hi = int(n)
            do while (lo < hi)
                mid = lo + (hi - lo) / 2
                if (nhypergeom_sf_real(shift + real(mid, dp), m, n, r, shift) <= probability) then
                    hi = mid
                else
                    lo = mid + 1
                end if
            end do
            y = shift + real(lo, dp)
        end if
    end function nhypergeom_isf

    pure subroutine tails(k, m, n, r, cdf, sf, logcdf, logsf)
        real(dp), intent(in) :: k !! shifted integer tail split
        real(dp), intent(in) :: m !! population size
        real(dp), intent(in) :: n !! success count
        real(dp), intent(in) :: r !! failures required
        real(dp), intent(out) :: cdf !! lower-tail probability
        real(dp), intent(out) :: sf !! upper-tail probability
        real(dp), intent(out) :: logcdf !! logarithm of lower tail
        real(dp), intent(out) :: logsf !! logarithm of upper tail
        integer :: ik, in

        if (.not. valid_parameters(m, n, r) .or. ieee_is_nan(k)) then
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
        ik = int(k); in = int(n)
        logcdf = logsum_range(0, ik, m, n, r)
        logsf = logsum_range(ik + 1, in, m, n, r)
        cdf = exp(logcdf)
        sf = exp(logsf)
    end subroutine tails

    pure function logsum_range(lo, hi, m, n, r) result(value)
        integer, intent(in) :: lo !! inclusive lower support index
        integer, intent(in) :: hi !! inclusive upper support index
        real(dp), intent(in) :: m !! population size
        real(dp), intent(in) :: n !! success count
        real(dp), intent(in) :: r !! failures required
        real(dp) :: value, term
        integer :: j
        value = negative_infinity(1.0_dp)
        do j = lo, hi
            term = nhypergeom_logpmf_real(real(j, dp), m, n, r)
            value = logadd(value, term)
        end do
    end function logsum_range

    pure elemental function log_choose(total, chosen) result(y)
        real(dp), intent(in) :: total !! nonnegative integer upper argument
        real(dp), intent(in) :: chosen !! integer lower argument
        real(dp) :: y
        if (chosen < 0.0_dp .or. chosen > total) then
            y = negative_infinity(total)
        else
            y = log_gamma(total + 1.0_dp) - log_gamma(chosen + 1.0_dp) - log_gamma(total - chosen + 1.0_dp)
        end if
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

    pure elemental logical function valid_parameters(m, n, r)
        real(dp), intent(in) :: m !! population size
        real(dp), intent(in) :: n !! success count
        real(dp), intent(in) :: r !! failures required
        valid_parameters = ieee_is_finite(m) .and. ieee_is_finite(n) .and. ieee_is_finite(r) .and. &
            m >= 0.0_dp .and. n >= 0.0_dp .and. n <= m .and. r >= 0.0_dp .and. r <= m - n .and. &
            m == aint(m) .and. n == aint(n) .and. r == aint(r) .and. m <= real(huge(0), dp)
    end function valid_parameters

    pure elemental logical function valid_probability(p)
        real(dp), intent(in) :: p !! probability candidate
        valid_probability = ieee_is_finite(p) .and. p >= 0.0_dp .and. p <= 1.0_dp
    end function valid_probability

end module scifort_nhypergeom
