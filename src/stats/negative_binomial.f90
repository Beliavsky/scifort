! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Negative binomial distribution matching scipy.stats.nbinom. X counts the
! failures observed before n successes, where n > 0 need not be an integer.
! Its mass is Gamma(k+n)/(Gamma(n) Gamma(k+1)) p**n (1-p)**k.
!
! DLMF 8.17.5 gives the cumulative relation
! P(X <= k) = I_p(n, k + 1). The complementary beta tail is evaluated in the
! same call, so survival probabilities are not obtained by subtracting from 1.

module scifort_negative_binomial
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_incomplete_beta, only : incomplete_beta_xy, log_beta_kernel
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, quiet_nan
    use scifort_normal, only : normal_isf, normal_ppf
    implicit none
    private

    public :: negative_binomial_cdf
    public :: negative_binomial_isf
    public :: negative_binomial_logcdf
    public :: negative_binomial_logpmf
    public :: negative_binomial_logsf
    public :: negative_binomial_pmf
    public :: negative_binomial_ppf
    public :: negative_binomial_sf

    interface negative_binomial_pmf
        module procedure negative_binomial_pmf_real
        module procedure negative_binomial_pmf_int
    end interface negative_binomial_pmf

    interface negative_binomial_logpmf
        module procedure negative_binomial_logpmf_real
        module procedure negative_binomial_logpmf_int
    end interface negative_binomial_logpmf

    interface negative_binomial_cdf
        module procedure negative_binomial_cdf_real
        module procedure negative_binomial_cdf_int
    end interface negative_binomial_cdf

    interface negative_binomial_sf
        module procedure negative_binomial_sf_real
        module procedure negative_binomial_sf_int
    end interface negative_binomial_sf

    interface negative_binomial_logcdf
        module procedure negative_binomial_logcdf_real
        module procedure negative_binomial_logcdf_int
    end interface negative_binomial_logcdf

    interface negative_binomial_logsf
        module procedure negative_binomial_logsf_real
        module procedure negative_binomial_logsf_int
    end interface negative_binomial_logsf

    integer, parameter :: max_search_steps = 200

contains

    pure elemental function negative_binomial_pmf_real(k, n, p, loc) result(y)
        real(dp), intent(in) :: k !! failure count; non-integer values have probability zero
        real(dp), intent(in) :: n !! required successes, finite and > 0
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        y = exp(negative_binomial_logpmf_real(k, n, p, loc))
    end function negative_binomial_pmf_real

    pure elemental function negative_binomial_pmf_int(k, n, p, loc) result(y)
        integer, intent(in) :: k !! failure count
        real(dp), intent(in) :: n !! required successes, finite and > 0
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        y = negative_binomial_pmf_real(real(k, dp), n, p, loc)
    end function negative_binomial_pmf_int

    pure elemental function negative_binomial_logpmf_real(k, n, p, loc) result(y)
        real(dp), intent(in) :: k !! failure count; non-integer values have probability zero
        real(dp), intent(in) :: n !! required successes, finite and > 0
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        real(dp) :: count

        count = shifted_count(k, loc)
        if (.not. valid_parameters(n, p) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (count < 0.0_dp .or. .not. ieee_is_finite(count)) then
            y = negative_infinity(k)
        else if (count /= aint(count)) then
            y = negative_infinity(k)
        else if (p == 1.0_dp) then
            if (count == 0.0_dp) then
                y = 0.0_dp
            else
                y = negative_infinity(k)
            end if
        else if (count == 0.0_dp) then
            y = n * log(p)
        else
            y = log_beta_kernel(n, count + 1.0_dp, p, 1.0_dp - p) - &
                log_sum_positive(n, count) - log1p_safe(-p)
        end if
    end function negative_binomial_logpmf_real

    pure elemental function negative_binomial_logpmf_int(k, n, p, loc) result(y)
        integer, intent(in) :: k !! failure count
        real(dp), intent(in) :: n !! required successes, finite and > 0
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        y = negative_binomial_logpmf_real(real(k, dp), n, p, loc)
    end function negative_binomial_logpmf_int

    pure elemental function negative_binomial_cdf_real(k, n, p, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of P(X <= k); the floor is used
        real(dp), intent(in) :: n !! required successes, finite and > 0
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        real(dp) :: logcdf
        real(dp) :: logsf
        real(dp) :: sf

        call tails(floor_count(k, loc), n, p, y, sf, logcdf, logsf)
    end function negative_binomial_cdf_real

    pure elemental function negative_binomial_cdf_int(k, n, p, loc) result(y)
        integer, intent(in) :: k !! upper limit of P(X <= k)
        real(dp), intent(in) :: n !! required successes, finite and > 0
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        y = negative_binomial_cdf_real(real(k, dp), n, p, loc)
    end function negative_binomial_cdf_int

    pure elemental function negative_binomial_sf_real(k, n, p, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of P(X > k); the floor is used
        real(dp), intent(in) :: n !! required successes, finite and > 0
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        real(dp) :: cdf
        real(dp) :: logcdf
        real(dp) :: logsf

        call tails(floor_count(k, loc), n, p, cdf, y, logcdf, logsf)
    end function negative_binomial_sf_real

    pure elemental function negative_binomial_sf_int(k, n, p, loc) result(y)
        integer, intent(in) :: k !! lower limit of P(X > k)
        real(dp), intent(in) :: n !! required successes, finite and > 0
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        y = negative_binomial_sf_real(real(k, dp), n, p, loc)
    end function negative_binomial_sf_int

    pure elemental function negative_binomial_logcdf_real(k, n, p, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of log(P(X <= k)); the floor is used
        real(dp), intent(in) :: n !! required successes, finite and > 0
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        real(dp) :: cdf
        real(dp) :: logsf
        real(dp) :: sf

        call tails(floor_count(k, loc), n, p, cdf, sf, y, logsf)
    end function negative_binomial_logcdf_real

    pure elemental function negative_binomial_logcdf_int(k, n, p, loc) result(y)
        integer, intent(in) :: k !! upper limit of log(P(X <= k))
        real(dp), intent(in) :: n !! required successes, finite and > 0
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        y = negative_binomial_logcdf_real(real(k, dp), n, p, loc)
    end function negative_binomial_logcdf_int

    pure elemental function negative_binomial_logsf_real(k, n, p, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of log(P(X > k)); the floor is used
        real(dp), intent(in) :: n !! required successes, finite and > 0
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        real(dp) :: cdf
        real(dp) :: logcdf
        real(dp) :: sf

        call tails(floor_count(k, loc), n, p, cdf, sf, logcdf, y)
    end function negative_binomial_logsf_real

    pure elemental function negative_binomial_logsf_int(k, n, p, loc) result(y)
        integer, intent(in) :: k !! lower limit of log(P(X > k))
        real(dp), intent(in) :: n !! required successes, finite and > 0
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        y = negative_binomial_logsf_real(real(k, dp), n, p, loc)
    end function negative_binomial_logsf_int

    pure elemental function negative_binomial_ppf(probability, n, p, loc) result(y)
        real(dp), intent(in) :: probability !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: n !! required successes, finite and > 0
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        real(dp) :: shift

        shift = 0.0_dp
        if (present(loc)) shift = loc
        if (.not. valid_parameters(n, p) .or. .not. ieee_is_finite(shift)) then
            y = quiet_nan(probability)
        else if (.not. (probability >= 0.0_dp .and. probability <= 1.0_dp)) then
            y = quiet_nan(probability)
        else if (probability == 0.0_dp) then
            y = shift - 1.0_dp
        else if (probability >= 1.0_dp) then
            y = positive_infinity(probability)
        else if (p == 1.0_dp) then
            y = shift
        else
            y = shift + search_quantile(probability, n, p, .true.)
        end if
    end function negative_binomial_ppf

    pure elemental function negative_binomial_isf(probability, n, p, loc) result(y)
        real(dp), intent(in) :: probability !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: n !! required successes, finite and > 0
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        real(dp) :: shift

        shift = 0.0_dp
        if (present(loc)) shift = loc
        if (.not. valid_parameters(n, p) .or. .not. ieee_is_finite(shift)) then
            y = quiet_nan(probability)
        else if (.not. (probability >= 0.0_dp .and. probability <= 1.0_dp)) then
            y = quiet_nan(probability)
        else if (probability >= 1.0_dp) then
            y = shift - 1.0_dp
        else if (probability == 0.0_dp) then
            y = positive_infinity(probability)
        else if (p == 1.0_dp) then
            y = shift
        else
            y = shift + search_quantile(probability, n, p, .false.)
        end if
    end function negative_binomial_isf

    pure elemental subroutine tails(count, n, p, cdf, sf, logcdf, logsf)
        real(dp), intent(in) :: count !! floor(k - loc)
        real(dp), intent(in) :: n !! required successes, finite and > 0
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(out) :: cdf !! P(X <= count)
        real(dp), intent(out) :: sf !! P(X > count)
        real(dp), intent(out) :: logcdf !! log(cdf)
        real(dp), intent(out) :: logsf !! log(sf)

        if (.not. valid_parameters(n, p) .or. ieee_is_nan(count)) then
            cdf = quiet_nan(p)
            sf = cdf
            logcdf = cdf
            logsf = cdf
        else if (count < 0.0_dp) then
            cdf = 0.0_dp
            sf = 1.0_dp
            logcdf = negative_infinity(p)
            logsf = 0.0_dp
        else if (.not. ieee_is_finite(count)) then
            cdf = 1.0_dp
            sf = 0.0_dp
            logcdf = 0.0_dp
            logsf = negative_infinity(p)
        else if (p == 1.0_dp) then
            cdf = 1.0_dp
            sf = 0.0_dp
            logcdf = 0.0_dp
            logsf = negative_infinity(p)
        else
            call incomplete_beta_xy(n, count + 1.0_dp, p, 1.0_dp - p, cdf, sf, &
                logcdf, logsf)
        end if
    end subroutine tails

    pure function search_quantile(target, n, p, lower) result(count)
        real(dp), intent(in) :: target !! target probability, strictly between 0 and 1
        real(dp), intent(in) :: n !! required successes, finite and > 0
        real(dp), intent(in) :: p !! success probability, strictly between 0 and 1
        logical, intent(in) :: lower !! true for the lower tail, false for the upper tail
        real(dp) :: count

        integer :: step
        real(dp) :: approximation
        real(dp) :: hi
        real(dp) :: lo
        real(dp) :: mean
        real(dp) :: middle
        real(dp) :: sd
        real(dp) :: width
        real(dp) :: z

        mean = n * (1.0_dp - p) / p
        sd = sqrt(n * (1.0_dp - p)) / p
        if (.not. ieee_is_finite(mean) .or. .not. ieee_is_finite(sd)) then
            count = positive_infinity(target)
            return
        end if

        if (lower) then
            z = normal_ppf(target)
        else
            z = normal_isf(target)
        end if
        approximation = mean + z * sd + 0.5_dp
        if (ieee_is_finite(approximation)) then
            count = max(0.0_dp, aint(approximation))
        else if (approximation < 0.0_dp) then
            count = 0.0_dp
        else
            count = positive_infinity(target)
            return
        end if

        if (satisfies(count, target, n, p, lower)) then
            hi = count
            width = 1.0_dp
            lo = max(-1.0_dp, count - width)
            do step = 1, max_search_steps
                if (lo < 0.0_dp) exit
                if (.not. satisfies(lo, target, n, p, lower)) exit
                hi = lo
                width = 2.0_dp * width
                lo = max(-1.0_dp, lo - width)
            end do
        else
            lo = count
            width = 1.0_dp
            hi = count + width
            do step = 1, max_search_steps
                if (.not. ieee_is_finite(hi)) then
                    count = positive_infinity(target)
                    return
                end if
                if (satisfies(hi, target, n, p, lower)) exit
                lo = hi
                width = 2.0_dp * width
                hi = hi + width
            end do
        end if

        do step = 1, max_search_steps
            if (hi - lo <= 1.0_dp) exit
            middle = aint(lo + 0.5_dp * (hi - lo))
            if (satisfies(middle, target, n, p, lower)) then
                hi = middle
            else
                lo = middle
            end if
        end do
        count = hi
    end function search_quantile

    pure function satisfies(count, target, n, p, lower) result(ok)
        real(dp), intent(in) :: count !! candidate failure count
        real(dp), intent(in) :: target !! target probability
        real(dp), intent(in) :: n !! required successes, finite and > 0
        real(dp), intent(in) :: p !! success probability, strictly between 0 and 1
        logical, intent(in) :: lower !! true for the lower tail, false for the upper tail
        logical :: ok

        real(dp) :: cdf
        real(dp) :: logcdf
        real(dp) :: logsf
        real(dp) :: sf

        if (count < 0.0_dp) then
            ok = .false.
            return
        end if
        call tails(count, n, p, cdf, sf, logcdf, logsf)
        if (lower) then
            if (target <= 0.5_dp) then
                ok = cdf >= target
            else
                ok = sf <= 1.0_dp - target
            end if
        else
            if (target <= 0.5_dp) then
                ok = sf <= target
            else
                ok = cdf >= 1.0_dp - target
            end if
        end if
    end function satisfies

    pure elemental function log_sum_positive(a, b) result(y)
        real(dp), intent(in) :: a !! first positive term
        real(dp), intent(in) :: b !! second nonnegative term
        real(dp) :: y

        if (a >= b) then
            y = log(a) + log1p_safe(b / a)
        else
            y = log(b) + log1p_safe(a / b)
        end if
    end function log_sum_positive

    pure elemental function shifted_count(k, loc) result(count)
        real(dp), intent(in) :: k !! argument in the shifted support
        real(dp), intent(in), optional :: loc !! shift of the support
        real(dp) :: count

        if (present(loc)) then
            count = k - loc
        else
            count = k
        end if
    end function shifted_count

    pure elemental function floor_count(k, loc) result(count)
        real(dp), intent(in) :: k !! argument in the shifted support
        real(dp), intent(in), optional :: loc !! shift of the support
        real(dp) :: count

        count = shifted_count(k, loc)
        if (ieee_is_finite(count)) count = floor(count)
    end function floor_count

    pure elemental logical function valid_parameters(n, p) result(valid)
        real(dp), intent(in) :: n !! required successes to check
        real(dp), intent(in) :: p !! success probability to check

        valid = ieee_is_finite(n) .and. n > 0.0_dp .and. p > 0.0_dp .and. p <= 1.0_dp
    end function valid_parameters

end module scifort_negative_binomial
