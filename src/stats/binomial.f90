! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Binomial distribution for n independent trials with success probability p
! and integer location loc, matching scipy.stats.binom. The support is
! loc, loc + 1, ..., loc + n and a non-integer argument has probability zero.
!
! DLMF 8.17.5 states I_x(m, n - m + 1) = sum over j = m..n of
! C(n, j) x**j (1 - x)**(n - j), so the upper tail is
! P(X > k) = I_p(k + 1, n - k) and the lower tail is its complement. Both are
! taken from one call to scifort_incomplete_beta, which forms each tail
! directly rather than by subtraction. The mass function uses the stable
! kernel x**a y**b / B(a, b) of the same module, since
! C(n, k) p**k q**(n - k) = p**(k+1) q**(n-k+1) /
! ((n + 1) p q B(k + 1, n - k + 1)).

module scifort_binomial
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_incomplete_beta, only : incomplete_beta_xy, log_beta_kernel
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, quiet_nan
    use scifort_normal, only : normal_isf, normal_ppf
    implicit none
    private

    public :: binomial_cdf
    public :: binomial_isf
    public :: binomial_logcdf
    public :: binomial_logpmf
    public :: binomial_logsf
    public :: binomial_pmf
    public :: binomial_ppf
    public :: binomial_sf

    !> Probability mass function.
    interface binomial_pmf
        module procedure binomial_pmf_real
        module procedure binomial_pmf_int
    end interface binomial_pmf

    !> Logarithm of the probability mass function.
    interface binomial_logpmf
        module procedure binomial_logpmf_real
        module procedure binomial_logpmf_int
    end interface binomial_logpmf

    !> Cumulative distribution function P(X <= k).
    interface binomial_cdf
        module procedure binomial_cdf_real
        module procedure binomial_cdf_int
    end interface binomial_cdf

    !> Survival function P(X > k).
    interface binomial_sf
        module procedure binomial_sf_real
        module procedure binomial_sf_int
    end interface binomial_sf

    !> Logarithm of the cumulative distribution function.
    interface binomial_logcdf
        module procedure binomial_logcdf_real
        module procedure binomial_logcdf_int
    end interface binomial_logcdf

    !> Logarithm of the survival function.
    interface binomial_logsf
        module procedure binomial_logsf_real
        module procedure binomial_logsf_int
    end interface binomial_logsf

    integer, parameter :: max_search_steps = 200

contains

    pure elemental function binomial_pmf_real(k, n, p, loc) result(y)
        real(dp), intent(in) :: k !! number of successes; non-integer values have probability zero
        real(dp), intent(in) :: n !! number of trials, a nonnegative integer
        real(dp), intent(in) :: p !! success probability in [0, 1]
        real(dp), intent(in), optional :: loc !! integer shift of the support (default 0)
        real(dp) :: y

        y = exp(binomial_logpmf_real(k, n, p, loc))
    end function binomial_pmf_real

    pure elemental function binomial_pmf_int(k, n, p, loc) result(y)
        integer, intent(in) :: k !! number of successes
        integer, intent(in) :: n !! number of trials, >= 0
        real(dp), intent(in) :: p !! success probability in [0, 1]
        real(dp), intent(in), optional :: loc !! integer shift of the support (default 0)
        real(dp) :: y

        y = binomial_pmf_real(real(k, dp), real(n, dp), p, loc)
    end function binomial_pmf_int

    pure elemental function binomial_logpmf_real(k, n, p, loc) result(y)
        real(dp), intent(in) :: k !! number of successes; non-integer values have probability zero
        real(dp), intent(in) :: n !! number of trials, a nonnegative integer
        real(dp), intent(in) :: p !! success probability in [0, 1]
        real(dp), intent(in), optional :: loc !! integer shift of the support (default 0)
        real(dp) :: y

        real(dp) :: count
        real(dp) :: q

        count = shifted_count(k, loc)
        q = 1.0_dp - p
        if (.not. valid_parameters(n, p) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (count < 0.0_dp .or. count > n .or. .not. ieee_is_finite(count)) then
            y = negative_infinity(k)
        else if (count /= aint(count)) then
            y = negative_infinity(k)
        else if (p == 0.0_dp) then
            y = log_indicator(count == 0.0_dp, k)
        else if (p == 1.0_dp) then
            y = log_indicator(count == n, k)
        else if (count == 0.0_dp) then
            y = n * log1p_safe(-p)
        else if (count == n) then
            y = n * log(p)
        else
            y = log_beta_kernel(count + 1.0_dp, n - count + 1.0_dp, p, q) - &
                log(n + 1.0_dp) - log(p) - log(q)
        end if
    end function binomial_logpmf_real

    pure elemental function binomial_logpmf_int(k, n, p, loc) result(y)
        integer, intent(in) :: k !! number of successes
        integer, intent(in) :: n !! number of trials, >= 0
        real(dp), intent(in) :: p !! success probability in [0, 1]
        real(dp), intent(in), optional :: loc !! integer shift of the support (default 0)
        real(dp) :: y

        y = binomial_logpmf_real(real(k, dp), real(n, dp), p, loc)
    end function binomial_logpmf_int

    pure elemental function binomial_cdf_real(k, n, p, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of P(X <= k); the floor is used
        real(dp), intent(in) :: n !! number of trials, a nonnegative integer
        real(dp), intent(in) :: p !! success probability in [0, 1]
        real(dp), intent(in), optional :: loc !! integer shift of the support (default 0)
        real(dp) :: y

        real(dp) :: logcdf
        real(dp) :: logsf
        real(dp) :: sf

        call tails(floor_count(k, loc), n, p, y, sf, logcdf, logsf)
    end function binomial_cdf_real

    pure elemental function binomial_cdf_int(k, n, p, loc) result(y)
        integer, intent(in) :: k !! upper limit of P(X <= k)
        integer, intent(in) :: n !! number of trials, >= 0
        real(dp), intent(in) :: p !! success probability in [0, 1]
        real(dp), intent(in), optional :: loc !! integer shift of the support (default 0)
        real(dp) :: y

        y = binomial_cdf_real(real(k, dp), real(n, dp), p, loc)
    end function binomial_cdf_int

    pure elemental function binomial_sf_real(k, n, p, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of P(X > k); the floor is used
        real(dp), intent(in) :: n !! number of trials, a nonnegative integer
        real(dp), intent(in) :: p !! success probability in [0, 1]
        real(dp), intent(in), optional :: loc !! integer shift of the support (default 0)
        real(dp) :: y

        real(dp) :: cdf
        real(dp) :: logcdf
        real(dp) :: logsf

        call tails(floor_count(k, loc), n, p, cdf, y, logcdf, logsf)
    end function binomial_sf_real

    pure elemental function binomial_sf_int(k, n, p, loc) result(y)
        integer, intent(in) :: k !! lower limit of P(X > k)
        integer, intent(in) :: n !! number of trials, >= 0
        real(dp), intent(in) :: p !! success probability in [0, 1]
        real(dp), intent(in), optional :: loc !! integer shift of the support (default 0)
        real(dp) :: y

        y = binomial_sf_real(real(k, dp), real(n, dp), p, loc)
    end function binomial_sf_int

    pure elemental function binomial_logcdf_real(k, n, p, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of log(P(X <= k)); the floor is used
        real(dp), intent(in) :: n !! number of trials, a nonnegative integer
        real(dp), intent(in) :: p !! success probability in [0, 1]
        real(dp), intent(in), optional :: loc !! integer shift of the support (default 0)
        real(dp) :: y

        real(dp) :: cdf
        real(dp) :: logsf
        real(dp) :: sf

        call tails(floor_count(k, loc), n, p, cdf, sf, y, logsf)
    end function binomial_logcdf_real

    pure elemental function binomial_logcdf_int(k, n, p, loc) result(y)
        integer, intent(in) :: k !! upper limit of log(P(X <= k))
        integer, intent(in) :: n !! number of trials, >= 0
        real(dp), intent(in) :: p !! success probability in [0, 1]
        real(dp), intent(in), optional :: loc !! integer shift of the support (default 0)
        real(dp) :: y

        y = binomial_logcdf_real(real(k, dp), real(n, dp), p, loc)
    end function binomial_logcdf_int

    pure elemental function binomial_logsf_real(k, n, p, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of log(P(X > k)); the floor is used
        real(dp), intent(in) :: n !! number of trials, a nonnegative integer
        real(dp), intent(in) :: p !! success probability in [0, 1]
        real(dp), intent(in), optional :: loc !! integer shift of the support (default 0)
        real(dp) :: y

        real(dp) :: cdf
        real(dp) :: logcdf
        real(dp) :: sf

        call tails(floor_count(k, loc), n, p, cdf, sf, logcdf, y)
    end function binomial_logsf_real

    pure elemental function binomial_logsf_int(k, n, p, loc) result(y)
        integer, intent(in) :: k !! lower limit of log(P(X > k))
        integer, intent(in) :: n !! number of trials, >= 0
        real(dp), intent(in) :: p !! success probability in [0, 1]
        real(dp), intent(in), optional :: loc !! integer shift of the support (default 0)
        real(dp) :: y

        y = binomial_logsf_real(real(k, dp), real(n, dp), p, loc)
    end function binomial_logsf_int

    ! Smallest k in the support with P(X <= k) >= probability. As in
    ! scipy.stats, a probability of zero gives loc - 1.
    pure elemental function binomial_ppf(probability, n, p, loc) result(y)
        real(dp), intent(in) :: probability !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: n !! number of trials, a nonnegative integer
        real(dp), intent(in) :: p !! success probability in [0, 1]
        real(dp), intent(in), optional :: loc !! integer shift of the support (default 0)
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
            y = shift + n
        else
            y = shift + search_quantile(probability, n, p, .true.)
        end if
    end function binomial_ppf

    ! Smallest k in the support with P(X > k) <= probability.
    pure elemental function binomial_isf(probability, n, p, loc) result(y)
        real(dp), intent(in) :: probability !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: n !! number of trials, a nonnegative integer
        real(dp), intent(in) :: p !! success probability in [0, 1]
        real(dp), intent(in), optional :: loc !! integer shift of the support (default 0)
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
            y = shift + n
        else
            y = shift + search_quantile(probability, n, p, .false.)
        end if
    end function binomial_isf

    ! Both tails and their logarithms at a floored count.
    pure elemental subroutine tails(count, n, p, cdf, sf, logcdf, logsf)
        real(dp), intent(in) :: count !! floor(k - loc)
        real(dp), intent(in) :: n !! number of trials, a nonnegative integer
        real(dp), intent(in) :: p !! success probability in [0, 1]
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
            call set_tails(0.0_dp, 1.0_dp, negative_infinity(p), 0.0_dp, cdf, sf, &
                logcdf, logsf)
        else if (count >= n) then
            call set_tails(1.0_dp, 0.0_dp, 0.0_dp, negative_infinity(p), cdf, sf, &
                logcdf, logsf)
        else if (p == 0.0_dp) then
            call set_tails(1.0_dp, 0.0_dp, 0.0_dp, negative_infinity(p), cdf, sf, &
                logcdf, logsf)
        else if (p == 1.0_dp) then
            call set_tails(0.0_dp, 1.0_dp, negative_infinity(p), 0.0_dp, cdf, sf, &
                logcdf, logsf)
        else
            ! P(X > count) = I_p(count + 1, n - count) by DLMF 8.17.5.
            call incomplete_beta_xy(count + 1.0_dp, n - count, p, 1.0_dp - p, sf, cdf, &
                logsf, logcdf)
        end if
    end subroutine tails

    pure elemental subroutine set_tails(cdf_in, sf_in, logcdf_in, logsf_in, cdf, sf, &
            logcdf, logsf)
        real(dp), intent(in) :: cdf_in !! value for the lower tail
        real(dp), intent(in) :: sf_in !! value for the upper tail
        real(dp), intent(in) :: logcdf_in !! value for the log lower tail
        real(dp), intent(in) :: logsf_in !! value for the log upper tail
        real(dp), intent(out) :: cdf !! P(X <= count)
        real(dp), intent(out) :: sf !! P(X > count)
        real(dp), intent(out) :: logcdf !! log(cdf)
        real(dp), intent(out) :: logsf !! log(sf)

        cdf = cdf_in
        sf = sf_in
        logcdf = logcdf_in
        logsf = logsf_in
    end subroutine set_tails

    ! Smallest count in [0, n] with P(X <= k) >= target (lower) or
    ! P(X > k) <= target (upper). A normal approximation with a continuity
    ! correction starts the search, a bracket is found by doubling, and the
    ! bracket is then halved.
    pure function search_quantile(target, n, p, lower) result(count)
        real(dp), intent(in) :: target !! target probability, 0 < target < 1
        real(dp), intent(in) :: n !! number of trials, a nonnegative integer
        real(dp), intent(in) :: p !! success probability, 0 < p < 1
        logical, intent(in) :: lower !! .true. for the lower tail, .false. for the upper
        real(dp) :: count

        integer :: step
        real(dp) :: hi
        real(dp) :: lo
        real(dp) :: middle
        real(dp) :: width
        real(dp) :: z

        if (lower) then
            z = normal_ppf(target)
        else
            z = normal_isf(target)
        end if
        count = min(n, max(0.0_dp, aint(n * p + z * sqrt(n * p * (1.0_dp - p)) + 0.5_dp)))

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
            hi = min(n, count + width)
            do step = 1, max_search_steps
                if (satisfies(hi, target, n, p, lower)) exit
                lo = hi
                width = 2.0_dp * width
                hi = min(n, hi + width)
            end do
        end if

        ! Invariant: lo does not satisfy the condition, hi does.
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
        real(dp), intent(in) :: count !! candidate count
        real(dp), intent(in) :: target !! target probability
        real(dp), intent(in) :: n !! number of trials, a nonnegative integer
        real(dp), intent(in) :: p !! success probability, 0 < p < 1
        logical, intent(in) :: lower !! .true. for the lower tail, .false. for the upper
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
            ok = cdf >= target
        else
            ok = sf <= target
        end if
    end function satisfies

    ! k - loc, with loc defaulting to zero.
    pure elemental function shifted_count(k, loc) result(count)
        real(dp), intent(in) :: k !! argument in the shifted support
        real(dp), intent(in), optional :: loc !! integer shift of the support
        real(dp) :: count

        if (present(loc)) then
            count = k - loc
        else
            count = k
        end if
    end function shifted_count

    ! floor(k - loc), which is where a discrete probability changes value.
    pure elemental function floor_count(k, loc) result(count)
        real(dp), intent(in) :: k !! argument in the shifted support
        real(dp), intent(in), optional :: loc !! integer shift of the support
        real(dp) :: count

        count = shifted_count(k, loc)
        if (ieee_is_finite(count)) count = floor(count)
    end function floor_count

    ! 0 for a true condition and -infinity for a false one.
    pure elemental function log_indicator(condition, x) result(y)
        logical, intent(in) :: condition !! whether the outcome has probability one
        real(dp), intent(in) :: x !! argument whose real kind selects the kind of the result
        real(dp) :: y

        if (condition) then
            y = 0.0_dp
        else
            y = negative_infinity(x)
        end if
    end function log_indicator

    pure elemental logical function valid_parameters(n, p) result(valid)
        real(dp), intent(in) :: n !! number of trials to check
        real(dp), intent(in) :: p !! success probability to check

        valid = ieee_is_finite(n) .and. n >= 0.0_dp .and. n == aint(n) .and. &
            p >= 0.0_dp .and. p <= 1.0_dp
    end function valid_parameters

end module scifort_binomial
