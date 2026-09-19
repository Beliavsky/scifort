! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Poisson distribution with mean mu >= 0 and integer location loc, matching
! scipy.stats.poisson. The support is loc, loc + 1, ... and a non-integer
! argument has probability zero.
!
! DLMF 8.4.8 with 8.4.11 gives Gamma(n + 1, z) = n! exp(-z) sum over
! j = 0..n of z**j / j!, so the regularized upper function Q(n + 1, z) is the
! Poisson cumulative probability. The tails are therefore evaluated as
! P(X <= k) = Q(k + 1, mu) and P(X > k) = P(k + 1, mu) without forming any
! sum, and the mass function reuses the stable kernel
! mu**k exp(-mu) / Gamma(k + 1) of scifort_incomplete_gamma.

module scifort_poisson
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_incomplete_gamma, only : gammainc, gammaincc, log_gamma_kernel, &
        log_gammainc, log_gammaincc
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan
    use scifort_normal, only : normal_isf, normal_ppf
    implicit none
    private

    public :: poisson_cdf
    public :: poisson_isf
    public :: poisson_logcdf
    public :: poisson_logpmf
    public :: poisson_logsf
    public :: poisson_pmf
    public :: poisson_ppf
    public :: poisson_sf

    !> Probability mass function.
    interface poisson_pmf
        module procedure poisson_pmf_real
        module procedure poisson_pmf_int
    end interface poisson_pmf

    !> Logarithm of the probability mass function.
    interface poisson_logpmf
        module procedure poisson_logpmf_real
        module procedure poisson_logpmf_int
    end interface poisson_logpmf

    !> Cumulative distribution function P(X <= k).
    interface poisson_cdf
        module procedure poisson_cdf_real
        module procedure poisson_cdf_int
    end interface poisson_cdf

    !> Survival function P(X > k).
    interface poisson_sf
        module procedure poisson_sf_real
        module procedure poisson_sf_int
    end interface poisson_sf

    !> Logarithm of the cumulative distribution function.
    interface poisson_logcdf
        module procedure poisson_logcdf_real
        module procedure poisson_logcdf_int
    end interface poisson_logcdf

    !> Logarithm of the survival function.
    interface poisson_logsf
        module procedure poisson_logsf_real
        module procedure poisson_logsf_int
    end interface poisson_logsf

    integer, parameter :: max_search_steps = 200

contains

    pure elemental function poisson_pmf_real(k, mu, loc) result(y)
        real(dp), intent(in) :: k !! number of events; non-integer values have probability zero
        real(dp), intent(in) :: mu !! mean, >= 0
        real(dp), intent(in), optional :: loc !! integer shift of the support (default 0)
        real(dp) :: y

        real(dp) :: count

        count = shifted_count(k, loc)
        if (.not. valid_mean(mu) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (count < 0.0_dp .or. .not. ieee_is_finite(count)) then
            y = 0.0_dp
        else if (count /= aint(count)) then
            y = 0.0_dp
        else if (mu == 0.0_dp) then
            y = merge(1.0_dp, 0.0_dp, count == 0.0_dp)
        else
            y = exp(log_gamma_kernel(count, mu))
        end if
    end function poisson_pmf_real

    pure elemental function poisson_pmf_int(k, mu, loc) result(y)
        integer, intent(in) :: k !! number of events
        real(dp), intent(in) :: mu !! mean, >= 0
        real(dp), intent(in), optional :: loc !! integer shift of the support (default 0)
        real(dp) :: y

        y = poisson_pmf_real(real(k, dp), mu, loc)
    end function poisson_pmf_int

    pure elemental function poisson_logpmf_real(k, mu, loc) result(y)
        real(dp), intent(in) :: k !! number of events; non-integer values have probability zero
        real(dp), intent(in) :: mu !! mean, >= 0
        real(dp), intent(in), optional :: loc !! integer shift of the support (default 0)
        real(dp) :: y

        real(dp) :: count

        count = shifted_count(k, loc)
        if (.not. valid_mean(mu) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (count < 0.0_dp .or. .not. ieee_is_finite(count)) then
            y = negative_infinity(k)
        else if (count /= aint(count)) then
            y = negative_infinity(k)
        else if (mu == 0.0_dp) then
            if (count == 0.0_dp) then
                y = 0.0_dp
            else
                y = negative_infinity(k)
            end if
        else
            y = log_gamma_kernel(count, mu)
        end if
    end function poisson_logpmf_real

    pure elemental function poisson_logpmf_int(k, mu, loc) result(y)
        integer, intent(in) :: k !! number of events
        real(dp), intent(in) :: mu !! mean, >= 0
        real(dp), intent(in), optional :: loc !! integer shift of the support (default 0)
        real(dp) :: y

        y = poisson_logpmf_real(real(k, dp), mu, loc)
    end function poisson_logpmf_int

    pure elemental function poisson_cdf_real(k, mu, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of P(X <= k); the floor is used
        real(dp), intent(in) :: mu !! mean, >= 0
        real(dp), intent(in), optional :: loc !! integer shift of the support (default 0)
        real(dp) :: y

        real(dp) :: count

        count = floor_count(k, loc)
        if (.not. valid_mean(mu) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (count < 0.0_dp) then
            y = 0.0_dp
        else if (mu == 0.0_dp) then
            y = 1.0_dp
        else
            y = gammaincc(count + 1.0_dp, mu)
        end if
    end function poisson_cdf_real

    pure elemental function poisson_cdf_int(k, mu, loc) result(y)
        integer, intent(in) :: k !! upper limit of P(X <= k)
        real(dp), intent(in) :: mu !! mean, >= 0
        real(dp), intent(in), optional :: loc !! integer shift of the support (default 0)
        real(dp) :: y

        y = poisson_cdf_real(real(k, dp), mu, loc)
    end function poisson_cdf_int

    pure elemental function poisson_sf_real(k, mu, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of P(X > k); the floor is used
        real(dp), intent(in) :: mu !! mean, >= 0
        real(dp), intent(in), optional :: loc !! integer shift of the support (default 0)
        real(dp) :: y

        real(dp) :: count

        count = floor_count(k, loc)
        if (.not. valid_mean(mu) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (count < 0.0_dp) then
            y = 1.0_dp
        else if (mu == 0.0_dp) then
            y = 0.0_dp
        else
            y = gammainc(count + 1.0_dp, mu)
        end if
    end function poisson_sf_real

    pure elemental function poisson_sf_int(k, mu, loc) result(y)
        integer, intent(in) :: k !! lower limit of P(X > k)
        real(dp), intent(in) :: mu !! mean, >= 0
        real(dp), intent(in), optional :: loc !! integer shift of the support (default 0)
        real(dp) :: y

        y = poisson_sf_real(real(k, dp), mu, loc)
    end function poisson_sf_int

    pure elemental function poisson_logcdf_real(k, mu, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of log(P(X <= k)); the floor is used
        real(dp), intent(in) :: mu !! mean, >= 0
        real(dp), intent(in), optional :: loc !! integer shift of the support (default 0)
        real(dp) :: y

        real(dp) :: count

        count = floor_count(k, loc)
        if (.not. valid_mean(mu) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (count < 0.0_dp) then
            y = negative_infinity(k)
        else if (mu == 0.0_dp) then
            y = 0.0_dp
        else
            y = log_gammaincc(count + 1.0_dp, mu)
        end if
    end function poisson_logcdf_real

    pure elemental function poisson_logcdf_int(k, mu, loc) result(y)
        integer, intent(in) :: k !! upper limit of log(P(X <= k))
        real(dp), intent(in) :: mu !! mean, >= 0
        real(dp), intent(in), optional :: loc !! integer shift of the support (default 0)
        real(dp) :: y

        y = poisson_logcdf_real(real(k, dp), mu, loc)
    end function poisson_logcdf_int

    pure elemental function poisson_logsf_real(k, mu, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of log(P(X > k)); the floor is used
        real(dp), intent(in) :: mu !! mean, >= 0
        real(dp), intent(in), optional :: loc !! integer shift of the support (default 0)
        real(dp) :: y

        real(dp) :: count

        count = floor_count(k, loc)
        if (.not. valid_mean(mu) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (count < 0.0_dp) then
            y = 0.0_dp
        else if (mu == 0.0_dp) then
            y = negative_infinity(k)
        else
            y = log_gammainc(count + 1.0_dp, mu)
        end if
    end function poisson_logsf_real

    pure elemental function poisson_logsf_int(k, mu, loc) result(y)
        integer, intent(in) :: k !! lower limit of log(P(X > k))
        real(dp), intent(in) :: mu !! mean, >= 0
        real(dp), intent(in), optional :: loc !! integer shift of the support (default 0)
        real(dp) :: y

        y = poisson_logsf_real(real(k, dp), mu, loc)
    end function poisson_logsf_int

    ! Smallest k in the support with P(X <= k) >= p. As in scipy.stats, a
    ! probability of zero gives loc - 1, the value just below the support.
    pure elemental function poisson_ppf(p, mu, loc) result(y)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: mu !! mean, >= 0
        real(dp), intent(in), optional :: loc !! integer shift of the support (default 0)
        real(dp) :: y

        real(dp) :: shift

        shift = 0.0_dp
        if (present(loc)) shift = loc
        if (.not. valid_mean(mu) .or. .not. ieee_is_finite(shift)) then
            y = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            y = quiet_nan(p)
        else if (p == 0.0_dp) then
            y = shift - 1.0_dp
        else if (p >= 1.0_dp) then
            y = positive_infinity(p)
        else if (mu == 0.0_dp) then
            y = shift
        else
            y = shift + search_quantile(p, mu, .true.)
        end if
    end function poisson_ppf

    ! Smallest k in the support with P(X > k) <= p.
    pure elemental function poisson_isf(p, mu, loc) result(y)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: mu !! mean, >= 0
        real(dp), intent(in), optional :: loc !! integer shift of the support (default 0)
        real(dp) :: y

        real(dp) :: shift

        shift = 0.0_dp
        if (present(loc)) shift = loc
        if (.not. valid_mean(mu) .or. .not. ieee_is_finite(shift)) then
            y = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            y = quiet_nan(p)
        else if (p >= 1.0_dp) then
            y = shift - 1.0_dp
        else if (p == 0.0_dp) then
            y = positive_infinity(p)
        else if (mu == 0.0_dp) then
            y = shift
        else
            y = shift + search_quantile(p, mu, .false.)
        end if
    end function poisson_isf

    ! Smallest count with P(X <= k) >= p (lower) or P(X > k) <= p (upper).
    ! A normal approximation with a continuity correction starts the search,
    ! a bracket is found by doubling, and the bracket is then halved.
    pure function search_quantile(p, mu, lower) result(count)
        real(dp), intent(in) :: p !! target probability, 0 < p < 1
        real(dp), intent(in) :: mu !! mean, > 0 and finite
        logical, intent(in) :: lower !! .true. for the lower tail, .false. for the upper
        real(dp) :: count

        integer :: step
        real(dp) :: hi
        real(dp) :: lo
        real(dp) :: middle
        real(dp) :: width
        real(dp) :: z

        if (.not. ieee_is_finite(mu)) then
            count = positive_infinity(p)
            return
        end if

        if (lower) then
            z = normal_ppf(p)
        else
            z = normal_isf(p)
        end if
        count = max(0.0_dp, aint(mu + z * sqrt(mu) + 0.5_dp))

        if (satisfies(count, p, mu, lower)) then
            ! Walk down to the first count that still satisfies the condition.
            hi = count
            width = 1.0_dp
            lo = max(-1.0_dp, count - width)
            do step = 1, max_search_steps
                if (lo < 0.0_dp) exit
                if (.not. satisfies(lo, p, mu, lower)) exit
                hi = lo
                width = 2.0_dp * width
                lo = max(-1.0_dp, lo - width)
            end do
        else
            lo = count
            width = 1.0_dp
            hi = count + width
            do step = 1, max_search_steps
                if (satisfies(hi, p, mu, lower)) exit
                lo = hi
                width = 2.0_dp * width
                hi = hi + width
            end do
        end if

        ! Invariant: lo does not satisfy the condition, hi does.
        do step = 1, max_search_steps
            if (hi - lo <= 1.0_dp) exit
            middle = aint(lo + 0.5_dp * (hi - lo))
            if (satisfies(middle, p, mu, lower)) then
                hi = middle
            else
                lo = middle
            end if
        end do
        count = hi
    end function search_quantile

    pure function satisfies(count, p, mu, lower) result(ok)
        real(dp), intent(in) :: count !! candidate count, >= 0
        real(dp), intent(in) :: p !! target probability
        real(dp), intent(in) :: mu !! mean, > 0 and finite
        logical, intent(in) :: lower !! .true. for the lower tail, .false. for the upper
        logical :: ok

        if (count < 0.0_dp) then
            ok = .false.
        else if (lower) then
            ok = gammaincc(count + 1.0_dp, mu) >= p
        else
            ok = gammainc(count + 1.0_dp, mu) <= p
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

    pure elemental logical function valid_mean(mu) result(valid)
        real(dp), intent(in) :: mu !! mean to check

        valid = mu >= 0.0_dp
    end function valid_mean

end module scifort_poisson
