! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Boltzmann (truncated discrete exponential) distribution matching
! scipy.stats.boltzmann. Stable log-tail formulas use expm1 near lambda=0.

module scifort_boltzmann
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, quiet_nan
    implicit none
    private

    public :: boltzmann_cdf, boltzmann_isf, boltzmann_logcdf, boltzmann_logpmf
    public :: boltzmann_logsf, boltzmann_pmf, boltzmann_ppf, boltzmann_sf

    interface boltzmann_pmf
        module procedure boltzmann_pmf_real
        module procedure boltzmann_pmf_int
    end interface boltzmann_pmf
    interface boltzmann_logpmf
        module procedure boltzmann_logpmf_real
        module procedure boltzmann_logpmf_int
    end interface boltzmann_logpmf
    interface boltzmann_cdf
        module procedure boltzmann_cdf_real
        module procedure boltzmann_cdf_int
    end interface boltzmann_cdf
    interface boltzmann_sf
        module procedure boltzmann_sf_real
        module procedure boltzmann_sf_int
    end interface boltzmann_sf
    interface boltzmann_logcdf
        module procedure boltzmann_logcdf_real
        module procedure boltzmann_logcdf_int
    end interface boltzmann_logcdf
    interface boltzmann_logsf
        module procedure boltzmann_logsf_real
        module procedure boltzmann_logsf_int
    end interface boltzmann_logsf

contains

    pure elemental function boltzmann_logpmf_real(k, lambda, n, loc) result(y)
        real(dp), intent(in) :: k !! lattice observation
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in) :: n !! positive integer support size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: count

        count = shifted_count(k, loc)
        if (.not. valid_parameters(lambda, n) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (.not. ieee_is_finite(count) .or. count /= aint(count) .or. count < 0.0_dp .or. count >= n) then
            y = negative_infinity(k)
        else
            y = log1mexp(-lambda) - lambda * count - log1mexp(-lambda * n)
        end if
    end function boltzmann_logpmf_real

    pure elemental function boltzmann_logpmf_int(k, lambda, n, loc) result(y)
        integer, intent(in) :: k !! integer observation
        real(dp), intent(in) :: lambda !! positive exponential rate
        integer, intent(in) :: n !! positive support size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = boltzmann_logpmf_real(real(k, dp), lambda, real(n, dp), loc)
    end function boltzmann_logpmf_int

    pure elemental function boltzmann_pmf_real(k, lambda, n, loc) result(y)
        real(dp), intent(in) :: k !! lattice observation
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in) :: n !! positive integer support size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = exp(boltzmann_logpmf_real(k, lambda, n, loc))
    end function boltzmann_pmf_real

    pure elemental function boltzmann_pmf_int(k, lambda, n, loc) result(y)
        integer, intent(in) :: k !! integer observation
        real(dp), intent(in) :: lambda !! positive exponential rate
        integer, intent(in) :: n !! positive support size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = boltzmann_pmf_real(real(k, dp), lambda, real(n, dp), loc)
    end function boltzmann_pmf_int

    pure elemental function boltzmann_cdf_real(k, lambda, n, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of P(X <= k)
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in) :: n !! positive integer support size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: logvalue
        logvalue = boltzmann_logcdf_real(k, lambda, n, loc)
        y = exp(logvalue)
    end function boltzmann_cdf_real

    pure elemental function boltzmann_cdf_int(k, lambda, n, loc) result(y)
        integer, intent(in) :: k !! upper integer limit
        real(dp), intent(in) :: lambda !! positive exponential rate
        integer, intent(in) :: n !! positive support size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = boltzmann_cdf_real(real(k, dp), lambda, real(n, dp), loc)
    end function boltzmann_cdf_int

    pure elemental function boltzmann_sf_real(k, lambda, n, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of P(X > k)
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in) :: n !! positive integer support size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: logvalue
        logvalue = boltzmann_logsf_real(k, lambda, n, loc)
        y = exp(logvalue)
    end function boltzmann_sf_real

    pure elemental function boltzmann_sf_int(k, lambda, n, loc) result(y)
        integer, intent(in) :: k !! lower integer limit
        real(dp), intent(in) :: lambda !! positive exponential rate
        integer, intent(in) :: n !! positive support size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = boltzmann_sf_real(real(k, dp), lambda, real(n, dp), loc)
    end function boltzmann_sf_int

    pure elemental function boltzmann_logcdf_real(k, lambda, n, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of log(P(X <= k))
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in) :: n !! positive integer support size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: count

        count = floor_count(k, loc)
        if (.not. valid_parameters(lambda, n) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (count < 0.0_dp) then
            y = negative_infinity(k)
        else if (.not. ieee_is_finite(count) .or. count >= n - 1.0_dp) then
            y = 0.0_dp
        else
            y = log1mexp(-lambda * (count + 1.0_dp)) - log1mexp(-lambda * n)
        end if
    end function boltzmann_logcdf_real

    pure elemental function boltzmann_logcdf_int(k, lambda, n, loc) result(y)
        integer, intent(in) :: k !! upper integer limit
        real(dp), intent(in) :: lambda !! positive exponential rate
        integer, intent(in) :: n !! positive support size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = boltzmann_logcdf_real(real(k, dp), lambda, real(n, dp), loc)
    end function boltzmann_logcdf_int

    pure elemental function boltzmann_logsf_real(k, lambda, n, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of log(P(X > k))
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in) :: n !! positive integer support size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: count

        count = floor_count(k, loc)
        if (.not. valid_parameters(lambda, n) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (count < 0.0_dp) then
            y = 0.0_dp
        else if (.not. ieee_is_finite(count) .or. count >= n - 1.0_dp) then
            y = negative_infinity(k)
        else
            y = -lambda * (count + 1.0_dp) + log1mexp(-lambda * (n - count - 1.0_dp)) - &
                log1mexp(-lambda * n)
        end if
    end function boltzmann_logsf_real

    pure elemental function boltzmann_logsf_int(k, lambda, n, loc) result(y)
        integer, intent(in) :: k !! lower integer limit
        real(dp), intent(in) :: lambda !! positive exponential rate
        integer, intent(in) :: n !! positive support size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = boltzmann_logsf_real(real(k, dp), lambda, real(n, dp), loc)
    end function boltzmann_logsf_int

    pure elemental function boltzmann_ppf(probability, lambda, n, loc) result(y)
        real(dp), intent(in) :: probability !! lower-tail probability in [0,1]
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in) :: n !! positive integer support size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: shift
        integer :: lo, hi, mid

        shift = optional_loc(loc)
        if (.not. valid_parameters(lambda, n) .or. .not. ieee_is_finite(shift) .or. &
            .not. valid_probability(probability)) then
            y = quiet_nan(probability)
        else if (probability <= 0.0_dp) then
            y = shift - 1.0_dp
        else if (probability >= 1.0_dp) then
            y = shift + n - 1.0_dp
        else
            lo = 0; hi = int(n) - 1
            do while (lo < hi)
                mid = lo + (hi - lo) / 2
                if (boltzmann_cdf_real(shift + real(mid, dp), lambda, n, shift) >= probability) then
                    hi = mid
                else
                    lo = mid + 1
                end if
            end do
            y = shift + real(lo, dp)
        end if
    end function boltzmann_ppf

    pure elemental function boltzmann_isf(probability, lambda, n, loc) result(y)
        real(dp), intent(in) :: probability !! upper-tail probability in [0,1]
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in) :: n !! positive integer support size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: shift
        integer :: lo, hi, mid

        shift = optional_loc(loc)
        if (.not. valid_parameters(lambda, n) .or. .not. ieee_is_finite(shift) .or. &
            .not. valid_probability(probability)) then
            y = quiet_nan(probability)
        else if (probability >= 1.0_dp) then
            y = shift - 1.0_dp
        else if (probability <= 0.0_dp) then
            y = shift + n - 1.0_dp
        else
            lo = 0; hi = int(n) - 1
            do while (lo < hi)
                mid = lo + (hi - lo) / 2
                if (boltzmann_sf_real(shift + real(mid, dp), lambda, n, shift) <= probability) then
                    hi = mid
                else
                    lo = mid + 1
                end if
            end do
            y = shift + real(lo, dp)
        end if
    end function boltzmann_isf

    pure elemental function log1mexp(x) result(y)
        real(dp), intent(in) :: x !! nonpositive logarithm
        real(dp) :: y
        if (x > -0.693147180559945309417232121458_dp) then
            y = log(-expm1_safe(x))
        else
            y = log1p_safe(-exp(x))
        end if
    end function log1mexp

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

    pure elemental logical function valid_parameters(lambda, n)
        real(dp), intent(in) :: lambda !! exponential rate
        real(dp), intent(in) :: n !! support size
        valid_parameters = ieee_is_finite(lambda) .and. lambda > 0.0_dp .and. ieee_is_finite(n) .and. &
            n > 0.0_dp .and. n == aint(n) .and. n <= real(huge(0), dp)
    end function valid_parameters

    pure elemental logical function valid_probability(p)
        real(dp), intent(in) :: p !! probability candidate
        valid_probability = ieee_is_finite(p) .and. p >= 0.0_dp .and. p <= 1.0_dp
    end function valid_probability

end module scifort_boltzmann
