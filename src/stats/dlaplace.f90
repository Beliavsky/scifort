! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Discrete Laplace distribution matching scipy.stats.dlaplace.
! P(X=k) = tanh(a/2)*exp(-a*abs(k)), k integer before loc.

module scifort_dlaplace
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, &
        positive_infinity, quiet_nan
    implicit none
    private

    public :: dlaplace_cdf, dlaplace_isf, dlaplace_logcdf, dlaplace_logpmf
    public :: dlaplace_logsf, dlaplace_pmf, dlaplace_ppf, dlaplace_sf

    interface dlaplace_pmf
        module procedure dlaplace_pmf_real
        module procedure dlaplace_pmf_int
    end interface dlaplace_pmf

    interface dlaplace_logpmf
        module procedure dlaplace_logpmf_real
        module procedure dlaplace_logpmf_int
    end interface dlaplace_logpmf

    interface dlaplace_cdf
        module procedure dlaplace_cdf_real
        module procedure dlaplace_cdf_int
    end interface dlaplace_cdf

    interface dlaplace_sf
        module procedure dlaplace_sf_real
        module procedure dlaplace_sf_int
    end interface dlaplace_sf

    interface dlaplace_logcdf
        module procedure dlaplace_logcdf_real
        module procedure dlaplace_logcdf_int
    end interface dlaplace_logcdf

    interface dlaplace_logsf
        module procedure dlaplace_logsf_real
        module procedure dlaplace_logsf_int
    end interface dlaplace_logsf

contains

    pure elemental function dlaplace_logpmf_real(k, a, loc) result(y)
        real(dp), intent(in) :: k !! lattice point; non-lattice values have log-mass -Inf
        real(dp), intent(in) :: a !! positive concentration/rate parameter
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: count

        count = shifted_count(k, loc)
        if (.not. valid_a(a) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (.not. ieee_is_finite(count) .or. count /= aint(count)) then
            y = negative_infinity(k)
        else
            y = log_tanh_half(a) - a * abs(count)
        end if
    end function dlaplace_logpmf_real

    pure elemental function dlaplace_logpmf_int(k, a, loc) result(y)
        integer, intent(in) :: k !! integer evaluation point
        real(dp), intent(in) :: a !! positive concentration/rate parameter
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = dlaplace_logpmf_real(real(k, dp), a, loc)
    end function dlaplace_logpmf_int

    pure elemental function dlaplace_pmf_real(k, a, loc) result(y)
        real(dp), intent(in) :: k !! lattice point
        real(dp), intent(in) :: a !! positive concentration/rate parameter
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = exp(dlaplace_logpmf_real(k, a, loc))
    end function dlaplace_pmf_real

    pure elemental function dlaplace_pmf_int(k, a, loc) result(y)
        integer, intent(in) :: k !! integer evaluation point
        real(dp), intent(in) :: a !! positive concentration/rate parameter
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = dlaplace_pmf_real(real(k, dp), a, loc)
    end function dlaplace_pmf_int

    pure elemental function dlaplace_logcdf_real(k, a, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of log(P(X <= k)); shifted floor is used
        real(dp), intent(in) :: a !! positive concentration/rate parameter
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: count, logden, logtail

        count = floor_count(k, loc)
        if (.not. valid_a(a) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (count == negative_infinity(count)) then
            y = negative_infinity(k)
        else if (count == positive_infinity(count)) then
            y = 0.0_dp
        else
            logden = log1p_safe(exp(-a))
            if (count < 0.0_dp) then
                y = a * count - logden
            else
                logtail = -a * (count + 1.0_dp) - logden
                y = log(-expm1_safe(logtail))
            end if
        end if
    end function dlaplace_logcdf_real

    pure elemental function dlaplace_logcdf_int(k, a, loc) result(y)
        integer, intent(in) :: k !! upper integer limit
        real(dp), intent(in) :: a !! positive concentration/rate parameter
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = dlaplace_logcdf_real(real(k, dp), a, loc)
    end function dlaplace_logcdf_int

    pure elemental function dlaplace_cdf_real(k, a, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of P(X <= k)
        real(dp), intent(in) :: a !! positive concentration/rate parameter
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: ly
        ly = dlaplace_logcdf_real(k, a, loc)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function dlaplace_cdf_real

    pure elemental function dlaplace_cdf_int(k, a, loc) result(y)
        integer, intent(in) :: k !! upper integer limit
        real(dp), intent(in) :: a !! positive concentration/rate parameter
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = dlaplace_cdf_real(real(k, dp), a, loc)
    end function dlaplace_cdf_int

    pure elemental function dlaplace_logsf_real(k, a, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of log(P(X > k)); shifted floor is used
        real(dp), intent(in) :: a !! positive concentration/rate parameter
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: count, logden, loglower

        count = floor_count(k, loc)
        if (.not. valid_a(a) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (count == negative_infinity(count)) then
            y = 0.0_dp
        else if (count == positive_infinity(count)) then
            y = negative_infinity(k)
        else
            logden = log1p_safe(exp(-a))
            if (count >= 0.0_dp) then
                y = -a * (count + 1.0_dp) - logden
            else
                loglower = a * count - logden
                y = log(-expm1_safe(loglower))
            end if
        end if
    end function dlaplace_logsf_real

    pure elemental function dlaplace_logsf_int(k, a, loc) result(y)
        integer, intent(in) :: k !! lower integer limit
        real(dp), intent(in) :: a !! positive concentration/rate parameter
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = dlaplace_logsf_real(real(k, dp), a, loc)
    end function dlaplace_logsf_int

    pure elemental function dlaplace_sf_real(k, a, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of P(X > k)
        real(dp), intent(in) :: a !! positive concentration/rate parameter
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: ly
        ly = dlaplace_logsf_real(k, a, loc)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function dlaplace_sf_real

    pure elemental function dlaplace_sf_int(k, a, loc) result(y)
        integer, intent(in) :: k !! lower integer limit
        real(dp), intent(in) :: a !! positive concentration/rate parameter
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = dlaplace_sf_real(real(k, dp), a, loc)
    end function dlaplace_sf_int

    pure elemental function dlaplace_ppf(probability, a, loc) result(y)
        real(dp), intent(in) :: probability !! lower-tail probability in [0,1]
        real(dp), intent(in) :: a !! positive concentration/rate parameter
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: shift, logden, threshold, count

        shift = optional_loc(loc)
        if (.not. valid_a(a) .or. .not. ieee_is_finite(shift) .or. &
            .not. (probability >= 0.0_dp .and. probability <= 1.0_dp)) then
            y = quiet_nan(probability)
        else if (probability <= 0.0_dp) then
            y = negative_infinity(probability)
        else if (probability >= 1.0_dp) then
            y = positive_infinity(probability)
        else
            logden = a + log1p_safe(exp(-a))
            threshold = 1.0_dp / (1.0_dp + exp(-a))
            if (probability < threshold) then
                count = ceiling((log(probability) + logden) / a - 1.0_dp)
            else
                count = ceiling(-(log1p_safe(-probability) + logden) / a)
            end if
            if (dlaplace_cdf_real(shift + count - 1.0_dp, a, shift) >= probability) &
                count = count - 1.0_dp
            y = shift + count
        end if
    end function dlaplace_ppf

    pure elemental function dlaplace_isf(probability, a, loc) result(y)
        real(dp), intent(in) :: probability !! upper-tail probability in [0,1]
        real(dp), intent(in) :: a !! positive concentration/rate parameter
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: shift, lower_q

        shift = optional_loc(loc)
        if (.not. valid_a(a) .or. .not. ieee_is_finite(shift) .or. &
            .not. (probability >= 0.0_dp .and. probability <= 1.0_dp)) then
            y = quiet_nan(probability)
        else if (probability <= 0.0_dp) then
            y = positive_infinity(probability)
        else if (probability >= 1.0_dp) then
            y = negative_infinity(probability)
        else
            lower_q = dlaplace_ppf(probability, a, 0.0_dp)
            y = shift - lower_q
        end if
    end function dlaplace_isf

    pure elemental function log_tanh_half(a) result(y)
        real(dp), intent(in) :: a !! positive rate
        real(dp) :: y
        y = log(-expm1_safe(-a)) - log1p_safe(exp(-a))
    end function log_tanh_half

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

    pure elemental logical function valid_a(a)
        real(dp), intent(in) :: a !! rate parameter
        valid_a = ieee_is_finite(a) .and. a > 0.0_dp
    end function valid_a

end module scifort_dlaplace
