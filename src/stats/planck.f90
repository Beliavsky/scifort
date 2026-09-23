! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Planck discrete exponential distribution matching scipy.stats.planck.
! P(X=k) = (1-exp(-lambda))*exp(-lambda*k), k = 0,1,... before loc.

module scifort_planck
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, &
        positive_infinity, quiet_nan
    implicit none
    private

    public :: planck_cdf, planck_isf, planck_logcdf, planck_logpmf
    public :: planck_logsf, planck_pmf, planck_ppf, planck_sf

    interface planck_pmf
        module procedure planck_pmf_real
        module procedure planck_pmf_int
    end interface planck_pmf

    interface planck_logpmf
        module procedure planck_logpmf_real
        module procedure planck_logpmf_int
    end interface planck_logpmf

    interface planck_cdf
        module procedure planck_cdf_real
        module procedure planck_cdf_int
    end interface planck_cdf

    interface planck_sf
        module procedure planck_sf_real
        module procedure planck_sf_int
    end interface planck_sf

    interface planck_logcdf
        module procedure planck_logcdf_real
        module procedure planck_logcdf_int
    end interface planck_logcdf

    interface planck_logsf
        module procedure planck_logsf_real
        module procedure planck_logsf_int
    end interface planck_logsf

contains

    pure elemental function planck_logpmf_real(k, lambda, loc) result(y)
        real(dp), intent(in) :: k !! lattice point; non-lattice values have log-mass -Inf
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: count

        count = shifted_count(k, loc)
        if (.not. valid_lambda(lambda) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (.not. ieee_is_finite(count) .or. count < 0.0_dp .or. count /= aint(count)) then
            y = negative_infinity(k)
        else
            y = log(-expm1_safe(-lambda)) - lambda * count
        end if
    end function planck_logpmf_real

    pure elemental function planck_logpmf_int(k, lambda, loc) result(y)
        integer, intent(in) :: k !! integer evaluation point
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = planck_logpmf_real(real(k, dp), lambda, loc)
    end function planck_logpmf_int

    pure elemental function planck_pmf_real(k, lambda, loc) result(y)
        real(dp), intent(in) :: k !! lattice point
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = exp(planck_logpmf_real(k, lambda, loc))
    end function planck_pmf_real

    pure elemental function planck_pmf_int(k, lambda, loc) result(y)
        integer, intent(in) :: k !! integer evaluation point
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = planck_pmf_real(real(k, dp), lambda, loc)
    end function planck_pmf_int

    pure elemental function planck_logsf_real(k, lambda, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of log(P(X > k)); shifted floor is used
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: count

        count = floor_count(k, loc)
        if (.not. valid_lambda(lambda) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (count < 0.0_dp) then
            y = 0.0_dp
        else if (.not. ieee_is_finite(count)) then
            y = negative_infinity(k)
        else
            y = -lambda * (count + 1.0_dp)
        end if
    end function planck_logsf_real

    pure elemental function planck_logsf_int(k, lambda, loc) result(y)
        integer, intent(in) :: k !! lower integer limit
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = planck_logsf_real(real(k, dp), lambda, loc)
    end function planck_logsf_int

    pure elemental function planck_sf_real(k, lambda, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of P(X > k)
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: ly
        ly = planck_logsf_real(k, lambda, loc)
        if (ieee_is_nan(ly)) then
            y = ly
        else
            y = exp(ly)
        end if
    end function planck_sf_real

    pure elemental function planck_sf_int(k, lambda, loc) result(y)
        integer, intent(in) :: k !! lower integer limit
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = planck_sf_real(real(k, dp), lambda, loc)
    end function planck_sf_int

    pure elemental function planck_logcdf_real(k, lambda, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of log(P(X <= k)); shifted floor is used
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: count, logsf

        count = floor_count(k, loc)
        if (.not. valid_lambda(lambda) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (count < 0.0_dp) then
            y = negative_infinity(k)
        else if (.not. ieee_is_finite(count)) then
            y = 0.0_dp
        else
            logsf = -lambda * (count + 1.0_dp)
            y = log(-expm1_safe(logsf))
        end if
    end function planck_logcdf_real

    pure elemental function planck_logcdf_int(k, lambda, loc) result(y)
        integer, intent(in) :: k !! upper integer limit
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = planck_logcdf_real(real(k, dp), lambda, loc)
    end function planck_logcdf_int

    pure elemental function planck_cdf_real(k, lambda, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of P(X <= k)
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: ly
        ly = planck_logcdf_real(k, lambda, loc)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function planck_cdf_real

    pure elemental function planck_cdf_int(k, lambda, loc) result(y)
        integer, intent(in) :: k !! upper integer limit
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = planck_cdf_real(real(k, dp), lambda, loc)
    end function planck_cdf_int

    pure elemental function planck_ppf(probability, lambda, loc) result(y)
        real(dp), intent(in) :: probability !! lower-tail probability in [0,1]
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: shift, count

        shift = optional_loc(loc)
        if (.not. valid_lambda(lambda) .or. .not. ieee_is_finite(shift) .or. &
            .not. (probability >= 0.0_dp .and. probability <= 1.0_dp)) then
            y = quiet_nan(probability)
        else if (probability <= 0.0_dp) then
            y = shift - 1.0_dp
        else if (probability >= 1.0_dp) then
            y = positive_infinity(probability)
        else
            count = ceiling(-log1p_safe(-probability) / lambda - 1.0_dp)
            count = max(0.0_dp, count)
            if (count > 0.0_dp) then
                if (planck_cdf_real(shift + count - 1.0_dp, lambda, shift) >= probability) &
                    count = count - 1.0_dp
            end if
            y = shift + count
        end if
    end function planck_ppf

    pure elemental function planck_isf(probability, lambda, loc) result(y)
        real(dp), intent(in) :: probability !! upper-tail probability in [0,1]
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: shift, count

        shift = optional_loc(loc)
        if (.not. valid_lambda(lambda) .or. .not. ieee_is_finite(shift) .or. &
            .not. (probability >= 0.0_dp .and. probability <= 1.0_dp)) then
            y = quiet_nan(probability)
        else if (probability >= 1.0_dp) then
            y = shift - 1.0_dp
        else if (probability <= 0.0_dp) then
            y = positive_infinity(probability)
        else
            count = ceiling(-log(probability) / lambda - 1.0_dp)
            count = max(0.0_dp, count)
            if (count > 0.0_dp) then
                if (planck_sf_real(shift + count - 1.0_dp, lambda, shift) <= probability) &
                    count = count - 1.0_dp
            end if
            y = shift + count
        end if
    end function planck_isf

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

    pure elemental logical function valid_lambda(lambda)
        real(dp), intent(in) :: lambda !! rate parameter
        valid_lambda = ieee_is_finite(lambda) .and. lambda > 0.0_dp
    end function valid_lambda

end module scifort_planck
