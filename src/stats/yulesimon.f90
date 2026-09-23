! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Yule-Simon distribution matching scipy.stats.yulesimon.
module scifort_yulesimon
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_log_gamma, only : log_beta
    use scifort_math, only : expm1_safe, negative_infinity, positive_infinity, quiet_nan
    implicit none
    private

    public :: yulesimon_cdf, yulesimon_isf, yulesimon_logcdf, yulesimon_logpmf
    public :: yulesimon_logsf, yulesimon_pmf, yulesimon_ppf, yulesimon_sf

    interface yulesimon_pmf
        module procedure yulesimon_pmf_real
        module procedure yulesimon_pmf_int
    end interface
    interface yulesimon_logpmf
        module procedure yulesimon_logpmf_real
        module procedure yulesimon_logpmf_int
    end interface
    interface yulesimon_cdf
        module procedure yulesimon_cdf_real
        module procedure yulesimon_cdf_int
    end interface
    interface yulesimon_sf
        module procedure yulesimon_sf_real
        module procedure yulesimon_sf_int
    end interface
    interface yulesimon_logcdf
        module procedure yulesimon_logcdf_real
        module procedure yulesimon_logcdf_int
    end interface
    interface yulesimon_logsf
        module procedure yulesimon_logsf_real
        module procedure yulesimon_logsf_int
    end interface

contains

    pure elemental function yulesimon_logpmf_real(k, alpha, loc) result(y)
        real(dp), intent(in) :: k !! lattice observation
        real(dp), intent(in) :: alpha !! positive Yule-Simon shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, count
        count = shifted_count(k, loc)
        if (.not. valid_alpha(alpha) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (.not. ieee_is_finite(count) .or. count < 1.0_dp .or. count /= aint(count)) then
            y = negative_infinity(k)
        else
            y = log(alpha) + log_beta(count, alpha + 1.0_dp)
        end if
    end function

    pure elemental function yulesimon_logpmf_int(k, alpha, loc) result(y)
        integer, intent(in) :: k !! integer observation
        real(dp), intent(in) :: alpha !! positive Yule-Simon shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = yulesimon_logpmf_real(real(k, dp), alpha, loc)
    end function

    pure elemental function yulesimon_pmf_real(k, alpha, loc) result(y)
        real(dp), intent(in) :: k !! lattice observation
        real(dp), intent(in) :: alpha !! positive Yule-Simon shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = exp(yulesimon_logpmf_real(k, alpha, loc))
    end function

    pure elemental function yulesimon_pmf_int(k, alpha, loc) result(y)
        integer, intent(in) :: k !! integer observation
        real(dp), intent(in) :: alpha !! positive Yule-Simon shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = yulesimon_pmf_real(real(k, dp), alpha, loc)
    end function

    pure elemental function yulesimon_logsf_real(k, alpha, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of log(P(X > k)); shifted floor is used
        real(dp), intent(in) :: alpha !! positive Yule-Simon shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, count
        count = floor_count(k, loc)
        if (.not. valid_alpha(alpha) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (count < 1.0_dp) then
            y = 0.0_dp
        else if (.not. ieee_is_finite(count)) then
            y = negative_infinity(k)
        else
            y = log(count) + log_beta(count, alpha + 1.0_dp)
        end if
    end function

    pure elemental function yulesimon_logsf_int(k, alpha, loc) result(y)
        integer, intent(in) :: k !! lower integer limit
        real(dp), intent(in) :: alpha !! positive Yule-Simon shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = yulesimon_logsf_real(real(k, dp), alpha, loc)
    end function

    pure elemental function yulesimon_sf_real(k, alpha, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of P(X > k)
        real(dp), intent(in) :: alpha !! positive Yule-Simon shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, ly
        ly = yulesimon_logsf_real(k, alpha, loc)
        if (ieee_is_nan(ly)) then
            y = ly
        else
            y = exp(ly)
        end if
    end function

    pure elemental function yulesimon_sf_int(k, alpha, loc) result(y)
        integer, intent(in) :: k !! lower integer limit
        real(dp), intent(in) :: alpha !! positive Yule-Simon shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = yulesimon_sf_real(real(k, dp), alpha, loc)
    end function

    pure elemental function yulesimon_logcdf_real(k, alpha, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of log(P(X <= k)); shifted floor is used
        real(dp), intent(in) :: alpha !! positive Yule-Simon shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, lsf
        lsf = yulesimon_logsf_real(k, alpha, loc)
        if (ieee_is_nan(lsf)) then
            y = lsf
        else if (lsf == 0.0_dp) then
            y = negative_infinity(k)
        else if (lsf == negative_infinity(lsf)) then
            y = 0.0_dp
        else
            y = log(-expm1_safe(lsf))
        end if
    end function

    pure elemental function yulesimon_logcdf_int(k, alpha, loc) result(y)
        integer, intent(in) :: k !! upper integer limit
        real(dp), intent(in) :: alpha !! positive Yule-Simon shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = yulesimon_logcdf_real(real(k, dp), alpha, loc)
    end function

    pure elemental function yulesimon_cdf_real(k, alpha, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of P(X <= k)
        real(dp), intent(in) :: alpha !! positive Yule-Simon shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, lsf
        lsf = yulesimon_logsf_real(k, alpha, loc)
        if (ieee_is_nan(lsf)) then
            y = lsf
        else
            y = -expm1_safe(lsf)
        end if
    end function

    pure elemental function yulesimon_cdf_int(k, alpha, loc) result(y)
        integer, intent(in) :: k !! upper integer limit
        real(dp), intent(in) :: alpha !! positive Yule-Simon shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = yulesimon_cdf_real(real(k, dp), alpha, loc)
    end function

    pure elemental function yulesimon_ppf(probability, alpha, loc) result(y)
        real(dp), intent(in) :: probability !! lower-tail probability in [0,1]
        real(dp), intent(in) :: alpha !! positive Yule-Simon shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, shift
        integer :: lo, hi, mid
        shift = optional_loc(loc)
        if (.not. valid_alpha(alpha) .or. .not. ieee_is_finite(shift) .or. .not. valid_probability(probability)) then
            y = quiet_nan(probability)
        else if (probability <= 0.0_dp) then
            y = shift
        else if (probability >= 1.0_dp) then
            y = positive_infinity(probability)
        else
            lo = 1; hi = 1
            do while (yulesimon_cdf_real(shift + real(hi, dp), alpha, shift) < probability .and. hi < huge(hi)/2)
                hi = hi * 2
            end do
            do while (lo < hi)
                mid = lo + (hi - lo)/2
                if (yulesimon_cdf_real(shift + real(mid, dp), alpha, shift) >= probability) then
                    hi = mid
                else
                    lo = mid + 1
                end if
            end do
            y = shift + real(lo, dp)
        end if
    end function

    pure elemental function yulesimon_isf(probability, alpha, loc) result(y)
        real(dp), intent(in) :: probability !! upper-tail probability in [0,1]
        real(dp), intent(in) :: alpha !! positive Yule-Simon shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, shift
        integer :: lo, hi, mid
        shift = optional_loc(loc)
        if (.not. valid_alpha(alpha) .or. .not. ieee_is_finite(shift) .or. .not. valid_probability(probability)) then
            y = quiet_nan(probability)
        else if (probability >= 1.0_dp) then
            y = shift
        else if (probability <= 0.0_dp) then
            y = positive_infinity(probability)
        else
            lo = 1; hi = 1
            do while (yulesimon_sf_real(shift + real(hi, dp), alpha, shift) > probability .and. hi < huge(hi)/2)
                hi = hi * 2
            end do
            do while (lo < hi)
                mid = lo + (hi - lo)/2
                if (yulesimon_sf_real(shift + real(mid, dp), alpha, shift) <= probability) then
                    hi = mid
                else
                    lo = mid + 1
                end if
            end do
            y = shift + real(lo, dp)
        end if
    end function

    pure elemental function shifted_count(k, loc) result(count)
        real(dp), intent(in) :: k !! observation
        real(dp), intent(in), optional :: loc !! support shift
        real(dp) :: count
        count = k - optional_loc(loc)
    end function

    pure elemental function floor_count(k, loc) result(count)
        real(dp), intent(in) :: k !! observation
        real(dp), intent(in), optional :: loc !! support shift
        real(dp) :: count
        count = floor(k - optional_loc(loc))
    end function

    pure elemental function optional_loc(loc) result(value)
        real(dp), intent(in), optional :: loc !! support shift
        real(dp) :: value
        value = 0.0_dp
        if (present(loc)) value = loc
    end function

    pure elemental logical function valid_alpha(alpha)
        real(dp), intent(in) :: alpha !! shape to validate
        valid_alpha = ieee_is_finite(alpha) .and. alpha > 0.0_dp
    end function

    pure elemental logical function valid_probability(p)
        real(dp), intent(in) :: p !! probability to validate
        valid_probability = p >= 0.0_dp .and. p <= 1.0_dp
    end function
end module scifort_yulesimon
