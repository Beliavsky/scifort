! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Geometric distribution matching scipy.stats.geom: X is the number of trials
! through the first success, with support 1, 2, ... before applying loc.
! The implementation follows the defining pmf p (1-p)**(k-1). CDF and SF use
! expm1/log1p forms so both tails remain accurate when p or the tail is small.

module scifort_geometric
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, &
        positive_infinity, quiet_nan
    implicit none
    private

    public :: geometric_cdf
    public :: geometric_isf
    public :: geometric_logcdf
    public :: geometric_logpmf
    public :: geometric_logsf
    public :: geometric_pmf
    public :: geometric_ppf
    public :: geometric_sf

    interface geometric_pmf
        module procedure geometric_pmf_real
        module procedure geometric_pmf_int
    end interface geometric_pmf

    interface geometric_logpmf
        module procedure geometric_logpmf_real
        module procedure geometric_logpmf_int
    end interface geometric_logpmf

    interface geometric_cdf
        module procedure geometric_cdf_real
        module procedure geometric_cdf_int
    end interface geometric_cdf

    interface geometric_sf
        module procedure geometric_sf_real
        module procedure geometric_sf_int
    end interface geometric_sf

    interface geometric_logcdf
        module procedure geometric_logcdf_real
        module procedure geometric_logcdf_int
    end interface geometric_logcdf

    interface geometric_logsf
        module procedure geometric_logsf_real
        module procedure geometric_logsf_int
    end interface geometric_logsf

contains

    pure elemental function geometric_pmf_real(k, p, loc) result(y)
        real(dp), intent(in) :: k !! trial count; non-integer values have probability zero
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        y = exp(geometric_logpmf_real(k, p, loc))
    end function geometric_pmf_real

    pure elemental function geometric_pmf_int(k, p, loc) result(y)
        integer, intent(in) :: k !! trial count through the first success
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        y = geometric_pmf_real(real(k, dp), p, loc)
    end function geometric_pmf_int

    pure elemental function geometric_logpmf_real(k, p, loc) result(y)
        real(dp), intent(in) :: k !! trial count; non-integer values have probability zero
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        real(dp) :: count

        count = shifted_count(k, loc)
        if (.not. valid_probability(p) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (count < 1.0_dp .or. .not. ieee_is_finite(count)) then
            y = negative_infinity(k)
        else if (count /= aint(count)) then
            y = negative_infinity(k)
        else if (p == 1.0_dp) then
            if (count == 1.0_dp) then
                y = 0.0_dp
            else
                y = negative_infinity(k)
            end if
        else
            y = log(p) + (count - 1.0_dp) * log1p_safe(-p)
        end if
    end function geometric_logpmf_real

    pure elemental function geometric_logpmf_int(k, p, loc) result(y)
        integer, intent(in) :: k !! trial count through the first success
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        y = geometric_logpmf_real(real(k, dp), p, loc)
    end function geometric_logpmf_int

    pure elemental function geometric_cdf_real(k, p, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of P(X <= k); the floor is used
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        real(dp) :: logcdf
        real(dp) :: logsf
        real(dp) :: sf

        call tails(floor_count(k, loc), p, y, sf, logcdf, logsf)
    end function geometric_cdf_real

    pure elemental function geometric_cdf_int(k, p, loc) result(y)
        integer, intent(in) :: k !! upper limit of P(X <= k)
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        y = geometric_cdf_real(real(k, dp), p, loc)
    end function geometric_cdf_int

    pure elemental function geometric_sf_real(k, p, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of P(X > k); the floor is used
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        real(dp) :: cdf
        real(dp) :: logcdf
        real(dp) :: logsf

        call tails(floor_count(k, loc), p, cdf, y, logcdf, logsf)
    end function geometric_sf_real

    pure elemental function geometric_sf_int(k, p, loc) result(y)
        integer, intent(in) :: k !! lower limit of P(X > k)
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        y = geometric_sf_real(real(k, dp), p, loc)
    end function geometric_sf_int

    pure elemental function geometric_logcdf_real(k, p, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of log(P(X <= k)); the floor is used
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        real(dp) :: cdf
        real(dp) :: logsf
        real(dp) :: sf

        call tails(floor_count(k, loc), p, cdf, sf, y, logsf)
    end function geometric_logcdf_real

    pure elemental function geometric_logcdf_int(k, p, loc) result(y)
        integer, intent(in) :: k !! upper limit of log(P(X <= k))
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        y = geometric_logcdf_real(real(k, dp), p, loc)
    end function geometric_logcdf_int

    pure elemental function geometric_logsf_real(k, p, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of log(P(X > k)); the floor is used
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        real(dp) :: cdf
        real(dp) :: logcdf
        real(dp) :: sf

        call tails(floor_count(k, loc), p, cdf, sf, logcdf, y)
    end function geometric_logsf_real

    pure elemental function geometric_logsf_int(k, p, loc) result(y)
        integer, intent(in) :: k !! lower limit of log(P(X > k))
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        y = geometric_logsf_real(real(k, dp), p, loc)
    end function geometric_logsf_int

    pure elemental function geometric_ppf(probability, p, loc) result(y)
        real(dp), intent(in) :: probability !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        real(dp) :: count
        real(dp) :: shift

        shift = 0.0_dp
        if (present(loc)) shift = loc
        if (.not. valid_probability(p) .or. .not. ieee_is_finite(shift)) then
            y = quiet_nan(probability)
        else if (.not. (probability >= 0.0_dp .and. probability <= 1.0_dp)) then
            y = quiet_nan(probability)
        else if (probability == 0.0_dp) then
            y = shift
        else if (probability >= 1.0_dp) then
            y = positive_infinity(probability)
        else if (p == 1.0_dp) then
            y = shift + 1.0_dp
        else
            count = ceil_positive(log1p_safe(-probability) / log1p_safe(-p))
            count = max(1.0_dp, count)
            if (count > 1.0_dp) then
                if (geometric_cdf_real(count - 1.0_dp, p) >= probability) &
                    count = count - 1.0_dp
            end if
            y = shift + count
        end if
    end function geometric_ppf

    pure elemental function geometric_isf(probability, p, loc) result(y)
        real(dp), intent(in) :: probability !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        real(dp) :: count
        real(dp) :: shift

        shift = 0.0_dp
        if (present(loc)) shift = loc
        if (.not. valid_probability(p) .or. .not. ieee_is_finite(shift)) then
            y = quiet_nan(probability)
        else if (.not. (probability >= 0.0_dp .and. probability <= 1.0_dp)) then
            y = quiet_nan(probability)
        else if (probability >= 1.0_dp) then
            y = shift
        else if (probability == 0.0_dp) then
            y = positive_infinity(probability)
        else if (p == 1.0_dp) then
            y = shift + 1.0_dp
        else
            count = ceil_positive(log(probability) / log1p_safe(-p))
            count = max(1.0_dp, count)
            if (count > 1.0_dp) then
                if (geometric_sf_real(count - 1.0_dp, p) <= probability) &
                    count = count - 1.0_dp
            end if
            y = shift + count
        end if
    end function geometric_isf

    pure elemental subroutine tails(count, p, cdf, sf, logcdf, logsf)
        real(dp), intent(in) :: count !! floor(k - loc)
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(out) :: cdf !! P(X <= count)
        real(dp), intent(out) :: sf !! P(X > count)
        real(dp), intent(out) :: logcdf !! log(cdf)
        real(dp), intent(out) :: logsf !! log(sf)

        if (.not. valid_probability(p) .or. ieee_is_nan(count)) then
            cdf = quiet_nan(p)
            sf = cdf
            logcdf = cdf
            logsf = cdf
        else if (count < 1.0_dp) then
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
            logsf = count * log1p_safe(-p)
            sf = exp(logsf)
            cdf = -expm1_safe(logsf)
            logcdf = log(cdf)
        end if
    end subroutine tails

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

    pure elemental function ceil_positive(x) result(y)
        real(dp), intent(in) :: x !! nonnegative value to round upward without integer conversion
        real(dp) :: y

        y = aint(x)
        if (y < x) y = y + 1.0_dp
    end function ceil_positive

    pure elemental logical function valid_probability(p) result(valid)
        real(dp), intent(in) :: p !! success probability to check

        valid = p > 0.0_dp .and. p <= 1.0_dp
    end function valid_probability

end module scifort_geometric
