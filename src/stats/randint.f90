! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Discrete uniform distribution matching scipy.stats.randint. The standard
! support is low, ..., high - 1, with integer low < high. A real-valued loc
! shifts the lattice exactly as in scipy.stats.rv_discrete.

module scifort_randint
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : floor_real, negative_infinity, quiet_nan
    implicit none
    private

    public :: randint_cdf, randint_isf, randint_logcdf, randint_logpmf
    public :: randint_logsf, randint_pmf, randint_ppf, randint_sf

    interface randint_pmf
        module procedure randint_pmf_real
        module procedure randint_pmf_int
    end interface randint_pmf

    interface randint_logpmf
        module procedure randint_logpmf_real
        module procedure randint_logpmf_int
    end interface randint_logpmf

    interface randint_cdf
        module procedure randint_cdf_real
        module procedure randint_cdf_int
    end interface randint_cdf

    interface randint_sf
        module procedure randint_sf_real
        module procedure randint_sf_int
    end interface randint_sf

    interface randint_logcdf
        module procedure randint_logcdf_real
        module procedure randint_logcdf_int
    end interface randint_logcdf

    interface randint_logsf
        module procedure randint_logsf_real
        module procedure randint_logsf_int
    end interface randint_logsf

contains

    pure elemental function randint_pmf_real(k, low, high, loc) result(y)
        real(dp), intent(in) :: k !! lattice point; non-lattice values have mass zero
        real(dp), intent(in) :: low !! integer lower shape bound, inclusive
        real(dp), intent(in) :: high !! integer upper shape bound, exclusive and > low
        real(dp), intent(in), optional :: loc !! real support shift (default 0)
        real(dp) :: y
        real(dp) :: count

        count = shifted_count(k, loc)
        if (.not. valid_parameters(low, high) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (.not. ieee_is_finite(count) .or. count /= aint(count)) then
            y = 0.0_dp
        else if (count < low .or. count >= high) then
            y = 0.0_dp
        else
            y = 1.0_dp / (high - low)
        end if
    end function randint_pmf_real

    pure elemental function randint_pmf_int(k, low, high, loc) result(y)
        integer, intent(in) :: k !! integer evaluation point
        integer, intent(in) :: low !! lower shape bound, inclusive
        integer, intent(in) :: high !! upper shape bound, exclusive and > low
        real(dp), intent(in), optional :: loc !! real support shift (default 0)
        real(dp) :: y
        y = randint_pmf_real(real(k, dp), real(low, dp), real(high, dp), loc)
    end function randint_pmf_int

    pure elemental function randint_logpmf_real(k, low, high, loc) result(y)
        real(dp), intent(in) :: k !! lattice point; non-lattice values have log-mass -Inf
        real(dp), intent(in) :: low !! integer lower shape bound, inclusive
        real(dp), intent(in) :: high !! integer upper shape bound, exclusive and > low
        real(dp), intent(in), optional :: loc !! real support shift (default 0)
        real(dp) :: y
        real(dp) :: count

        count = shifted_count(k, loc)
        if (.not. valid_parameters(low, high) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (.not. ieee_is_finite(count) .or. count /= aint(count)) then
            y = negative_infinity(k)
        else if (count < low .or. count >= high) then
            y = negative_infinity(k)
        else
            y = -log(high - low)
        end if
    end function randint_logpmf_real

    pure elemental function randint_logpmf_int(k, low, high, loc) result(y)
        integer, intent(in) :: k !! integer evaluation point
        integer, intent(in) :: low !! lower shape bound, inclusive
        integer, intent(in) :: high !! upper shape bound, exclusive and > low
        real(dp), intent(in), optional :: loc !! real support shift (default 0)
        real(dp) :: y
        y = randint_logpmf_real(real(k, dp), real(low, dp), real(high, dp), loc)
    end function randint_logpmf_int

    pure elemental function randint_cdf_real(k, low, high, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of P(X <= k); shifted lattice floor is used
        real(dp), intent(in) :: low !! integer lower shape bound, inclusive
        real(dp), intent(in) :: high !! integer upper shape bound, exclusive and > low
        real(dp), intent(in), optional :: loc !! real support shift (default 0)
        real(dp) :: y
        real(dp) :: count

        count = floor_count(k, loc)
        if (.not. valid_parameters(low, high) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (count < low) then
            y = 0.0_dp
        else if (count >= high - 1.0_dp) then
            y = 1.0_dp
        else
            y = (count - low + 1.0_dp) / (high - low)
        end if
    end function randint_cdf_real

    pure elemental function randint_cdf_int(k, low, high, loc) result(y)
        integer, intent(in) :: k !! upper integer limit
        integer, intent(in) :: low !! lower shape bound, inclusive
        integer, intent(in) :: high !! upper shape bound, exclusive and > low
        real(dp), intent(in), optional :: loc !! real support shift (default 0)
        real(dp) :: y
        y = randint_cdf_real(real(k, dp), real(low, dp), real(high, dp), loc)
    end function randint_cdf_int

    pure elemental function randint_sf_real(k, low, high, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of P(X > k); shifted lattice floor is used
        real(dp), intent(in) :: low !! integer lower shape bound, inclusive
        real(dp), intent(in) :: high !! integer upper shape bound, exclusive and > low
        real(dp), intent(in), optional :: loc !! real support shift (default 0)
        real(dp) :: y
        real(dp) :: count

        count = floor_count(k, loc)
        if (.not. valid_parameters(low, high) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (count < low) then
            y = 1.0_dp
        else if (count >= high - 1.0_dp) then
            y = 0.0_dp
        else
            y = (high - count - 1.0_dp) / (high - low)
        end if
    end function randint_sf_real

    pure elemental function randint_sf_int(k, low, high, loc) result(y)
        integer, intent(in) :: k !! lower integer limit
        integer, intent(in) :: low !! lower shape bound, inclusive
        integer, intent(in) :: high !! upper shape bound, exclusive and > low
        real(dp), intent(in), optional :: loc !! real support shift (default 0)
        real(dp) :: y
        y = randint_sf_real(real(k, dp), real(low, dp), real(high, dp), loc)
    end function randint_sf_int

    pure elemental function randint_logcdf_real(k, low, high, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of log(P(X <= k))
        real(dp), intent(in) :: low !! integer lower shape bound, inclusive
        real(dp), intent(in) :: high !! integer upper shape bound, exclusive and > low
        real(dp), intent(in), optional :: loc !! real support shift (default 0)
        real(dp) :: y
        real(dp) :: p
        p = randint_cdf_real(k, low, high, loc)
        if (ieee_is_nan(p)) then
            y = p
        else if (p <= 0.0_dp) then
            y = negative_infinity(k)
        else
            y = log(p)
        end if
    end function randint_logcdf_real

    pure elemental function randint_logcdf_int(k, low, high, loc) result(y)
        integer, intent(in) :: k !! upper integer limit
        integer, intent(in) :: low !! lower shape bound, inclusive
        integer, intent(in) :: high !! upper shape bound, exclusive and > low
        real(dp), intent(in), optional :: loc !! real support shift (default 0)
        real(dp) :: y
        y = randint_logcdf_real(real(k, dp), real(low, dp), real(high, dp), loc)
    end function randint_logcdf_int

    pure elemental function randint_logsf_real(k, low, high, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of log(P(X > k))
        real(dp), intent(in) :: low !! integer lower shape bound, inclusive
        real(dp), intent(in) :: high !! integer upper shape bound, exclusive and > low
        real(dp), intent(in), optional :: loc !! real support shift (default 0)
        real(dp) :: y
        real(dp) :: p
        p = randint_sf_real(k, low, high, loc)
        if (ieee_is_nan(p)) then
            y = p
        else if (p <= 0.0_dp) then
            y = negative_infinity(k)
        else
            y = log(p)
        end if
    end function randint_logsf_real

    pure elemental function randint_logsf_int(k, low, high, loc) result(y)
        integer, intent(in) :: k !! lower integer limit
        integer, intent(in) :: low !! lower shape bound, inclusive
        integer, intent(in) :: high !! upper shape bound, exclusive and > low
        real(dp), intent(in), optional :: loc !! real support shift (default 0)
        real(dp) :: y
        y = randint_logsf_real(real(k, dp), real(low, dp), real(high, dp), loc)
    end function randint_logsf_int

    pure elemental function randint_ppf(probability, low, high, loc) result(y)
        real(dp), intent(in) :: probability !! lower-tail probability in [0,1]
        real(dp), intent(in) :: low !! integer lower shape bound, inclusive
        real(dp), intent(in) :: high !! integer upper shape bound, exclusive and > low
        real(dp), intent(in), optional :: loc !! real support shift (default 0)
        real(dp) :: y
        real(dp) :: shift, width, count

        shift = optional_loc(loc)
        if (.not. valid_parameters(low, high) .or. .not. ieee_is_finite(shift) .or. &
            .not. (probability >= 0.0_dp .and. probability <= 1.0_dp)) then
            y = quiet_nan(probability)
        else if (probability <= 0.0_dp) then
            y = shift + low - 1.0_dp
        else if (probability >= 1.0_dp) then
            y = shift + high - 1.0_dp
        else
            width = high - low
            count = -floor_real(-(probability * width + low)) - 1.0_dp
            count = max(low, min(high - 1.0_dp, count))
            if (count > low) then
                if (randint_cdf_real(shift + count - 1.0_dp, low, high, shift) >= probability) &
                    count = count - 1.0_dp
            end if
            y = shift + count
        end if
    end function randint_ppf

    pure elemental function randint_isf(probability, low, high, loc) result(y)
        real(dp), intent(in) :: probability !! upper-tail probability in [0,1]
        real(dp), intent(in) :: low !! integer lower shape bound, inclusive
        real(dp), intent(in) :: high !! integer upper shape bound, exclusive and > low
        real(dp), intent(in), optional :: loc !! real support shift (default 0)
        real(dp) :: y
        real(dp) :: shift, width, count

        shift = optional_loc(loc)
        if (.not. valid_parameters(low, high) .or. .not. ieee_is_finite(shift) .or. &
            .not. (probability >= 0.0_dp .and. probability <= 1.0_dp)) then
            y = quiet_nan(probability)
        else if (probability <= 0.0_dp) then
            y = shift + high - 1.0_dp
        else if (probability >= 1.0_dp) then
            y = shift + low - 1.0_dp
        else
            width = high - low
            count = floor_real(probability * width)
            count = high - 1.0_dp - count
            count = max(low, min(high - 1.0_dp, count))
            do while (count > low .and. randint_sf_real(shift + count - 1.0_dp, low, high, shift) <= probability)
                count = count - 1.0_dp
            end do
            do while (count < high - 1.0_dp .and. randint_sf_real(shift + count, low, high, shift) > probability)
                count = count + 1.0_dp
            end do
            y = shift + count
        end if
    end function randint_isf

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

    pure elemental logical function valid_parameters(low, high)
        real(dp), intent(in) :: low !! lower shape bound
        real(dp), intent(in) :: high !! upper shape bound
        valid_parameters = ieee_is_finite(low) .and. ieee_is_finite(high) .and. &
            low == aint(low) .and. high == aint(high) .and. high > low
    end function valid_parameters

end module scifort_randint
