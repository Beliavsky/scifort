! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Zipf/zeta distribution matching scipy.stats.zipf.
module scifort_zipf
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : floor_real, expm1_safe, negative_infinity, positive_infinity, quiet_nan
    use scifort_zeta, only : hurwitz_zeta
    implicit none
    private

    public :: zipf_cdf, zipf_isf, zipf_logcdf, zipf_logpmf
    public :: zipf_logsf, zipf_pmf, zipf_ppf, zipf_sf

    interface zipf_pmf
        module procedure zipf_pmf_real
        module procedure zipf_pmf_int
    end interface
    interface zipf_logpmf
        module procedure zipf_logpmf_real
        module procedure zipf_logpmf_int
    end interface
    interface zipf_cdf
        module procedure zipf_cdf_real
        module procedure zipf_cdf_int
    end interface
    interface zipf_sf
        module procedure zipf_sf_real
        module procedure zipf_sf_int
    end interface
    interface zipf_logcdf
        module procedure zipf_logcdf_real
        module procedure zipf_logcdf_int
    end interface
    interface zipf_logsf
        module procedure zipf_logsf_real
        module procedure zipf_logsf_int
    end interface
contains
    pure elemental function zipf_logpmf_real(k, a, loc) result(y)
        real(dp), intent(in) :: k !! lattice observation
        real(dp), intent(in) :: a !! power exponent, strictly greater than one
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, count
        count = k - optional_loc(loc)
        if (.not. valid_a(a) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (.not. ieee_is_finite(count) .or. count < 1.0_dp .or. count /= aint(count)) then
            y = negative_infinity(k)
        else
            y = -a * log(count) - log(hurwitz_zeta(a, 1.0_dp))
        end if
    end function
    pure elemental function zipf_logpmf_int(k, a, loc) result(y)
        integer, intent(in) :: k !! integer observation
        real(dp), intent(in) :: a !! power exponent, strictly greater than one
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = zipf_logpmf_real(real(k, dp), a, loc)
    end function
    pure elemental function zipf_pmf_real(k, a, loc) result(y)
        real(dp), intent(in) :: k !! lattice observation
        real(dp), intent(in) :: a !! power exponent, strictly greater than one
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = exp(zipf_logpmf_real(k, a, loc))
    end function
    pure elemental function zipf_pmf_int(k, a, loc) result(y)
        integer, intent(in) :: k !! integer observation
        real(dp), intent(in) :: a !! power exponent, strictly greater than one
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = zipf_pmf_real(real(k, dp), a, loc)
    end function
    pure elemental function zipf_logsf_real(k, a, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of log(P(X > k)); shifted floor is used
        real(dp), intent(in) :: a !! power exponent, strictly greater than one
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, count
        count = floor_real(k - optional_loc(loc))
        if (.not. valid_a(a) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (count < 1.0_dp) then
            y = 0.0_dp
        else if (.not. ieee_is_finite(count)) then
            y = negative_infinity(k)
        else
            y = log(hurwitz_zeta(a, count + 1.0_dp)) - log(hurwitz_zeta(a, 1.0_dp))
        end if
    end function
    pure elemental function zipf_logsf_int(k, a, loc) result(y)
        integer, intent(in) :: k !! lower integer limit
        real(dp), intent(in) :: a !! power exponent, strictly greater than one
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = zipf_logsf_real(real(k, dp), a, loc)
    end function
    pure elemental function zipf_sf_real(k, a, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of P(X > k)
        real(dp), intent(in) :: a !! power exponent, strictly greater than one
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, ly
        ly = zipf_logsf_real(k, a, loc)
        if (ieee_is_nan(ly)) then
            y = ly
        else
            y = exp(ly)
        end if
    end function
    pure elemental function zipf_sf_int(k, a, loc) result(y)
        integer, intent(in) :: k !! lower integer limit
        real(dp), intent(in) :: a !! power exponent, strictly greater than one
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = zipf_sf_real(real(k, dp), a, loc)
    end function
    pure elemental function zipf_logcdf_real(k, a, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of log(P(X <= k)); shifted floor is used
        real(dp), intent(in) :: a !! power exponent, strictly greater than one
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, lsf
        lsf = zipf_logsf_real(k, a, loc)
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
    pure elemental function zipf_logcdf_int(k, a, loc) result(y)
        integer, intent(in) :: k !! upper integer limit
        real(dp), intent(in) :: a !! power exponent, strictly greater than one
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = zipf_logcdf_real(real(k, dp), a, loc)
    end function
    pure elemental function zipf_cdf_real(k, a, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of P(X <= k)
        real(dp), intent(in) :: a !! power exponent, strictly greater than one
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, lsf
        lsf = zipf_logsf_real(k, a, loc)
        if (ieee_is_nan(lsf)) then
            y = lsf
        else
            y = -expm1_safe(lsf)
        end if
    end function
    pure elemental function zipf_cdf_int(k, a, loc) result(y)
        integer, intent(in) :: k !! upper integer limit
        real(dp), intent(in) :: a !! power exponent, strictly greater than one
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = zipf_cdf_real(real(k, dp), a, loc)
    end function
    pure elemental function zipf_ppf(probability, a, loc) result(y)
        real(dp), intent(in) :: probability !! lower-tail probability in [0,1]
        real(dp), intent(in) :: a !! power exponent, strictly greater than one
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, shift
        integer :: lo, hi, mid
        shift = optional_loc(loc)
        if (.not. valid_a(a) .or. .not. ieee_is_finite(shift) .or. .not. valid_probability(probability)) then
            y = quiet_nan(probability)
        else if (probability <= 0.0_dp) then
            y = shift
        else if (probability >= 1.0_dp) then
            y = positive_infinity(probability)
        else
            lo = 1; hi = 1
            do while (zipf_cdf_real(shift + real(hi, dp), a, shift) < probability .and. hi < huge(hi)/2)
                hi = hi * 2
            end do
            do while (lo < hi)
                mid = lo + (hi-lo)/2
                if (zipf_cdf_real(shift + real(mid, dp), a, shift) >= probability) then
                    hi = mid
                else
                    lo = mid + 1
                end if
            end do
            y = shift + real(lo, dp)
        end if
    end function
    pure elemental function zipf_isf(probability, a, loc) result(y)
        real(dp), intent(in) :: probability !! upper-tail probability in [0,1]
        real(dp), intent(in) :: a !! power exponent, strictly greater than one
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, shift
        integer :: lo, hi, mid
        shift = optional_loc(loc)
        if (.not. valid_a(a) .or. .not. ieee_is_finite(shift) .or. .not. valid_probability(probability)) then
            y = quiet_nan(probability)
        else if (probability >= 1.0_dp) then
            y = shift
        else if (probability <= 0.0_dp) then
            y = positive_infinity(probability)
        else
            lo = 1; hi = 1
            do while (zipf_sf_real(shift + real(hi, dp), a, shift) > probability .and. hi < huge(hi)/2)
                hi = hi * 2
            end do
            do while (lo < hi)
                mid = lo + (hi-lo)/2
                if (zipf_sf_real(shift + real(mid, dp), a, shift) <= probability) then
                    hi = mid
                else
                    lo = mid + 1
                end if
            end do
            y = shift + real(lo, dp)
        end if
    end function
    pure elemental function optional_loc(loc) result(value)
        real(dp), intent(in), optional :: loc !! support shift
        real(dp) :: value
        value = 0.0_dp
        if (present(loc)) value = loc
    end function
    pure elemental logical function valid_a(a)
        real(dp), intent(in) :: a !! exponent to validate
        valid_a = ieee_is_finite(a) .and. a > 1.0_dp
    end function
    pure elemental logical function valid_probability(p)
        real(dp), intent(in) :: p !! probability to validate
        valid_probability = p >= 0.0_dp .and. p <= 1.0_dp
    end function
end module scifort_zipf
