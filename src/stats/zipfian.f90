! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Finite Zipfian distribution matching scipy.stats.zipfian.
module scifort_zipfian
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, quiet_nan
    implicit none
    private

    public :: zipfian_cdf, zipfian_isf, zipfian_logcdf, zipfian_logpmf
    public :: zipfian_logsf, zipfian_pmf, zipfian_ppf, zipfian_sf

    interface zipfian_pmf
        module procedure zipfian_pmf_real
        module procedure zipfian_pmf_int
    end interface
    interface zipfian_logpmf
        module procedure zipfian_logpmf_real
        module procedure zipfian_logpmf_int
    end interface
    interface zipfian_cdf
        module procedure zipfian_cdf_real
        module procedure zipfian_cdf_int
    end interface
    interface zipfian_sf
        module procedure zipfian_sf_real
        module procedure zipfian_sf_int
    end interface
    interface zipfian_logcdf
        module procedure zipfian_logcdf_real
        module procedure zipfian_logcdf_int
    end interface
    interface zipfian_logsf
        module procedure zipfian_logsf_real
        module procedure zipfian_logsf_int
    end interface
contains
    pure elemental function zipfian_logpmf_real(k, a, n, loc) result(y)
        real(dp), intent(in) :: k !! lattice observation
        real(dp), intent(in) :: a !! nonnegative power exponent
        real(dp), intent(in) :: n !! positive integer upper support bound
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, count, h
        count = k - optional_loc(loc)
        if (.not. valid_parameters(a, n) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (.not. ieee_is_finite(count) .or. count < 1.0_dp .or. count > n .or. count /= aint(count)) then
            y = negative_infinity(k)
        else
            h = harmonic_sum(int(n), a)
            y = -a * log(count) - log(h)
        end if
    end function
    pure elemental function zipfian_logpmf_int(k, a, n, loc) result(y)
        integer, intent(in) :: k !! integer observation
        real(dp), intent(in) :: a !! nonnegative power exponent
        integer, intent(in) :: n !! positive integer upper support bound
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = zipfian_logpmf_real(real(k, dp), a, real(n, dp), loc)
    end function
    pure elemental function zipfian_pmf_real(k, a, n, loc) result(y)
        real(dp), intent(in) :: k !! lattice observation
        real(dp), intent(in) :: a !! nonnegative power exponent
        real(dp), intent(in) :: n !! positive integer upper support bound
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = exp(zipfian_logpmf_real(k, a, n, loc))
    end function
    pure elemental function zipfian_pmf_int(k, a, n, loc) result(y)
        integer, intent(in) :: k !! integer observation
        real(dp), intent(in) :: a !! nonnegative power exponent
        integer, intent(in) :: n !! positive integer upper support bound
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = zipfian_pmf_real(real(k, dp), a, real(n, dp), loc)
    end function
    pure elemental function zipfian_cdf_real(k, a, n, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of P(X <= k); shifted floor is used
        real(dp), intent(in) :: a !! nonnegative power exponent
        real(dp), intent(in) :: n !! positive integer upper support bound
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, count, total
        count = floor(k - optional_loc(loc))
        if (.not. valid_parameters(a, n) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (count < 1.0_dp) then
            y = 0.0_dp
        else if (.not. ieee_is_finite(count) .or. count >= n) then
            y = 1.0_dp
        else
            total = harmonic_sum(int(n), a)
            y = harmonic_sum(int(count), a) / total
            y = min(1.0_dp, max(0.0_dp, y))
        end if
    end function
    pure elemental function zipfian_cdf_int(k, a, n, loc) result(y)
        integer, intent(in) :: k !! upper integer limit
        real(dp), intent(in) :: a !! nonnegative power exponent
        integer, intent(in) :: n !! positive integer upper support bound
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = zipfian_cdf_real(real(k, dp), a, real(n, dp), loc)
    end function
    pure elemental function zipfian_sf_real(k, a, n, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of P(X > k); shifted floor is used
        real(dp), intent(in) :: a !! nonnegative power exponent
        real(dp), intent(in) :: n !! positive integer upper support bound
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, count, total
        count = floor(k - optional_loc(loc))
        if (.not. valid_parameters(a, n) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (count < 1.0_dp) then
            y = 1.0_dp
        else if (.not. ieee_is_finite(count) .or. count >= n) then
            y = 0.0_dp
        else
            total = harmonic_sum(int(n), a)
            y = harmonic_range(int(count) + 1, int(n), a) / total
            y = min(1.0_dp, max(0.0_dp, y))
        end if
    end function
    pure elemental function zipfian_sf_int(k, a, n, loc) result(y)
        integer, intent(in) :: k !! lower integer limit
        real(dp), intent(in) :: a !! nonnegative power exponent
        integer, intent(in) :: n !! positive integer upper support bound
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = zipfian_sf_real(real(k, dp), a, real(n, dp), loc)
    end function
    pure elemental function zipfian_logcdf_real(k, a, n, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of log(P(X <= k))
        real(dp), intent(in) :: a !! nonnegative power exponent
        real(dp), intent(in) :: n !! positive integer upper support bound
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, p
        p = zipfian_cdf_real(k, a, n, loc)
        if (ieee_is_nan(p)) then
            y = p
        else if (p <= 0.0_dp) then
            y = negative_infinity(k)
        else
            y = log(p)
        end if
    end function
    pure elemental function zipfian_logcdf_int(k, a, n, loc) result(y)
        integer, intent(in) :: k !! upper integer limit
        real(dp), intent(in) :: a !! nonnegative power exponent
        integer, intent(in) :: n !! positive integer upper support bound
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = zipfian_logcdf_real(real(k, dp), a, real(n, dp), loc)
    end function
    pure elemental function zipfian_logsf_real(k, a, n, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of log(P(X > k))
        real(dp), intent(in) :: a !! nonnegative power exponent
        real(dp), intent(in) :: n !! positive integer upper support bound
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, p
        p = zipfian_sf_real(k, a, n, loc)
        if (ieee_is_nan(p)) then
            y = p
        else if (p <= 0.0_dp) then
            y = negative_infinity(k)
        else
            y = log(p)
        end if
    end function
    pure elemental function zipfian_logsf_int(k, a, n, loc) result(y)
        integer, intent(in) :: k !! lower integer limit
        real(dp), intent(in) :: a !! nonnegative power exponent
        integer, intent(in) :: n !! positive integer upper support bound
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = zipfian_logsf_real(real(k, dp), a, real(n, dp), loc)
    end function
    pure elemental function zipfian_ppf(probability, a, n, loc) result(y)
        real(dp), intent(in) :: probability !! lower-tail probability in [0,1]
        real(dp), intent(in) :: a !! nonnegative power exponent
        real(dp), intent(in) :: n !! positive integer upper support bound
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, shift
        integer :: lo, hi, mid
        shift = optional_loc(loc)
        if (.not. valid_parameters(a, n) .or. .not. ieee_is_finite(shift) .or. .not. valid_probability(probability)) then
            y = quiet_nan(probability)
        else if (probability <= 0.0_dp) then
            y = shift
        else if (probability >= 1.0_dp) then
            y = shift + n
        else
            lo = 1; hi = int(n)
            do while (lo < hi)
                mid = lo + (hi-lo)/2
                if (zipfian_cdf_real(shift + real(mid, dp), a, n, shift) >= probability) then
                    hi = mid
                else
                    lo = mid + 1
                end if
            end do
            y = shift + real(lo, dp)
        end if
    end function
    pure elemental function zipfian_isf(probability, a, n, loc) result(y)
        real(dp), intent(in) :: probability !! upper-tail probability in [0,1]
        real(dp), intent(in) :: a !! nonnegative power exponent
        real(dp), intent(in) :: n !! positive integer upper support bound
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, shift
        integer :: lo, hi, mid
        shift = optional_loc(loc)
        if (.not. valid_parameters(a, n) .or. .not. ieee_is_finite(shift) .or. .not. valid_probability(probability)) then
            y = quiet_nan(probability)
        else if (probability >= 1.0_dp) then
            y = shift
        else if (probability <= 0.0_dp) then
            y = shift + n
        else
            lo = 1; hi = int(n)
            do while (lo < hi)
                mid = lo + (hi-lo)/2
                if (zipfian_sf_real(shift + real(mid, dp), a, n, shift) <= probability) then
                    hi = mid
                else
                    lo = mid + 1
                end if
            end do
            y = shift + real(lo, dp)
        end if
    end function
    pure function harmonic_sum(n, a) result(value)
        integer, intent(in) :: n !! upper summation index
        real(dp), intent(in) :: a !! nonnegative exponent
        real(dp) :: value
        integer :: j
        value = 0.0_dp
        if (a == 0.0_dp) then
            value = real(n, dp)
            return
        end if
        do j = 1, n
            value = value + exp(-a * log(real(j, dp)))
        end do
    end function
    pure function harmonic_range(first, last, a) result(value)
        integer, intent(in) :: first !! first included index
        integer, intent(in) :: last !! last included index
        real(dp), intent(in) :: a !! nonnegative exponent
        real(dp) :: value
        integer :: j
        value = 0.0_dp
        if (first > last) return
        if (a == 0.0_dp) then
            value = real(last-first+1, dp)
            return
        end if
        do j = first, last
            value = value + exp(-a * log(real(j, dp)))
        end do
    end function
    pure elemental function optional_loc(loc) result(value)
        real(dp), intent(in), optional :: loc !! support shift
        real(dp) :: value
        value = 0.0_dp
        if (present(loc)) value = loc
    end function
    pure elemental logical function valid_parameters(a, n)
        real(dp), intent(in) :: a !! exponent to validate
        real(dp), intent(in) :: n !! upper support bound to validate
        valid_parameters = ieee_is_finite(a) .and. ieee_is_finite(n) .and. &
            a >= 0.0_dp .and. n >= 1.0_dp .and. n == aint(n) .and. n <= real(huge(0), dp)
    end function
    pure elemental logical function valid_probability(p)
        real(dp), intent(in) :: p !! probability to validate
        valid_probability = p >= 0.0_dp .and. p <= 1.0_dp
    end function
end module scifort_zipfian
