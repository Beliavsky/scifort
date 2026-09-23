! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Logarithmic-series distribution matching scipy.stats.logser.
! P(X=k) = -p**k/(k*log(1-p)), k = 1,2,... before loc.
! Finite lower-tail sums and direct upper-tail recurrence avoid subtracting a
! nearly unit CDF from one. The upper-tail recurrence is exact in binary64 up
! to its convergence tolerance and has a generous deterministic term cap.

module scifort_logser
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, quiet_nan
    implicit none
    private

    integer, parameter :: max_tail_terms = 2000000

    public :: logser_cdf, logser_isf, logser_logcdf, logser_logpmf
    public :: logser_logsf, logser_pmf, logser_ppf, logser_sf

    interface logser_pmf
        module procedure logser_pmf_real
        module procedure logser_pmf_int
    end interface logser_pmf

    interface logser_logpmf
        module procedure logser_logpmf_real
        module procedure logser_logpmf_int
    end interface logser_logpmf

    interface logser_cdf
        module procedure logser_cdf_real
        module procedure logser_cdf_int
    end interface logser_cdf

    interface logser_sf
        module procedure logser_sf_real
        module procedure logser_sf_int
    end interface logser_sf

    interface logser_logcdf
        module procedure logser_logcdf_real
        module procedure logser_logcdf_int
    end interface logser_logcdf

    interface logser_logsf
        module procedure logser_logsf_real
        module procedure logser_logsf_int
    end interface logser_logsf

contains

    pure elemental function logser_logpmf_real(k, p, loc) result(y)
        real(dp), intent(in) :: k !! lattice point; non-lattice values have log-mass -Inf
        real(dp), intent(in) :: p !! probability parameter strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: count

        count = shifted_count(k, loc)
        if (.not. valid_p(p) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (.not. ieee_is_finite(count) .or. count < 1.0_dp .or. count /= aint(count)) then
            y = negative_infinity(k)
        else
            y = count * log(p) - log(count) - log(-log1p_safe(-p))
        end if
    end function logser_logpmf_real

    pure elemental function logser_logpmf_int(k, p, loc) result(y)
        integer, intent(in) :: k !! integer evaluation point
        real(dp), intent(in) :: p !! probability parameter strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = logser_logpmf_real(real(k, dp), p, loc)
    end function logser_logpmf_int

    pure elemental function logser_pmf_real(k, p, loc) result(y)
        real(dp), intent(in) :: k !! lattice point
        real(dp), intent(in) :: p !! probability parameter strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = exp(logser_logpmf_real(k, p, loc))
    end function logser_pmf_real

    pure elemental function logser_pmf_int(k, p, loc) result(y)
        integer, intent(in) :: k !! integer evaluation point
        real(dp), intent(in) :: p !! probability parameter strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = logser_pmf_real(real(k, dp), p, loc)
    end function logser_pmf_int

    pure elemental function logser_cdf_real(k, p, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of P(X <= k); shifted floor is used
        real(dp), intent(in) :: p !! probability parameter strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: count, cdf, sf

        count = floor_count(k, loc)
        if (.not. valid_p(p) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (count < 1.0_dp) then
            y = 0.0_dp
        else if (.not. ieee_is_finite(count)) then
            y = 1.0_dp
        else
            call tails(int(count), p, cdf, sf)
            y = cdf
        end if
    end function logser_cdf_real

    pure elemental function logser_cdf_int(k, p, loc) result(y)
        integer, intent(in) :: k !! upper integer limit
        real(dp), intent(in) :: p !! probability parameter strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = logser_cdf_real(real(k, dp), p, loc)
    end function logser_cdf_int

    pure elemental function logser_sf_real(k, p, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of P(X > k); shifted floor is used
        real(dp), intent(in) :: p !! probability parameter strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: count, cdf, sf

        count = floor_count(k, loc)
        if (.not. valid_p(p) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (count < 1.0_dp) then
            y = 1.0_dp
        else if (.not. ieee_is_finite(count)) then
            y = 0.0_dp
        else
            call tails(int(count), p, cdf, sf)
            y = sf
        end if
    end function logser_sf_real

    pure elemental function logser_sf_int(k, p, loc) result(y)
        integer, intent(in) :: k !! lower integer limit
        real(dp), intent(in) :: p !! probability parameter strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = logser_sf_real(real(k, dp), p, loc)
    end function logser_sf_int

    pure elemental function logser_logcdf_real(k, p, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of log(P(X <= k))
        real(dp), intent(in) :: p !! probability parameter strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: value
        value = logser_cdf_real(k, p, loc)
        if (ieee_is_nan(value)) then
            y = value
        else if (value <= 0.0_dp) then
            y = negative_infinity(k)
        else
            y = log(value)
        end if
    end function logser_logcdf_real

    pure elemental function logser_logcdf_int(k, p, loc) result(y)
        integer, intent(in) :: k !! upper integer limit
        real(dp), intent(in) :: p !! probability parameter strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = logser_logcdf_real(real(k, dp), p, loc)
    end function logser_logcdf_int

    pure elemental function logser_logsf_real(k, p, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of log(P(X > k))
        real(dp), intent(in) :: p !! probability parameter strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: value
        value = logser_sf_real(k, p, loc)
        if (ieee_is_nan(value)) then
            y = value
        else if (value <= 0.0_dp) then
            y = negative_infinity(k)
        else
            y = log(value)
        end if
    end function logser_logsf_real

    pure elemental function logser_logsf_int(k, p, loc) result(y)
        integer, intent(in) :: k !! lower integer limit
        real(dp), intent(in) :: p !! probability parameter strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = logser_logsf_real(real(k, dp), p, loc)
    end function logser_logsf_int

    pure elemental function logser_ppf(probability, p, loc) result(y)
        real(dp), intent(in) :: probability !! lower-tail probability in [0,1]
        real(dp), intent(in) :: p !! probability parameter strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: shift
        integer :: lo, hi, mid

        shift = optional_loc(loc)
        if (.not. valid_p(p) .or. .not. ieee_is_finite(shift) .or. &
            .not. (probability >= 0.0_dp .and. probability <= 1.0_dp)) then
            y = quiet_nan(probability)
        else if (probability <= 0.0_dp) then
            y = shift
        else if (probability >= 1.0_dp) then
            y = positive_infinity(probability)
        else
            lo = 1
            hi = 1
            do while (logser_cdf_real(shift + real(hi, dp), p, shift) < probability)
                lo = hi + 1
                if (hi > huge(hi) / 2) then
                    y = positive_infinity(probability)
                    return
                end if
                hi = 2 * hi
            end do
            do while (lo < hi)
                mid = lo + (hi - lo) / 2
                if (logser_cdf_real(shift + real(mid, dp), p, shift) >= probability) then
                    hi = mid
                else
                    lo = mid + 1
                end if
            end do
            y = shift + real(lo, dp)
        end if
    end function logser_ppf

    pure elemental function logser_isf(probability, p, loc) result(y)
        real(dp), intent(in) :: probability !! upper-tail probability in [0,1]
        real(dp), intent(in) :: p !! probability parameter strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        real(dp) :: shift
        integer :: lo, hi, mid

        shift = optional_loc(loc)
        if (.not. valid_p(p) .or. .not. ieee_is_finite(shift) .or. &
            .not. (probability >= 0.0_dp .and. probability <= 1.0_dp)) then
            y = quiet_nan(probability)
        else if (probability >= 1.0_dp) then
            y = shift
        else if (probability <= 0.0_dp) then
            y = positive_infinity(probability)
        else
            lo = 1
            hi = 1
            do while (logser_sf_real(shift + real(hi, dp), p, shift) > probability)
                lo = hi + 1
                if (hi > huge(hi) / 2) then
                    y = positive_infinity(probability)
                    return
                end if
                hi = 2 * hi
            end do
            do while (lo < hi)
                mid = lo + (hi - lo) / 2
                if (logser_sf_real(shift + real(mid, dp), p, shift) <= probability) then
                    hi = mid
                else
                    lo = mid + 1
                end if
            end do
            y = shift + real(lo, dp)
        end if
    end function logser_isf

    pure subroutine tails(k, p, cdf, sf)
        integer, intent(in) :: k !! positive integer tail split point
        real(dp), intent(in) :: p !! probability parameter in (0,1)
        real(dp), intent(out) :: cdf !! P(X <= k)
        real(dp), intent(out) :: sf !! P(X > k)

        real(dp) :: norm, term, lower_sum, upper_sum
        integer :: j, n

        norm = -log1p_safe(-p)
        if (k <= max(32, int(4.0_dp * p / max(1.0e-12_dp, 1.0_dp - p)))) then
            term = p
            lower_sum = term
            do j = 2, k
                term = term * p * real(j - 1, dp) / real(j, dp)
                lower_sum = lower_sum + term
            end do
            cdf = min(1.0_dp, lower_sum / norm)
            if (cdf <= 0.5_dp) then
                sf = 1.0_dp - cdf
                return
            end if
        end if

        n = k + 1
        term = exp(real(n, dp) * log(p)) / real(n, dp)
        upper_sum = 0.0_dp
        do j = 0, max_tail_terms - 1
            upper_sum = upper_sum + term
            if (term <= epsilon(1.0_dp) * max(upper_sum, tiny(1.0_dp))) exit
            n = n + 1
            term = term * p * real(n - 1, dp) / real(n, dp)
            if (term == 0.0_dp) exit
        end do
        sf = max(0.0_dp, min(1.0_dp, upper_sum / norm))
        cdf = 1.0_dp - sf
        if (cdf < 0.0_dp) cdf = 0.0_dp
    end subroutine tails

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

    pure elemental logical function valid_p(p)
        real(dp), intent(in) :: p !! probability parameter
        valid_p = ieee_is_finite(p) .and. p > 0.0_dp .and. p < 1.0_dp
    end function valid_p

end module scifort_logser
