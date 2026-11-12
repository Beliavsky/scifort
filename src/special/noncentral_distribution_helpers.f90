! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! SciPy-compatible legacy noncentral-distribution helpers.  The forward CDFs
! and variate inverses reuse SciFort's independently validated distribution
! kernels; parameter inverses use safeguarded scalar root searches around the
! same kernels.
module scifort_special_noncentral_distribution_helpers
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan
    use scifort_ncf, only : ncf_cdf, ncf_ppf
    use scifort_nct, only : nct_cdf, nct_ppf
    use scifort_ncx2, only : ncx2_cdf, ncx2_ppf
    use scifort_special_elementary, only : ndtr
    implicit none
    private

    public :: chndtr
    public :: chndtridf
    public :: chndtrinc
    public :: chndtrix
    public :: ncfdtr
    public :: ncfdtri
    public :: ncfdtridfd
    public :: ncfdtridfn
    public :: ncfdtrinc
    public :: nctdtr
    public :: nctdtridf
    public :: nctdtrinc
    public :: nctdtrit

    integer, parameter :: family_chi_df = 1
    integer, parameter :: family_chi_nc = 2
    integer, parameter :: family_ncf_dfd = 3
    integer, parameter :: family_ncf_dfn = 4
    integer, parameter :: family_ncf_nc = 5
    integer, parameter :: family_nct_df = 6

contains

    pure elemental function chndtr(x, df, nc) result(value)
        real(dp), intent(in) :: x !! nonnegative chi-square variate; +infinity allowed
        real(dp), intent(in) :: df !! positive finite degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative finite noncentrality parameter
        real(dp) :: value

        if (.not. valid_positive_finite(df) .or. .not. valid_nonnegative_finite(nc) .or. &
                ieee_is_nan(x) .or. x < 0.0_dp) then
            value = quiet_nan(x)
        else
            value = ncx2_cdf(x, df, nc)
        end if
    end function chndtr

    pure elemental function chndtrix(p, df, nc) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: df !! positive finite degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative finite noncentrality parameter
        real(dp) :: x

        if (.not. valid_probability(p) .or. .not. valid_positive_finite(df) .or. &
                .not. valid_nonnegative_finite(nc)) then
            x = quiet_nan(p)
        else
            x = ncx2_ppf(p, df, nc)
        end if
    end function chndtrix

    pure elemental function chndtridf(x, p, nc) result(df)
        real(dp), intent(in) :: x !! positive finite chi-square variate
        real(dp), intent(in) :: p !! target lower-tail probability in [0, 1]
        real(dp), intent(in) :: nc !! nonnegative finite noncentrality parameter
        real(dp) :: df

        real(dp) :: near_zero

        if (.not. valid_positive_finite(x) .or. .not. valid_probability(p) .or. &
                .not. valid_nonnegative_finite(nc)) then
            df = quiet_nan(p)
            return
        end if
        if (p <= 0.0_dp) then
            df = positive_infinity(p)
            return
        end if
        near_zero = ncx2_cdf(x, 1.0e-12_dp, nc)
        if (p > near_zero) then
            df = quiet_nan(p)
        else if (abs(p - near_zero) <= 32.0_dp * epsilon(1.0_dp) * max(1.0_dp, p)) then
            df = 0.0_dp
        else
            df = solve_positive_monotone(p, family_chi_df, x, nc, 0.0_dp, .false.)
        end if
    end function chndtridf

    pure elemental function chndtrinc(x, df, p) result(nc)
        real(dp), intent(in) :: x !! positive finite chi-square variate
        real(dp), intent(in) :: df !! positive finite degrees of freedom
        real(dp), intent(in) :: p !! target lower-tail probability in [0, 1]
        real(dp) :: nc

        real(dp) :: central

        if (.not. valid_positive_finite(x) .or. .not. valid_positive_finite(df) .or. &
                .not. valid_probability(p)) then
            nc = quiet_nan(p)
            return
        end if
        central = ncx2_cdf(x, df, 0.0_dp)
        if (abs(p - central) <= probability_tolerance(central)) then
            nc = 0.0_dp
        else if (p > central) then
            nc = quiet_nan(p)
        else if (p <= 0.0_dp) then
            nc = positive_infinity(p)
        else
            nc = solve_nonnegative_monotone(p, family_chi_nc, x, df, 0.0_dp, .false.)
        end if
    end function chndtrinc

    pure elemental function ncfdtr(dfn, dfd, nc, f) result(value)
        real(dp), intent(in) :: dfn !! positive finite numerator degrees of freedom
        real(dp), intent(in) :: dfd !! positive finite denominator degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative finite noncentrality parameter
        real(dp), intent(in) :: f !! nonnegative F variate; +infinity allowed
        real(dp) :: value

        if (.not. valid_positive_finite(dfn) .or. .not. valid_positive_finite(dfd) .or. &
                .not. valid_nonnegative_finite(nc) .or. ieee_is_nan(f) .or. f < 0.0_dp) then
            value = quiet_nan(f)
        else
            value = ncf_cdf(f, dfn, dfd, nc)
        end if
    end function ncfdtr

    pure elemental function ncfdtri(dfn, dfd, nc, p) result(f)
        real(dp), intent(in) :: dfn !! positive finite numerator degrees of freedom
        real(dp), intent(in) :: dfd !! positive finite denominator degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative finite noncentrality parameter
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp) :: f

        if (.not. valid_positive_finite(dfn) .or. .not. valid_positive_finite(dfd) .or. &
                .not. valid_nonnegative_finite(nc) .or. .not. valid_probability(p)) then
            f = quiet_nan(p)
        else
            f = ncf_ppf(p, dfn, dfd, nc)
        end if
    end function ncfdtri

    pure elemental function ncfdtridfd(dfn, p, nc, f) result(dfd)
        real(dp), intent(in) :: dfn !! positive finite numerator degrees of freedom
        real(dp), intent(in) :: p !! target lower-tail probability in (0, 1)
        real(dp), intent(in) :: nc !! nonnegative finite noncentrality parameter
        real(dp), intent(in) :: f !! positive finite F variate
        real(dp) :: dfd

        if (.not. valid_positive_finite(dfn) .or. .not. valid_open_probability(p) .or. &
                .not. valid_nonnegative_finite(nc) .or. .not. valid_positive_finite(f)) then
            dfd = quiet_nan(p)
        else
            dfd = solve_first_positive_root(p, family_ncf_dfd, dfn, nc, f)
        end if
    end function ncfdtridfd

    pure elemental function ncfdtridfn(p, dfd, nc, f) result(dfn)
        real(dp), intent(in) :: p !! target lower-tail probability in (0, 1)
        real(dp), intent(in) :: dfd !! positive finite denominator degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative finite noncentrality parameter
        real(dp), intent(in) :: f !! positive finite F variate
        real(dp) :: dfn

        if (.not. valid_open_probability(p) .or. .not. valid_positive_finite(dfd) .or. &
                .not. valid_nonnegative_finite(nc) .or. .not. valid_positive_finite(f)) then
            dfn = quiet_nan(p)
        else
            dfn = solve_first_positive_root(p, family_ncf_dfn, dfd, nc, f)
        end if
    end function ncfdtridfn

    pure elemental function ncfdtrinc(dfn, dfd, p, f) result(nc)
        real(dp), intent(in) :: dfn !! positive finite numerator degrees of freedom
        real(dp), intent(in) :: dfd !! positive finite denominator degrees of freedom
        real(dp), intent(in) :: p !! target lower-tail probability in [0, 1]
        real(dp), intent(in) :: f !! positive finite F variate
        real(dp) :: nc

        real(dp) :: central

        if (.not. valid_positive_finite(dfn) .or. .not. valid_positive_finite(dfd) .or. &
                .not. valid_probability(p) .or. .not. valid_positive_finite(f)) then
            nc = quiet_nan(p)
            return
        end if
        central = ncf_cdf(f, dfn, dfd, 0.0_dp)
        if (abs(p - central) <= probability_tolerance(central)) then
            nc = 0.0_dp
        else if (p > central) then
            nc = quiet_nan(p)
        else if (p <= 0.0_dp) then
            nc = positive_infinity(p)
        else
            nc = solve_nonnegative_monotone(p, family_ncf_nc, dfn, dfd, f, .false.)
        end if
    end function ncfdtrinc

    pure elemental function nctdtr(df, nc, t) result(value)
        real(dp), intent(in) :: df !! positive finite degrees of freedom
        real(dp), intent(in) :: nc !! finite noncentrality parameter
        real(dp), intent(in) :: t !! Student-t variate; infinities allowed
        real(dp) :: value

        if (.not. valid_positive_finite(df) .or. .not. ieee_is_finite(nc) .or. &
                ieee_is_nan(t)) then
            value = quiet_nan(t)
        else
            value = nct_cdf(t, df, nc)
        end if
    end function nctdtr

    pure elemental function nctdtrit(df, nc, p) result(t)
        real(dp), intent(in) :: df !! positive finite degrees of freedom
        real(dp), intent(in) :: nc !! finite noncentrality parameter
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp) :: t

        if (.not. valid_positive_finite(df) .or. .not. ieee_is_finite(nc) .or. &
                .not. valid_probability(p)) then
            t = quiet_nan(p)
        else if (p <= 0.0_dp) then
            t = negative_infinity(p)
        else if (p >= 1.0_dp) then
            t = positive_infinity(p)
        else
            t = nct_ppf(p, df, nc)
        end if
    end function nctdtrit

    pure elemental function nctdtrinc(df, p, t) result(nc)
        real(dp), intent(in) :: df !! positive finite degrees of freedom
        real(dp), intent(in) :: p !! target lower-tail probability in [0, 1]
        real(dp), intent(in) :: t !! finite Student-t variate
        real(dp) :: nc

        if (.not. valid_positive_finite(df) .or. .not. valid_probability(p) .or. &
                .not. ieee_is_finite(t)) then
            nc = quiet_nan(p)
        else if (p <= 0.0_dp) then
            nc = positive_infinity(p)
        else if (p >= 1.0_dp) then
            nc = negative_infinity(p)
        else
            nc = solve_real_nct_nc(df, p, t)
        end if
    end function nctdtrinc

    pure elemental function nctdtridf(p, nc, t) result(df)
        real(dp), intent(in) :: p !! target lower-tail probability in (0, 1)
        real(dp), intent(in) :: nc !! finite noncentrality parameter
        real(dp), intent(in) :: t !! finite Student-t variate
        real(dp) :: df

        real(dp) :: limit_probability

        if (.not. valid_open_probability(p) .or. .not. ieee_is_finite(nc) .or. &
                .not. ieee_is_finite(t)) then
            df = quiet_nan(p)
            return
        end if
        limit_probability = ndtr(t - nc)
        if (abs(p - limit_probability) <= 16.0_dp * epsilon(1.0_dp) * &
                max(1.0_dp, abs(limit_probability))) then
            df = positive_infinity(p)
        else
            df = solve_first_positive_root(p, family_nct_df, nc, t, 0.0_dp)
        end if
    end function nctdtridf

    pure function solve_positive_monotone(target, family, arg1, arg2, arg3, increasing) result(root)
        real(dp), intent(in) :: target !! target probability strictly inside the attainable range
        integer, intent(in) :: family !! equation selector
        real(dp), intent(in) :: arg1 !! first fixed equation argument
        real(dp), intent(in) :: arg2 !! second fixed equation argument
        real(dp), intent(in) :: arg3 !! third fixed equation argument
        logical, intent(in) :: increasing !! true when the probability increases with the parameter
        real(dp) :: root

        integer :: iteration
        real(dp) :: hi
        real(dp) :: lo
        real(dp) :: value

        root = 1.0_dp
        value = family_probability(root, family, arg1, arg2, arg3)
        if (too_low(value, target, increasing)) then
            lo = root
            hi = 2.0_dp
            do iteration = 1, 1024
                value = family_probability(hi, family, arg1, arg2, arg3)
                if (.not. too_low(value, target, increasing)) exit
                lo = hi
                if (hi > huge(1.0_dp) / 2.0_dp) then
                    root = positive_infinity(target)
                    return
                end if
                hi = 2.0_dp * hi
            end do
        else
            hi = root
            lo = 0.5_dp
            do iteration = 1, 1024
                value = family_probability(lo, family, arg1, arg2, arg3)
                if (too_low(value, target, increasing)) exit
                hi = lo
                if (lo <= 2.0_dp * tiny(1.0_dp)) then
                    root = 0.0_dp
                    return
                end if
                lo = 0.5_dp * lo
            end do
        end if

        root = bisect_positive(target, family, arg1, arg2, arg3, lo, hi, increasing)
    end function solve_positive_monotone

    pure function solve_nonnegative_monotone(target, family, arg1, arg2, arg3, increasing) result(root)
        real(dp), intent(in) :: target !! target probability strictly inside the attainable range
        integer, intent(in) :: family !! equation selector
        real(dp), intent(in) :: arg1 !! first fixed equation argument
        real(dp), intent(in) :: arg2 !! second fixed equation argument
        real(dp), intent(in) :: arg3 !! third fixed equation argument
        logical, intent(in) :: increasing !! true when the probability increases with the parameter
        real(dp) :: root

        integer :: iteration
        real(dp) :: hi
        real(dp) :: lo
        real(dp) :: value

        lo = 0.0_dp
        hi = 1.0_dp
        do iteration = 1, 1024
            value = family_probability(hi, family, arg1, arg2, arg3)
            if (.not. too_low(value, target, increasing)) exit
            lo = hi
            if (hi > huge(1.0_dp) / 2.0_dp) then
                root = positive_infinity(target)
                return
            end if
            hi = 2.0_dp * hi
        end do
        root = bisect_arithmetic(target, family, arg1, arg2, arg3, lo, hi, increasing)
    end function solve_nonnegative_monotone

    pure function solve_first_positive_root(target, family, arg1, arg2, arg3) result(root)
        real(dp), intent(in) :: target !! target probability strictly between zero and one
        integer, intent(in) :: family !! equation selector
        real(dp), intent(in) :: arg1 !! first fixed equation argument
        real(dp), intent(in) :: arg2 !! second fixed equation argument
        real(dp), intent(in) :: arg3 !! third fixed equation argument
        real(dp) :: root

        integer :: iteration
        real(dp), parameter :: start = 9.0949470177292823792e-13_dp ! 2**(-40)
        real(dp) :: current
        real(dp) :: f_current
        real(dp) :: f_previous
        real(dp) :: previous

        previous = start
        f_previous = family_probability(previous, family, arg1, arg2, arg3) - target
        if (abs(f_previous) <= probability_tolerance(target)) then
            root = previous
            return
        end if

        current = previous
        do iteration = 1, 160
            if (current > huge(1.0_dp) / 2.0_dp) exit
            current = 2.0_dp * current
            f_current = family_probability(current, family, arg1, arg2, arg3) - target
            if (abs(f_current) <= probability_tolerance(target)) then
                root = current
                return
            end if
            if (opposite_signs(f_previous, f_current)) then
                root = bisect_sign_change(target, family, arg1, arg2, arg3, previous, current)
                return
            end if
            previous = current
            f_previous = f_current
        end do
        root = quiet_nan(target)
    end function solve_first_positive_root

    pure function solve_real_nct_nc(df, target, t) result(root)
        real(dp), intent(in) :: df !! positive finite degrees of freedom
        real(dp), intent(in) :: target !! target probability strictly between zero and one
        real(dp), intent(in) :: t !! finite Student-t variate
        real(dp) :: root

        integer :: iteration
        real(dp) :: hi
        real(dp) :: lo
        real(dp) :: mid
        real(dp) :: value

        lo = -1.0_dp
        hi = 1.0_dp
        do iteration = 1, 1024
            if (nct_cdf(t, df, lo) >= target) exit
            if (lo < -huge(1.0_dp) / 2.0_dp) then
                root = negative_infinity(target)
                return
            end if
            lo = 2.0_dp * lo
        end do
        do iteration = 1, 1024
            if (nct_cdf(t, df, hi) <= target) exit
            if (hi > huge(1.0_dp) / 2.0_dp) then
                root = positive_infinity(target)
                return
            end if
            hi = 2.0_dp * hi
        end do

        do iteration = 1, 120
            mid = lo + 0.5_dp * (hi - lo)
            value = nct_cdf(t, df, mid)
            if (value > target) then
                lo = mid
            else
                hi = mid
            end if
            if (hi - lo <= 32.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(mid))) exit
        end do
        root = lo + 0.5_dp * (hi - lo)
    end function solve_real_nct_nc

    pure function bisect_positive(target, family, arg1, arg2, arg3, lo_in, hi_in, increasing) result(root)
        real(dp), intent(in) :: target !! target probability
        integer, intent(in) :: family !! equation selector
        real(dp), intent(in) :: arg1 !! first fixed equation argument
        real(dp), intent(in) :: arg2 !! second fixed equation argument
        real(dp), intent(in) :: arg3 !! third fixed equation argument
        real(dp), intent(in) :: lo_in !! positive lower bracket endpoint
        real(dp), intent(in) :: hi_in !! positive upper bracket endpoint
        logical, intent(in) :: increasing !! true when probability rises with parameter
        real(dp) :: root

        integer :: iteration
        real(dp) :: hi
        real(dp) :: lo
        real(dp) :: mid
        real(dp) :: value

        lo = lo_in
        hi = hi_in
        do iteration = 1, 120
            mid = sqrt(lo) * sqrt(hi)
            value = family_probability(mid, family, arg1, arg2, arg3)
            if (abs(value - target) <= probability_tolerance(target)) exit
            if (too_low(value, target, increasing)) then
                lo = mid
            else
                hi = mid
            end if
            if (hi <= lo * (1.0_dp + 32.0_dp * epsilon(1.0_dp))) exit
        end do
        root = mid
    end function bisect_positive

    pure function bisect_arithmetic(target, family, arg1, arg2, arg3, lo_in, hi_in, increasing) result(root)
        real(dp), intent(in) :: target !! target probability
        integer, intent(in) :: family !! equation selector
        real(dp), intent(in) :: arg1 !! first fixed equation argument
        real(dp), intent(in) :: arg2 !! second fixed equation argument
        real(dp), intent(in) :: arg3 !! third fixed equation argument
        real(dp), intent(in) :: lo_in !! nonnegative lower bracket endpoint
        real(dp), intent(in) :: hi_in !! positive upper bracket endpoint
        logical, intent(in) :: increasing !! true when probability rises with parameter
        real(dp) :: root

        integer :: iteration
        real(dp) :: hi
        real(dp) :: lo
        real(dp) :: mid
        real(dp) :: value

        lo = lo_in
        hi = hi_in
        do iteration = 1, 120
            mid = lo + 0.5_dp * (hi - lo)
            value = family_probability(mid, family, arg1, arg2, arg3)
            if (abs(value - target) <= probability_tolerance(target)) exit
            if (too_low(value, target, increasing)) then
                lo = mid
            else
                hi = mid
            end if
            if (hi - lo <= 32.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(mid))) exit
        end do
        root = mid
    end function bisect_arithmetic

    pure function bisect_sign_change(target, family, arg1, arg2, arg3, lo_in, hi_in) result(root)
        real(dp), intent(in) :: target !! target probability
        integer, intent(in) :: family !! equation selector
        real(dp), intent(in) :: arg1 !! first fixed equation argument
        real(dp), intent(in) :: arg2 !! second fixed equation argument
        real(dp), intent(in) :: arg3 !! third fixed equation argument
        real(dp), intent(in) :: lo_in !! positive lower bracket endpoint
        real(dp), intent(in) :: hi_in !! positive upper bracket endpoint
        real(dp) :: root

        integer :: iteration
        real(dp) :: f_lo
        real(dp) :: f_mid
        real(dp) :: hi
        real(dp) :: lo
        real(dp) :: mid

        lo = lo_in
        hi = hi_in
        f_lo = family_probability(lo, family, arg1, arg2, arg3) - target
        do iteration = 1, 120
            mid = sqrt(lo) * sqrt(hi)
            f_mid = family_probability(mid, family, arg1, arg2, arg3) - target
            if (abs(f_mid) <= probability_tolerance(target)) exit
            if (opposite_signs(f_lo, f_mid)) then
                hi = mid
            else
                lo = mid
                f_lo = f_mid
            end if
            if (hi <= lo * (1.0_dp + 32.0_dp * epsilon(1.0_dp))) exit
        end do
        root = mid
    end function bisect_sign_change

    pure function family_probability(parameter, family, arg1, arg2, arg3) result(value)
        real(dp), intent(in) :: parameter !! parameter being inverted
        integer, intent(in) :: family !! equation selector
        real(dp), intent(in) :: arg1 !! first fixed equation argument
        real(dp), intent(in) :: arg2 !! second fixed equation argument
        real(dp), intent(in) :: arg3 !! third fixed equation argument
        real(dp) :: value

        select case (family)
        case (family_chi_df)
            value = ncx2_cdf(arg1, parameter, arg2)
        case (family_chi_nc)
            value = ncx2_cdf(arg1, arg2, parameter)
        case (family_ncf_dfd)
            value = ncf_cdf(arg3, arg1, parameter, arg2)
        case (family_ncf_dfn)
            value = ncf_cdf(arg3, parameter, arg1, arg2)
        case (family_ncf_nc)
            value = ncf_cdf(arg3, arg1, arg2, parameter)
        case (family_nct_df)
            value = nct_cdf(arg2, parameter, arg1)
        case default
            value = quiet_nan(parameter)
        end select
    end function family_probability

    pure elemental logical function valid_probability(p) result(valid)
        real(dp), intent(in) :: p !! candidate probability

        valid = ieee_is_finite(p) .and. p >= 0.0_dp .and. p <= 1.0_dp
    end function valid_probability

    pure elemental logical function valid_open_probability(p) result(valid)
        real(dp), intent(in) :: p !! candidate probability

        valid = ieee_is_finite(p) .and. p > 0.0_dp .and. p < 1.0_dp
    end function valid_open_probability

    pure elemental logical function valid_positive_finite(x) result(valid)
        real(dp), intent(in) :: x !! positive finite parameter candidate

        valid = ieee_is_finite(x) .and. x > 0.0_dp
    end function valid_positive_finite

    pure elemental logical function valid_nonnegative_finite(x) result(valid)
        real(dp), intent(in) :: x !! nonnegative finite parameter candidate

        valid = ieee_is_finite(x) .and. x >= 0.0_dp
    end function valid_nonnegative_finite

    pure elemental logical function too_low(value, target, increasing) result(answer)
        real(dp), intent(in) :: value !! probability at current parameter
        real(dp), intent(in) :: target !! target probability
        logical, intent(in) :: increasing !! true when probability rises with parameter

        if (increasing) then
            answer = value < target
        else
            answer = value > target
        end if
    end function too_low

    pure elemental logical function opposite_signs(a, b) result(answer)
        real(dp), intent(in) :: a !! first finite residual
        real(dp), intent(in) :: b !! second finite residual

        answer = (a < 0.0_dp .and. b > 0.0_dp) .or. (a > 0.0_dp .and. b < 0.0_dp)
    end function opposite_signs

    pure elemental function probability_tolerance(target) result(tol)
        real(dp), intent(in) :: target !! probability target
        real(dp) :: tol

        tol = 32.0_dp * epsilon(1.0_dp) * max(tiny(1.0_dp), abs(target))
    end function probability_tolerance

end module scifort_special_noncentral_distribution_helpers
