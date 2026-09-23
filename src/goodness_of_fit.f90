! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2001-2002 Enthought, Inc. 2003-2026, SciPy Developers.
! Copyright (c) 2026 SciFort contributors
!
! Goodness-of-fit and normality tests matching scipy.stats 1.17 semantics.
! The Shapiro-Wilk implementation is adapted from SciPy's double-precision
! translation of Royston's AS R94 / Algorithm AS 181 implementation.
! The Cramer-von Mises finite-sample correction follows scipy.stats._hypotests.
! See THIRD_PARTY_LICENSES.md and CODE_PROVENANCE.md.

module scifort_goodness_of_fit
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_bessel_k, only : log_besselk
    use scifort_chi2, only : chi2_sf
    use scifort_constants, only : scifort_pi
    use scifort_descriptive, only : central_moment, mean, rankdata, standard_deviation
    use scifort_exponential, only : exponential_logcdf, exponential_logsf
    use scifort_gumbel_l, only : gumbel_l_logcdf, gumbel_l_logsf
    use scifort_gumbel_r, only : gumbel_r_logcdf, gumbel_r_logsf
    use scifort_kinds, only : dp
    use scifort_ksone, only : ksone_sf
    use scifort_kstwo, only : kstwo_sf
    use scifort_kstwobign, only : kstwobign_sf
    use scifort_logistic, only : logistic_logcdf, logistic_logsf
    use scifort_math, only : quiet_nan
    use scifort_normal, only : normal_cdf, normal_logcdf, normal_logsf, normal_sf
    implicit none
    private

    abstract interface
        function cdf_callback(x) result(p)
            import :: dp
            real(dp), intent(in) :: x !! point at which the reference CDF is evaluated
            real(dp) :: p
        end function cdf_callback
    end interface

    type, public :: significance_result
        real(dp) :: statistic !! test statistic
        real(dp) :: pvalue !! p-value
    end type significance_result

    type, public :: kstest_result
        real(dp) :: statistic !! Kolmogorov-Smirnov distance
        real(dp) :: pvalue !! exact or asymptotic p-value
        real(dp) :: statistic_location !! observation at which the extremum occurs
        integer :: statistic_sign !! +1 for D+, -1 for D-
    end type kstest_result

    type, public :: shapiro_result
        real(dp) :: statistic !! Shapiro-Wilk W statistic
        real(dp) :: pvalue !! Royston approximation to the null p-value
        integer :: status !! 0 success, 2 n > 5000 warning, 6 zero-range data
    end type shapiro_result

    type, public :: anderson_result
        real(dp) :: statistic !! Anderson-Darling A2 statistic
        real(dp) :: pvalue !! interpolated p-value when requested, NaN otherwise
        real(dp), allocatable :: critical_values(:) !! legacy SciPy critical values
        real(dp), allocatable :: significance_level(:) !! corresponding percentages
        real(dp), allocatable :: fit_params(:) !! fitted distribution parameters
        integer :: status !! 0 success; nonzero indicates unsupported/invalid input
    end type anderson_result

    interface kstest
        module procedure kstest_one_sample
        module procedure kstest_two_sample
    end interface kstest

    public :: anderson
    public :: cramervonmises
    public :: cramervonmises_2samp
    public :: jarque_bera
    public :: ks_1samp
    public :: ks_2samp
    public :: kstest
    public :: kurtosistest
    public :: normaltest
    public :: shapiro
    public :: skewtest

contains

    function skewtest(a, alternative) result(res)
        real(dp), intent(in) :: a(:) !! sample, at least eight finite observations
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        type(significance_result) :: res

        character(len=9) :: alt
        integer :: n
        real(dp) :: alpha
        real(dp) :: b2
        real(dp) :: beta2
        real(dp) :: delta
        real(dp) :: w2
        real(dp) :: y
        real(dp) :: z

        alt = select_alternative(alternative)
        n = size(a)
        if (n < 8 .or. .not. finite_sample(a) .or. .not. valid_alternative(alt)) then
            call invalidate_significance(res)
            return
        end if
        b2 = sample_skew(a)
        if (.not. ieee_is_finite(b2)) then
            call invalidate_significance(res)
            return
        end if
        y = b2 * sqrt(real((n + 1) * (n + 3), dp) / (6.0_dp * real(n - 2, dp)))
        beta2 = 3.0_dp * real(n * n + 27 * n - 70, dp) * real((n + 1) * (n + 3), dp) / &
            (real(n - 2, dp) * real(n + 5, dp) * real(n + 7, dp) * real(n + 9, dp))
        w2 = -1.0_dp + sqrt(2.0_dp * (beta2 - 1.0_dp))
        delta = 1.0_dp / sqrt(0.5_dp * log(w2))
        alpha = sqrt(2.0_dp / (w2 - 1.0_dp))
        if (y == 0.0_dp) y = 1.0_dp
        z = delta * log(y / alpha + sqrt((y / alpha) ** 2 + 1.0_dp))
        res%statistic = z
        res%pvalue = normal_pvalue(z, alt)
    end function skewtest

    function kurtosistest(a, alternative) result(res)
        real(dp), intent(in) :: a(:) !! sample, at least five finite observations
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        type(significance_result) :: res

        character(len=9) :: alt
        integer :: n
        real(dp) :: aa
        real(dp) :: b2
        real(dp) :: denom
        real(dp) :: e
        real(dp) :: sqrtbeta1
        real(dp) :: term1
        real(dp) :: term2
        real(dp) :: varb2
        real(dp) :: x
        real(dp) :: z

        alt = select_alternative(alternative)
        n = size(a)
        if (n < 5 .or. .not. finite_sample(a) .or. .not. valid_alternative(alt)) then
            call invalidate_significance(res)
            return
        end if
        b2 = sample_kurtosis_pearson(a)
        if (.not. ieee_is_finite(b2)) then
            call invalidate_significance(res)
            return
        end if
        e = 3.0_dp * real(n - 1, dp) / real(n + 1, dp)
        varb2 = 24.0_dp * real(n * (n - 2) * (n - 3), dp) / &
            (real((n + 1) * (n + 1), dp) * real(n + 3, dp) * real(n + 5, dp))
        x = (b2 - e) / sqrt(varb2)
        sqrtbeta1 = 6.0_dp * real(n * n - 5 * n + 2, dp) / &
            real((n + 7) * (n + 9), dp) * &
            sqrt(6.0_dp * real((n + 3) * (n + 5), dp) / &
            real(n * (n - 2) * (n - 3), dp))
        aa = 6.0_dp + 8.0_dp / sqrtbeta1 * &
            (2.0_dp / sqrtbeta1 + sqrt(1.0_dp + 4.0_dp / (sqrtbeta1 * sqrtbeta1)))
        term1 = 1.0_dp - 2.0_dp / (9.0_dp * aa)
        denom = 1.0_dp + x * sqrt(2.0_dp / (aa - 4.0_dp))
        if (denom == 0.0_dp) then
            call invalidate_significance(res)
            return
        end if
        term2 = sign(1.0_dp, denom) * ((1.0_dp - 2.0_dp / aa) / abs(denom)) ** (1.0_dp / 3.0_dp)
        z = (term1 - term2) / sqrt(2.0_dp / (9.0_dp * aa))
        res%statistic = z
        res%pvalue = normal_pvalue(z, alt)
    end function kurtosistest

    function normaltest(a) result(res)
        real(dp), intent(in) :: a(:) !! sample, at least eight finite observations
        type(significance_result) :: res

        type(significance_result) :: kr
        type(significance_result) :: sk

        sk = skewtest(a)
        kr = kurtosistest(a)
        if (ieee_is_nan(sk%statistic) .or. ieee_is_nan(kr%statistic)) then
            call invalidate_significance(res)
            return
        end if
        res%statistic = sk%statistic ** 2 + kr%statistic ** 2
        res%pvalue = chi2_sf(res%statistic, 2.0_dp)
    end function normaltest

    function jarque_bera(x) result(res)
        real(dp), intent(in) :: x(:) !! finite sample observations
        type(significance_result) :: res

        real(dp) :: k
        real(dp) :: s

        if (size(x) < 2 .or. .not. finite_sample(x)) then
            call invalidate_significance(res)
            return
        end if
        s = sample_skew(x)
        k = sample_kurtosis_pearson(x) - 3.0_dp
        if (.not. ieee_is_finite(s) .or. .not. ieee_is_finite(k)) then
            call invalidate_significance(res)
            return
        end if
        res%statistic = real(size(x), dp) / 6.0_dp * (s * s + 0.25_dp * k * k)
        res%pvalue = chi2_sf(res%statistic, 2.0_dp)
    end function jarque_bera

    function ks_1samp(x, cdf, alternative, method) result(res)
        real(dp), intent(in) :: x(:) !! observations from a continuous distribution
        procedure(cdf_callback) :: cdf !! reference cumulative distribution function
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        character(len=*), intent(in), optional :: method !! auto, exact, approx, or asymp
        type(kstest_result) :: res

        character(len=9) :: alt
        character(len=8) :: selected_method
        integer :: i
        integer :: iminus
        integer :: iplus
        integer :: n
        real(dp) :: d
        real(dp) :: dminus
        real(dp) :: dplus
        real(dp), allocatable :: sx(:)
        real(dp) :: f
        real(dp) :: minus_value
        real(dp) :: plus_value

        alt = select_alternative(alternative)
        selected_method = select_method(method, 'exact')
        n = size(x)
        if (n < 1 .or. .not. finite_sample(x) .or. .not. valid_alternative(alt)) then
            call invalidate_ks(res)
            return
        end if
        if (.not. valid_ks1_method(selected_method)) then
            call invalidate_ks(res)
            return
        end if
        sx = x
        call sort_real(sx)
        dplus = -1.0_dp
        dminus = -1.0_dp
        iplus = 1
        iminus = 1
        do i = 1, n
            f = cdf(sx(i))
            if (.not. ieee_is_finite(f)) then
                call invalidate_ks(res)
                return
            end if
            f = min(1.0_dp, max(0.0_dp, f))
            plus_value = real(i, dp) / real(n, dp) - f
            minus_value = f - real(i - 1, dp) / real(n, dp)
            if (plus_value > dplus) then
                dplus = plus_value
                iplus = i
            end if
            if (minus_value > dminus) then
                dminus = minus_value
                iminus = i
            end if
        end do

        if (alt == 'greater') then
            res%statistic = max(0.0_dp, dplus)
            res%statistic_location = sx(iplus)
            res%statistic_sign = 1
            res%pvalue = ksone_sf(res%statistic, real(n, dp))
            return
        else if (alt == 'less') then
            res%statistic = max(0.0_dp, dminus)
            res%statistic_location = sx(iminus)
            res%statistic_sign = -1
            res%pvalue = ksone_sf(res%statistic, real(n, dp))
            return
        end if

        if (dplus > dminus) then
            d = dplus
            res%statistic_location = sx(iplus)
            res%statistic_sign = 1
        else
            d = dminus
            res%statistic_location = sx(iminus)
            res%statistic_sign = -1
        end if
        res%statistic = max(0.0_dp, d)
        if (selected_method == 'auto' .or. selected_method == 'exact') then
            res%pvalue = kstwo_sf(res%statistic, real(n, dp))
        else if (selected_method == 'asymp') then
            res%pvalue = kstwobign_sf(res%statistic * sqrt(real(n, dp)))
        else
            res%pvalue = min(1.0_dp, 2.0_dp * ksone_sf(res%statistic, real(n, dp)))
        end if
        res%pvalue = min(1.0_dp, max(0.0_dp, res%pvalue))
    end function ks_1samp

    function ks_2samp(data1, data2, alternative, method) result(res)
        real(dp), intent(in) :: data1(:) !! first independent sample
        real(dp), intent(in) :: data2(:) !! second independent sample
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        character(len=*), intent(in), optional :: method !! auto, exact, or asymp
        type(kstest_result) :: res

        character(len=9) :: alt
        character(len=8) :: selected_method
        integer :: i
        integer :: j
        integer :: n1
        integer :: n2
        real(dp) :: d
        real(dp) :: diff
        real(dp) :: maxs
        real(dp) :: mins
        real(dp) :: loc_max
        real(dp) :: loc_min
        real(dp) :: m
        real(dp) :: n
        real(dp) :: en
        real(dp) :: z
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: y(:)

        alt = select_alternative(alternative)
        selected_method = select_method(method, 'auto')
        n1 = size(data1)
        n2 = size(data2)
        if (n1 < 1 .or. n2 < 1 .or. .not. finite_sample(data1) .or. &
                .not. finite_sample(data2) .or. .not. valid_alternative(alt)) then
            call invalidate_ks(res)
            return
        end if
        if (.not. (selected_method == 'auto' .or. selected_method == 'exact' .or. &
                selected_method == 'asymp')) then
            call invalidate_ks(res)
            return
        end if
        x = data1
        y = data2
        call sort_real(x)
        call sort_real(y)
        i = 0
        j = 0
        maxs = -1.0_dp
        mins = -1.0_dp
        loc_max = min(x(1), y(1))
        loc_min = loc_max
        do while (i < n1 .or. j < n2)
            if (j >= n2) then
                d = x(i + 1)
            else if (i >= n1) then
                d = y(j + 1)
            else if (x(i + 1) <= y(j + 1)) then
                d = x(i + 1)
            else
                d = y(j + 1)
            end if
            do while (i < n1)
                if (x(i + 1) > d) exit
                i = i + 1
            end do
            do while (j < n2)
                if (y(j + 1) > d) exit
                j = j + 1
            end do
            diff = real(i, dp) / real(n1, dp) - real(j, dp) / real(n2, dp)
            if (diff > maxs) then
                maxs = diff
                loc_max = d
            end if
            if (-diff > mins) then
                mins = -diff
                loc_min = d
            end if
        end do
        maxs = max(0.0_dp, maxs)
        mins = max(0.0_dp, mins)
        if (alt == 'less' .or. (alt == 'two-sided' .and. mins > maxs)) then
            res%statistic = mins
            res%statistic_location = loc_min
            res%statistic_sign = -1
        else
            res%statistic = maxs
            res%statistic_location = loc_max
            res%statistic_sign = 1
        end if

        if (selected_method == 'auto') then
            if (max(n1, n2) <= 10000 .and. int(n1, kind=8) * int(n2, kind=8) <= 25000000_8) then
                selected_method = 'exact'
            else
                selected_method = 'asymp'
            end if
        end if
        if (selected_method == 'exact' .and. &
                int(n1, kind=8) * int(n2, kind=8) <= 25000000_8) then
            res%pvalue = ks2_exact_probability(n1, n2, res%statistic, alt)
        else
            m = real(max(n1, n2), dp)
            n = real(min(n1, n2), dp)
            en = m * n / (m + n)
            if (alt == 'two-sided') then
                res%pvalue = kstwo_sf(res%statistic, real(nint(en), dp))
            else
                z = sqrt(en) * res%statistic
                res%pvalue = exp(-2.0_dp * z * z - &
                    2.0_dp * z * (m + 2.0_dp * n) / sqrt(m * n * (m + n)) / 3.0_dp)
            end if
        end if
        res%pvalue = min(1.0_dp, max(0.0_dp, res%pvalue))
    end function ks_2samp

    function kstest_one_sample(x, cdf, alternative, method) result(res)
        real(dp), intent(in) :: x(:) !! one sample
        procedure(cdf_callback) :: cdf !! reference cumulative distribution function
        character(len=*), intent(in), optional :: alternative !! alternative hypothesis
        character(len=*), intent(in), optional :: method !! p-value method
        type(kstest_result) :: res

        res = ks_1samp(x, cdf, alternative, method)
    end function kstest_one_sample

    function kstest_two_sample(x, y, alternative, method) result(res)
        real(dp), intent(in) :: x(:) !! first sample
        real(dp), intent(in) :: y(:) !! second sample
        character(len=*), intent(in), optional :: alternative !! alternative hypothesis
        character(len=*), intent(in), optional :: method !! p-value method
        type(kstest_result) :: res

        res = ks_2samp(x, y, alternative, method)
    end function kstest_two_sample

    function cramervonmises(rvs, cdf) result(res)
        real(dp), intent(in) :: rvs(:) !! observations, at least two
        procedure(cdf_callback) :: cdf !! reference cumulative distribution function
        type(significance_result) :: res

        integer :: i
        integer :: n
        real(dp) :: f
        real(dp) :: u
        real(dp) :: w
        real(dp), allocatable :: x(:)

        n = size(rvs)
        if (n <= 1 .or. .not. finite_sample(rvs)) then
            call invalidate_significance(res)
            return
        end if
        x = rvs
        call sort_real(x)
        w = 1.0_dp / (12.0_dp * real(n, dp))
        do i = 1, n
            f = cdf(x(i))
            if (.not. ieee_is_finite(f)) then
                call invalidate_significance(res)
                return
            end if
            u = real(2 * i - 1, dp) / (2.0_dp * real(n, dp))
            w = w + (u - f) ** 2
        end do
        res%statistic = w
        res%pvalue = max(0.0_dp, min(1.0_dp, 1.0_dp - cvm_cdf(w, n)))
    end function cramervonmises

    function cramervonmises_2samp(x, y, method) result(res)
        real(dp), intent(in) :: x(:) !! first sample, at least two observations
        real(dp), intent(in) :: y(:) !! second sample, at least two observations
        character(len=*), intent(in), optional :: method !! auto, exact, or asymptotic
        type(significance_result) :: res

        character(len=10) :: selected_method
        integer :: i
        integer :: k
        integer :: n_total
        integer :: nx
        integer :: ny
        real(dp) :: et
        real(dp) :: t
        real(dp) :: tn
        real(dp) :: u
        real(dp) :: vt
        real(dp), allocatable :: pooled(:)
        real(dp), allocatable :: ranks(:)

        nx = size(x)
        ny = size(y)
        if (nx <= 1 .or. ny <= 1 .or. .not. finite_sample(x) .or. .not. finite_sample(y)) then
            call invalidate_significance(res)
            return
        end if
        selected_method = 'auto'
        if (present(method)) selected_method = lower_string(adjustl(method))
        if (.not. (trim(selected_method) == 'auto' .or. trim(selected_method) == 'exact' .or. &
                trim(selected_method) == 'asymptotic' .or. trim(selected_method) == 'asymp')) then
            call invalidate_significance(res)
            return
        end if
        n_total = nx + ny
        allocate(pooled(n_total))
        pooled(1:nx) = x
        pooled(nx + 1:) = y
        ranks = rankdata(pooled)
        u = 0.0_dp
        do i = 1, nx
            u = u + real(nx, dp) * (ranks(i) - real(i, dp)) ** 2
        end do
        do i = 1, ny
            u = u + real(ny, dp) * (ranks(nx + i) - real(i, dp)) ** 2
        end do
        k = nx * ny
        t = u / real(k * n_total, dp) - real(4 * k - 1, dp) / real(6 * n_total, dp)
        res%statistic = t

        if (trim(selected_method) == 'auto') then
            if (max(nx, ny) <= 20 .and. combination_leq(nx + ny, nx, 2000000_8)) then
                selected_method = 'exact'
            else
                selected_method = 'asymptotic'
            end if
        end if
        if (trim(selected_method) == 'exact' .and. combination_leq(nx + ny, nx, 2000000_8)) then
            res%pvalue = cvm2_exact_enumeration(u, nx, ny)
        else
            et = (1.0_dp + 1.0_dp / real(n_total, dp)) / 6.0_dp
            vt = real(n_total + 1, dp) * &
                real(4 * k * n_total - 3 * (nx * nx + ny * ny) - 2 * k, dp)
            vt = vt / (45.0_dp * real(n_total * n_total, dp) * 4.0_dp * real(k, dp))
            tn = 1.0_dp / 6.0_dp + (t - et) / sqrt(45.0_dp * vt)
            if (tn < 0.003_dp) then
                res%pvalue = 1.0_dp
            else
                res%pvalue = max(0.0_dp, min(1.0_dp, 1.0_dp - cvm_cdf_inf(tn)))
            end if
        end if
    end function cramervonmises_2samp

    function shapiro(x) result(res)
        real(dp), intent(in) :: x(:) !! sample, at least three finite observations
        type(shapiro_result) :: res

        integer :: n
        real(dp), allocatable :: a(:)
        real(dp), allocatable :: y(:)

        n = size(x)
        if (n < 3 .or. .not. finite_sample(x)) then
            res%statistic = quiet_nan(0.0_dp)
            res%pvalue = quiet_nan(0.0_dp)
            res%status = 1
            return
        end if
        y = x
        call sort_real(y)
        y = y - x(n / 2 + 1)
        allocate(a(n / 2))
        call swilk(y, a, res%statistic, res%pvalue, res%status)
    end function shapiro

    function anderson(x, dist, interpolate_pvalue) result(res)
        real(dp), intent(in) :: x(:) !! finite sample observations
        character(len=*), intent(in), optional :: dist !! norm, expon, logistic, gumbel_l, or gumbel_r
        logical, intent(in), optional :: interpolate_pvalue !! compute table-interpolated p-value
        type(anderson_result) :: res

        character(len=16) :: selected_dist
        integer :: i
        integer :: n
        logical :: want_pvalue
        real(dp) :: loc
        real(dp) :: scale
        real(dp), allocatable :: logcdf(:)
        real(dp), allocatable :: logsf(:)
        real(dp), allocatable :: y(:)

        selected_dist = 'norm'
        if (present(dist)) selected_dist = lower_string(adjustl(dist))
        if (trim(selected_dist) == 'gumbel' .or. trim(selected_dist) == 'extreme1') then
            selected_dist = 'gumbel_l'
        end if
        want_pvalue = .false.
        if (present(interpolate_pvalue)) want_pvalue = interpolate_pvalue
        n = size(x)
        if (n < 2 .or. .not. finite_sample(x)) then
            call invalidate_anderson(res)
            return
        end if
        y = x
        call sort_real(y)
        allocate(logcdf(n), logsf(n))

        select case (trim(selected_dist))
        case ('norm')
            loc = mean(x)
            scale = standard_deviation(x, 1)
            if (.not. ieee_is_finite(scale) .or. scale <= 0.0_dp) then
                call invalidate_anderson(res)
                return
            end if
            allocate(res%fit_params(2), res%critical_values(5), res%significance_level(5))
            res%fit_params = [loc, scale]
            res%significance_level = [15.0_dp, 10.0_dp, 5.0_dp, 2.5_dp, 1.0_dp]
            res%critical_values = nint(1000.0_dp * &
                [0.561_dp, 0.631_dp, 0.752_dp, 0.873_dp, 1.035_dp] / &
                (1.0_dp + 0.75_dp / real(n, dp) + 2.25_dp / real(n * n, dp))) / 1000.0_dp
            do i = 1, n
                logcdf(i) = normal_logcdf((y(i) - loc) / scale)
                logsf(i) = normal_logsf((y(i) - loc) / scale)
            end do
        case ('expon')
            loc = 0.0_dp
            scale = mean(x)
            if (.not. ieee_is_finite(scale) .or. scale <= 0.0_dp) then
                call invalidate_anderson(res)
                return
            end if
            allocate(res%fit_params(2), res%critical_values(5), res%significance_level(5))
            res%fit_params = [loc, scale]
            res%significance_level = [15.0_dp, 10.0_dp, 5.0_dp, 2.5_dp, 1.0_dp]
            res%critical_values = nint(1000.0_dp * &
                [0.916_dp, 1.062_dp, 1.321_dp, 1.591_dp, 1.959_dp] / &
                (1.0_dp + 0.6_dp / real(n, dp))) / 1000.0_dp
            do i = 1, n
                logcdf(i) = exponential_logcdf(y(i), 0.0_dp, scale)
                logsf(i) = exponential_logsf(y(i), 0.0_dp, scale)
            end do
        case ('logistic')
            call logistic_mle(x, loc, scale)
            if (.not. ieee_is_finite(scale) .or. scale <= 0.0_dp) then
                call invalidate_anderson(res)
                return
            end if
            allocate(res%fit_params(2), res%critical_values(6), res%significance_level(6))
            res%fit_params = [loc, scale]
            res%significance_level = [25.0_dp, 10.0_dp, 5.0_dp, 2.5_dp, 1.0_dp, 0.5_dp]
            res%critical_values = nint(1000.0_dp * &
                [0.426_dp, 0.563_dp, 0.660_dp, 0.769_dp, 0.906_dp, 1.010_dp] / &
                (1.0_dp + 0.25_dp / real(n, dp))) / 1000.0_dp
            do i = 1, n
                logcdf(i) = logistic_logcdf(y(i), loc, scale)
                logsf(i) = logistic_logsf(y(i), loc, scale)
            end do
        case ('gumbel_r')
            call gumbel_r_mle(x, loc, scale)
            if (.not. ieee_is_finite(scale) .or. scale <= 0.0_dp) then
                call invalidate_anderson(res)
                return
            end if
            allocate(res%fit_params(2), res%critical_values(5), res%significance_level(5))
            res%fit_params = [loc, scale]
            res%significance_level = [25.0_dp, 10.0_dp, 5.0_dp, 2.5_dp, 1.0_dp]
            res%critical_values = nint(1000.0_dp * &
                [0.474_dp, 0.637_dp, 0.757_dp, 0.877_dp, 1.038_dp] / &
                (1.0_dp + 0.2_dp / sqrt(real(n, dp)))) / 1000.0_dp
            do i = 1, n
                logcdf(i) = gumbel_r_logcdf(y(i), loc, scale)
                logsf(i) = gumbel_r_logsf(y(i), loc, scale)
            end do
        case ('gumbel_l')
            call gumbel_l_mle(x, loc, scale)
            if (.not. ieee_is_finite(scale) .or. scale <= 0.0_dp) then
                call invalidate_anderson(res)
                return
            end if
            allocate(res%fit_params(2), res%critical_values(5), res%significance_level(5))
            res%fit_params = [loc, scale]
            res%significance_level = [25.0_dp, 10.0_dp, 5.0_dp, 2.5_dp, 1.0_dp]
            res%critical_values = nint(1000.0_dp * &
                [0.474_dp, 0.637_dp, 0.757_dp, 0.877_dp, 1.038_dp] / &
                (1.0_dp + 0.2_dp / sqrt(real(n, dp)))) / 1000.0_dp
            do i = 1, n
                logcdf(i) = gumbel_l_logcdf(y(i), loc, scale)
                logsf(i) = gumbel_l_logsf(y(i), loc, scale)
            end do
        case default
            call invalidate_anderson(res)
            res%status = 2
            return
        end select

        res%statistic = -real(n, dp)
        do i = 1, n
            res%statistic = res%statistic - real(2 * i - 1, dp) / real(n, dp) * &
                (logcdf(i) + logsf(n + 1 - i))
        end do
        res%status = 0
        res%pvalue = quiet_nan(0.0_dp)
        if (want_pvalue) then
            res%pvalue = interpolate_table(res%statistic, res%critical_values, &
                res%significance_level / 100.0_dp)
        end if
    end function anderson

    pure function sample_skew(x) result(value)
        real(dp), intent(in) :: x(:) !! finite sample
        real(dp) :: value
        real(dp) :: m2
        real(dp) :: m3

        m2 = central_moment(x, 2)
        m3 = central_moment(x, 3)
        if (m2 <= 0.0_dp .or. .not. ieee_is_finite(m2) .or. .not. ieee_is_finite(m3)) then
            value = quiet_nan(0.0_dp)
        else
            value = m3 / (m2 ** 1.5_dp)
        end if
    end function sample_skew

    pure function sample_kurtosis_pearson(x) result(value)
        real(dp), intent(in) :: x(:) !! finite sample
        real(dp) :: value
        real(dp) :: m2
        real(dp) :: m4

        m2 = central_moment(x, 2)
        m4 = central_moment(x, 4)
        if (m2 <= 0.0_dp .or. .not. ieee_is_finite(m2) .or. .not. ieee_is_finite(m4)) then
            value = quiet_nan(0.0_dp)
        else
            value = m4 / (m2 * m2)
        end if
    end function sample_kurtosis_pearson

    pure function normal_pvalue(z, alternative) result(p)
        real(dp), intent(in) :: z !! standard-normal test statistic
        character(len=*), intent(in) :: alternative !! selected alternative
        real(dp) :: p

        select case (trim(alternative))
        case ('less')
            p = normal_cdf(z)
        case ('greater')
            p = normal_sf(z)
        case default
            p = min(1.0_dp, 2.0_dp * normal_sf(abs(z)))
        end select
    end function normal_pvalue

    function ks2_exact_probability(n1, n2, d_observed, alternative) result(pout)
        integer, intent(in) :: n1 !! first sample size
        integer, intent(in) :: n2 !! second sample size
        real(dp), intent(in) :: d_observed !! observed ECDF distance
        character(len=*), intent(in) :: alternative !! test alternative
        real(dp) :: pout

        integer :: g
        integer :: h
        integer :: i
        integer :: j
        integer(kind=8) :: lcm
        real(dp) :: d
        real(dp) :: denom
        real(dp) :: diff
        real(dp), allocatable :: current(:)
        real(dp), allocatable :: previous(:)

        g = gcd_int(n1, n2)
        lcm = int(n1 / g, kind=8) * int(n2, kind=8)
        h = nint(d_observed * real(lcm, dp))
        if (h <= 0) then
            pout = 1.0_dp
            return
        end if
        d = real(h, dp) / real(lcm, dp)
        allocate(previous(0:n2), current(0:n2))
        previous = 0.0_dp
        previous(0) = 1.0_dp
        do i = 0, n1
            current = 0.0_dp
            do j = 0, n2
                if (i == 0 .and. j == 0) then
                    current(j) = 1.0_dp
                    cycle
                end if
                diff = real(i, dp) / real(n1, dp) - real(j, dp) / real(n2, dp)
                if (.not. ks_inside(diff, d, alternative)) cycle
                denom = real(n1 + n2 - (i + j) + 1, dp)
                if (i > 0) then
                    current(j) = current(j) + previous(j) * &
                        real(n1 - i + 1, dp) / denom
                end if
                if (j > 0) then
                    current(j) = current(j) + current(j - 1) * &
                        real(n2 - j + 1, dp) / denom
                end if
            end do
            previous = current
        end do
        pout = min(1.0_dp, max(0.0_dp, 1.0_dp - current(n2)))
    end function ks2_exact_probability

    pure function ks_inside(diff, d, alternative) result(inside)
        real(dp), intent(in) :: diff !! F1-F2 lattice difference
        real(dp), intent(in) :: d !! observed lattice distance
        character(len=*), intent(in) :: alternative !! alternative hypothesis
        logical :: inside
        real(dp), parameter :: tol = 32.0_dp * epsilon(1.0_dp)

        select case (trim(alternative))
        case ('greater')
            inside = diff < d - tol
        case ('less')
            inside = -diff < d - tol
        case default
            inside = abs(diff) < d - tol
        end select
    end function ks_inside

    function cvm_cdf(x, n) result(value)
        real(dp), intent(in) :: x !! Cramer-von Mises statistic
        integer, intent(in) :: n !! sample size
        real(dp) :: value
        real(dp) :: v

        if (x <= 1.0_dp / (12.0_dp * real(n, dp))) then
            value = 0.0_dp
        else if (x >= real(n, dp) / 3.0_dp) then
            value = 1.0_dp
        else
            v = cvm_cdf_inf(x)
            value = v * (1.0_dp + 1.0_dp / (12.0_dp * real(n, dp))) + &
                cvm_psi1_mod(x) / real(n, dp)
        end if
    end function cvm_cdf

    function cvm_cdf_inf(x) result(value)
        real(dp), intent(in) :: x !! positive asymptotic statistic
        real(dp) :: value

        integer :: k
        real(dp) :: logterm
        real(dp) :: q
        real(dp) :: term
        real(dp) :: y

        if (x <= 0.0_dp) then
            value = 0.0_dp
            return
        end if
        value = 0.0_dp
        do k = 0, 1000
            y = real(4 * k + 1, dp)
            q = y * y / (16.0_dp * x)
            logterm = log_gamma(real(k, dp) + 0.5_dp) - log_gamma(real(k + 1, dp)) - &
                1.5_dp * log(scifort_pi) - 0.5_dp * log(x) + 0.5_dp * log(y) - q + &
                log_besselk(0.25_dp, q)
            term = exp(logterm)
            value = value + term
            if (abs(term) < 1.0e-7_dp) exit
        end do
    end function cvm_cdf_inf

    function cvm_psi1_mod(x) result(value)
        real(dp), intent(in) :: x !! positive asymptotic statistic
        real(dp) :: value

        integer :: k
        real(dp) :: ak
        real(dp) :: term

        if (x <= 0.0_dp) then
            value = 0.0_dp
            return
        end if
        value = 0.0_dp
        do k = 0, 1000
            ak = cvm_ak(k, x)
            term = -ak / (scifort_pi * gamma(real(k + 1, dp)))
            value = value + term
            if (abs(term) < 1.0e-7_dp) exit
        end do
    end function cvm_psi1_mod

    function cvm_ak(k, x) result(value)
        integer, intent(in) :: k !! nonnegative series index
        real(dp), intent(in) :: x !! positive statistic
        real(dp) :: value

        real(dp) :: e1
        real(dp) :: e2
        real(dp) :: e3
        real(dp) :: e4
        real(dp) :: e5
        real(dp) :: g1
        real(dp) :: g3
        real(dp) :: m
        real(dp) :: sx

        m = real(2 * k + 1, dp)
        sx = 2.0_dp * sqrt(x)
        g1 = gamma(real(k, dp) + 0.5_dp)
        g3 = gamma(real(k, dp) + 1.5_dp)
        e1 = m * g1 * cvm_ed2(real(4 * k + 3, dp) / sx) / (9.0_dp * x ** 0.75_dp)
        e2 = g1 * cvm_ed3(real(4 * k + 1, dp) / sx) / (72.0_dp * x ** 1.25_dp)
        e3 = 2.0_dp * (m + 2.0_dp) * g3 * &
            cvm_ed3(real(4 * k + 5, dp) / sx) / (12.0_dp * x ** 1.25_dp)
        e4 = 7.0_dp * m * g1 * cvm_ed2(real(4 * k + 1, dp) / sx) / (144.0_dp * x ** 0.75_dp)
        e5 = 7.0_dp * m * g1 * cvm_ed2(real(4 * k + 5, dp) / sx) / (144.0_dp * x ** 0.75_dp)
        value = e1 + e2 + e3 + e4 + e5
    end function cvm_ak

    function cvm_ed2(y) result(value)
        real(dp), intent(in) :: y !! positive auxiliary argument
        real(dp) :: value
        real(dp) :: z
        real(dp) :: l1
        real(dp) :: l2
        real(dp) :: mlog

        z = y * y / 4.0_dp
        l1 = log_besselk(0.25_dp, z)
        l2 = log_besselk(0.75_dp, z)
        mlog = max(l1, l2)
        value = exp(-z + 1.5_dp * log(y / 2.0_dp) - 0.5_dp * log(scifort_pi) + mlog) * &
            (exp(l1 - mlog) + exp(l2 - mlog))
    end function cvm_ed2

    function cvm_ed3(y) result(value)
        real(dp), intent(in) :: y !! positive auxiliary argument
        real(dp) :: value
        real(dp) :: z
        real(dp) :: k1
        real(dp) :: k3
        real(dp) :: k5
        real(dp) :: mlog

        z = y * y / 4.0_dp
        mlog = max(log_besselk(0.25_dp, z), &
            max(log_besselk(0.75_dp, z), log_besselk(1.25_dp, z)))
        k1 = exp(log_besselk(0.25_dp, z) - mlog)
        k3 = exp(log_besselk(0.75_dp, z) - mlog)
        k5 = exp(log_besselk(1.25_dp, z) - mlog)
        value = exp(-z + 2.5_dp * log(y / 2.0_dp) - 0.5_dp * log(scifort_pi) + mlog) * &
            (2.0_dp * k1 + 3.0_dp * k3 - k5)
    end function cvm_ed3

    function cvm2_exact_enumeration(observed_u, nx, ny) result(pvalue)
        real(dp), intent(in) :: observed_u !! observed rank statistic U
        integer, intent(in) :: nx !! first sample size
        integer, intent(in) :: ny !! second sample size
        real(dp) :: pvalue

        integer :: count_extreme
        integer :: total
        integer, allocatable :: choice(:)

        allocate(choice(nx))
        count_extreme = 0
        total = 0
        call enumerate_cvm_choices(1, 0, choice, nx, ny, observed_u, count_extreme, total)
        if (total <= 0) then
            pvalue = quiet_nan(0.0_dp)
        else
            pvalue = real(count_extreme, dp) / real(total, dp)
        end if
    end function cvm2_exact_enumeration

    recursive subroutine enumerate_cvm_choices(start, depth, choice, nx, ny, observed_u, &
            count_extreme, total)
        integer, intent(in) :: start !! smallest pooled rank available at this recursion level
        integer, intent(in) :: depth !! number of selected first-sample ranks
        integer, intent(inout) :: choice(:) !! selected pooled ranks
        integer, intent(in) :: nx !! first sample size
        integer, intent(in) :: ny !! second sample size
        real(dp), intent(in) :: observed_u !! observed U threshold
        integer, intent(inout) :: count_extreme !! number of allocations at least as extreme
        integer, intent(inout) :: total !! total allocations visited

        integer :: i
        integer :: n_total
        real(dp) :: u

        n_total = nx + ny
        if (depth == nx) then
            total = total + 1
            u = cvm_u_from_choice(choice, nx, ny)
            if (u >= observed_u - 64.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(observed_u))) then
                count_extreme = count_extreme + 1
            end if
            return
        end if
        do i = start, n_total - (nx - depth) + 1
            choice(depth + 1) = i
            call enumerate_cvm_choices(i + 1, depth + 1, choice, nx, ny, observed_u, &
                count_extreme, total)
        end do
    end subroutine enumerate_cvm_choices

    pure function cvm_u_from_choice(choice, nx, ny) result(u)
        integer, intent(in) :: choice(:) !! sorted pooled ranks assigned to first sample
        integer, intent(in) :: nx !! first sample size
        integer, intent(in) :: ny !! second sample size
        real(dp) :: u

        integer :: i
        integer :: ix
        integer :: iy
        integer :: rank_value

        u = 0.0_dp
        do i = 1, nx
            u = u + real(nx, dp) * (real(choice(i), dp) - real(i, dp)) ** 2
        end do
        iy = 0
        ix = 1
        do rank_value = 1, nx + ny
            if (ix <= nx) then
                if (choice(ix) == rank_value) then
                    ix = ix + 1
                    cycle
                end if
            end if
            iy = iy + 1
            u = u + real(ny, dp) * (real(rank_value, dp) - real(iy, dp)) ** 2
        end do
    end function cvm_u_from_choice

    subroutine swilk(x, a, w, pw, ifault)
        real(dp), intent(in) :: x(:) !! sorted and median-shifted sample
        real(dp), intent(out) :: a(:) !! Shapiro-Wilk coefficients, length floor(n/2)
        real(dp), intent(out) :: w !! Shapiro-Wilk statistic
        real(dp), intent(out) :: pw !! approximate p-value
        integer, intent(out) :: ifault !! AS R94 status code

        real(dp), parameter :: c1(6) = [0.0_dp, 0.221157_dp, -0.147981_dp, &
            -2.07119_dp, 4.434685_dp, -2.706056_dp]
        real(dp), parameter :: c2(6) = [0.0_dp, 0.042981_dp, -0.293762_dp, &
            -1.752461_dp, 5.682633_dp, -3.582633_dp]
        real(dp), parameter :: c3(4) = [0.5440_dp, -0.39978_dp, 0.025054_dp, -0.0006714_dp]
        real(dp), parameter :: c4(4) = [1.3822_dp, -0.77857_dp, 0.062767_dp, -0.0020322_dp]
        real(dp), parameter :: c5(4) = [-1.5861_dp, -0.31082_dp, -0.083751_dp, 0.0038915_dp]
        real(dp), parameter :: c6(3) = [-0.4803_dp, -0.082676_dp, 0.0030302_dp]
        real(dp), parameter :: g(2) = [-2.273_dp, 0.459_dp]
        real(dp), parameter :: small = 1.0e-19_dp
        integer :: i
        integer :: ind2
        integer :: n
        integer :: n2
        real(dp) :: a1
        real(dp) :: a2
        real(dp) :: an
        real(dp) :: an25
        real(dp) :: asa
        real(dp) :: fac
        real(dp) :: gamma_value
        real(dp) :: m
        real(dp) :: range_value
        real(dp) :: rsn
        real(dp) :: sa
        real(dp) :: sax
        real(dp) :: ssa
        real(dp) :: ssassx
        real(dp) :: ssumm2
        real(dp) :: ssx
        real(dp) :: summ2
        real(dp) :: sx
        real(dp) :: xi
        real(dp) :: xsx
        real(dp) :: xx
        real(dp) :: y
        real(dp) :: s
        real(dp) :: w1

        n = size(x)
        n2 = n / 2
        w = 1.0_dp
        pw = 1.0_dp
        ifault = 0
        an = real(n, dp)
        if (size(a) < n2) then
            ifault = 3
            return
        end if
        if (n < 3) then
            ifault = 1
            return
        end if
        if (n == 3) then
            a(1) = sqrt(0.5_dp)
        else
            an25 = an + 0.25_dp
            summ2 = 0.0_dp
            do i = 1, n2
                a(i) = shapiro_ppnd((real(i, dp) - 0.375_dp) / an25)
                summ2 = summ2 + a(i) ** 2
            end do
            summ2 = 2.0_dp * summ2
            ssumm2 = sqrt(summ2)
            rsn = 1.0_dp / sqrt(an)
            a1 = shapiro_poly(c1, rsn) - a(1) / ssumm2
            if (n > 5) then
                a2 = -a(2) / ssumm2 + shapiro_poly(c2, rsn)
                fac = sqrt((summ2 - 2.0_dp * a(1) ** 2 - 2.0_dp * a(2) ** 2) / &
                    (1.0_dp - 2.0_dp * a1 ** 2 - 2.0_dp * a2 ** 2))
                a(2) = a2
                do i = 3, n2
                    a(i) = -a(i) / fac
                end do
            else
                fac = sqrt((summ2 - 2.0_dp * a(1) ** 2) / (1.0_dp - 2.0_dp * a1 ** 2))
                do i = 2, n2
                    a(i) = -a(i) / fac
                end do
            end if
            a(1) = a1
        end if
        range_value = x(n) - x(1)
        if (range_value < small) then
            ifault = 6
            return
        end if
        xx = x(1) / range_value
        sx = xx
        sa = -a(1)
        ind2 = n - 2
        do i = 2, n
            xi = x(i) / range_value
            sx = sx + xi
            if (i - 1 /= ind2) then
                if (i - 1 < ind2) then
                    sa = sa - a(min(i - 1, ind2) + 1)
                else
                    sa = sa + a(min(i - 1, ind2) + 1)
                end if
            end if
            ind2 = ind2 - 1
        end do
        if (n > 5000) ifault = 2
        sa = sa / an
        sx = sx / an
        ssa = 0.0_dp
        ssx = 0.0_dp
        sax = 0.0_dp
        ind2 = n - 1
        do i = 1, n
            if (i - 1 /= ind2) then
                if (i - 1 < ind2) then
                    asa = -a(min(i - 1, ind2) + 1) - sa
                else
                    asa = a(min(i - 1, ind2) + 1) - sa
                end if
            else
                asa = -sa
            end if
            xsx = x(i) / range_value - sx
            ssa = ssa + asa * asa
            ssx = ssx + xsx * xsx
            sax = sax + asa * xsx
            ind2 = ind2 - 1
        end do
        ssassx = sqrt(ssa * ssx)
        w1 = (ssassx - sax) * (ssassx + sax) / (ssa * ssx)
        w = 1.0_dp - w1
        if (n == 3) then
            if (w < 0.75_dp) then
                w = 0.75_dp
                pw = 0.0_dp
            else
                pw = 1.0_dp - (6.0_dp / scifort_pi) * acos(sqrt(w))
            end if
            return
        end if
        y = log(w1)
        xx = log(an)
        if (n <= 11) then
            gamma_value = shapiro_poly(g, an)
            if (y >= gamma_value) then
                pw = small
                return
            end if
            y = -log(gamma_value - y)
            m = shapiro_poly(c3, an)
            s = exp(shapiro_poly(c4, an))
        else
            m = shapiro_poly(c5, xx)
            s = exp(shapiro_poly(c6, xx))
        end if
        pw = normal_sf((y - m) / s)
    end subroutine swilk

    pure function shapiro_ppnd(p) result(value)
        real(dp), intent(in) :: p !! probability in (0,1)
        real(dp) :: value

        real(dp), parameter :: a0 = 2.50662823884_dp
        real(dp), parameter :: a1 = -18.61500062529_dp
        real(dp), parameter :: a2 = 41.39119773534_dp
        real(dp), parameter :: a3 = -25.44106049637_dp
        real(dp), parameter :: b1 = -8.47351093090_dp
        real(dp), parameter :: b2 = 23.08336743743_dp
        real(dp), parameter :: b3 = -21.06224101826_dp
        real(dp), parameter :: b4 = 3.13082909833_dp
        real(dp), parameter :: c0 = -2.78718931138_dp
        real(dp), parameter :: c1 = -2.29796479134_dp
        real(dp), parameter :: c2 = 4.85014127135_dp
        real(dp), parameter :: c3 = 2.32121276858_dp
        real(dp), parameter :: d1 = 3.54388924762_dp
        real(dp), parameter :: d2 = 1.63706781897_dp
        real(dp) :: q
        real(dp) :: r
        real(dp) :: temp

        q = p - 0.5_dp
        if (abs(q) <= 0.42_dp) then
            r = q * q
            temp = q * (((a3 * r + a2) * r + a1) * r + a0)
            temp = temp / ((((b4 * r + b3) * r + b2) * r + b1) * r + 1.0_dp)
            value = temp
            return
        end if
        r = p
        if (q > 0.0_dp) r = 1.0_dp - p
        if (r <= 0.0_dp) then
            value = 0.0_dp
            return
        end if
        r = sqrt(-log(r))
        temp = ((c3 * r + c2) * r + c1) * r + c0
        temp = temp / ((d2 * r + d1) * r + 1.0_dp)
        if (q < 0.0_dp) then
            value = -temp
        else
            value = temp
        end if
    end function shapiro_ppnd

    pure function shapiro_poly(c, x) result(value)
        real(dp), intent(in) :: c(:) !! AS R94 polynomial coefficients
        real(dp), intent(in) :: x !! polynomial argument
        real(dp) :: value

        integer :: i
        real(dp) :: p

        value = c(1)
        if (size(c) == 1) return
        p = x * c(size(c))
        if (size(c) == 2) then
            value = value + p
            return
        end if
        do i = size(c) - 1, 2, -1
            p = (p + c(i)) * x
        end do
        value = value + p
    end function shapiro_poly

    subroutine logistic_mle(x, loc, scale)
        real(dp), intent(in) :: x(:) !! finite observations
        real(dp), intent(out) :: loc !! fitted logistic location
        real(dp), intent(out) :: scale !! fitted positive logistic scale

        integer :: iter
        real(dp) :: det
        real(dp) :: epsa
        real(dp) :: epsb
        real(dp) :: f1
        real(dp) :: f2
        real(dp) :: f1a
        real(dp) :: f1b
        real(dp) :: f2a
        real(dp) :: f2b
        real(dp) :: da
        real(dp) :: db

        loc = mean(x)
        scale = standard_deviation(x, 1)
        if (.not. ieee_is_finite(scale) .or. scale <= 0.0_dp) return
        scale = max(scale, sqrt(tiny(1.0_dp)))
        do iter = 1, 100
            call logistic_score_equations(x, loc, scale, f1, f2)
            if (max(abs(f1), abs(f2)) < 1.0e-10_dp * real(size(x), dp)) exit
            epsa = sqrt(epsilon(1.0_dp)) * max(1.0_dp, abs(loc))
            epsb = sqrt(epsilon(1.0_dp)) * max(1.0_dp, scale)
            call logistic_score_equations(x, loc + epsa, scale, f1a, f2a)
            call logistic_score_equations(x, loc, scale + epsb, f1b, f2b)
            f1a = (f1a - f1) / epsa
            f2a = (f2a - f2) / epsa
            f1b = (f1b - f1) / epsb
            f2b = (f2b - f2) / epsb
            det = f1a * f2b - f1b * f2a
            if (abs(det) <= tiny(1.0_dp)) exit
            da = (-f1 * f2b + f1b * f2) / det
            db = (-f1a * f2 + f1 * f2a) / det
            if (scale + db <= 0.0_dp) db = -0.5_dp * scale
            loc = loc + da
            scale = scale + db
            if (max(abs(da), abs(db)) < 1.0e-12_dp * max(1.0_dp, scale)) exit
        end do
    end subroutine logistic_mle

    subroutine logistic_score_equations(x, loc, scale, f1, f2)
        real(dp), intent(in) :: x(:) !! finite observations
        real(dp), intent(in) :: loc !! trial location
        real(dp), intent(in) :: scale !! trial positive scale
        real(dp), intent(out) :: f1 !! first likelihood equation
        real(dp), intent(out) :: f2 !! second likelihood equation

        integer :: i
        real(dp) :: e
        real(dp) :: t
        real(dp) :: inv

        f1 = -0.5_dp * real(size(x), dp)
        f2 = real(size(x), dp)
        do i = 1, size(x)
            t = (x(i) - loc) / scale
            if (t > 40.0_dp) then
                inv = exp(-t)
                e = 1.0_dp / (1.0_dp + inv)
            else if (t < -40.0_dp) then
                e = exp(t)
            else
                e = 1.0_dp / (1.0_dp + exp(-t))
            end if
            f1 = f1 + (1.0_dp - e)
            f2 = f2 + t * (1.0_dp - 2.0_dp * e)
        end do
    end subroutine logistic_score_equations

    subroutine gumbel_r_mle(x, loc, scale)
        real(dp), intent(in) :: x(:) !! finite observations
        real(dp), intent(out) :: loc !! fitted Gumbel-right location
        real(dp), intent(out) :: scale !! fitted positive scale

        integer :: iter
        real(dp) :: fl
        real(dp) :: fr
        real(dp) :: fm
        real(dp) :: left
        real(dp) :: mid
        real(dp) :: right

        left = 0.5_dp
        right = 2.0_dp
        call gumbel_scale_equation(x, left, fl)
        call gumbel_scale_equation(x, right, fr)
        do while (sign(1.0_dp, fl) == sign(1.0_dp, fr))
            left = left / 2.0_dp
            right = right * 2.0_dp
            call gumbel_scale_equation(x, left, fl)
            call gumbel_scale_equation(x, right, fr)
            if (right > huge(1.0_dp) ** 0.25_dp) exit
        end do
        do iter = 1, 200
            mid = 0.5_dp * (left + right)
            call gumbel_scale_equation(x, mid, fm)
            if (abs(fm) <= 1.0e-13_dp * max(1.0_dp, abs(mean(x)))) exit
            if (sign(1.0_dp, fm) == sign(1.0_dp, fl)) then
                left = mid
                fl = fm
            else
                right = mid
                fr = fm
            end if
        end do
        scale = 0.5_dp * (left + right)
        loc = gumbel_loc_from_scale(x, scale)
    end subroutine gumbel_r_mle

    subroutine gumbel_l_mle(x, loc, scale)
        real(dp), intent(in) :: x(:) !! finite observations
        real(dp), intent(out) :: loc !! fitted Gumbel-left location
        real(dp), intent(out) :: scale !! fitted positive scale

        real(dp), allocatable :: negx(:)
        real(dp) :: rloc

        negx = -x
        call gumbel_r_mle(negx, rloc, scale)
        loc = -rloc
    end subroutine gumbel_l_mle

    subroutine gumbel_scale_equation(x, scale, value)
        real(dp), intent(in) :: x(:) !! finite observations
        real(dp), intent(in) :: scale !! positive trial scale
        real(dp), intent(out) :: value !! scale likelihood equation

        integer :: i
        real(dp) :: max_logw
        real(dp) :: sumw
        real(dp) :: sumwx

        max_logw = maxval(-x / scale)
        sumw = 0.0_dp
        sumwx = 0.0_dp
        do i = 1, size(x)
            sumw = sumw + exp(-x(i) / scale - max_logw)
            sumwx = sumwx + x(i) * exp(-x(i) / scale - max_logw)
        end do
        value = mean(x) - sumwx / sumw - scale
    end subroutine gumbel_scale_equation

    function gumbel_loc_from_scale(x, scale) result(loc)
        real(dp), intent(in) :: x(:) !! finite observations
        real(dp), intent(in) :: scale !! fitted positive scale
        real(dp) :: loc

        integer :: i
        real(dp) :: max_logw
        real(dp) :: sumw

        max_logw = maxval(-x / scale)
        sumw = 0.0_dp
        do i = 1, size(x)
            sumw = sumw + exp(-x(i) / scale - max_logw)
        end do
        loc = -scale * (max_logw + log(sumw) - log(real(size(x), dp)))
    end function gumbel_loc_from_scale

    function interpolate_table(x, grid_x, grid_y) result(value)
        real(dp), intent(in) :: x !! statistic to interpolate
        real(dp), intent(in) :: grid_x(:) !! increasing critical values
        real(dp), intent(in) :: grid_y(:) !! corresponding p-values
        real(dp) :: value

        integer :: i

        if (x <= grid_x(1)) then
            value = grid_y(1)
            return
        else if (x >= grid_x(size(grid_x))) then
            value = grid_y(size(grid_y))
            return
        end if
        do i = 1, size(grid_x) - 1
            if (x <= grid_x(i + 1)) then
                value = grid_y(i) + (grid_y(i + 1) - grid_y(i)) * &
                    (x - grid_x(i)) / (grid_x(i + 1) - grid_x(i))
                return
            end if
        end do
        value = grid_y(size(grid_y))
    end function interpolate_table

    pure function finite_sample(x) result(ok)
        real(dp), intent(in) :: x(:) !! candidate finite observations
        logical :: ok

        integer :: i

        ok = .true.
        do i = 1, size(x)
            if (.not. ieee_is_finite(x(i))) then
                ok = .false.
                return
            end if
        end do
    end function finite_sample

    pure function valid_alternative(alt) result(ok)
        character(len=*), intent(in) :: alt !! normalized alternative name
        logical :: ok

        ok = trim(alt) == 'two-sided' .or. trim(alt) == 'less' .or. trim(alt) == 'greater'
    end function valid_alternative

    pure function select_alternative(alternative) result(alt)
        character(len=*), intent(in), optional :: alternative !! requested alternative
        character(len=9) :: alt

        alt = 'two-sided'
        if (present(alternative)) alt = lower_string(adjustl(alternative))
        if (alt(1:1) == 't') alt = 'two-sided'
        if (alt(1:1) == 'g') alt = 'greater'
        if (alt(1:1) == 'l') alt = 'less'
    end function select_alternative

    pure function select_method(method, default_value) result(selected)
        character(len=*), intent(in), optional :: method !! requested method
        character(len=*), intent(in) :: default_value !! default method
        character(len=8) :: selected

        selected = default_value
        if (present(method)) selected = lower_string(adjustl(method))
    end function select_method

    pure function valid_ks1_method(method) result(ok)
        character(len=*), intent(in) :: method !! normalized method
        logical :: ok

        ok = trim(method) == 'auto' .or. trim(method) == 'exact' .or. &
            trim(method) == 'approx' .or. trim(method) == 'asymp'
    end function valid_ks1_method

    pure function lower_string(text) result(lower)
        character(len=*), intent(in) :: text !! text converted to lower case
        character(len=len(text)) :: lower

        integer :: c
        integer :: i

        lower = text
        do i = 1, len(text)
            c = iachar(lower(i:i))
            if (c >= iachar('A') .and. c <= iachar('Z')) lower(i:i) = achar(c + 32)
        end do
    end function lower_string

    pure function gcd_int(a, b) result(g)
        integer, intent(in) :: a !! positive integer
        integer, intent(in) :: b !! positive integer
        integer :: g

        integer :: x
        integer :: y
        integer :: t

        x = abs(a)
        y = abs(b)
        do while (y /= 0)
            t = modulo(x, y)
            x = y
            y = t
        end do
        g = x
    end function gcd_int

    pure function combination_leq(n, k, limit) result(ok)
        integer, intent(in) :: n !! total items
        integer, intent(in) :: k !! selected items
        integer(kind=8), intent(in) :: limit !! maximum accepted combination count
        logical :: ok

        integer :: i
        integer :: kk
        real(dp) :: value

        kk = min(k, n - k)
        value = 1.0_dp
        do i = 1, kk
            value = value * real(n - kk + i, dp) / real(i, dp)
            if (value > real(limit, dp)) then
                ok = .false.
                return
            end if
        end do
        ok = .true.
    end function combination_leq

    subroutine invalidate_significance(res)
        type(significance_result), intent(out) :: res !! result filled with NaNs

        res%statistic = quiet_nan(0.0_dp)
        res%pvalue = quiet_nan(0.0_dp)
    end subroutine invalidate_significance

    subroutine invalidate_ks(res)
        type(kstest_result), intent(out) :: res !! result filled with invalid sentinels

        res%statistic = quiet_nan(0.0_dp)
        res%pvalue = quiet_nan(0.0_dp)
        res%statistic_location = quiet_nan(0.0_dp)
        res%statistic_sign = 0
    end subroutine invalidate_ks

    subroutine invalidate_anderson(res)
        type(anderson_result), intent(out) :: res !! result filled with invalid sentinels

        res%statistic = quiet_nan(0.0_dp)
        res%pvalue = quiet_nan(0.0_dp)
        res%status = 1
        if (allocated(res%critical_values)) deallocate(res%critical_values)
        if (allocated(res%significance_level)) deallocate(res%significance_level)
        if (allocated(res%fit_params)) deallocate(res%fit_params)
    end subroutine invalidate_anderson

    subroutine sort_real(values)
        real(dp), intent(inout) :: values(:) !! values sorted into increasing order

        if (size(values) > 1) call quicksort_real(values, 1, size(values))
    end subroutine sort_real

    recursive subroutine quicksort_real(values, left, right)
        real(dp), intent(inout) :: values(:) !! array being sorted
        integer, intent(in) :: left !! first index of this partition
        integer, intent(in) :: right !! last index of this partition

        integer :: i
        integer :: j
        real(dp) :: pivot
        real(dp) :: tmp

        i = left
        j = right
        pivot = values((left + right) / 2)
        do
            do while (values(i) < pivot)
                i = i + 1
            end do
            do while (values(j) > pivot)
                j = j - 1
            end do
            if (i <= j) then
                tmp = values(i)
                values(i) = values(j)
                values(j) = tmp
                i = i + 1
                j = j - 1
            end if
            if (i > j) exit
        end do
        if (left < j) call quicksort_real(values, left, j)
        if (i < right) call quicksort_real(values, i, right)
    end subroutine quicksort_real

end module scifort_goodness_of_fit
