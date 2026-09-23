! SPDX-License-Identifier: MIT AND BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
! Adapted portions: SciPy 1.17.0 scipy/stats/_mstats_basic.py (Kendall exact recurrence).
! Retained BSD terms: THIRD_PARTY_LICENSES.md; details: CODE_PROVENANCE.md.
!
! Association, robust-regression, and ordered-trend procedures following
! scipy.stats 1.17.0 scalar semantics where practical.

module scifort_association_extended
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_descriptive, only : covariance, mean, median, pearson_correlation, &
        rankdata, variance
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    use scifort_normal, only : normal_cdf, normal_ppf, normal_sf
    use scifort_student_t, only : t_cdf, t_sf
    implicit none
    private

    type, public :: association_result
        real(dp) :: statistic !! association or test statistic
        real(dp) :: pvalue !! p-value for the requested alternative
    end type association_result

    type, public :: linregress_result
        real(dp) :: slope !! least-squares slope
        real(dp) :: intercept !! least-squares intercept
        real(dp) :: rvalue !! Pearson correlation coefficient
        real(dp) :: pvalue !! p-value for a zero-slope test
        real(dp) :: stderr !! standard error of the slope
        real(dp) :: intercept_stderr !! standard error of the intercept
    end type linregress_result

    type, public :: theilslopes_result
        real(dp) :: slope !! Theil-Sen median slope
        real(dp) :: intercept !! selected robust intercept estimate
        real(dp) :: low_slope !! lower confidence bound for the slope
        real(dp) :: high_slope !! upper confidence bound for the slope
    end type theilslopes_result

    type, public :: siegelslopes_result
        real(dp) :: slope !! Siegel repeated-median slope
        real(dp) :: intercept !! selected robust intercept estimate
    end type siegelslopes_result

    type, public :: page_trend_result
        real(dp) :: statistic !! Page L statistic
        real(dp) :: pvalue !! exact or asymptotic one-sided p-value
        character(len=10) :: method !! exact or asymptotic
    end type page_trend_result

    public :: brunnermunzel
    public :: kendalltau
    public :: linregress
    public :: page_trend_test
    public :: pointbiserialr
    public :: siegelslopes
    public :: theilslopes

contains

    function pointbiserialr(x, y) result(result_value)
        real(dp), intent(in) :: x(:) !! binary-coded observations; SciPy treats these as numeric
        real(dp), intent(in) :: y(:) !! continuous observations, same size as x
        type(association_result) :: result_value

        real(dp) :: r
        real(dp) :: tstat
        real(dp) :: df

        call invalidate_association(result_value)
        if (size(x) /= size(y) .or. size(x) < 2) return
        if (any(ieee_is_nan(x)) .or. any(ieee_is_nan(y))) return

        r = pearson_correlation(x, y)
        result_value%statistic = r
        if (ieee_is_nan(r)) return
        if (size(x) == 2) then
            if (abs(r) >= 1.0_dp) then
                result_value%pvalue = 1.0_dp
            else
                result_value%pvalue = quiet_nan(0.0_dp)
            end if
            return
        end if

        df = real(size(x) - 2, dp)
        if (abs(r) >= 1.0_dp) then
            result_value%pvalue = 0.0_dp
        else
            tstat = r * sqrt(df / max(tiny(1.0_dp), (1.0_dp - r) * (1.0_dp + r)))
            result_value%pvalue = min(1.0_dp, 2.0_dp * t_sf(abs(tstat), df))
        end if
    end function pointbiserialr

    function kendalltau(x, y, method, variant, alternative) result(result_value)
        real(dp), intent(in) :: x(:) !! first ordinal observations
        real(dp), intent(in) :: y(:) !! second ordinal observations, same size as x
        character(len=*), intent(in), optional :: method !! auto, exact, or asymptotic
        character(len=*), intent(in), optional :: variant !! b or c
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        type(association_result) :: result_value

        character(len=10) :: selected_method
        character(len=9) :: selected_alternative
        character(len=1) :: selected_variant
        integer(int64) :: concordant
        integer(int64) :: discordant
        integer :: i
        integer :: j
        integer :: minclasses
        integer :: n
        integer(int64) :: ntie
        integer(int64) :: tot
        integer(int64) :: x0
        integer(int64) :: x1
        integer(int64) :: xtie
        integer(int64) :: y0
        integer(int64) :: y1
        integer(int64) :: ytie
        real(dp) :: con_minus_dis
        real(dp) :: m
        real(dp) :: var_s
        real(dp) :: z

        call invalidate_association(result_value)
        if (size(x) /= size(y) .or. size(x) < 2) return
        if (any(ieee_is_nan(x)) .or. any(ieee_is_nan(y))) return

        selected_method = 'auto'
        if (present(method)) selected_method = adjustl(method)
        selected_variant = 'b'
        if (present(variant)) selected_variant = variant(1:1)
        selected_alternative = 'two-sided'
        if (present(alternative)) selected_alternative = adjustl(alternative)
        if (.not. valid_alternative(selected_alternative)) return
        if (selected_method /= 'auto' .and. selected_method /= 'exact' .and. &
                selected_method /= 'asymptotic') return
        if (selected_variant /= 'b' .and. selected_variant /= 'c') return

        n = size(x)
        concordant = 0
        discordant = 0
        ntie = 0
        do i = 1, n - 1
            do j = i + 1, n
                if (same_value(x(i), x(j)) .and. same_value(y(i), y(j))) then
                    ntie = ntie + 1
                else if (.not. same_value(x(i), x(j)) .and. &
                        .not. same_value(y(i), y(j))) then
                    if ((x(i) - x(j)) * (y(i) - y(j)) > 0.0_dp) then
                        concordant = concordant + 1
                    else
                        discordant = discordant + 1
                    end if
                end if
            end do
        end do

        call tie_statistics(x, xtie, x0, x1, minclasses)
        call tie_statistics(y, ytie, y0, y1, i)
        minclasses = min(minclasses, i)
        tot = int(n, int64) * int(n - 1, int64) / 2_int64
        if (xtie == tot .or. ytie == tot) return

        con_minus_dis = real(concordant - discordant, dp)
        if (selected_variant == 'b') then
            result_value%statistic = con_minus_dis / &
                sqrt(real(tot - xtie, dp) * real(tot - ytie, dp))
        else
            if (minclasses <= 1) return
            result_value%statistic = 2.0_dp * con_minus_dis / &
                (real(n, dp) ** 2 * real(minclasses - 1, dp) / real(minclasses, dp))
        end if
        result_value%statistic = max(-1.0_dp, min(1.0_dp, result_value%statistic))

        if (selected_method == 'exact' .and. (xtie /= 0 .or. ytie /= 0)) return
        if (selected_method == 'auto') then
            if (xtie == 0 .and. ytie == 0 .and. &
                    (n <= 33 .or. min(discordant, tot - discordant) <= 1_int64)) then
                selected_method = 'exact'
            else
                selected_method = 'asymptotic'
            end if
        end if

        if (selected_method == 'exact') then
            if (concordant > int(huge(0), int64)) return
            result_value%pvalue = kendall_p_exact(n, int(concordant), selected_alternative)
        else
            m = real(n, dp) * real(n - 1, dp)
            var_s = (m * real(2 * n + 5, dp) - real(x1 + y1, dp)) / 18.0_dp
            var_s = var_s + 2.0_dp * real(xtie, dp) * real(ytie, dp) / m
            if (n > 2) then
                var_s = var_s + real(x0, dp) * real(y0, dp) / &
                    (9.0_dp * m * real(n - 2, dp))
            end if
            if (var_s <= 0.0_dp) then
                result_value%pvalue = quiet_nan(0.0_dp)
            else
                z = con_minus_dis / sqrt(var_s)
                result_value%pvalue = normal_pvalue(z, selected_alternative)
            end if
        end if
    end function kendalltau

    function linregress(x, y, alternative) result(result_value)
        real(dp), intent(in) :: x(:) !! independent variable
        real(dp), intent(in) :: y(:) !! dependent variable, same size as x
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        type(linregress_result) :: result_value

        character(len=9) :: selected_alternative
        integer :: n
        real(dp), parameter :: tiny_r = 1.0e-20_dp
        real(dp) :: df
        real(dp) :: r
        real(dp) :: ssxm
        real(dp) :: ssxym
        real(dp) :: ssym
        real(dp) :: tstat
        real(dp) :: xmean
        real(dp) :: ymean

        call invalidate_linregress(result_value)
        selected_alternative = 'two-sided'
        if (present(alternative)) selected_alternative = adjustl(alternative)
        if (.not. valid_alternative(selected_alternative)) return
        if (size(x) /= size(y) .or. size(x) < 2) return
        if (any(ieee_is_nan(x)) .or. any(ieee_is_nan(y))) return
        if (maxval(x) == minval(x)) return

        n = size(x)
        xmean = mean(x)
        ymean = mean(y)
        ssxm = variance(x, 0)
        ssym = variance(y, 0)
        ssxym = covariance(x, y, 0)
        if (ssxm <= 0.0_dp) return

        if (ssxm == 0.0_dp .or. ssym == 0.0_dp) then
            if (ssxym == 0.0_dp) then
                r = quiet_nan(0.0_dp)
            else
                r = 0.0_dp
            end if
        else
            r = ssxym / sqrt(ssxm * ssym)
            r = max(-1.0_dp, min(1.0_dp, r))
        end if

        result_value%slope = ssxym / ssxm
        result_value%intercept = ymean - result_value%slope * xmean
        result_value%rvalue = r

        if (n == 2) then
            if (y(1) == y(2)) then
                result_value%pvalue = 1.0_dp
            else
                result_value%pvalue = 0.0_dp
            end if
            result_value%stderr = 0.0_dp
            result_value%intercept_stderr = 0.0_dp
            return
        end if

        df = real(n - 2, dp)
        if (ieee_is_nan(r)) then
            result_value%pvalue = quiet_nan(0.0_dp)
            result_value%stderr = quiet_nan(0.0_dp)
            result_value%intercept_stderr = quiet_nan(0.0_dp)
            return
        end if
        tstat = r * sqrt(df / ((1.0_dp - r + tiny_r) * (1.0_dp + r + tiny_r)))
        result_value%pvalue = student_t_pvalue(tstat, df, selected_alternative)
        result_value%stderr = sqrt(max(0.0_dp, (1.0_dp - r ** 2) * ssym / ssxm / df))
        result_value%intercept_stderr = result_value%stderr * sqrt(ssxm + xmean ** 2)
    end function linregress

    function theilslopes(y, x, alpha, method) result(result_value)
        real(dp), intent(in) :: y(:) !! dependent variable
        real(dp), intent(in), optional :: x(:) !! independent variable; default 0,1,...
        real(dp), intent(in), optional :: alpha !! confidence degree, default 0.95
        character(len=*), intent(in), optional :: method !! separate or joint
        type(theilslopes_result) :: result_value

        character(len=8) :: selected_method
        integer :: i
        integer :: j
        integer :: low_index
        integer :: n
        integer :: nslopes
        integer :: position
        integer :: upper_index
        real(dp) :: a
        real(dp), allocatable :: intercept_data(:)
        real(dp), allocatable :: slopes(:)
        real(dp), allocatable :: xv(:)
        real(dp) :: sigsq
        real(dp) :: sigma
        real(dp) :: z

        call invalidate_theil(result_value)
        n = size(y)
        if (n < 2) return
        if (any(ieee_is_nan(y))) return
        selected_method = 'separate'
        if (present(method)) selected_method = adjustl(method)
        if (selected_method /= 'separate' .and. selected_method /= 'joint') return

        allocate(xv(n))
        if (present(x)) then
            if (size(x) /= n .or. any(ieee_is_nan(x))) return
            xv = x
        else
            do i = 1, n
                xv(i) = real(i - 1, dp)
            end do
        end if

        nslopes = 0
        do i = 1, n
            do j = 1, n
                if (xv(i) > xv(j)) nslopes = nslopes + 1
            end do
        end do
        if (nslopes <= 0) return
        allocate(slopes(nslopes))
        position = 0
        do i = 1, n
            do j = 1, n
                if (xv(i) > xv(j)) then
                    position = position + 1
                    slopes(position) = (y(i) - y(j)) / (xv(i) - xv(j))
                end if
            end do
        end do
        call sort_real(slopes)
        result_value%slope = median_sorted(slopes)

        if (selected_method == 'joint') then
            allocate(intercept_data(n))
            intercept_data = y - result_value%slope * xv
            result_value%intercept = median(intercept_data)
        else
            result_value%intercept = median(y) - result_value%slope * median(xv)
        end if

        a = 0.95_dp
        if (present(alpha)) a = alpha
        if (.not. ieee_is_finite(a) .or. a <= 0.0_dp .or. a >= 1.0_dp) return
        if (a > 0.5_dp) a = 1.0_dp - a
        z = normal_ppf(a / 2.0_dp)
        sigsq = real(n * (n - 1) * (2 * n + 5), dp)
        sigsq = sigsq - tie_cubic5_sum(xv) - tie_cubic5_sum(y)
        sigsq = sigsq / 18.0_dp
        if (sigsq < 0.0_dp) return
        sigma = sqrt(sigsq)
        upper_index = min(round_nearest_even((real(nslopes, dp) - z * sigma) / 2.0_dp), &
            nslopes - 1) + 1
        low_index = max(round_nearest_even((real(nslopes, dp) + z * sigma) / 2.0_dp) - 1, 0) + 1
        result_value%low_slope = slopes(low_index)
        result_value%high_slope = slopes(upper_index)
    end function theilslopes

    function siegelslopes(y, x, method) result(result_value)
        real(dp), intent(in) :: y(:) !! dependent variable
        real(dp), intent(in), optional :: x(:) !! independent variable; default 0,1,...
        character(len=*), intent(in), optional :: method !! hierarchical or separate
        type(siegelslopes_result) :: result_value

        character(len=12) :: selected_method
        integer :: count_nonzero
        integer :: i
        integer :: j
        integer :: n
        integer :: position
        real(dp), allocatable :: intercepts(:)
        real(dp), allocatable :: local_values(:)
        real(dp), allocatable :: slopes(:)
        real(dp), allocatable :: xv(:)

        call invalidate_siegel(result_value)
        n = size(y)
        if (n < 2 .or. any(ieee_is_nan(y))) return
        selected_method = 'hierarchical'
        if (present(method)) selected_method = adjustl(method)
        if (selected_method /= 'hierarchical' .and. selected_method /= 'separate') return

        allocate(xv(n))
        if (present(x)) then
            if (size(x) /= n .or. any(ieee_is_nan(x))) return
            xv = x
        else
            do i = 1, n
                xv(i) = real(i - 1, dp)
            end do
        end if

        allocate(slopes(n))
        if (selected_method == 'separate') allocate(intercepts(n))
        do i = 1, n
            count_nonzero = count(xv(i) /= xv)
            if (count_nonzero <= 0) return
            allocate(local_values(count_nonzero))
            position = 0
            do j = 1, n
                if (xv(i) /= xv(j)) then
                    position = position + 1
                    local_values(position) = (y(i) - y(j)) / (xv(i) - xv(j))
                end if
            end do
            slopes(i) = median(local_values)
            if (selected_method == 'separate') then
                position = 0
                do j = 1, n
                    if (xv(i) /= xv(j)) then
                        position = position + 1
                        local_values(position) = &
                            (y(j) * xv(i) - y(i) * xv(j)) / (xv(i) - xv(j))
                    end if
                end do
                intercepts(i) = median(local_values)
            end if
            deallocate(local_values)
        end do

        result_value%slope = median(slopes)
        if (selected_method == 'separate') then
            result_value%intercept = median(intercepts)
        else
            allocate(local_values(n))
            local_values = y - result_value%slope * xv
            result_value%intercept = median(local_values)
        end if
    end function siegelslopes

    function brunnermunzel(x, y, alternative, distribution) result(result_value)
        real(dp), intent(in) :: x(:) !! first independent sample
        real(dp), intent(in) :: y(:) !! second independent sample
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        character(len=*), intent(in), optional :: distribution !! t or normal
        type(association_result) :: result_value

        character(len=9) :: selected_alternative
        character(len=6) :: selected_distribution
        integer :: nx
        integer :: ny
        real(dp) :: df
        real(dp) :: df_denom
        real(dp) :: df_numer
        real(dp), allocatable :: pooled(:)
        real(dp), allocatable :: rankc(:)
        real(dp), allocatable :: rankcx(:)
        real(dp), allocatable :: rankcy(:)
        real(dp), allocatable :: rankx(:)
        real(dp), allocatable :: ranky(:)
        real(dp) :: rankcx_mean
        real(dp) :: rankcy_mean
        real(dp) :: rankx_mean
        real(dp) :: ranky_mean
        real(dp) :: sx
        real(dp) :: sy
        real(dp) :: w

        call invalidate_association(result_value)
        nx = size(x)
        ny = size(y)
        if (nx < 2 .or. ny < 2 .or. any(ieee_is_nan(x)) .or. any(ieee_is_nan(y))) return
        selected_alternative = 'two-sided'
        if (present(alternative)) selected_alternative = adjustl(alternative)
        selected_distribution = 't'
        if (present(distribution)) selected_distribution = adjustl(distribution)
        if (.not. valid_alternative(selected_alternative)) return
        if (selected_distribution /= 't' .and. selected_distribution /= 'normal') return

        allocate(pooled(nx + ny))
        pooled(1:nx) = x
        pooled(nx + 1:) = y
        rankc = rankdata(pooled)
        rankcx = rankc(1:nx)
        rankcy = rankc(nx + 1:)
        rankx = rankdata(x)
        ranky = rankdata(y)
        rankcx_mean = mean(rankcx)
        rankcy_mean = mean(rankcy)
        rankx_mean = mean(rankx)
        ranky_mean = mean(ranky)
        sx = sum((rankcx - rankx - rankcx_mean + rankx_mean) ** 2) / real(nx - 1, dp)
        sy = sum((rankcy - ranky - rankcy_mean + ranky_mean) ** 2) / real(ny - 1, dp)
        if (real(nx, dp) * sx + real(ny, dp) * sy <= 0.0_dp) return

        w = real(nx * ny, dp) * (rankcy_mean - rankcx_mean) / &
            (real(nx + ny, dp) * sqrt(real(nx, dp) * sx + real(ny, dp) * sy))
        result_value%statistic = w
        if (selected_distribution == 't') then
            df_numer = (real(nx, dp) * sx + real(ny, dp) * sy) ** 2
            df_denom = (real(nx, dp) * sx) ** 2 / real(nx - 1, dp) + &
                (real(ny, dp) * sy) ** 2 / real(ny - 1, dp)
            if (df_denom <= 0.0_dp) return
            df = df_numer / df_denom
            result_value%pvalue = student_t_pvalue(-w, df, selected_alternative)
        else
            result_value%pvalue = normal_pvalue(-w, selected_alternative)
        end if
    end function brunnermunzel

    function page_trend_test(data, ranked, predicted_ranks, method) result(result_value)
        real(dp), intent(in) :: data(:, :) !! subjects by treatments
        logical, intent(in), optional :: ranked !! data already contain row ranks
        integer, intent(in), optional :: predicted_ranks(:) !! permutation of 1..n treatments
        character(len=*), intent(in), optional :: method !! auto, exact, or asymptotic
        type(page_trend_result) :: result_value

        logical :: already_ranked
        character(len=10) :: selected_method
        integer :: i
        integer :: j
        integer :: m
        integer :: n
        integer, allocatable :: expected_order(:)
        real(dp) :: e0
        real(dp) :: lambda
        real(dp), allocatable :: ranks(:, :)
        real(dp) :: v0

        call invalidate_page(result_value)
        m = size(data, 1)
        n = size(data, 2)
        if (m < 2 .or. n < 3 .or. any(ieee_is_nan(data))) return
        already_ranked = .false.
        if (present(ranked)) already_ranked = ranked
        selected_method = 'auto'
        if (present(method)) selected_method = adjustl(method)
        if (selected_method /= 'auto' .and. selected_method /= 'exact' .and. &
                selected_method /= 'asymptotic') return

        allocate(expected_order(n))
        if (present(predicted_ranks)) then
            if (size(predicted_ranks) /= n) return
            expected_order = predicted_ranks
            if (.not. valid_permutation(expected_order)) return
        else
            do j = 1, n
                expected_order(j) = j
            end do
        end if

        allocate(ranks(m, n))
        if (already_ranked) then
            if (minval(data) < 1.0_dp .or. maxval(data) > real(n, dp)) return
            ranks = data
        else
            do i = 1, m
                ranks(i, :) = rankdata(data(i, :))
            end do
        end if

        result_value%statistic = 0.0_dp
        do j = 1, n
            result_value%statistic = result_value%statistic + &
                real(expected_order(j), dp) * sum(ranks(:, j))
        end do

        if (selected_method == 'auto') then
            if (n > 8 .or. (m > 12 .and. n > 3) .or. m > 20) then
                selected_method = 'asymptotic'
            else
                selected_method = 'exact'
            end if
        end if
        result_value%method = selected_method

        if (selected_method == 'asymptotic') then
            e0 = real(m * n * (n + 1) ** 2, dp) / 4.0_dp
            v0 = real(m * n ** 2 * (n + 1) * (n ** 2 - 1), dp) / 144.0_dp
            if (v0 <= 0.0_dp) return
            lambda = (result_value%statistic - e0) / sqrt(v0)
            result_value%pvalue = normal_sf(lambda)
        else
            result_value%pvalue = page_exact_sf(int(result_value%statistic), m, n)
        end if
    end function page_trend_test

    function kendall_p_exact(n, c_input, alternative) result(probability)
        integer, intent(in) :: n !! sample size
        integer, intent(in) :: c_input !! number of concordant pairs
        character(len=*), intent(in) :: alternative !! requested alternative
        real(dp) :: probability

        integer :: c
        integer :: i
        integer :: j
        integer :: total
        logical :: in_right_tail
        logical :: alternative_greater
        real(dp) :: fact
        real(dp) :: p_mass
        real(dp), allocatable :: new(:)
        real(dp), allocatable :: temp(:)
        real(dp), allocatable :: cumulative(:)

        probability = quiet_nan(0.0_dp)
        if (n <= 0) return
        total = n * (n - 1) / 2
        c = min(c_input, total - c_input)
        if (c < 0 .or. 2 * c > total) return
        in_right_tail = c_input >= total - c_input
        alternative_greater = alternative == 'greater'

        if (n == 1) then
            probability = 1.0_dp
            p_mass = 1.0_dp
        else if (n == 2) then
            probability = 1.0_dp
            p_mass = 0.5_dp
        else if (c == 0) then
            if (n < 171) then
                fact = factorial_real(n)
                probability = 2.0_dp / fact
                p_mass = probability / 2.0_dp
            else
                probability = 0.0_dp
                p_mass = 0.0_dp
            end if
        else if (c == 1) then
            if (n < 172) then
                fact = factorial_real(n)
                probability = 2.0_dp * real(n, dp) / fact
                p_mass = real(n - 1, dp) / fact
            else
                probability = 0.0_dp
                p_mass = 0.0_dp
            end if
        else if (2 * c == total .and. alternative == 'two-sided') then
            probability = 1.0_dp
            return
        else
            allocate(new(0:c), temp(0:c), cumulative(0:c))
            new = 0.0_dp
            new(0:min(1, c)) = 1.0_dp
            if (n < 171) then
                do j = 3, n
                    temp = 0.0_dp
                    temp(0) = new(0)
                    do i = 1, c
                        temp(i) = temp(i - 1) + new(i)
                    end do
                    if (j <= c) then
                        cumulative = temp
                        do i = j, c
                            temp(i) = temp(i) - cumulative(i - j)
                        end do
                    end if
                    new = temp
                end do
                fact = factorial_real(n)
                probability = 2.0_dp * sum(new) / fact
                p_mass = new(c) / fact
            else
                do j = 3, n
                    temp = 0.0_dp
                    temp(0) = new(0) / real(j, dp)
                    do i = 1, c
                        temp(i) = (temp(i - 1) * real(j, dp) + new(i)) / real(j, dp)
                    end do
                    if (j <= c) then
                        cumulative = temp
                        do i = j, c
                            temp(i) = temp(i) - cumulative(i - j)
                        end do
                    end if
                    new = temp
                end do
                probability = sum(new)
                p_mass = new(c) / 2.0_dp
            end if
        end if

        if (alternative /= 'two-sided') then
            if (in_right_tail .eqv. alternative_greater) then
                probability = probability / 2.0_dp
            else
                probability = 1.0_dp - probability / 2.0_dp + p_mass
            end if
        end if
        probability = max(0.0_dp, min(1.0_dp, probability))
    end function kendall_p_exact

    subroutine tie_statistics(values, pair_ties, cubic0, cubic1, nclasses)
        real(dp), intent(in) :: values(:) !! values without NaNs
        integer(int64), intent(out) :: pair_ties !! sum choose(group_size, 2)
        integer(int64), intent(out) :: cubic0 !! sum t(t-1)(t-2)
        integer(int64), intent(out) :: cubic1 !! sum t(t-1)(2t+5)
        integer, intent(out) :: nclasses !! number of distinct values

        integer :: group_end
        integer :: group_start
        integer(int64) :: t
        real(dp), allocatable :: sorted(:)

        pair_ties = 0
        cubic0 = 0
        cubic1 = 0
        nclasses = 0
        if (size(values) == 0) return
        sorted = values
        call sort_real(sorted)
        group_start = 1
        do while (group_start <= size(sorted))
            group_end = group_start
            do while (group_end < size(sorted))
                if (.not. same_value(sorted(group_end), sorted(group_end + 1))) exit
                group_end = group_end + 1
            end do
            t = int(group_end - group_start + 1, int64)
            nclasses = nclasses + 1
            if (t > 1) then
                pair_ties = pair_ties + t * (t - 1_int64) / 2_int64
                cubic0 = cubic0 + t * (t - 1_int64) * (t - 2_int64)
                cubic1 = cubic1 + t * (t - 1_int64) * (2_int64 * t + 5_int64)
            end if
            group_start = group_end + 1
        end do
    end subroutine tie_statistics

    function tie_cubic5_sum(values) result(value)
        real(dp), intent(in) :: values(:) !! values without NaNs
        real(dp) :: value

        integer(int64) :: dummy0
        integer(int64) :: dummy_pairs
        integer :: dummy_classes
        integer(int64) :: cubic1

        call tie_statistics(values, dummy_pairs, dummy0, cubic1, dummy_classes)
        value = real(cubic1, dp)
    end function tie_cubic5_sum

    function page_exact_sf(l_observed, m, k) result(probability)
        integer, intent(in) :: l_observed !! observed L truncated toward zero as in SciPy
        integer, intent(in) :: m !! number of rows/subjects
        integer, intent(in) :: k !! number of treatments
        real(dp) :: probability

        integer :: a
        integer :: b
        integer :: i
        integer :: l
        integer :: max_sum
        integer :: min_sum
        integer :: row
        real(dp), allocatable :: current(:)
        real(dp), allocatable :: next(:)
        real(dp), allocatable :: single(:)

        probability = quiet_nan(0.0_dp)
        if (m < 1 .or. k < 1 .or. k > 10) return
        a = k * (k + 1) * (k + 2) / 6
        b = k * (k + 1) * (2 * k + 1) / 6
        allocate(single(a:b))
        single = 0.0_dp
        call page_single_row_pmf(k, a, b, single)

        min_sum = a
        max_sum = b
        allocate(current(min_sum:max_sum))
        current = single
        do row = 2, m
            allocate(next(row * a:row * b))
            next = 0.0_dp
            do i = min_sum, max_sum
                if (current(i) == 0.0_dp) cycle
                do l = a, b
                    if (single(l) == 0.0_dp) cycle
                    next(i + l) = next(i + l) + current(i) * single(l)
                end do
            end do
            call move_alloc(next, current)
            min_sum = row * a
            max_sum = row * b
        end do

        if (l_observed <= min_sum) then
            probability = 1.0_dp
        else if (l_observed > max_sum) then
            probability = 0.0_dp
        else
            probability = sum(current(l_observed:max_sum))
        end if
        probability = max(0.0_dp, min(1.0_dp, probability))
    end function page_exact_sf

    subroutine page_single_row_pmf(k, a, b, pmf)
        integer, intent(in) :: k !! number of treatments
        integer, intent(in) :: a !! minimum one-row L
        integer, intent(in) :: b !! maximum one-row L
        real(dp), intent(out) :: pmf(a:b) !! one-row null PMF

        integer :: perm(k)
        logical :: used(k)
        real(dp) :: fact

        pmf = 0.0_dp
        used = .false.
        call enumerate_page_permutations(1, k, a, perm, used, pmf)
        fact = factorial_real(k)
        pmf = pmf / fact
    end subroutine page_single_row_pmf

    recursive subroutine enumerate_page_permutations(position, k, a, perm, used, pmf)
        integer, intent(in) :: position !! current permutation position
        integer, intent(in) :: k !! permutation length
        integer, intent(in) :: a !! minimum one-row L used as PMF offset
        integer, intent(inout) :: perm(k) !! current permutation
        logical, intent(inout) :: used(k) !! value-used flags
        real(dp), intent(inout) :: pmf(:) !! offset-indexed by declared lower bound

        integer :: j
        integer :: lvalue

        if (position > k) then
            lvalue = 0
            do j = 1, k
                lvalue = lvalue + j * perm(j)
            end do
            pmf(lvalue - a + 1) = pmf(lvalue - a + 1) + 1.0_dp
            return
        end if
        do j = 1, k
            if (.not. used(j)) then
                used(j) = .true.
                perm(position) = j
                call enumerate_page_permutations(position + 1, k, a, perm, used, pmf)
                used(j) = .false.
            end if
        end do
    end subroutine enumerate_page_permutations

    function factorial_real(n) result(value)
        integer, intent(in) :: n !! nonnegative integer
        real(dp) :: value
        integer :: i

        value = 1.0_dp
        do i = 2, n
            value = value * real(i, dp)
        end do
    end function factorial_real

    function median_sorted(values) result(value)
        real(dp), intent(in) :: values(:) !! nonempty sorted values
        real(dp) :: value
        integer :: n

        n = size(values)
        if (n <= 0) then
            value = quiet_nan(0.0_dp)
        else if (mod(n, 2) == 1) then
            value = values((n + 1) / 2)
        else
            value = 0.5_dp * (values(n / 2) + values(n / 2 + 1))
        end if
    end function median_sorted

    subroutine sort_real(values)
        real(dp), intent(inout) :: values(:) !! sorted ascending in place
        real(dp), allocatable :: workspace(:)

        if (size(values) <= 1) return
        allocate(workspace(size(values)))
        call merge_sort_real(values, workspace, 1, size(values))
    end subroutine sort_real

    recursive subroutine merge_sort_real(values, workspace, left, right)
        real(dp), intent(inout) :: values(:) !! values being sorted
        real(dp), intent(inout) :: workspace(:) !! merge workspace
        integer, intent(in) :: left !! first index
        integer, intent(in) :: right !! last index

        integer :: i
        integer :: j
        integer :: k
        integer :: mid

        if (left >= right) return
        mid = left + (right - left) / 2
        call merge_sort_real(values, workspace, left, mid)
        call merge_sort_real(values, workspace, mid + 1, right)
        i = left
        j = mid + 1
        do k = left, right
            if (i > mid) then
                workspace(k) = values(j)
                j = j + 1
            else if (j > right) then
                workspace(k) = values(i)
                i = i + 1
            else if (values(i) <= values(j)) then
                workspace(k) = values(i)
                i = i + 1
            else
                workspace(k) = values(j)
                j = j + 1
            end if
        end do
        values(left:right) = workspace(left:right)
    end subroutine merge_sort_real

    integer function round_nearest_even(value) result(rounded)
        real(dp), intent(in) :: value !! finite nonnegative value
        integer :: lower
        real(dp) :: fraction

        lower = int(floor(value))
        fraction = value - real(lower, dp)
        if (fraction < 0.5_dp) then
            rounded = lower
        else if (fraction > 0.5_dp) then
            rounded = lower + 1
        else if (mod(lower, 2) == 0) then
            rounded = lower
        else
            rounded = lower + 1
        end if
    end function round_nearest_even

    pure logical function valid_permutation(values) result(valid)
        integer, intent(in) :: values(:) !! candidate permutation
        logical, allocatable :: seen(:)
        integer :: i

        valid = .false.
        if (size(values) == 0) return
        allocate(seen(size(values)))
        seen = .false.
        do i = 1, size(values)
            if (values(i) < 1 .or. values(i) > size(values)) return
            if (seen(values(i))) return
            seen(values(i)) = .true.
        end do
        valid = .true.
    end function valid_permutation

    pure logical function same_value(x, y) result(same)
        real(dp), intent(in) :: x !! first non-NaN value
        real(dp), intent(in) :: y !! second non-NaN value

        same = x <= y .and. y <= x
    end function same_value

    pure function normal_pvalue(statistic, alternative) result(pvalue)
        real(dp), intent(in) :: statistic !! standard-normal test statistic
        character(len=*), intent(in) :: alternative !! valid alternative
        real(dp) :: pvalue

        select case (alternative)
        case ('less')
            pvalue = normal_cdf(statistic)
        case ('greater')
            pvalue = normal_sf(statistic)
        case default
            pvalue = min(1.0_dp, 2.0_dp * normal_sf(abs(statistic)))
        end select
    end function normal_pvalue

    pure function student_t_pvalue(statistic, df, alternative) result(pvalue)
        real(dp), intent(in) :: statistic !! Student t statistic
        real(dp), intent(in) :: df !! positive degrees of freedom
        character(len=*), intent(in) :: alternative !! valid alternative
        real(dp) :: pvalue

        select case (alternative)
        case ('less')
            pvalue = t_cdf(statistic, df)
        case ('greater')
            pvalue = t_sf(statistic, df)
        case default
            pvalue = min(1.0_dp, 2.0_dp * t_sf(abs(statistic), df))
        end select
    end function student_t_pvalue

    pure logical function valid_alternative(alternative) result(valid)
        character(len=*), intent(in) :: alternative !! alternative name

        valid = alternative == 'two-sided' .or. alternative == 'less' .or. &
            alternative == 'greater'
    end function valid_alternative

    subroutine invalidate_association(result_value)
        type(association_result), intent(out) :: result_value !! result set to NaN

        result_value%statistic = quiet_nan(0.0_dp)
        result_value%pvalue = quiet_nan(0.0_dp)
    end subroutine invalidate_association

    subroutine invalidate_linregress(result_value)
        type(linregress_result), intent(out) :: result_value !! result set to NaN

        result_value%slope = quiet_nan(0.0_dp)
        result_value%intercept = quiet_nan(0.0_dp)
        result_value%rvalue = quiet_nan(0.0_dp)
        result_value%pvalue = quiet_nan(0.0_dp)
        result_value%stderr = quiet_nan(0.0_dp)
        result_value%intercept_stderr = quiet_nan(0.0_dp)
    end subroutine invalidate_linregress

    subroutine invalidate_theil(result_value)
        type(theilslopes_result), intent(out) :: result_value !! result set to NaN

        result_value%slope = quiet_nan(0.0_dp)
        result_value%intercept = quiet_nan(0.0_dp)
        result_value%low_slope = quiet_nan(0.0_dp)
        result_value%high_slope = quiet_nan(0.0_dp)
    end subroutine invalidate_theil

    subroutine invalidate_siegel(result_value)
        type(siegelslopes_result), intent(out) :: result_value !! result set to NaN

        result_value%slope = quiet_nan(0.0_dp)
        result_value%intercept = quiet_nan(0.0_dp)
    end subroutine invalidate_siegel

    subroutine invalidate_page(result_value)
        type(page_trend_result), intent(out) :: result_value !! result set to NaN

        result_value%statistic = quiet_nan(0.0_dp)
        result_value%pvalue = quiet_nan(0.0_dp)
        result_value%method = ''
    end subroutine invalidate_page

end module scifort_association_extended
