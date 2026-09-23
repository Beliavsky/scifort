! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Additional nonparametric and multi-sample tests modeled after scipy.stats
! 1.17.0.  The formulas are independently implemented from the published
! definitions and checked numerically against SciPy reference values.

module scifort_nonparametric_extended
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_chi2, only : chi2_sf
    use scifort_descriptive, only : median, rankdata
    use scifort_goodness_of_fit, only : significance_result
    use scifort_hypothesis_extended, only : chi2_contingency, chi2_contingency_result, sample_group
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, symmetric_eigen_jacobi
    use scifort_math, only : quiet_nan
    use scifort_normal, only : normal_cdf, normal_sf
    implicit none
    private

    type, public :: anderson_ksamp_result
        real(dp) :: statistic !! normalized k-sample Anderson-Darling statistic
        real(dp) :: pvalue !! interpolated/capped/floored p-value
        real(dp), allocatable :: critical_values(:) !! critical values at standard significance levels
    end type anderson_ksamp_result

    type, public :: median_test_result
        real(dp) :: statistic !! contingency-table power-divergence statistic
        real(dp) :: pvalue !! chi-square reference p-value
        real(dp) :: median !! grand median of all observations
        integer, allocatable :: table(:, :) !! 2 by k table: above then below median
    end type median_test_result

    public :: anderson_ksamp
    public :: ansari
    public :: epps_singleton_2samp
    public :: median_test
    public :: median_test_named
    public :: mood

contains

    function anderson_ksamp(groups, variant) result(result_value)
        type(sample_group), intent(in) :: groups(:) !! two or more independent samples
        character(len=*), intent(in), optional :: variant !! midrank, right, or continuous
        type(anderson_ksamp_result) :: result_value

        real(dp), parameter :: b0(7) = [0.675_dp, 1.281_dp, 1.645_dp, 1.96_dp, &
            2.326_dp, 2.573_dp, 3.085_dp]
        real(dp), parameter :: b1(7) = [-0.245_dp, 0.25_dp, 0.678_dp, 1.149_dp, &
            1.822_dp, 2.364_dp, 3.615_dp]
        real(dp), parameter :: b2(7) = [-0.105_dp, -0.305_dp, -0.362_dp, -0.391_dp, &
            -0.396_dp, -0.345_dp, -0.154_dp]
        real(dp), parameter :: sig(7) = [0.25_dp, 0.1_dp, 0.05_dp, 0.025_dp, &
            0.01_dp, 0.005_dp, 0.001_dp]

        real(dp) :: acoef
        real(dp) :: a2
        real(dp) :: a2kn
        real(dp) :: bcoef
        real(dp) :: ccoef
        character(len=10) :: choice
        real(dp) :: cum
        real(dp) :: dcoef
        real(dp) :: g
        real(dp) :: h
        real(dp) :: hsum
        integer :: i
        integer :: k
        integer :: n_total
        integer, allocatable :: sizes(:)
        real(dp) :: sigma_sq
        real(dp), allocatable :: pooled(:)
        real(dp), allocatable :: unique_values(:)

        allocate(result_value%critical_values(7))
        result_value%critical_values = quiet_nan(0.0_dp)
        call invalidate_pair(result_value%statistic, result_value%pvalue)

        k = size(groups)
        if (k < 2) return
        allocate(sizes(k))
        do i = 1, k
            if (.not. allocated(groups(i)%values)) return
            sizes(i) = size(groups(i)%values)
            if (sizes(i) < 1) return
            if (any(.not. ieee_is_finite(groups(i)%values))) return
        end do
        n_total = sum(sizes)
        if (n_total < 4) return

        call pool_groups(groups, pooled)
        call sort_values(pooled)
        unique_values = unique_sorted(pooled)
        if (size(unique_values) < 2) return

        choice = 'midrank'
        if (present(variant)) choice = trim(variant)
        select case (trim(choice))
        case ('midrank')
            a2kn = anderson_midrank_statistic(groups, pooled, unique_values, sizes)
        case ('right')
            a2kn = anderson_right_statistic(groups, pooled, unique_values, sizes)
        case ('continuous')
            a2kn = anderson_continuous_statistic(groups, pooled, sizes)
        case default
            return
        end select
        if (.not. ieee_is_finite(a2kn)) return

        hsum = sum(1.0_dp / real(sizes, dp))
        cum = 0.0_dp
        g = 0.0_dp
        do i = 1, n_total - 2
            cum = cum + 1.0_dp / real(n_total - i, dp)
            g = g + cum / real(i + 1, dp)
        end do
        h = cum + 1.0_dp

        acoef = (4.0_dp * g - 6.0_dp) * real(k - 1, dp) + (10.0_dp - 6.0_dp * g) * hsum
        bcoef = (2.0_dp * g - 4.0_dp) * real(k * k, dp) + 8.0_dp * h * real(k, dp) + &
            (2.0_dp * g - 14.0_dp * h - 4.0_dp) * hsum - 8.0_dp * h + 4.0_dp * g - 6.0_dp
        ccoef = (6.0_dp * h + 2.0_dp * g - 2.0_dp) * real(k * k, dp) + &
            (4.0_dp * h - 4.0_dp * g + 6.0_dp) * real(k, dp) + &
            (2.0_dp * h - 6.0_dp) * hsum + 4.0_dp * h
        dcoef = (2.0_dp * h + 6.0_dp) * real(k * k, dp) - 4.0_dp * h * real(k, dp)
        sigma_sq = (acoef * real(n_total, dp) ** 3 + bcoef * real(n_total, dp) ** 2 + &
            ccoef * real(n_total, dp) + dcoef) / &
            real((n_total - 1) * (n_total - 2) * (n_total - 3), dp)
        if (.not. ieee_is_finite(sigma_sq) .or. sigma_sq <= 0.0_dp) return

        a2 = (a2kn - real(k - 1, dp)) / sqrt(sigma_sq)
        result_value%critical_values = b0 + b1 / sqrt(real(k - 1, dp)) + b2 / real(k - 1, dp)
        result_value%statistic = a2
        if (a2 < minval(result_value%critical_values)) then
            result_value%pvalue = maxval(sig)
        else if (a2 > maxval(result_value%critical_values)) then
            result_value%pvalue = minval(sig)
        else
            result_value%pvalue = exp(quadratic_least_squares_eval(&
                result_value%critical_values, log(sig), a2))
        end if
    end function anderson_ksamp

    function ansari(x, y, alternative) result(result_value)
        real(dp), intent(in) :: x(:) !! first independent sample
        real(dp), intent(in) :: y(:) !! second independent sample
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        type(significance_result) :: result_value

        real(dp) :: ab
        character(len=9) :: alt
        real(dp) :: fac
        integer :: m
        real(dp) :: mn_ab
        integer :: n
        integer :: n_total
        real(dp), allocatable :: pooled(:)
        real(dp), allocatable :: ranks(:)
        logical :: repeats
        real(dp), allocatable :: symrank(:)
        real(dp) :: var_ab
        real(dp) :: z

        call invalidate_pair(result_value%statistic, result_value%pvalue)
        alt = select_alternative(alternative)
        if (.not. valid_alternative(alt)) return
        n = size(x)
        m = size(y)
        if (n < 1 .or. m < 1) return
        if (any(.not. ieee_is_finite(x)) .or. any(.not. ieee_is_finite(y))) return
        n_total = n + m
        allocate(pooled(n_total))
        pooled(1:n) = x
        pooled(n + 1:) = y
        ranks = rankdata(pooled)
        allocate(symrank(n_total))
        symrank = min(ranks, real(n_total + 1, dp) - ranks)
        ab = sum(symrank(1:n))
        repeats = has_duplicate(pooled)
        result_value%statistic = ab

        if (n < 55 .and. m < 55 .and. .not. repeats) then
            result_value%pvalue = ansari_exact_pvalue(n, m, nint(ab), alt)
            return
        end if

        if (mod(n_total, 2) == 1) then
            mn_ab = real(n, dp) * real(n_total + 1, dp) ** 2 / (4.0_dp * real(n_total, dp))
        else
            mn_ab = real(n * (n_total + 2), dp) / 4.0_dp
        end if

        if (repeats) then
            fac = sum(symrank ** 2)
            if (mod(n_total, 2) == 1) then
                var_ab = real(m * n, dp) * &
                    (16.0_dp * real(n_total, dp) * fac - real(n_total + 1, dp) ** 4) / &
                    (16.0_dp * real(n_total * n_total * (n_total - 1), dp))
            else
                var_ab = real(m * n, dp) * &
                    (16.0_dp * fac - real(n_total * (n_total + 2) ** 2, dp)) / &
                    (16.0_dp * real(n_total * (n_total - 1), dp))
            end if
        else if (mod(n_total, 2) == 1) then
            var_ab = real(n * m * (n_total + 1) * (3 + n_total * n_total), dp) / &
                (48.0_dp * real(n_total * n_total, dp))
        else
            var_ab = real(m * n * (n_total + 2), dp) * real(n_total - 2, dp) / &
                (48.0_dp * real(n_total - 1, dp))
        end if
        if (.not. ieee_is_finite(var_ab) .or. var_ab <= 0.0_dp) then
            result_value%pvalue = quiet_nan(0.0_dp)
            return
        end if
        z = (mn_ab - ab) / sqrt(var_ab)
        result_value%pvalue = normal_pvalue(z, alt)
    end function ansari

    function mood(x, y, alternative) result(result_value)
        real(dp), intent(in) :: x(:) !! first independent sample
        real(dp), intent(in) :: y(:) !! second independent sample
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        type(significance_result) :: result_value

        character(len=9) :: alt
        integer :: m
        integer :: n
        integer :: n_total
        real(dp), allocatable :: pooled(:)
        real(dp), allocatable :: ranks(:)
        real(dp) :: z

        call invalidate_pair(result_value%statistic, result_value%pvalue)
        alt = select_alternative(alternative)
        if (.not. valid_alternative(alt)) return
        m = size(x)
        n = size(y)
        n_total = m + n
        if (m < 1 .or. n < 1 .or. n_total < 3) return
        if (any(.not. ieee_is_finite(x)) .or. any(.not. ieee_is_finite(y))) return
        allocate(pooled(n_total))
        pooled(1:m) = x
        pooled(m + 1:) = y
        ranks = rankdata(pooled)
        if (has_duplicate(pooled)) then
            z = mood_tied_statistic(x, y)
        else
            z = mood_untied_statistic(ranks, m, n)
        end if
        result_value%statistic = z
        if (.not. ieee_is_finite(z)) then
            result_value%pvalue = quiet_nan(0.0_dp)
        else
            result_value%pvalue = normal_pvalue(z, alt)
        end if
    end function mood

    function epps_singleton_2samp(x, y, t) result(result_value)
        real(dp), intent(in) :: x(:) !! first independent sample, at least five observations
        real(dp), intent(in) :: y(:) !! second independent sample, at least five observations
        real(dp), intent(in), optional :: t(:) !! positive distinct characteristic-function points
        type(significance_result) :: result_value

        real(dp), parameter :: default_t(2) = [0.4_dp, 0.8_dp]
        real(dp), allocatable :: cov_x(:, :)
        real(dp), allocatable :: cov_y(:, :)
        real(dp), allocatable :: eigvals(:)
        real(dp), allocatable :: eigvecs(:, :)
        real(dp), allocatable :: gx(:, :)
        real(dp), allocatable :: gy(:, :)
        real(dp), allocatable :: gdiff(:)
        integer :: i
        integer :: j
        integer :: nx
        integer :: ny
        integer :: n_total
        integer :: p
        integer :: rank_cov
        real(dp) :: corr
        real(dp) :: max_eig
        real(dp), allocatable :: pooled(:)
        real(dp) :: projection
        real(dp) :: sigma
        integer :: status
        real(dp), allocatable :: tvals(:)
        real(dp) :: tol
        real(dp) :: w

        call invalidate_pair(result_value%statistic, result_value%pvalue)
        nx = size(x)
        ny = size(y)
        if (nx < 5 .or. ny < 5) return
        if (any(.not. ieee_is_finite(x)) .or. any(.not. ieee_is_finite(y))) return
        if (present(t)) then
            if (size(t) < 1 .or. any(.not. ieee_is_finite(t)) .or. any(t <= 0.0_dp)) return
            if (has_duplicate(t)) return
            tvals = t
        else
            tvals = default_t
        end if

        n_total = nx + ny
        allocate(pooled(n_total))
        pooled(1:nx) = x
        pooled(nx + 1:) = y
        call sort_values(pooled)
        sigma = 0.5_dp * (linear_quantile_sorted(pooled, 0.75_dp) - &
            linear_quantile_sorted(pooled, 0.25_dp))
        if (.not. ieee_is_finite(sigma) .or. sigma <= 0.0_dp) return

        p = 2 * size(tvals)
        allocate(gx(p, nx), gy(p, ny))
        do j = 1, size(tvals)
            do i = 1, nx
                gx(j, i) = cos(tvals(j) * x(i) / sigma)
                gx(size(tvals) + j, i) = sin(tvals(j) * x(i) / sigma)
            end do
            do i = 1, ny
                gy(j, i) = cos(tvals(j) * y(i) / sigma)
                gy(size(tvals) + j, i) = sin(tvals(j) * y(i) / sigma)
            end do
        end do

        cov_x = biased_covariance(gx)
        cov_y = biased_covariance(gy)
        cov_x = real(n_total, dp) / real(nx, dp) * cov_x + &
            real(n_total, dp) / real(ny, dp) * cov_y
        allocate(eigvals(p), eigvecs(p, p), gdiff(p))
        call symmetric_eigen_jacobi(cov_x, eigvals, eigvecs, status)
        if (status /= linalg_status_success) return
        max_eig = max(0.0_dp, maxval(eigvals))
        tol = real(p, dp) * epsilon(1.0_dp) * max_eig
        rank_cov = count(eigvals > tol)
        if (rank_cov < 1) return

        gdiff = sum(gx, dim=2) / real(nx, dp) - sum(gy, dim=2) / real(ny, dp)
        w = 0.0_dp
        do j = 1, p
            if (eigvals(j) > tol) then
                projection = dot_product(eigvecs(:, j), gdiff)
                w = w + projection * projection / eigvals(j)
            end if
        end do
        w = real(n_total, dp) * w
        if (max(nx, ny) < 25) then
            corr = 1.0_dp / (1.0_dp + real(n_total, dp) ** (-0.45_dp) + &
                10.1_dp * (real(nx, dp) ** (-1.7_dp) + real(ny, dp) ** (-1.7_dp)))
            w = w * corr
        end if
        result_value%statistic = w
        result_value%pvalue = chi2_sf(w, real(rank_cov, dp))
    end function epps_singleton_2samp

    function median_test(groups, ties, correction, lambda_) result(result_value)
        type(sample_group), intent(in) :: groups(:) !! two or more independent samples
        character(len=*), intent(in), optional :: ties !! below, above, or ignore
        logical, intent(in), optional :: correction !! Yates correction for two samples
        real(dp), intent(in), optional :: lambda_ !! Cressie-Read lambda, default one
        type(median_test_result) :: result_value

        real(dp) :: lambda_value

        lambda_value = 1.0_dp
        if (present(lambda_)) lambda_value = lambda_
        result_value = median_test_core(groups, ties, correction, lambda_value)
    end function median_test

    function median_test_named(groups, ties, correction, lambda_name) result(result_value)
        type(sample_group), intent(in) :: groups(:) !! two or more independent samples
        character(len=*), intent(in), optional :: ties !! below, above, or ignore
        logical, intent(in), optional :: correction !! Yates correction for two samples
        character(len=*), intent(in) :: lambda_name !! SciPy Cressie-Read alias
        type(median_test_result) :: result_value

        real(dp) :: lambda_value

        select case (trim(lambda_name))
        case ('pearson')
            lambda_value = 1.0_dp
        case ('log-likelihood')
            lambda_value = 0.0_dp
        case ('freeman-tukey')
            lambda_value = -0.5_dp
        case ('mod-log-likelihood')
            lambda_value = -1.0_dp
        case ('neyman')
            lambda_value = -2.0_dp
        case ('cressie-read')
            lambda_value = 2.0_dp / 3.0_dp
        case default
            call invalidate_median_result(result_value, size(groups))
            return
        end select
        result_value = median_test_core(groups, ties, correction, lambda_value)
    end function median_test_named

    function median_test_core(groups, ties, correction, lambda_value) result(result_value)
        type(sample_group), intent(in) :: groups(:) !! two or more independent samples
        character(len=*), intent(in), optional :: ties !! below, above, or ignore
        logical, intent(in), optional :: correction !! Yates correction for two samples
        real(dp), intent(in) :: lambda_value !! validated Cressie-Read lambda
        type(median_test_result) :: result_value

        integer :: above
        integer :: below
        character(len=6) :: choice
        type(chi2_contingency_result) :: contingency
        logical :: correct
        integer :: equal
        integer :: i
        integer :: k
        real(dp), allocatable :: observed(:, :)
        real(dp), allocatable :: pooled(:)

        k = size(groups)
        call invalidate_median_result(result_value, k)
        if (k < 2) return
        choice = 'below'
        if (present(ties)) choice = trim(ties)
        if (choice /= 'below' .and. choice /= 'above' .and. choice /= 'ignore') return
        do i = 1, k
            if (.not. allocated(groups(i)%values)) return
            if (size(groups(i)%values) < 1) return
            if (any(.not. ieee_is_finite(groups(i)%values))) return
        end do
        call pool_groups(groups, pooled)
        result_value%median = median(pooled)
        result_value%table = 0
        do i = 1, k
            above = count(groups(i)%values > result_value%median)
            below = count(groups(i)%values < result_value%median)
            equal = size(groups(i)%values) - above - below
            result_value%table(1, i) = above
            result_value%table(2, i) = below
            if (choice == 'below') result_value%table(2, i) = result_value%table(2, i) + equal
            if (choice == 'above') result_value%table(1, i) = result_value%table(1, i) + equal
        end do
        if (sum(result_value%table(1, :)) == 0 .or. sum(result_value%table(2, :)) == 0) return
        if (choice == 'ignore') then
            do i = 1, k
                if (sum(result_value%table(:, i)) == 0) return
            end do
        end if

        allocate(observed(2, k))
        observed = real(result_value%table, dp)
        correct = .true.
        if (present(correction)) correct = correction
        contingency = chi2_contingency(observed, correction=correct, lambda_=lambda_value)
        result_value%statistic = contingency%statistic
        result_value%pvalue = contingency%pvalue
    end function median_test_core

    function anderson_midrank_statistic(groups, pooled, unique_values, sizes) result(statistic)
        type(sample_group), intent(in) :: groups(:) !! validated groups
        real(dp), intent(in) :: pooled(:) !! sorted pooled observations
        real(dp), intent(in) :: unique_values(:) !! sorted unique pooled values
        integer, intent(in) :: sizes(:) !! sample sizes
        real(dp) :: statistic

        real(dp) :: bj
        real(dp) :: denominator
        real(dp) :: fij
        integer :: i
        integer :: j
        real(dp) :: lj
        real(dp) :: mij
        integer :: n_total

        n_total = size(pooled)
        statistic = 0.0_dp
        do j = 1, size(unique_values)
            lj = real(count_equal_sorted(pooled, unique_values(j)), dp)
            bj = real(count_less_sorted(pooled, unique_values(j)), dp) + 0.5_dp * lj
            denominator = bj * (real(n_total, dp) - bj) - real(n_total, dp) * lj / 4.0_dp
            if (denominator <= 0.0_dp) cycle
            do i = 1, size(groups)
                fij = real(count(groups(i)%values == unique_values(j)), dp)
                mij = real(count(groups(i)%values <= unique_values(j)), dp) - 0.5_dp * fij
                statistic = statistic + lj / real(n_total, dp) * &
                    (real(n_total, dp) * mij - bj * real(sizes(i), dp)) ** 2 / &
                    denominator / real(sizes(i), dp)
            end do
        end do
        statistic = statistic * real(n_total - 1, dp) / real(n_total, dp)
    end function anderson_midrank_statistic

    function anderson_right_statistic(groups, pooled, unique_values, sizes) result(statistic)
        type(sample_group), intent(in) :: groups(:) !! validated groups
        real(dp), intent(in) :: pooled(:) !! sorted pooled observations
        real(dp), intent(in) :: unique_values(:) !! sorted unique pooled values
        integer, intent(in) :: sizes(:) !! sample sizes
        real(dp) :: statistic

        real(dp) :: bj
        integer :: i
        integer :: j
        real(dp) :: lj
        real(dp) :: mij
        integer :: n_total

        n_total = size(pooled)
        statistic = 0.0_dp
        bj = 0.0_dp
        do j = 1, size(unique_values) - 1
            lj = real(count_equal_sorted(pooled, unique_values(j)), dp)
            bj = bj + lj
            do i = 1, size(groups)
                mij = real(count(groups(i)%values <= unique_values(j)), dp)
                statistic = statistic + lj / real(n_total, dp) * &
                    (real(n_total, dp) * mij - bj * real(sizes(i), dp)) ** 2 / &
                    (bj * (real(n_total, dp) - bj)) / real(sizes(i), dp)
            end do
        end do
    end function anderson_right_statistic

    function anderson_continuous_statistic(groups, pooled, sizes) result(statistic)
        type(sample_group), intent(in) :: groups(:) !! validated groups
        real(dp), intent(in) :: pooled(:) !! sorted pooled observations
        integer, intent(in) :: sizes(:) !! sample sizes
        real(dp) :: statistic

        integer :: i
        integer :: j
        real(dp) :: mij
        integer :: n_total

        n_total = size(pooled)
        statistic = 0.0_dp
        do j = 1, n_total - 1
            do i = 1, size(groups)
                mij = real(count(groups(i)%values <= pooled(j)), dp)
                statistic = statistic + &
                    (real(n_total, dp) * mij - real(j * sizes(i), dp)) ** 2 / &
                    (real(j * (n_total - j) * sizes(i), dp) * real(n_total, dp))
            end do
        end do
    end function anderson_continuous_statistic

    function ansari_exact_pvalue(n, m, observed, alternative) result(pvalue)
        integer, intent(in) :: n !! size of first sample
        integer, intent(in) :: m !! size of second sample
        integer, intent(in) :: observed !! observed integer Ansari-Bradley statistic
        character(len=*), intent(in) :: alternative !! validated alternative
        real(dp) :: pvalue

        real(dp) :: cdf
        real(dp), allocatable :: counts(:, :)
        integer :: i
        integer :: k
        integer :: max_sum
        integer :: n_total
        integer :: s
        integer :: score
        real(dp) :: sf
        real(dp) :: total

        n_total = n + m
        max_sum = n * ((n_total + 1) / 2)
        allocate(counts(0:n, 0:max_sum))
        counts = 0.0_dp
        counts(0, 0) = 1.0_dp
        do i = 1, n_total
            score = min(i, n_total - i + 1)
            do k = min(i, n), 1, -1
                do s = max_sum, score, -1
                    counts(k, s) = counts(k, s) + counts(k - 1, s - score)
                end do
            end do
        end do
        total = sum(counts(n, :))
        if (total <= 0.0_dp) then
            pvalue = quiet_nan(0.0_dp)
            return
        end if
        cdf = sum(counts(n, 0:min(max_sum, observed))) / total
        sf = sum(counts(n, max(0, observed):max_sum)) / total
        select case (alternative)
        case ('greater')
            pvalue = cdf
        case ('less')
            pvalue = sf
        case default
            pvalue = min(1.0_dp, 2.0_dp * min(cdf, sf))
        end select
    end function ansari_exact_pvalue

    function mood_untied_statistic(ranks, m, n) result(z)
        real(dp), intent(in) :: ranks(:) !! pooled ranks with no ties
        integer, intent(in) :: m !! first sample size
        integer, intent(in) :: n !! second sample size
        real(dp) :: z

        real(dp) :: expected
        real(dp) :: moment
        integer :: n_total
        real(dp) :: variance_m

        n_total = m + n
        moment = sum((ranks(1:m) - real(n_total + 1, dp) / 2.0_dp) ** 2)
        expected = real(m * (n_total * n_total - 1), dp) / 12.0_dp
        variance_m = real(m * n * (n_total + 1) * (n_total + 2) * (n_total - 2), dp) / 180.0_dp
        if (variance_m <= 0.0_dp) then
            z = quiet_nan(0.0_dp)
        else
            z = (moment - expected) / sqrt(variance_m)
        end if
    end function mood_untied_statistic

    function mood_tied_statistic(x, y) result(z)
        real(dp), intent(in) :: x(:) !! first independent sample
        real(dp), intent(in) :: y(:) !! second independent sample
        real(dp) :: z

        integer :: a_count
        real(dp) :: c
        integer :: i
        integer :: j
        integer :: m
        integer :: n
        integer :: n_total
        real(dp), allocatable :: pooled(:)
        integer :: s_i
        integer :: s_prev
        real(dp) :: sum_i
        real(dp) :: sum_i2
        real(dp) :: term
        integer :: tie_count
        real(dp) :: tstat
        real(dp), allocatable :: unique_values(:)
        real(dp) :: variance_m
        real(dp) :: expected

        m = size(x)
        n = size(y)
        n_total = m + n
        allocate(pooled(n_total))
        pooled(1:m) = x
        pooled(m + 1:) = y
        call sort_values(pooled)
        unique_values = unique_sorted(pooled)
        c = real(n_total + 1, dp) / 2.0_dp
        s_prev = 0
        tstat = 0.0_dp
        variance_m = real(m * n * (n_total + 1) * (n_total * n_total - 4), dp) / 180.0_dp
        do j = 1, size(unique_values)
            tie_count = count_equal_sorted(pooled, unique_values(j))
            s_i = s_prev + tie_count
            sum_i = 0.0_dp
            sum_i2 = 0.0_dp
            do i = s_prev + 1, s_i
                sum_i = sum_i + real(i, dp)
                sum_i2 = sum_i2 + real(i * i, dp)
            end do
            a_count = count(x == unique_values(j))
            tstat = tstat + real(a_count, dp) * &
                (sum_i2 - 2.0_dp * c * sum_i + real(tie_count, dp) * c * c) / &
                real(tie_count, dp)
            term = real(tie_count * (tie_count * tie_count - 1), dp) * &
                (real(tie_count * tie_count - 4, dp) + &
                15.0_dp * real(n_total - s_i - s_prev, dp) ** 2)
            variance_m = variance_m - real(m * n, dp) * term / &
                (180.0_dp * real(n_total * (n_total - 1), dp))
            s_prev = s_i
        end do
        expected = real(m * (n_total * n_total - 1), dp) / 12.0_dp
        if (variance_m <= 0.0_dp) then
            z = quiet_nan(0.0_dp)
        else
            z = (tstat - expected) / sqrt(variance_m)
        end if
    end function mood_tied_statistic

    function biased_covariance(values) result(covariance)
        real(dp), intent(in) :: values(:, :) !! variables by observations
        real(dp), allocatable :: covariance(:, :)

        integer :: i
        integer :: j
        integer :: nobs
        real(dp), allocatable :: means(:)

        nobs = size(values, 2)
        allocate(covariance(size(values, 1), size(values, 1)))
        allocate(means(size(values, 1)))
        means = sum(values, dim=2) / real(nobs, dp)
        do j = 1, size(values, 1)
            do i = 1, size(values, 1)
                covariance(i, j) = dot_product(values(i, :) - means(i), values(j, :) - means(j)) / &
                    real(nobs, dp)
            end do
        end do
    end function biased_covariance

    function quadratic_least_squares_eval(x, y, x0) result(value)
        real(dp), intent(in) :: x(:) !! abscissae
        real(dp), intent(in) :: y(:) !! ordinates
        real(dp), intent(in) :: x0 !! evaluation point
        real(dp) :: value

        real(dp) :: a(3, 3)
        real(dp) :: b(3)
        real(dp) :: coef(3)
        integer :: i
        integer :: status

        a = 0.0_dp
        b = 0.0_dp
        do i = 1, size(x)
            a(1, 1) = a(1, 1) + 1.0_dp
            a(1, 2) = a(1, 2) + x(i)
            a(1, 3) = a(1, 3) + x(i) * x(i)
            a(2, 2) = a(2, 2) + x(i) * x(i)
            a(2, 3) = a(2, 3) + x(i) ** 3
            a(3, 3) = a(3, 3) + x(i) ** 4
            b(1) = b(1) + y(i)
            b(2) = b(2) + x(i) * y(i)
            b(3) = b(3) + x(i) * x(i) * y(i)
        end do
        a(2, 1) = a(1, 2)
        a(3, 1) = a(1, 3)
        a(3, 2) = a(2, 3)
        call solve_3x3(a, b, coef, status)
        if (status /= 0) then
            value = quiet_nan(0.0_dp)
        else
            value = coef(1) + coef(2) * x0 + coef(3) * x0 * x0
        end if
    end function quadratic_least_squares_eval

    subroutine solve_3x3(a, b, x, status)
        real(dp), intent(in) :: a(3, 3) !! coefficient matrix
        real(dp), intent(in) :: b(3) !! right-hand side
        real(dp), intent(out) :: x(3) !! solution
        integer, intent(out) :: status !! zero on success

        real(dp) :: aug(3, 4)
        real(dp) :: factor
        integer :: i
        integer :: j
        integer :: pivot
        real(dp) :: row(4)

        aug(:, 1:3) = a
        aug(:, 4) = b
        status = 0
        do i = 1, 3
            pivot = i
            do j = i + 1, 3
                if (abs(aug(j, i)) > abs(aug(pivot, i))) pivot = j
            end do
            if (abs(aug(pivot, i)) <= tiny(1.0_dp)) then
                status = 1
                x = quiet_nan(0.0_dp)
                return
            end if
            if (pivot /= i) then
                row = aug(i, :)
                aug(i, :) = aug(pivot, :)
                aug(pivot, :) = row
            end if
            aug(i, :) = aug(i, :) / aug(i, i)
            do j = 1, 3
                if (j == i) cycle
                factor = aug(j, i)
                aug(j, :) = aug(j, :) - factor * aug(i, :)
            end do
        end do
        x = aug(:, 4)
    end subroutine solve_3x3

    function normal_pvalue(z, alternative) result(pvalue)
        real(dp), intent(in) :: z !! standard-normal statistic
        character(len=*), intent(in) :: alternative !! validated alternative
        real(dp) :: pvalue

        select case (alternative)
        case ('less')
            pvalue = normal_cdf(z)
        case ('greater')
            pvalue = normal_sf(z)
        case default
            pvalue = min(1.0_dp, 2.0_dp * normal_sf(abs(z)))
        end select
    end function normal_pvalue

    function select_alternative(alternative) result(choice)
        character(len=*), intent(in), optional :: alternative !! optional alternative name
        character(len=9) :: choice

        choice = 'two-sided'
        if (present(alternative)) choice = trim(alternative)
    end function select_alternative

    pure function valid_alternative(alternative) result(valid)
        character(len=*), intent(in) :: alternative !! alternative name
        logical :: valid

        valid = alternative == 'two-sided' .or. alternative == 'less' .or. alternative == 'greater'
    end function valid_alternative

    subroutine invalidate_pair(statistic, pvalue)
        real(dp), intent(out) :: statistic !! statistic set to NaN
        real(dp), intent(out) :: pvalue !! p-value set to NaN

        statistic = quiet_nan(0.0_dp)
        pvalue = quiet_nan(0.0_dp)
    end subroutine invalidate_pair

    subroutine invalidate_median_result(result_value, k)
        type(median_test_result), intent(out) :: result_value !! result initialized to invalid values
        integer, intent(in) :: k !! requested number of sample columns

        result_value%statistic = quiet_nan(0.0_dp)
        result_value%pvalue = quiet_nan(0.0_dp)
        result_value%median = quiet_nan(0.0_dp)
        allocate(result_value%table(2, max(0, k)))
        result_value%table = 0
    end subroutine invalidate_median_result

    subroutine pool_groups(groups, pooled)
        type(sample_group), intent(in) :: groups(:) !! allocated groups
        real(dp), allocatable, intent(out) :: pooled(:) !! concatenated observations

        integer :: first
        integer :: i
        integer :: last
        integer :: n_total

        n_total = 0
        do i = 1, size(groups)
            n_total = n_total + size(groups(i)%values)
        end do
        allocate(pooled(n_total))
        first = 1
        do i = 1, size(groups)
            last = first + size(groups(i)%values) - 1
            pooled(first:last) = groups(i)%values
            first = last + 1
        end do
    end subroutine pool_groups

    function unique_sorted(sorted) result(values)
        real(dp), intent(in) :: sorted(:) !! sorted observations
        real(dp), allocatable :: values(:)

        integer :: i
        integer :: n_unique

        if (size(sorted) == 0) then
            allocate(values(0))
            return
        end if
        n_unique = 1
        do i = 2, size(sorted)
            if (sorted(i) /= sorted(i - 1)) n_unique = n_unique + 1
        end do
        allocate(values(n_unique))
        values(1) = sorted(1)
        n_unique = 1
        do i = 2, size(sorted)
            if (sorted(i) /= sorted(i - 1)) then
                n_unique = n_unique + 1
                values(n_unique) = sorted(i)
            end if
        end do
    end function unique_sorted

    pure integer function count_equal_sorted(sorted, value) result(number)
        real(dp), intent(in) :: sorted(:) !! sorted observations
        real(dp), intent(in) :: value !! value to count

        number = count(sorted == value)
    end function count_equal_sorted

    pure integer function count_less_sorted(sorted, value) result(number)
        real(dp), intent(in) :: sorted(:) !! sorted observations
        real(dp), intent(in) :: value !! threshold

        number = count(sorted < value)
    end function count_less_sorted

    function has_duplicate(values) result(found)
        real(dp), intent(in) :: values(:) !! finite observations
        logical :: found

        integer :: i
        real(dp), allocatable :: work(:)

        found = .false.
        if (size(values) < 2) return
        work = values
        call sort_values(work)
        do i = 2, size(work)
            if (work(i) == work(i - 1)) then
                found = .true.
                return
            end if
        end do
    end function has_duplicate

    function linear_quantile_sorted(sorted, probability) result(value)
        real(dp), intent(in) :: sorted(:) !! nonempty sorted observations
        real(dp), intent(in) :: probability !! probability in [0, 1]
        real(dp) :: value

        real(dp) :: fraction
        integer :: lower
        real(dp) :: position
        integer :: upper

        position = probability * real(size(sorted) - 1, dp)
        lower = int(floor(position)) + 1
        upper = min(size(sorted), lower + 1)
        fraction = position - floor(position)
        value = (1.0_dp - fraction) * sorted(lower) + fraction * sorted(upper)
    end function linear_quantile_sorted

    subroutine sort_values(values)
        real(dp), intent(inout) :: values(:) !! values sorted ascending in place

        real(dp), allocatable :: temporary(:)

        if (size(values) <= 1) return
        allocate(temporary(size(values)))
        call merge_sort_values(values, temporary, 1, size(values))
    end subroutine sort_values

    recursive subroutine merge_sort_values(values, temporary, left, right)
        real(dp), intent(inout) :: values(:) !! values being sorted
        real(dp), intent(inout) :: temporary(:) !! merge workspace
        integer, intent(in) :: left !! first index of inclusive range
        integer, intent(in) :: right !! last index of inclusive range

        integer :: i
        integer :: j
        integer :: k
        integer :: middle

        if (left >= right) return
        middle = left + (right - left) / 2
        call merge_sort_values(values, temporary, left, middle)
        call merge_sort_values(values, temporary, middle + 1, right)
        i = left
        j = middle + 1
        do k = left, right
            if (i > middle) then
                temporary(k) = values(j)
                j = j + 1
            else if (j > right) then
                temporary(k) = values(i)
                i = i + 1
            else if (values(i) <= values(j)) then
                temporary(k) = values(i)
                i = i + 1
            else
                temporary(k) = values(j)
                j = j + 1
            end if
        end do
        values(left:right) = temporary(left:right)
    end subroutine merge_sort_values

end module scifort_nonparametric_extended
