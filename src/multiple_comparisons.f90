! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Multiple-comparison and heterogeneous-mean procedures corresponding to
! scipy.stats 1.17.0.  Tukey/Games-Howell reuse the SciFort studentized-range
! distribution, Dunnett reuses the multivariate Student-t box integrator, and
! the Poisson means E-test reuses the translated Poisson distribution.
module scifort_multiple_comparisons
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_chi2, only : chi2_sf
    use scifort_descriptive, only : mean, variance
    use scifort_hypothesis_extended, only : sample_group
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan
    use scifort_multivariate_t, only : multivariate_t_cdf
    use scifort_poisson, only : poisson_pmf, poisson_ppf
    use scifort_studentized_range, only : studentized_range_ppf, studentized_range_sf
    implicit none
    private

    type, public :: alexandergovern_result
        real(dp) :: statistic !! Alexander-Govern A statistic
        real(dp) :: pvalue !! upper-tail chi-square p-value
    end type alexandergovern_result

    type, public :: matrix_confidence_interval
        real(dp), allocatable :: low(:, :) !! lower pairwise confidence limits
        real(dp), allocatable :: high(:, :) !! upper pairwise confidence limits
    end type matrix_confidence_interval

    type, public :: vector_confidence_interval
        real(dp), allocatable :: low(:) !! lower simultaneous confidence limits
        real(dp), allocatable :: high(:) !! upper simultaneous confidence limits
    end type vector_confidence_interval

    type, public :: tukey_hsd_result
        real(dp), allocatable :: statistic(:, :) !! pairwise mean differences i-j
        real(dp), allocatable :: pvalue(:, :) !! studentized-range adjusted p-values
        real(dp), allocatable :: df(:, :) !! pairwise reference degrees of freedom
        real(dp), allocatable :: stand_err(:, :) !! SciPy studentized-range standard error
        integer :: ntreatments = 0 !! number of groups
    contains
        procedure :: confidence_interval => tukey_hsd_confidence_interval
    end type tukey_hsd_result

    type, public :: dunnett_result
        real(dp), allocatable :: statistic(:) !! treatment-vs-control t statistics
        real(dp), allocatable :: pvalue(:) !! single-step Dunnett p-values
        real(dp), allocatable :: rho(:, :) !! correlation matrix of comparisons
        real(dp), allocatable :: mean_samples(:) !! treatment sample means
        integer, allocatable :: n_samples(:) !! treatment sample sizes
        integer :: df = 0 !! pooled residual degrees of freedom
        integer :: n_control = 0 !! control sample size
        real(dp) :: std = 0.0_dp !! pooled residual standard deviation
        real(dp) :: mean_control = 0.0_dp !! control sample mean
        character(len=9) :: alternative = 'two-sided' !! test alternative
    contains
        procedure :: confidence_interval => dunnett_confidence_interval
    end type dunnett_result

    type, public :: poisson_means_test_result
        real(dp) :: statistic !! observed E-test pivot statistic
        real(dp) :: pvalue !! E-test p-value
    end type poisson_means_test_result

    public :: alexandergovern
    public :: dunnett
    public :: poisson_means_test
    public :: tukey_hsd

contains

    function alexandergovern(groups) result(result_value)
        type(sample_group), intent(in) :: groups(:) !! two or more independent samples, each size >= 2
        type(alexandergovern_result) :: result_value

        integer :: i
        integer :: k
        real(dp) :: a
        real(dp) :: b
        real(dp) :: c
        real(dp) :: eps
        real(dp), allocatable :: means(:)
        real(dp), allocatable :: se2(:)
        real(dp), allocatable :: se(:)
        real(dp), allocatable :: tstat(:)
        real(dp), allocatable :: weights(:)
        real(dp), allocatable :: z(:)
        real(dp) :: v
        real(dp) :: weighted_mean

        call invalidate_pair(result_value%statistic, result_value%pvalue)
        k = size(groups)
        if (.not. valid_groups(groups, 2)) return

        allocate(means(k), se2(k), se(k), tstat(k), weights(k), z(k))
        do i = 1, k
            means(i) = mean(groups(i)%values)
            se2(i) = variance(groups(i)%values) / real(size(groups(i)%values), dp)
            se(i) = sqrt(max(0.0_dp, se2(i)))
        end do
        eps = epsilon(1.0_dp)
        do i = 1, k
            if (se(i) <= abs(eps * means(i))) return
            if (se2(i) <= 0.0_dp .or. .not. ieee_is_finite(se2(i))) return
        end do

        weights = 1.0_dp / se2
        weights = weights / sum(weights)
        weighted_mean = sum(weights * means)
        tstat = (means - weighted_mean) / se

        do i = 1, k
            v = real(size(groups(i)%values) - 1, dp)
            a = v - 0.5_dp
            b = 48.0_dp * a * a
            c = sqrt(max(0.0_dp, a * log(1.0_dp + tstat(i) * tstat(i) / v)))
            z(i) = c + (c**3 + 3.0_dp * c) / b - &
                (4.0_dp * c**7 + 33.0_dp * c**5 + 240.0_dp * c**3 + 855.0_dp * c) / &
                (10.0_dp * b * b + 8.0_dp * b * c**4 + 1000.0_dp * b)
        end do
        result_value%statistic = sum(z * z)
        result_value%pvalue = chi2_sf(result_value%statistic, real(k - 1, dp))
    end function alexandergovern

    function tukey_hsd(groups, equal_var) result(result_value)
        type(sample_group), intent(in) :: groups(:) !! two or more independent samples of size > 1
        logical, intent(in), optional :: equal_var !! true: Tukey/Tukey-Kramer; false: Games-Howell
        type(tukey_hsd_result) :: result_value

        logical :: pooled
        logical :: equal_sizes
        integer :: i
        integer :: j
        integer :: k
        integer :: nobs
        integer, allocatable :: ns(:)
        real(dp), allocatable :: means(:)
        real(dp), allocatable :: vars(:)
        real(dp) :: denom
        real(dp) :: mse
        real(dp) :: s_i
        real(dp) :: s_j

        pooled = .true.
        if (present(equal_var)) pooled = equal_var
        k = size(groups)
        if (.not. valid_groups(groups, 2)) then
            call allocate_invalid_tukey(result_value, max(0, k))
            return
        end if

        allocate(ns(k), means(k), vars(k))
        nobs = 0
        do i = 1, k
            ns(i) = size(groups(i)%values)
            means(i) = mean(groups(i)%values)
            vars(i) = variance(groups(i)%values)
            nobs = nobs + ns(i)
        end do
        if (any(.not. ieee_is_finite(means)) .or. any(.not. ieee_is_finite(vars))) then
            call allocate_invalid_tukey(result_value, k)
            return
        end if

        result_value%ntreatments = k
        allocate(result_value%statistic(k, k), result_value%pvalue(k, k), &
            result_value%df(k, k), result_value%stand_err(k, k))
        equal_sizes = all(ns == ns(1))

        if (pooled) then
            denom = real(nobs - k, dp)
            if (denom <= 0.0_dp) then
                call fill_invalid_tukey(result_value)
                return
            end if
            mse = sum(vars * real(ns - 1, dp)) / denom
            do j = 1, k
                do i = 1, k
                    if (equal_sizes) then
                        result_value%stand_err(i, j) = sqrt(mse / real(ns(1), dp))
                    else
                        result_value%stand_err(i, j) = sqrt( &
                            (1.0_dp / real(ns(i), dp) + 1.0_dp / real(ns(j), dp)) * mse / 2.0_dp)
                    end if
                    result_value%df(i, j) = real(nobs - k, dp)
                end do
            end do
        else
            do j = 1, k
                do i = 1, k
                    s_i = vars(i) / real(ns(i), dp)
                    s_j = vars(j) / real(ns(j), dp)
                    result_value%stand_err(i, j) = sqrt(s_i + s_j) / sqrt(2.0_dp)
                    denom = s_i * s_i / real(ns(i) - 1, dp) + &
                        s_j * s_j / real(ns(j) - 1, dp)
                    if (denom > 0.0_dp) then
                        result_value%df(i, j) = (s_i + s_j) ** 2 / denom
                    else
                        result_value%df(i, j) = quiet_nan(0.0_dp)
                    end if
                end do
            end do
        end if

        do j = 1, k
            do i = 1, k
                result_value%statistic(i, j) = means(i) - means(j)
                if (result_value%stand_err(i, j) > 0.0_dp .and. &
                        result_value%df(i, j) > 0.0_dp) then
                    result_value%pvalue(i, j) = studentized_range_sf( &
                        abs(result_value%statistic(i, j)) / result_value%stand_err(i, j), &
                        real(k, dp), result_value%df(i, j))
                else if (i == j .and. result_value%stand_err(i, j) == 0.0_dp) then
                    result_value%pvalue(i, j) = 1.0_dp
                else
                    result_value%pvalue(i, j) = quiet_nan(0.0_dp)
                end if
            end do
        end do
    end function tukey_hsd

    function tukey_hsd_confidence_interval(self, confidence_level) result(ci)
        class(tukey_hsd_result), intent(in) :: self !! result returned by tukey_hsd
        real(dp), intent(in), optional :: confidence_level !! simultaneous confidence level, default 0.95
        type(matrix_confidence_interval) :: ci

        integer :: i
        integer :: j
        integer :: k
        real(dp) :: level
        real(dp) :: critical

        level = 0.95_dp
        if (present(confidence_level)) level = confidence_level
        k = self%ntreatments
        allocate(ci%low(k, k), ci%high(k, k))
        if (k < 2 .or. level <= 0.0_dp .or. level >= 1.0_dp .or. &
                .not. allocated(self%statistic)) then
            ci%low = quiet_nan(0.0_dp)
            ci%high = quiet_nan(0.0_dp)
            return
        end if

        do j = 1, k
            do i = 1, k
                critical = studentized_range_ppf(level, real(k, dp), self%df(i, j))
                ci%low(i, j) = self%statistic(i, j) - critical * self%stand_err(i, j)
                ci%high(i, j) = self%statistic(i, j) + critical * self%stand_err(i, j)
            end do
        end do
    end function tukey_hsd_confidence_interval

    function dunnett(groups, control, alternative, maxpts) result(result_value)
        type(sample_group), intent(in) :: groups(:) !! one or more treatment samples
        real(dp), intent(in) :: control(:) !! control sample
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        integer, intent(in), optional :: maxpts !! multivariate-t integration point budget
        type(dunnett_result) :: result_value

        integer :: i
        integer :: j
        integer :: k
        integer :: n_total
        integer :: points
        real(dp) :: pooled_ss
        real(dp) :: scale
        character(len=9) :: alt

        alt = choose_alternative(alternative)
        k = size(groups)
        call allocate_dunnett(result_value, max(0, k))
        result_value%alternative = alt
        if (k < 1 .or. size(control) < 1 .or. .not. valid_alternative(alt)) then
            call invalidate_dunnett(result_value)
            return
        end if
        if (any(.not. ieee_is_finite(control))) then
            call invalidate_dunnett(result_value)
            return
        end if
        do i = 1, k
            if (.not. allocated(groups(i)%values) .or. size(groups(i)%values) < 1) then
                call invalidate_dunnett(result_value)
                return
            end if
            if (any(.not. ieee_is_finite(groups(i)%values))) then
                call invalidate_dunnett(result_value)
                return
            end if
        end do

        result_value%n_control = size(control)
        result_value%mean_control = mean(control)
        n_total = size(control)
        do i = 1, k
            result_value%n_samples(i) = size(groups(i)%values)
            result_value%mean_samples(i) = mean(groups(i)%values)
            n_total = n_total + result_value%n_samples(i)
        end do
        result_value%df = n_total - k - 1
        if (result_value%df <= 0) then
            call invalidate_dunnett(result_value)
            return
        end if

        pooled_ss = sum((control - result_value%mean_control) ** 2)
        do i = 1, k
            pooled_ss = pooled_ss + sum((groups(i)%values - result_value%mean_samples(i)) ** 2)
        end do
        result_value%std = sqrt(pooled_ss / real(result_value%df, dp))
        if (result_value%std <= 0.0_dp .or. .not. ieee_is_finite(result_value%std)) then
            call invalidate_dunnett(result_value)
            return
        end if

        do j = 1, k
            do i = 1, k
                if (i == j) then
                    result_value%rho(i, j) = 1.0_dp
                else
                    result_value%rho(i, j) = 1.0_dp / sqrt( &
                        (real(result_value%n_control, dp) / real(result_value%n_samples(i), dp) + 1.0_dp) * &
                        (real(result_value%n_control, dp) / real(result_value%n_samples(j), dp) + 1.0_dp))
                end if
            end do
        end do
        do i = 1, k
            scale = sqrt(1.0_dp / real(result_value%n_samples(i), dp) + &
                1.0_dp / real(result_value%n_control, dp))
            result_value%statistic(i) = (result_value%mean_samples(i) - result_value%mean_control) / &
                (result_value%std * scale)
        end do

        points = max(8000, 4000 * k)
        if (present(maxpts)) then
            if (maxpts > 0) points = maxpts
        end if
        do i = 1, k
            result_value%pvalue(i) = dunnett_adjusted_pvalue(result_value%rho, &
                real(result_value%df, dp), result_value%statistic(i), alt, points)
        end do
    end function dunnett

    function dunnett_confidence_interval(self, confidence_level, maxpts) result(ci)
        class(dunnett_result), intent(in) :: self !! result returned by dunnett
        real(dp), intent(in), optional :: confidence_level !! simultaneous confidence level, default 0.95
        integer, intent(in), optional :: maxpts !! multivariate-t integration point budget
        type(vector_confidence_interval) :: ci

        integer :: i
        integer :: k
        integer :: points
        real(dp) :: allowance
        real(dp) :: critical
        real(dp) :: diff
        real(dp) :: level
        real(dp) :: scale

        level = 0.95_dp
        if (present(confidence_level)) level = confidence_level
        k = size(self%statistic)
        allocate(ci%low(k), ci%high(k))
        if (k < 1 .or. level <= 0.0_dp .or. level >= 1.0_dp .or. self%df <= 0 .or. &
                self%std <= 0.0_dp) then
            ci%low = quiet_nan(0.0_dp)
            ci%high = quiet_nan(0.0_dp)
            return
        end if
        points = max(10000, 5000 * k)
        if (present(maxpts)) then
            if (maxpts > 0) points = maxpts
        end if

        critical = dunnett_critical_value(self%rho, real(self%df, dp), &
            self%alternative, 1.0_dp - level, points)
        do i = 1, k
            scale = sqrt(1.0_dp / real(self%n_samples(i), dp) + &
                1.0_dp / real(self%n_control, dp))
            allowance = abs(critical) * self%std * scale
            diff = self%mean_samples(i) - self%mean_control
            ci%low(i) = diff - allowance
            ci%high(i) = diff + allowance
            if (trim(self%alternative) == 'greater') ci%high(i) = positive_infinity(0.0_dp)
            if (trim(self%alternative) == 'less') ci%low(i) = negative_infinity(0.0_dp)
        end do
    end function dunnett_confidence_interval

    function poisson_means_test(k1, n1, k2, n2, diff, alternative) result(result_value)
        integer, intent(in) :: k1 !! observed count for sample one, >= 0
        real(dp), intent(in) :: n1 !! exposure/sample size for sample one, > 0
        integer, intent(in) :: k2 !! observed count for sample two, >= 0
        real(dp), intent(in) :: n2 !! exposure/sample size for sample two, > 0
        real(dp), intent(in), optional :: diff !! null mean difference, >= 0, default 0
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        type(poisson_means_test_result) :: result_value

        character(len=9) :: alt
        integer :: x1
        integer :: x1_lb
        integer :: x1_ub
        integer :: x2
        integer :: x2_lb
        integer :: x2_ub
        real(dp) :: d
        real(dp) :: lambda_hat2
        real(dp) :: lambda_x1
        real(dp) :: lambda_x2
        real(dp) :: mu1
        real(dp) :: mu2
        real(dp) :: prob1
        real(dp) :: prob2
        real(dp) :: t_obs
        real(dp) :: t_value
        real(dp) :: var
        real(dp) :: var_x

        call invalidate_pair(result_value%statistic, result_value%pvalue)
        d = 0.0_dp
        if (present(diff)) d = diff
        alt = choose_alternative(alternative)
        if (k1 < 0 .or. k2 < 0 .or. n1 <= 0.0_dp .or. n2 <= 0.0_dp .or. &
                d < 0.0_dp .or. .not. valid_alternative(alt)) return
        if (.not. ieee_is_finite(n1) .or. .not. ieee_is_finite(n2) .or. &
                .not. ieee_is_finite(d)) return

        lambda_hat2 = (real(k1 + k2, dp) - d * n1) / (n1 + n2)
        if (lambda_hat2 <= 0.0_dp) then
            result_value%statistic = 0.0_dp
            result_value%pvalue = 1.0_dp
            return
        end if

        var = real(k1, dp) / (n1 * n1) + real(k2, dp) / (n2 * n2)
        if (var == 0.0_dp) then
            if (d == 0.0_dp) then
                t_obs = quiet_nan(0.0_dp)
            else
                t_obs = -positive_infinity(0.0_dp)
            end if
        else
            t_obs = (real(k1, dp) / n1 - real(k2, dp) / n2 - d) / sqrt(var)
        end if
        result_value%statistic = t_obs

        mu1 = n1 * (lambda_hat2 + d)
        mu2 = n2 * lambda_hat2
        x1_lb = max(0, int(poisson_ppf(1.0e-10_dp, mu1)))
        x2_lb = max(0, int(poisson_ppf(1.0e-10_dp, mu2)))
        x1_ub = max(x1_lb, int(poisson_ppf(1.0_dp - 1.0e-16_dp, mu1)))
        x2_ub = max(x2_lb, int(poisson_ppf(1.0_dp - 1.0e-16_dp, mu2)))

        result_value%pvalue = 0.0_dp
        do x2 = x2_lb, x2_ub
            prob2 = poisson_pmf(x2, mu2)
            lambda_x2 = real(x2, dp) / n2
            do x1 = x1_lb, x1_ub
                lambda_x1 = real(x1, dp) / n1
                var_x = lambda_x1 / n1 + lambda_x2 / n2
                if (var_x == 0.0_dp) then
                    if (lambda_x1 - lambda_x2 - d == 0.0_dp) then
                        cycle
                    else
                        t_value = sign(positive_infinity(0.0_dp), lambda_x1 - lambda_x2 - d)
                    end if
                else
                    t_value = (lambda_x1 - lambda_x2 - d) / sqrt(var_x)
                end if
                if (poisson_indicator(t_value, t_obs, alt)) then
                    prob1 = poisson_pmf(x1, mu1)
                    result_value%pvalue = result_value%pvalue + prob1 * prob2
                end if
            end do
        end do
        result_value%pvalue = min(1.0_dp, max(0.0_dp, result_value%pvalue))
    end function poisson_means_test

    pure logical function poisson_indicator(value, observed, alternative) result(include)
        real(dp), intent(in) :: value !! candidate pivot statistic
        real(dp), intent(in) :: observed !! observed pivot statistic
        character(len=*), intent(in) :: alternative !! validated alternative

        if (ieee_is_nan(value) .or. ieee_is_nan(observed)) then
            include = .false.
        else if (trim(alternative) == 'less') then
            include = value <= observed
        else if (trim(alternative) == 'greater') then
            include = value >= observed
        else
            include = abs(value) >= abs(observed)
        end if
    end function poisson_indicator

    function dunnett_adjusted_pvalue(rho, df, statistic, alternative, maxpts) result(pvalue)
        real(dp), intent(in) :: rho(:, :) !! comparison correlation matrix
        real(dp), intent(in) :: df !! residual degrees of freedom
        real(dp), intent(in) :: statistic !! scalar comparison statistic
        character(len=*), intent(in) :: alternative !! validated alternative
        integer, intent(in) :: maxpts !! multivariate-t integration point budget
        real(dp) :: pvalue

        integer :: k
        real(dp), allocatable :: high(:)
        real(dp), allocatable :: low(:)
        real(dp) :: c

        k = size(rho, 1)
        allocate(high(k), low(k))
        if (trim(alternative) == 'two-sided') then
            c = abs(statistic)
            high = c
            low = -c
        else if (trim(alternative) == 'greater') then
            high = statistic
            low = negative_infinity(0.0_dp)
        else
            high = positive_infinity(0.0_dp)
            low = statistic
        end if
        pvalue = 1.0_dp - multivariate_t_cdf(high, shape=rho, df=df, &
            allow_singular=.false., maxpts=maxpts, lower_limit=low)
        pvalue = min(1.0_dp, max(0.0_dp, pvalue))
    end function dunnett_adjusted_pvalue

    function dunnett_critical_value(rho, df, alternative, alpha, maxpts) result(critical)
        real(dp), intent(in) :: rho(:, :) !! comparison correlation matrix
        real(dp), intent(in) :: df !! residual degrees of freedom
        character(len=*), intent(in) :: alternative !! validated alternative
        real(dp), intent(in) :: alpha !! target familywise tail probability
        integer, intent(in) :: maxpts !! multivariate-t integration point budget
        real(dp) :: critical

        integer :: iter
        real(dp) :: hi
        real(dp) :: lo
        real(dp) :: mid
        real(dp) :: p
        real(dp) :: sign_value

        sign_value = 1.0_dp
        if (trim(alternative) == 'less') sign_value = -1.0_dp
        lo = 0.0_dp
        hi = 1.0_dp
        do iter = 1, 64
            p = dunnett_adjusted_pvalue(rho, df, sign_value * hi, alternative, maxpts)
            if (p <= alpha) exit
            hi = 2.0_dp * hi
        end do
        do iter = 1, 24
            mid = 0.5_dp * (lo + hi)
            p = dunnett_adjusted_pvalue(rho, df, sign_value * mid, alternative, maxpts)
            if (p > alpha) then
                lo = mid
            else
                hi = mid
            end if
        end do
        critical = sign_value * 0.5_dp * (lo + hi)
    end function dunnett_critical_value

    subroutine allocate_invalid_tukey(result_value, k)
        type(tukey_hsd_result), intent(out) :: result_value !! result to allocate and invalidate
        integer, intent(in) :: k !! result matrix dimension

        result_value%ntreatments = k
        allocate(result_value%statistic(k, k), result_value%pvalue(k, k), &
            result_value%df(k, k), result_value%stand_err(k, k))
        call fill_invalid_tukey(result_value)
    end subroutine allocate_invalid_tukey

    subroutine fill_invalid_tukey(result_value)
        type(tukey_hsd_result), intent(inout) :: result_value !! allocated result to invalidate

        result_value%statistic = quiet_nan(0.0_dp)
        result_value%pvalue = quiet_nan(0.0_dp)
        result_value%df = quiet_nan(0.0_dp)
        result_value%stand_err = quiet_nan(0.0_dp)
    end subroutine fill_invalid_tukey

    subroutine allocate_dunnett(result_value, k)
        type(dunnett_result), intent(out) :: result_value !! result to allocate
        integer, intent(in) :: k !! number of treatment samples

        allocate(result_value%statistic(k), result_value%pvalue(k), result_value%rho(k, k), &
            result_value%mean_samples(k), result_value%n_samples(k))
        result_value%statistic = quiet_nan(0.0_dp)
        result_value%pvalue = quiet_nan(0.0_dp)
        result_value%rho = quiet_nan(0.0_dp)
        result_value%mean_samples = quiet_nan(0.0_dp)
        result_value%n_samples = 0
    end subroutine allocate_dunnett

    subroutine invalidate_dunnett(result_value)
        type(dunnett_result), intent(inout) :: result_value !! allocated result to invalidate

        result_value%statistic = quiet_nan(0.0_dp)
        result_value%pvalue = quiet_nan(0.0_dp)
        result_value%rho = quiet_nan(0.0_dp)
        result_value%mean_samples = quiet_nan(0.0_dp)
        result_value%n_samples = 0
        result_value%df = 0
        result_value%n_control = 0
        result_value%std = quiet_nan(0.0_dp)
        result_value%mean_control = quiet_nan(0.0_dp)
    end subroutine invalidate_dunnett

    pure logical function valid_groups(groups, minimum_groups) result(ok)
        type(sample_group), intent(in) :: groups(:) !! sample groups to validate
        integer, intent(in) :: minimum_groups !! required minimum number of groups
        integer :: i

        ok = size(groups) >= minimum_groups
        if (.not. ok) return
        do i = 1, size(groups)
            if (.not. allocated(groups(i)%values)) then
                ok = .false.
                return
            end if
            if (size(groups(i)%values) < 2) then
                ok = .false.
                return
            end if
            if (any(.not. ieee_is_finite(groups(i)%values))) then
                ok = .false.
                return
            end if
        end do
    end function valid_groups

    pure function choose_alternative(alternative) result(selected)
        character(len=*), intent(in), optional :: alternative !! optional user selection
        character(len=9) :: selected

        selected = 'two-sided'
        if (present(alternative)) selected = adjustl(alternative)
    end function choose_alternative

    pure logical function valid_alternative(alternative) result(ok)
        character(len=*), intent(in) :: alternative !! alternative name

        ok = trim(alternative) == 'two-sided' .or. trim(alternative) == 'less' .or. &
            trim(alternative) == 'greater'
    end function valid_alternative

    pure subroutine invalidate_pair(statistic, pvalue)
        real(dp), intent(out) :: statistic !! statistic set to NaN
        real(dp), intent(out) :: pvalue !! p-value set to NaN

        statistic = quiet_nan(0.0_dp)
        pvalue = quiet_nan(0.0_dp)
    end subroutine invalidate_pair

end module scifort_multiple_comparisons
