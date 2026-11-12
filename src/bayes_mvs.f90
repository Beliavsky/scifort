! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Bayesian posterior distributions and intervals for a normal sample with
! unknown mean and variance, matching scipy.stats.mvsdist and bayes_mvs.
!
! SciPy 1.17.0 uses Jeffreys' prior and represents the posterior mean,
! variance, and standard deviation by Student-t, inverse-gamma, and
! generalized-gamma distributions for n <= 1000, respectively.  For larger
! samples it uses normal approximations.  The formulas and n > 1000 switch are
! retained here; the Fortran API exposes small frozen-posterior objects instead
! of Python rv_frozen objects.

module scifort_bayes_mvs
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_gengamma, only : gengamma_cdf, gengamma_isf, gengamma_logcdf, &
        gengamma_logpdf, gengamma_logsf, gengamma_pdf, gengamma_ppf, gengamma_sf
    use scifort_invgamma, only : invgamma_cdf, invgamma_isf, invgamma_logcdf, &
        invgamma_logpdf, invgamma_logsf, invgamma_pdf, invgamma_ppf, invgamma_sf
    use scifort_kinds, only : dp
    use scifort_math, only : positive_infinity, quiet_nan
    use scifort_normal, only : normal_cdf, normal_isf, normal_logcdf, normal_logpdf, &
        normal_logsf, normal_pdf, normal_ppf, normal_sf
    use scifort_student_t, only : t_cdf, t_isf, t_logcdf, t_logpdf, t_logsf, &
        t_pdf, t_ppf, t_sf
    implicit none
    private

    integer, parameter, public :: mvs_status_success = 0
    integer, parameter, public :: mvs_status_invalid_data = 1
    integer, parameter, public :: mvs_status_invalid_alpha = 2

    integer, parameter :: posterior_mean = 1
    integer, parameter :: posterior_variance = 2
    integer, parameter :: posterior_std = 3

    type, public :: mvs_posterior
        private
        integer :: posterior_kind = 0
        integer :: n = 0
        real(dp) :: xbar = 0.0_dp
        real(dp) :: sample_variance = 0.0_dp
    contains
        procedure, public :: cdf => mvs_posterior_cdf
        procedure, public :: interval => mvs_posterior_interval
        procedure, public :: isf => mvs_posterior_isf
        procedure, public :: logcdf => mvs_posterior_logcdf
        procedure, public :: logpdf => mvs_posterior_logpdf
        procedure, public :: logsf => mvs_posterior_logsf
        procedure, public :: mean => mvs_posterior_mean
        procedure, public :: pdf => mvs_posterior_pdf
        procedure, public :: ppf => mvs_posterior_ppf
        procedure, public :: sf => mvs_posterior_sf
        procedure, public :: std => mvs_posterior_std
        procedure, public :: var => mvs_posterior_var
    end type mvs_posterior

    type, public :: mvsdist_result
        type(mvs_posterior) :: mean !! posterior distribution of the population mean
        type(mvs_posterior) :: variance !! posterior distribution of the population variance
        type(mvs_posterior) :: std_dev !! posterior distribution of the population standard deviation
    end type mvsdist_result

    type, public :: bayes_mvs_estimate
        real(dp) :: statistic !! posterior mean of the parameter
        real(dp) :: minmax(2) !! central probability interval about the posterior median
    end type bayes_mvs_estimate

    type, public :: bayes_mvs_result
        type(bayes_mvs_estimate) :: mean !! posterior summary for the population mean
        type(bayes_mvs_estimate) :: variance !! posterior summary for the population variance
        type(bayes_mvs_estimate) :: std_dev !! posterior summary for the population standard deviation
    end type bayes_mvs_result

    public :: bayes_mvs
    public :: mvsdist

contains

    function mvsdist(data, status) result(result_value)
        real(dp), intent(in) :: data(..) !! sample; any rank is flattened in array element order
        integer, intent(out), optional :: status !! zero on success, otherwise mvs_status_invalid_data
        type(mvsdist_result) :: result_value

        real(dp), allocatable :: flat(:)

        call flatten_data(data, flat)
        if (size(flat) < 2) then
            call invalidate_mvsdist(result_value)
            if (present(status)) status = mvs_status_invalid_data
            return
        end if

        call initialize_mvsdist(flat, result_value)
        if (present(status)) status = mvs_status_success
    end function mvsdist

    function bayes_mvs(data, alpha, status) result(result_value)
        real(dp), intent(in) :: data(..) !! sample; any rank is flattened in array element order
        real(dp), intent(in), optional :: alpha !! interval probability, default 0.9; must be in (0,1)
        integer, intent(out), optional :: status !! mvs_status_* result code
        type(bayes_mvs_result) :: result_value

        real(dp) :: probability
        integer :: local_status
        type(mvsdist_result) :: distributions

        probability = 0.9_dp
        if (present(alpha)) probability = alpha
        if (.not. ieee_is_finite(probability)) then
            call invalidate_bayes_result(result_value)
            if (present(status)) status = mvs_status_invalid_alpha
            return
        end if
        if (probability <= 0.0_dp .or. probability >= 1.0_dp) then
            call invalidate_bayes_result(result_value)
            if (present(status)) status = mvs_status_invalid_alpha
            return
        end if

        distributions = mvsdist(data, local_status)
        if (local_status /= mvs_status_success) then
            call invalidate_bayes_result(result_value)
            if (present(status)) status = local_status
            return
        end if

        result_value%mean%statistic = distributions%mean%mean()
        result_value%mean%minmax = distributions%mean%interval(probability)
        result_value%variance%statistic = distributions%variance%mean()
        result_value%variance%minmax = distributions%variance%interval(probability)
        result_value%std_dev%statistic = distributions%std_dev%mean()
        result_value%std_dev%minmax = distributions%std_dev%interval(probability)
        if (present(status)) status = mvs_status_success
    end function bayes_mvs

    subroutine initialize_mvsdist(data, result_value)
        real(dp), intent(in) :: data(:) !! flattened sample with at least two observations
        type(mvsdist_result), intent(out) :: result_value !! initialized posterior distributions

        integer :: n
        real(dp) :: c
        real(dp) :: xbar

        n = size(data)
        if (any(.not. ieee_is_finite(data))) then
            xbar = quiet_nan(0.0_dp)
            c = quiet_nan(0.0_dp)
        else
            xbar = sum(data) / real(n, dp)
            c = sum((data - xbar) ** 2) / real(n, dp)
        end if

        result_value%mean%posterior_kind = posterior_mean
        result_value%variance%posterior_kind = posterior_variance
        result_value%std_dev%posterior_kind = posterior_std
        result_value%mean%n = n
        result_value%variance%n = n
        result_value%std_dev%n = n
        result_value%mean%xbar = xbar
        result_value%variance%xbar = xbar
        result_value%std_dev%xbar = xbar
        result_value%mean%sample_variance = c
        result_value%variance%sample_variance = c
        result_value%std_dev%sample_variance = c
    end subroutine initialize_mvsdist

    subroutine invalidate_mvsdist(result_value)
        type(mvsdist_result), intent(out) :: result_value !! result set to invalid posterior objects

        result_value%mean%posterior_kind = 0
        result_value%variance%posterior_kind = 0
        result_value%std_dev%posterior_kind = 0
    end subroutine invalidate_mvsdist

    subroutine invalidate_bayes_result(result_value)
        type(bayes_mvs_result), intent(out) :: result_value !! result set to NaNs

        real(dp) :: nan_value

        nan_value = quiet_nan(0.0_dp)
        result_value%mean%statistic = nan_value
        result_value%mean%minmax = nan_value
        result_value%variance%statistic = nan_value
        result_value%variance%minmax = nan_value
        result_value%std_dev%statistic = nan_value
        result_value%std_dev%minmax = nan_value
    end subroutine invalidate_bayes_result

    pure function mvs_posterior_pdf(self, x) result(value)
        class(mvs_posterior), intent(in) :: self !! frozen posterior distribution
        real(dp), intent(in) :: x !! evaluation point
        real(dp) :: value

        real(dp) :: a
        real(dp) :: df
        real(dp) :: loc
        real(dp) :: scale

        call posterior_parameters(self, a, df, loc, scale)
        if (.not. valid_posterior(self)) then
            value = quiet_nan(x)
        else if (.not. (scale > 0.0_dp)) then
            value = quiet_nan(x)
        else if (self%n > 1000) then
            value = normal_pdf(x, loc, scale)
        else
            select case (self%posterior_kind)
            case (posterior_mean)
                value = t_pdf(x, df, loc, scale)
            case (posterior_variance)
                value = invgamma_pdf(x, a, 0.0_dp, scale)
            case (posterior_std)
                value = gengamma_pdf(x, a, -2.0_dp, 0.0_dp, scale)
            case default
                value = quiet_nan(x)
            end select
        end if
    end function mvs_posterior_pdf

    pure function mvs_posterior_logpdf(self, x) result(value)
        class(mvs_posterior), intent(in) :: self !! frozen posterior distribution
        real(dp), intent(in) :: x !! evaluation point
        real(dp) :: value

        real(dp) :: a
        real(dp) :: df
        real(dp) :: loc
        real(dp) :: scale

        call posterior_parameters(self, a, df, loc, scale)
        if (.not. valid_posterior(self)) then
            value = quiet_nan(x)
        else if (.not. (scale > 0.0_dp)) then
            value = quiet_nan(x)
        else if (self%n > 1000) then
            value = normal_logpdf(x, loc, scale)
        else
            select case (self%posterior_kind)
            case (posterior_mean)
                value = t_logpdf(x, df, loc, scale)
            case (posterior_variance)
                value = invgamma_logpdf(x, a, 0.0_dp, scale)
            case (posterior_std)
                value = gengamma_logpdf(x, a, -2.0_dp, 0.0_dp, scale)
            case default
                value = quiet_nan(x)
            end select
        end if
    end function mvs_posterior_logpdf

    pure function mvs_posterior_cdf(self, x) result(value)
        class(mvs_posterior), intent(in) :: self !! frozen posterior distribution
        real(dp), intent(in) :: x !! evaluation point
        real(dp) :: value

        real(dp) :: a
        real(dp) :: df
        real(dp) :: loc
        real(dp) :: scale

        call posterior_parameters(self, a, df, loc, scale)
        if (.not. valid_posterior(self)) then
            value = quiet_nan(x)
        else if (.not. (scale > 0.0_dp)) then
            value = quiet_nan(x)
        else if (self%n > 1000) then
            value = normal_cdf(x, loc, scale)
        else
            select case (self%posterior_kind)
            case (posterior_mean)
                value = t_cdf(x, df, loc, scale)
            case (posterior_variance)
                value = invgamma_cdf(x, a, 0.0_dp, scale)
            case (posterior_std)
                value = gengamma_cdf(x, a, -2.0_dp, 0.0_dp, scale)
            case default
                value = quiet_nan(x)
            end select
        end if
    end function mvs_posterior_cdf

    pure function mvs_posterior_logcdf(self, x) result(value)
        class(mvs_posterior), intent(in) :: self !! frozen posterior distribution
        real(dp), intent(in) :: x !! evaluation point
        real(dp) :: value

        real(dp) :: a
        real(dp) :: df
        real(dp) :: loc
        real(dp) :: scale

        call posterior_parameters(self, a, df, loc, scale)
        if (.not. valid_posterior(self)) then
            value = quiet_nan(x)
        else if (.not. (scale > 0.0_dp)) then
            value = quiet_nan(x)
        else if (self%n > 1000) then
            value = normal_logcdf(x, loc, scale)
        else
            select case (self%posterior_kind)
            case (posterior_mean)
                value = t_logcdf(x, df, loc, scale)
            case (posterior_variance)
                value = invgamma_logcdf(x, a, 0.0_dp, scale)
            case (posterior_std)
                value = gengamma_logcdf(x, a, -2.0_dp, 0.0_dp, scale)
            case default
                value = quiet_nan(x)
            end select
        end if
    end function mvs_posterior_logcdf

    pure function mvs_posterior_sf(self, x) result(value)
        class(mvs_posterior), intent(in) :: self !! frozen posterior distribution
        real(dp), intent(in) :: x !! evaluation point
        real(dp) :: value

        real(dp) :: a
        real(dp) :: df
        real(dp) :: loc
        real(dp) :: scale

        call posterior_parameters(self, a, df, loc, scale)
        if (.not. valid_posterior(self)) then
            value = quiet_nan(x)
        else if (.not. (scale > 0.0_dp)) then
            value = quiet_nan(x)
        else if (self%n > 1000) then
            value = normal_sf(x, loc, scale)
        else
            select case (self%posterior_kind)
            case (posterior_mean)
                value = t_sf(x, df, loc, scale)
            case (posterior_variance)
                value = invgamma_sf(x, a, 0.0_dp, scale)
            case (posterior_std)
                value = gengamma_sf(x, a, -2.0_dp, 0.0_dp, scale)
            case default
                value = quiet_nan(x)
            end select
        end if
    end function mvs_posterior_sf

    pure function mvs_posterior_logsf(self, x) result(value)
        class(mvs_posterior), intent(in) :: self !! frozen posterior distribution
        real(dp), intent(in) :: x !! evaluation point
        real(dp) :: value

        real(dp) :: a
        real(dp) :: df
        real(dp) :: loc
        real(dp) :: scale

        call posterior_parameters(self, a, df, loc, scale)
        if (.not. valid_posterior(self)) then
            value = quiet_nan(x)
        else if (.not. (scale > 0.0_dp)) then
            value = quiet_nan(x)
        else if (self%n > 1000) then
            value = normal_logsf(x, loc, scale)
        else
            select case (self%posterior_kind)
            case (posterior_mean)
                value = t_logsf(x, df, loc, scale)
            case (posterior_variance)
                value = invgamma_logsf(x, a, 0.0_dp, scale)
            case (posterior_std)
                value = gengamma_logsf(x, a, -2.0_dp, 0.0_dp, scale)
            case default
                value = quiet_nan(x)
            end select
        end if
    end function mvs_posterior_logsf

    pure function mvs_posterior_ppf(self, probability) result(value)
        class(mvs_posterior), intent(in) :: self !! frozen posterior distribution
        real(dp), intent(in) :: probability !! lower-tail probability in [0,1]
        real(dp) :: value

        real(dp) :: a
        real(dp) :: df
        real(dp) :: loc
        real(dp) :: scale

        if (.not. ieee_is_finite(probability)) then
            value = quiet_nan(probability)
            return
        end if
        call posterior_parameters(self, a, df, loc, scale)
        if (.not. valid_posterior(self)) then
            value = quiet_nan(probability)
        else if (.not. (scale > 0.0_dp)) then
            value = quiet_nan(probability)
        else if (self%n > 1000) then
            value = normal_ppf(probability, loc, scale)
        else
            select case (self%posterior_kind)
            case (posterior_mean)
                value = t_ppf(probability, df, loc, scale)
            case (posterior_variance)
                value = invgamma_ppf(probability, a, 0.0_dp, scale)
            case (posterior_std)
                value = gengamma_ppf(probability, a, -2.0_dp, 0.0_dp, scale)
            case default
                value = quiet_nan(probability)
            end select
        end if
    end function mvs_posterior_ppf

    pure function mvs_posterior_isf(self, probability) result(value)
        class(mvs_posterior), intent(in) :: self !! frozen posterior distribution
        real(dp), intent(in) :: probability !! upper-tail probability in [0,1]
        real(dp) :: value

        real(dp) :: a
        real(dp) :: df
        real(dp) :: loc
        real(dp) :: scale

        if (.not. ieee_is_finite(probability)) then
            value = quiet_nan(probability)
            return
        end if
        call posterior_parameters(self, a, df, loc, scale)
        if (.not. valid_posterior(self)) then
            value = quiet_nan(probability)
        else if (.not. (scale > 0.0_dp)) then
            value = quiet_nan(probability)
        else if (self%n > 1000) then
            value = normal_isf(probability, loc, scale)
        else
            select case (self%posterior_kind)
            case (posterior_mean)
                value = t_isf(probability, df, loc, scale)
            case (posterior_variance)
                value = invgamma_isf(probability, a, 0.0_dp, scale)
            case (posterior_std)
                value = gengamma_isf(probability, a, -2.0_dp, 0.0_dp, scale)
            case default
                value = quiet_nan(probability)
            end select
        end if
    end function mvs_posterior_isf

    pure function mvs_posterior_interval(self, alpha) result(limits)
        class(mvs_posterior), intent(in) :: self !! frozen posterior distribution
        real(dp), intent(in) :: alpha !! central probability in [0,1]
        real(dp) :: limits(2)

        real(dp) :: lower_probability

        if (.not. ieee_is_finite(alpha)) then
            limits = quiet_nan(alpha)
            return
        end if
        if (alpha < 0.0_dp .or. alpha > 1.0_dp) then
            limits = quiet_nan(alpha)
            return
        end if
        lower_probability = 0.5_dp * (1.0_dp - alpha)
        limits(1) = self%ppf(lower_probability)
        limits(2) = self%ppf(1.0_dp - lower_probability)
    end function mvs_posterior_interval

    pure function mvs_posterior_mean(self) result(value)
        class(mvs_posterior), intent(in) :: self !! frozen posterior distribution
        real(dp) :: value

        real(dp) :: a
        real(dp) :: df
        real(dp) :: loc
        real(dp) :: scale

        call posterior_parameters(self, a, df, loc, scale)
        if (.not. valid_posterior(self)) then
            value = quiet_nan(0.0_dp)
            return
        end if
        if (.not. (scale > 0.0_dp)) then
            value = quiet_nan(0.0_dp)
            return
        end if
        if (self%n > 1000) then
            value = loc
            return
        end if

        select case (self%posterior_kind)
        case (posterior_mean)
            if (df > 1.0_dp) then
                value = loc
            else
                value = positive_infinity(0.0_dp)
            end if
        case (posterior_variance)
            if (a > 1.0_dp) then
                value = scale / (a - 1.0_dp)
            else
                value = positive_infinity(0.0_dp)
            end if
        case (posterior_std)
            if (a > 0.5_dp) then
                value = scale * exp(log_gamma(a - 0.5_dp) - log_gamma(a))
            else
                value = positive_infinity(0.0_dp)
            end if
        case default
            value = quiet_nan(0.0_dp)
        end select
    end function mvs_posterior_mean

    pure function mvs_posterior_var(self) result(value)
        class(mvs_posterior), intent(in) :: self !! frozen posterior distribution
        real(dp) :: value

        real(dp) :: a
        real(dp) :: df
        real(dp) :: first_moment
        real(dp) :: loc
        real(dp) :: scale
        real(dp) :: second_moment

        call posterior_parameters(self, a, df, loc, scale)
        if (.not. valid_posterior(self)) then
            value = quiet_nan(0.0_dp)
            return
        end if
        if (.not. (scale > 0.0_dp)) then
            value = quiet_nan(0.0_dp)
            return
        end if
        if (self%n > 1000) then
            value = scale * scale
            return
        end if

        select case (self%posterior_kind)
        case (posterior_mean)
            if (df > 2.0_dp) then
                value = scale * scale * df / (df - 2.0_dp)
            else if (df > 1.0_dp) then
                value = positive_infinity(0.0_dp)
            else
                value = quiet_nan(0.0_dp)
            end if
        case (posterior_variance)
            if (a > 2.0_dp) then
                value = scale * scale / ((a - 1.0_dp) ** 2 * (a - 2.0_dp))
            else
                value = positive_infinity(0.0_dp)
            end if
        case (posterior_std)
            if (a > 1.0_dp) then
                first_moment = scale * exp(log_gamma(a - 0.5_dp) - log_gamma(a))
                second_moment = scale * scale / (a - 1.0_dp)
                value = max(0.0_dp, second_moment - first_moment * first_moment)
            else
                value = positive_infinity(0.0_dp)
            end if
        case default
            value = quiet_nan(0.0_dp)
        end select
    end function mvs_posterior_var

    pure function mvs_posterior_std(self) result(value)
        class(mvs_posterior), intent(in) :: self !! frozen posterior distribution
        real(dp) :: value

        value = sqrt(self%var())
    end function mvs_posterior_std

    pure subroutine posterior_parameters(self, a, df, loc, scale)
        class(mvs_posterior), intent(in) :: self !! frozen posterior distribution
        real(dp), intent(out) :: a !! inverse/generalized-gamma shape
        real(dp), intent(out) :: df !! Student-t degrees of freedom
        real(dp), intent(out) :: loc !! location of selected posterior
        real(dp), intent(out) :: scale !! scale of selected posterior

        real(dp) :: c
        real(dp) :: n_real

        a = quiet_nan(0.0_dp)
        df = quiet_nan(0.0_dp)
        loc = quiet_nan(0.0_dp)
        scale = quiet_nan(0.0_dp)
        if (.not. valid_posterior(self)) return

        c = self%sample_variance
        n_real = real(self%n, dp)
        if (self%n > 1000) then
            select case (self%posterior_kind)
            case (posterior_mean)
                loc = self%xbar
                scale = sqrt(c / n_real)
            case (posterior_variance)
                loc = c
                scale = sqrt(2.0_dp / n_real) * c
            case (posterior_std)
                loc = sqrt(c)
                scale = sqrt(c / (2.0_dp * n_real))
            end select
        else
            df = real(self%n - 1, dp)
            a = 0.5_dp * df
            select case (self%posterior_kind)
            case (posterior_mean)
                loc = self%xbar
                scale = sqrt(c / df)
            case (posterior_variance)
                loc = 0.0_dp
                scale = 0.5_dp * n_real * c
            case (posterior_std)
                loc = 0.0_dp
                scale = sqrt(0.5_dp * n_real * c)
            end select
        end if
    end subroutine posterior_parameters

    pure logical function valid_posterior(self)
        class(mvs_posterior), intent(in) :: self !! frozen posterior distribution

        valid_posterior = .false.
        if (self%n < 2) return
        if (self%posterior_kind < posterior_mean .or. self%posterior_kind > posterior_std) return
        if (.not. ieee_is_finite(self%xbar)) return
        if (.not. ieee_is_finite(self%sample_variance)) return
        if (self%sample_variance < 0.0_dp) return
        valid_posterior = .true.
    end function valid_posterior

    subroutine flatten_data(data, flat)
        real(dp), intent(in) :: data(..) !! assumed-rank source array
        real(dp), allocatable, intent(out) :: flat(:) !! rank-one copy

        select rank (data)
        rank (0)
            flat = [data]
        rank (1)
            flat = data
        rank (2)
            flat = reshape(data, [size(data)])
        rank (3)
            flat = reshape(data, [size(data)])
        rank (4)
            flat = reshape(data, [size(data)])
        rank (5)
            flat = reshape(data, [size(data)])
        rank (6)
            flat = reshape(data, [size(data)])
        rank (7)
            flat = reshape(data, [size(data)])
        rank (8)
            flat = reshape(data, [size(data)])
        rank (9)
            flat = reshape(data, [size(data)])
        rank (10)
            flat = reshape(data, [size(data)])
        rank (11)
            flat = reshape(data, [size(data)])
        rank (12)
            flat = reshape(data, [size(data)])
        rank (13)
            flat = reshape(data, [size(data)])
        rank (14)
            flat = reshape(data, [size(data)])
        rank (15)
            flat = reshape(data, [size(data)])
        rank default
            allocate(flat(0))
        end select
    end subroutine flatten_data

end module scifort_bayes_mvs
