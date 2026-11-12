! SPDX-License-Identifier: MIT AND BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright (c) 2001-2002 Enthought, Inc. 2003-2026, SciPy Developers.
! Adapted from SciPy 1.17.0 scipy/stats/_distribution_infrastructure.py.
! Retained BSD terms: THIRD_PARTY_LICENSES.md; details: CODE_PROVENANCE.md.
!
! First-class scalar continuous-distribution adapters and truncation helpers.
! SciPy's public make_distribution/truncate semantics and truncation formulas are
! adapted to an idiomatic Fortran callback interface and explicit RNG state.

module scifort_distribution_infrastructure
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, quiet_nan
    use scifort_random, only : rng_state, rng_uniform
    implicit none
    private

    integer, parameter, public :: distribution_status_success = 0
    integer, parameter, public :: distribution_status_invalid_support = 1
    integer, parameter, public :: distribution_status_invalid_bounds = 2
    integer, parameter, public :: distribution_status_zero_mass = 3

    abstract interface
        pure function distribution_scalar_callback(x, parameters) result(y)
            import dp
            real(dp), intent(in) :: x !! scalar argument or probability
            real(dp), intent(in) :: parameters(:) !! frozen distribution parameters
            real(dp) :: y
        end function distribution_scalar_callback
    end interface

    type, public :: continuous_distribution
        private
        real(dp), allocatable :: parameters(:)
        real(dp) :: lower = 0.0_dp
        real(dp) :: upper = 0.0_dp
        logical :: inclusive(2) = [.false., .false.]
        logical :: valid = .false.
        procedure(distribution_scalar_callback), pointer, nopass :: pdf_callback => null()
        procedure(distribution_scalar_callback), pointer, nopass :: cdf_callback => null()
        procedure(distribution_scalar_callback), pointer, nopass :: ppf_callback => null()
        procedure(distribution_scalar_callback), pointer, nopass :: logpdf_callback => null()
        procedure(distribution_scalar_callback), pointer, nopass :: logcdf_callback => null()
        procedure(distribution_scalar_callback), pointer, nopass :: sf_callback => null()
        procedure(distribution_scalar_callback), pointer, nopass :: logsf_callback => null()
        procedure(distribution_scalar_callback), pointer, nopass :: isf_callback => null()
    contains
        procedure, public :: ccdf => continuous_sf
        procedure, public :: cdf => continuous_cdf
        procedure, public :: iccdf => continuous_isf
        procedure, public :: icdf => continuous_ppf
        procedure, public :: ilogccdf => continuous_ilogsf
        procedure, public :: ilogcdf => continuous_ilogcdf
        procedure, public :: isf => continuous_isf
        procedure, public :: logccdf => continuous_logsf
        procedure, public :: logcdf => continuous_logcdf
        procedure, public :: logpdf => continuous_logpdf
        procedure, public :: logsf => continuous_logsf
        procedure, public :: median => continuous_median
        procedure, public :: pdf => continuous_pdf
        procedure, public :: ppf => continuous_ppf
        procedure, public :: sample => continuous_sample
        procedure, public :: sf => continuous_sf
        procedure, public :: support => continuous_support
    end type continuous_distribution

    type, public :: truncated_distribution
        private
        type(continuous_distribution) :: base
        real(dp) :: lb = 0.0_dp
        real(dp) :: ub = 0.0_dp
        real(dp) :: lower = 0.0_dp
        real(dp) :: upper = 0.0_dp
        real(dp) :: logmass = 0.0_dp
        logical :: valid = .false.
    contains
        procedure, public :: ccdf => truncated_sf
        procedure, public :: cdf => truncated_cdf
        procedure, public :: iccdf => truncated_isf
        procedure, public :: icdf => truncated_ppf
        procedure, public :: ilogccdf => truncated_ilogsf
        procedure, public :: ilogcdf => truncated_ilogcdf
        procedure, public :: isf => truncated_isf
        procedure, public :: logccdf => truncated_logsf
        procedure, public :: logcdf => truncated_logcdf
        procedure, public :: logpdf => truncated_logpdf
        procedure, public :: logsf => truncated_logsf
        procedure, public :: median => truncated_median
        procedure, public :: pdf => truncated_pdf
        procedure, public :: ppf => truncated_ppf
        procedure, public :: sample => truncated_sample
        procedure, public :: sf => truncated_sf
        procedure, public :: support => truncated_support
    end type truncated_distribution

    public :: distribution_scalar_callback
    public :: make_distribution
    public :: truncate

contains

    function make_distribution(pdf, cdf, ppf, support, parameters, inclusive, &
            logpdf, logcdf, sf, logsf, isf, status) result(dist)
        procedure(distribution_scalar_callback) :: pdf !! base PDF callback
        procedure(distribution_scalar_callback) :: cdf !! base CDF callback
        procedure(distribution_scalar_callback) :: ppf !! base inverse-CDF callback
        real(dp), intent(in) :: support(2) !! lower and upper support endpoints
        real(dp), intent(in), optional :: parameters(:) !! frozen callback parameters
        logical, intent(in), optional :: inclusive(2) !! support endpoint inclusion flags
        procedure(distribution_scalar_callback), optional :: logpdf !! stable log-PDF callback
        procedure(distribution_scalar_callback), optional :: logcdf !! stable log-CDF callback
        procedure(distribution_scalar_callback), optional :: sf !! stable survival callback
        procedure(distribution_scalar_callback), optional :: logsf !! stable log-survival callback
        procedure(distribution_scalar_callback), optional :: isf !! stable inverse-survival callback
        integer, intent(out), optional :: status !! distribution_status_* result code
        type(continuous_distribution) :: dist

        logical :: finite_or_infinite

        finite_or_infinite = .not. ieee_is_nan(support(1)) .and. .not. ieee_is_nan(support(2))
        if (.not. finite_or_infinite .or. support(1) >= support(2)) then
            if (present(status)) status = distribution_status_invalid_support
            return
        end if

        dist%lower = support(1)
        dist%upper = support(2)
        if (present(inclusive)) dist%inclusive = inclusive
        if (present(parameters)) then
            allocate(dist%parameters(size(parameters)))
            dist%parameters = parameters
        else
            allocate(dist%parameters(0))
        end if

        dist%pdf_callback => pdf
        dist%cdf_callback => cdf
        dist%ppf_callback => ppf
        if (present(logpdf)) dist%logpdf_callback => logpdf
        if (present(logcdf)) dist%logcdf_callback => logcdf
        if (present(sf)) dist%sf_callback => sf
        if (present(logsf)) dist%logsf_callback => logsf
        if (present(isf)) dist%isf_callback => isf
        dist%valid = .true.
        if (present(status)) status = distribution_status_success
    end function make_distribution

    function truncate(base, lb, ub, status) result(dist)
        type(continuous_distribution), intent(in) :: base !! distribution to truncate
        real(dp), intent(in), optional :: lb !! lower truncation point, default -infinity
        real(dp), intent(in), optional :: ub !! upper truncation point, default +infinity
        integer, intent(out), optional :: status !! distribution_status_* result code
        type(truncated_distribution) :: dist

        real(dp) :: fa
        real(dp) :: lower_bound
        real(dp) :: sb
        real(dp) :: upper_bound

        lower_bound = negative_infinity(0.0_dp)
        upper_bound = positive_infinity(0.0_dp)
        if (present(lb)) lower_bound = lb
        if (present(ub)) upper_bound = ub

        if (.not. base%valid) then
            if (present(status)) status = distribution_status_invalid_support
            return
        end if
        if (ieee_is_nan(lower_bound) .or. ieee_is_nan(upper_bound)) then
            if (present(status)) status = distribution_status_invalid_bounds
            return
        end if
        if (lower_bound >= upper_bound) then
            if (present(status)) status = distribution_status_invalid_bounds
            return
        end if

        dist%base = base
        dist%lb = lower_bound
        dist%ub = upper_bound
        dist%lower = max(base%lower, lower_bound)
        dist%upper = min(base%upper, upper_bound)
        if (dist%lower >= dist%upper) then
            if (present(status)) status = distribution_status_zero_mass
            return
        end if

        fa = base%cdf(dist%lower)
        sb = base%sf(dist%upper)
        if (ieee_is_nan(fa) .or. ieee_is_nan(sb)) then
            if (present(status)) status = distribution_status_zero_mass
            return
        end if

        if (fa <= sb) then
            dist%logmass = logdiffexp_pair(base%logcdf(dist%upper), &
                base%logcdf(dist%lower))
        else
            dist%logmass = logdiffexp_pair(base%logsf(dist%lower), &
                base%logsf(dist%upper))
        end if
        if (.not. ieee_is_finite(dist%logmass)) then
            if (present(status)) status = distribution_status_zero_mass
            return
        end if
        dist%logmass = min(0.0_dp, dist%logmass)
        dist%valid = .true.
        if (present(status)) status = distribution_status_success
    end function truncate

    pure elemental function continuous_pdf(self, x) result(y)
        class(continuous_distribution), intent(in) :: self !! frozen distribution
        real(dp), intent(in) :: x !! evaluation point
        real(dp) :: y

        if (.not. self%valid .or. ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. in_density_support(x, self%lower, self%upper, self%inclusive)) then
            y = 0.0_dp
        else
            y = self%pdf_callback(x, self%parameters)
        end if
    end function continuous_pdf

    pure elemental function continuous_logpdf(self, x) result(y)
        class(continuous_distribution), intent(in) :: self !! frozen distribution
        real(dp), intent(in) :: x !! evaluation point
        real(dp) :: y

        if (.not. self%valid .or. ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. in_density_support(x, self%lower, self%upper, self%inclusive)) then
            y = negative_infinity(x)
        else if (associated(self%logpdf_callback)) then
            y = self%logpdf_callback(x, self%parameters)
        else
            y = safe_log_nonnegative(self%pdf_callback(x, self%parameters), x)
        end if
    end function continuous_logpdf

    pure elemental function continuous_cdf(self, x) result(y)
        class(continuous_distribution), intent(in) :: self !! frozen distribution
        real(dp), intent(in) :: x !! evaluation point
        real(dp) :: y

        if (.not. self%valid .or. ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= self%lower) then
            y = 0.0_dp
        else if (x >= self%upper) then
            y = 1.0_dp
        else
            y = clamp_probability(self%cdf_callback(x, self%parameters))
        end if
    end function continuous_cdf

    pure elemental function continuous_sf(self, x) result(y)
        class(continuous_distribution), intent(in) :: self !! frozen distribution
        real(dp), intent(in) :: x !! evaluation point
        real(dp) :: y

        if (.not. self%valid .or. ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= self%lower) then
            y = 1.0_dp
        else if (x >= self%upper) then
            y = 0.0_dp
        else if (associated(self%sf_callback)) then
            y = clamp_probability(self%sf_callback(x, self%parameters))
        else if (associated(self%logsf_callback)) then
            y = safe_exp_log(self%logsf_callback(x, self%parameters), x)
        else
            y = max(0.0_dp, 1.0_dp - self%cdf_callback(x, self%parameters))
        end if
    end function continuous_sf

    pure elemental function continuous_logcdf(self, x) result(y)
        class(continuous_distribution), intent(in) :: self !! frozen distribution
        real(dp), intent(in) :: x !! evaluation point
        real(dp) :: y

        if (.not. self%valid .or. ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= self%lower) then
            y = negative_infinity(x)
        else if (x >= self%upper) then
            y = 0.0_dp
        else if (associated(self%logcdf_callback)) then
            y = self%logcdf_callback(x, self%parameters)
        else
            y = safe_log_probability(self%cdf_callback(x, self%parameters), x)
        end if
    end function continuous_logcdf

    pure elemental function continuous_logsf(self, x) result(y)
        class(continuous_distribution), intent(in) :: self !! frozen distribution
        real(dp), intent(in) :: x !! evaluation point
        real(dp) :: y

        if (.not. self%valid .or. ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= self%lower) then
            y = 0.0_dp
        else if (x >= self%upper) then
            y = negative_infinity(x)
        else if (associated(self%logsf_callback)) then
            y = self%logsf_callback(x, self%parameters)
        else if (associated(self%sf_callback)) then
            y = safe_log_probability(self%sf_callback(x, self%parameters), x)
        else
            y = safe_log_probability(max(0.0_dp, &
                1.0_dp - self%cdf_callback(x, self%parameters)), x)
        end if
    end function continuous_logsf

    pure elemental function continuous_ppf(self, p) result(x)
        class(continuous_distribution), intent(in) :: self !! frozen distribution
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp) :: x

        if (.not. self%valid .or. ieee_is_nan(p)) then
            x = quiet_nan(p)
        else if (p < 0.0_dp .or. p > 1.0_dp) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = self%lower
        else if (p == 1.0_dp) then
            x = self%upper
        else
            x = self%ppf_callback(p, self%parameters)
        end if
    end function continuous_ppf

    pure elemental function continuous_isf(self, p) result(x)
        class(continuous_distribution), intent(in) :: self !! frozen distribution
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp) :: x

        if (.not. self%valid .or. ieee_is_nan(p)) then
            x = quiet_nan(p)
        else if (p < 0.0_dp .or. p > 1.0_dp) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = self%upper
        else if (p == 1.0_dp) then
            x = self%lower
        else if (associated(self%isf_callback)) then
            x = self%isf_callback(p, self%parameters)
        else
            x = self%ppf_callback(1.0_dp - p, self%parameters)
        end if
    end function continuous_isf

    pure elemental function continuous_ilogcdf(self, logp) result(x)
        class(continuous_distribution), intent(in) :: self !! frozen distribution
        real(dp), intent(in) :: logp !! logarithm of a lower-tail probability
        real(dp) :: x

        x = self%ppf(probability_from_log(logp))
    end function continuous_ilogcdf

    pure elemental function continuous_ilogsf(self, logp) result(x)
        class(continuous_distribution), intent(in) :: self !! frozen distribution
        real(dp), intent(in) :: logp !! logarithm of an upper-tail probability
        real(dp) :: x

        x = self%isf(probability_from_log(logp))
    end function continuous_ilogsf

    pure elemental function continuous_median(self) result(x)
        class(continuous_distribution), intent(in) :: self !! frozen distribution
        real(dp) :: x

        x = self%ppf(0.5_dp)
    end function continuous_median

    pure function continuous_support(self) result(bounds)
        class(continuous_distribution), intent(in) :: self !! frozen distribution
        real(dp) :: bounds(2)

        if (self%valid) then
            bounds = [self%lower, self%upper]
        else
            bounds = quiet_nan(0.0_dp)
        end if
    end function continuous_support

    function continuous_sample(self, n, state) result(values)
        class(continuous_distribution), intent(in) :: self !! frozen distribution
        integer, intent(in) :: n !! number of variates; negative values produce an empty array
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by sampling
        real(dp), allocatable :: values(:)

        integer :: i

        allocate(values(max(0, n)))
        do i = 1, size(values)
            values(i) = self%ppf(rng_uniform(state))
        end do
    end function continuous_sample

    pure elemental function truncated_pdf(self, x) result(y)
        class(truncated_distribution), intent(in) :: self !! truncated distribution
        real(dp), intent(in) :: x !! evaluation point
        real(dp) :: y
        real(dp) :: logvalue

        if (.not. self%valid .or. ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. in_density_support(x, self%lower, self%upper, self%base%inclusive)) then
            y = 0.0_dp
        else
            logvalue = self%base%logpdf(x) - self%logmass
            y = safe_exp_log(logvalue, x)
        end if
    end function truncated_pdf

    pure elemental function truncated_logpdf(self, x) result(y)
        class(truncated_distribution), intent(in) :: self !! truncated distribution
        real(dp), intent(in) :: x !! evaluation point
        real(dp) :: y

        if (.not. self%valid .or. ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. in_density_support(x, self%lower, self%upper, self%base%inclusive)) then
            y = negative_infinity(x)
        else
            y = self%base%logpdf(x) - self%logmass
        end if
    end function truncated_logpdf

    pure elemental function truncated_cdf(self, x) result(y)
        class(truncated_distribution), intent(in) :: self !! truncated distribution
        real(dp), intent(in) :: x !! evaluation point
        real(dp) :: y
        real(dp) :: lognum

        if (.not. self%valid .or. ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= self%lower) then
            y = 0.0_dp
        else if (x >= self%upper) then
            y = 1.0_dp
        else
            lognum = logdiffexp_pair(self%base%logcdf(x), self%base%logcdf(self%lower))
            y = min(1.0_dp, safe_exp_log(lognum - self%logmass, x))
        end if
    end function truncated_cdf

    pure elemental function truncated_sf(self, x) result(y)
        class(truncated_distribution), intent(in) :: self !! truncated distribution
        real(dp), intent(in) :: x !! evaluation point
        real(dp) :: y
        real(dp) :: lognum

        if (.not. self%valid .or. ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= self%lower) then
            y = 1.0_dp
        else if (x >= self%upper) then
            y = 0.0_dp
        else
            lognum = logdiffexp_pair(self%base%logsf(x), self%base%logsf(self%upper))
            y = min(1.0_dp, safe_exp_log(lognum - self%logmass, x))
        end if
    end function truncated_sf

    pure elemental function truncated_logcdf(self, x) result(y)
        class(truncated_distribution), intent(in) :: self !! truncated distribution
        real(dp), intent(in) :: x !! evaluation point
        real(dp) :: y

        if (.not. self%valid .or. ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= self%lower) then
            y = negative_infinity(x)
        else if (x >= self%upper) then
            y = 0.0_dp
        else
            y = logdiffexp_pair(self%base%logcdf(x), self%base%logcdf(self%lower)) - &
                self%logmass
            y = min(0.0_dp, y)
        end if
    end function truncated_logcdf

    pure elemental function truncated_logsf(self, x) result(y)
        class(truncated_distribution), intent(in) :: self !! truncated distribution
        real(dp), intent(in) :: x !! evaluation point
        real(dp) :: y

        if (.not. self%valid .or. ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= self%lower) then
            y = 0.0_dp
        else if (x >= self%upper) then
            y = negative_infinity(x)
        else
            y = logdiffexp_pair(self%base%logsf(x), self%base%logsf(self%upper)) - &
                self%logmass
            y = min(0.0_dp, y)
        end if
    end function truncated_logsf

    pure elemental function truncated_ppf(self, p) result(x)
        class(truncated_distribution), intent(in) :: self !! truncated distribution
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp) :: x

        real(dp) :: fa
        real(dp) :: mass
        real(dp) :: q
        real(dp) :: sb

        if (.not. self%valid .or. ieee_is_nan(p)) then
            x = quiet_nan(p)
        else if (p < 0.0_dp .or. p > 1.0_dp) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = self%lower
        else if (p == 1.0_dp) then
            x = self%upper
        else
            fa = self%base%cdf(self%lower)
            sb = self%base%sf(self%upper)
            mass = safe_exp_log(self%logmass, p)
            if (fa <= sb) then
                q = min(1.0_dp, fa + p * mass)
                x = self%base%ppf(q)
            else
                q = min(1.0_dp, sb + (1.0_dp - p) * mass)
                x = self%base%isf(q)
            end if
        end if
    end function truncated_ppf

    pure elemental function truncated_isf(self, p) result(x)
        class(truncated_distribution), intent(in) :: self !! truncated distribution
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp) :: x

        real(dp) :: fa
        real(dp) :: mass
        real(dp) :: q
        real(dp) :: sb

        if (.not. self%valid .or. ieee_is_nan(p)) then
            x = quiet_nan(p)
        else if (p < 0.0_dp .or. p > 1.0_dp) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = self%upper
        else if (p == 1.0_dp) then
            x = self%lower
        else
            fa = self%base%cdf(self%lower)
            sb = self%base%sf(self%upper)
            mass = safe_exp_log(self%logmass, p)
            if (sb < fa) then
                q = min(1.0_dp, sb + p * mass)
                x = self%base%isf(q)
            else
                q = min(1.0_dp, fa + (1.0_dp - p) * mass)
                x = self%base%ppf(q)
            end if
        end if
    end function truncated_isf

    pure elemental function truncated_ilogcdf(self, logp) result(x)
        class(truncated_distribution), intent(in) :: self !! truncated distribution
        real(dp), intent(in) :: logp !! logarithm of a lower-tail probability
        real(dp) :: x

        x = self%ppf(probability_from_log(logp))
    end function truncated_ilogcdf

    pure elemental function truncated_ilogsf(self, logp) result(x)
        class(truncated_distribution), intent(in) :: self !! truncated distribution
        real(dp), intent(in) :: logp !! logarithm of an upper-tail probability
        real(dp) :: x

        x = self%isf(probability_from_log(logp))
    end function truncated_ilogsf

    pure elemental function truncated_median(self) result(x)
        class(truncated_distribution), intent(in) :: self !! truncated distribution
        real(dp) :: x

        x = self%ppf(0.5_dp)
    end function truncated_median

    pure function truncated_support(self) result(bounds)
        class(truncated_distribution), intent(in) :: self !! truncated distribution
        real(dp) :: bounds(2)

        if (self%valid) then
            bounds = [self%lower, self%upper]
        else
            bounds = quiet_nan(0.0_dp)
        end if
    end function truncated_support

    function truncated_sample(self, n, state) result(values)
        class(truncated_distribution), intent(in) :: self !! truncated distribution
        integer, intent(in) :: n !! number of variates; negative values produce an empty array
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by sampling
        real(dp), allocatable :: values(:)

        integer :: i

        allocate(values(max(0, n)))
        do i = 1, size(values)
            values(i) = self%ppf(rng_uniform(state))
        end do
    end function truncated_sample

    pure logical function in_density_support(x, lower, upper, inclusive) result(inside)
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: lower !! lower support endpoint
        real(dp), intent(in) :: upper !! upper support endpoint
        logical, intent(in) :: inclusive(2) !! endpoint inclusion flags

        logical :: left
        logical :: right

        if (inclusive(1)) then
            left = x >= lower
        else
            left = x > lower
        end if
        if (inclusive(2)) then
            right = x <= upper
        else
            right = x < upper
        end if
        inside = left .and. right
    end function in_density_support

    pure elemental function safe_log_nonnegative(value, seed) result(y)
        real(dp), intent(in) :: value !! nonnegative quantity
        real(dp), intent(in) :: seed !! NaN/infinity construction seed
        real(dp) :: y

        if (ieee_is_nan(value) .or. value < 0.0_dp) then
            y = quiet_nan(seed)
        else if (value == 0.0_dp) then
            y = negative_infinity(seed)
        else
            y = log(value)
        end if
    end function safe_log_nonnegative

    pure elemental function safe_log_probability(value, seed) result(y)
        real(dp), intent(in) :: value !! probability-like quantity
        real(dp), intent(in) :: seed !! NaN/infinity construction seed
        real(dp) :: y

        if (ieee_is_nan(value) .or. value < 0.0_dp .or. value > 1.0_dp) then
            y = quiet_nan(seed)
        else if (value == 0.0_dp) then
            y = negative_infinity(seed)
        else
            y = log(value)
        end if
    end function safe_log_probability

    pure elemental function safe_exp_log(logvalue, seed) result(value)
        real(dp), intent(in) :: logvalue !! logarithm of a nonnegative quantity
        real(dp), intent(in) :: seed !! NaN/infinity construction seed
        real(dp) :: value

        if (ieee_is_nan(logvalue)) then
            value = quiet_nan(seed)
        else if (logvalue == negative_infinity(logvalue)) then
            value = 0.0_dp
        else if (logvalue > log(huge(1.0_dp))) then
            value = positive_infinity(seed)
        else
            value = exp(logvalue)
        end if
    end function safe_exp_log

    pure elemental function clamp_probability(value) result(y)
        real(dp), intent(in) :: value !! probability-like value
        real(dp) :: y

        if (ieee_is_nan(value)) then
            y = value
        else
            y = min(1.0_dp, max(0.0_dp, value))
        end if
    end function clamp_probability

    pure elemental function logdiffexp_pair(a, b) result(y)
        real(dp), intent(in) :: a !! logarithm of the leading nonnegative term
        real(dp), intent(in) :: b !! logarithm of the subtracted nonnegative term
        real(dp) :: y

        if (ieee_is_nan(a) .or. ieee_is_nan(b)) then
            y = quiet_nan(a + b)
        else if (a == negative_infinity(a)) then
            y = negative_infinity(a)
        else if (b == negative_infinity(b)) then
            y = a
        else if (b >= a) then
            y = negative_infinity(a)
        else
            y = a + log1p_safe(-exp(b - a))
        end if
    end function logdiffexp_pair

    pure elemental function probability_from_log(logp) result(p)
        real(dp), intent(in) :: logp !! logarithm of a probability
        real(dp) :: p

        if (ieee_is_nan(logp) .or. logp > 0.0_dp) then
            p = quiet_nan(logp)
        else if (logp == negative_infinity(logp)) then
            p = 0.0_dp
        else
            p = exp(logp)
        end if
    end function probability_from_log

end module scifort_distribution_infrastructure
