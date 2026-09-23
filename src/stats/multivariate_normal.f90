! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_multivariate_normal
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_log_sqrt_two_pi
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, psd_decompose
    use scifort_multivariate_integration, only : qmv_status_success, qmvn_box
    use scifort_math, only : expm1_safe, negative_infinity, quiet_nan
    use scifort_normal, only : normal_logcdf, normal_logsf, normal_ppf
    use scifort_random, only : rng_state, rng_uniform
    implicit none
    private

    integer, parameter, public :: multivariate_status_success = 0
    integer, parameter, public :: multivariate_status_invalid_shape = 1
    integer, parameter, public :: multivariate_status_invalid_parameter = 2
    integer, parameter, public :: multivariate_status_linalg_failure = 3

    public :: multivariate_normal_cdf
    public :: multivariate_normal_entropy
    public :: multivariate_normal_fit
    public :: multivariate_normal_logcdf
    public :: multivariate_normal_logpdf
    public :: multivariate_normal_pdf
    public :: multivariate_normal_marginal
    public :: multivariate_normal_rvs
    public :: multivariate_normal_rvs_array

contains

    function multivariate_normal_logpdf(x, mean, cov, allow_singular) result(y)
        real(dp), intent(in) :: x(:) !! point at which the density is evaluated
        real(dp), intent(in), optional :: mean(:) !! mean vector; zero vector when absent
        real(dp), intent(in), optional :: cov(:, :) !! covariance matrix; identity when absent
        logical, intent(in), optional :: allow_singular !! permit positive-semidefinite covariance (default false)
        real(dp) :: y

        integer :: i
        integer :: n
        integer :: rank
        integer :: status
        logical :: singular_ok
        real(dp), allocatable :: covariance(:, :)
        real(dp), allocatable :: dev(:)
        real(dp), allocatable :: eigenvalues(:)
        real(dp), allocatable :: eigenvectors(:, :)
        real(dp), allocatable :: mu(:)
        real(dp) :: cutoff
        real(dp) :: log_pdet
        real(dp) :: maha
        real(dp) :: projection
        real(dp) :: residual2
        real(dp) :: support_eps

        n = size(x)
        if (n < 1) then
            y = quiet_nan(0.0_dp)
            return
        end if
        allocate(mu(n), covariance(n, n), dev(n), eigenvalues(n), eigenvectors(n, n))
        call process_mean_cov(n, mean, cov, mu, covariance, status)
        if (status /= multivariate_status_success) then
            y = quiet_nan(x(1))
            return
        end if

        singular_ok = .false.
        if (present(allow_singular)) singular_ok = allow_singular
        call psd_decompose(covariance, eigenvalues, eigenvectors, rank, log_pdet, cutoff, &
            status, singular_ok)
        if (status /= linalg_status_success) then
            y = quiet_nan(x(1))
            return
        end if

        dev = x - mu
        if (any(ieee_is_nan(dev))) then
            y = quiet_nan(x(1))
            return
        end if

        maha = 0.0_dp
        residual2 = 0.0_dp
        do i = 1, n
            projection = dot_product(dev, eigenvectors(:, i))
            if (eigenvalues(i) > cutoff) then
                maha = maha + projection * projection / eigenvalues(i)
            else
                residual2 = residual2 + projection * projection
            end if
        end do

        if (rank < n) then
            support_eps = 1.0e3_dp * cutoff
            if (sqrt(residual2) >= support_eps) then
                y = negative_infinity(x(1))
                return
            end if
        end if

        y = -0.5_dp * (real(rank, dp) * (2.0_dp * scifort_log_sqrt_two_pi) + &
            log_pdet + maha)
    end function multivariate_normal_logpdf

    function multivariate_normal_pdf(x, mean, cov, allow_singular) result(y)
        real(dp), intent(in) :: x(:) !! point at which the density is evaluated
        real(dp), intent(in), optional :: mean(:) !! mean vector; zero vector when absent
        real(dp), intent(in), optional :: cov(:, :) !! covariance matrix; identity when absent
        logical, intent(in), optional :: allow_singular !! permit positive-semidefinite covariance (default false)
        real(dp) :: y

        real(dp) :: logp

        if (present(mean)) then
            if (present(cov)) then
                if (present(allow_singular)) then
                    logp = multivariate_normal_logpdf(x, mean, cov, allow_singular)
                else
                    logp = multivariate_normal_logpdf(x, mean, cov)
                end if
            else if (present(allow_singular)) then
                logp = multivariate_normal_logpdf(x, mean=mean, allow_singular=allow_singular)
            else
                logp = multivariate_normal_logpdf(x, mean=mean)
            end if
        else if (present(cov)) then
            if (present(allow_singular)) then
                logp = multivariate_normal_logpdf(x, cov=cov, allow_singular=allow_singular)
            else
                logp = multivariate_normal_logpdf(x, cov=cov)
            end if
        else
            logp = multivariate_normal_logpdf(x)
        end if
        y = exp(logp)
    end function multivariate_normal_pdf

    function multivariate_normal_cdf(x, mean, cov, allow_singular, maxpts, abseps, releps, &
            lower_limit, state) result(y)
        real(dp), intent(in) :: x(:) !! upper integration limits; last-axis point in SciPy terminology
        real(dp), intent(in), optional :: mean(:) !! mean vector; zero vector when absent
        real(dp), intent(in), optional :: cov(:, :) !! covariance matrix; identity when absent
        logical, intent(in), optional :: allow_singular !! permit positive-semidefinite covariance (default false)
        integer, intent(in), optional :: maxpts !! approximate maximum QMC points (default 1000000*dimension)
        real(dp), intent(in), optional :: abseps !! requested absolute integration error (default 1e-5)
        real(dp), intent(in), optional :: releps !! SciPy-compatible relative tolerance argument; currently informational
        real(dp), intent(in), optional :: lower_limit(:) !! lower integration limits; -infinity when absent
        type(rng_state), intent(inout), optional :: state !! optional explicit RNG state for randomized QMC shifts
        real(dp) :: y

        integer :: i
        integer :: info
        integer :: n
        integer :: nsamples
        integer :: rank
        integer :: sign_value
        logical :: singular_ok
        real(dp), allocatable :: covariance(:, :)
        real(dp), allocatable :: eigenvalues(:)
        real(dp), allocatable :: eigenvectors(:, :)
        real(dp), allocatable :: high(:)
        real(dp), allocatable :: low(:)
        real(dp), allocatable :: mu(:)
        real(dp) :: cutoff
        real(dp) :: estimated_error
        real(dp) :: log_pdet
        real(dp) :: tmp

        n = size(x)
        if (n < 1 .or. any(ieee_is_nan(x))) then
            y = quiet_nan(0.0_dp)
            return
        end if
        if (present(releps)) then
            if (ieee_is_nan(releps) .or. releps <= 0.0_dp) then
                y = quiet_nan(0.0_dp)
                return
            end if
        end if

        allocate(mu(n), covariance(n, n), eigenvalues(n), eigenvectors(n, n), high(n), low(n))
        call process_mean_cov(n, mean, cov, mu, covariance, info)
        if (info /= multivariate_status_success) then
            y = quiet_nan(0.0_dp)
            return
        end if

        singular_ok = .false.
        if (present(allow_singular)) singular_ok = allow_singular
        call psd_decompose(covariance, eigenvalues, eigenvectors, rank, log_pdet, cutoff, &
            info, singular_ok)
        if (info /= linalg_status_success) then
            y = quiet_nan(0.0_dp)
            return
        end if

        high = x - mu
        if (present(lower_limit)) then
            if (size(lower_limit) /= n .or. any(ieee_is_nan(lower_limit))) then
                y = quiet_nan(0.0_dp)
                return
            end if
            low = lower_limit - mu
        else
            low = negative_infinity(0.0_dp)
        end if

        sign_value = 1
        do i = 1, n
            if (high(i) < low(i)) then
                tmp = high(i)
                high(i) = low(i)
                low(i) = tmp
                sign_value = -sign_value
            end if
        end do

        if (present(state)) then
            call qmvn_box(covariance, low, high, y, estimated_error, nsamples, info, &
                state=state, maxpts=maxpts, abseps=abseps)
        else
            call qmvn_box(covariance, low, high, y, estimated_error, nsamples, info, &
                maxpts=maxpts, abseps=abseps)
        end if
        if (info /= qmv_status_success) then
            y = quiet_nan(0.0_dp)
        else
            y = real(sign_value, dp) * y
        end if
    end function multivariate_normal_cdf

    function multivariate_normal_logcdf(x, mean, cov, allow_singular, maxpts, abseps, releps, &
            lower_limit, state) result(y)
        real(dp), intent(in) :: x(:) !! upper integration limits
        real(dp), intent(in), optional :: mean(:) !! mean vector; zero vector when absent
        real(dp), intent(in), optional :: cov(:, :) !! covariance matrix; identity when absent
        logical, intent(in), optional :: allow_singular !! permit positive-semidefinite covariance (default false)
        integer, intent(in), optional :: maxpts !! approximate maximum QMC points
        real(dp), intent(in), optional :: abseps !! requested absolute integration error
        real(dp), intent(in), optional :: releps !! SciPy-compatible relative tolerance argument
        real(dp), intent(in), optional :: lower_limit(:) !! lower integration limits
        type(rng_state), intent(inout), optional :: state !! optional explicit RNG state for randomized QMC shifts
        real(dp) :: y

        real(dp) :: p, cutoff, log_pdet, lower, upper, sd, tmp
        real(dp), allocatable :: mu(:), covariance(:, :), eigenvalues(:), eigenvectors(:, :)
        integer :: n, i, j, info, rank, sign_value
        logical :: diagonal, singular_ok

        n = size(x)
        y = quiet_nan(0.0_dp)
        if (n < 1 .or. any(ieee_is_nan(x))) return
        if (present(abseps)) then
            if (.not. ieee_is_finite(abseps) .or. abseps <= 0.0_dp) return
        end if
        if (present(releps)) then
            if (ieee_is_nan(releps) .or. releps <= 0.0_dp) return
        end if
        if (present(lower_limit)) then
            if (size(lower_limit) /= n .or. any(ieee_is_nan(lower_limit))) return
        end if
        allocate(mu(n), covariance(n,n), eigenvalues(n), eigenvectors(n,n))
        call process_mean_cov(n, mean, cov, mu, covariance, info)
        if (info /= multivariate_status_success) return
        singular_ok = .false.
        if (present(allow_singular)) singular_ok = allow_singular
        call psd_decompose(covariance, eigenvalues, eigenvectors, rank, log_pdet, cutoff, &
            info, singular_ok)
        if (info /= linalg_status_success) return
        diagonal = .true.
        do i = 1, n
            do j = 1, i - 1
                if (covariance(i,j) /= 0.0_dp) diagonal = .false.
            end do
        end do

        if (diagonal) then
            ! Independent coordinates factor exactly; never form the tiny product.
            y = 0.0_dp
            sign_value = 1
            do i = 1, n
                upper = x(i) - mu(i)
                lower = negative_infinity(0.0_dp)
                if (present(lower_limit)) lower = lower_limit(i) - mu(i)
                if (upper < lower) then
                    tmp = lower
                    lower = upper
                    upper = tmp
                    sign_value = -sign_value
                end if
                if (covariance(i,i) <= 0.0_dp) then
                    if (lower > 0.0_dp .or. upper < 0.0_dp) y = negative_infinity(0.0_dp)
                else
                    sd = sqrt(covariance(i,i))
                    y = y + normal_log_interval(lower / sd, upper / sd)
                end if
            end do
            if (sign_value < 0 .and. y /= negative_infinity(0.0_dp)) y = quiet_nan(0.0_dp)
            return
        end if

        ! Correlated boxes retain the existing probability-domain QMC integration.
        p = multivariate_normal_cdf(x, mean, cov, allow_singular, maxpts, abseps, releps, &
            lower_limit, state)
        if (ieee_is_nan(p) .or. p < 0.0_dp) then
            y = quiet_nan(0.0_dp)
        else if (p == 0.0_dp) then
            y = negative_infinity(0.0_dp)
        else
            y = log(p)
        end if
    end function multivariate_normal_logcdf

    pure function normal_log_interval(lower, upper) result(value)
        implicit none
        real(dp), intent(in) :: lower !! standardized lower endpoint, <= upper
        real(dp), intent(in) :: upper !! standardized upper endpoint
        real(dp) :: value, larger, smaller

        if (lower == upper) then
            value = negative_infinity(0.0_dp)
            return
        end if
        if (lower >= 0.0_dp) then
            larger = normal_logsf(lower)
            smaller = normal_logsf(upper)
        else
            larger = normal_logcdf(upper)
            smaller = normal_logcdf(lower)
        end if
        if (larger == negative_infinity(0.0_dp)) then
            value = larger
        else
            value = larger + log(-expm1_safe(smaller - larger))
        end if
    end function normal_log_interval

    subroutine multivariate_normal_marginal(dimensions, mean, cov, marginal_mean, marginal_cov, status)
        integer, intent(in) :: dimensions(:) !! one-based indices of retained variables; must be unique
        real(dp), intent(in) :: mean(:) !! original mean vector
        real(dp), intent(in) :: cov(:, :) !! original covariance matrix; lower triangle is authoritative
        real(dp), intent(out) :: marginal_mean(:) !! selected mean vector, in requested order
        real(dp), intent(out) :: marginal_cov(:, :) !! selected covariance matrix, in requested order
        integer, intent(out) :: status !! zero on success; nonzero for invalid dimensions or shapes

        integer :: i
        integer :: j
        integer :: n
        integer :: m

        n = size(mean)
        m = size(dimensions)
        if (n < 1 .or. m < 1 .or. size(cov, 1) /= n .or. size(cov, 2) /= n .or. &
                size(marginal_mean) /= m .or. size(marginal_cov, 1) /= m .or. &
                size(marginal_cov, 2) /= m) then
            status = multivariate_status_invalid_shape
            return
        end if
        if (any(.not. ieee_is_finite(mean)) .or. any(.not. ieee_is_finite(cov)) .or. &
                any(dimensions < 1) .or. any(dimensions > n)) then
            status = multivariate_status_invalid_parameter
            return
        end if
        do i = 1, m
            do j = 1, i - 1
                if (dimensions(i) == dimensions(j)) then
                    status = multivariate_status_invalid_parameter
                    return
                end if
            end do
        end do

        do i = 1, m
            marginal_mean(i) = mean(dimensions(i))
            do j = 1, m
                if (dimensions(i) >= dimensions(j)) then
                    marginal_cov(i, j) = cov(dimensions(i), dimensions(j))
                else
                    marginal_cov(i, j) = cov(dimensions(j), dimensions(i))
                end if
            end do
        end do
        status = multivariate_status_success
    end subroutine multivariate_normal_marginal

    function multivariate_normal_entropy(mean, cov) result(h)
        real(dp), intent(in), optional :: mean(:) !! mean vector used only to infer dimension
        real(dp), intent(in), optional :: cov(:, :) !! covariance matrix; identity when absent
        real(dp) :: h

        integer :: n
        integer :: rank
        integer :: status
        real(dp), allocatable :: covariance(:, :)
        real(dp), allocatable :: eigenvalues(:)
        real(dp), allocatable :: eigenvectors(:, :)
        real(dp), allocatable :: mu(:)
        real(dp) :: cutoff
        real(dp) :: log_pdet

        if (present(mean)) then
            n = size(mean)
        else if (present(cov)) then
            n = size(cov, 1)
        else
            n = 1
        end if
        if (n < 1) then
            h = quiet_nan(0.0_dp)
            return
        end if
        allocate(mu(n), covariance(n, n), eigenvalues(n), eigenvectors(n, n))
        call process_mean_cov(n, mean, cov, mu, covariance, status)
        if (status /= multivariate_status_success) then
            h = quiet_nan(0.0_dp)
            return
        end if

        call psd_decompose(covariance, eigenvalues, eigenvectors, rank, log_pdet, cutoff, &
            status, .true.)
        if (status /= linalg_status_success) then
            h = quiet_nan(0.0_dp)
            return
        end if
        h = 0.5_dp * (real(rank, dp) * (2.0_dp * scifort_log_sqrt_two_pi + 1.0_dp) + &
            log_pdet)
    end function multivariate_normal_entropy

    subroutine multivariate_normal_rvs(state, sample, mean, cov, status)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by the draw
        real(dp), intent(out) :: sample(:) !! generated vector; its length defines the dimension
        real(dp), intent(in), optional :: mean(:) !! mean vector; zero vector when absent
        real(dp), intent(in), optional :: cov(:, :) !! covariance matrix; identity when absent
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid covariance or shape

        integer :: i
        integer :: info
        integer :: j
        integer :: n
        integer :: rank
        real(dp), allocatable :: covariance(:, :)
        real(dp), allocatable :: eigenvalues(:)
        real(dp), allocatable :: eigenvectors(:, :)
        real(dp), allocatable :: mu(:)
        real(dp), allocatable :: z(:)
        real(dp) :: cutoff
        real(dp) :: log_pdet
        real(dp) :: u

        n = size(sample)
        if (n < 1) then
            if (present(status)) status = multivariate_status_invalid_shape
            return
        end if
        allocate(mu(n), covariance(n, n), eigenvalues(n), eigenvectors(n, n), z(n))
        call process_mean_cov(n, mean, cov, mu, covariance, info)
        if (info /= multivariate_status_success) then
            sample = quiet_nan(0.0_dp)
            if (present(status)) status = info
            return
        end if
        call psd_decompose(covariance, eigenvalues, eigenvectors, rank, log_pdet, cutoff, &
            info, .true.)
        if (info /= linalg_status_success) then
            sample = quiet_nan(0.0_dp)
            if (present(status)) status = multivariate_status_linalg_failure
            return
        end if

        do i = 1, n
            u = rng_uniform(state)
            z(i) = normal_ppf(u)
        end do
        sample = mu
        do j = 1, n
            if (eigenvalues(j) > 0.0_dp) then
                sample = sample + eigenvectors(:, j) * sqrt(eigenvalues(j)) * z(j)
            end if
        end do
        if (present(status)) status = multivariate_status_success
    end subroutine multivariate_normal_rvs

    subroutine multivariate_normal_rvs_array(state, samples, mean, cov, status)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by all draws
        real(dp), intent(out) :: samples(:, :) !! samples by column, shape (dimension, n_samples)
        real(dp), intent(in), optional :: mean(:) !! mean vector; zero vector when absent
        real(dp), intent(in), optional :: cov(:, :) !! covariance matrix; identity when absent
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid covariance or shape

        integer :: i
        integer :: info
        integer :: j
        integer :: k
        integer :: n
        integer :: rank
        real(dp), allocatable :: covariance(:, :)
        real(dp), allocatable :: eigenvalues(:)
        real(dp), allocatable :: eigenvectors(:, :)
        real(dp), allocatable :: mu(:)
        real(dp), allocatable :: z(:)
        real(dp) :: cutoff
        real(dp) :: log_pdet

        n = size(samples, 1)
        if (n < 1) then
            if (present(status)) status = multivariate_status_invalid_shape
            return
        end if
        allocate(mu(n), covariance(n,n), eigenvalues(n), eigenvectors(n,n), z(n))
        call process_mean_cov(n, mean, cov, mu, covariance, info)
        if (info /= multivariate_status_success) then
            samples = quiet_nan(0.0_dp)
            if (present(status)) status = info
            return
        end if
        call psd_decompose(covariance, eigenvalues, eigenvectors, rank, log_pdet, cutoff, &
            info, .true.)
        if (info /= linalg_status_success) then
            samples = quiet_nan(0.0_dp)
            if (present(status)) status = multivariate_status_linalg_failure
            return
        end if

        do k = 1, size(samples, 2)
            do i = 1, n
                z(i) = normal_ppf(rng_uniform(state))
            end do
            samples(:,k) = mu
            do j = 1, n
                if (eigenvalues(j) > 0.0_dp) then
                    samples(:,k) = samples(:,k) + eigenvectors(:,j) * sqrt(eigenvalues(j)) * z(j)
                end if
            end do
        end do
        if (present(status)) status = multivariate_status_success
    end subroutine multivariate_normal_rvs_array

    subroutine multivariate_normal_fit(data, mean, cov, status, fixed_mean, fixed_cov)
        real(dp), intent(in) :: data(:, :) !! observations by column, shape (dimension, n_observations)
        real(dp), intent(out) :: mean(:) !! fitted or fixed mean vector
        real(dp), intent(out) :: cov(:, :) !! fitted or fixed maximum-likelihood covariance matrix
        integer, intent(out) :: status !! zero on success; nonzero for invalid input
        real(dp), intent(in), optional :: fixed_mean(:) !! fixed mean vector, if supplied
        real(dp), intent(in), optional :: fixed_cov(:, :) !! fixed PSD covariance matrix, if supplied

        integer :: i
        integer :: j
        integer :: n
        integer :: nobs
        integer :: rank
        integer :: lstatus
        real(dp), allocatable :: centered(:)
        real(dp), allocatable :: eigenvalues(:)
        real(dp), allocatable :: eigenvectors(:, :)
        real(dp) :: cutoff
        real(dp) :: log_pdet

        n = size(data, 1)
        nobs = size(data, 2)
        if (n < 1 .or. nobs < 1 .or. size(mean) /= n .or. size(cov, 1) /= n .or. &
                size(cov, 2) /= n .or. any(.not. ieee_is_finite(data))) then
            mean = quiet_nan(0.0_dp)
            cov = quiet_nan(0.0_dp)
            status = multivariate_status_invalid_shape
            return
        end if

        if (present(fixed_mean)) then
            if (size(fixed_mean) /= n .or. any(.not. ieee_is_finite(fixed_mean))) then
                mean = quiet_nan(0.0_dp)
                cov = quiet_nan(0.0_dp)
                status = multivariate_status_invalid_parameter
                return
            end if
            mean = fixed_mean
        else
            mean = sum(data, dim=2) / real(nobs, dp)
        end if

        if (present(fixed_cov)) then
            if (size(fixed_cov, 1) /= n .or. size(fixed_cov, 2) /= n .or. &
                    any(.not. ieee_is_finite(fixed_cov))) then
                cov = quiet_nan(0.0_dp)
                status = multivariate_status_invalid_parameter
                return
            end if
            allocate(eigenvalues(n), eigenvectors(n, n))
            call psd_decompose(fixed_cov, eigenvalues, eigenvectors, rank, log_pdet, cutoff, &
                lstatus, .true.)
            if (lstatus /= linalg_status_success) then
                cov = quiet_nan(0.0_dp)
                status = multivariate_status_invalid_parameter
                return
            end if
            call copy_lower_symmetric(fixed_cov, cov)
        else
            allocate(centered(n))
            cov = 0.0_dp
            do j = 1, nobs
                centered = data(:, j) - mean
                do i = 1, n
                    cov(i, :) = cov(i, :) + centered(i) * centered
                end do
            end do
            cov = cov / real(nobs, dp)
        end if
        status = multivariate_status_success
    end subroutine multivariate_normal_fit

    subroutine process_mean_cov(n, mean, cov, mu, covariance, status)
        integer, intent(in) :: n !! requested dimension
        real(dp), intent(in), optional :: mean(:) !! optional mean vector
        real(dp), intent(in), optional :: cov(:, :) !! optional covariance matrix
        real(dp), intent(out) :: mu(:) !! expanded mean vector
        real(dp), intent(out) :: covariance(:, :) !! expanded symmetric covariance matrix
        integer, intent(out) :: status !! zero on success; nonzero on invalid shape/value

        integer :: i

        if (size(mu) /= n .or. size(covariance, 1) /= n .or. size(covariance, 2) /= n) then
            status = multivariate_status_invalid_shape
            return
        end if

        if (present(mean)) then
            if (size(mean) /= n .or. any(.not. ieee_is_finite(mean))) then
                status = multivariate_status_invalid_parameter
                return
            end if
            mu = mean
        else
            mu = 0.0_dp
        end if

        if (present(cov)) then
            if (size(cov, 1) /= n .or. size(cov, 2) /= n) then
                status = multivariate_status_invalid_shape
                return
            end if
            if (any(.not. ieee_is_finite(cov))) then
                status = multivariate_status_invalid_parameter
                return
            end if
            call copy_lower_symmetric(cov, covariance)
        else
            covariance = 0.0_dp
            do i = 1, n
                covariance(i, i) = 1.0_dp
            end do
        end if
        status = multivariate_status_success
    end subroutine process_mean_cov

    pure subroutine copy_lower_symmetric(source, target)
        real(dp), intent(in) :: source(:, :) !! source whose lower triangle is authoritative
        real(dp), intent(out) :: target(:, :) !! symmetric copy of the lower triangle

        integer :: i
        integer :: j

        target = 0.0_dp
        do i = 1, size(target, 1)
            target(i, i) = source(i, i)
            do j = 1, i - 1
                target(i, j) = source(i, j)
                target(j, i) = source(i, j)
            end do
        end do
    end subroutine copy_lower_symmetric

end module scifort_multivariate_normal
