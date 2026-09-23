! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Multivariate Student t distribution corresponding to scipy.stats.multivariate_t.
! Dense shape matrices use their lower triangle, matching SciPy's covariance
! convention.  The shape matrix is not the covariance matrix; when df > 2 the
! covariance is df/(df-2) times shape.

module scifort_multivariate_t
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_log_sqrt_two_pi, scifort_pi
    use scifort_digamma, only : digamma_positive
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, psd_decompose
    use scifort_math, only : negative_infinity, quiet_nan
    use scifort_multivariate_integration, only : qmv_status_success, qmvt_box
    use scifort_multivariate_normal, only : multivariate_normal_entropy, &
        multivariate_normal_logpdf, multivariate_normal_rvs
    use scifort_random, only : rng_state, rng_uniform
    use scifort_chi2, only : chi2_ppf
    implicit none
    private

    integer, parameter, public :: multivariate_t_status_success = 0
    integer, parameter, public :: multivariate_t_status_invalid_shape = 1
    integer, parameter, public :: multivariate_t_status_invalid_parameter = 2
    integer, parameter, public :: multivariate_t_status_linalg_failure = 3

    public :: multivariate_t_cdf
    public :: multivariate_t_entropy
    public :: multivariate_t_logpdf
    public :: multivariate_t_marginal
    public :: multivariate_t_pdf
    public :: multivariate_t_rvs
    public :: multivariate_t_rvs_array

contains

    function multivariate_t_logpdf(x, loc, shape, df) result(y)
        real(dp), intent(in) :: x(:) !! point at which the log density is evaluated
        real(dp), intent(in), optional :: loc(:) !! location vector; zero when absent
        real(dp), intent(in), optional :: shape(:, :) !! positive-semidefinite shape matrix; identity when absent
        real(dp), intent(in), optional :: df !! degrees of freedom, > 0; default 1; +infinity gives normal
        real(dp) :: y

        integer :: info
        integer :: n
        integer :: rank
        real(dp), allocatable :: matrix(:, :)
        real(dp), allocatable :: eigenvalues(:)
        real(dp), allocatable :: eigenvectors(:, :)
        real(dp), allocatable :: mu(:)
        real(dp) :: cutoff
        real(dp) :: log_pdet
        real(dp) :: nu

        n = size(x)
        if (n < 1 .or. any(ieee_is_nan(x))) then
            y = quiet_nan(0.0_dp)
            return
        end if
        allocate(mu(n), matrix(n, n), eigenvalues(n), eigenvectors(n, n))
        call process_loc_shape(n, loc, shape, df, mu, matrix, nu, info)
        if (info /= multivariate_t_status_success) then
            y = quiet_nan(0.0_dp)
            return
        end if

        if (.not. ieee_is_finite(nu)) then
            y = multivariate_normal_logpdf(x, mean=mu, cov=matrix, allow_singular=.true.)
            return
        end if

        ! scipy.stats.multivariate_t.logpdf constructs _PSD with its default
        ! allow_singular=True, unlike pdf whose public allow_singular default is false.
        call psd_decompose(matrix, eigenvalues, eigenvectors, rank, log_pdet, cutoff, &
            info, .true.)
        if (info /= linalg_status_success) then
            y = quiet_nan(0.0_dp)
            return
        end if
        y = logpdf_core(x, mu, eigenvalues, eigenvectors, log_pdet, cutoff, nu)
    end function multivariate_t_logpdf

    function multivariate_t_pdf(x, loc, shape, df, allow_singular) result(y)
        real(dp), intent(in) :: x(:) !! point at which the density is evaluated
        real(dp), intent(in), optional :: loc(:) !! location vector; zero when absent
        real(dp), intent(in), optional :: shape(:, :) !! positive-semidefinite shape matrix; identity when absent
        real(dp), intent(in), optional :: df !! degrees of freedom, > 0; default 1; +infinity gives normal
        logical, intent(in), optional :: allow_singular !! permit singular shape matrix (default false)
        real(dp) :: y

        integer :: info
        integer :: n
        integer :: rank
        logical :: singular_ok
        real(dp), allocatable :: matrix(:, :)
        real(dp), allocatable :: eigenvalues(:)
        real(dp), allocatable :: eigenvectors(:, :)
        real(dp), allocatable :: mu(:)
        real(dp) :: cutoff
        real(dp) :: log_pdet
        real(dp) :: nu

        n = size(x)
        if (n < 1 .or. any(ieee_is_nan(x))) then
            y = quiet_nan(0.0_dp)
            return
        end if
        allocate(mu(n), matrix(n, n), eigenvalues(n), eigenvectors(n, n))
        call process_loc_shape(n, loc, shape, df, mu, matrix, nu, info)
        if (info /= multivariate_t_status_success) then
            y = quiet_nan(0.0_dp)
            return
        end if

        singular_ok = .false.
        if (present(allow_singular)) singular_ok = allow_singular
        if (.not. ieee_is_finite(nu)) then
            y = exp(multivariate_normal_logpdf(x, mean=mu, cov=matrix, allow_singular=singular_ok))
            return
        end if
        call psd_decompose(matrix, eigenvalues, eigenvectors, rank, log_pdet, cutoff, &
            info, singular_ok)
        if (info /= linalg_status_success) then
            y = quiet_nan(0.0_dp)
            return
        end if
        y = exp(logpdf_core(x, mu, eigenvalues, eigenvectors, log_pdet, cutoff, nu))
    end function multivariate_t_pdf

    function multivariate_t_cdf(x, loc, shape, df, allow_singular, maxpts, lower_limit, state) result(y)
        real(dp), intent(in) :: x(:) !! upper integration limits
        real(dp), intent(in), optional :: loc(:) !! location vector; zero when absent
        real(dp), intent(in), optional :: shape(:, :) !! positive-semidefinite shape matrix; identity when absent
        real(dp), intent(in), optional :: df !! degrees of freedom, > 0; default 1; +infinity gives normal
        logical, intent(in), optional :: allow_singular !! permit singular shape matrix (default false)
        integer, intent(in), optional :: maxpts !! approximate QMC point limit (default 1000*dimension)
        real(dp), intent(in), optional :: lower_limit(:) !! lower limits; -infinity when absent
        type(rng_state), intent(inout), optional :: state !! optional explicit RNG for randomized shifts
        real(dp) :: y

        integer :: i
        integer :: info
        integer :: limit
        integer :: n
        integer :: nsamples
        integer :: rank
        integer :: sign_value
        logical :: singular_ok
        real(dp), allocatable :: matrix(:, :)
        real(dp), allocatable :: eigenvalues(:)
        real(dp), allocatable :: eigenvectors(:, :)
        real(dp), allocatable :: high(:)
        real(dp), allocatable :: low(:)
        real(dp), allocatable :: mu(:)
        real(dp) :: cutoff
        real(dp) :: estimated_error
        real(dp) :: log_pdet
        real(dp) :: nu
        real(dp) :: tmp

        n = size(x)
        if (n < 1 .or. any(ieee_is_nan(x))) then
            y = quiet_nan(0.0_dp)
            return
        end if
        allocate(mu(n), matrix(n, n), eigenvalues(n), eigenvectors(n, n), high(n), low(n))
        call process_loc_shape(n, loc, shape, df, mu, matrix, nu, info)
        if (info /= multivariate_t_status_success) then
            y = quiet_nan(0.0_dp)
            return
        end if

        singular_ok = .false.
        if (present(allow_singular)) singular_ok = allow_singular
        call psd_decompose(matrix, eigenvalues, eigenvectors, rank, log_pdet, cutoff, &
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

        if (1000.0_dp * real(n, dp) > real(huge(1), dp)) then
            limit = huge(1)
        else
            limit = 1000 * n
        end if
        if (present(maxpts)) then
            if (maxpts > 0) limit = maxpts
        end if

        if (present(state)) then
            call qmvt_box(matrix, low, high, nu, y, estimated_error, nsamples, info, &
                state=state, maxpts=limit)
        else
            call qmvt_box(matrix, low, high, nu, y, estimated_error, nsamples, info, maxpts=limit)
        end if
        if (info /= qmv_status_success) then
            y = quiet_nan(0.0_dp)
        else
            y = real(sign_value, dp) * y
        end if
    end function multivariate_t_cdf

    function multivariate_t_entropy(loc, shape, df) result(h)
        real(dp), intent(in), optional :: loc(:) !! location vector, used to infer dimension
        real(dp), intent(in), optional :: shape(:, :) !! positive-semidefinite shape matrix; identity when absent
        real(dp), intent(in), optional :: df !! degrees of freedom, > 0; default 1; +infinity gives normal
        real(dp) :: h

        integer :: info
        integer :: n
        integer :: rank
        real(dp), allocatable :: matrix(:, :)
        real(dp), allocatable :: eigenvalues(:)
        real(dp), allocatable :: eigenvectors(:, :)
        real(dp), allocatable :: mu(:)
        real(dp) :: cutoff
        real(dp) :: d
        real(dp) :: half_df
        real(dp) :: half_sum
        real(dp) :: log_pdet
        real(dp) :: nu
        real(dp) :: threshold

        if (present(loc)) then
            n = size(loc)
        else if (present(shape)) then
            n = size(shape, 1)
        else
            n = 1
        end if
        if (n < 1) then
            h = quiet_nan(0.0_dp)
            return
        end if
        allocate(mu(n), matrix(n, n), eigenvalues(n), eigenvectors(n, n))
        call process_loc_shape(n, loc, shape, df, mu, matrix, nu, info)
        if (info /= multivariate_t_status_success) then
            h = quiet_nan(0.0_dp)
            return
        end if
        if (.not. ieee_is_finite(nu)) then
            h = multivariate_normal_entropy(mu, matrix)
            return
        end if

        call psd_decompose(matrix, eigenvalues, eigenvectors, rank, log_pdet, cutoff, info, .true.)
        if (info /= linalg_status_success) then
            h = quiet_nan(0.0_dp)
            return
        end if

        d = real(n, dp)
        threshold = d * 400.0_dp / (log(d) + 1.0_dp)
        if (nu >= threshold) then
            h = d * (scifort_log_sqrt_two_pi + 0.5_dp) + d / nu &
                - d * (d - 2.0_dp) / (4.0_dp * nu ** 2) &
                + d * d * (d - 2.0_dp) / (6.0_dp * nu ** 3) &
                + d * (-3.0_dp * d ** 3 + 8.0_dp * d ** 2 - 8.0_dp) / (24.0_dp * nu ** 4) &
                + d * d * (3.0_dp * d ** 3 - 10.0_dp * d ** 2 + 16.0_dp) / (30.0_dp * nu ** 5) &
                + 0.5_dp * log_pdet
        else
            half_df = 0.5_dp * nu
            half_sum = 0.5_dp * (real(n, dp) + nu)
            h = -log_gamma(half_sum) + log_gamma(half_df) &
                + 0.5_dp * real(n, dp) * log(nu * scifort_pi) &
                + half_sum * (digamma_positive(half_sum) - digamma_positive(half_df)) &
                + 0.5_dp * log_pdet
        end if
    end function multivariate_t_entropy

    subroutine multivariate_t_rvs(state, sample, loc, shape, df, status)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by the draw
        real(dp), intent(out) :: sample(:) !! generated vector; length defines dimension
        real(dp), intent(in), optional :: loc(:) !! location vector; zero when absent
        real(dp), intent(in), optional :: shape(:, :) !! positive-semidefinite shape matrix; identity when absent
        real(dp), intent(in), optional :: df !! degrees of freedom, > 0; default 1; +infinity gives normal
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid input

        integer :: info
        integer :: n
        real(dp), allocatable :: matrix(:, :)
        real(dp), allocatable :: mu(:)
        real(dp), allocatable :: zero(:)
        real(dp), allocatable :: z(:)
        real(dp) :: nu
        real(dp) :: scale_factor

        n = size(sample)
        if (n < 1) then
            if (present(status)) status = multivariate_t_status_invalid_shape
            return
        end if
        allocate(mu(n), matrix(n, n), zero(n), z(n))
        call process_loc_shape(n, loc, shape, df, mu, matrix, nu, info)
        if (info /= multivariate_t_status_success) then
            sample = quiet_nan(0.0_dp)
            if (present(status)) status = info
            return
        end if
        zero = 0.0_dp
        call multivariate_normal_rvs(state, z, zero, matrix, info)
        if (info /= 0) then
            sample = quiet_nan(0.0_dp)
            if (present(status)) status = multivariate_t_status_linalg_failure
            return
        end if

        if (.not. ieee_is_finite(nu)) then
            scale_factor = 1.0_dp
        else
            scale_factor = chi2_ppf(rng_uniform(state), nu) / nu
            if (.not. ieee_is_finite(scale_factor) .or. scale_factor <= 0.0_dp) then
                sample = quiet_nan(0.0_dp)
                if (present(status)) status = multivariate_t_status_invalid_parameter
                return
            end if
        end if
        sample = mu + z / sqrt(scale_factor)
        if (present(status)) status = multivariate_t_status_success
    end subroutine multivariate_t_rvs

    subroutine multivariate_t_rvs_array(state, samples, loc, shape, df, status)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by all draws
        real(dp), intent(out) :: samples(:, :) !! samples by column, shape (dimension, n_samples)
        real(dp), intent(in), optional :: loc(:) !! location vector; zero when absent
        real(dp), intent(in), optional :: shape(:, :) !! positive-semidefinite shape matrix; identity when absent
        real(dp), intent(in), optional :: df !! degrees of freedom, > 0; default 1; +infinity gives normal
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid input

        integer :: info
        integer :: j

        do j = 1, size(samples, 2)
            call multivariate_t_rvs(state, samples(:, j), loc, shape, df, info)
            if (info /= multivariate_t_status_success) then
                samples = quiet_nan(0.0_dp)
                if (present(status)) status = info
                return
            end if
        end do
        if (present(status)) status = multivariate_t_status_success
    end subroutine multivariate_t_rvs_array

    subroutine multivariate_t_marginal(dimensions, loc, shape, marginal_loc, marginal_shape, status)
        integer, intent(in) :: dimensions(:) !! one-based unique indices of retained variables
        real(dp), intent(in) :: loc(:) !! original location vector
        real(dp), intent(in) :: shape(:, :) !! original shape matrix; lower triangle is authoritative
        real(dp), intent(out) :: marginal_loc(:) !! selected location vector
        real(dp), intent(out) :: marginal_shape(:, :) !! selected shape matrix
        integer, intent(out) :: status !! zero on success; nonzero for invalid dimensions/shapes

        integer :: i
        integer :: j
        integer :: m
        integer :: n

        n = size(loc)
        m = size(dimensions)
        if (n < 1 .or. m < 1 .or. size(shape, 1) /= n .or. size(shape, 2) /= n .or. &
                size(marginal_loc) /= m .or. size(marginal_shape, 1) /= m .or. &
                size(marginal_shape, 2) /= m) then
            status = multivariate_t_status_invalid_shape
            return
        end if
        if (any(.not. ieee_is_finite(loc)) .or. any(.not. ieee_is_finite(shape)) .or. &
                any(dimensions < 1) .or. any(dimensions > n)) then
            status = multivariate_t_status_invalid_parameter
            return
        end if
        do i = 1, m
            do j = 1, i - 1
                if (dimensions(i) == dimensions(j)) then
                    status = multivariate_t_status_invalid_parameter
                    return
                end if
            end do
        end do

        do i = 1, m
            marginal_loc(i) = loc(dimensions(i))
            do j = 1, m
                if (dimensions(i) >= dimensions(j)) then
                    marginal_shape(i, j) = shape(dimensions(i), dimensions(j))
                else
                    marginal_shape(i, j) = shape(dimensions(j), dimensions(i))
                end if
            end do
        end do
        status = multivariate_t_status_success
    end subroutine multivariate_t_marginal

    function logpdf_core(x, mu, eigenvalues, eigenvectors, log_pdet, cutoff, df) result(y)
        real(dp), intent(in) :: x(:) !! evaluation point
        real(dp), intent(in) :: mu(:) !! location vector
        real(dp), intent(in) :: eigenvalues(:) !! ascending shape eigenvalues
        real(dp), intent(in) :: eigenvectors(:, :) !! corresponding eigenvectors
        real(dp), intent(in) :: log_pdet !! log pseudo-determinant of shape
        real(dp), intent(in) :: cutoff !! numerical rank cutoff
        real(dp), intent(in) :: df !! finite positive degrees of freedom
        real(dp) :: y

        integer :: i
        integer :: n
        real(dp) :: maha
        real(dp) :: projection
        real(dp) :: t
        real(dp), allocatable :: dev(:)

        n = size(x)
        allocate(dev(n))
        dev = x - mu
        maha = 0.0_dp
        do i = 1, n
            if (eigenvalues(i) > cutoff) then
                projection = dot_product(dev, eigenvectors(:, i))
                maha = maha + projection * projection / eigenvalues(i)
            end if
        end do

        t = 0.5_dp * (df + real(n, dp))
        y = log_gamma(t) - log_gamma(0.5_dp * df) &
            - 0.5_dp * real(n, dp) * log(df * scifort_pi) &
            - 0.5_dp * log_pdet - t * log(1.0_dp + maha / df)
    end function logpdf_core

    subroutine process_loc_shape(n, loc, shape, df, mu, matrix, nu, status)
        integer, intent(in) :: n !! requested dimension
        real(dp), intent(in), optional :: loc(:) !! optional location vector
        real(dp), intent(in), optional :: shape(:, :) !! optional shape matrix
        real(dp), intent(in), optional :: df !! optional degrees of freedom
        real(dp), intent(out) :: mu(:) !! expanded location vector
        real(dp), intent(out) :: matrix(:, :) !! expanded symmetric shape matrix
        real(dp), intent(out) :: nu !! processed degrees of freedom
        integer, intent(out) :: status !! zero on success

        integer :: i

        if (size(mu) /= n .or. size(matrix, 1) /= n .or. size(matrix, 2) /= n) then
            status = multivariate_t_status_invalid_shape
            return
        end if
        if (present(loc)) then
            if (size(loc) /= n .or. any(.not. ieee_is_finite(loc))) then
                status = multivariate_t_status_invalid_parameter
                return
            end if
            mu = loc
        else
            mu = 0.0_dp
        end if

        if (present(shape)) then
            if (size(shape, 1) /= n .or. size(shape, 2) /= n) then
                status = multivariate_t_status_invalid_shape
                return
            end if
            if (any(.not. ieee_is_finite(shape))) then
                status = multivariate_t_status_invalid_parameter
                return
            end if
            call copy_lower_symmetric(shape, matrix)
        else
            matrix = 0.0_dp
            do i = 1, n
                matrix(i, i) = 1.0_dp
            end do
        end if

        nu = 1.0_dp
        if (present(df)) nu = df
        if (ieee_is_nan(nu) .or. nu <= 0.0_dp) then
            status = multivariate_t_status_invalid_parameter
            return
        end if
        status = multivariate_t_status_success
    end subroutine process_loc_shape

    pure subroutine copy_lower_symmetric(source, target)
        real(dp), intent(in) :: source(:, :) !! source whose lower triangle is authoritative
        real(dp), intent(out) :: target(:, :) !! symmetric lower-triangle copy

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

end module scifort_multivariate_t
