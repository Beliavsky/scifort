! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Gaussian kernel density estimation compatible with the core numerical
! behavior of scipy.stats.gaussian_kde. The estimator stores observations by
! column, uses unbiased weighted sample covariance, and supports Scott,
! Silverman, and constant bandwidth factors.

module scifort_gaussian_kde
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_linalg, only : cholesky_lower, linalg_status_success
    use scifort_matrix_distribution_helpers, only : matrix_status_success, &
        solve_lower_matrix, solve_spd_matrix
    use scifort_multivariate_normal, only : multivariate_normal_cdf
    use scifort_math, only : quiet_nan
    use scifort_normal, only : normal_cdf, normal_ppf
    use scifort_random, only : rng_state, rng_uniform
    implicit none
    private

    integer, parameter, public :: gaussian_kde_status_success = 0
    integer, parameter, public :: gaussian_kde_status_invalid_shape = 1
    integer, parameter, public :: gaussian_kde_status_invalid_parameter = 2
    integer, parameter, public :: gaussian_kde_status_singular = 3
    integer, parameter, public :: gaussian_kde_status_linalg_failure = 4

    integer, parameter :: bandwidth_scott = 1
    integer, parameter :: bandwidth_silverman = 2
    integer, parameter :: bandwidth_constant = 3

    type, public :: gaussian_kde
        integer :: d = 0 !! number of data dimensions
        integer :: n = 0 !! number of observations
        real(dp) :: neff = 0.0_dp !! effective sample size, 1/sum(weights**2)
        real(dp) :: factor = 0.0_dp !! current scalar bandwidth factor
        real(dp), allocatable :: dataset(:, :) !! observations by column, shape (d,n)
        real(dp), allocatable :: weights(:) !! normalized nonnegative observation weights
        real(dp), allocatable :: covariance(:, :) !! Gaussian-kernel covariance matrix
        real(dp), allocatable :: cho_cov(:, :) !! lower Cholesky factor of covariance
        real(dp), allocatable :: data_covariance(:, :) !! unbiased weighted sample covariance
        real(dp), allocatable :: data_cho_cov(:, :) !! lower Cholesky factor of data covariance
        integer, private :: bandwidth_mode = bandwidth_scott
        real(dp), private :: bandwidth_constant_value = 0.0_dp
    contains
        procedure, private :: evaluate_point => gaussian_kde_evaluate_point
        procedure, private :: evaluate_points => gaussian_kde_evaluate_points
        generic, public :: evaluate => evaluate_point, evaluate_points
        generic, public :: pdf => evaluate_point, evaluate_points
        procedure, private :: logpdf_point => gaussian_kde_logpdf_point
        procedure, private :: logpdf_points => gaussian_kde_logpdf_points
        generic, public :: logpdf => logpdf_point, logpdf_points
        procedure, public :: scotts_factor => gaussian_kde_scotts_factor
        procedure, public :: silverman_factor => gaussian_kde_silverman_factor
        procedure, private :: set_bandwidth_name => gaussian_kde_set_bandwidth_name
        procedure, private :: set_bandwidth_factor => gaussian_kde_set_bandwidth_factor
        generic, public :: set_bandwidth => set_bandwidth_name, set_bandwidth_factor
        procedure, public :: integrate_gaussian => gaussian_kde_integrate_gaussian
        procedure, public :: integrate_box_1d => gaussian_kde_integrate_box_1d
        procedure, public :: integrate_box => gaussian_kde_integrate_box
        procedure, public :: integrate_kde => gaussian_kde_integrate_kde
        procedure, public :: resample => gaussian_kde_resample
        procedure, public :: marginal => gaussian_kde_marginal
        procedure, public :: inv_cov => gaussian_kde_inv_cov
    end type gaussian_kde

    interface gaussian_kde_init
        module procedure gaussian_kde_init_1d
        module procedure gaussian_kde_init_nd
    end interface gaussian_kde_init

    public :: gaussian_kde_init

contains

    subroutine gaussian_kde_init_1d(kde, dataset, status, weights, bw_method, bw_factor)
        type(gaussian_kde), intent(out) :: kde !! initialized KDE object
        real(dp), intent(in) :: dataset(:) !! univariate observations
        integer, intent(out) :: status !! gaussian_kde_status_* result
        real(dp), intent(in), optional :: weights(:) !! optional nonnegative observation weights
        character(len=*), intent(in), optional :: bw_method !! 'scott' or 'silverman'
        real(dp), intent(in), optional :: bw_factor !! positive constant bandwidth factor

        real(dp), allocatable :: data2(:, :)

        allocate(data2(1, size(dataset)))
        data2(1, :) = dataset
        call gaussian_kde_init_nd(kde, data2, status, weights, bw_method, bw_factor)
    end subroutine gaussian_kde_init_1d

    subroutine gaussian_kde_init_nd(kde, dataset, status, weights, bw_method, bw_factor)
        type(gaussian_kde), intent(out) :: kde !! initialized KDE object
        real(dp), intent(in) :: dataset(:, :) !! observations by column, shape (dimension,n)
        integer, intent(out) :: status !! gaussian_kde_status_* result
        real(dp), intent(in), optional :: weights(:) !! optional nonnegative observation weights
        character(len=*), intent(in), optional :: bw_method !! 'scott' or 'silverman'
        real(dp), intent(in), optional :: bw_factor !! positive constant bandwidth factor

        integer :: i
        integer :: info
        real(dp), allocatable :: dev(:)
        real(dp), allocatable :: mean(:)
        real(dp) :: denom
        real(dp) :: weight_sum

        status = gaussian_kde_status_success
        kde%d = size(dataset, 1)
        kde%n = size(dataset, 2)
        if (kde%d < 1 .or. kde%n < 2 .or. kde%d > kde%n .or. .not. all(ieee_is_finite(dataset))) then
            status = gaussian_kde_status_invalid_shape
            return
        end if
        if (present(bw_method) .and. present(bw_factor)) then
            status = gaussian_kde_status_invalid_parameter
            return
        end if

        allocate(kde%dataset(kde%d,kde%n), kde%weights(kde%n))
        allocate(kde%data_covariance(kde%d,kde%d), kde%data_cho_cov(kde%d,kde%d))
        allocate(kde%covariance(kde%d,kde%d), kde%cho_cov(kde%d,kde%d))
        kde%dataset = dataset

        if (present(weights)) then
            if (size(weights) /= kde%n .or. .not. all(ieee_is_finite(weights)) .or. &
                    any(weights < 0.0_dp)) then
                status = gaussian_kde_status_invalid_parameter
                return
            end if
            weight_sum = sum(weights)
            if (.not. ieee_is_finite(weight_sum) .or. weight_sum <= 0.0_dp) then
                status = gaussian_kde_status_invalid_parameter
                return
            end if
            kde%weights = weights / weight_sum
        else
            kde%weights = 1.0_dp / real(kde%n, dp)
        end if

        kde%neff = 1.0_dp / sum(kde%weights * kde%weights)
        denom = 1.0_dp - sum(kde%weights * kde%weights)
        if (denom <= 0.0_dp .or. .not. ieee_is_finite(denom)) then
            status = gaussian_kde_status_singular
            return
        end if

        allocate(mean(kde%d), dev(kde%d))
        mean = matmul(kde%dataset, kde%weights)
        kde%data_covariance = 0.0_dp
        do i = 1, kde%n
            dev = kde%dataset(:,i) - mean
            kde%data_covariance = kde%data_covariance + kde%weights(i) * outer_product(dev, dev)
        end do
        kde%data_covariance = kde%data_covariance / denom
        call cholesky_lower(kde%data_covariance, kde%data_cho_cov, info)
        if (info /= linalg_status_success) then
            status = gaussian_kde_status_singular
            return
        end if

        kde%bandwidth_mode = bandwidth_scott
        if (present(bw_method)) then
            select case (trim(lower_ascii(bw_method)))
            case ('scott')
                kde%bandwidth_mode = bandwidth_scott
            case ('silverman')
                kde%bandwidth_mode = bandwidth_silverman
            case default
                status = gaussian_kde_status_invalid_parameter
                return
            end select
        else if (present(bw_factor)) then
            if (.not. ieee_is_finite(bw_factor) .or. bw_factor <= 0.0_dp) then
                status = gaussian_kde_status_invalid_parameter
                return
            end if
            kde%bandwidth_mode = bandwidth_constant
            kde%bandwidth_constant_value = bw_factor
        end if
        call compute_covariance(kde, status)
    end subroutine gaussian_kde_init_nd

    function gaussian_kde_scotts_factor(self) result(value)
        class(gaussian_kde), intent(in) :: self !! initialized KDE object
        real(dp) :: value

        value = self%neff ** (-1.0_dp / real(self%d + 4, dp))
    end function gaussian_kde_scotts_factor

    function gaussian_kde_silverman_factor(self) result(value)
        class(gaussian_kde), intent(in) :: self !! initialized KDE object
        real(dp) :: value

        value = (self%neff * real(self%d + 2, dp) / 4.0_dp) ** &
            (-1.0_dp / real(self%d + 4, dp))
    end function gaussian_kde_silverman_factor

    subroutine gaussian_kde_set_bandwidth_name(self, bw_method, status)
        class(gaussian_kde), intent(inout) :: self !! KDE object whose bandwidth is changed
        character(len=*), intent(in) :: bw_method !! 'scott' or 'silverman'
        integer, intent(out), optional :: status !! gaussian_kde_status_* result

        integer :: info

        select case (trim(lower_ascii(bw_method)))
        case ('scott')
            self%bandwidth_mode = bandwidth_scott
        case ('silverman')
            self%bandwidth_mode = bandwidth_silverman
        case default
            if (present(status)) status = gaussian_kde_status_invalid_parameter
            return
        end select
        call compute_covariance(self, info)
        if (present(status)) status = info
    end subroutine gaussian_kde_set_bandwidth_name

    subroutine gaussian_kde_set_bandwidth_factor(self, bw_factor, status)
        class(gaussian_kde), intent(inout) :: self !! KDE object whose bandwidth is changed
        real(dp), intent(in) :: bw_factor !! positive constant bandwidth factor
        integer, intent(out), optional :: status !! gaussian_kde_status_* result

        integer :: info

        if (.not. ieee_is_finite(bw_factor) .or. bw_factor <= 0.0_dp) then
            if (present(status)) status = gaussian_kde_status_invalid_parameter
            return
        end if
        self%bandwidth_mode = bandwidth_constant
        self%bandwidth_constant_value = bw_factor
        call compute_covariance(self, info)
        if (present(status)) status = info
    end subroutine gaussian_kde_set_bandwidth_factor

    function gaussian_kde_evaluate_point(self, point, status) result(value)
        class(gaussian_kde), intent(in) :: self !! initialized KDE object
        real(dp), intent(in) :: point(:) !! one evaluation point, length d
        integer, intent(out), optional :: status !! gaussian_kde_status_* result
        real(dp) :: value

        real(dp), allocatable :: points(:, :)
        real(dp), allocatable :: values(:)
        integer :: info

        allocate(points(size(point),1))
        points(:,1) = point
        values = evaluate_many(self, points, .false., info)
        if (size(values) == 1) then
            value = values(1)
        else
            value = quiet_nan(0.0_dp)
        end if
        if (present(status)) status = info
    end function gaussian_kde_evaluate_point

    function gaussian_kde_evaluate_points(self, points, status) result(values)
        class(gaussian_kde), intent(in) :: self !! initialized KDE object
        real(dp), intent(in) :: points(:, :) !! evaluation points by column, shape (d,m)
        integer, intent(out), optional :: status !! gaussian_kde_status_* result
        real(dp), allocatable :: values(:)

        integer :: info

        values = evaluate_many(self, points, .false., info)
        if (present(status)) status = info
    end function gaussian_kde_evaluate_points

    function gaussian_kde_logpdf_point(self, point, status) result(value)
        class(gaussian_kde), intent(in) :: self !! initialized KDE object
        real(dp), intent(in) :: point(:) !! one evaluation point, length d
        integer, intent(out), optional :: status !! gaussian_kde_status_* result
        real(dp) :: value

        real(dp), allocatable :: points(:, :)
        real(dp), allocatable :: values(:)
        integer :: info

        allocate(points(size(point),1))
        points(:,1) = point
        values = evaluate_many(self, points, .true., info)
        if (size(values) == 1) then
            value = values(1)
        else
            value = quiet_nan(0.0_dp)
        end if
        if (present(status)) status = info
    end function gaussian_kde_logpdf_point

    function gaussian_kde_logpdf_points(self, points, status) result(values)
        class(gaussian_kde), intent(in) :: self !! initialized KDE object
        real(dp), intent(in) :: points(:, :) !! evaluation points by column, shape (d,m)
        integer, intent(out), optional :: status !! gaussian_kde_status_* result
        real(dp), allocatable :: values(:)

        integer :: info

        values = evaluate_many(self, points, .true., info)
        if (present(status)) status = info
    end function gaussian_kde_logpdf_points

    function gaussian_kde_integrate_gaussian(self, mean, cov, status) result(value)
        class(gaussian_kde), intent(in) :: self !! initialized KDE object
        real(dp), intent(in) :: mean(:) !! Gaussian mean vector, length d
        real(dp), intent(in) :: cov(:, :) !! Gaussian covariance matrix, shape (d,d)
        integer, intent(out), optional :: status !! gaussian_kde_status_* result
        real(dp) :: value

        integer :: i
        integer :: info
        real(dp), allocatable :: diff(:, :)
        real(dp), allocatable :: l(:, :)
        real(dp), allocatable :: solved(:, :)
        real(dp), allocatable :: sum_cov(:, :)
        real(dp) :: energy
        real(dp) :: norm_const

        value = 0.0_dp
        if (size(mean) /= self%d .or. size(cov,1) /= self%d .or. size(cov,2) /= self%d .or. &
                .not. all(ieee_is_finite(mean)) .or. .not. all(ieee_is_finite(cov))) then
            if (present(status)) status = gaussian_kde_status_invalid_shape
            return
        end if
        allocate(sum_cov(self%d,self%d), l(self%d,self%d))
        allocate(diff(self%d,self%n), solved(self%d,self%n))
        sum_cov = self%covariance + cov
        call cholesky_lower(sum_cov, l, info)
        if (info /= linalg_status_success) then
            if (present(status)) status = gaussian_kde_status_linalg_failure
            return
        end if
        do i = 1, self%n
            diff(:,i) = self%dataset(:,i) - mean
        end do
        call solve_lower_matrix(l, diff, solved, info)
        if (info /= matrix_status_success) then
            if (present(status)) status = gaussian_kde_status_linalg_failure
            return
        end if
        norm_const = exp(gaussian_log_norm_from_cholesky(l))
        do i = 1, self%n
            energy = 0.5_dp * sum(solved(:,i) * solved(:,i))
            value = value + self%weights(i) * exp(-energy)
        end do
        value = value / norm_const
        if (present(status)) status = gaussian_kde_status_success
    end function gaussian_kde_integrate_gaussian

    function gaussian_kde_integrate_box_1d(self, low, high, status) result(value)
        class(gaussian_kde), intent(in) :: self !! initialized one-dimensional KDE object
        real(dp), intent(in) :: low !! lower integration bound
        real(dp), intent(in) :: high !! upper integration bound
        integer, intent(out), optional :: status !! gaussian_kde_status_* result
        real(dp) :: value

        integer :: i
        real(dp) :: stdev

        value = 0.0_dp
        if (self%d /= 1) then
            if (present(status)) status = gaussian_kde_status_invalid_shape
            return
        end if
        stdev = sqrt(self%covariance(1,1))
        do i = 1, self%n
            value = value + self%weights(i) * &
                (normal_cdf(high, self%dataset(1,i), stdev) - normal_cdf(low, self%dataset(1,i), stdev))
        end do
        if (present(status)) status = gaussian_kde_status_success
    end function gaussian_kde_integrate_box_1d

    function gaussian_kde_integrate_box(self, low, high, maxpts, state, status) result(value)
        class(gaussian_kde), intent(in) :: self !! initialized KDE object
        real(dp), intent(in) :: low(:) !! lower rectangular integration bounds
        real(dp), intent(in) :: high(:) !! upper rectangular integration bounds
        integer, intent(in), optional :: maxpts !! approximate maximum QMC points per Gaussian component
        type(rng_state), intent(inout), optional :: state !! explicit randomized-QMC state
        integer, intent(out), optional :: status !! gaussian_kde_status_* result
        real(dp) :: value

        integer :: i
        real(dp), allocatable :: hi_dev(:)
        real(dp), allocatable :: lo_dev(:)
        real(dp) :: p

        value = 0.0_dp
        if (size(low) /= self%d .or. size(high) /= self%d) then
            if (present(status)) status = gaussian_kde_status_invalid_shape
            return
        end if
        allocate(hi_dev(self%d), lo_dev(self%d))
        do i = 1, self%n
            hi_dev = high - self%dataset(:,i)
            lo_dev = low - self%dataset(:,i)
            if (present(state)) then
                p = multivariate_normal_cdf(hi_dev, cov=self%covariance, lower_limit=lo_dev, &
                    maxpts=maxpts, state=state)
            else
                p = multivariate_normal_cdf(hi_dev, cov=self%covariance, lower_limit=lo_dev, &
                    maxpts=maxpts)
            end if
            if (.not. ieee_is_finite(p)) then
                if (present(status)) status = gaussian_kde_status_linalg_failure
                return
            end if
            value = value + self%weights(i) * p
        end do
        if (present(status)) status = gaussian_kde_status_success
    end function gaussian_kde_integrate_box

    function gaussian_kde_integrate_kde(self, other, status) result(value)
        class(gaussian_kde), intent(in) :: self !! first KDE object
        class(gaussian_kde), intent(in) :: other !! second KDE object of equal dimension
        integer, intent(out), optional :: status !! gaussian_kde_status_* result
        real(dp) :: value

        integer :: i
        integer :: info
        integer :: j
        real(dp), allocatable :: diff(:, :)
        real(dp), allocatable :: l(:, :)
        real(dp), allocatable :: solved(:, :)
        real(dp), allocatable :: sum_cov(:, :)
        real(dp) :: norm_const

        value = 0.0_dp
        if (self%d /= other%d .or. self%d < 1) then
            if (present(status)) status = gaussian_kde_status_invalid_shape
            return
        end if
        allocate(sum_cov(self%d,self%d), l(self%d,self%d))
        allocate(diff(self%d,other%n), solved(self%d,other%n))
        sum_cov = self%covariance + other%covariance
        call cholesky_lower(sum_cov, l, info)
        if (info /= linalg_status_success) then
            if (present(status)) status = gaussian_kde_status_linalg_failure
            return
        end if
        norm_const = exp(gaussian_log_norm_from_cholesky(l))
        do i = 1, self%n
            do j = 1, other%n
                diff(:,j) = other%dataset(:,j) - self%dataset(:,i)
            end do
            call solve_lower_matrix(l, diff, solved, info)
            if (info /= matrix_status_success) then
                if (present(status)) status = gaussian_kde_status_linalg_failure
                return
            end if
            do j = 1, other%n
                value = value + self%weights(i) * other%weights(j) * &
                    exp(-0.5_dp * sum(solved(:,j) * solved(:,j)))
            end do
        end do
        value = value / norm_const
        if (present(status)) status = gaussian_kde_status_success
    end function gaussian_kde_integrate_kde

    subroutine gaussian_kde_resample(self, state, samples, status)
        class(gaussian_kde), intent(in) :: self !! initialized KDE object
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by all draws
        real(dp), intent(out) :: samples(:, :) !! generated samples by column, shape (d,m)
        integer, intent(out), optional :: status !! gaussian_kde_status_* result

        integer :: i
        integer :: index
        integer :: j
        real(dp), allocatable :: z(:)
        real(dp) :: cumulative
        real(dp) :: u

        if (size(samples,1) /= self%d) then
            if (present(status)) status = gaussian_kde_status_invalid_shape
            return
        end if
        allocate(z(self%d))
        do j = 1, size(samples,2)
            u = rng_uniform(state)
            cumulative = 0.0_dp
            index = self%n
            do i = 1, self%n
                cumulative = cumulative + self%weights(i)
                if (u <= cumulative) then
                    index = i
                    exit
                end if
            end do
            do i = 1, self%d
                z(i) = normal_ppf(rng_uniform(state))
            end do
            samples(:,j) = self%dataset(:,index) + matmul(self%cho_cov, z)
        end do
        if (present(status)) status = gaussian_kde_status_success
    end subroutine gaussian_kde_resample

    function gaussian_kde_marginal(self, dimensions, status) result(marginal)
        class(gaussian_kde), intent(in) :: self !! initialized multivariate KDE object
        integer, intent(in) :: dimensions(:) !! unique one-based dimensions to retain
        integer, intent(out), optional :: status !! gaussian_kde_status_* result
        type(gaussian_kde) :: marginal

        integer :: i
        integer :: info
        integer :: j
        real(dp), allocatable :: selected(:, :)

        if (size(dimensions) < 1 .or. any(dimensions < 1) .or. any(dimensions > self%d)) then
            if (present(status)) status = gaussian_kde_status_invalid_parameter
            return
        end if
        do i = 1, size(dimensions)
            do j = 1, i - 1
                if (dimensions(i) == dimensions(j)) then
                    if (present(status)) status = gaussian_kde_status_invalid_parameter
                    return
                end if
            end do
        end do
        allocate(selected(size(dimensions),self%n))
        do i = 1, size(dimensions)
            selected(i,:) = self%dataset(dimensions(i),:)
        end do
        call gaussian_kde_init_nd(marginal, selected, info, weights=self%weights, bw_factor=self%factor)
        if (present(status)) status = info
    end function gaussian_kde_marginal

    subroutine gaussian_kde_inv_cov(self, inv_cov, status)
        class(gaussian_kde), intent(in) :: self !! initialized KDE object
        real(dp), intent(out) :: inv_cov(:, :) !! inverse of current kernel covariance
        integer, intent(out), optional :: status !! gaussian_kde_status_* result

        integer :: i
        integer :: info
        real(dp), allocatable :: ident(:, :)

        if (size(inv_cov,1) /= self%d .or. size(inv_cov,2) /= self%d) then
            if (present(status)) status = gaussian_kde_status_invalid_shape
            return
        end if
        allocate(ident(self%d,self%d))
        ident = 0.0_dp
        do i = 1, self%d
            ident(i,i) = 1.0_dp
        end do
        call solve_spd_matrix(self%cho_cov, ident, inv_cov, info)
        if (info /= matrix_status_success) then
            if (present(status)) status = gaussian_kde_status_linalg_failure
        else if (present(status)) then
            status = gaussian_kde_status_success
        end if
    end subroutine gaussian_kde_inv_cov

    subroutine compute_covariance(self, status)
        class(gaussian_kde), intent(inout) :: self !! KDE object with cached data covariance
        integer, intent(out) :: status !! gaussian_kde_status_* result

        select case (self%bandwidth_mode)
        case (bandwidth_scott)
            self%factor = self%scotts_factor()
        case (bandwidth_silverman)
            self%factor = self%silverman_factor()
        case (bandwidth_constant)
            self%factor = self%bandwidth_constant_value
        case default
            status = gaussian_kde_status_invalid_parameter
            return
        end select
        if (.not. ieee_is_finite(self%factor) .or. self%factor <= 0.0_dp) then
            status = gaussian_kde_status_invalid_parameter
            return
        end if
        self%covariance = self%data_covariance * self%factor * self%factor
        self%cho_cov = self%data_cho_cov * self%factor
        status = gaussian_kde_status_success
    end subroutine compute_covariance

    function evaluate_many(self, points, return_log, status) result(values)
        class(gaussian_kde), intent(in) :: self !! initialized KDE object
        real(dp), intent(in) :: points(:, :) !! points by column, shape (d,m)
        logical, intent(in) :: return_log !! return log density when true
        integer, intent(out) :: status !! gaussian_kde_status_* result
        real(dp), allocatable :: values(:)

        integer :: i
        integer :: info
        integer :: j
        real(dp), allocatable :: diff(:, :)
        real(dp), allocatable :: solved(:, :)
        real(dp), allocatable :: terms(:)
        real(dp) :: log_norm
        real(dp) :: max_term
        real(dp) :: sum_exp

        allocate(values(size(points,2)))
        values = 0.0_dp
        if (size(points,1) /= self%d .or. .not. all(ieee_is_finite(points))) then
            status = gaussian_kde_status_invalid_shape
            return
        end if
        allocate(diff(self%d,self%n), solved(self%d,self%n), terms(self%n))
        log_norm = gaussian_log_norm_from_cholesky(self%cho_cov)
        do j = 1, size(points,2)
            do i = 1, self%n
                diff(:,i) = points(:,j) - self%dataset(:,i)
            end do
            call solve_lower_matrix(self%cho_cov, diff, solved, info)
            if (info /= matrix_status_success) then
                status = gaussian_kde_status_linalg_failure
                return
            end if
            max_term = -huge(1.0_dp)
            do i = 1, self%n
                if (self%weights(i) > 0.0_dp) then
                    terms(i) = log(self%weights(i)) - 0.5_dp * sum(solved(:,i) * solved(:,i))
                    max_term = max(max_term, terms(i))
                else
                    terms(i) = -huge(1.0_dp)
                end if
            end do
            sum_exp = 0.0_dp
            do i = 1, self%n
                if (self%weights(i) > 0.0_dp) sum_exp = sum_exp + exp(terms(i) - max_term)
            end do
            if (return_log) then
                values(j) = max_term + log(sum_exp) - log_norm
            else
                values(j) = exp(max_term + log(sum_exp) - log_norm)
            end if
        end do
        status = gaussian_kde_status_success
    end function evaluate_many

    pure function gaussian_log_norm_from_cholesky(l) result(value)
        real(dp), intent(in) :: l(:, :) !! lower Cholesky factor of covariance
        real(dp) :: value

        integer :: i

        value = 0.5_dp * real(size(l,1),dp) * log(2.0_dp * acos(-1.0_dp))
        do i = 1, size(l,1)
            value = value + log(l(i,i))
        end do
    end function gaussian_log_norm_from_cholesky

    pure function outer_product(x, y) result(a)
        real(dp), intent(in) :: x(:) !! left vector
        real(dp), intent(in) :: y(:) !! right vector
        real(dp) :: a(size(x),size(y))

        integer :: i

        do i = 1, size(x)
            a(i,:) = x(i) * y
        end do
    end function outer_product

    pure function lower_ascii(text) result(lower)
        character(len=*), intent(in) :: text !! ASCII text to convert to lowercase
        character(len=len(text)) :: lower

        integer :: code
        integer :: i

        lower = text
        do i = 1, len(text)
            code = iachar(text(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
        end do
    end function lower_ascii

end module scifort_gaussian_kde
