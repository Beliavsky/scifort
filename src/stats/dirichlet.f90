! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_dirichlet
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_gamma, only : gamma_ppf
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    use scifort_random, only : rng_state, rng_uniform
    use scifort_special_elementary, only : digamma, gammaln, xlogy
    implicit none
    private

    integer, parameter, public :: dirichlet_status_success = 0
    integer, parameter, public :: dirichlet_status_invalid_parameter = 1
    integer, parameter, public :: dirichlet_status_invalid_shape = 2

    public :: dirichlet_cov
    public :: dirichlet_entropy
    public :: dirichlet_logpdf
    public :: dirichlet_mean
    public :: dirichlet_pdf
    public :: dirichlet_rvs
    public :: dirichlet_rvs_array
    public :: dirichlet_var

contains

    function dirichlet_logpdf(x, alpha) result(y)
        real(dp), intent(in) :: x(:) !! simplex point with K or K-1 entries
        real(dp), intent(in) :: alpha(:) !! positive concentration parameters
        real(dp) :: y

        integer :: k
        real(dp), allocatable :: xx(:)
        real(dp) :: alpha0
        real(dp) :: log_beta

        k = size(alpha)
        if (.not. valid_alpha(alpha) .or. (size(x) /= k .and. size(x) /= k - 1)) then
            y = quiet_nan(0.0_dp)
            return
        end if

        allocate(xx(k))
        if (size(x) == k) then
            xx = x
        else
            if (k <= 1) then
                y = quiet_nan(0.0_dp)
                return
            end if
            xx(1:k - 1) = x
            xx(k) = 1.0_dp - sum(x)
        end if

        if (any(.not. ieee_is_finite(xx)) .or. any(xx < 0.0_dp) .or. any(xx > 1.0_dp)) then
            y = quiet_nan(0.0_dp)
            return
        end if
        if (any((xx == 0.0_dp) .and. (alpha < 1.0_dp))) then
            y = quiet_nan(0.0_dp)
            return
        end if
        if (abs(sum(xx) - 1.0_dp) > 1.0e-9_dp) then
            y = quiet_nan(0.0_dp)
            return
        end if

        alpha0 = sum(alpha)
        log_beta = sum(gammaln(alpha)) - gammaln(alpha0)
        y = -log_beta + sum(xlogy(alpha - 1.0_dp, xx))
    end function dirichlet_logpdf

    function dirichlet_pdf(x, alpha) result(y)
        real(dp), intent(in) :: x(:) !! simplex point with K or K-1 entries
        real(dp), intent(in) :: alpha(:) !! positive concentration parameters
        real(dp) :: y

        y = exp(dirichlet_logpdf(x, alpha))
    end function dirichlet_pdf

    function dirichlet_mean(alpha) result(mu)
        real(dp), intent(in) :: alpha(:) !! positive concentration parameters
        real(dp) :: mu(size(alpha))

        if (.not. valid_alpha(alpha)) then
            mu = quiet_nan(0.0_dp)
        else
            mu = alpha / sum(alpha)
        end if
    end function dirichlet_mean

    function dirichlet_var(alpha) result(v)
        real(dp), intent(in) :: alpha(:) !! positive concentration parameters
        real(dp) :: v(size(alpha))

        real(dp) :: alpha0

        if (.not. valid_alpha(alpha)) then
            v = quiet_nan(0.0_dp)
            return
        end if
        alpha0 = sum(alpha)
        v = alpha * (alpha0 - alpha) / (alpha0 * alpha0 * (alpha0 + 1.0_dp))
    end function dirichlet_var

    function dirichlet_cov(alpha) result(cov)
        real(dp), intent(in) :: alpha(:) !! positive concentration parameters
        real(dp) :: cov(size(alpha), size(alpha))

        integer :: i
        integer :: j
        real(dp) :: alpha0
        real(dp) :: a(size(alpha))

        if (.not. valid_alpha(alpha)) then
            cov = quiet_nan(0.0_dp)
            return
        end if
        alpha0 = sum(alpha)
        a = alpha / alpha0
        do i = 1, size(alpha)
            do j = 1, size(alpha)
                cov(i, j) = -a(i) * a(j) / (alpha0 + 1.0_dp)
            end do
            cov(i, i) = cov(i, i) + a(i) / (alpha0 + 1.0_dp)
        end do
    end function dirichlet_cov

    function dirichlet_entropy(alpha) result(h)
        real(dp), intent(in) :: alpha(:) !! positive concentration parameters
        real(dp) :: h

        integer :: k
        real(dp) :: alpha0
        real(dp) :: log_beta

        if (.not. valid_alpha(alpha)) then
            h = quiet_nan(0.0_dp)
            return
        end if
        k = size(alpha)
        alpha0 = sum(alpha)
        log_beta = sum(gammaln(alpha)) - gammaln(alpha0)
        h = log_beta + (alpha0 - real(k, dp)) * digamma(alpha0) - &
            sum((alpha - 1.0_dp) * digamma(alpha))
    end function dirichlet_entropy

    subroutine dirichlet_rvs(state, sample, alpha, status)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by the draw
        real(dp), intent(out) :: sample(:) !! generated simplex vector
        real(dp), intent(in) :: alpha(:) !! positive concentration parameters
        integer, intent(out), optional :: status !! zero on success; nonzero on invalid input

        integer :: i
        real(dp) :: total

        if (size(sample) /= size(alpha)) then
            sample = quiet_nan(0.0_dp)
            if (present(status)) status = dirichlet_status_invalid_shape
            return
        end if
        if (.not. valid_alpha(alpha)) then
            sample = quiet_nan(0.0_dp)
            if (present(status)) status = dirichlet_status_invalid_parameter
            return
        end if

        do i = 1, size(alpha)
            sample(i) = gamma_ppf(rng_uniform(state), alpha(i))
        end do
        total = sum(sample)
        if (.not. ieee_is_finite(total) .or. total <= 0.0_dp) then
            sample = quiet_nan(0.0_dp)
            if (present(status)) status = dirichlet_status_invalid_parameter
            return
        end if
        sample = sample / total
        if (present(status)) status = dirichlet_status_success
    end subroutine dirichlet_rvs

    subroutine dirichlet_rvs_array(state, samples, alpha, status)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by all draws
        real(dp), intent(out) :: samples(:, :) !! samples by column, shape (K, n_samples)
        real(dp), intent(in) :: alpha(:) !! positive concentration parameters
        integer, intent(out), optional :: status !! zero on success; nonzero on invalid input

        integer :: info
        integer :: j

        if (size(samples, 1) /= size(alpha)) then
            samples = quiet_nan(0.0_dp)
            if (present(status)) status = dirichlet_status_invalid_shape
            return
        end if
        do j = 1, size(samples, 2)
            call dirichlet_rvs(state, samples(:, j), alpha, info)
            if (info /= dirichlet_status_success) then
                if (j < size(samples, 2)) samples(:, j + 1:) = quiet_nan(0.0_dp)
                if (present(status)) status = info
                return
            end if
        end do
        if (present(status)) status = dirichlet_status_success
    end subroutine dirichlet_rvs_array

    pure logical function valid_alpha(alpha) result(valid)
        real(dp), intent(in) :: alpha(:) !! candidate concentration vector

        valid = size(alpha) >= 1 .and. all(ieee_is_finite(alpha)) .and. all(alpha > 0.0_dp)
    end function valid_alpha

end module scifort_dirichlet
