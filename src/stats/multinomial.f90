! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_multinomial
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_binomial, only : binomial_pmf, binomial_ppf
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, quiet_nan
    use scifort_random, only : rng_state, rng_uniform
    use scifort_special_elementary, only : entr, gammaln, xlogy
    implicit none
    private

    integer, parameter, public :: multinomial_status_success = 0
    integer, parameter, public :: multinomial_status_invalid_parameter = 1
    integer, parameter, public :: multinomial_status_invalid_shape = 2

    public :: multinomial_cov
    public :: multinomial_entropy
    public :: multinomial_logpmf
    public :: multinomial_mean
    public :: multinomial_pmf
    public :: multinomial_rvs
    public :: multinomial_rvs_array

contains

    function multinomial_logpmf(x, n, p) result(y)
        real(dp), intent(in) :: x(:) !! category counts; each value must be a nonnegative integer and sum to n
        integer, intent(in) :: n !! total number of trials, >= 0
        real(dp), intent(in) :: p(:) !! category probabilities; SciPy-compatible final-entry normalization is applied
        real(dp) :: y

        real(dp), allocatable :: pp(:)
        integer :: status

        call process_probabilities(p, pp, status)
        if (n < 0 .or. status /= multinomial_status_success) then
            y = quiet_nan(0.0_dp)
            return
        end if
        if (size(x) /= size(pp)) then
            y = quiet_nan(0.0_dp)
            return
        end if
        if (.not. valid_counts(x, n)) then
            y = negative_infinity(0.0_dp)
            return
        end if

        y = gammaln(real(n + 1, dp)) + &
            sum(xlogy(x, pp) - gammaln(x + 1.0_dp))
    end function multinomial_logpmf

    function multinomial_pmf(x, n, p) result(y)
        real(dp), intent(in) :: x(:) !! category counts; each value must be a nonnegative integer and sum to n
        integer, intent(in) :: n !! total number of trials, >= 0
        real(dp), intent(in) :: p(:) !! category probabilities; SciPy-compatible final-entry normalization is applied
        real(dp) :: y

        y = exp(multinomial_logpmf(x, n, p))
    end function multinomial_pmf

    function multinomial_mean(n, p) result(mu)
        integer, intent(in) :: n !! total number of trials, >= 0
        real(dp), intent(in) :: p(:) !! category probabilities; SciPy-compatible final-entry normalization is applied
        real(dp) :: mu(size(p))

        real(dp), allocatable :: pp(:)
        integer :: status

        call process_probabilities(p, pp, status)
        if (n < 0 .or. status /= multinomial_status_success) then
            mu = quiet_nan(0.0_dp)
        else
            mu = real(n, dp) * pp
        end if
    end function multinomial_mean

    function multinomial_cov(n, p) result(cov)
        integer, intent(in) :: n !! total number of trials, >= 0
        real(dp), intent(in) :: p(:) !! category probabilities; SciPy-compatible final-entry normalization is applied
        real(dp) :: cov(size(p), size(p))

        real(dp), allocatable :: pp(:)
        integer :: i
        integer :: j
        integer :: status

        call process_probabilities(p, pp, status)
        if (n < 0 .or. status /= multinomial_status_success) then
            cov = quiet_nan(0.0_dp)
            return
        end if

        do i = 1, size(pp)
            do j = 1, size(pp)
                cov(i, j) = -real(n, dp) * pp(i) * pp(j)
            end do
            cov(i, i) = cov(i, i) + real(n, dp) * pp(i)
        end do
    end function multinomial_cov

    function multinomial_entropy(n, p) result(h)
        integer, intent(in) :: n !! total number of trials, >= 0
        real(dp), intent(in) :: p(:) !! category probabilities; SciPy-compatible final-entry normalization is applied
        real(dp) :: h

        real(dp), allocatable :: pp(:)
        integer :: i
        integer :: k
        integer :: status
        real(dp) :: term2

        call process_probabilities(p, pp, status)
        if (n < 0 .or. status /= multinomial_status_success) then
            h = quiet_nan(0.0_dp)
            return
        end if
        if (n == 0) then
            h = 0.0_dp
            return
        end if

        h = real(n, dp) * sum(entr(pp)) - gammaln(real(n + 1, dp))
        term2 = 0.0_dp
        do i = 1, size(pp)
            do k = 1, n
                term2 = term2 + binomial_pmf(k, n, pp(i)) * gammaln(real(k + 1, dp))
            end do
        end do
        h = h + term2
    end function multinomial_entropy

    subroutine multinomial_rvs(state, sample, n, p, status)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by the draw
        integer, intent(out) :: sample(:) !! generated category counts; sum equals n on success
        integer, intent(in) :: n !! total number of trials, >= 0
        real(dp), intent(in) :: p(:) !! category probabilities; SciPy-compatible final-entry normalization is applied
        integer, intent(out), optional :: status !! zero on success; nonzero on invalid input

        real(dp), allocatable :: pp(:)
        integer :: i
        integer :: info
        integer :: remaining
        real(dp) :: remaining_probability
        real(dp) :: conditional_probability

        call process_probabilities(p, pp, info)
        if (size(sample) /= size(p)) then
            sample = 0
            if (present(status)) status = multinomial_status_invalid_shape
            return
        end if
        if (n < 0 .or. info /= multinomial_status_success) then
            sample = 0
            if (present(status)) status = multinomial_status_invalid_parameter
            return
        end if

        sample = 0
        remaining = n
        remaining_probability = 1.0_dp
        do i = 1, size(pp) - 1
            if (remaining <= 0) exit
            if (remaining_probability <= 0.0_dp) exit
            conditional_probability = pp(i) / remaining_probability
            conditional_probability = max(0.0_dp, min(1.0_dp, conditional_probability))
            sample(i) = int(binomial_ppf(rng_uniform(state), real(remaining, dp), conditional_probability))
            remaining = remaining - sample(i)
            remaining_probability = remaining_probability - pp(i)
        end do
        sample(size(pp)) = remaining
        if (present(status)) status = multinomial_status_success
    end subroutine multinomial_rvs

    subroutine multinomial_rvs_array(state, samples, n, p, status)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by all draws
        integer, intent(out) :: samples(:, :) !! samples by column, shape (K, n_samples)
        integer, intent(in) :: n !! total number of trials, >= 0
        real(dp), intent(in) :: p(:) !! category probabilities; SciPy-compatible final-entry normalization is applied
        integer, intent(out), optional :: status !! zero on success; nonzero on invalid input

        integer :: info
        integer :: j

        if (size(samples, 1) /= size(p)) then
            samples = 0
            if (present(status)) status = multinomial_status_invalid_shape
            return
        end if
        do j = 1, size(samples, 2)
            call multinomial_rvs(state, samples(:, j), n, p, info)
            if (info /= multinomial_status_success) then
                if (j < size(samples, 2)) samples(:, j + 1:) = 0
                if (present(status)) status = info
                return
            end if
        end do
        if (present(status)) status = multinomial_status_success
    end subroutine multinomial_rvs_array

    subroutine process_probabilities(p, pp, status)
        real(dp), intent(in) :: p(:) !! input category probabilities
        real(dp), allocatable, intent(out) :: pp(:) !! adjusted category probabilities
        integer, intent(out) :: status !! zero on success; nonzero for invalid probabilities

        real(dp), parameter :: probability_tolerance = 10.0_dp * epsilon(1.0_dp)
        real(dp) :: p_sum

        allocate(pp(size(p)))
        if (size(p) < 1) then
            status = multinomial_status_invalid_shape
            return
        end if
        if (any(.not. ieee_is_finite(p))) then
            pp = quiet_nan(0.0_dp)
            status = multinomial_status_invalid_parameter
            return
        end if

        pp = p
        p_sum = sum(pp)
        if (abs(1.0_dp - p_sum) > probability_tolerance) then
            if (size(pp) == 1) then
                pp(1) = 1.0_dp
            else
                pp(size(pp)) = 1.0_dp - sum(pp(1:size(pp) - 1))
            end if
        end if
        if (any(pp < 0.0_dp) .or. any(pp > 1.0_dp)) then
            pp = quiet_nan(0.0_dp)
            status = multinomial_status_invalid_parameter
            return
        end if
        status = multinomial_status_success
    end subroutine process_probabilities

    pure logical function valid_counts(x, n) result(valid)
        real(dp), intent(in) :: x(:) !! candidate category counts
        integer, intent(in) :: n !! required total count

        valid = all(ieee_is_finite(x)) .and. all(x >= 0.0_dp) .and. &
            all(x == aint(x)) .and. sum(x) == real(n, dp)
    end function valid_counts

end module scifort_multinomial
