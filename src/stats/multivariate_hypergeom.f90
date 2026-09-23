! SPDX-License-Identifier: MIT AND BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
! Adapted portions: SciPy 1.17.0 scipy/stats/_multivariate.py (sequential sampling).
! Retained BSD terms: THIRD_PARTY_LICENSES.md; details: CODE_PROVENANCE.md.
!
! Multivariate hypergeometric distribution. Formulas and the sequential
! hypergeometric sampling decomposition were checked against SciPy 1.17.0.

module scifort_multivariate_hypergeom
    use scifort_hypergeom, only : hypergeom_ppf
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, quiet_nan
    use scifort_random, only : rng_state, rng_uniform
    use scifort_special_elementary, only : gammaln
    implicit none
    private

    integer, parameter, public :: multivariate_hypergeom_status_success = 0
    integer, parameter, public :: multivariate_hypergeom_status_invalid_parameter = 1
    integer, parameter, public :: multivariate_hypergeom_status_invalid_shape = 2

    public :: multivariate_hypergeom_cov
    public :: multivariate_hypergeom_logpmf
    public :: multivariate_hypergeom_mean
    public :: multivariate_hypergeom_pmf
    public :: multivariate_hypergeom_rvs
    public :: multivariate_hypergeom_rvs_array
    public :: multivariate_hypergeom_var

contains

    function multivariate_hypergeom_logpmf(x, m, n) result(y)
        integer, intent(in) :: x(:) !! observed counts by category
        integer, intent(in) :: m(:) !! population counts by category
        integer, intent(in) :: n !! number of objects sampled without replacement
        real(dp) :: y

        integer :: i
        integer :: total

        if (size(m) < 1 .or. size(x) /= size(m)) then
            y = quiet_nan(0.0_dp)
            return
        end if
        if (any(m < 0) .or. n < 0) then
            y = quiet_nan(0.0_dp)
            return
        end if
        total = sum(m)
        if (n > total) then
            y = quiet_nan(0.0_dp)
            return
        end if
        if (any(x < 0) .or. any(x > m) .or. sum(x) /= n) then
            y = negative_infinity(0.0_dp)
            return
        end if

        y = -log_binomial(total, n)
        do i = 1, size(m)
            y = y + log_binomial(m(i), x(i))
        end do
    end function multivariate_hypergeom_logpmf

    function multivariate_hypergeom_pmf(x, m, n) result(y)
        integer, intent(in) :: x(:) !! observed counts by category
        integer, intent(in) :: m(:) !! population counts by category
        integer, intent(in) :: n !! number of objects sampled without replacement
        real(dp) :: y

        y = exp(multivariate_hypergeom_logpmf(x, m, n))
    end function multivariate_hypergeom_pmf

    function multivariate_hypergeom_mean(m, n) result(mu)
        integer, intent(in) :: m(:) !! population counts by category
        integer, intent(in) :: n !! number of objects sampled without replacement
        real(dp) :: mu(size(m))

        integer :: total

        if (.not. valid_parameters(m, n)) then
            mu = quiet_nan(0.0_dp)
            return
        end if
        total = sum(m)
        if (total == 0) then
            mu = 0.0_dp
        else
            mu = real(n, dp) * real(m, dp) / real(total, dp)
        end if
    end function multivariate_hypergeom_mean

    function multivariate_hypergeom_var(m, n) result(v)
        integer, intent(in) :: m(:) !! population counts by category
        integer, intent(in) :: n !! number of objects sampled without replacement
        real(dp) :: v(size(m))

        real(dp) :: nn
        real(dp) :: total
        real(dp) :: mi(size(m))

        if (.not. valid_parameters(m, n)) then
            v = quiet_nan(0.0_dp)
            return
        end if
        total = real(sum(m), dp)
        if (total <= 1.0_dp) then
            v = 0.0_dp
            return
        end if
        nn = real(n, dp)
        mi = real(m, dp)
        v = nn * (mi / total) * ((total - mi) / total) * &
            ((total - nn) / (total - 1.0_dp))
    end function multivariate_hypergeom_var

    function multivariate_hypergeom_cov(m, n) result(cov)
        integer, intent(in) :: m(:) !! population counts by category
        integer, intent(in) :: n !! number of objects sampled without replacement
        real(dp) :: cov(size(m), size(m))

        integer :: i
        integer :: j
        real(dp) :: factor
        real(dp) :: mi(size(m))
        real(dp) :: nn
        real(dp) :: total

        if (.not. valid_parameters(m, n)) then
            cov = quiet_nan(0.0_dp)
            return
        end if
        total = real(sum(m), dp)
        if (total <= 1.0_dp) then
            cov = 0.0_dp
            return
        end if
        nn = real(n, dp)
        mi = real(m, dp)
        factor = nn * (total - nn) / ((total - 1.0_dp) * total * total)
        do j = 1, size(m)
            do i = 1, size(m)
                cov(i, j) = -factor * mi(i) * mi(j)
            end do
            cov(j, j) = factor * mi(j) * (total - mi(j))
        end do
    end function multivariate_hypergeom_cov

    subroutine multivariate_hypergeom_rvs(state, m, n, sample, status)
        type(rng_state), intent(inout) :: state !! explicit random-number generator state
        integer, intent(in) :: m(:) !! population counts by category
        integer, intent(in) :: n !! number of objects sampled without replacement
        integer, intent(out) :: sample(:) !! sampled category counts
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid input

        integer :: draw
        integer :: i
        integer :: remaining_draws
        integer :: remaining_population
        integer :: total

        sample = 0
        if (size(sample) /= size(m) .or. size(m) < 1) then
            if (present(status)) status = multivariate_hypergeom_status_invalid_shape
            return
        end if
        if (.not. valid_parameters(m, n)) then
            if (present(status)) status = multivariate_hypergeom_status_invalid_parameter
            return
        end if

        total = sum(m)
        remaining_population = total
        remaining_draws = n
        do i = 1, size(m) - 1
            if (remaining_draws == 0) then
                draw = 0
            else
                draw = int(hypergeom_ppf(rng_uniform(state), real(remaining_population, dp), &
                    real(m(i), dp), real(remaining_draws, dp)))
            end if
            sample(i) = draw
            remaining_draws = remaining_draws - draw
            remaining_population = remaining_population - m(i)
        end do
        sample(size(m)) = remaining_draws
        if (present(status)) status = multivariate_hypergeom_status_success
    end subroutine multivariate_hypergeom_rvs

    subroutine multivariate_hypergeom_rvs_array(state, m, n, samples, status)
        type(rng_state), intent(inout) :: state !! explicit random-number generator state
        integer, intent(in) :: m(:) !! population counts by category
        integer, intent(in) :: n !! number of objects sampled without replacement
        integer, intent(out) :: samples(:, :) !! samples by column, shape (categories, nsamples)
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid input

        integer :: j
        integer :: local_status

        samples = 0
        if (size(samples, 1) /= size(m) .or. size(m) < 1) then
            if (present(status)) status = multivariate_hypergeom_status_invalid_shape
            return
        end if
        if (.not. valid_parameters(m, n)) then
            if (present(status)) status = multivariate_hypergeom_status_invalid_parameter
            return
        end if
        do j = 1, size(samples, 2)
            call multivariate_hypergeom_rvs(state, m, n, samples(:, j), local_status)
            if (local_status /= multivariate_hypergeom_status_success) then
                if (present(status)) status = local_status
                return
            end if
        end do
        if (present(status)) status = multivariate_hypergeom_status_success
    end subroutine multivariate_hypergeom_rvs_array

    pure logical function valid_parameters(m, n) result(valid)
        integer, intent(in) :: m(:) !! candidate population counts
        integer, intent(in) :: n !! candidate sample size

        valid = size(m) >= 1 .and. all(m >= 0) .and. n >= 0
        if (valid) valid = n <= sum(m)
    end function valid_parameters

    pure function log_binomial(n, k) result(y)
        integer, intent(in) :: n !! nonnegative total count
        integer, intent(in) :: k !! selected count in [0,n]
        real(dp) :: y

        y = gammaln(real(n + 1, dp)) - gammaln(real(k + 1, dp)) - &
            gammaln(real(n - k + 1, dp))
    end function log_binomial

end module scifort_multivariate_hypergeom
