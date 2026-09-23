! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Distribution of contingency tables with fixed row and column marginals.
! The exact sampler draws each row from the corresponding multivariate
! hypergeometric conditional distribution, which is distributionally
! equivalent to the Boyett random-permutation construction.

module scifort_random_table
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, quiet_nan
    use scifort_multivariate_hypergeom, only : multivariate_hypergeom_rvs, &
        multivariate_hypergeom_status_success
    use scifort_random, only : rng_state
    implicit none
    private

    integer, parameter, public :: random_table_status_success = 0
    integer, parameter, public :: random_table_status_invalid_shape = 1
    integer, parameter, public :: random_table_status_invalid_parameter = 2
    integer, parameter, public :: random_table_status_numerical_failure = 3

    public :: random_table_logpmf
    public :: random_table_mean
    public :: random_table_pmf
    public :: random_table_rvs
    public :: random_table_rvs_array

contains

    function random_table_logpmf(x, row, col) result(y)
        integer, intent(in) :: x(:, :) !! nonnegative contingency table
        integer, intent(in) :: row(:) !! prescribed row sums, nonnegative integers
        integer, intent(in) :: col(:) !! prescribed column sums, nonnegative integers
        real(dp) :: y

        integer :: i
        integer :: j
        integer :: n

        if (.not. valid_marginals(row, col)) then
            y = quiet_nan(0.0_dp)
            return
        end if
        if (size(x,1) /= size(row) .or. size(x,2) /= size(col) .or. any(x < 0)) then
            y = quiet_nan(0.0_dp)
            return
        end if
        do i = 1, size(row)
            if (sum(x(i,:)) /= row(i)) then
                y = negative_infinity(0.0_dp)
                return
            end if
        end do
        do j = 1, size(col)
            if (sum(x(:,j)) /= col(j)) then
                y = negative_infinity(0.0_dp)
                return
            end if
        end do

        n = sum(row)
        y = -log_gamma(real(n + 1, dp))
        do i = 1, size(row)
            y = y + log_gamma(real(row(i) + 1, dp))
        end do
        do j = 1, size(col)
            y = y + log_gamma(real(col(j) + 1, dp))
        end do
        do j = 1, size(col)
            do i = 1, size(row)
                y = y - log_gamma(real(x(i,j) + 1, dp))
            end do
        end do
    end function random_table_logpmf

    function random_table_pmf(x, row, col) result(y)
        integer, intent(in) :: x(:, :) !! nonnegative contingency table
        integer, intent(in) :: row(:) !! prescribed row sums, nonnegative integers
        integer, intent(in) :: col(:) !! prescribed column sums, nonnegative integers
        real(dp) :: y

        real(dp) :: logp

        logp = random_table_logpmf(x, row, col)
        if (ieee_is_nan(logp)) then
            y = logp
        else
            y = exp(logp)
        end if
    end function random_table_pmf

    function random_table_mean(row, col) result(mu)
        integer, intent(in) :: row(:) !! prescribed row sums, nonnegative integers
        integer, intent(in) :: col(:) !! prescribed column sums, nonnegative integers
        real(dp) :: mu(size(row), size(col))

        integer :: i
        integer :: j
        integer :: n

        if (.not. valid_marginals(row, col)) then
            mu = quiet_nan(0.0_dp)
            return
        end if
        n = sum(row)
        if (n == 0) then
            mu = quiet_nan(0.0_dp)
            return
        end if
        do j = 1, size(col)
            do i = 1, size(row)
                mu(i,j) = real(row(i),dp) * real(col(j),dp) / real(n,dp)
            end do
        end do
    end function random_table_mean

    subroutine random_table_rvs(state, row, col, sample, status)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by conditional draws
        integer, intent(in) :: row(:) !! prescribed row sums, nonnegative integers
        integer, intent(in) :: col(:) !! prescribed column sums, nonnegative integers
        integer, intent(out) :: sample(:, :) !! sampled table with shape (size(row), size(col))
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid/numerical failure

        integer :: i
        integer :: info
        integer, allocatable :: remaining_col(:)
        integer, allocatable :: draw(:)

        sample = 0
        if (size(sample,1) /= size(row) .or. size(sample,2) /= size(col) .or. &
                size(row) < 1 .or. size(col) < 1) then
            if (present(status)) status = random_table_status_invalid_shape
            return
        end if
        if (.not. valid_marginals(row, col)) then
            if (present(status)) status = random_table_status_invalid_parameter
            return
        end if

        allocate(remaining_col(size(col)), draw(size(col)))
        remaining_col = col
        do i = 1, size(row) - 1
            call multivariate_hypergeom_rvs(state, remaining_col, row(i), draw, info)
            if (info /= multivariate_hypergeom_status_success) then
                if (present(status)) status = random_table_status_numerical_failure
                return
            end if
            sample(i,:) = draw
            remaining_col = remaining_col - draw
        end do
        sample(size(row),:) = remaining_col
        if (present(status)) status = random_table_status_success
    end subroutine random_table_rvs

    subroutine random_table_rvs_array(state, row, col, samples, status)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by all table draws
        integer, intent(in) :: row(:) !! prescribed row sums, nonnegative integers
        integer, intent(in) :: col(:) !! prescribed column sums, nonnegative integers
        integer, intent(out) :: samples(:, :, :) !! tables indexed by third dimension
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid/numerical failure

        integer :: info
        integer :: k

        samples = 0
        if (size(samples,1) /= size(row) .or. size(samples,2) /= size(col)) then
            if (present(status)) status = random_table_status_invalid_shape
            return
        end if
        do k = 1, size(samples,3)
            call random_table_rvs(state, row, col, samples(:,:,k), info)
            if (info /= random_table_status_success) then
                if (present(status)) status = info
                return
            end if
        end do
        if (present(status)) status = random_table_status_success
    end subroutine random_table_rvs_array

    pure logical function valid_marginals(row, col) result(valid)
        integer, intent(in) :: row(:) !! candidate row sums
        integer, intent(in) :: col(:) !! candidate column sums

        valid = size(row) >= 1 .and. size(col) >= 1
        if (.not. valid) return
        valid = all(row >= 0) .and. all(col >= 0) .and. sum(row) == sum(col)
    end function valid_marginals

end module scifort_random_table
