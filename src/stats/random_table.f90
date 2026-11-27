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
        !! Evaluate the log probability of a table with fixed margins.
        !! @dimension n_rows = size(row)
        !! @dimension n_cols = size(col)
        !! @require n_rows >= 1
        !! @require n_cols >= 1
        !! Nonempty margins must be nonnegative with equal totals <= huge(0)-1.
        integer, intent(in) :: x(:, :) !! shape (n_rows,n_cols); nonnegative contingency table
        integer, intent(in) :: row(:) !! shape (n_rows); nonnegative prescribed row sums
        integer, intent(in) :: col(:) !! shape (n_cols); nonnegative prescribed column sums
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
            if (bounded_count_total(x(i,:)) /= row(i)) then
                y = negative_infinity(0.0_dp)
                return
            end if
        end do
        do j = 1, size(col)
            if (bounded_count_total(x(:,j)) /= col(j)) then
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
        !! Evaluate the probability of a table with fixed margins.
        !! @dimension n_rows = size(row)
        !! @dimension n_cols = size(col)
        !! @require n_rows >= 1
        !! @require n_cols >= 1
        !! Nonempty margins must be nonnegative with equal totals <= huge(0)-1.
        integer, intent(in) :: x(:, :) !! shape (n_rows,n_cols); nonnegative contingency table
        integer, intent(in) :: row(:) !! shape (n_rows); nonnegative prescribed row sums
        integer, intent(in) :: col(:) !! shape (n_cols); nonnegative prescribed column sums
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
        !! Return expected counts for a table with fixed margins.
        !! @dimension n_rows = size(row)
        !! @dimension n_cols = size(col)
        !! @require n_rows >= 1
        !! @require n_cols >= 1
        !! Nonempty margins must be nonnegative with equal totals <= huge(0)-1.
        !! A zero total returns NaN in every entry.
        integer, intent(in) :: row(:) !! shape (n_rows); nonnegative prescribed row sums
        integer, intent(in) :: col(:) !! shape (n_cols); nonnegative prescribed column sums
        real(dp) :: mu(size(row), size(col)) !! shape (n_rows,n_cols); expected table counts

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
        !! Draw a fixed-margin table using explicit RNG state.
        !! @dimension n_rows = size(row)
        !! @dimension n_cols = size(col)
        !! @require n_rows >= 1
        !! @require n_cols >= 1
        !! Nonempty margins must be nonnegative with equal totals <= huge(0)-1.
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by conditional draws
        integer, intent(in) :: row(:) !! shape (n_rows); nonnegative prescribed row sums
        integer, intent(in) :: col(:) !! shape (n_cols); nonnegative prescribed column sums
        integer, intent(out) :: sample(:, :) !! shape (n_rows,n_cols); sampled table
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
        !! Draw fixed-margin tables using explicit RNG state.
        !! @dimension n_rows = size(row)
        !! @dimension n_cols = size(col)
        !! @require n_rows >= 0
        !! @require n_cols >= 0
        !! @dimension n_samples = size(samples,3)
        !! @require n_samples >= 0
        !! An empty batch skips native marginal validation and leaves state unchanged.
        !! Nonempty margins must be nonnegative with equal totals <= huge(0)-1.
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by all table draws
        integer, intent(in) :: row(:) !! shape (n_rows); nonnegative prescribed row sums
        integer, intent(in) :: col(:) !! shape (n_cols); nonnegative prescribed column sums
        integer, intent(out) :: samples(:, :, :) !! shape (n_rows,n_cols,n_samples); tables along the third axis
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
        !! Validate margins without overflowing default-integer sums.
        !! n_rows = size(row), n_cols = size(col); both must be >= 1 for valid margins.
        !! Nonnegative margins must have equal totals <= huge(0)-1.
        integer, intent(in) :: row(:) !! shape (n_rows); candidate row sums
        integer, intent(in) :: col(:) !! shape (n_cols); candidate column sums

        integer :: row_total
        integer :: col_total

        valid = .false.
        if (size(row) < 1 .or. size(col) < 1) return
        row_total = bounded_count_total(row)
        col_total = bounded_count_total(col)
        valid = row_total >= 0 .and. col_total >= 0 .and. row_total == col_total
    end function valid_marginals

    pure integer function bounded_count_total(counts) result(total)
        !! Sum nonnegative counts up to huge(0)-1; return -1 for invalid or excessive counts.
        !! The vector length is unrestricted, including zero.
        integer, intent(in) :: counts(:) !! nonnegative counts to sum; empty input has total zero

        integer :: i
        integer, parameter :: maximum = huge(0)-1

        total = 0
        do i = 1, size(counts)
            if (counts(i) < 0) then
                total = -1
                return
            end if
            ! Check the remaining capacity before adding, so the addition cannot overflow.
            if (counts(i) > maximum-total) then
                total = -1
                return
            end if
            total = total + counts(i)
        end do
    end function bounded_count_total

end module scifort_random_table
