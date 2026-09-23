! SPDX-License-Identifier: MIT AND BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
! Adapted portions: SciPy 1.17.0 scipy/stats/_multivariate.py (sampling decomposition).
! Retained BSD terms: THIRD_PARTY_LICENSES.md; details: CODE_PROVENANCE.md.
!
! Matrix-variate Student t distribution following scipy.stats.matrix_t.

module scifort_matrix_t
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_invwishart, only : invwishart_rvs, invwishart_status_success
    use scifort_kinds, only : dp
    use scifort_matrix_distribution_helpers, only : log_multivariate_gamma, matrix_status_success, &
        solve_lower_matrix, spd_factor
    use scifort_math, only : quiet_nan
    use scifort_normal, only : normal_ppf
    use scifort_random, only : rng_state, rng_uniform
    implicit none
    private

    integer, parameter, public :: matrix_t_status_success = 0
    integer, parameter, public :: matrix_t_status_invalid_shape = 1
    integer, parameter, public :: matrix_t_status_invalid_parameter = 2
    integer, parameter, public :: matrix_t_status_linalg_failure = 3

    public :: matrix_t_logpdf
    public :: matrix_t_pdf
    public :: matrix_t_rvs
    public :: matrix_t_rvs_array

contains

    function matrix_t_logpdf(x, mean, row_spread, col_spread, df) result(y)
        real(dp), intent(in) :: x(:, :) !! matrix variate
        real(dp), intent(in), optional :: mean(:, :) !! mean matrix; zero when absent
        real(dp), intent(in), optional :: row_spread(:, :) !! SPD row spread; identity when absent
        real(dp), intent(in), optional :: col_spread(:, :) !! SPD column spread; identity when absent
        real(dp), intent(in), optional :: df !! positive degrees of freedom; default 1
        real(dp) :: y

        integer :: i
        integer :: info
        integer :: m
        integer :: n
        real(dp), allocatable :: ccol(:, :)
        real(dp), allocatable :: crow(:, :)
        real(dp), allocatable :: col(:, :)
        real(dp), allocatable :: dev(:, :)
        real(dp), allocatable :: whitened_row(:, :)
        real(dp), allocatable :: whitened_both(:, :)
        real(dp), allocatable :: mean0(:, :)
        real(dp), allocatable :: row(:, :)
        real(dp), allocatable :: detarg(:, :)
        real(dp), allocatable :: ident(:, :)
        real(dp) :: dof
        real(dp) :: logdet_arg
        real(dp) :: logdet_col
        real(dp) :: logdet_row

        m = size(x,1)
        n = size(x,2)
        if (m < 1 .or. n < 1 .or. .not. all(ieee_is_finite(x))) then
            y = quiet_nan(0.0_dp)
            return
        end if
        dof = 1.0_dp
        if (present(df)) dof = df
        if (.not. ieee_is_finite(dof) .or. dof <= 0.0_dp) then
            y = quiet_nan(0.0_dp)
            return
        end if
        allocate(mean0(m,n), row(m,m), col(n,n), crow(m,m), ccol(n,n), dev(m,n), &
            whitened_row(m,n), whitened_both(n,m), detarg(n,n), ident(n,n))
        call process_parameters(m, n, mean, row_spread, col_spread, mean0, row, col, info)
        if (info /= matrix_t_status_success) then
            y = quiet_nan(0.0_dp)
            return
        end if
        call spd_factor(row, crow, logdet_row, info)
        if (info /= matrix_status_success) then
            y = quiet_nan(0.0_dp)
            return
        end if
        call spd_factor(col, ccol, logdet_col, info)
        if (info /= matrix_status_success) then
            y = quiet_nan(0.0_dp)
            return
        end if

        dev = x - mean0
        call solve_lower_matrix(crow, dev, whitened_row, info)
        if (info /= matrix_status_success) then
            y = quiet_nan(0.0_dp)
            return
        end if
        call solve_lower_matrix(ccol, transpose(whitened_row), whitened_both, info)
        if (info /= matrix_status_success) then
            y = quiet_nan(0.0_dp)
            return
        end if
        ident = 0.0_dp
        do i = 1, n
            ident(i,i) = 1.0_dp
        end do
        ! Similarity transform of SciPy's I + D^T R^-1 D C^-1 into
        ! the symmetric positive-definite I + B B^T.
        detarg = ident + matmul(whitened_both, transpose(whitened_both))
        call spd_factor(detarg, ccol, logdet_arg, info)
        if (info /= matrix_status_success) then
            y = quiet_nan(0.0_dp)
            return
        end if

        y = -0.5_dp * (dof + real(m + n - 1,dp)) * logdet_arg + &
            log_multivariate_gamma(0.5_dp * (dof + real(m + n - 1,dp)), n) - &
            log_multivariate_gamma(0.5_dp * (dof + real(n - 1,dp)), n) - &
            0.5_dp * real(m*n,dp) * log(acos(-1.0_dp)) - &
            0.5_dp * real(n,dp) * logdet_row - 0.5_dp * real(m,dp) * logdet_col
    end function matrix_t_logpdf

    function matrix_t_pdf(x, mean, row_spread, col_spread, df) result(y)
        real(dp), intent(in) :: x(:, :) !! matrix variate
        real(dp), intent(in), optional :: mean(:, :) !! mean matrix
        real(dp), intent(in), optional :: row_spread(:, :) !! row spread matrix
        real(dp), intent(in), optional :: col_spread(:, :) !! column spread matrix
        real(dp), intent(in), optional :: df !! degrees of freedom
        real(dp) :: y

        y = exp(matrix_t_logpdf(x, mean, row_spread, col_spread, df))
    end function matrix_t_pdf

    subroutine matrix_t_rvs(state, sample, mean, row_spread, col_spread, df, status)
        type(rng_state), intent(inout) :: state !! explicit random-number generator state
        real(dp), intent(out) :: sample(:, :) !! random matrix variate
        real(dp), intent(in), optional :: mean(:, :) !! mean matrix; zero when absent
        real(dp), intent(in), optional :: row_spread(:, :) !! SPD row spread; identity when absent
        real(dp), intent(in), optional :: col_spread(:, :) !! SPD column spread; identity when absent
        real(dp), intent(in), optional :: df !! positive degrees of freedom; default 1
        integer, intent(out), optional :: status !! matrix_t_status_* result

        integer :: i
        integer :: info
        integer :: j
        integer :: m
        integer :: n
        real(dp), allocatable :: ccol(:, :)
        real(dp), allocatable :: crow(:, :)
        real(dp), allocatable :: col(:, :)
        real(dp), allocatable :: iw(:, :)
        real(dp), allocatable :: mean0(:, :)
        real(dp), allocatable :: row(:, :)
        real(dp), allocatable :: z(:, :)
        real(dp) :: dof
        real(dp) :: dummy

        m = size(sample,1)
        n = size(sample,2)
        if (m < 1 .or. n < 1) then
            if (present(status)) status = matrix_t_status_invalid_shape
            return
        end if
        dof = 1.0_dp
        if (present(df)) dof = df
        if (.not. ieee_is_finite(dof) .or. dof <= 0.0_dp) then
            sample = quiet_nan(0.0_dp)
            if (present(status)) status = matrix_t_status_invalid_parameter
            return
        end if
        allocate(mean0(m,n), row(m,m), col(n,n), crow(m,m), ccol(n,n), iw(m,m), z(m,n))
        call process_parameters(m, n, mean, row_spread, col_spread, mean0, row, col, info)
        if (info /= matrix_t_status_success) then
            sample = quiet_nan(0.0_dp)
            if (present(status)) status = info
            return
        end if

        ! Matrix-t mixture: N(M, W^{-1}, col_spread),
        ! W^{-1} ~ InvWishart(df + m - 1, row_spread).
        call invwishart_rvs(state, dof + real(m - 1,dp), row, iw, info)
        if (info /= invwishart_status_success) then
            sample = quiet_nan(0.0_dp)
            if (present(status)) status = matrix_t_status_linalg_failure
            return
        end if
        call spd_factor(iw, crow, dummy, info)
        if (info /= matrix_status_success) then
            sample = quiet_nan(0.0_dp)
            if (present(status)) status = matrix_t_status_linalg_failure
            return
        end if
        call spd_factor(col, ccol, dummy, info)
        if (info /= matrix_status_success) then
            sample = quiet_nan(0.0_dp)
            if (present(status)) status = matrix_t_status_linalg_failure
            return
        end if
        do j = 1, n
            do i = 1, m
                z(i,j) = normal_ppf(rng_uniform(state))
            end do
        end do
        sample = mean0 + matmul(matmul(crow, z), transpose(ccol))
        if (present(status)) status = matrix_t_status_success
    end subroutine matrix_t_rvs

    subroutine matrix_t_rvs_array(state, samples, mean, row_spread, col_spread, df, status)
        type(rng_state), intent(inout) :: state !! explicit random-number generator state
        real(dp), intent(out) :: samples(:, :, :) !! samples indexed by third dimension
        real(dp), intent(in), optional :: mean(:, :) !! mean matrix
        real(dp), intent(in), optional :: row_spread(:, :) !! row spread matrix
        real(dp), intent(in), optional :: col_spread(:, :) !! column spread matrix
        real(dp), intent(in), optional :: df !! degrees of freedom
        integer, intent(out), optional :: status !! matrix_t_status_* result

        integer :: info
        integer :: k

        do k = 1, size(samples,3)
            call matrix_t_rvs(state, samples(:,:,k), mean, row_spread, col_spread, df, info)
            if (info /= matrix_t_status_success) then
                if (present(status)) status = info
                return
            end if
        end do
        if (present(status)) status = matrix_t_status_success
    end subroutine matrix_t_rvs_array

    pure subroutine process_parameters(m, n, mean, row_spread, col_spread, mean0, row, col, status)
        integer, intent(in) :: m !! matrix row dimension
        integer, intent(in) :: n !! matrix column dimension
        real(dp), intent(in), optional :: mean(:, :) !! optional mean matrix
        real(dp), intent(in), optional :: row_spread(:, :) !! optional row spread matrix
        real(dp), intent(in), optional :: col_spread(:, :) !! optional column spread matrix
        real(dp), intent(out) :: mean0(:, :) !! resolved mean
        real(dp), intent(out) :: row(:, :) !! resolved row spread
        real(dp), intent(out) :: col(:, :) !! resolved column spread
        integer, intent(out) :: status !! matrix_t_status_* result

        integer :: i

        if (m < 1 .or. n < 1 .or. size(mean0,1) /= m .or. size(mean0,2) /= n .or. &
                size(row,1) /= m .or. size(row,2) /= m .or. size(col,1) /= n .or. size(col,2) /= n) then
            status = matrix_t_status_invalid_shape
            return
        end if
        mean0 = 0.0_dp
        if (present(mean)) then
            if (size(mean,1) /= m .or. size(mean,2) /= n) then
                status = matrix_t_status_invalid_shape
                return
            end if
            if (.not. all(ieee_is_finite(mean))) then
                status = matrix_t_status_invalid_parameter
                return
            end if
            mean0 = mean
        end if
        row = 0.0_dp
        if (present(row_spread)) then
            if (size(row_spread,1) /= m .or. size(row_spread,2) /= m) then
                status = matrix_t_status_invalid_shape
                return
            end if
            if (.not. all(ieee_is_finite(row_spread))) then
                status = matrix_t_status_invalid_parameter
                return
            end if
            row = row_spread
        else
            do i = 1, m
                row(i,i) = 1.0_dp
            end do
        end if
        col = 0.0_dp
        if (present(col_spread)) then
            if (size(col_spread,1) /= n .or. size(col_spread,2) /= n) then
                status = matrix_t_status_invalid_shape
                return
            end if
            if (.not. all(ieee_is_finite(col_spread))) then
                status = matrix_t_status_invalid_parameter
                return
            end if
            col = col_spread
        else
            do i = 1, n
                col(i,i) = 1.0_dp
            end do
        end if
        status = matrix_t_status_success
    end subroutine process_parameters

end module scifort_matrix_t
