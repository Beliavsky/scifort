! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_matrix_normal
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_matrix_distribution_helpers, only : matrix_status_success, solve_lower_matrix, spd_factor
    use scifort_math, only : quiet_nan
    use scifort_normal, only : normal_ppf
    use scifort_random, only : rng_state, rng_uniform
    implicit none
    private

    integer, parameter, public :: matrix_normal_status_success = 0
    integer, parameter, public :: matrix_normal_status_invalid_shape = 1
    integer, parameter, public :: matrix_normal_status_invalid_parameter = 2
    integer, parameter, public :: matrix_normal_status_linalg_failure = 3

    public :: matrix_normal_entropy
    public :: matrix_normal_logpdf
    public :: matrix_normal_pdf
    public :: matrix_normal_rvs
    public :: matrix_normal_rvs_array

contains

    function matrix_normal_logpdf(x, mean, rowcov, colcov) result(y)
        !! Evaluate the log density at one rectangular matrix.
        !! @dimension m = size(x,1)
        !! @dimension n = size(x,2)
        !! @require m >= 1
        !! @require n >= 1
        !! Mean defaults to zero; row and column matrices default to identities.
        !! All entries must be finite; covariance/spread factors use lower triangles.
        real(dp), intent(in) :: x(:, :) !! shape (m,n); finite evaluation matrix
        real(dp), intent(in), optional :: mean(:, :) !! shape (m,n); finite mean; zero when absent
        real(dp), intent(in), optional :: rowcov(:, :) !! shape (m,m); finite SPD row covariance; identity when absent
        real(dp), intent(in), optional :: colcov(:, :) !! shape (n,n); finite SPD column covariance; identity when absent
        real(dp) :: y

        integer :: info
        integer :: m
        integer :: n
        real(dp), allocatable :: a(:, :)
        real(dp), allocatable :: b(:, :)
        real(dp), allocatable :: cm(:, :)
        real(dp), allocatable :: dev(:, :)
        real(dp), allocatable :: lr(:, :)
        real(dp), allocatable :: lc(:, :)
        real(dp), allocatable :: mu(:, :)
        real(dp), allocatable :: rm(:, :)
        real(dp) :: logdet_c
        real(dp) :: logdet_r
        real(dp) :: maha

        m = size(x, 1)
        n = size(x, 2)
        if (m < 1 .or. n < 1 .or. .not. all(ieee_is_finite(x))) then
            y = quiet_nan(0.0_dp)
            return
        end if
        allocate(mu(m,n), rm(m,m), cm(n,n), lr(m,m), lc(n,n), dev(m,n), a(m,n), b(n,m))
        call process_parameters(m, n, mean, rowcov, colcov, mu, rm, cm, info)
        if (info /= matrix_normal_status_success) then
            y = quiet_nan(0.0_dp)
            return
        end if
        call spd_factor(rm, lr, logdet_r, info)
        if (info /= matrix_status_success) then
            y = quiet_nan(0.0_dp)
            return
        end if
        call spd_factor(cm, lc, logdet_c, info)
        if (info /= matrix_status_success) then
            y = quiet_nan(0.0_dp)
            return
        end if

        dev = x - mu
        call solve_lower_matrix(lr, dev, a, info)
        if (info /= matrix_status_success) then
            y = quiet_nan(0.0_dp)
            return
        end if
        call solve_lower_matrix(lc, transpose(a), b, info)
        if (info /= matrix_status_success) then
            y = quiet_nan(0.0_dp)
            return
        end if
        maha = sum(b * b)
        y = -0.5_dp * (real(m*n,dp) * log(2.0_dp * acos(-1.0_dp)) + &
            real(n,dp) * logdet_r + real(m,dp) * logdet_c + maha)
    end function matrix_normal_logpdf

    function matrix_normal_pdf(x, mean, rowcov, colcov) result(y)
        !! Evaluate the density at one rectangular matrix.
        !! @dimension m = size(x,1)
        !! @dimension n = size(x,2)
        !! @require m >= 1
        !! @require n >= 1
        !! Mean defaults to zero; row and column matrices default to identities.
        !! All entries must be finite; covariance/spread factors use lower triangles.
        real(dp), intent(in) :: x(:, :) !! shape (m,n); finite evaluation matrix
        real(dp), intent(in), optional :: mean(:, :) !! shape (m,n); finite mean; zero when absent
        real(dp), intent(in), optional :: rowcov(:, :) !! shape (m,m); finite SPD row covariance; identity when absent
        real(dp), intent(in), optional :: colcov(:, :) !! shape (n,n); finite SPD column covariance; identity when absent
        real(dp) :: y

        y = exp(matrix_normal_logpdf(x, mean, rowcov, colcov))
    end function matrix_normal_pdf

    function matrix_normal_entropy(rowcov, colcov) result(y)
        !! Return the matrix-normal differential entropy.
        !! @dimension m = size(rowcov,1)
        !! @dimension n = size(colcov,1)
        !! @require m >= 1
        !! @require n >= 1
        !! When rowcov is absent m = 1; when colcov is absent n = 1.
        !! Each absent covariance independently defaults to a 1-by-1 identity.
        real(dp), intent(in), optional :: rowcov(:, :) !! shape (m,m); finite SPD row covariance; identity when absent
        real(dp), intent(in), optional :: colcov(:, :) !! shape (n,n); finite SPD column covariance; identity when absent
        real(dp) :: y

        integer :: info
        integer :: m
        integer :: n
        real(dp), allocatable :: cm(:, :)
        real(dp), allocatable :: lc(:, :)
        real(dp), allocatable :: lr(:, :)
        real(dp), allocatable :: rm(:, :)
        real(dp) :: logdet_c
        real(dp) :: logdet_r

        m = 1
        n = 1
        if (present(rowcov)) m = size(rowcov, 1)
        if (present(colcov)) n = size(colcov, 1)
        if (m < 1 .or. n < 1) then
            y = quiet_nan(0.0_dp)
            return
        end if
        allocate(rm(m,m), cm(n,n), lr(m,m), lc(n,n))
        call covariance_or_identity(rowcov, m, rm, info)
        if (info /= matrix_normal_status_success) then
            y = quiet_nan(0.0_dp)
            return
        end if
        call covariance_or_identity(colcov, n, cm, info)
        if (info /= matrix_normal_status_success) then
            y = quiet_nan(0.0_dp)
            return
        end if
        call spd_factor(rm, lr, logdet_r, info)
        if (info /= matrix_status_success) then
            y = quiet_nan(0.0_dp)
            return
        end if
        call spd_factor(cm, lc, logdet_c, info)
        if (info /= matrix_status_success) then
            y = quiet_nan(0.0_dp)
            return
        end if
        y = 0.5_dp * real(m*n,dp) * (1.0_dp + log(2.0_dp * acos(-1.0_dp))) + &
            0.5_dp * real(n,dp) * logdet_r + 0.5_dp * real(m,dp) * logdet_c
    end function matrix_normal_entropy

    subroutine matrix_normal_rvs(state, sample, mean, rowcov, colcov, status)
        !! Draw one matrix using explicit RNG state.
        !! m = size(sample,1), n = size(sample,2); m >= 1 and n >= 1.
        type(rng_state), intent(inout) :: state !! explicit random-number generator state
        real(dp), intent(out) :: sample(:, :) !! shape (m,n); output matrix variate
        real(dp), intent(in), optional :: mean(:, :) !! shape (m,n); finite mean; zero when absent
        real(dp), intent(in), optional :: rowcov(:, :) !! shape (m,m); finite SPD row covariance; identity when absent
        real(dp), intent(in), optional :: colcov(:, :) !! shape (n,n); finite SPD column covariance; identity when absent
        integer, intent(out), optional :: status !! matrix_normal_status_* result

        integer :: i
        integer :: info
        integer :: j
        integer :: m
        integer :: n
        real(dp), allocatable :: cm(:, :)
        real(dp), allocatable :: lc(:, :)
        real(dp), allocatable :: lr(:, :)
        real(dp), allocatable :: mu(:, :)
        real(dp), allocatable :: rm(:, :)
        real(dp), allocatable :: z(:, :)
        real(dp) :: dummy

        m = size(sample, 1)
        n = size(sample, 2)
        if (m < 1 .or. n < 1) then
            if (present(status)) status = matrix_normal_status_invalid_shape
            return
        end if
        allocate(mu(m,n), rm(m,m), cm(n,n), lr(m,m), lc(n,n), z(m,n))
        call process_parameters(m, n, mean, rowcov, colcov, mu, rm, cm, info)
        if (info /= matrix_normal_status_success) then
            sample = quiet_nan(0.0_dp)
            if (present(status)) status = info
            return
        end if
        call spd_factor(rm, lr, dummy, info)
        if (info /= matrix_status_success) then
            sample = quiet_nan(0.0_dp)
            if (present(status)) status = matrix_normal_status_linalg_failure
            return
        end if
        call spd_factor(cm, lc, dummy, info)
        if (info /= matrix_status_success) then
            sample = quiet_nan(0.0_dp)
            if (present(status)) status = matrix_normal_status_linalg_failure
            return
        end if
        do j = 1, n
            do i = 1, m
                z(i,j) = normal_ppf(rng_uniform(state))
            end do
        end do
        sample = mu + matmul(matmul(lr, z), transpose(lc))
        if (present(status)) status = matrix_normal_status_success
    end subroutine matrix_normal_rvs

    subroutine matrix_normal_rvs_array(state, samples, mean, rowcov, colcov, status)
        !! Draw matrix samples along the third axis using explicit RNG state.
        !! @dimension m = size(samples,1)
        !! @dimension n = size(samples,2)
        !! @dimension n_samples = size(samples,3)
        !! @require m >= 0
        !! @require n >= 0
        !! @require n_samples >= 0
        !! @when n_samples == 0: All native shape and distribution validation is skipped.
        !! m=size(samples,1), n=size(samples,2), n_samples=size(samples,3); all may be zero.
        !! Successful nonempty draws require positive m and n. Present mean has shape
        !! (m,n), row matrix (m,m) and column matrix (n,n), with finite entries.
        !! Row/column matrices must be SPD; only their lower triangles are used.
        !! Absent mean is zero and absent row/column matrices are identities.
        !! On failure, completed matrices are retained and all matrices after the
        !! failing matrix are NaN. State retains all draws consumed before failure.
        implicit none
        type(rng_state), intent(inout) :: state !! explicit random-number generator state
        real(dp), intent(out) :: samples(:, :, :) !! shape (m,n,n_samples); matrices along the third axis
        real(dp), intent(in), optional :: mean(:, :) !! shape (m,n); finite mean; zero when absent
        real(dp), intent(in), optional :: rowcov(:, :) !! shape (m,m); finite SPD row covariance; identity when absent
        real(dp), intent(in), optional :: colcov(:, :) !! shape (n,n); finite SPD column covariance; identity when absent
        integer, intent(out), optional :: status !! matrix_normal_status_* result

        integer :: info
        integer :: k

        do k = 1, size(samples, 3)
            call matrix_normal_rvs(state, samples(:,:,k), mean, rowcov, colcov, info)
            if (info /= matrix_normal_status_success) then
                if (k < size(samples,3)) samples(:,:,k+1:) = quiet_nan(0.0_dp)
                if (present(status)) status = info
                return
            end if
        end do
        if (present(status)) status = matrix_normal_status_success
    end subroutine matrix_normal_rvs_array

    pure subroutine process_parameters(m, n, mean, rowcov, colcov, mu, rm, cm, status)
        !! Resolve mean and row/column matrices without testing positive-definiteness.
        !! m >= 1 and n >= 1; omitted inputs become zero mean or identities.
        !! Output shapes are checked before use; mismatches return invalid-shape status.
        integer, intent(in) :: m !! matrix row dimension
        integer, intent(in) :: n !! matrix column dimension
        real(dp), intent(in), optional :: mean(:, :) !! shape (m,n); finite mean; zero when absent
        real(dp), intent(in), optional :: rowcov(:, :) !! shape (m,m); finite SPD row covariance; identity when absent
        real(dp), intent(in), optional :: colcov(:, :) !! shape (n,n); finite SPD column covariance; identity when absent
        real(dp), intent(out) :: mu(:, :) !! shape (m,n); resolved finite mean
        real(dp), intent(out) :: rm(:, :) !! shape (m,m); resolved row covariance
        real(dp), intent(out) :: cm(:, :) !! shape (n,n); resolved column covariance
        integer, intent(out) :: status !! matrix_normal_status_* result

        if (m < 1 .or. n < 1 .or. size(mu,1) /= m .or. size(mu,2) /= n .or. &
                size(rm,1) /= m .or. size(rm,2) /= m .or. &
                size(cm,1) /= n .or. size(cm,2) /= n) then
            status = matrix_normal_status_invalid_shape
            return
        end if
        mu = 0.0_dp
        if (present(mean)) then
            if (size(mean,1) /= m .or. size(mean,2) /= n .or. .not. all(ieee_is_finite(mean))) then
                status = matrix_normal_status_invalid_shape
                return
            end if
            mu = mean
        end if
        call covariance_or_identity(rowcov, m, rm, status)
        if (status /= matrix_normal_status_success) return
        call covariance_or_identity(colcov, n, cm, status)
    end subroutine process_parameters

    pure subroutine covariance_or_identity(cov, n, out, status)
        !! Copy a finite covariance or construct an identity without SPD validation.
        !! n >= 0; shape (n,n) is required for out and any present cov.
        !! Zero extents are permitted by this helper; the caller enforces n >= 1.
        real(dp), intent(in), optional :: cov(:, :) !! shape (n,n); finite covariance; identity when absent
        integer, intent(in) :: n !! required dimension
        real(dp), intent(out) :: out(:, :) !! shape (n,n); resolved covariance
        integer, intent(out) :: status !! matrix_normal_status_* result

        integer :: i

        out = 0.0_dp
        if (size(out,1) /= n .or. size(out,2) /= n) then
            status = matrix_normal_status_invalid_shape
            return
        end if
        if (present(cov)) then
            if (size(cov,1) /= n .or. size(cov,2) /= n) then
                status = matrix_normal_status_invalid_shape
                return
            end if
            if (.not. all(ieee_is_finite(cov))) then
                status = matrix_normal_status_invalid_parameter
                return
            end if
            out = cov
        else
            do i = 1, n
                out(i,i) = 1.0_dp
            end do
        end if
        status = matrix_normal_status_success
    end subroutine covariance_or_identity

end module scifort_matrix_normal
