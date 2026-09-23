! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Shared dense-matrix helpers for matrix-valued probability distributions.

module scifort_matrix_distribution_helpers
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_linalg, only : cholesky_lower, linalg_status_success
    use scifort_special_elementary, only : digamma, gammaln
    implicit none
    private

    integer, parameter, public :: matrix_status_success = 0
    integer, parameter, public :: matrix_status_invalid_shape = 1
    integer, parameter, public :: matrix_status_invalid_parameter = 2
    integer, parameter, public :: matrix_status_not_positive_definite = 3

    public :: invert_lower
    public :: log_multivariate_gamma
    public :: multivariate_digamma_sum
    public :: solve_lower_matrix
    public :: solve_spd_matrix
    public :: spd_factor

contains

    pure subroutine spd_factor(a, l, logdet, status)
        real(dp), intent(in) :: a(:, :) !! square SPD matrix; lower triangle defines symmetry
        real(dp), intent(out) :: l(:, :) !! lower Cholesky factor
        real(dp), intent(out) :: logdet !! logarithm of determinant
        integer, intent(out) :: status !! matrix_status_* result

        integer :: i
        integer :: info
        integer :: n

        l = 0.0_dp
        logdet = 0.0_dp
        n = size(a, 1)
        if (n < 1 .or. size(a, 2) /= n .or. size(l, 1) /= n .or. size(l, 2) /= n) then
            status = matrix_status_invalid_shape
            return
        end if
        if (.not. all(ieee_is_finite(a))) then
            status = matrix_status_invalid_parameter
            return
        end if
        call cholesky_lower(a, l, info)
        if (info /= linalg_status_success) then
            status = matrix_status_not_positive_definite
            return
        end if
        logdet = 0.0_dp
        do i = 1, n
            logdet = logdet + 2.0_dp * log(l(i, i))
        end do
        status = matrix_status_success
    end subroutine spd_factor

    pure subroutine solve_lower_matrix(l, b, x, status)
        real(dp), intent(in) :: l(:, :) !! nonsingular lower triangular matrix
        real(dp), intent(in) :: b(:, :) !! right-hand sides by column
        real(dp), intent(out) :: x(:, :) !! solution of L X = B
        integer, intent(out) :: status !! matrix_status_* result

        integer :: i
        integer :: j
        integer :: n
        real(dp) :: rhs

        x = 0.0_dp
        n = size(l, 1)
        if (n < 1 .or. size(l, 2) /= n .or. size(b, 1) /= n .or. &
                size(x, 1) /= n .or. size(x, 2) /= size(b, 2)) then
            status = matrix_status_invalid_shape
            return
        end if
        do j = 1, size(b, 2)
            do i = 1, n
                if (.not. ieee_is_finite(l(i, i)) .or. l(i, i) == 0.0_dp) then
                    status = matrix_status_not_positive_definite
                    return
                end if
                rhs = b(i, j)
                if (i > 1) rhs = rhs - dot_product(l(i, 1:i-1), x(1:i-1, j))
                x(i, j) = rhs / l(i, i)
            end do
        end do
        status = matrix_status_success
    end subroutine solve_lower_matrix

    pure subroutine solve_spd_matrix(l, b, x, status)
        real(dp), intent(in) :: l(:, :) !! lower Cholesky factor of SPD coefficient matrix
        real(dp), intent(in) :: b(:, :) !! right-hand sides by column
        real(dp), intent(out) :: x(:, :) !! solution of L L^T X = B
        integer, intent(out) :: status !! matrix_status_* result

        integer :: i
        integer :: j
        integer :: n
        real(dp), allocatable :: y(:, :)
        real(dp) :: rhs

        x = 0.0_dp
        n = size(l, 1)
        if (n < 1 .or. size(l, 2) /= n .or. size(b, 1) /= n .or. &
                size(x, 1) /= n .or. size(x, 2) /= size(b, 2)) then
            status = matrix_status_invalid_shape
            return
        end if
        allocate(y(n, size(b, 2)))
        call solve_lower_matrix(l, b, y, status)
        if (status /= matrix_status_success) return

        do j = 1, size(b, 2)
            do i = n, 1, -1
                rhs = y(i, j)
                if (i < n) rhs = rhs - dot_product(l(i+1:n, i), x(i+1:n, j))
                x(i, j) = rhs / l(i, i)
            end do
        end do
        status = matrix_status_success
    end subroutine solve_spd_matrix

    pure subroutine invert_lower(l, linv, status)
        real(dp), intent(in) :: l(:, :) !! nonsingular lower triangular matrix
        real(dp), intent(out) :: linv(:, :) !! inverse lower triangular matrix
        integer, intent(out) :: status !! matrix_status_* result

        integer :: i
        integer :: n
        real(dp), allocatable :: ident(:, :)

        n = size(l, 1)
        linv = 0.0_dp
        if (n < 1 .or. size(l, 2) /= n .or. size(linv, 1) /= n .or. size(linv, 2) /= n) then
            status = matrix_status_invalid_shape
            return
        end if
        allocate(ident(n, n))
        ident = 0.0_dp
        do i = 1, n
            ident(i, i) = 1.0_dp
        end do
        call solve_lower_matrix(l, ident, linv, status)
    end subroutine invert_lower

    pure function log_multivariate_gamma(a, p) result(y)
        real(dp), intent(in) :: a !! scalar argument, must exceed (p-1)/2
        integer, intent(in) :: p !! positive matrix dimension
        real(dp) :: y

        integer :: i

        if (p < 1 .or. a <= 0.5_dp * real(p - 1, dp) .or. .not. ieee_is_finite(a)) then
            y = huge(1.0_dp)
            return
        end if
        y = 0.25_dp * real(p * (p - 1), dp) * log(acos(-1.0_dp))
        do i = 1, p
            y = y + gammaln(a + 0.5_dp * real(1 - i, dp))
        end do
    end function log_multivariate_gamma

    pure function multivariate_digamma_sum(a, p) result(y)
        real(dp), intent(in) :: a !! base argument to multivariate gamma derivative
        integer, intent(in) :: p !! positive matrix dimension
        real(dp) :: y

        integer :: i

        y = 0.0_dp
        do i = 1, p
            y = y + digamma(a + 0.5_dp * real(1 - i, dp))
        end do
    end function multivariate_digamma_sum

end module scifort_matrix_distribution_helpers
