! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Shared dense random-matrix helpers. Haar orthogonal/unitary matrices are
! generated from Gaussian columns by modified Gram-Schmidt with positive
! real diagonal convention, equivalent in distribution to the QR construction
! described by Mezzadri (2007).

module scifort_random_matrix_helpers
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_normal, only : normal_ppf
    use scifort_random, only : rng_state, rng_uniform
    implicit none
    private

    integer, parameter, public :: random_matrix_status_success = 0
    integer, parameter, public :: random_matrix_status_invalid_shape = 1
    integer, parameter, public :: random_matrix_status_numerical_failure = 2

    public :: haar_orthogonal
    public :: haar_unitary
    public :: real_determinant

contains

    subroutine haar_orthogonal(state, q, status)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by the Gaussian draws
        real(dp), intent(out) :: q(:, :) !! square Haar-distributed orthogonal matrix
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid/numerical failure

        integer, parameter :: max_attempts = 4
        integer :: attempt
        integer :: i
        integer :: j
        integer :: n
        real(dp) :: coeff
        real(dp) :: normv
        real(dp), allocatable :: v(:)
        real(dp), allocatable :: z(:, :)

        n = size(q, 1)
        q = 0.0_dp
        if (size(q, 2) /= n) then
            if (present(status)) status = random_matrix_status_invalid_shape
            return
        end if
        if (n == 0) then
            if (present(status)) status = random_matrix_status_success
            return
        end if

        allocate(v(n), z(n,n))
        do attempt = 1, max_attempts
            do j = 1, n
                do i = 1, n
                    z(i,j) = normal_ppf(rng_uniform(state))
                end do
            end do
            q = 0.0_dp
            do j = 1, n
                v = z(:,j)
                do i = 1, j - 1
                    coeff = dot_product(q(:,i), v)
                    v = v - coeff * q(:,i)
                end do
                ! One re-orthogonalization pass improves loss of orthogonality
                ! without changing the Gaussian-QR distribution.
                do i = 1, j - 1
                    coeff = dot_product(q(:,i), v)
                    v = v - coeff * q(:,i)
                end do
                normv = sqrt(max(0.0_dp, dot_product(v, v)))
                if (.not. ieee_is_finite(normv) .or. normv <= sqrt(tiny(1.0_dp))) exit
                q(:,j) = v / normv
            end do
            if (j > n) then
                if (present(status)) status = random_matrix_status_success
                return
            end if
        end do

        q = 0.0_dp
        if (present(status)) status = random_matrix_status_numerical_failure
    end subroutine haar_orthogonal

    subroutine haar_unitary(state, q, status)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by the complex Gaussian draws
        complex(dp), intent(out) :: q(:, :) !! square Haar-distributed unitary matrix
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid/numerical failure

        integer, parameter :: max_attempts = 4
        integer :: attempt
        integer :: i
        integer :: j
        integer :: n
        complex(dp) :: coeff
        real(dp) :: normv
        real(dp), parameter :: inv_sqrt_two = 0.707106781186547524400844362104849039_dp
        complex(dp), allocatable :: v(:)
        complex(dp), allocatable :: z(:, :)

        n = size(q, 1)
        q = cmplx(0.0_dp, 0.0_dp, kind=dp)
        if (size(q, 2) /= n) then
            if (present(status)) status = random_matrix_status_invalid_shape
            return
        end if
        if (n == 0) then
            if (present(status)) status = random_matrix_status_success
            return
        end if

        allocate(v(n), z(n,n))
        do attempt = 1, max_attempts
            do j = 1, n
                do i = 1, n
                    z(i,j) = inv_sqrt_two * cmplx(normal_ppf(rng_uniform(state)), &
                        normal_ppf(rng_uniform(state)), kind=dp)
                end do
            end do
            q = cmplx(0.0_dp, 0.0_dp, kind=dp)
            do j = 1, n
                v = z(:,j)
                do i = 1, j - 1
                    coeff = dot_product(q(:,i), v)
                    v = v - coeff * q(:,i)
                end do
                do i = 1, j - 1
                    coeff = dot_product(q(:,i), v)
                    v = v - coeff * q(:,i)
                end do
                normv = sqrt(max(0.0_dp, real(dot_product(v, v), dp)))
                if (.not. ieee_is_finite(normv) .or. normv <= sqrt(tiny(1.0_dp))) exit
                q(:,j) = v / normv
            end do
            if (j > n) then
                if (present(status)) status = random_matrix_status_success
                return
            end if
        end do

        q = cmplx(0.0_dp, 0.0_dp, kind=dp)
        if (present(status)) status = random_matrix_status_numerical_failure
    end subroutine haar_unitary

    pure function real_determinant(a) result(det)
        real(dp), intent(in) :: a(:, :) !! square real matrix
        real(dp) :: det

        integer :: i
        integer :: j
        integer :: k
        integer :: n
        integer :: pivot
        real(dp), allocatable :: b(:, :)
        real(dp) :: factor
        real(dp) :: largest
        real(dp) :: temp_row(size(a,2))

        n = size(a, 1)
        if (size(a, 2) /= n) then
            det = 0.0_dp
            return
        end if
        if (n == 0) then
            det = 1.0_dp
            return
        end if

        allocate(b(n,n))
        b = a
        det = 1.0_dp
        do k = 1, n
            pivot = k
            largest = abs(b(k,k))
            do i = k + 1, n
                if (abs(b(i,k)) > largest) then
                    largest = abs(b(i,k))
                    pivot = i
                end if
            end do
            if (largest <= tiny(1.0_dp)) then
                det = 0.0_dp
                return
            end if
            if (pivot /= k) then
                temp_row = b(k,:)
                b(k,:) = b(pivot,:)
                b(pivot,:) = temp_row
                det = -det
            end if
            det = det * b(k,k)
            do i = k + 1, n
                factor = b(i,k) / b(k,k)
                do j = k + 1, n
                    b(i,j) = b(i,j) - factor * b(k,j)
                end do
            end do
        end do
    end function real_determinant

end module scifort_random_matrix_helpers
