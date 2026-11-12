! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Dense matrix exponential, scipy.linalg.expm style.
!
! An independent implementation of the standard 13/13 Pade approximation
! with scaling and squaring (Higham, Functions of Matrices, 2008).
! This fixed-degree method is simpler than SciPy's variable-degree scheme.
! It uses SciFort's LU backend and requires no external BLAS or LAPACK.
module scifort_linalg_expm
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_dense, only : lu_factor, lu_solve
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    ! Higham's 1-norm backward-error threshold for [13/13] Pade.
    real(dp), parameter :: theta_13 = 5.371920351148152_dp
    real(dp), parameter :: log_two = log(2.0_dp)
    integer, parameter :: max_scalings = 1022
    real(dp), parameter :: coeff(0:13) = [ &
        64764752532480000.0_dp, 32382376266240000.0_dp, &
        7771770303897600.0_dp, 1187353796428800.0_dp, &
        129060195264000.0_dp, 10559470521600.0_dp, &
        670442572800.0_dp, 33522128640.0_dp, &
        1323241920.0_dp, 40840800.0_dp, &
        960960.0_dp, 16380.0_dp, 182.0_dp, 1.0_dp]

    public :: expm

    interface expm
        module procedure expm_real
        module procedure expm_complex
        module procedure expm_real_batch
        module procedure expm_complex_batch
    end interface expm

contains

    subroutine expm_real(a, x, status, scalings)
        real(dp), intent(in) :: a(:, :) !! square real matrix of order n
        real(dp), intent(out) :: x(:, :) !! output exp(A), shape (n,n), zero on failure
        integer, intent(out) :: status !! linalg_status_*; zero on success
        integer, intent(out), optional :: scalings !! number of squarings; zero on failure

        real(dp), allocatable :: b(:, :), b2(:, :), b4(:, :), b6(:, :)
        real(dp), allocatable :: eye(:, :), u(:, :), v(:, :), factors(:, :), rhs(:, :)
        integer, allocatable :: pivots(:)
        real(dp) :: norm_a, scale_factor
        integer :: n, i, s, code

        x = 0.0_dp
        status = linalg_status_invalid_shape
        if (present(scalings)) scalings = 0
        n = size(a, 1)
        if (size(a, 2) /= n .or. size(x, 1) /= n .or. size(x, 2) /= n) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(a))) return
        if (n == 0) then
            status = linalg_status_success
            return
        end if
        ! Diagonal and scalar inputs avoid avoidable cancellation during squaring.
        if (is_diagonal_real(a)) then
            do i = 1, n
                x(i, i) = exp(a(i, i))
            end do
            if (.not. all(ieee_is_finite(x))) then
                x = 0.0_dp
                return
            end if
            status = linalg_status_success
            return
        end if

        norm_a = maxval(sum(abs(a), dim=1))
        if (.not. ieee_is_finite(norm_a)) return
        call choose_scalings(norm_a, s, code)
        if (code /= linalg_status_success) return
        scale_factor = scale(1.0_dp, -s)

        allocate(b(n,n), b2(n,n), b4(n,n), b6(n,n), &
                 eye(n,n), u(n,n), v(n,n), factors(n,n), rhs(n,n), pivots(n))
        eye = 0.0_dp
        do i = 1, n
            eye(i, i) = 1.0_dp
        end do
        b = scale_factor*a
        b2 = matmul(b,b)
        b4 = matmul(b2,b2)
        b6 = matmul(b2,b4)
        u = matmul(b6, coeff(13)*b6 + coeff(11)*b4 + coeff(9)*b2)
        u = matmul(b, u + coeff(7)*b6 + coeff(5)*b4 + &
            coeff(3)*b2 + coeff(1)*eye)
        v = matmul(b6, coeff(12)*b6 + coeff(10)*b4 + coeff(8)*b2) + &
            coeff(6)*b6 + coeff(4)*b4 + coeff(2)*b2 + coeff(0)*eye
        if (.not. all(ieee_is_finite(u)) .or. .not. all(ieee_is_finite(v))) return
        rhs = v + u
        call lu_factor(v-u, factors, pivots, code)
        if (code /= linalg_status_success) return
        call lu_solve(factors, pivots, rhs, x, code)
        if (code /= linalg_status_success) then
            x = 0.0_dp
            return
        end if
        do i = 1, s
            x = matmul(x,x)
            if (.not. all(ieee_is_finite(x))) then
                x = 0.0_dp
                return
            end if
        end do
        status = linalg_status_success
        if (present(scalings)) scalings = s
    end subroutine expm_real

    subroutine expm_complex(a, x, status, scalings)
        complex(dp), intent(in) :: a(:, :) !! square complex matrix of order n
        complex(dp), intent(out) :: x(:, :) !! output exp(A), shape (n,n), zero on failure
        integer, intent(out) :: status !! linalg_status_*; zero on success
        integer, intent(out), optional :: scalings !! number of squarings; zero on failure

        complex(dp), allocatable :: b(:, :), b2(:, :), b4(:, :), b6(:, :)
        complex(dp), allocatable :: eye(:, :), u(:, :), v(:, :), factors(:, :), rhs(:, :)
        integer, allocatable :: pivots(:)
        real(dp) :: norm_a, scale_factor
        integer :: n, i, s, code

        x = cmplx(0.0_dp, 0.0_dp, dp)
        status = linalg_status_invalid_shape
        if (present(scalings)) scalings = 0
        n = size(a, 1)
        if (size(a, 2) /= n .or. size(x, 1) /= n .or. size(x, 2) /= n) return
        status = linalg_status_no_convergence
        if (.not. finite_complex(a)) return
        if (n == 0) then
            status = linalg_status_success
            return
        end if
        if (is_diagonal_complex(a)) then
            do i = 1, n
                x(i, i) = exp(a(i, i))
            end do
            if (.not. finite_complex(x)) then
                x = cmplx(0.0_dp, 0.0_dp, dp)
                return
            end if
            status = linalg_status_success
            return
        end if

        norm_a = maxval(sum(abs(a), dim=1))
        if (.not. ieee_is_finite(norm_a)) return
        call choose_scalings(norm_a, s, code)
        if (code /= linalg_status_success) return
        scale_factor = scale(1.0_dp, -s)

        allocate(b(n,n), b2(n,n), b4(n,n), b6(n,n), &
                 eye(n,n), u(n,n), v(n,n), factors(n,n), rhs(n,n), pivots(n))
        eye = cmplx(0.0_dp, 0.0_dp, dp)
        do i = 1, n
            eye(i, i) = cmplx(1.0_dp, 0.0_dp, dp)
        end do
        b = scale_factor*a
        b2 = matmul(b,b)
        b4 = matmul(b2,b2)
        b6 = matmul(b2,b4)
        u = matmul(b6, coeff(13)*b6 + coeff(11)*b4 + coeff(9)*b2)
        u = matmul(b, u + coeff(7)*b6 + coeff(5)*b4 + &
            coeff(3)*b2 + coeff(1)*eye)
        v = matmul(b6, coeff(12)*b6 + coeff(10)*b4 + coeff(8)*b2) + &
            coeff(6)*b6 + coeff(4)*b4 + coeff(2)*b2 + coeff(0)*eye
        if (.not. finite_complex(u) .or. .not. finite_complex(v)) return
        rhs = v + u
        call lu_factor(v-u, factors, pivots, code)
        if (code /= linalg_status_success) return
        call lu_solve(factors, pivots, rhs, x, code)
        if (code /= linalg_status_success) then
            x = cmplx(0.0_dp, 0.0_dp, dp)
            return
        end if
        do i = 1, s
            x = matmul(x,x)
            if (.not. finite_complex(x)) then
                x = cmplx(0.0_dp, 0.0_dp, dp)
                return
            end if
        end do
        status = linalg_status_success
        if (present(scalings)) scalings = s
    end subroutine expm_complex

    subroutine expm_real_batch(a, x, status)
        real(dp), intent(in) :: a(:, :, :) !! batch of square matrices, (batch,n,n)
        real(dp), intent(out) :: x(:, :, :) !! batched matrix exponentials, shape of A
        integer, intent(out) :: status !! zero on success, otherwise first error
        integer :: n, j, code

        x = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a, 2)
        if (size(a,3) /= n .or. size(x,1) /= size(a,1) .or. &
            size(x,2) /= n .or. size(x,3) /= n) return
        do j = 1, size(a,1)
            call expm_real(a(j,:,:), x(j,:,:), code)
            if (code /= linalg_status_success) then
                x = 0.0_dp
                status = code
                return
            end if
        end do
        status = linalg_status_success
    end subroutine expm_real_batch

    subroutine expm_complex_batch(a, x, status)
        complex(dp), intent(in) :: a(:, :, :) !! batch of square matrices, (batch,n,n)
        complex(dp), intent(out) :: x(:, :, :) !! batched matrix exponentials, shape of A
        integer, intent(out) :: status !! zero on success, otherwise first error
        integer :: n, j, code

        x = cmplx(0.0_dp, 0.0_dp, dp)
        status = linalg_status_invalid_shape
        n = size(a, 2)
        if (size(a,3) /= n .or. size(x,1) /= size(a,1) .or. &
            size(x,2) /= n .or. size(x,3) /= n) return
        do j = 1, size(a,1)
            call expm_complex(a(j,:,:), x(j,:,:), code)
            if (code /= linalg_status_success) then
                x = cmplx(0.0_dp, 0.0_dp, dp)
                status = code
                return
            end if
        end do
        status = linalg_status_success
    end subroutine expm_complex_batch

    pure subroutine choose_scalings(norm_a, count, status)
        real(dp), intent(in) :: norm_a !! finite matrix 1-norm
        integer, intent(out) :: count !! scaling exponent
        integer, intent(out) :: status !! zero if exponent is safe
        real(dp) :: estimated

        count = 0
        status = linalg_status_success
        if (norm_a <= theta_13) return
        estimated = (log(norm_a)-log(theta_13))/log_two
        if (.not. ieee_is_finite(estimated) .or. estimated >= real(max_scalings,dp)) then
            status = linalg_status_no_convergence
            return
        end if
        count = max(0, ceiling(estimated))
    end subroutine choose_scalings

    pure logical function is_diagonal_real(a)
        real(dp), intent(in) :: a(:, :) !! square matrix
        integer :: i, j
        is_diagonal_real = .false.
        do j = 1, size(a,2)
            do i = 1, size(a,1)
                if (i == j) cycle
                if (abs(a(i,j)) > 0.0_dp) return
            end do
        end do
        is_diagonal_real = .true.
    end function is_diagonal_real

    pure logical function is_diagonal_complex(a)
        complex(dp), intent(in) :: a(:, :) !! square matrix
        integer :: i, j
        is_diagonal_complex = .false.
        do j = 1, size(a,2)
            do i = 1, size(a,1)
                if (i == j) cycle
                if (abs(a(i,j)) > 0.0_dp) return
            end do
        end do
        is_diagonal_complex = .true.
    end function is_diagonal_complex

    pure logical function finite_complex(a)
        complex(dp), intent(in) :: a(:, :) !! input complex matrix
        finite_complex = all(ieee_is_finite(real(a,dp))) .and. &
            all(ieee_is_finite(aimag(a)))
    end function finite_complex

end module scifort_linalg_expm
