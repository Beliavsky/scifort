! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Relative Frobenius condition number of the matrix exponential.
! The norm of the Frechet derivative is the largest singular value of its
! linear map E -> L_exp(A,E). Its adjoint is E -> L_exp(A^H,E). A fully
! reorthogonalized Lanczos iteration estimates the largest eigenvalue of
! L_exp(A,.)^* L_exp(A,.) without forming an n^2-by-n^2 Kronecker matrix.
module scifort_linalg_expm_condition
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence, symmetric_eigen_jacobi
    use scifort_linalg_expm_frechet, only : expm_frechet, expm_frechet_max_order
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private
    real(dp), parameter :: default_tol = 1.0e-9_dp
    integer, parameter :: default_max_iterations = 80
    public :: expm_cond

    interface expm_cond
        module procedure expm_cond_real
        module procedure expm_cond_complex
    end interface expm_cond

contains

    subroutine expm_cond_real(a, condition, status, tol, max_iterations, iterations)
        real(dp), intent(in) :: a(:, :) !! square real matrix, 1 <= n <= 48
        real(dp), intent(out) :: condition !! relative Frobenius condition number, zero on failure
        integer, intent(out) :: status !! linalg_status_* result
        real(dp), intent(in), optional :: tol !! relative Ritz residual tolerance, default 1e-9
        integer, intent(in), optional :: max_iterations !! Lanczos steps, default 80
        integer, intent(out), optional :: iterations !! actual Lanczos steps executed

        call condition_core(cmplx(a, kind=dp), .true., condition, status, &
            tol, max_iterations, iterations)
    end subroutine expm_cond_real

    subroutine expm_cond_complex(a, condition, status, tol, max_iterations, iterations)
        complex(dp), intent(in) :: a(:, :) !! square complex matrix, 1 <= n <= 48
        real(dp), intent(out) :: condition !! relative Frobenius condition number, zero on failure
        integer, intent(out) :: status !! linalg_status_* result
        real(dp), intent(in), optional :: tol !! relative Ritz residual tolerance, default 1e-9
        integer, intent(in), optional :: max_iterations !! Lanczos steps, default 80
        integer, intent(out), optional :: iterations !! actual Lanczos steps executed

        call condition_core(a, .false., condition, status, tol, max_iterations, iterations)
    end subroutine expm_cond_complex

    subroutine condition_core(a, real_input, condition, status, tol, max_iterations, iterations)
        complex(dp), intent(in) :: a(:, :) !! input square matrix
        logical, intent(in) :: real_input !! use real perturbations for the real API
        real(dp), intent(out) :: condition !! relative condition number
        integer, intent(out) :: status !! linalg_status_* result
        real(dp), intent(in), optional :: tol !! positive residual tolerance
        integer, intent(in), optional :: max_iterations !! positive iteration limit
        integer, intent(out), optional :: iterations !! Lanczos steps actually completed
        complex(dp), allocatable :: basis(:, :, :), shifted(:, :), exp_a(:, :), work(:, :), &
            derivative(:, :), image(:, :), w(:, :)
        real(dp), allocatable :: diagonal(:), offdiag(:), tridiagonal(:, :), &
            eigenvalues(:), eigenvectors(:, :)
        real(dp) :: tolerance, anorm, xnorm, scale_norm, bnorm, beta, lambda, &
            residual, seed_norm
        integer :: n, k, j, pass, limit, checked, eig_status, i
        complex(dp) :: projection, shift

        condition = 0.0_dp
        status = linalg_status_invalid_shape
        if (present(iterations)) iterations = 0
        n = size(a, 1)
        if (n < 1 .or. n > expm_frechet_max_order .or. size(a, 2) /= n) return
        if (.not. all(ieee_is_finite(real(a, dp))) .or. &
                .not. all(ieee_is_finite(aimag(a)))) return
        tolerance = default_tol
        if (present(tol)) tolerance = tol
        limit = default_max_iterations
        if (present(max_iterations)) limit = max_iterations
        if (.not. ieee_is_finite(tolerance) .or. tolerance <= 0.0_dp .or. &
                limit < 1) return
        limit = min(limit, n*n)
        anorm = frobenius(a)
        if (.not. ieee_is_finite(anorm)) return
        if (anorm == 0.0_dp) then
            status = linalg_status_success
            return
        end if

        allocate(basis(n,n,limit+1), shifted(n,n), exp_a(n,n), work(n,n), &
            derivative(n,n), image(n,n), w(n,n), diagonal(limit), &
            offdiag(limit), tridiagonal(limit,limit), &
            eigenvalues(limit), eigenvectors(limit,limit))
        ! A scalar shift changes exp(A) and L_exp(A,.) by the same
        ! factor exp(shift), which cancels from the relative condition
        ! number. This avoids overflow/underflow for large scalar offsets.
        shift = cmplx(maxval([(real(a(i,i),dp), i=1,n)]), &
            sum([(aimag(a(i,i)), i=1,n)])/real(n,dp), dp)
        shifted = a
        do i = 1, n
            shifted(i,i) = shifted(i,i)-shift
        end do
        ! For diagonal matrices the Frechet derivative acts entrywise
        ! through exponential divided differences; compute its exact norm.
        if (diagonal_matrix(shifted)) then
            call diagonal_condition(shifted, anorm, condition, status)
            return
        end if
        basis = cmplx(0.0_dp, 0.0_dp, dp)
        work = cmplx(0.0_dp, 0.0_dp, dp)
        call expm_frechet(shifted, work, exp_a, derivative, status)
        if (status /= linalg_status_success) return
        xnorm = frobenius(exp_a)
        if (.not. ieee_is_finite(xnorm) .or. xnorm <= tiny(1.0_dp)) then
            status = linalg_status_no_convergence
            return
        end if

        ! A reproducible non-symmetric seed, with all entries nonzero.
        ! Use a purely real seed when the input and perturbations are real.
        do j = 1, n
            do i = 1, n
                basis(i,j,1) = cmplx(sin(real(13*i+17*j,dp)), &
                    merge(cos(real(19*i+11*j,dp)), 0.0_dp, .not. real_input), dp)
            end do
        end do
        seed_norm = frobenius(basis(:,:,1))
        basis(:,:,1) = basis(:,:,1)/seed_norm
        diagonal = 0.0_dp
        offdiag = 0.0_dp
        tridiagonal = 0.0_dp
        status = linalg_status_no_convergence

        do k = 1, limit
            if (present(iterations)) iterations = k
            ! For exp: L(A,.)^*(E) = L(A^H,E) under the Frobenius product.
            call expm_frechet(shifted, basis(:,:,k), work, derivative, checked)
            if (checked /= linalg_status_success) exit
            call expm_frechet(transpose(conjg(shifted)), derivative, work, image, checked)
            if (checked /= linalg_status_success) exit
            bnorm = frobenius(image)
            if (.not. ieee_is_finite(bnorm)) exit
            projection = sum(conjg(basis(:,:,k))*image)
            diagonal(k) = real(projection, dp)
            w = image - diagonal(k)*basis(:,:,k)
            if (k > 1) w = w - offdiag(k-1)*basis(:,:,k-1)
            ! Reorthogonalize twice to keep the Lanczos Ritz values meaningful.
            do pass = 1, 2
                do j = 1, k
                    projection = sum(conjg(basis(:,:,j))*w)
                    w = w - projection*basis(:,:,j)
                end do
            end do
            beta = frobenius(w)
            if (.not. ieee_is_finite(beta)) exit
            offdiag(k) = beta

            tridiagonal(k,k) = diagonal(k)
            if (k > 1) then
                tridiagonal(k-1,k) = offdiag(k-1)
                tridiagonal(k,k-1) = offdiag(k-1)
            end if
            call symmetric_eigen_jacobi(tridiagonal(:k,:k), &
                eigenvalues(:k), eigenvectors(:k,:k), eig_status)
            if (eig_status /= linalg_status_success) exit
            lambda = max(0.0_dp, eigenvalues(k))
            residual = beta*abs(eigenvectors(k,k))
            if (lambda <= 0.0_dp) exit
            scale_norm = sqrt(lambda)
            condition = scale_norm*anorm/xnorm
            if (.not. ieee_is_finite(condition)) exit
            if (beta <= 32.0_dp*epsilon(1.0_dp)*max(bnorm,lambda) .or. &
                    residual <= tolerance*lambda) then
                status = linalg_status_success
                return
            end if
            if (k < limit) basis(:,:,k+1) = w/beta
        end do
        condition = 0.0_dp
        status = linalg_status_no_convergence
    end subroutine condition_core

    pure logical function diagonal_matrix(a) result(is_diagonal)
        complex(dp), intent(in) :: a(:, :) !! square matrix to inspect
        integer :: i, j
        is_diagonal = .false.
        do j = 1, size(a,2)
            do i = 1, size(a,1)
                if (i /= j .and. a(i,j) /= cmplx(0.0_dp, 0.0_dp, dp)) return
            end do
        end do
        is_diagonal = .true.
    end function diagonal_matrix

    subroutine diagonal_condition(a, anorm, condition, status)
        complex(dp), intent(in) :: a(:, :) !! shifted diagonal input matrix
        real(dp), intent(in) :: anorm !! Frobenius norm of unshifted matrix
        real(dp), intent(out) :: condition !! exact relative Frobenius condition number
        integer, intent(out) :: status !! success or no convergence
        complex(dp) :: lambda_i, lambda_j, delta, divided, term
        real(dp) :: norm_exp, max_divided
        integer :: i, j, k, n
        n = size(a,1)
        norm_exp = 0.0_dp
        max_divided = 0.0_dp
        do i = 1, n
            norm_exp = norm_exp + abs(exp(a(i,i)))**2
        end do
        norm_exp = sqrt(norm_exp)
        if (norm_exp <= tiny(1.0_dp)) then
            status = linalg_status_no_convergence
            return
        end if
        do j = 1, n
            lambda_j = a(j,j)
            do i = 1, n
                lambda_i = a(i,i)
                delta = lambda_i-lambda_j
                if (abs(delta) < 1.0e-2_dp) then
                    ! exp(lambda_j)*expm1(delta)/delta, evaluated as a
                    ! short series to retain accuracy for close eigenvalues.
                    divided = cmplx(1.0_dp, 0.0_dp, dp)
                    term = divided
                    do k = 1, 14
                        term = term*delta/real(k+1,dp)
                        divided = divided+term
                    end do
                    divided = exp(lambda_j)*divided
                else
                    divided = (exp(lambda_i)-exp(lambda_j))/delta
                end if
                max_divided = max(max_divided,abs(divided))
            end do
        end do
        condition = anorm*max_divided/norm_exp
        status = linalg_status_success
        if (.not. ieee_is_finite(condition)) then
            condition = 0.0_dp
            status = linalg_status_no_convergence
        end if
    end subroutine diagonal_condition

    pure real(dp) function frobenius(a) result(value)
        complex(dp), intent(in) :: a(:, :) !! matrix whose Frobenius norm is returned
        value = sqrt(sum(abs(a)**2))
    end function frobenius
end module scifort_linalg_expm_condition
