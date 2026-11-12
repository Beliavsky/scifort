! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Ordinary nonsymmetric matrix eigenproblems for modest dense matrices.
! Schur reduction reuses SciFort's MIT-licensed dense_complex_schur;
! left and right eigenvectors are recovered by triangular substitution.
module scifort_linalg_general_eigen
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_matrix_functions, only : dense_complex_schur, matrix_function_max_order
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    public :: eig, eigvals

    interface eig
        module procedure eig_real
        module procedure eig_complex
    end interface eig

    interface eigvals
        module procedure eigvals_real
        module procedure eigvals_complex
    end interface eigvals

contains

    subroutine eig_real(a, w, status, vr, vl)
        real(dp), intent(in) :: a(:, :) !! real square coefficient matrix (n,n), 0 <= n <= 40
        complex(dp), intent(out) :: w(:) !! n eigenvalues (order unspecified)
        integer, intent(out) :: status !! success, invalid shape, or numerical nonconvergence
        complex(dp), intent(out), optional :: vr(:, :) !! normalized right eigenvectors (n,n)
        complex(dp), intent(out), optional :: vl(:, :) !! normalized left eigenvectors (n,n)

        w = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(vr)) vr = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(vl)) vl = cmplx(0.0_dp, 0.0_dp, dp)
        status = linalg_status_invalid_shape
        if (.not. valid_shapes(size(a,1),size(a,2),size(w),vr,vl)) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(a))) return
        call eig_core(cmplx(a, kind=dp), w, status, vr, vl)
    end subroutine eig_real

    subroutine eig_complex(a, w, status, vr, vl)
        complex(dp), intent(in) :: a(:, :) !! complex square coefficient matrix (n,n), 0 <= n <= 40
        complex(dp), intent(out) :: w(:) !! n eigenvalues (order unspecified)
        integer, intent(out) :: status !! success, invalid shape, or numerical nonconvergence
        complex(dp), intent(out), optional :: vr(:, :) !! normalized right eigenvectors (n,n)
        complex(dp), intent(out), optional :: vl(:, :) !! normalized left eigenvectors (n,n)

        w = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(vr)) vr = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(vl)) vl = cmplx(0.0_dp, 0.0_dp, dp)
        status = linalg_status_invalid_shape
        if (.not. valid_shapes(size(a,1),size(a,2),size(w),vr,vl)) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(real(a,dp))) .or. &
            .not. all(ieee_is_finite(aimag(a)))) return
        call eig_core(a, w, status, vr, vl)
    end subroutine eig_complex

    subroutine eigvals_real(a, w, status)
        real(dp), intent(in) :: a(:, :) !! square real matrix, order 0..40
        complex(dp), intent(out) :: w(:) !! n eigenvalues, order unspecified
        integer, intent(out) :: status !! linalg_status_* outcome

        call eig_real(a, w, status)
    end subroutine eigvals_real

    subroutine eigvals_complex(a, w, status)
        complex(dp), intent(in) :: a(:, :) !! square complex matrix, order 0..40
        complex(dp), intent(out) :: w(:) !! n eigenvalues, order unspecified
        integer, intent(out) :: status !! linalg_status_* outcome

        call eig_complex(a, w, status)
    end subroutine eigvals_complex

    logical function valid_shapes(n, m, nw, vr, vl)
        integer, intent(in) :: n !! number of rows
        integer, intent(in) :: m !! number of columns
        integer, intent(in) :: nw !! eigenvalue output capacity
        complex(dp), intent(in), optional :: vr(:, :) !! optional right eigenvector output (shape only)
        complex(dp), intent(in), optional :: vl(:, :) !! optional left eigenvector output (shape only)

        valid_shapes = .false.
        if (n /= m .or. n > matrix_function_max_order .or. nw /= n) return
        if (present(vr)) then
            if (size(vr,1) /= n .or. size(vr,2) /= n) return
        end if
        if (present(vl)) then
            if (size(vl,1) /= n .or. size(vl,2) /= n) return
        end if
        valid_shapes = .true.
    end function valid_shapes

    subroutine eig_core(a, w, status, vr, vl)
        complex(dp), intent(in) :: a(:, :) !! finite general square matrix
        complex(dp), intent(out) :: w(:) !! eigenvalues, in Schur diagonal order
        integer, intent(out) :: status !! success or no convergence
        complex(dp), intent(out), optional :: vr(:, :) !! right eigenvectors
        complex(dp), intent(out), optional :: vl(:, :) !! left eigenvectors

        complex(dp), allocatable :: t(:, :), q(:, :), y(:), v(:)
        real(dp) :: magnitude, residual, tolerance
        integer :: n, i, j, code
        logical :: upper_triangular, lower_triangular

        status = linalg_status_no_convergence
        n = size(a,1)
        if (n == 0) then
            status = linalg_status_success
            return
        end if
        magnitude = max(maxval(abs(real(a,dp))),maxval(abs(aimag(a))))
        if (magnitude == 0.0_dp) then
            if (present(vr)) then
                do i = 1, n
                    vr(i,i) = cmplx(1.0_dp,0.0_dp,dp)
                end do
            end if
            if (present(vl)) then
                do i = 1, n
                    vl(i,i) = cmplx(1.0_dp,0.0_dp,dp)
                end do
            end if
            status = linalg_status_success
            return
        end if
        allocate(t(n,n),q(n,n),y(n),v(n))
        ! A triangular matrix already gives its exact diagonal eigenvalues.
        ! Reverse a lower triangular input to upper triangular form. Avoiding
        ! iterative QR here also protects highly nonnormal triangular cases.
        upper_triangular = .true.
        lower_triangular = .true.
        do j = 1, n-1
            if (any(abs(a(j+1:n,j)) > 0.0_dp)) upper_triangular = .false.
            if (any(abs(a(j,j+1:n)) > 0.0_dp)) lower_triangular = .false.
        end do
        q = cmplx(0.0_dp,0.0_dp,dp)
        if (upper_triangular) then
            t = a/magnitude
            do i = 1, n
                q(i,i) = cmplx(1.0_dp,0.0_dp,dp)
            end do
        else if (lower_triangular) then
            t = a(n:1:-1,n:1:-1)/magnitude
            do i = 1, n
                q(i,n+1-i) = cmplx(1.0_dp,0.0_dp,dp)
            end do
        else
            call dense_complex_schur(a/magnitude,t,q,code)
            if (code /= linalg_status_success) return
        end if
        do j = 1, n
            w(j) = t(j,j)*magnitude
        end do
        if (.not. all(ieee_is_finite(real(w,dp))) .or. &
            .not. all(ieee_is_finite(aimag(w)))) return
        if (.not. present(vr) .and. .not. present(vl)) then
            status = linalg_status_success
            return
        end if
        tolerance = 5.0e-9_dp
        do j = 1, n
            if (present(vr)) then
                call triangular_eigenvector(t,j,.false.,y)
                v = matmul(q,y)
                call normalize_eigenvector(v,code)
                if (code /= linalg_status_success) return
                residual = maxval(abs(matmul(t,y)-t(j,j)*y)) / &
                    max(1.0_dp,maxval(abs(t)))
                if (.not. ieee_is_finite(residual) .or. residual > tolerance) return
                vr(:,j) = v
            end if
            if (present(vl)) then
                call triangular_eigenvector(t,j,.true.,y)
                v = matmul(q,y)
                call normalize_eigenvector(v,code)
                if (code /= linalg_status_success) return
                residual = maxval(abs(matmul(transpose(conjg(t)),y)-conjg(t(j,j))*y)) / &
                    max(1.0_dp,maxval(abs(t)))
                if (.not. ieee_is_finite(residual) .or. residual > tolerance) return
                vl(:,j) = v
            end if
        end do
        status = linalg_status_success
    end subroutine eig_core

    subroutine triangular_eigenvector(t, j, left, y)
        complex(dp), intent(in) :: t(:, :) !! upper triangular complex Schur matrix
        integer, intent(in) :: j !! target diagonal index, 1..n
        logical, intent(in) :: left !! true: eigenvector of T^H
        complex(dp), intent(out) :: y(:) !! unnormalized Schur-coordinate eigenvector

        complex(dp) :: delta, subtotal
        real(dp) :: separation, magnitude
        integer :: i, n

        n = size(t,1)
        separation = 32.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(t)))
        y = cmplx(0.0_dp,0.0_dp,dp)
        y(j) = cmplx(1.0_dp,0.0_dp,dp)
        if (left) then
            do i = j+1, n
                subtotal = sum(conjg(t(j:i-1,i))*y(j:i-1))
                delta = conjg(t(i,i)-t(j,j))
                if (abs(delta) < separation) then
                    if (abs(delta) > 0.0_dp) then
                        delta = separation*(delta/abs(delta))
                    else
                        delta = cmplx(separation,0.0_dp,dp)
                    end if
                end if
                y(i) = -subtotal/delta
                magnitude = maxval(abs(y))
                if (magnitude > 1.0e50_dp) y = y/magnitude
            end do
        else
            do i = j-1, 1, -1
                subtotal = sum(t(i,i+1:j)*y(i+1:j))
                delta = t(i,i)-t(j,j)
                if (abs(delta) < separation) then
                    if (abs(delta) > 0.0_dp) then
                        delta = separation*(delta/abs(delta))
                    else
                        delta = cmplx(separation,0.0_dp,dp)
                    end if
                end if
                y(i) = -subtotal/delta
                magnitude = maxval(abs(y))
                if (magnitude > 1.0e50_dp) y = y/magnitude
            end do
        end if
        magnitude = maxval(abs(y))
        if (magnitude > 0.0_dp) y = y/magnitude
    end subroutine triangular_eigenvector

    subroutine normalize_eigenvector(v, status)
        complex(dp), intent(inout) :: v(:) !! vector to normalize to Euclidean norm one
        integer, intent(out) :: status !! success or numerical failure
        real(dp) :: magnitude, norm

        status = linalg_status_no_convergence
        magnitude = maxval(abs(v))
        if (.not. ieee_is_finite(magnitude) .or. magnitude <= 0.0_dp) return
        v = v/magnitude
        norm = sqrt(sum(real(conjg(v)*v,dp)))
        if (.not. ieee_is_finite(norm) .or. norm <= 0.0_dp) return
        v = v/norm
        status = linalg_status_success
    end subroutine normalize_eigenvector

end module scifort_linalg_general_eigen
