! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Dense matrix functions for modest-size real and complex matrices.
! The complex Schur/QR engine is adapted from SciFort's MIT-licensed
! src/sparse_linalg_complex_eigs.f90 (not from a third-party library).
!
module scifort_linalg_matrix_functions
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_singular, linalg_status_no_convergence
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    integer, parameter, public :: matrix_function_max_order = 40
    integer, parameter :: default_sign_iterations = 100
    real(dp), parameter :: default_tolerance = 1.0e-11_dp
    public :: sqrtm, sqrtm_real, signm
    public :: dense_inverse ! Reused by other small dense matrix-function implementations.
    public :: dense_complex_schur ! Reused by generic dense matrix functions.

    interface sqrtm
        module procedure sqrtm_real_input
        module procedure sqrtm_complex_input
    end interface
    interface signm
        module procedure signm_real_input
        module procedure signm_complex_input
    end interface

contains

    subroutine sqrtm_real_input(a, x, status, residual)
        real(dp), intent(in) :: a(:, :) !! real square input matrix (n,n)
        complex(dp), intent(out) :: x(:, :) !! principal complex square root (n,n)
        integer, intent(out) :: status !! linalg_status_* result
        real(dp), intent(out), optional :: residual !! relative Frobenius residual

        call sqrtm_complex_input(cmplx(a, kind=dp), x, status, residual, real_input=.true.)
    end subroutine sqrtm_real_input

    subroutine sqrtm_real(a, x, status, residual)
        real(dp), intent(in) :: a(:, :) !! real square input matrix (n,n)
        real(dp), intent(out) :: x(:, :) !! real-valued principal square root when it exists
        integer, intent(out) :: status !! success only when the complex root is numerically real
        real(dp), intent(out), optional :: residual !! relative Frobenius residual
        complex(dp), allocatable :: xc(:, :)
        real(dp) :: scale

        x = 0.0_dp
        allocate(xc(size(x, 1), size(x, 2)))
        call sqrtm_real_input(a, xc, status, residual)
        if (status /= linalg_status_success) return
        scale = max(1.0_dp, maxval(abs(xc)))
        if (maxval(abs(aimag(xc))) > 1.0e-9_dp*scale) then
            status = linalg_status_no_convergence
            return
        end if
        x = real(xc, dp)
    end subroutine sqrtm_real

    subroutine sqrtm_complex_input(a, x, status, residual, real_input)
        complex(dp), intent(in) :: a(:, :) !! complex square input matrix (n,n)
        complex(dp), intent(out) :: x(:, :) !! principal complex square root (n,n)
        integer, intent(out) :: status !! linalg_status_* result
        real(dp), intent(out), optional :: residual !! relative Frobenius residual
        logical, intent(in), optional :: real_input !! force +imaginary branch for negative real eigenvalues
        complex(dp), allocatable :: t(:, :), z(:, :), f(:, :), work(:, :)
        complex(dp) :: rhs, denominator
        integer :: n, i, j, k
        real(dp) :: anorm, root_norm, threshold, err
        logical :: snap_real

        x = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(residual)) residual = huge(1.0_dp)
        n = size(a, 1)
        if (n < 1 .or. n > matrix_function_max_order .or. size(a, 2) /= n .or. &
                size(x, 1) /= n .or. size(x, 2) /= n) then
            status = linalg_status_invalid_shape
            return
        end if
        if (.not. finite_matrix(a)) then
            status = linalg_status_no_convergence
            return
        end if
        if (all(a == cmplx(0.0_dp, 0.0_dp, dp))) then
            status = linalg_status_success
            if (present(residual)) residual = 0.0_dp
            return
        end if
        snap_real = .false.
        if (present(real_input)) snap_real = real_input
        anorm = fro_norm(a)
        allocate(t(n, n), z(n, n), f(n, n), work(n, n))
        call dense_complex_schur(a, t, z, status)
        if (status /= linalg_status_success) return

        f = cmplx(0.0_dp, 0.0_dp, dp)
        do i = 1, n
            if (snap_real .and. real(t(i, i), dp) < 0.0_dp .and. &
                    abs(aimag(t(i, i))) < 128.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(t(i, i)))) then
                f(i, i) = cmplx(0.0_dp, sqrt(-real(t(i, i), dp)), dp)
            else
                f(i, i) = sqrt(t(i, i))
            end if
        end do
        ! R_ij = (T_ij - sum_k R_ik R_kj)/(R_ii + R_jj).
        do j = 2, n
            do i = j-1, 1, -1
                rhs = t(i, j)
                do k = i+1, j-1
                    rhs = rhs-f(i, k)*f(k, j)
                end do
                denominator = f(i, i)+f(j, j)
                threshold = 16.0_dp*epsilon(1.0_dp)* &
                    max(1.0_dp, abs(f(i, i)), abs(f(j, j)))
                if (abs(denominator) <= threshold) then
                    status = linalg_status_no_convergence
                    return
                end if
                f(i, j) = rhs/denominator
            end do
        end do
        if (.not. finite_matrix(f)) then
            status = linalg_status_no_convergence
            return
        end if
        work = matmul(z, f)
        x = matmul(work, transpose(conjg(z)))
        root_norm = fro_norm(x)
        err = fro_norm(matmul(x, x)-a)/max(1.0_dp, anorm, root_norm*root_norm)
        if (present(residual)) residual = err
        if (.not. finite_matrix(x) .or. .not. ieee_is_finite(err) .or. err > 1.0e-8_dp) then
            x = cmplx(0.0_dp, 0.0_dp, dp)
            status = linalg_status_no_convergence
            return
        end if
        status = linalg_status_success
    end subroutine sqrtm_complex_input

    subroutine signm_real_input(a, x, status, tol, max_iter, iterations, residual)
        real(dp), intent(in) :: a(:, :) !! real square matrix (n,n)
        real(dp), intent(out) :: x(:, :) !! real matrix sign (n,n)
        integer, intent(out) :: status !! linalg_status_* result
        real(dp), intent(in), optional :: tol !! relative tolerance, default 1e-11
        integer, intent(in), optional :: max_iter !! maximum scaled-Newton iterations, default 100
        integer, intent(out), optional :: iterations !! Newton iterations performed
        real(dp), intent(out), optional :: residual !! relative involution residual
        complex(dp), allocatable :: xc(:, :)
        real(dp) :: scale

        x = 0.0_dp
        allocate(xc(size(x, 1), size(x, 2)))
        call signm_complex_input(cmplx(a, kind=dp), xc, status, tol, max_iter, iterations, residual)
        if (status /= linalg_status_success) return
        scale = max(1.0_dp, maxval(abs(xc)))
        if (maxval(abs(aimag(xc))) > 1.0e-8_dp*scale) then
            status = linalg_status_no_convergence
            return
        end if
        x = real(xc, dp)
    end subroutine signm_real_input

    subroutine signm_complex_input(a, x, status, tol, max_iter, iterations, residual)
        complex(dp), intent(in) :: a(:, :) !! square matrix without imaginary-axis eigenvalues
        complex(dp), intent(out) :: x(:, :) !! complex matrix sign (n,n)
        integer, intent(out) :: status !! linalg_status_* result
        real(dp), intent(in), optional :: tol !! relative tolerance, default 1e-11
        integer, intent(in), optional :: max_iter !! maximum scaled-Newton iterations, default 100
        integer, intent(out), optional :: iterations !! Newton iterations performed
        real(dp), intent(out), optional :: residual !! relative involution residual
        complex(dp), allocatable :: b(:, :), binv(:, :), next(:, :), identity(:, :)
        real(dp) :: tolerance, scale, invscale, mu, change, err, commutator, norma
        integer :: n, it, limit, code, i

        x = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(iterations)) iterations = 0
        if (present(residual)) residual = huge(1.0_dp)
        n = size(a, 1)
        if (n < 1 .or. n > matrix_function_max_order .or. size(a, 2) /= n .or. &
                size(x, 1) /= n .or. size(x, 2) /= n) then
            status = linalg_status_invalid_shape
            return
        end if
        tolerance = default_tolerance
        if (present(tol)) tolerance = tol
        limit = default_sign_iterations
        if (present(max_iter)) limit = max_iter
        if (limit < 1 .or. .not. ieee_is_finite(tolerance) .or. tolerance <= 0.0_dp) then
            status = linalg_status_invalid_shape
            return
        end if
        if (.not. finite_matrix(a)) then
            status = linalg_status_no_convergence
            return
        end if
        if (all(a == cmplx(0.0_dp, 0.0_dp, dp))) then
            status = linalg_status_success
            if (present(residual)) residual = 0.0_dp
            return
        end if
        norma = fro_norm(a)
        allocate(b(n, n), binv(n, n), next(n, n), identity(n, n))
        identity = cmplx(0.0_dp, 0.0_dp, dp)
        do i = 1, n
            identity(i, i) = cmplx(1.0_dp, 0.0_dp, dp)
        end do
        b = a/max(1.0_dp, norma)
        status = linalg_status_no_convergence
        do it = 1, limit
            call dense_inverse(b, binv, code)
            if (code /= linalg_status_success) return
            scale = fro_norm(b)
            invscale = fro_norm(binv)
            if (scale <= 0.0_dp .or. invscale <= 0.0_dp) return
            mu = sqrt(invscale/scale)
            if (.not. ieee_is_finite(mu) .or. mu <= 0.0_dp) return
            next = 0.5_dp*(mu*b+binv/mu)
            if (.not. finite_matrix(next)) return
            change = fro_norm(next-b)/max(1.0_dp, fro_norm(next))
            b = next
            if (present(iterations)) iterations = it
            err = fro_norm(matmul(b, b)-identity)/max(1.0_dp, real(n, dp))
            if (present(residual)) residual = err
            if (change <= tolerance .and. err <= max(1.0e-9_dp, 10.0_dp*tolerance)) then
                commutator = fro_norm(matmul(a, b)-matmul(b, a))/ &
                    max(1.0_dp, norma*fro_norm(b))
                if (.not. ieee_is_finite(commutator) .or. commutator > 1.0e-8_dp) return
                x = b
                status = linalg_status_success
                return
            end if
        end do
    end subroutine signm_complex_input

    subroutine dense_inverse(a, inverse, status)
        complex(dp), intent(in) :: a(:, :) !! complex square matrix (n,n)
        complex(dp), intent(out) :: inverse(:, :) !! matrix inverse on success
        integer, intent(out) :: status !! success or singular/invalid
        complex(dp), allocatable :: b(:, :), trow(:)
        complex(dp) :: factor
        real(dp) :: pivot_scale
        integer :: n, i, k, pivot

        n = size(a, 1)
        inverse = cmplx(0.0_dp, 0.0_dp, dp)
        if (n < 1 .or. size(a, 2) /= n .or. size(inverse, 1) /= n .or. &
                size(inverse, 2) /= n) then
            status = linalg_status_invalid_shape
            return
        end if
        if (.not. finite_matrix(a)) then
            status = linalg_status_no_convergence
            return
        end if
        allocate(b(n, n), trow(n))
        b = a
        do i = 1, n
            inverse(i, i) = cmplx(1.0_dp, 0.0_dp, dp)
        end do
        pivot_scale = maxval(abs(b))
        do k = 1, n
            pivot = k-1+maxloc(abs(b(k:n, k)), dim=1)
            if (abs(b(pivot, k)) < max(32.0_dp*tiny(1.0_dp), &
                    8.0_dp*epsilon(1.0_dp)*pivot_scale)) then
                status = linalg_status_singular
                return
            end if
            if (pivot /= k) then
                trow = b(k, :)
                b(k, :) = b(pivot, :)
                b(pivot, :) = trow
                trow = inverse(k, :)
                inverse(k, :) = inverse(pivot, :)
                inverse(pivot, :) = trow
            end if
            factor = b(k, k)
            b(k, :) = b(k, :)/factor
            inverse(k, :) = inverse(k, :)/factor
            do i = 1, n
                if (i == k) cycle
                factor = b(i, k)
                b(i, :) = b(i, :)-factor*b(k, :)
                inverse(i, :) = inverse(i, :)-factor*inverse(k, :)
            end do
        end do
        if (.not. finite_matrix(inverse)) then
            status = linalg_status_no_convergence
        else
            status = linalg_status_success
        end if
    end subroutine dense_inverse

    ! Unitary Householder QR and Hessenberg shifted QR, adapted from SciFort's
    ! scifort_sparse_linalg_complex_eigs (MIT).  Returns A = Z*T*Z^H.
    subroutine dense_complex_schur(a, t, z, status)
        complex(dp), intent(in) :: a(:, :) !! dense square matrix
        complex(dp), intent(out) :: t(:, :) !! upper triangular Schur factor
        complex(dp), intent(out) :: z(:, :) !! unitary Schur basis
        integer, intent(out) :: status !! zero on successful convergence
        complex(dp), allocatable :: u(:, :), r(:, :), block(:, :), temp(:, :), v(:)
        complex(dp) :: scalar, alpha, beta, gamma, delta, mean, disc, mu
        integer :: n, i, j, l, m, it
        real(dp) :: threshold, err

        n = size(a, 1)
        status = linalg_status_no_convergence
        allocate(u(n, n), r(n, n), block(n, n), temp(n, n), v(n))
        t = a
        z = cmplx(0.0_dp, 0.0_dp, dp)
        do i = 1, n
            z(i, i) = cmplx(1.0_dp, 0.0_dp, dp)
        end do
        do j = 1, n-2
            call reflector(t(j+1:n, j), v(j+1:n))
            if (vector_norm(v(j+1:n)) == 0.0_dp) cycle
            do l = j, n
                scalar = 2.0_dp*dot_product(v(j+1:n), t(j+1:n, l))
                t(j+1:n, l) = t(j+1:n, l)-v(j+1:n)*scalar
            end do
            do i = 1, n
                scalar = 2.0_dp*sum(t(i, j+1:n)*v(j+1:n))
                t(i, j+1:n) = t(i, j+1:n)-scalar*conjg(v(j+1:n))
                scalar = 2.0_dp*sum(z(i, j+1:n)*v(j+1:n))
                z(i, j+1:n) = z(i, j+1:n)-scalar*conjg(v(j+1:n))
            end do
            t(j+2:n, j) = cmplx(0.0_dp, 0.0_dp, dp)
        end do
        m = n
        it = 0
        do while (m > 1)
            threshold = 32.0_dp*epsilon(1.0_dp)* &
                (abs(t(m-1, m-1))+abs(t(m, m))+1.0_dp)
            if (abs(t(m, m-1)) <= threshold) then
                t(m, m-1) = cmplx(0.0_dp, 0.0_dp, dp)
                m = m-1
                cycle
            end if
            it = it+1
            if (it > max(600, 250*n)) return
            alpha = t(m-1, m-1)
            beta = t(m-1, m)
            gamma = t(m, m-1)
            delta = t(m, m)
            mean = 0.5_dp*(alpha+delta)
            disc = sqrt((0.5_dp*(alpha-delta))**2+beta*gamma)
            mu = mean+disc
            if (abs(mean-disc-delta) < abs(mu-delta)) mu = mean-disc
            block(:m, :m) = t(:m, :m)
            do i = 1, m
                block(i, i) = block(i, i)-mu
            end do
            call complex_qr(block(:m, :m), u(:m, :m), r(:m, :m))
            t(:m, :m) = matmul(r(:m, :m), u(:m, :m))
            do i = 1, m
                t(i, i) = t(i, i)+mu
            end do
            if (m < n) then
                temp(:m, m+1:n) = matmul(transpose(conjg(u(:m, :m))), t(:m, m+1:n))
                t(:m, m+1:n) = temp(:m, m+1:n)
            end if
            temp(:, 1:m) = matmul(z(:, 1:m), u(:m, :m))
            z(:, 1:m) = temp(:, 1:m)
            do j = 1, m-2
                t(j+2:m, j) = cmplx(0.0_dp, 0.0_dp, dp)
            end do
            if (.not. finite_matrix(t)) return
        end do
        if (.not. finite_matrix(z)) return
        err = fro_norm(matmul(matmul(z, t), transpose(conjg(z)))-a)/max(1.0_dp, fro_norm(a))
        if (err > 1.0e-9_dp .or. .not. ieee_is_finite(err)) return
        status = linalg_status_success
    end subroutine dense_complex_schur

    subroutine complex_qr(b, u, r)
        complex(dp), intent(in) :: b(:, :) !! square input matrix
        complex(dp), intent(out) :: u(:, :) !! unitary Q
        complex(dp), intent(out) :: r(:, :) !! upper-triangular R
        complex(dp) :: v(size(b, 1)), scalar
        integer :: n, j, col, row

        n = size(b, 1)
        u = cmplx(0.0_dp, 0.0_dp, dp)
        do j = 1, n
            u(j, j) = cmplx(1.0_dp, 0.0_dp, dp)
        end do
        r = b
        do j = 1, n-1
            call reflector(r(j:n, j), v(j:n))
            if (vector_norm(v(j:n)) == 0.0_dp) cycle
            do col = j, n
                scalar = 2.0_dp*dot_product(v(j:n), r(j:n, col))
                r(j:n, col) = r(j:n, col)-v(j:n)*scalar
            end do
            do row = 1, n
                scalar = 2.0_dp*sum(u(row, j:n)*v(j:n))
                u(row, j:n) = u(row, j:n)-scalar*conjg(v(j:n))
            end do
        end do
    end subroutine complex_qr

    subroutine reflector(x, v)
        complex(dp), intent(in) :: x(:) !! column segment
        complex(dp), intent(out) :: v(:) !! normalized Householder vector
        complex(dp) :: phase
        real(dp) :: magnitude, denominator

        v = cmplx(0.0_dp, 0.0_dp, dp)
        magnitude = vector_norm(x)
        if (magnitude <= tiny(magnitude)) return
        phase = cmplx(1.0_dp, 0.0_dp, dp)
        if (abs(x(1)) > 0.0_dp) phase = x(1)/abs(x(1))
        v = x
        v(1) = v(1)+phase*magnitude
        denominator = vector_norm(v)
        if (denominator > 0.0_dp) v = v/denominator
    end subroutine reflector

    pure real(dp) function fro_norm(a) result(value)
        complex(dp), intent(in) :: a(:, :) !! matrix
        value = sqrt(sum(abs(a)**2))
    end function fro_norm

    pure real(dp) function vector_norm(a) result(value)
        complex(dp), intent(in) :: a(:) !! vector
        value = sqrt(sum(abs(a)**2))
    end function vector_norm

    pure logical function finite_matrix(a)
        complex(dp), intent(in) :: a(:, :) !! matrix
        finite_matrix = all(ieee_is_finite(real(a, dp))) .and. &
            all(ieee_is_finite(aimag(a)))
    end function finite_matrix

end module scifort_linalg_matrix_functions
