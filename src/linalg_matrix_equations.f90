! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Dense matrix equations: original pivoted Kronecker elimination for small
! problems, and Schur/Bartels-Stewart triangular recurrences for larger ones.
! The latter use SciFort's own complex Schur decomposition, not LAPACK.
!
module scifort_linalg_matrix_equations
    use scifort_kinds, only : dp
    use scifort_linalg_matrix_functions, only : dense_complex_schur, matrix_function_max_order
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_singular, linalg_status_no_convergence
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    integer, parameter, public :: matrix_equation_max_unknowns = 400
    integer, parameter, public :: matrix_equation_schur_max_order = matrix_function_max_order

    public :: solve_sylvester
    public :: solve_continuous_lyapunov
    public :: solve_discrete_lyapunov

    interface solve_sylvester
        module procedure solve_sylvester_real
        module procedure solve_sylvester_complex
    end interface

    interface solve_continuous_lyapunov
        module procedure solve_continuous_lyapunov_real
        module procedure solve_continuous_lyapunov_complex
    end interface

    interface solve_discrete_lyapunov
        module procedure solve_discrete_lyapunov_real
        module procedure solve_discrete_lyapunov_complex
    end interface

contains

    subroutine solve_sylvester_complex(a, b, q, x, status, method)
        complex(dp), intent(in) :: a(:, :) !! leading square matrix (m,m)
        complex(dp), intent(in) :: b(:, :) !! trailing square matrix (n,n)
        complex(dp), intent(in) :: q(:, :) !! right-hand side (m,n)
        complex(dp), intent(out) :: x(:, :) !! solution satisfying A*X + X*B = Q
        integer, intent(out) :: status !! linalg_status_* result
        character(*), intent(in), optional :: method !! 'auto' (default), 'direct', or 'schur'
        character(len=8) :: chosen

        call choose_method(method, size(a,1), size(b,1), chosen, status)
        x = cmplx(0.0_dp, 0.0_dp, dp)
        if (status /= linalg_status_success) return
        select case (chosen)
        case ('direct')
            call solve_sylvester_direct_complex(a, b, q, x, status)
        case ('schur')
            call solve_sylvester_schur_complex(a, b, q, x, status)
        end select
    end subroutine solve_sylvester_complex

    subroutine solve_discrete_lyapunov_complex(a, q, x, status, method)
        complex(dp), intent(in) :: a(:, :) !! square dynamics matrix (n,n)
        complex(dp), intent(in) :: q(:, :) !! right-hand side (n,n)
        complex(dp), intent(out) :: x(:, :) !! solution X - A*X*A^H = Q
        integer, intent(out) :: status !! linalg_status_* result
        character(*), intent(in), optional :: method !! 'auto' (default), 'direct', or 'schur'
        character(len=8) :: chosen

        call choose_method(method, size(a,1), size(a,1), chosen, status)
        x = cmplx(0.0_dp, 0.0_dp, dp)
        if (status /= linalg_status_success) return
        select case (chosen)
        case ('direct')
            call solve_discrete_lyapunov_direct_complex(a, q, x, status)
        case ('schur')
            call solve_lyapunov_schur_complex(a, q, x, status, .true.)
        end select
    end subroutine solve_discrete_lyapunov_complex

    subroutine choose_method(method, m, n, chosen, status)
        character(*), intent(in), optional :: method !! optional method name
        integer, intent(in) :: m !! matrix A order
        integer, intent(in) :: n !! matrix B order
        character(len=8), intent(out) :: chosen !! selected algorithm
        integer, intent(out) :: status !! invalid method or unsupported order

        chosen = 'auto'
        status = linalg_status_success
        if (present(method)) chosen = method
        if (m < 1 .or. n < 1) then
            status = linalg_status_invalid_shape
            return
        end if
        select case (trim(chosen))
        case ('auto')
            if (m <= matrix_equation_max_unknowns/n) then
                chosen = 'direct'
            else
                chosen = 'schur'
            end if
        case ('direct', 'schur')
            continue
        case default
            status = linalg_status_invalid_shape
            return
        end select
        if (m < 1 .or. n < 1) then
            status = linalg_status_invalid_shape
        else if (chosen == 'direct') then
            if (m > matrix_equation_max_unknowns/n) status = linalg_status_invalid_shape
        else
            if (m > matrix_equation_schur_max_order .or. n > matrix_equation_schur_max_order) &
                status = linalg_status_invalid_shape
        end if
    end subroutine choose_method

    subroutine solve_sylvester_direct_complex(a, b, q, x, status)
        complex(dp), intent(in) :: a(:, :) !! leading square matrix (m,m)
        complex(dp), intent(in) :: b(:, :) !! trailing square matrix (n,n)
        complex(dp), intent(in) :: q(:, :) !! right-hand side (m,n)
        complex(dp), intent(out) :: x(:, :) !! result satisfying A*X + X*B = Q
        integer, intent(out) :: status !! success, invalid shape, singular, or nonconvergence
        complex(dp), allocatable :: k(:, :)
        complex(dp), allocatable :: rhs(:)
        complex(dp), allocatable :: sol(:)
        integer :: m, n, i, j, p, r, row, col, nvar

        x = cmplx(0.0_dp, 0.0_dp, kind=dp)
        m = size(a, 1)
        n = size(b, 1)
        if (.not. valid_sylvester_shapes(m, n, size(a, 2), size(b, 2), &
                size(q, 1), size(q, 2), size(x, 1), size(x, 2))) then
            status = linalg_status_invalid_shape
            return
        end if
        if (.not. (finite_matrix(a) .and. finite_matrix(b) .and. finite_matrix(q))) then
            status = linalg_status_no_convergence
            return
        end if
        nvar = m * n
        allocate(k(nvar, nvar), rhs(nvar), sol(nvar))
        k = cmplx(0.0_dp, 0.0_dp, kind=dp)
        do j = 1, n
            do i = 1, m
                row = i + (j - 1) * m
                rhs(row) = q(i, j)
                do p = 1, m
                    col = p + (j - 1) * m
                    k(row, col) = k(row, col) + a(i, p)
                end do
                do r = 1, n
                    col = i + (r - 1) * m
                    k(row, col) = k(row, col) + b(r, j)
                end do
            end do
        end do
        call solve_pivoted_complex(k, rhs, sol, status)
        if (status == linalg_status_success) x = reshape(sol, [m, n])
    end subroutine solve_sylvester_direct_complex

    subroutine solve_continuous_lyapunov_complex(a, q, x, status, method)
        complex(dp), intent(in) :: a(:, :) !! square dynamics matrix (n,n)
        complex(dp), intent(in) :: q(:, :) !! right-hand side (n,n)
        complex(dp), intent(out) :: x(:, :) !! solution A*X + X*A^H = Q
        integer, intent(out) :: status !! success, invalid shape, singular, or nonconvergence
        character(*), intent(in), optional :: method !! 'auto' (default), 'direct', or 'schur'
        character(len=8) :: chosen

        call choose_method(method, size(a,1), size(a,1), chosen, status)
        x = cmplx(0.0_dp, 0.0_dp, dp)
        if (status /= linalg_status_success) return
        select case (chosen)
        case ('direct')
            call solve_sylvester_direct_complex(a, transpose(conjg(a)), q, x, status)
        case ('schur')
            call solve_lyapunov_schur_complex(a, q, x, status, .false.)
        end select
    end subroutine solve_continuous_lyapunov_complex

    subroutine solve_discrete_lyapunov_direct_complex(a, q, x, status)
        complex(dp), intent(in) :: a(:, :) !! square dynamics matrix (n,n)
        complex(dp), intent(in) :: q(:, :) !! right-hand side (n,n)
        complex(dp), intent(out) :: x(:, :) !! solution X - A*X*A^H = Q
        integer, intent(out) :: status !! success, invalid shape, singular, or nonconvergence
        complex(dp), allocatable :: k(:, :)
        complex(dp), allocatable :: rhs(:)
        complex(dp), allocatable :: sol(:)
        integer :: n, i, j, p, r, row, col, nvar

        x = cmplx(0.0_dp, 0.0_dp, kind=dp)
        n = size(a, 1)
        if (.not. valid_sylvester_shapes(n, n, size(a, 2), size(a, 2), &
                size(q, 1), size(q, 2), size(x, 1), size(x, 2))) then
            status = linalg_status_invalid_shape
            return
        end if
        if (.not. (finite_matrix(a) .and. finite_matrix(q))) then
            status = linalg_status_no_convergence
            return
        end if
        nvar = n * n
        allocate(k(nvar, nvar), rhs(nvar), sol(nvar))
        k = cmplx(0.0_dp, 0.0_dp, kind=dp)
        do j = 1, n
            do i = 1, n
                row = i + (j - 1) * n
                rhs(row) = q(i, j)
                k(row, row) = cmplx(1.0_dp, 0.0_dp, kind=dp)
                do r = 1, n
                    do p = 1, n
                        col = p + (r - 1) * n
                        k(row, col) = k(row, col) - a(i, p) * conjg(a(j, r))
                    end do
                end do
            end do
        end do
        call solve_pivoted_complex(k, rhs, sol, status)
        if (status == linalg_status_success) x = reshape(sol, [n, n])
    end subroutine solve_discrete_lyapunov_direct_complex

    ! Classical Bartels-Stewart transformation for the generalized Sylvester
    ! equation. T and S are triangular, and the scalar recurrences require
    ! only O(m*n*(m+n)) workspace operations after the two Schur factorizations.
    subroutine solve_sylvester_schur_complex(a, b, q, x, status)
        complex(dp), intent(in) :: a(:, :) !! square matrix of order m
        complex(dp), intent(in) :: b(:, :) !! square matrix of order n
        complex(dp), intent(in) :: q(:, :) !! right-hand side (m,n)
        complex(dp), intent(out) :: x(:, :) !! solution A*X + X*B = Q
        integer, intent(out) :: status !! linalg_status_* result
        complex(dp), allocatable :: t(:, :), u(:, :), s(:, :), v(:, :)
        complex(dp), allocatable :: c(:, :), y(:, :)
        complex(dp) :: rhs, denominator
        integer :: m, n, i, j, k
        real(dp) :: backward_error

        x = cmplx(0.0_dp, 0.0_dp, dp)
        m = size(a, 1)
        n = size(b, 1)
        if (m < 1 .or. n < 1 .or. m > matrix_equation_schur_max_order .or. &
                n > matrix_equation_schur_max_order .or. size(a,2) /= m .or. &
                size(b,2) /= n .or. size(q,1) /= m .or. size(q,2) /= n .or. &
                size(x,1) /= m .or. size(x,2) /= n) then
            status = linalg_status_invalid_shape
            return
        end if
        if (.not. (finite_matrix(a) .and. finite_matrix(b) .and. finite_matrix(q))) then
            status = linalg_status_no_convergence
            return
        end if
        allocate(t(m,m), u(m,m), s(n,n), v(n,n), c(m,n), y(m,n))
        call dense_complex_schur(a, t, u, status)
        if (status /= linalg_status_success) return
        call dense_complex_schur(b, s, v, status)
        if (status /= linalg_status_success) return
        c = matmul(matmul(transpose(conjg(u)), q), v)
        y = cmplx(0.0_dp, 0.0_dp, dp)
        do j = 1, n
            do i = m, 1, -1
                denominator = t(i,i) + s(j,j)
                if (singular_coefficient(denominator, t(i,i), s(j,j))) then
                    status = linalg_status_singular
                    return
                end if
                rhs = c(i,j)
                do k = i+1, m
                    rhs = rhs-t(i,k)*y(k,j)
                end do
                do k = 1, j-1
                    rhs = rhs-y(i,k)*s(k,j)
                end do
                y(i,j) = rhs/denominator
            end do
        end do
        if (.not. finite_matrix(y)) then
            status = linalg_status_no_convergence
            return
        end if
        x = matmul(matmul(u, y), transpose(conjg(v)))
        if (.not. finite_matrix(x)) then
            x = cmplx(0.0_dp, 0.0_dp, dp)
            status = linalg_status_no_convergence
            return
        end if
        backward_error = fro_norm(matmul(a,x)+matmul(x,b)-q) / &
            max(tiny(1.0_dp), (fro_norm(a)+fro_norm(b))*fro_norm(x)+fro_norm(q))
        if (.not. ieee_is_finite(backward_error) .or. backward_error > 1.0e-9_dp) then
            x = cmplx(0.0_dp, 0.0_dp, dp)
            status = linalg_status_no_convergence
            return
        end if
        status = linalg_status_success
    end subroutine solve_sylvester_schur_complex

    ! A single Schur factorization suffices for both continuous and discrete
    ! Lyapunov equations. The upper-triangular recurrences visit columns and
    ! rows backwards; entries needed by each step have already been computed.
    subroutine solve_lyapunov_schur_complex(a, q, x, status, discrete)
        complex(dp), intent(in) :: a(:, :) !! square matrix of order n
        complex(dp), intent(in) :: q(:, :) !! right-hand side (n,n)
        complex(dp), intent(out) :: x(:, :) !! Lyapunov solution (n,n)
        integer, intent(out) :: status !! linalg_status_* result
        logical, intent(in) :: discrete !! true: X-A*X*A^H=Q; false: A*X+X*A^H=Q
        complex(dp), allocatable :: t(:, :), u(:, :), c(:, :), y(:, :)
        complex(dp) :: rhs, denominator, self_term
        integer :: n, i, j, p, r
        real(dp) :: backward_error

        x = cmplx(0.0_dp, 0.0_dp, dp)
        n = size(a, 1)
        if (n < 1 .or. n > matrix_equation_schur_max_order .or. size(a,2) /= n .or. &
                size(q,1) /= n .or. size(q,2) /= n .or. &
                size(x,1) /= n .or. size(x,2) /= n) then
            status = linalg_status_invalid_shape
            return
        end if
        if (.not. (finite_matrix(a) .and. finite_matrix(q))) then
            status = linalg_status_no_convergence
            return
        end if
        allocate(t(n,n), u(n,n), c(n,n), y(n,n))
        call dense_complex_schur(a, t, u, status)
        if (status /= linalg_status_success) return
        c = matmul(matmul(transpose(conjg(u)),q),u)
        y = cmplx(0.0_dp, 0.0_dp, dp)
        do j = n, 1, -1
            do i = n, 1, -1
                self_term = t(i,i)*conjg(t(j,j))
                if (discrete) then
                    denominator = cmplx(1.0_dp,0.0_dp,dp)-self_term
                    if (singular_coefficient(denominator, cmplx(1.0_dp,0.0_dp,dp), -self_term)) then
                        status = linalg_status_singular
                        return
                    end if
                    rhs = c(i,j)
                    do r = j, n
                        do p = i, n
                            if (p == i .and. r == j) cycle
                            rhs = rhs+t(i,p)*y(p,r)*conjg(t(j,r))
                        end do
                    end do
                else
                    denominator = t(i,i)+conjg(t(j,j))
                    if (singular_coefficient(denominator, t(i,i), conjg(t(j,j)))) then
                        status = linalg_status_singular
                        return
                    end if
                    rhs = c(i,j)
                    do p = i+1, n
                        rhs = rhs-t(i,p)*y(p,j)
                    end do
                    do r = j+1, n
                        rhs = rhs-y(i,r)*conjg(t(j,r))
                    end do
                end if
                y(i,j) = rhs/denominator
            end do
        end do
        if (.not. finite_matrix(y)) then
            status = linalg_status_no_convergence
            return
        end if
        x = matmul(matmul(u,y),transpose(conjg(u)))
        if (.not. finite_matrix(x)) then
            x = cmplx(0.0_dp, 0.0_dp, dp)
            status = linalg_status_no_convergence
            return
        end if
        if (discrete) then
            backward_error = fro_norm(x-matmul(matmul(a,x),transpose(conjg(a)))-q) / &
                max(tiny(1.0_dp), (1.0_dp+fro_norm(a)**2)*fro_norm(x)+fro_norm(q))
        else
            backward_error = fro_norm(matmul(a,x)+matmul(x,transpose(conjg(a)))-q) / &
                max(tiny(1.0_dp), 2.0_dp*fro_norm(a)*fro_norm(x)+fro_norm(q))
        end if
        if (.not. ieee_is_finite(backward_error) .or. backward_error > 1.0e-9_dp) then
            x = cmplx(0.0_dp, 0.0_dp, dp)
            status = linalg_status_no_convergence
            return
        end if
        status = linalg_status_success
    end subroutine solve_lyapunov_schur_complex

    pure logical function singular_coefficient(denom, term_a, term_b) result(singular)
        complex(dp), intent(in) :: denom !! triangular scalar denominator
        complex(dp), intent(in) :: term_a !! first summand
        complex(dp), intent(in) :: term_b !! second summand
        real(dp) :: tolerance

        tolerance = 32.0_dp*epsilon(1.0_dp)* &
            max(abs(term_a), abs(term_b), tiny(1.0_dp))
        singular = abs(denom) <= tolerance
    end function singular_coefficient

    subroutine solve_sylvester_real(a, b, q, x, status, method)
        real(dp), intent(in) :: a(:, :) !! leading square matrix (m,m)
        real(dp), intent(in) :: b(:, :) !! trailing square matrix (n,n)
        real(dp), intent(in) :: q(:, :) !! real right-hand side (m,n)
        real(dp), intent(out) :: x(:, :) !! real solution (m,n)
        integer, intent(out) :: status !! success, invalid shape, singular, or nonconvergence
        character(*), intent(in), optional :: method !! 'auto', 'direct', or 'schur'
        complex(dp), allocatable :: xc(:, :)

        allocate(xc(size(x, 1), size(x, 2)))
        call solve_sylvester_complex(cmplx(a, kind=dp), cmplx(b, kind=dp), &
            cmplx(q, kind=dp), xc, status, method=method)
        x = real(xc, kind=dp)
    end subroutine solve_sylvester_real

    subroutine solve_continuous_lyapunov_real(a, q, x, status, method)
        real(dp), intent(in) :: a(:, :) !! real square dynamics matrix (n,n)
        real(dp), intent(in) :: q(:, :) !! real right-hand side (n,n)
        real(dp), intent(out) :: x(:, :) !! solution A*X + X*A^T = Q
        integer, intent(out) :: status !! success, invalid shape, singular, or nonconvergence
        character(*), intent(in), optional :: method !! 'auto', 'direct', or 'schur'
        complex(dp), allocatable :: xc(:, :)

        allocate(xc(size(x, 1), size(x, 2)))
        call solve_continuous_lyapunov_complex(cmplx(a, kind=dp), &
            cmplx(q, kind=dp), xc, status, method=method)
        x = real(xc, kind=dp)
    end subroutine solve_continuous_lyapunov_real

    subroutine solve_discrete_lyapunov_real(a, q, x, status, method)
        real(dp), intent(in) :: a(:, :) !! real square dynamics matrix (n,n)
        real(dp), intent(in) :: q(:, :) !! real right-hand side (n,n)
        real(dp), intent(out) :: x(:, :) !! solution X - A*X*A^T = Q
        integer, intent(out) :: status !! success, invalid shape, singular, or nonconvergence
        character(*), intent(in), optional :: method !! 'auto', 'direct', or 'schur'
        complex(dp), allocatable :: xc(:, :)

        allocate(xc(size(x, 1), size(x, 2)))
        call solve_discrete_lyapunov_complex(cmplx(a, kind=dp), &
            cmplx(q, kind=dp), xc, status, method=method)
        x = real(xc, kind=dp)
    end subroutine solve_discrete_lyapunov_real

    pure logical function valid_sylvester_shapes(m, n, am, bn, qm, qn, xm, xn)
        integer, intent(in) :: m !! number of rows of X
        integer, intent(in) :: n !! number of columns of X
        integer, intent(in) :: am !! number of columns of A
        integer, intent(in) :: bn !! number of columns of B
        integer, intent(in) :: qm !! number of rows of Q
        integer, intent(in) :: qn !! number of columns of Q
        integer, intent(in) :: xm !! number of rows of X
        integer, intent(in) :: xn !! number of columns of X

        valid_sylvester_shapes = .false.
        if (m <= 0 .or. n <= 0) return
        if (m > matrix_equation_max_unknowns / n) return
        if (am /= m .or. bn /= n .or. qm /= m .or. qn /= n) return
        if (xm /= m .or. xn /= n) return
        valid_sylvester_shapes = .true.
    end function valid_sylvester_shapes

    pure real(dp) function fro_norm(a) result(value)
        complex(dp), intent(in) :: a(:, :) !! input matrix
        value = sqrt(sum(abs(a)**2))
    end function fro_norm

    pure logical function finite_matrix(a)
        complex(dp), intent(in) :: a(:, :) !! complex input matrix
        finite_matrix = all(ieee_is_finite(real(a, dp))) .and. &
            all(ieee_is_finite(aimag(a)))
    end function finite_matrix

    subroutine solve_pivoted_complex(a, b, x, status)
        complex(dp), intent(inout) :: a(:, :) !! nonsingular coefficient matrix, overwritten by elimination
        complex(dp), intent(inout) :: b(:) !! right-hand side, overwritten by elimination
        complex(dp), intent(out) :: x(:) !! solution vector
        integer, intent(out) :: status !! zero if successful; singular/nonfinite otherwise
        complex(dp) :: factor, temp
        real(dp) :: scale, pivot_tol
        integer :: n, i, j, k, pivot

        x = cmplx(0.0_dp, 0.0_dp, kind=dp)
        n = size(b)
        if (n < 1 .or. size(a, 1) /= n .or. size(a, 2) /= n .or. size(x) /= n) then
            status = linalg_status_invalid_shape
            return
        end if
        if (.not. finite_matrix(a(:, :)) .or. .not. all(ieee_is_finite(real(b, dp))) .or. &
                .not. all(ieee_is_finite(aimag(b)))) then
            status = linalg_status_no_convergence
            return
        end if
        scale = maxval(abs(a))
        pivot_tol = max(32.0_dp * tiny(1.0_dp), 8.0_dp * epsilon(1.0_dp) * scale)
        do k = 1, n
            pivot = k - 1 + maxloc(abs(a(k:n, k)), dim=1)
            if (abs(a(pivot, k)) <= pivot_tol) then
                status = linalg_status_singular
                return
            end if
            if (pivot /= k) then
                do j = k, n
                    temp = a(k, j)
                    a(k, j) = a(pivot, j)
                    a(pivot, j) = temp
                end do
                temp = b(k)
                b(k) = b(pivot)
                b(pivot) = temp
            end if
            do i = k + 1, n
                factor = a(i, k) / a(k, k)
                a(i, k) = cmplx(0.0_dp, 0.0_dp, kind=dp)
                a(i, k + 1:n) = a(i, k + 1:n) - factor * a(k, k + 1:n)
                b(i) = b(i) - factor * b(k)
            end do
        end do
        do i = n, 1, -1
            x(i) = (b(i) - sum(a(i, i + 1:n) * x(i + 1:n))) / a(i, i)
        end do
        if (.not. all(ieee_is_finite(real(x, dp))) .or. &
                .not. all(ieee_is_finite(aimag(x)))) then
            status = linalg_status_no_convergence
            x = cmplx(0.0_dp, 0.0_dp, kind=dp)
            return
        end if
        status = linalg_status_success
    end subroutine solve_pivoted_complex

end module scifort_linalg_matrix_equations
