program test_integrate_sparse_lu
    use scifort_kinds, only : dp
    use scifort_integrate_sparse_jacobian, only : ivp_sparse_matrix, &
        ivp_complex_sparse_matrix, ivp_csc_matrix
    use scifort_integrate_sparse_lu, only : ivp_sparse_lu_real, &
        ivp_sparse_lu_complex, sparse_lu_factor
    implicit none

    call test_real_tridiagonal()
    call test_real_pivoting()
    call test_real_to_complex_shift()
    call test_complex_factor()
    call test_singular_detection()
    print *, 'test_integrate_sparse_lu: PASS'

contains

    subroutine test_real_tridiagonal()
        integer, parameter :: n = 40
        integer, parameter :: nnz = 3 * n - 2
        integer :: indptr(n + 1)
        integer :: indices(nnz)
        real(dp) :: data(nnz)
        real(dp) :: b(n)
        real(dp) :: expected(n)
        real(dp) :: x(n)
        type(ivp_sparse_matrix) :: matrix
        type(ivp_sparse_lu_real) :: lu
        integer :: i
        integer :: info
        integer :: k

        call make_real_tridiagonal(n, indptr, indices, data)
        matrix = ivp_csc_matrix(n, n, indptr, indices, data)
        expected = [(sin(real(i, dp)), i = 1, n)]
        call real_sparse_matvec(matrix, expected, b)
        call sparse_lu_factor(matrix, 0.0_dp, 1.0_dp, lu, info)
        call require(info == 0, 'real tridiagonal factorization')
        call lu%solve(b, x, info)
        call require(info == 0, 'real tridiagonal solve')
        call require(maxval(abs(x - expected)) < 2.0e-13_dp, &
            'real tridiagonal solution accuracy')
        k = lu%nnz()
        call require(k < n * n / 4, 'real tridiagonal sparse fill')
    end subroutine test_real_tridiagonal

    subroutine test_real_pivoting()
        type(ivp_sparse_matrix) :: matrix
        type(ivp_sparse_lu_real) :: lu
        real(dp) :: b(3)
        real(dp) :: expected(3)
        real(dp) :: x(3)
        integer :: info

        matrix = ivp_csc_matrix(3, 3, [1, 2, 5, 7], [2, 1, 2, 3, 2, 3], &
            [1.0_dp, 2.0_dp, 3.0_dp, 5.0_dp, 4.0_dp, 6.0_dp])
        expected = [1.0_dp, -2.0_dp, 0.5_dp]
        call real_sparse_matvec(matrix, expected, b)
        call sparse_lu_factor(matrix, 0.0_dp, 1.0_dp, lu, info)
        call require(info == 0, 'real pivoting factorization')
        call require(lu%piv(1) == 2, 'real pivoting row selection')
        call lu%solve(b, x, info)
        call require(info == 0, 'real pivoting solve')
        call require(maxval(abs(x - expected)) < 2.0e-13_dp, &
            'real pivoting solution accuracy')
    end subroutine test_real_pivoting

    subroutine test_real_to_complex_shift()
        type(ivp_sparse_matrix) :: jac
        type(ivp_sparse_lu_complex) :: lu
        complex(dp) :: b(3)
        complex(dp) :: expected(3)
        complex(dp) :: x(3)
        complex(dp), parameter :: diagonal = cmplx(2.5_dp, -0.75_dp, kind=dp)
        complex(dp), parameter :: jac_scale = cmplx(-1.0_dp, 0.0_dp, kind=dp)
        integer :: info

        jac = ivp_csc_matrix(3, 3, [1, 3, 6, 8], [1, 2, 1, 2, 3, 2, 3], &
            [-2.0_dp, 1.0_dp, 1.0_dp, -3.0_dp, 1.0_dp, 1.0_dp, -2.0_dp])
        expected = [cmplx(1.0_dp, 0.5_dp, dp), cmplx(-0.5_dp, 1.0_dp, dp), &
            cmplx(2.0_dp, -0.25_dp, dp)]
        call shifted_real_sparse_matvec(jac, diagonal, jac_scale, expected, b)
        call sparse_lu_factor(jac, diagonal, jac_scale, lu, info)
        call require(info == 0, 'real-to-complex shifted factorization')
        call lu%solve(b, x, info)
        call require(info == 0, 'real-to-complex shifted solve')
        call require(maxval(abs(x - expected)) < 3.0e-13_dp, &
            'real-to-complex shifted solution accuracy')
    end subroutine test_real_to_complex_shift

    subroutine test_complex_factor()
        type(ivp_complex_sparse_matrix) :: jac
        type(ivp_sparse_lu_complex) :: lu
        complex(dp) :: b(3)
        complex(dp) :: expected(3)
        complex(dp) :: x(3)
        complex(dp), parameter :: diagonal = cmplx(1.0_dp, 0.0_dp, kind=dp)
        complex(dp), parameter :: jac_scale = cmplx(-0.25_dp, 0.0_dp, kind=dp)
        integer :: info

        jac = ivp_csc_matrix(3, 3, [1, 3, 6, 8], [1, 2, 1, 2, 3, 2, 3], &
            [cmplx(-3.0_dp, 1.0_dp, dp), cmplx(1.0_dp, -0.5_dp, dp), &
             cmplx(2.0_dp, 0.25_dp, dp), cmplx(-4.0_dp, -1.0_dp, dp), &
             cmplx(1.0_dp, 0.0_dp, dp), cmplx(1.0_dp, 0.75_dp, dp), &
             cmplx(-2.0_dp, 0.5_dp, dp)])
        expected = [cmplx(0.5_dp, 1.0_dp, dp), cmplx(-1.0_dp, 0.25_dp, dp), &
            cmplx(1.5_dp, -0.75_dp, dp)]
        call shifted_complex_sparse_matvec(jac, diagonal, jac_scale, expected, b)
        call sparse_lu_factor(jac, diagonal, jac_scale, lu, info)
        call require(info == 0, 'complex sparse factorization')
        call lu%solve(b, x, info)
        call require(info == 0, 'complex sparse solve')
        call require(maxval(abs(x - expected)) < 3.0e-13_dp, &
            'complex sparse solution accuracy')
    end subroutine test_complex_factor

    subroutine test_singular_detection()
        type(ivp_sparse_matrix) :: matrix
        type(ivp_sparse_lu_real) :: lu
        integer :: info

        matrix = ivp_csc_matrix(2, 2, [1, 1, 1], [integer ::], [real(dp) ::])
        call sparse_lu_factor(matrix, 0.0_dp, 1.0_dp, lu, info)
        call require(info == 1, 'sparse singular detection')
    end subroutine test_singular_detection

    subroutine make_real_tridiagonal(n, indptr, indices, data)
        integer, intent(in) :: n !! square tridiagonal matrix order, at least 2
        integer, intent(out) :: indptr(n + 1) !! one-based CSC column pointers
        integer, intent(out) :: indices(3 * n - 2) !! one-based CSC row indices
        real(dp), intent(out) :: data(3 * n - 2) !! tridiagonal matrix values

        integer :: j
        integer :: k

        k = 1
        indptr(1) = 1
        do j = 1, n
            if (j > 1) then
                indices(k) = j - 1
                data(k) = -1.0_dp
                k = k + 1
            end if
            indices(k) = j
            data(k) = 4.0_dp
            k = k + 1
            if (j < n) then
                indices(k) = j + 1
                data(k) = -1.0_dp
                k = k + 1
            end if
            indptr(j + 1) = k
        end do
    end subroutine make_real_tridiagonal

    subroutine real_sparse_matvec(matrix, x, y)
        type(ivp_sparse_matrix), intent(in) :: matrix !! real CSC matrix
        real(dp), intent(in) :: x(:) !! input vector, length ncols
        real(dp), intent(out) :: y(:) !! matrix-vector product, length nrows

        integer :: j
        integer :: k

        y = 0.0_dp
        do j = 1, matrix%ncols
            do k = matrix%indptr(j), matrix%indptr(j + 1) - 1
                y(matrix%indices(k)) = y(matrix%indices(k)) + matrix%data(k) * x(j)
            end do
        end do
    end subroutine real_sparse_matvec

    subroutine shifted_real_sparse_matvec(matrix, diagonal, jac_scale, x, y)
        type(ivp_sparse_matrix), intent(in) :: matrix !! real CSC Jacobian
        complex(dp), intent(in) :: diagonal !! diagonal shift applied to the identity
        complex(dp), intent(in) :: jac_scale !! multiplier applied to Jacobian entries
        complex(dp), intent(in) :: x(:) !! complex input vector
        complex(dp), intent(out) :: y(:) !! shifted sparse matrix-vector product

        integer :: j
        integer :: k

        y = diagonal * x
        do j = 1, matrix%ncols
            do k = matrix%indptr(j), matrix%indptr(j + 1) - 1
                y(matrix%indices(k)) = y(matrix%indices(k)) + &
                    jac_scale * matrix%data(k) * x(j)
            end do
        end do
    end subroutine shifted_real_sparse_matvec

    subroutine shifted_complex_sparse_matvec(matrix, diagonal, jac_scale, x, y)
        type(ivp_complex_sparse_matrix), intent(in) :: matrix !! complex CSC Jacobian
        complex(dp), intent(in) :: diagonal !! diagonal shift applied to the identity
        complex(dp), intent(in) :: jac_scale !! multiplier applied to Jacobian entries
        complex(dp), intent(in) :: x(:) !! complex input vector
        complex(dp), intent(out) :: y(:) !! shifted sparse matrix-vector product

        integer :: j
        integer :: k

        y = diagonal * x
        do j = 1, matrix%ncols
            do k = matrix%indptr(j), matrix%indptr(j + 1) - 1
                y(matrix%indices(k)) = y(matrix%indices(k)) + &
                    jac_scale * matrix%data(k) * x(j)
            end do
        end do
    end subroutine shifted_complex_sparse_matvec

    subroutine require(condition, label)
        logical, intent(in) :: condition !! assertion condition
        character(len=*), intent(in) :: label !! assertion label for failure output

        if (.not. condition) then
            write (*, '(a)') 'FAIL '//trim(label)
            error stop 1
        end if
    end subroutine require

end program test_integrate_sparse_lu
