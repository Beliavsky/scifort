! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Small dense linear-algebra kernels used by the multivariate statistics code.
! The Cholesky interface and failure semantics were checked against the
! BSD-3-Clause fortran-lapack project supplied with the SciFort translation
! work. The implementation below is an independent compact unblocked kernel.

module scifort_linalg
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: linalg_status_success = 0
    integer, parameter, public :: linalg_status_invalid_shape = 1
    integer, parameter, public :: linalg_status_not_positive_definite = 2
    integer, parameter, public :: linalg_status_not_positive_semidefinite = 3
    integer, parameter, public :: linalg_status_no_convergence = 4
    integer, parameter, public :: linalg_status_singular = 5

    public :: cholesky_lower
    public :: psd_decompose
    public :: solve_lower
    public :: symmetric_eigen_jacobi

contains

    pure subroutine cholesky_lower(a, l, status)
        real(dp), intent(in) :: a(:, :) !! square symmetric matrix; only the lower triangle is read
        real(dp), intent(out) :: l(:, :) !! lower Cholesky factor, A = L L^T
        integer, intent(out) :: status !! zero on success; nonzero for invalid/non-SPD input

        integer :: i
        integer :: j
        integer :: n
        real(dp) :: diagonal
        real(dp) :: work

        l = 0.0_dp
        n = size(a, 1)
        if (n < 1 .or. size(a, 2) /= n .or. size(l, 1) /= n .or. size(l, 2) /= n) then
            status = linalg_status_invalid_shape
            return
        end if

        do j = 1, n
            if (.not. ieee_is_finite(a(j, j))) then
                status = linalg_status_not_positive_definite
                return
            end if
            diagonal = a(j, j)
            if (j > 1) diagonal = diagonal - dot_product(l(j, 1:j - 1), l(j, 1:j - 1))
            if (.not. ieee_is_finite(diagonal) .or. diagonal <= 0.0_dp) then
                status = linalg_status_not_positive_definite
                return
            end if
            l(j, j) = sqrt(diagonal)

            do i = j + 1, n
                if (.not. ieee_is_finite(a(i, j))) then
                    status = linalg_status_not_positive_definite
                    return
                end if
                work = a(i, j)
                if (j > 1) work = work - dot_product(l(i, 1:j - 1), l(j, 1:j - 1))
                l(i, j) = work / l(j, j)
            end do
        end do
        status = linalg_status_success
    end subroutine cholesky_lower

    pure subroutine solve_lower(l, b, x, status)
        real(dp), intent(in) :: l(:, :) !! nonsingular lower-triangular coefficient matrix
        real(dp), intent(in) :: b(:) !! right-hand side vector
        real(dp), intent(out) :: x(:) !! solution of L x = b
        integer, intent(out) :: status !! zero on success; nonzero for invalid/singular input

        integer :: i
        integer :: n
        real(dp) :: rhs

        x = 0.0_dp
        n = size(b)
        if (n < 1 .or. size(x) /= n .or. size(l, 1) /= n .or. size(l, 2) /= n) then
            status = linalg_status_invalid_shape
            return
        end if

        do i = 1, n
            if (.not. ieee_is_finite(l(i, i)) .or. l(i, i) == 0.0_dp) then
                status = linalg_status_singular
                return
            end if
            rhs = b(i)
            if (i > 1) rhs = rhs - dot_product(l(i, 1:i - 1), x(1:i - 1))
            x(i) = rhs / l(i, i)
        end do
        status = linalg_status_success
    end subroutine solve_lower

    pure subroutine symmetric_eigen_jacobi(a, eigenvalues, eigenvectors, status)
        real(dp), intent(in) :: a(:, :) !! square symmetric matrix; only the lower triangle is read
        real(dp), intent(out) :: eigenvalues(:) !! eigenvalues in ascending order
        real(dp), intent(out) :: eigenvectors(:, :) !! corresponding orthonormal eigenvectors by column
        integer, intent(out) :: status !! zero on convergence; nonzero otherwise

        integer, parameter :: max_sweeps = 100
        integer :: i
        integer :: k
        integer :: n
        integer :: p
        integer :: q
        integer :: sweep
        real(dp), allocatable :: b(:, :)
        real(dp) :: app
        real(dp) :: apq
        real(dp) :: aqq
        real(dp) :: brp
        real(dp) :: brq
        real(dp) :: c
        real(dp) :: max_off
        real(dp) :: norm_a
        real(dp) :: s
        real(dp) :: t
        real(dp) :: tau
        real(dp) :: tol
        real(dp) :: vrp
        real(dp) :: vrq

        n = size(a, 1)
        eigenvalues = 0.0_dp
        eigenvectors = 0.0_dp
        if (n < 1 .or. size(a, 2) /= n .or. size(eigenvalues) /= n .or. &
                size(eigenvectors, 1) /= n .or. size(eigenvectors, 2) /= n) then
            status = linalg_status_invalid_shape
            return
        end if

        allocate(b(n, n))
        do i = 1, n
            if (.not. ieee_is_finite(a(i, i))) then
                status = linalg_status_no_convergence
                return
            end if
            b(i, i) = a(i, i)
        end do
        do p = 2, n
            do q = 1, p - 1
                if (.not. ieee_is_finite(a(p, q))) then
                    status = linalg_status_no_convergence
                    return
                end if
                b(p, q) = a(p, q)
                b(q, p) = a(p, q)
            end do
        end do

        do i = 1, n
            eigenvectors(i, i) = 1.0_dp
        end do
        norm_a = max(1.0_dp, maxval(abs(b)))
        tol = 32.0_dp * epsilon(1.0_dp) * norm_a

        do sweep = 1, max_sweeps
            max_off = 0.0_dp
            do p = 1, n - 1
                do q = p + 1, n
                    apq = b(p, q)
                    max_off = max(max_off, abs(apq))
                    if (abs(apq) <= tol) cycle

                    app = b(p, p)
                    aqq = b(q, q)
                    tau = (aqq - app) / (2.0_dp * apq)
                    if (tau >= 0.0_dp) then
                        t = 1.0_dp / (tau + sqrt(1.0_dp + tau * tau))
                    else
                        t = -1.0_dp / (-tau + sqrt(1.0_dp + tau * tau))
                    end if
                    c = 1.0_dp / sqrt(1.0_dp + t * t)
                    s = t * c

                    b(p, p) = app - t * apq
                    b(q, q) = aqq + t * apq
                    b(p, q) = 0.0_dp
                    b(q, p) = 0.0_dp

                    do k = 1, n
                        if (k /= p .and. k /= q) then
                            brp = b(k, p)
                            brq = b(k, q)
                            b(k, p) = c * brp - s * brq
                            b(p, k) = b(k, p)
                            b(k, q) = s * brp + c * brq
                            b(q, k) = b(k, q)
                        end if

                        vrp = eigenvectors(k, p)
                        vrq = eigenvectors(k, q)
                        eigenvectors(k, p) = c * vrp - s * vrq
                        eigenvectors(k, q) = s * vrp + c * vrq
                    end do
                end do
            end do
            if (max_off <= tol) exit
        end do

        if (max_off > tol) then
            status = linalg_status_no_convergence
            return
        end if

        do i = 1, n
            eigenvalues(i) = b(i, i)
        end do
        call sort_eigenpairs(eigenvalues, eigenvectors)
        status = linalg_status_success
    end subroutine symmetric_eigen_jacobi

    pure subroutine psd_decompose(a, eigenvalues, eigenvectors, rank, log_pdet, cutoff, &
            status, allow_singular)
        real(dp), intent(in) :: a(:, :) !! symmetric covariance-like matrix; lower triangle is authoritative
        real(dp), intent(out) :: eigenvalues(:) !! ascending eigenvalues
        real(dp), intent(out) :: eigenvectors(:, :) !! orthonormal eigenvectors by column
        integer, intent(out) :: rank !! numerical rank using SciPy-compatible binary64 cutoff
        real(dp), intent(out) :: log_pdet !! logarithm of the pseudo-determinant
        real(dp), intent(out) :: cutoff !! eigenvalue magnitude cutoff for numerical rank
        integer, intent(out) :: status !! zero on success; nonzero for non-PSD/singular input
        logical, intent(in), optional :: allow_singular !! permit PSD rank deficiency (default true)

        integer :: i
        logical :: singular_ok
        real(dp) :: largest

        singular_ok = .true.
        if (present(allow_singular)) singular_ok = allow_singular

        call symmetric_eigen_jacobi(a, eigenvalues, eigenvectors, status)
        if (status /= linalg_status_success) then
            rank = 0
            log_pdet = 0.0_dp
            cutoff = 0.0_dp
            return
        end if

        largest = maxval(abs(eigenvalues))
        cutoff = 1.0e6_dp * epsilon(1.0_dp) * largest
        if (minval(eigenvalues) < -cutoff) then
            rank = 0
            log_pdet = 0.0_dp
            status = linalg_status_not_positive_semidefinite
            return
        end if

        rank = count(eigenvalues > cutoff)
        if (rank < size(eigenvalues) .and. .not. singular_ok) then
            log_pdet = 0.0_dp
            status = linalg_status_singular
            return
        end if

        log_pdet = 0.0_dp
        do i = 1, size(eigenvalues)
            if (eigenvalues(i) > cutoff) log_pdet = log_pdet + log(eigenvalues(i))
        end do
        status = linalg_status_success
    end subroutine psd_decompose

    pure subroutine sort_eigenpairs(eigenvalues, eigenvectors)
        real(dp), intent(inout) :: eigenvalues(:) !! values sorted in ascending order in place
        real(dp), intent(inout) :: eigenvectors(:, :) !! columns permuted with eigenvalues

        integer :: i
        integer :: j
        integer :: k
        real(dp) :: tmp
        real(dp) :: tmp_col(size(eigenvectors, 1))

        do i = 1, size(eigenvalues) - 1
            k = i
            do j = i + 1, size(eigenvalues)
                if (eigenvalues(j) < eigenvalues(k)) k = j
            end do
            if (k /= i) then
                tmp = eigenvalues(i)
                eigenvalues(i) = eigenvalues(k)
                eigenvalues(k) = tmp
                tmp_col = eigenvectors(:, i)
                eigenvectors(:, i) = eigenvectors(:, k)
                eigenvectors(:, k) = tmp_col
            end if
        end do
    end subroutine sort_eigenpairs

end module scifort_linalg
