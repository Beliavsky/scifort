! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
!
! Rank-revealing small dense least squares for augmented Krylov systems.
! Classical one-sided Jacobi orthogonalizes the columns of the input matrix,
! without forming A**T*A, which squares the spectral condition number.
!
module scifort_optimize_krylov_rank_ls
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private
    public :: krylov_rank_least_squares

contains

    subroutine krylov_rank_least_squares(a, rhs, solution, rank, ok)
        real(dp), intent(in) :: a(:, :) !! finite m-by-n least-squares matrix, m >= n
        real(dp), intent(in) :: rhs(:) !! finite right-hand side, length m
        real(dp), intent(out) :: solution(:) !! truncated-SVD minimum-norm solution, length n
        integer, intent(out) :: rank !! estimated numerical rank
        logical, intent(out) :: ok !! true when finite least-squares coefficients obtained

        real(dp), allocatable :: work(:, :)
        real(dp), allocatable :: right(:, :)
        real(dp), allocatable :: pcol(:)
        real(dp), allocatable :: qcol(:)
        real(dp) :: c
        real(dp) :: coeff
        real(dp) :: cosine
        real(dp) :: norm_p
        real(dp) :: norm_q
        real(dp) :: rank_tol
        real(dp) :: scale
        real(dp) :: s
        real(dp) :: sigma
        real(dp) :: sigma_max
        real(dp) :: t
        real(dp) :: tau
        real(dp) :: tmp
        integer :: i
        integer :: m
        integer :: n
        integer :: p
        integer :: q
        integer :: sweep
        logical :: changed
        logical :: converged

        m = size(a, 1)
        n = size(a, 2)
        solution = 0.0_dp
        rank = 0
        ok = .false.
        if (m < n .or. n < 1) return
        if (size(rhs) /= m .or. size(solution) /= n) return
        if (.not. all(ieee_is_finite(a)) .or. .not. all(ieee_is_finite(rhs))) return

        scale = maxval(abs(a))
        if (scale == 0.0_dp) then
            ok = .true.
            return
        end if
        ! Uniform scaling bounds the matrix entries and permits rank decisions
        ! from singular values rather than squared singular values.
        allocate(work(m, n), right(n, n), pcol(m), qcol(m))
        work = a / scale
        right = 0.0_dp
        do i = 1, n
            right(i, i) = 1.0_dp
        end do
        sigma_max = 0.0_dp
        do i = 1, n
            sigma_max = max(sigma_max, column_norm(work(:, i)))
        end do
        rank_tol = real(max(m, n), dp) * epsilon(1.0_dp) * sigma_max

        ! Each plane rotation orthogonalizes a pair of columns and applies the
        ! same rotation to right singular vectors. We skip columns already
        ! below the eventual numerical-rank cutoff.
        converged = .false.
        do sweep = 1, max(30, min(100, 3 * n))
            changed = .false.
            do p = 1, n - 1
                norm_p = column_norm(work(:, p))
                if (norm_p <= rank_tol) cycle
                do q = p + 1, n
                    norm_q = column_norm(work(:, q))
                    if (norm_q <= rank_tol) cycle
                    cosine = dot_product(work(:, p) / norm_p, work(:, q) / norm_q)
                    if (.not. ieee_is_finite(cosine)) return
                    if (abs(cosine) <= 8.0_dp * epsilon(1.0_dp)) cycle

                    ! Avoid squaring column norms when calculating the Jacobi
                    ! angle; the norm ratio is bounded for retained columns.
                    tau = 0.5_dp * (norm_q / norm_p - norm_p / norm_q) / cosine
                    if (.not. ieee_is_finite(tau)) return
                    t = sign(1.0_dp, tau) / (abs(tau) + hypot(1.0_dp, tau))
                    c = 1.0_dp / hypot(1.0_dp, t)
                    s = t * c
                    pcol = work(:, p)
                    qcol = work(:, q)
                    work(:, p) = c * pcol - s * qcol
                    work(:, q) = s * pcol + c * qcol
                    pcol(1:n) = right(:, p)
                    qcol(1:n) = right(:, q)
                    right(:, p) = c * pcol(1:n) - s * qcol(1:n)
                    right(:, q) = s * pcol(1:n) + c * qcol(1:n)
                    changed = .true.
                end do
            end do
            if (.not. changed) then
                converged = .true.
                exit
            end if
        end do
        if (.not. converged) return

        sigma_max = 0.0_dp
        do i = 1, n
            sigma_max = max(sigma_max, column_norm(work(:, i)))
        end do
        rank_tol = real(max(m, n), dp) * epsilon(1.0_dp) * sigma_max
        do i = 1, n
            sigma = column_norm(work(:, i))
            if (sigma <= rank_tol) cycle
            rank = rank + 1
            ! Project onto the normalized left singular vector. Division by
            ! sigma * scale applies the inverse singular value only once.
            tmp = dot_product(work(:, i) / sigma, rhs)
            if (.not. ieee_is_finite(tmp)) return
            ! Divide by the larger denominator first. This avoids an
            ! overflowing intermediate for matrices with huge common scale,
            ! even when the final coefficient is representable.
            if (scale >= 1.0_dp) then
                tmp = tmp / scale
                if (sigma < 1.0_dp) then
                    if (abs(tmp) > huge(1.0_dp) * sigma) return
                end if
                coeff = tmp / sigma
            else
                if (sigma < 1.0_dp) then
                    if (abs(tmp) > huge(1.0_dp) * sigma) return
                end if
                tmp = tmp / sigma
                if (abs(tmp) > huge(1.0_dp) * scale) return
                coeff = tmp / scale
            end if
            if (.not. ieee_is_finite(coeff)) return
            solution = solution + coeff * right(:, i)
            if (.not. all(ieee_is_finite(solution))) return
        end do
        ok = .true.
    end subroutine krylov_rank_least_squares

    pure function column_norm(values) result(value)
        real(dp), intent(in) :: values(:) !! real vector to measure
        real(dp) :: value !! overflow-safe Euclidean norm
        real(dp) :: largest

        largest = maxval(abs(values))
        if (largest == 0.0_dp) then
            value = 0.0_dp
        else
            value = largest * sqrt(sum((values / largest)**2))
        end if
    end function column_norm

end module scifort_optimize_krylov_rank_ls
