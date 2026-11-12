! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Rectangular real/complex polar decomposition, scipy.linalg.polar style.
! Independent one-sided complex Jacobi method, without BLAS/LAPACK.
module scifort_linalg_polar
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private
    integer, parameter, public :: polar_max_min_dimension = 64
    integer, parameter, public :: polar_max_max_dimension = 256
    integer, parameter :: default_sweeps = 120
    public :: polar

    interface polar
        module procedure polar_real
        module procedure polar_complex
    end interface polar

contains

    subroutine polar_real(a, u, p, status, side, residual, orthogonality, max_sweeps)
        real(dp), intent(in) :: a(:, :) !! input real matrix (m,n)
        real(dp), intent(out) :: u(:, :) !! partial orthogonal factor (m,n)
        real(dp), intent(out) :: p(:, :) !! symmetric PSD factor (n,n) or (m,m)
        integer, intent(out) :: status !! linalg_status_* result
        character(*), intent(in), optional :: side !! 'right' (default) or 'left'
        real(dp), intent(out), optional :: residual !! relative reconstruction error
        real(dp), intent(out), optional :: orthogonality !! orthogonal-factor error
        integer, intent(in), optional :: max_sweeps !! number of Jacobi sweeps, default 120
        complex(dp), allocatable :: uc(:, :), pc(:, :)
        real(dp) :: scale

        u = 0.0_dp
        p = 0.0_dp
        allocate(uc(size(u, 1), size(u, 2)), pc(size(p, 1), size(p, 2)))
        call polar_complex(cmplx(a, kind=dp), uc, pc, status, side, residual, orthogonality, max_sweeps)
        if (status /= linalg_status_success) return
        scale = max(1.0_dp, maxval(abs(uc)), maxval(abs(pc)))
        if (max(maxval(abs(aimag(uc))), maxval(abs(aimag(pc)))) > 2.0e-10_dp*scale) then
            status = linalg_status_no_convergence
            return
        end if
        u = real(uc, dp)
        p = real(pc, dp)
    end subroutine polar_real

    subroutine polar_complex(a, u, p, status, side, residual, orthogonality, max_sweeps)
        complex(dp), intent(in) :: a(:, :) !! input complex matrix (m,n)
        complex(dp), intent(out) :: u(:, :) !! partial unitary factor (m,n)
        complex(dp), intent(out) :: p(:, :) !! Hermitian PSD factor (n,n) or (m,m)
        integer, intent(out) :: status !! linalg_status_* result
        character(*), intent(in), optional :: side !! 'right' (default) or 'left'
        real(dp), intent(out), optional :: residual !! relative reconstruction error
        real(dp), intent(out), optional :: orthogonality !! partial-unitarity error
        integer, intent(in), optional :: max_sweeps !! number of Jacobi sweeps, default 120
        complex(dp), allocatable :: adjoint(:, :), uadj(:, :), right_p(:, :), projector(:, :)
        real(dp) :: reconstruct_err, orth_err, denom, root_scale
        integer :: m, n, i, k, sweeps
        logical :: left_side

        u = cmplx(0.0_dp, 0.0_dp, dp)
        p = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(residual)) residual = huge(1.0_dp)
        if (present(orthogonality)) orthogonality = huge(1.0_dp)
        status = linalg_status_invalid_shape
        m = size(a, 1)
        n = size(a, 2)
        if (m < 1 .or. n < 1 .or. min(m, n) > polar_max_min_dimension .or. &
                max(m, n) > polar_max_max_dimension) return
        if (size(u, 1) /= m .or. size(u, 2) /= n) return
        left_side = .false.
        if (present(side)) then
            select case (side)
            case ('left')
                left_side = .true.
            case ('right')
                left_side = .false.
            case default
                return
            end select
        end if
        if (left_side) then
            if (size(p, 1) /= m .or. size(p, 2) /= m) return
        else
            if (size(p, 1) /= n .or. size(p, 2) /= n) return
        end if
        sweeps = default_sweeps
        if (present(max_sweeps)) sweeps = max_sweeps
        if (sweeps < 1) return
        if (.not. all(ieee_is_finite(real(a, dp))) .or. &
                .not. all(ieee_is_finite(aimag(a)))) return

        if (m >= n) then
            allocate(right_p(n, n))
            call tall_polar(a, u, right_p, status, sweeps)
        else
            allocate(adjoint(n, m), uadj(n, m), right_p(m, m))
            adjoint = transpose(conjg(a))
            call tall_polar(adjoint, uadj, right_p, status, sweeps)
            if (status == linalg_status_success) u = transpose(conjg(uadj))
        end if
        if (status /= linalg_status_success) return

        if (left_side) then
            ! A U^H = U (A^H A)^(1/2) U^H, a Hermitian PSD matrix.
            p = matmul(a, transpose(conjg(u)))
        else
            p = matmul(transpose(conjg(u)), a)
        end if
        p = 0.5_dp*(p+transpose(conjg(p)))

        if (left_side) then
            reconstruct_err = norm_fro(matmul(p, u)-a)
        else
            reconstruct_err = norm_fro(matmul(u, p)-a)
        end if
        denom = max(tiny(1.0_dp), norm_fro(a))
        reconstruct_err = reconstruct_err/denom
        k = min(m, n)
        allocate(projector(k, k))
        if (m >= n) then
            projector = matmul(transpose(conjg(u)), u)
        else
            projector = matmul(u, transpose(conjg(u)))
        end if
        do i = 1, k
            projector(i, i) = projector(i, i)-cmplx(1.0_dp, 0.0_dp, dp)
        end do
        orth_err = norm_fro(projector)/sqrt(real(k, dp))
        if (present(residual)) residual = reconstruct_err
        if (present(orthogonality)) orthogonality = orth_err
        root_scale = max(1.0_dp, maxval(abs(p)))
        if (.not. ieee_is_finite(reconstruct_err) .or. .not. ieee_is_finite(orth_err) .or. &
                .not. all(ieee_is_finite(real(p/root_scale, dp))) .or. &
                .not. all(ieee_is_finite(aimag(p/root_scale))) .or. &
                reconstruct_err > 1.0e-8_dp .or. orth_err > 1.0e-8_dp) then
            status = linalg_status_no_convergence
            u = cmplx(0.0_dp, 0.0_dp, dp)
            p = cmplx(0.0_dp, 0.0_dp, dp)
        end if
    end subroutine polar_complex

    subroutine tall_polar(a, u, p, status, limit)
        complex(dp), intent(in) :: a(:, :) !! tall input matrix, m>=n
        complex(dp), intent(out) :: u(:, :) !! partial unitary polar factor (m,n)
        complex(dp), intent(out) :: p(:, :) !! Hermitian PSD factor (n,n)
        integer, intent(out) :: status !! linalg_status_* result
        integer, intent(in) :: limit !! maximum Jacobi sweeps
        complex(dp), allocatable :: columns(:, :), right(:, :), left(:, :), tmp(:, :)
        complex(dp), allocatable :: column_work(:), right_work(:), best_work(:)
        real(dp), allocatable :: sigma(:)
        complex(dp) :: gamma, phase
        real(dp) :: a2, b2, g, angle, cosine, sine, scale, threshold, rank_tol
        real(dp) :: magnitude, best_magnitude, temp_sigma
        integer :: m, n, i, j, k, sweep, pass, best, pivot
        logical :: changed, converged

        m = size(a, 1)
        n = size(a, 2)
        status = linalg_status_no_convergence
        u = cmplx(0.0_dp, 0.0_dp, dp)
        p = cmplx(0.0_dp, 0.0_dp, dp)
        scale = maxval(abs(a))
        if (scale == 0.0_dp) then
            do i = 1, n
                u(i, i) = cmplx(1.0_dp, 0.0_dp, dp)
            end do
            status = linalg_status_success
            return
        end if
        allocate(columns(m, n), right(n, n), left(m, n), tmp(n, n), &
            column_work(m), right_work(n), best_work(m), sigma(n))
        columns = a/scale
        right = cmplx(0.0_dp, 0.0_dp, dp)
        do i = 1, n
            right(i, i) = cmplx(1.0_dp, 0.0_dp, dp)
        end do
        threshold = 64.0_dp*epsilon(1.0_dp)
        converged = .false.
        do sweep = 1, limit
            changed = .false.
            do i = 1, n-1
                do j = i+1, n
                    a2 = sum(abs(columns(:, i))**2)
                    b2 = sum(abs(columns(:, j))**2)
                    if (a2 == 0.0_dp .or. b2 == 0.0_dp) cycle
                    gamma = dot_product(columns(:, i), columns(:, j))
                    g = abs(gamma)
                    if (g <= threshold*sqrt(a2)*sqrt(b2)) cycle
                    phase = conjg(gamma)/g
                    if (b2 >= a2) then
                        angle = 0.5_dp*atan2(2.0_dp*g, b2-a2)
                    else
                        angle = -0.5_dp*atan2(2.0_dp*g, a2-b2)
                    end if
                    if (abs(angle) <= epsilon(1.0_dp)) cycle
                    cosine = cos(angle)
                    sine = sin(angle)
                    column_work = cosine*columns(:, i)-sine*phase*columns(:, j)
                    columns(:, j) = sine*columns(:, i)+cosine*phase*columns(:, j)
                    columns(:, i) = column_work
                    right_work = cosine*right(:, i)-sine*phase*right(:, j)
                    right(:, j) = sine*right(:, i)+cosine*phase*right(:, j)
                    right(:, i) = right_work
                    changed = .true.
                end do
            end do
            if (.not. changed) then
                converged = .true.
                exit
            end if
        end do
        if (.not. converged) return
        do j = 1, n
            sigma(j) = sqrt(sum(abs(columns(:, j))**2))
        end do
        ! Complete the numerical null space only after the nonzero directions.
        do j = 1, n-1
            best = j-1+maxloc(sigma(j:n), dim=1)
            if (best == j) cycle
            temp_sigma = sigma(j)
            sigma(j) = sigma(best)
            sigma(best) = temp_sigma
            column_work = columns(:, j)
            columns(:, j) = columns(:, best)
            columns(:, best) = column_work
            right_work = right(:, j)
            right(:, j) = right(:, best)
            right(:, best) = right_work
        end do
        left = cmplx(0.0_dp, 0.0_dp, dp)
        rank_tol = 64.0_dp*epsilon(1.0_dp)*real(max(m, n), dp)*maxval(sigma)
        do j = 1, n
            if (sigma(j) > rank_tol) then
                column_work = columns(:, j)/sigma(j)
            else
                ! Deterministic orthonormal completion, including zero A.
                best_magnitude = -1.0_dp
                do pivot = 1, m
                    column_work = cmplx(0.0_dp, 0.0_dp, dp)
                    column_work(pivot) = cmplx(1.0_dp, 0.0_dp, dp)
                    do pass = 1, 2
                        do k = 1, j-1
                            column_work = column_work-dot_product(left(:, k), column_work)*left(:, k)
                        end do
                    end do
                    magnitude = sqrt(sum(abs(column_work)**2))
                    if (magnitude > best_magnitude) then
                        best_magnitude = magnitude
                        best_work = column_work
                    end if
                end do
                column_work = best_work
            end if
            do pass = 1, 2
                do k = 1, j-1
                    column_work = column_work-dot_product(left(:, k), column_work)*left(:, k)
                end do
            end do
            magnitude = sqrt(sum(abs(column_work)**2))
            if (magnitude <= 32.0_dp*epsilon(1.0_dp)) return
            left(:, j) = column_work/magnitude
        end do
        u = matmul(left, transpose(conjg(right)))
        do j = 1, n
            tmp(:, j) = sigma(j)*right(:, j)
        end do
        p = matmul(tmp, transpose(conjg(right)))*scale
        p = 0.5_dp*(p+transpose(conjg(p)))
        status = linalg_status_success
    end subroutine tall_polar

    pure real(dp) function norm_fro(a) result(value)
        complex(dp), intent(in) :: a(:, :) !! input matrix
        value = sqrt(sum(abs(a)**2))
    end function norm_fro

end module scifort_linalg_polar
