! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Matrix balancing by similarity permutations and binary-power diagonal scaling.
! Independent implementation of the Parlett-Reinsch balancing principles.
module scifort_linalg_balance
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    public :: matrix_balance

    interface matrix_balance
        module procedure matrix_balance_real
        module procedure matrix_balance_complex
    end interface matrix_balance

contains

    subroutine matrix_balance_real(a, b, status, t, permute, scale, scale_factors, permutation)
        real(dp), intent(in) :: a(:, :) !! finite square input matrix, including (0,0)
        real(dp), intent(out) :: b(:, :) !! balanced matrix B = inv(T) A T
        integer, intent(out) :: status !! success, invalid shape, or numerical failure
        real(dp), intent(out), optional :: t(:, :) !! similarity transformation (n,n)
        logical, intent(in), optional :: permute !! isolate zero-coupled rows/columns (default true)
        logical, intent(in), optional :: scale !! use binary-power diagonal scaling (default true)
        real(dp), intent(out), optional :: scale_factors(:) !! positive diagonal factors (n)
        integer, intent(out), optional :: permutation(:) !! one-based original row for each new column (n)

        integer :: n, j
        real(dp), allocatable :: factors(:)
        integer, allocatable :: perm(:)
        complex(dp), allocatable :: ac(:, :), bc(:, :)

        status = linalg_status_invalid_shape
        b = 0.0_dp
        if (present(t)) t = 0.0_dp
        if (present(scale_factors)) scale_factors = 0.0_dp
        if (present(permutation)) permutation = 0
        n = size(a, 1)
        if (size(a, 2) /= n .or. size(b, 1) /= n .or. size(b, 2) /= n) return
        if (present(t)) then
            if (size(t, 1) /= n .or. size(t, 2) /= n) return
        end if
        if (present(scale_factors)) then
            if (size(scale_factors) /= n) return
        end if
        if (present(permutation)) then
            if (size(permutation) /= n) return
        end if
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(a))) return

        allocate(ac(n,n), bc(n,n), factors(n), perm(n))
        ac = cmplx(a, 0.0_dp, dp)
        call balance_core(ac, bc, factors, perm, status, permute, scale)
        if (status /= linalg_status_success) return
        b = real(bc, dp)
        if (present(t)) then
            do j = 1, n
                t(perm(j), j) = factors(j)
            end do
        end if
        if (present(scale_factors)) scale_factors = factors
        if (present(permutation)) permutation = perm
    end subroutine matrix_balance_real

    subroutine matrix_balance_complex(a, b, status, t, permute, scale, scale_factors, permutation)
        complex(dp), intent(in) :: a(:, :) !! finite square input matrix, including (0,0)
        complex(dp), intent(out) :: b(:, :) !! balanced matrix B = inv(T) A T
        integer, intent(out) :: status !! success, invalid shape, or numerical failure
        complex(dp), intent(out), optional :: t(:, :) !! similarity transformation (n,n)
        logical, intent(in), optional :: permute !! isolate zero-coupled rows/columns (default true)
        logical, intent(in), optional :: scale !! use binary-power diagonal scaling (default true)
        real(dp), intent(out), optional :: scale_factors(:) !! positive diagonal factors (n)
        integer, intent(out), optional :: permutation(:) !! one-based original row for each new column (n)

        integer :: n, j
        real(dp), allocatable :: factors(:)
        integer, allocatable :: perm(:)

        status = linalg_status_invalid_shape
        b = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(t)) t = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(scale_factors)) scale_factors = 0.0_dp
        if (present(permutation)) permutation = 0
        n = size(a, 1)
        if (size(a, 2) /= n .or. size(b, 1) /= n .or. size(b, 2) /= n) return
        if (present(t)) then
            if (size(t, 1) /= n .or. size(t, 2) /= n) return
        end if
        if (present(scale_factors)) then
            if (size(scale_factors) /= n) return
        end if
        if (present(permutation)) then
            if (size(permutation) /= n) return
        end if
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(real(a, dp)))) return
        if (.not. all(ieee_is_finite(aimag(a)))) return

        allocate(factors(n), perm(n))
        call balance_core(a, b, factors, perm, status, permute, scale)
        if (status /= linalg_status_success) then
            b = cmplx(0.0_dp, 0.0_dp, dp)
            return
        end if
        if (present(t)) then
            do j = 1, n
                t(perm(j), j) = cmplx(factors(j), 0.0_dp, dp)
            end do
        end if
        if (present(scale_factors)) scale_factors = factors
        if (present(permutation)) permutation = perm
    end subroutine matrix_balance_complex

    subroutine balance_core(a, b, factors, perm, status, permute, scale)
        complex(dp), intent(in) :: a(:, :) !! finite square input
        complex(dp), intent(out) :: b(:, :) !! balanced output in original units
        real(dp), intent(out) :: factors(:) !! positive binary-power similarity factors
        integer, intent(out) :: perm(:) !! one-based permutation of original rows
        integer, intent(out) :: status !! numerical status
        logical, intent(in), optional :: permute !! permit similarity permutations (default true)
        logical, intent(in), optional :: scale !! permit diagonal scaling (default true)

        integer, parameter :: max_sweeps = 100
        integer, parameter :: max_binary_steps = 512
        real(dp), parameter :: improvement = 0.95_dp
        integer :: i, j, n, low, high, sweep, steps
        real(dp) :: amax, c, r, original_norm, factor, row_max, column_max, safe_limit
        complex(dp), allocatable :: tmp(:)
        logical :: do_permute, do_scale, isolated, changed

        n = size(a, 1)
        b = a
        factors = 1.0_dp
        do i = 1, n
            perm(i) = i
        end do
        status = linalg_status_success
        if (n == 0) return
        do_permute = .true.
        do_scale = .true.
        if (present(permute)) do_permute = permute
        if (present(scale)) do_scale = scale
        low = 1
        high = n
        if (do_permute) then
            allocate(tmp(n))
            ! Isolate uncoupled rows at the bottom; then uncoupled columns at the top.
            do while (low <= high)
                changed = .false.
                do j = high, low, -1
                    isolated = .true.
                    do i = low, high
                        if (i == j) cycle
                        if (b(j, i) /= cmplx(0.0_dp, 0.0_dp, dp)) isolated = .false.
                    end do
                    if (.not. isolated) cycle
                    call swap_similarity(b, perm, tmp, j, high)
                    high = high - 1
                    changed = .true.
                    exit
                end do
                if (changed) cycle
                do j = low, high
                    isolated = .true.
                    do i = low, high
                        if (i == j) cycle
                        if (b(i, j) /= cmplx(0.0_dp, 0.0_dp, dp)) isolated = .false.
                    end do
                    if (.not. isolated) cycle
                    call swap_similarity(b, perm, tmp, j, low)
                    low = low + 1
                    changed = .true.
                    exit
                end do
                if (.not. changed) exit
            end do
        end if
        if (.not. do_scale .or. low >= high) return

        ! Normalize globally before summing row and column magnitudes. This also
        ! handles matrices uniformly scaled to near underflow/overflow levels.
        amax = max(maxval(abs(real(b, dp))), maxval(abs(aimag(b))))
        if (amax == 0.0_dp) return
        b = b / amax
        safe_limit = sqrt(huge(1.0_dp)) * 0.25_dp
        do sweep = 1, max_sweeps
            changed = .false.
            do i = low, high
                r = 0.0_dp
                c = 0.0_dp
                do j = low, high
                    if (i == j) cycle
                    r = r + abs(real(b(i, j), dp)) + abs(aimag(b(i, j)))
                    c = c + abs(real(b(j, i), dp)) + abs(aimag(b(j, i)))
                end do
                if (r == 0.0_dp .or. c == 0.0_dp) cycle
                original_norm = c + r
                factor = 1.0_dp
                steps = 0
                do while (c < 0.5_dp * r .and. steps < max_binary_steps)
                    factor = factor * 2.0_dp
                    c = c * 2.0_dp
                    r = r * 0.5_dp
                    steps = steps + 1
                end do
                do while (r < 0.5_dp * c .and. steps < max_binary_steps)
                    factor = factor * 0.5_dp
                    c = c * 0.5_dp
                    r = r * 2.0_dp
                    steps = steps + 1
                end do
                if (c + r >= improvement * original_norm) cycle
                ! Do not create unrepresentable intermediates or factor values.
                if (factor > 1.0_dp) then
                    if (factors(i) > safe_limit / factor) cycle
                else if (factor < 1.0_dp) then
                    if (factors(i) < (1.0_dp / safe_limit) / factor) cycle
                end if
                row_max = max(maxval(abs(real(b(i, :), dp))), &
                    maxval(abs(aimag(b(i, :)))))
                column_max = max(maxval(abs(real(b(:, i), dp))), &
                    maxval(abs(aimag(b(:, i)))))
                if (factor > 1.0_dp) then
                    if (column_max > safe_limit / factor) cycle
                else if (factor < 1.0_dp) then
                    if (row_max > safe_limit * factor) cycle
                end if
                b(i, :) = b(i, :) / factor
                b(:, i) = b(:, i) * factor
                factors(i) = factors(i) * factor
                changed = .true.
            end do
            if (.not. changed) exit
        end do
        if (changed) then
            status = linalg_status_no_convergence
            return
        end if
        ! Rescale with original global magnitude. All similarity factors are
        ! exact powers of two; overflow is checked before claiming success.
        if (amax > 1.0_dp) then
            if (max(maxval(abs(real(b, dp))), maxval(abs(aimag(b)))) > &
                huge(1.0_dp) / amax) then
                status = linalg_status_no_convergence
                return
            end if
        end if
        b = b * amax
        if (.not. all(ieee_is_finite(real(b, dp))) .or. &
            .not. all(ieee_is_finite(aimag(b)))) status = linalg_status_no_convergence
    end subroutine balance_core

    subroutine swap_similarity(b, perm, tmp, p, q)
        complex(dp), intent(inout) :: b(:, :) !! similarity-permuted square matrix
        integer, intent(inout) :: perm(:) !! accumulated original row indices
        complex(dp), intent(inout) :: tmp(:) !! working row/column buffer
        integer, intent(in) :: p !! first index to swap
        integer, intent(in) :: q !! second index to swap
        integer :: old

        if (p == q) return
        tmp = b(p, :)
        b(p, :) = b(q, :)
        b(q, :) = tmp
        tmp = b(:, p)
        b(:, p) = b(:, q)
        b(:, q) = tmp
        old = perm(p)
        perm(p) = perm(q)
        perm(q) = old
    end subroutine swap_similarity
end module scifort_linalg_balance
