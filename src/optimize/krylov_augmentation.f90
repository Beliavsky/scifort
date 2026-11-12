! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
!
! Optional redundancy filtering for retained LGMRES solution directions.
! The directions themselves are left unchanged: only the independence test
! uses an orthonormal copy of the retained span.
!
module scifort_optimize_krylov_augmentation
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private

    public :: retain_lgmres_direction
    public :: retain_lgmres_image_direction
    public :: lgmres_pruning_workspace
    public :: reset_lgmres_image_pruning_stream
    public :: retain_lgmres_image_direction_stream
    public :: retain_lgmres_image_candidate_cached
    public :: reset_lgmres_direction_pruning_stream
    public :: retain_lgmres_direction_stream
    public :: retain_lgmres_direction_fifo

    ! Optional work arrays for repeated span tests.  A workspace belongs to
    ! one caller/solver, never to module SAVE state.  The traditional entry
    ! points still allocate their own scratch when no workspace is passed.
    type :: lgmres_pruning_workspace
        private
        real(dp), allocatable :: basis(:, :)
        real(dp), allocatable :: candidate(:)
        real(dp), allocatable :: work(:)
        integer :: stream_count = 0
        integer :: stream_basis_count = 0
        integer :: direction_stream_count = 0
        integer :: direction_basis_count = 0
        real(dp), allocatable :: direction_r(:, :)
        integer :: direction_fifo_updates = 0
        integer :: direction_fifo_count = 0
        logical :: direction_fifo_ready = .false.
        logical :: direction_fifo_disabled = .false.
    end type lgmres_pruning_workspace

contains

    subroutine ensure_pruning_workspace(workspace, n, ncols, need_basis)
        type(lgmres_pruning_workspace), intent(inout) :: workspace !! local reusable storage
        integer, intent(in) :: n !! vector length
        integer, intent(in) :: ncols !! active retained directions
        logical, intent(in), optional :: need_basis !! whether to allocate a separate basis
        logical :: own_basis

        own_basis = .true.
        if (present(need_basis)) own_basis = need_basis
        if (own_basis) then
            if (allocated(workspace%basis)) then
                if (size(workspace%basis, 1) /= n .or. &
                        size(workspace%basis, 2) < ncols) deallocate(workspace%basis)
            end if
            if (.not. allocated(workspace%basis)) allocate(workspace%basis(n, ncols))
        end if
        if (allocated(workspace%candidate)) then
            if (size(workspace%candidate) /= n) deallocate(workspace%candidate)
        end if
        if (.not. allocated(workspace%candidate)) allocate(workspace%candidate(n))
        if (allocated(workspace%work)) then
            if (size(workspace%work) /= n) deallocate(workspace%work)
        end if
        if (.not. allocated(workspace%work)) allocate(workspace%work(n))
    end subroutine ensure_pruning_workspace

    subroutine reset_lgmres_image_pruning_stream(workspace)
        type(lgmres_pruning_workspace), intent(inout) :: workspace !! work arrays and incremental state

        ! Jacobian images from a previous nonlinear iterate are not valid at
        ! the current point. Do not preserve the old orthonormal image span.
        workspace%stream_count = 0
        workspace%stream_basis_count = 0
    end subroutine reset_lgmres_image_pruning_stream

    subroutine retain_lgmres_image_direction_stream(stored, count, direction, &
            direction_image, tolerance, added, workspace, auto_threshold, basis_work, stream_limit)
        real(dp), intent(inout) :: stored(:, :) !! retained solution directions, append-only during stream
        integer, intent(inout) :: count !! retained direction count (initially zero)
        real(dp), intent(in) :: direction(:) !! candidate solution direction
        real(dp), intent(in) :: direction_image(:) !! candidate image under current Jacobian
        real(dp), intent(in) :: tolerance !! fixed dependence threshold or zero in automatic mode
        logical, intent(out) :: added !! whether the candidate is independent and retained
        type(lgmres_pruning_workspace), intent(inout) :: workspace !! private stream state and scratch
        logical, intent(in), optional :: auto_threshold !! use the projection-dependent cutoff
        real(dp), intent(inout), optional :: basis_work(:, :) !! scratch reused from retired Arnoldi basis
        integer, intent(in), optional :: stream_limit !! maximum active columns in append-only pass

        integer :: n
        integer :: capacity
        logical :: automatic

        added = .false.
        n = size(direction)
        automatic = .false.
        if (present(auto_threshold)) automatic = auto_threshold
        if (n == 0 .or. size(direction_image) /= n) return
        if (size(stored, 1) /= n .or. size(stored, 2) == 0) return
        capacity = size(stored, 2)
        if (present(stream_limit)) capacity = stream_limit
        if (capacity < 1 .or. capacity > size(stored, 2)) return
        ! A stream starts with count=0, only appends, and must be reset for
        ! each changed Jacobian. FIFO eviction requires a fresh QR factor.
        if (count /= workspace%stream_count .or. count >= capacity) return
        if (.not. ieee_is_finite(tolerance)) return
        if (automatic) then
            if (tolerance /= 0.0_dp) return
        else
            if (tolerance <= 0.0_dp .or. tolerance >= 1.0_dp) return
        end if
        if (.not. all(ieee_is_finite(direction))) return
        if (.not. all(ieee_is_finite(direction_image))) return
        if (maxval(abs(direction)) == 0.0_dp) return
        if (maxval(abs(direction_image)) == 0.0_dp) return

        ! Reserve space for the whole append-only stream once. Never grow or
        ! switch the backing basis midway through a sequence of candidates.
        if (present(basis_work)) then
            if (size(basis_work, 1) == n .and. &
                    size(basis_work, 2) >= capacity) then
                call ensure_pruning_workspace(workspace, n, 0, need_basis=.false.)
                call retain_image_stream_core(stored, count, direction, direction_image, &
                    tolerance, added, automatic, basis_work, workspace)
                return
            end if
        end if
        call ensure_pruning_workspace(workspace, n, capacity)
        call retain_image_stream_core(stored, count, direction, direction_image, &
            tolerance, added, automatic, workspace%basis, workspace)
    end subroutine retain_lgmres_image_direction_stream

    subroutine retain_image_stream_core(stored, count, direction, direction_image, &
            tolerance, added, automatic, basis, workspace)
        real(dp), intent(inout) :: stored(:, :) !! append-only retained solution directions
        integer, intent(inout) :: count !! retained count
        real(dp), intent(in) :: direction(:) !! candidate solution direction
        real(dp), intent(in) :: direction_image(:) !! current-Jacobian candidate image
        real(dp), intent(in) :: tolerance !! relative cutoff
        logical, intent(out) :: added !! retained or rejected
        logical, intent(in) :: automatic !! projection-dependent cutoff
        real(dp), intent(inout) :: basis(:, :) !! incrementally maintained image QR basis
        type(lgmres_pruning_workspace), intent(inout) :: workspace !! reusable scratch and state

        real(dp) :: scale
        real(dp) :: coefficient
        real(dp) :: distance
        real(dp) :: threshold
        integer :: i
        integer :: pass
        integer :: nbasis

        added = .false.
        nbasis = workspace%stream_basis_count
        scale = maxval(abs(direction_image))
        workspace%work = direction_image / scale
        workspace%work = workspace%work / norm2(workspace%work)
        threshold = tolerance
        do pass = 1, 2
            do i = 1, nbasis
                coefficient = dot_product(basis(:, i), workspace%work)
                workspace%work = workspace%work - coefficient * basis(:, i)
            end do
            if (automatic .and. pass == 1) &
                threshold = automatic_image_cutoff(norm2(workspace%work), size(direction))
        end do
        distance = norm2(workspace%work)
        if (distance <= threshold) return

        scale = maxval(abs(direction))
        workspace%candidate = direction / scale
        workspace%candidate = workspace%candidate / norm2(workspace%candidate)
        count = count + 1
        stored(:, count) = workspace%candidate
        workspace%stream_count = count
        added = .true.

        ! The ordinary span reconstruction skips previous images with
        ! orthogonal residual <= max(tolerance,32*epsilon). Apply that same
        ! floor when deciding which accepted images enlarge the basis.
        if (distance > max(threshold, 32.0_dp * epsilon(1.0_dp))) then
            nbasis = nbasis + 1
            basis(:, nbasis) = workspace%work / distance
            workspace%stream_basis_count = nbasis
        end if
    end subroutine retain_image_stream_core

    subroutine retain_lgmres_image_candidate_cached(stored, count, direction, &
            direction_image, tolerance, added, workspace, valid, auto_threshold, basis_work)
        ! Reuse the current-Jacobian image span built by the append-only
        ! filtering pass. Exactly one candidate is inserted after Arnoldi;
        ! a new nonlinear point must reset and rebuild the image basis.
        ! This remains valid at full capacity because we first test against
        ! the *old* retained span, then evict the oldest accepted column.
        real(dp), intent(inout) :: stored(:, :) !! retained solution directions
        integer, intent(inout) :: count !! active retained count
        real(dp), intent(in) :: direction(:) !! prospective Newton direction
        real(dp), intent(in) :: direction_image(:) !! J(x) times that direction
        real(dp), intent(in) :: tolerance !! fixed dependence cutoff, or zero for automatic
        logical, intent(out) :: added !! true if candidate was retained
        type(lgmres_pruning_workspace), intent(inout) :: workspace !! preserved current-Jacobian basis
        logical, intent(out) :: valid !! false if caller must use the reconstructing path
        logical, intent(in), optional :: auto_threshold !! automatic image cutoff
        real(dp), intent(in), optional :: basis_work(:, :) !! caller-owned current-Jacobian image basis

        integer :: n
        integer :: nbasis
        integer :: i
        integer :: pass
        real(dp) :: scale
        real(dp) :: coefficient
        real(dp) :: distance
        real(dp) :: threshold
        logical :: automatic

        added = .false.
        valid = .false.
        n = size(direction)
        automatic = .false.
        if (present(auto_threshold)) automatic = auto_threshold
        if (n == 0 .or. size(direction_image) /= n .or. &
                size(stored, 1) /= n .or. size(stored, 2) == 0) return
        if (count < 0 .or. count > size(stored, 2)) return
        if (workspace%stream_count /= count) return
        nbasis = workspace%stream_basis_count
        if (nbasis < 0 .or. nbasis > count) return
        if (present(basis_work)) then
            if (size(basis_work, 1) /= n .or. size(basis_work, 2) < nbasis) return
        else
            if (.not. allocated(workspace%basis)) return
            if (size(workspace%basis, 1) /= n .or. &
                    size(workspace%basis, 2) < nbasis) return
        end if
        if (.not. allocated(workspace%work)) return
        if (.not. allocated(workspace%candidate)) return
        if (size(workspace%work) /= n .or. size(workspace%candidate) /= n) return
        if (.not. ieee_is_finite(tolerance)) return
        if (automatic) then
            if (tolerance /= 0.0_dp) return
        else
            if (tolerance <= 0.0_dp .or. tolerance >= 1.0_dp) return
        end if
        if (.not. all(ieee_is_finite(direction))) return
        if (.not. all(ieee_is_finite(direction_image))) return
        if (count > 0) then
            if (.not. all(ieee_is_finite(stored(:, 1:count)))) return
        end if
        if (maxval(abs(direction)) == 0.0_dp .or. &
                maxval(abs(direction_image)) == 0.0_dp) then
            valid = .true.
            workspace%stream_count = -1
            return
        end if

        valid = .true.
        ! Invalidate the cached span after this one insertion. Even rejection
        ! must not permit another insertion without a fresh image pass.
        workspace%stream_count = -1
        scale = maxval(abs(direction_image))
        workspace%work = direction_image / scale
        workspace%work = workspace%work / norm2(workspace%work)
        threshold = tolerance
        do pass = 1, 2
            if (present(basis_work)) then
                do i = 1, nbasis
                    coefficient = dot_product(basis_work(:, i), workspace%work)
                    workspace%work = workspace%work - coefficient * basis_work(:, i)
                end do
            else
                do i = 1, nbasis
                    coefficient = dot_product(workspace%basis(:, i), workspace%work)
                    workspace%work = workspace%work - coefficient * workspace%basis(:, i)
                end do
            end if
            if (automatic .and. pass == 1) &
                threshold = automatic_image_cutoff(norm2(workspace%work), n)
        end do
        distance = norm2(workspace%work)
        if (distance <= threshold) return

        scale = maxval(abs(direction))
        workspace%candidate = direction / scale
        workspace%candidate = workspace%candidate / norm2(workspace%candidate)
        if (count < size(stored, 2)) then
            count = count + 1
        else if (count > 1) then
            stored(:, 1:count - 1) = stored(:, 2:count)
        end if
        stored(:, count) = workspace%candidate
        added = .true.
    end subroutine retain_lgmres_image_candidate_cached

    subroutine reset_lgmres_direction_pruning_stream(workspace)
        type(lgmres_pruning_workspace), intent(inout) :: workspace !! private basis and stream state

        ! After FIFO eviction (or a change of retained directions), rebuild
        ! instead of using an orthogonal span with the evicted vector in it.
        workspace%direction_stream_count = 0
        workspace%direction_basis_count = 0
        workspace%direction_fifo_ready = .false.
        workspace%direction_fifo_disabled = .false.
        workspace%direction_fifo_updates = 0
        workspace%direction_fifo_count = 0
    end subroutine reset_lgmres_direction_pruning_stream

    subroutine retain_lgmres_direction_stream(stored, count, direction, &
            tolerance, added, workspace)
        real(dp), intent(inout) :: stored(:, :) !! append-only retained solution directions
        integer, intent(inout) :: count !! current retained direction count
        real(dp), intent(in) :: direction(:) !! proposed unnormalized direction
        real(dp), intent(in) :: tolerance !! relative distance cutoff in (0, 1)
        logical, intent(out) :: added !! true if accepted and appended
        type(lgmres_pruning_workspace), intent(inout) :: workspace !! persistent orthogonal span

        integer :: n
        integer :: nbasis
        integer :: i
        integer :: pass
        real(dp) :: scale
        real(dp) :: coefficient
        real(dp) :: distance

        added = .false.
        n = size(direction)
        if (n == 0 .or. size(stored, 1) /= n .or. size(stored, 2) == 0) return
        if (count < 0 .or. count >= size(stored, 2)) return
        if (count /= workspace%direction_stream_count) return
        if (.not. ieee_is_finite(tolerance)) return
        if (tolerance <= 0.0_dp .or. tolerance >= 1.0_dp) return
        if (.not. all(ieee_is_finite(direction))) return
        scale = maxval(abs(direction))
        if (scale == 0.0_dp) return

        ! Reserve the complete append-only history so a resize never loses
        ! the existing QR basis. The caller can limit the stream capacity.
        call ensure_pruning_workspace(workspace, n, size(stored, 2))
        nbasis = workspace%direction_basis_count
        workspace%candidate = direction / scale
        workspace%candidate = workspace%candidate / norm2(workspace%candidate)
        workspace%work = workspace%candidate
        do pass = 1, 2
            do i = 1, nbasis
                coefficient = dot_product(workspace%basis(:, i), workspace%work)
                workspace%work = workspace%work - coefficient * workspace%basis(:, i)
            end do
        end do
        distance = norm2(workspace%work)
        if (distance <= tolerance) return

        count = count + 1
        stored(:, count) = workspace%candidate
        workspace%direction_stream_count = count
        added = .true.

        ! Reproduce the legacy reconstruction rule: normalize the *stored*
        ! column again before updating the basis, and skip columns below the
        ! max(tolerance, 32*epsilon) floor. A tiny accepted direction may not
        ! enlarge the numerical QR span, and this distinction matters for
        ! subsequent candidates when tolerance is below the roundoff floor.
        workspace%work = stored(:, count)
        scale = maxval(abs(workspace%work))
        workspace%work = workspace%work / scale
        workspace%work = workspace%work / norm2(workspace%work)
        do pass = 1, 2
            do i = 1, nbasis
                coefficient = dot_product(workspace%basis(:, i), workspace%work)
                workspace%work = workspace%work - coefficient * workspace%basis(:, i)
            end do
        end do
        distance = norm2(workspace%work)
        if (distance > max(tolerance, 32.0_dp * epsilon(1.0_dp))) then
            nbasis = nbasis + 1
            workspace%basis(:, nbasis) = workspace%work / distance
            workspace%direction_basis_count = nbasis
        end if
    end subroutine retain_lgmres_direction_stream

    subroutine retain_lgmres_direction_fifo(stored, count, direction, tolerance, &
            added, workspace)
        real(dp), intent(inout) :: stored(:, :) !! retained normalized directions, oldest first
        integer, intent(inout) :: count !! number of active columns, including FIFO history
        real(dp), intent(in) :: direction(:) !! proposed unnormalized direction
        real(dp), intent(in) :: tolerance !! relative dependence threshold in (0,1)
        logical, intent(out) :: added !! accepted candidate
        type(lgmres_pruning_workspace), intent(inout) :: workspace !! solver-local QR state

        integer :: n
        integer :: k
        integer :: i
        integer :: j
        integer :: pass
        real(dp) :: scale
        real(dp) :: distance
        real(dp) :: coefficient
        real(dp) :: a
        real(dp) :: b
        real(dp) :: c
        real(dp) :: sn
        real(dp) :: rotation_norm
        logical :: valid

        added = .false.
        n = size(direction)
        k = size(stored, 2)
        if (n == 0 .or. size(stored, 1) /= n .or. k < 2) return
        if (count /= k) return
        if (.not. ieee_is_finite(tolerance)) return
        if (tolerance <= 0.0_dp .or. tolerance >= 1.0_dp) return
        if (.not. all(ieee_is_finite(direction))) return
        if (.not. all(ieee_is_finite(stored(:, 1:k)))) return
        scale = maxval(abs(direction))
        if (scale == 0.0_dp) return

        ! A direct QR downdate is valid only while every retained direction
        ! enlarges a well-resolved, full-rank span. Near dependence uses the
        ! existing reconstruction path instead of amplifying roundoff.
        if (workspace%direction_fifo_disabled) then
            call retain_lgmres_direction(stored, count, direction, tolerance, added, workspace)
            return
        end if
        if (workspace%direction_fifo_ready) then
            if (workspace%direction_fifo_count /= k) workspace%direction_fifo_ready = .false.
        end if
        if (.not. workspace%direction_fifo_ready) then
            call ensure_pruning_workspace(workspace, n, k)
            if (allocated(workspace%direction_r)) then
                if (size(workspace%direction_r, 1) /= k .or. &
                        size(workspace%direction_r, 2) /= k) deallocate(workspace%direction_r)
            end if
            if (.not. allocated(workspace%direction_r)) allocate(workspace%direction_r(k, k))
            call build_direction_fifo_qr(stored(:, 1:k), tolerance, workspace, valid)
            if (.not. valid) then
                workspace%direction_fifo_disabled = .true.
                call retain_lgmres_direction(stored, count, direction, tolerance, added, workspace)
                return
            end if
            workspace%direction_fifo_ready = .true.
            workspace%direction_fifo_count = k
            workspace%direction_fifo_updates = 0
        end if

        workspace%candidate = direction / scale
        workspace%candidate = workspace%candidate / norm2(workspace%candidate)
        workspace%work = workspace%candidate
        do pass = 1, 2
            do i = 1, k
                coefficient = dot_product(workspace%basis(:, i), workspace%work)
                workspace%work = workspace%work - coefficient * workspace%basis(:, i)
            end do
        end do
        distance = norm2(workspace%work)
        if (distance <= tolerance) return

        ! Remove the first column of R, then chase the subdiagonal of the
        ! resulting k x (k-1) upper-Hessenberg factor with Givens rotations.
        ! Apply the transposed rotations to Q so that Q*R is unchanged.
        do j = 1, k - 1
            workspace%direction_r(:, j) = workspace%direction_r(:, j + 1)
        end do
        do i = 1, k - 1
            a = workspace%direction_r(i, i)
            b = workspace%direction_r(i + 1, i)
            rotation_norm = norm2([a, b])
            if (rotation_norm <= 64.0_dp * epsilon(1.0_dp)) then
                workspace%direction_fifo_disabled = .true.
                workspace%direction_fifo_ready = .false.
                call retain_lgmres_direction(stored, count, direction, tolerance, added, workspace)
                return
            end if
            c = a / rotation_norm
            sn = b / rotation_norm
            do j = i, k - 1
                a = workspace%direction_r(i, j)
                b = workspace%direction_r(i + 1, j)
                workspace%direction_r(i, j) = c * a + sn * b
                workspace%direction_r(i + 1, j) = -sn * a + c * b
            end do
            workspace%work = c * workspace%basis(:, i) + sn * workspace%basis(:, i + 1)
            workspace%basis(:, i + 1) = -sn * workspace%basis(:, i) + c * workspace%basis(:, i + 1)
            workspace%basis(:, i) = workspace%work
        end do

        ! Append the accepted normalized direction to the remaining k-1
        ! columns. This is O(n*k), versus O(n*k*k) to rebuild the QR span.
        workspace%work = workspace%candidate
        workspace%direction_r(:, k) = 0.0_dp
        do pass = 1, 2
            do i = 1, k - 1
                coefficient = dot_product(workspace%basis(:, i), workspace%work)
                workspace%direction_r(i, k) = workspace%direction_r(i, k) + coefficient
                workspace%work = workspace%work - coefficient * workspace%basis(:, i)
            end do
        end do
        distance = norm2(workspace%work)
        if (distance <= max(tolerance, 64.0_dp * epsilon(1.0_dp))) then
            ! The FIFO deletion exposed nearly dependent retained columns.
            ! Rebuild through the original algorithm for this and later steps.
            workspace%direction_fifo_disabled = .true.
            workspace%direction_fifo_ready = .false.
            call retain_lgmres_direction(stored, count, direction, tolerance, added, workspace)
            return
        end if
        workspace%direction_r(k, k) = distance
        workspace%basis(:, k) = workspace%work / distance
        stored(:, 1:k-1) = stored(:, 2:k)
        stored(:, k) = workspace%candidate
        added = .true.
        workspace%direction_fifo_updates = workspace%direction_fifo_updates + 1
        if (workspace%direction_fifo_updates >= 32) workspace%direction_fifo_ready = .false.
    end subroutine retain_lgmres_direction_fifo

    subroutine build_direction_fifo_qr(stored, tolerance, workspace, valid)
        real(dp), intent(in) :: stored(:, :) !! current complete augmentation history
        real(dp), intent(in) :: tolerance !! same numerical cutoff as normal pruning
        type(lgmres_pruning_workspace), intent(inout) :: workspace !! preallocated Q, R, and scratch
        logical, intent(out) :: valid !! full-rank, numerically reliable factorization

        integer :: i
        integer :: j
        integer :: pass
        real(dp) :: scale
        real(dp) :: coefficient
        real(dp) :: distance

        valid = .false.
        workspace%direction_r = 0.0_dp
        do j = 1, size(stored, 2)
            scale = maxval(abs(stored(:, j)))
            if (scale == 0.0_dp) return
            workspace%work = stored(:, j) / scale
            workspace%work = workspace%work / norm2(workspace%work)
            do pass = 1, 2
                do i = 1, j - 1
                    coefficient = dot_product(workspace%basis(:, i), workspace%work)
                    workspace%direction_r(i, j) = workspace%direction_r(i, j) + coefficient
                    workspace%work = workspace%work - coefficient * workspace%basis(:, i)
                end do
            end do
            distance = norm2(workspace%work)
            if (distance <= max(tolerance, 64.0_dp * epsilon(1.0_dp))) return
            workspace%direction_r(j, j) = distance
            workspace%basis(:, j) = workspace%work / distance
        end do
        valid = .true.
    end subroutine build_direction_fifo_qr

    subroutine retain_lgmres_direction(stored, count, direction, tolerance, added, &
            workspace, basis_work)
        real(dp), intent(inout) :: stored(:, :) !! retained normalized solution directions
        integer, intent(inout) :: count !! number of active retained columns
        real(dp), intent(in) :: direction(:) !! new unnormalized solution direction
        real(dp), intent(in) :: tolerance !! maximum relative distance from retained span
        logical, intent(out) :: added !! whether the direction was retained
        type(lgmres_pruning_workspace), intent(inout), optional :: workspace !! reusable MGS scratch
        real(dp), intent(inout), optional :: basis_work(:, :) !! caller-owned Arnoldi scratch basis

        real(dp), allocatable :: basis(:, :)
        real(dp), allocatable :: candidate(:)
        real(dp), allocatable :: work(:)
        real(dp) :: scale
        integer :: n

        added = .false.
        n = size(direction)
        if (n == 0 .or. size(stored, 1) /= n .or. size(stored, 2) == 0) return
        if (count < 0 .or. count > size(stored, 2)) return
        if (.not. ieee_is_finite(tolerance)) return
        if (tolerance <= 0.0_dp .or. tolerance >= 1.0_dp) return
        if (.not. all(ieee_is_finite(direction))) return
        if (count > 0) then
            if (.not. all(ieee_is_finite(stored(:, 1:count)))) return
        end if

        ! Scale before taking a norm, so huge and subnormal but finite
        ! Newton steps can be retained without squaring their magnitudes.
        scale = maxval(abs(direction))
        if (scale == 0.0_dp) return
        if (present(workspace)) then
            if (present(basis_work)) then
                if (size(basis_work, 1) == n .and. size(basis_work, 2) >= count) then
                    call ensure_pruning_workspace(workspace, n, count, need_basis=.false.)
                    call retain_direction_core(stored, count, direction, tolerance, added, &
                        basis_work, workspace%candidate, workspace%work)
                    return
                end if
            end if
            call ensure_pruning_workspace(workspace, n, count)
            call retain_direction_core(stored, count, direction, tolerance, added, &
                workspace%basis, workspace%candidate, workspace%work)
        else
            allocate(candidate(n), work(n), basis(n, count))
            call retain_direction_core(stored, count, direction, tolerance, added, &
                basis, candidate, work)
        end if
    end subroutine retain_lgmres_direction

    subroutine retain_direction_core(stored, count, direction, tolerance, added, &
            basis, candidate, work)
        real(dp), intent(inout) :: stored(:, :) !! retained vectors
        integer, intent(inout) :: count !! active retained vector count
        real(dp), intent(in) :: direction(:) !! new direction
        real(dp), intent(in) :: tolerance !! numerical dependence threshold
        logical, intent(out) :: added !! whether appended
        real(dp), intent(inout) :: basis(:, :) !! orthonormal basis scratch
        real(dp), intent(inout) :: candidate(:) !! candidate scratch
        real(dp), intent(inout) :: work(:) !! projection scratch
        real(dp) :: coefficient
        real(dp) :: distance
        real(dp) :: scale
        integer :: i
        integer :: j
        integer :: nbasis
        integer :: pass

        scale = maxval(abs(direction))
        candidate = direction / scale
        candidate = candidate / norm2(candidate)

        ! Construct a *temporary* orthonormal basis of the old directions.
        ! The original stored vectors must remain unchanged: changing them
        ! would change the augmentation operator in subsequent GMRES solves.
        nbasis = 0
        do j = 1, count
            work = stored(:, j)
            scale = maxval(abs(work))
            if (scale == 0.0_dp) cycle
            work = work / scale
            work = work / norm2(work)
            do pass = 1, 2
                do i = 1, nbasis
                    coefficient = dot_product(basis(:, i), work)
                    work = work - coefficient * basis(:, i)
                end do
            end do
            distance = norm2(work)
            if (distance <= max(tolerance, 32.0_dp * epsilon(1.0_dp))) cycle
            nbasis = nbasis + 1
            basis(:, nbasis) = work / distance
        end do

        work = candidate
        do pass = 1, 2
            do i = 1, nbasis
                coefficient = dot_product(basis(:, i), work)
                work = work - coefficient * basis(:, i)
            end do
        end do
        ! Include equality: when tolerance is met, the new direction does
        ! not add a numerically distinguishable component to the old span.
        if (norm2(work) <= tolerance) return

        if (count < size(stored, 2)) then
            count = count + 1
        else if (count > 1) then
            stored(:, 1:count-1) = stored(:, 2:count)
        end if
        stored(:, count) = candidate
        added = .true.
    end subroutine retain_direction_core


    subroutine retain_lgmres_image_direction(stored, count, direction, images, &
            direction_image, tolerance, added, auto_threshold, workspace, basis_work)
        real(dp), intent(inout) :: stored(:, :) !! retained normalized solution directions
        integer, intent(inout) :: count !! active columns in stored and images
        real(dp), intent(in) :: direction(:) !! new unnormalized solution direction
        real(dp), intent(in) :: images(:, :) !! images of retained columns at the current Jacobian
        real(dp), intent(in) :: direction_image(:) !! image of the proposed direction at the same Jacobian
        real(dp), intent(in) :: tolerance !! relative Jacobian-image dependence threshold in (0, 1)
        logical, intent(out) :: added !! true if the direction was retained
        logical, intent(in), optional :: auto_threshold !! derive cutoff from image projections
        type(lgmres_pruning_workspace), intent(inout), optional :: workspace !! reusable MGS scratch
        real(dp), intent(inout), optional :: basis_work(:, :) !! caller-owned Arnoldi scratch basis

        real(dp), allocatable :: basis(:, :)
        real(dp), allocatable :: candidate(:)
        real(dp), allocatable :: work(:)
        real(dp) :: scale
        integer :: n
        logical :: auto_use

        added = .false.
        auto_use = .false.
        if (present(auto_threshold)) auto_use = auto_threshold
        n = size(direction)
        if (n == 0 .or. size(stored, 1) /= n .or. size(stored, 2) == 0) return
        if (size(images, 1) /= n .or. size(direction_image) /= n) return
        if (count < 0 .or. count > min(size(stored, 2), size(images, 2))) return
        if (.not. ieee_is_finite(tolerance)) return
        if (auto_use) then
            if (tolerance /= 0.0_dp) return
        else
            if (tolerance <= 0.0_dp .or. tolerance >= 1.0_dp) return
        end if
        if (.not. all(ieee_is_finite(direction))) return
        if (.not. all(ieee_is_finite(direction_image))) return
        if (count > 0) then
            if (.not. all(ieee_is_finite(stored(:, 1:count)))) return
            if (.not. all(ieee_is_finite(images(:, 1:count)))) return
        end if

        ! Normalize independently in solution and Jacobian-image spaces.
        ! Filtering uses the images only, but the retained vector must stay
        ! in the original solution coordinates for the Arnoldi operator.
        scale = maxval(abs(direction))
        if (scale == 0.0_dp) return
        if (present(workspace)) then
            if (present(basis_work)) then
                if (size(basis_work, 1) == n .and. size(basis_work, 2) >= count) then
                    call ensure_pruning_workspace(workspace, n, count, need_basis=.false.)
                    call retain_image_core(stored, count, direction, images, direction_image, &
                        tolerance, added, auto_use, basis_work, workspace%candidate, &
                        workspace%work)
                    return
                end if
            end if
            call ensure_pruning_workspace(workspace, n, count)
            call retain_image_core(stored, count, direction, images, direction_image, &
                tolerance, added, auto_use, workspace%basis, workspace%candidate, &
                workspace%work)
        else
            allocate(candidate(n), work(n), basis(n, count))
            call retain_image_core(stored, count, direction, images, direction_image, &
                tolerance, added, auto_use, basis, candidate, work)
        end if
    end subroutine retain_lgmres_image_direction

    subroutine retain_image_core(stored, count, direction, images, direction_image, &
            tolerance, added, auto_use, basis, candidate, work)
        real(dp), intent(inout) :: stored(:, :) !! retained solution directions
        integer, intent(inout) :: count !! active retained vector count
        real(dp), intent(in) :: direction(:) !! proposed direction
        real(dp), intent(in) :: images(:, :) !! Jacobian images of retained directions
        real(dp), intent(in) :: direction_image(:) !! Jacobian image of proposal
        real(dp), intent(in) :: tolerance !! numerical dependence threshold
        logical, intent(out) :: added !! whether appended
        logical, intent(in) :: auto_use !! adaptive numerical cutoff enabled
        real(dp), intent(inout) :: basis(:, :) !! orthonormal image basis scratch
        real(dp), intent(inout) :: candidate(:) !! candidate scratch
        real(dp), intent(inout) :: work(:) !! projection scratch
        real(dp) :: coefficient
        real(dp) :: distance
        real(dp) :: scale
        real(dp) :: threshold
        integer :: i
        integer :: j
        integer :: n
        integer :: nbasis
        integer :: pass

        n = size(direction)
        scale = maxval(abs(direction))
        candidate = direction / scale
        candidate = candidate / norm2(candidate)

        scale = maxval(abs(direction_image))
        if (scale == 0.0_dp) return
        work = direction_image / scale
        work = work / norm2(work)
        nbasis = 0
        do j = 1, count
            scale = maxval(abs(images(:, j)))
            if (scale == 0.0_dp) cycle
            candidate = images(:, j) / scale
            candidate = candidate / norm2(candidate)
            threshold = tolerance
            do pass = 1, 2
                do i = 1, nbasis
                    coefficient = dot_product(basis(:, i), candidate)
                    candidate = candidate - coefficient * basis(:, i)
                end do
                if (auto_use .and. pass == 1) then
                    threshold = automatic_image_cutoff(norm2(candidate), n)
                end if
            end do
            distance = norm2(candidate)
            if (distance <= max(threshold, 32.0_dp * epsilon(1.0_dp))) cycle
            nbasis = nbasis + 1
            basis(:, nbasis) = candidate / distance
        end do

        threshold = tolerance
        do pass = 1, 2
            do i = 1, nbasis
                coefficient = dot_product(basis(:, i), work)
                work = work - coefficient * basis(:, i)
            end do
            if (auto_use .and. pass == 1) then
                threshold = automatic_image_cutoff(norm2(work), n)
            end if
        end do
        if (norm2(work) <= threshold) return

        ! Preserve the solution vector (not its image), and FIFO history.
        scale = maxval(abs(direction))
        candidate = direction / scale
        candidate = candidate / norm2(candidate)
        if (count < size(stored, 2)) then
            count = count + 1
        else if (count > 1) then
            stored(:, 1:count - 1) = stored(:, 2:count)
        end if
        stored(:, count) = candidate
        added = .true.
    end subroutine retain_image_core

    pure function automatic_image_cutoff(projected_norm, dimension) result(threshold)
        real(dp), intent(in) :: projected_norm !! normalized image norm after first MGS pass
        integer, intent(in) :: dimension !! number of image components
        real(dp) :: threshold !! automatic relative numerical-dependence cutoff

        real(dp) :: alignment
        real(dp) :: remaining
        real(dp) :: roundoff_floor

        ! The roundoff floor grows slowly with dimension. The additional
        ! sqrt(epsilon) term reflects the *observed* amount projected away:
        ! nearly parallel images warrant a larger cutoff than almost
        ! orthogonal images. This is a heuristic, not an error estimator for
        ! an arbitrary user-supplied (possibly noisy) Jacobian product.
        remaining = min(1.0_dp, max(0.0_dp, projected_norm))
        alignment = sqrt(max(0.0_dp, (1.0_dp - remaining) * (1.0_dp + remaining)))
        roundoff_floor = 32.0_dp * epsilon(1.0_dp) * &
            sqrt(real(max(1, dimension), dp))
        threshold = min(0.25_dp, max(roundoff_floor, &
            sqrt(epsilon(1.0_dp)) * alignment))
    end function automatic_image_cutoff

end module scifort_optimize_krylov_augmentation
