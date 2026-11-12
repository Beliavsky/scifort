! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Continuity tracking for conditional canonical spectral directions.
! Original SciFort extension; scipy.signal does not expose this operation.
! Match successive isolated, nonzero canonical directions with an optimal
! one-to-one assignment, then align the shared complex phase of each pair.
module scifort_signal_canonical_tracking
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_canonical_coherence, only : canonical_coherence_result, &
        canonical_coherence_tf_result
    implicit none
    private

    integer, parameter, public :: canonical_tracking_success = 0
    integer, parameter, public :: canonical_tracking_invalid_input = 1

    type, public :: canonical_mode_tracks_result
        real(dp), allocatable :: frequencies(:) !! frequency coordinates, unchanged
        real(dp), allocatable :: squared_modes(:,:) !! (frequency,track), not necessarily descending
        complex(dp), allocatable :: coefficients_a(:,:,:) !! aligned (frequency,A-channel,track)
        complex(dp), allocatable :: coefficients_b(:,:,:) !! aligned (frequency,B-channel,track)
        logical, allocatable :: directions_stable(:,:) !! isolated nonzero root flags, reordered
        integer, allocatable :: source_mode(:,:) !! 1-based descending root index, 0 at invalid bins
        real(dp), allocatable :: overlap(:,:) !! phase-invariant matched similarity, NaN at restarts
        logical, allocatable :: matched(:,:) !! successfully linked to previous frequency bin
        logical, allocatable :: segment_start(:) !! valid bin at which tracking restarts
        logical, allocatable :: valid(:) !! input spectral validity mask
        integer :: source_status = 0 !! status from the original canonical estimate
        integer :: status = canonical_tracking_invalid_input !! tracking input validation
    end type canonical_mode_tracks_result

    type, public :: canonical_mode_tracks_tf_result
        real(dp), allocatable :: frequencies(:) !! frequency coordinates, unchanged
        real(dp), allocatable :: times(:) !! time coordinates, unchanged
        real(dp), allocatable :: squared_modes(:,:,:) !! (frequency,time,track)
        complex(dp), allocatable :: coefficients_a(:,:,:,:) !! aligned (frequency,time,A,track)
        complex(dp), allocatable :: coefficients_b(:,:,:,:) !! aligned (frequency,time,B,track)
        logical, allocatable :: directions_stable(:,:,:) !! reordered isolated-root flags
        integer, allocatable :: source_mode(:,:,:) !! original sorted root index, 0 if invalid
        real(dp), allocatable :: overlap(:,:,:) !! similarity to predecessor; NaN at restarts
        logical, allocatable :: matched(:,:,:) !! successfully linked to predecessor
        logical, allocatable :: segment_start(:,:) !! (frequency,time), independent per time
        logical, allocatable :: valid(:,:) !! input spectral validity mask
        integer :: source_status = 0 !! original canonical estimate status
        integer :: status = canonical_tracking_invalid_input !! tracking input validation
    end type canonical_mode_tracks_tf_result

    interface track_canonical_modes
        module procedure track_canonical_modes_static
        module procedure track_canonical_modes_tf
    end interface track_canonical_modes
    public :: track_canonical_modes

contains

    ! The Hungarian shortest augmenting path algorithm finds the globally
    ! optimal assignment in O(n**3), avoiding greedy label swaps.
    pure subroutine maximum_assignment(scores,assignment)
        real(dp), intent(in) :: scores(:,:) !! finite similarity scores in [0,1], square
        integer, intent(out) :: assignment(:) !! selected column for each row, 1-based
        real(dp) :: u(0:size(scores,1)), v(0:size(scores,1))
        real(dp) :: minima(0:size(scores,1)), change, current
        logical :: used(0:size(scores,1))
        integer :: path(0:size(scores,1)), owner(0:size(scores,1))
        integer :: i, j, j0, j1, row, n

        n = size(scores,1)
        assignment = 0
        u = 0.0_dp
        v = 0.0_dp
        path = 0
        owner = 0
        do i = 1, n
            owner(0) = i
            j0 = 0
            used = .false.
            minima = huge(1.0_dp)
            do
                used(j0) = .true.
                row = owner(j0)
                change = huge(1.0_dp)
                j1 = 0
                do j = 1, n
                    if (used(j)) cycle
                    current = 1.0_dp-scores(row,j)-u(row)-v(j)
                    if (current < minima(j)) then
                        minima(j) = current
                        path(j) = j0
                    end if
                    if (minima(j) < change) then
                        change = minima(j)
                        j1 = j
                    end if
                end do
                do j = 0, n
                    if (used(j)) then
                        u(owner(j)) = u(owner(j))+change
                        v(j) = v(j)-change
                    else
                        minima(j) = minima(j)-change
                    end if
                end do
                j0 = j1
                if (owner(j0) == 0) exit
            end do
            do
                j1 = path(j0)
                owner(j0) = owner(j1)
                j0 = j1
                if (j0 == 0) exit
            end do
        end do
        do j = 1, n
            assignment(owner(j)) = j
        end do
    end subroutine maximum_assignment

    pure logical function usable_directions(a,b) result(ok)
        complex(dp), intent(in) :: a(:,:) !! A coefficients (channel,mode)
        complex(dp), intent(in) :: b(:,:) !! B coefficients (channel,mode)
        integer :: k

        ok = .false.
        if (.not. all(ieee_is_finite(real(a,dp)))) return
        if (.not. all(ieee_is_finite(aimag(a)))) return
        if (.not. all(ieee_is_finite(real(b,dp)))) return
        if (.not. all(ieee_is_finite(aimag(b)))) return
        do k = 1, size(a,2)
            if (maxval(abs(a(:,k))) <= tiny(1.0_dp)) return
            if (maxval(abs(b(:,k))) <= tiny(1.0_dp)) return
        end do
        ok = .true.
    end function usable_directions

    ! Both groups must receive the same phase to retain the canonical
    ! cross-spectrum convention. Normalize each Euclidean inner product
    ! separately to avoid artificially favoring one group by its units.
    pure subroutine pair_similarity(prev_a,prev_b,new_a,new_b,score,phase,ok)
        complex(dp), intent(in) :: prev_a(:) !! preceding A coefficient vector
        complex(dp), intent(in) :: prev_b(:) !! preceding B coefficient vector
        complex(dp), intent(in) :: new_a(:) !! current A coefficient vector
        complex(dp), intent(in) :: new_b(:) !! current B coefficient vector
        real(dp), intent(out) :: score !! magnitude of mean normalized overlap
        complex(dp), intent(out) :: phase !! common multiplier making mean overlap real
        logical, intent(out) :: ok !! finite nonzero Euclidean norms and overlap
        real(dp) :: na, nb, nc, nd
        complex(dp) :: overlap

        score = 0.0_dp
        phase = (1.0_dp,0.0_dp)
        ok = .false.
        na = sqrt(sum(abs(prev_a)**2))
        nb = sqrt(sum(abs(prev_b)**2))
        nc = sqrt(sum(abs(new_a)**2))
        nd = sqrt(sum(abs(new_b)**2))
        if (.not. all(ieee_is_finite([na,nb,nc,nd]))) return
        if (min(na,nb,nc,nd) <= tiny(1.0_dp)) return
        overlap = 0.5_dp*(sum(conjg(prev_a/na)*(new_a/nc)) + &
                          sum(conjg(prev_b/nb)*(new_b/nd)))
        if (.not. ieee_is_finite(real(overlap,dp))) return
        if (.not. ieee_is_finite(aimag(overlap))) return
        score = min(1.0_dp,abs(overlap))
        if (score <= 64.0_dp*epsilon(1.0_dp)) then
            ok = .true.
            return
        end if
        phase = conjg(overlap)/abs(overlap)
        ok = .true.
    end subroutine pair_similarity

    pure subroutine track_sequence(freq,values,a,b,stable,valid,min_overlap, &
            tracked_values,tracked_a,tracked_b,tracked_stable,source_mode,overlap,matched,starts)
        real(dp), intent(in) :: freq(:) !! frequency coordinates; only increasing steps link
        real(dp), intent(in) :: values(:,:) !! (frequency,sorted mode), nonnegative
        complex(dp), intent(in) :: a(:,:,:) !! (frequency,A-channel,sorted mode)
        complex(dp), intent(in) :: b(:,:,:) !! (frequency,B-channel,sorted mode)
        logical, intent(in) :: stable(:,:) !! isolated nonzero modes; all required to link
        logical, intent(in) :: valid(:) !! valid spectral cells
        real(dp), intent(in) :: min_overlap !! acceptance floor, [0,1]
        real(dp), intent(out) :: tracked_values(:,:) !! values reordered by continuity
        complex(dp), intent(out) :: tracked_a(:,:,:) !! reordered and phase-aligned A vectors
        complex(dp), intent(out) :: tracked_b(:,:,:) !! reordered and phase-aligned B vectors
        logical, intent(out) :: tracked_stable(:,:) !! stability flags in tracked order
        integer, intent(out) :: source_mode(:,:) !! sorted index at each bin, zero if invalid
        real(dp), intent(out) :: overlap(:,:) !! matched similarities, NaN at restarts
        logical, intent(out) :: matched(:,:) !! row matched from preceding bin
        logical, intent(out) :: starts(:) !! track boundary for each valid bin
        real(dp) :: scores(size(values,2),size(values,2)), nan
        complex(dp) :: phases(size(values,2),size(values,2)), gauge
        integer :: n, nf, f, i, j, permutation(size(values,2))
        logical :: good, active, ready, pair_ok

        nf = size(values,1)
        n = size(values,2)
        nan = ieee_value(0.0_dp,ieee_quiet_nan)
        tracked_values = values
        tracked_a = a
        tracked_b = b
        tracked_stable = stable
        source_mode = 0
        overlap = nan
        matched = .false.
        starts = .false.
        active = .false.
        do f = 1, nf
            if (.not. valid(f)) then
                active = .false.
                cycle
            end if
            do i = 1, n
                source_mode(f,i) = i
            end do
            ready = all(stable(f,:)) .and. usable_directions(a(f,:,:),b(f,:,:))
            if (ready .and. active) then
                if (.not. ieee_is_finite(freq(f))) active = .false.
                if (.not. ieee_is_finite(freq(f-1))) active = .false.
                if (freq(f) <= freq(f-1)) active = .false.
            end if
            if (ready .and. active) then
                good = .true.
                do i = 1, n
                    do j = 1, n
                        call pair_similarity(tracked_a(f-1,:,i),tracked_b(f-1,:,i), &
                            a(f,:,j),b(f,:,j),scores(i,j),phases(i,j),pair_ok)
                        if (.not. pair_ok) good = .false.
                    end do
                end do
                if (good) then
                    call maximum_assignment(scores,permutation)
                    do i = 1, n
                        if (scores(i,permutation(i)) <= 64.0_dp*epsilon(1.0_dp)) good = .false.
                        if (scores(i,permutation(i)) < min_overlap) good = .false.
                    end do
                end if
                if (good) then
                    do i = 1, n
                        j = permutation(i)
                        gauge = phases(i,j)
                        source_mode(f,i) = j
                        overlap(f,i) = scores(i,j)
                        matched(f,i) = .true.
                        tracked_values(f,i) = values(f,j)
                        tracked_stable(f,i) = stable(f,j)
                        tracked_a(f,:,i) = gauge*a(f,:,j)
                        tracked_b(f,:,i) = gauge*b(f,:,j)
                    end do
                else
                    starts(f) = .true.
                end if
            else
                starts(f) = .true.
            end if
            active = ready
        end do
    end subroutine track_sequence

    function track_canonical_modes_static(source,min_overlap) result(out)
        type(canonical_coherence_result), intent(in) :: source !! computed canonical modes, with vectors
        real(dp), intent(in), optional :: min_overlap !! matching threshold [0,1], default 0.5
        type(canonical_mode_tracks_result) :: out
        integer :: nf, nm, na, nb
        real(dp) :: minimum

        minimum = 0.5_dp
        if (present(min_overlap)) minimum = min_overlap
        if (.not. ieee_is_finite(minimum)) return
        if (minimum < 0.0_dp .or. minimum > 1.0_dp) return
        if (.not. allocated(source%frequencies)) return
        if (.not. allocated(source%squared_modes)) return
        if (.not. allocated(source%coefficients_a)) return
        if (.not. allocated(source%coefficients_b)) return
        if (.not. allocated(source%directions_stable)) return
        if (.not. allocated(source%valid)) return
        nf = size(source%frequencies)
        nm = size(source%squared_modes,2)
        na = size(source%coefficients_a,2)
        nb = size(source%coefficients_b,2)
        if (nf < 1 .or. nm < 1 .or. na < 1 .or. nb < 1) return
        if (any(shape(source%squared_modes) /= [nf,nm])) return
        if (any(shape(source%coefficients_a) /= [nf,na,nm])) return
        if (any(shape(source%coefficients_b) /= [nf,nb,nm])) return
        if (any(shape(source%directions_stable) /= [nf,nm])) return
        if (size(source%valid) /= nf) return
        allocate(out%frequencies(nf),out%squared_modes(nf,nm), &
            out%coefficients_a(nf,na,nm),out%coefficients_b(nf,nb,nm), &
            out%directions_stable(nf,nm),out%source_mode(nf,nm), &
            out%overlap(nf,nm),out%matched(nf,nm),out%segment_start(nf),out%valid(nf))
        out%frequencies = source%frequencies
        out%valid = source%valid
        out%source_status = source%status
        call track_sequence(source%frequencies,source%squared_modes,source%coefficients_a, &
            source%coefficients_b,source%directions_stable,source%valid,minimum, &
            out%squared_modes,out%coefficients_a,out%coefficients_b, &
            out%directions_stable,out%source_mode,out%overlap,out%matched,out%segment_start)
        out%status = canonical_tracking_success
    end function track_canonical_modes_static

    function track_canonical_modes_tf(source,min_overlap) result(out)
        type(canonical_coherence_tf_result), intent(in) :: source !! spectral modes with vectors
        real(dp), intent(in), optional :: min_overlap !! matching floor [0,1], default 0.5
        type(canonical_mode_tracks_tf_result) :: out
        integer :: nf, nt, nm, na, nb, t
        real(dp) :: minimum

        minimum = 0.5_dp
        if (present(min_overlap)) minimum = min_overlap
        if (.not. ieee_is_finite(minimum)) return
        if (minimum < 0.0_dp .or. minimum > 1.0_dp) return
        if (.not. allocated(source%frequencies)) return
        if (.not. allocated(source%times)) return
        if (.not. allocated(source%squared_modes)) return
        if (.not. allocated(source%coefficients_a)) return
        if (.not. allocated(source%coefficients_b)) return
        if (.not. allocated(source%directions_stable)) return
        if (.not. allocated(source%valid)) return
        nf = size(source%frequencies)
        nt = size(source%times)
        nm = size(source%squared_modes,3)
        na = size(source%coefficients_a,3)
        nb = size(source%coefficients_b,3)
        if (nf < 1 .or. nt < 1 .or. nm < 1 .or. na < 1 .or. nb < 1) return
        if (any(shape(source%squared_modes) /= [nf,nt,nm])) return
        if (any(shape(source%coefficients_a) /= [nf,nt,na,nm])) return
        if (any(shape(source%coefficients_b) /= [nf,nt,nb,nm])) return
        if (any(shape(source%directions_stable) /= [nf,nt,nm])) return
        if (any(shape(source%valid) /= [nf,nt])) return
        allocate(out%frequencies(nf),out%times(nt),out%squared_modes(nf,nt,nm), &
            out%coefficients_a(nf,nt,na,nm),out%coefficients_b(nf,nt,nb,nm), &
            out%directions_stable(nf,nt,nm),out%source_mode(nf,nt,nm), &
            out%overlap(nf,nt,nm),out%matched(nf,nt,nm), &
            out%segment_start(nf,nt),out%valid(nf,nt))
        out%frequencies = source%frequencies
        out%times = source%times
        out%valid = source%valid
        out%source_status = source%status
        do t = 1, nt
            call track_sequence(source%frequencies,source%squared_modes(:,t,:), &
                source%coefficients_a(:,t,:,:),source%coefficients_b(:,t,:,:), &
                source%directions_stable(:,t,:),source%valid(:,t),minimum, &
                out%squared_modes(:,t,:),out%coefficients_a(:,t,:,:), &
                out%coefficients_b(:,t,:,:),out%directions_stable(:,t,:), &
                out%source_mode(:,t,:),out%overlap(:,t,:),out%matched(:,t,:), &
                out%segment_start(:,t))
        end do
        out%status = canonical_tracking_success
    end function track_canonical_modes_tf
end module scifort_signal_canonical_tracking
