! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! One-dimensional peak detection, prominence and width measurements.
! Independent implementation from geometric definitions of peak prominence.
module scifort_signal_peaks
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: peak_status_success = 0
    integer, parameter, public :: peak_status_invalid_input = 1

    type, public :: peak_result
        integer, allocatable :: peaks(:) !! one-based local-maximum locations
        real(dp), allocatable :: peak_heights(:) !! amplitudes when height was requested
        real(dp), allocatable :: left_thresholds(:) !! differences to left neighbors
        real(dp), allocatable :: right_thresholds(:) !! differences to right neighbors
        real(dp), allocatable :: prominences(:) !! prominence when requested or width used
        integer, allocatable :: left_bases(:) !! left prominence bases, one-based
        integer, allocatable :: right_bases(:) !! right prominence bases, one-based
        real(dp), allocatable :: widths(:) !! widths measured in sample intervals
        real(dp), allocatable :: width_heights(:) !! amplitude at the width evaluation level
        real(dp), allocatable :: left_ips(:) !! fractional one-based left intersections
        real(dp), allocatable :: right_ips(:) !! fractional one-based right intersections
        integer, allocatable :: plateau_sizes(:) !! flat-top lengths in samples
        integer, allocatable :: left_edges(:) !! one-based plateau starts
        integer, allocatable :: right_edges(:) !! one-based plateau ends
        integer :: status = peak_status_invalid_input !! zero on success
    end type peak_result

    public :: find_peaks, peak_prominences, peak_widths

contains

    pure logical function valid_sample_bound(values, n) result(valid)
        real(dp), intent(in), optional :: values(:) !! per-sample interval endpoint
        integer, intent(in) :: n !! signal length

        valid = .true.
        if (.not. present(values)) return
        valid = size(values) == n
        if (.not. valid) return
        valid = all(ieee_is_finite(values))
    end function valid_sample_bound


    pure subroutine peak_prominences(x, peaks, prominences, left_bases, right_bases, &
            status, wlen)
        real(dp), intent(in) :: x(:) !! finite, one-dimensional signal
        integer, intent(in) :: peaks(:) !! one-based indices at which to measure prominence
        real(dp), allocatable, intent(out) :: prominences(:) !! vertical contour differences
        integer, allocatable, intent(out) :: left_bases(:) !! nearest left minimum on ties
        integer, allocatable, intent(out) :: right_bases(:) !! nearest right minimum on ties
        integer, intent(out) :: status !! zero on success, one for invalid input
        integer, intent(in), optional :: wlen !! centered prominence window, at least 2

        integer :: i, j, k, left_limit, right_limit, radius, n, p
        real(dp) :: h, left_min, right_min

        n = size(x)
        allocate(prominences(size(peaks)), left_bases(size(peaks)), right_bases(size(peaks)))
        prominences = 0.0_dp
        left_bases = 0
        right_bases = 0
        status = peak_status_invalid_input
        if (any(.not. ieee_is_finite(x))) return
        if (any(peaks < 1) .or. any(peaks > n)) return
        radius = n
        if (present(wlen)) then
            if (wlen < 2) return
            radius = wlen / 2
        end if

        do k = 1, size(peaks)
            p = peaks(k)
            h = x(p)
            left_min = h
            right_min = h
            left_bases(k) = p
            right_bases(k) = p
            left_limit = max(1, p - radius)
            right_limit = min(n, p + radius)
            ! Equality at the peak height must not stop either search:
            ! another peak of the same height does not bound prominence.
            do i = p, left_limit, -1
                if (x(i) > h) exit
                if (x(i) < left_min) then
                    left_min = x(i)
                    left_bases(k) = i
                end if
            end do
            do j = p, right_limit
                if (x(j) > h) exit
                if (x(j) < right_min) then
                    right_min = x(j)
                    right_bases(k) = j
                end if
            end do
            prominences(k) = h - max(left_min, right_min)
        end do
        status = peak_status_success
    end subroutine peak_prominences

    pure subroutine peak_widths(x, peaks, widths, width_heights, left_ips, right_ips, &
            status, rel_height, wlen, prominences, left_bases, right_bases)
        real(dp), intent(in) :: x(:) !! finite signal
        integer, intent(in) :: peaks(:) !! one-based peak indices
        real(dp), allocatable, intent(out) :: widths(:) !! widths in sample spacings
        real(dp), allocatable, intent(out) :: width_heights(:) !! height at each width
        real(dp), allocatable, intent(out) :: left_ips(:) !! fractional one-based positions
        real(dp), allocatable, intent(out) :: right_ips(:) !! fractional one-based positions
        integer, intent(out) :: status !! zero on success
        real(dp), intent(in), optional :: rel_height !! fraction of prominence, default 0.5; >= 0
        integer, intent(in), optional :: wlen !! optional window for prominence calculation
        real(dp), intent(in), optional :: prominences(:) !! precomputed peak prominences
        integer, intent(in), optional :: left_bases(:) !! precomputed one-based left bases
        integer, intent(in), optional :: right_bases(:) !! precomputed one-based right bases

        integer :: i, k, np, p, ierr
        integer, allocatable :: lb(:), rb(:)
        real(dp) :: fraction, h, li, ri
        real(dp), allocatable :: prom(:)

        np = size(peaks)
        allocate(widths(np), width_heights(np), left_ips(np), right_ips(np))
        widths = 0.0_dp
        width_heights = 0.0_dp
        left_ips = 0.0_dp
        right_ips = 0.0_dp
        status = peak_status_invalid_input
        fraction = 0.5_dp
        if (present(rel_height)) fraction = rel_height
        if (.not. ieee_is_finite(fraction)) return
        if (fraction < 0.0_dp) return
        if (any(.not. ieee_is_finite(x))) return
        if (any(peaks < 1) .or. any(peaks > size(x))) return

        if (present(prominences) .and. present(left_bases) .and. present(right_bases)) then
            if (size(prominences) /= np .or. size(left_bases) /= np .or. &
                    size(right_bases) /= np) return
            if (any(.not. ieee_is_finite(prominences))) return
            if (any(prominences < 0.0_dp)) return
            if (any(left_bases > peaks) .or. any(right_bases < peaks)) return
            if (any(left_bases < 1) .or. any(right_bases > size(x))) return
            prom = prominences
            lb = left_bases
            rb = right_bases
        else
            call peak_prominences(x, peaks, prom, lb, rb, ierr, wlen)
            if (ierr /= peak_status_success) return
        end if

        do k = 1, np
            p = peaks(k)
            h = x(p) - prom(k) * fraction
            i = p
            do while (i > lb(k))
                if (x(i) <= h) exit
                i = i - 1
            end do
            li = real(i, dp)
            if (x(i) < h .and. i < p) then
                li = li + (h - x(i)) / (x(i + 1) - x(i))
            end if

            i = p
            do while (i < rb(k))
                if (x(i) <= h) exit
                i = i + 1
            end do
            ri = real(i, dp)
            if (x(i) < h .and. i > p) then
                ri = ri - (h - x(i)) / (x(i - 1) - x(i))
            end if
            width_heights(k) = h
            left_ips(k) = li
            right_ips(k) = ri
            widths(k) = ri - li
        end do
        status = peak_status_success
    end subroutine peak_widths

    pure function find_peaks(x, height_min, height_max, threshold_min, threshold_max, &
            distance, prominence_min, prominence_max, width_min, width_max, &
            wlen, rel_height, plateau_size_min, plateau_size_max, &
            height_min_samples, height_max_samples, threshold_min_samples, &
            threshold_max_samples, prominence_min_samples, prominence_max_samples, &
            width_min_samples, width_max_samples, plateau_size_min_samples, &
            plateau_size_max_samples, height_open, threshold_open, &
            prominence_open, width_open, plateau_size_open) result(res)
        real(dp), intent(in) :: x(:) !! finite one-dimensional data
        real(dp), intent(in), optional :: height_min !! minimum peak height (inclusive)
        real(dp), intent(in), optional :: height_max !! maximum peak height (inclusive)
        real(dp), intent(in), optional :: threshold_min !! minimum drop to each immediate neighbor
        real(dp), intent(in), optional :: threshold_max !! maximum drop to each immediate neighbor
        integer, intent(in), optional :: distance !! minimum spacing between centers, >= 1
        real(dp), intent(in), optional :: prominence_min !! minimum peak prominence
        real(dp), intent(in), optional :: prominence_max !! maximum peak prominence
        real(dp), intent(in), optional :: width_min !! minimum peak width in samples
        real(dp), intent(in), optional :: width_max !! maximum peak width in samples
        integer, intent(in), optional :: wlen !! local window for prominence or width, >= 2
        real(dp), intent(in), optional :: rel_height !! width as fraction of prominence, >= 0
        integer, intent(in), optional :: plateau_size_min !! minimum plateau length
        integer, intent(in), optional :: plateau_size_max !! maximum plateau length
        real(dp), intent(in), optional :: height_min_samples(:) !! per-sample inclusive lower height
        real(dp), intent(in), optional :: height_max_samples(:) !! per-sample inclusive upper height
        real(dp), intent(in), optional :: threshold_min_samples(:) !! per-sample lower neighbor drop
        real(dp), intent(in), optional :: threshold_max_samples(:) !! per-sample upper neighbor drop
        real(dp), intent(in), optional :: prominence_min_samples(:) !! per-sample lower prominence
        real(dp), intent(in), optional :: prominence_max_samples(:) !! per-sample upper prominence
        real(dp), intent(in), optional :: width_min_samples(:) !! per-sample lower width
        real(dp), intent(in), optional :: width_max_samples(:) !! per-sample upper width
        real(dp), intent(in), optional :: plateau_size_min_samples(:) !! per-sample lower plateau length
        real(dp), intent(in), optional :: plateau_size_max_samples(:) !! per-sample upper plateau length
        ! An open interval (None, None) requests measurements but does not
        ! filter peak positions. A true flag is equivalent to SciPy's open
        ! interval for that property; ordinary scalar/array filters still apply.
        logical, intent(in), optional :: height_open !! return peak heights without filtering
        logical, intent(in), optional :: threshold_open !! return neighbor drops without filtering
        logical, intent(in), optional :: prominence_open !! return prominence and bases
        logical, intent(in), optional :: width_open !! return widths and prominence details
        logical, intent(in), optional :: plateau_size_open !! return flat-top lengths and edges
        type(peak_result) :: res

        integer :: i, j, k, n, np, first, last, mid, chosen, ierr, spacing
        integer, allocatable :: pks(:), left(:), right(:), plevel(:), lb(:), rb(:)
        logical, allocatable :: active(:), handled(:)
        logical :: want_height, want_threshold, want_plateau, want_prom, want_width
        real(dp) :: v, best
        real(dp), allocatable :: heights(:), lt(:), rt(:), prom(:), wid(:), wh(:), li(:), ri(:)

        allocate(res%peaks(0), res%peak_heights(0), res%left_thresholds(0), &
            res%right_thresholds(0), res%prominences(0), res%left_bases(0), &
            res%right_bases(0), res%widths(0), res%width_heights(0), &
            res%left_ips(0), res%right_ips(0), res%plateau_sizes(0), &
            res%left_edges(0), res%right_edges(0))
        if (any(.not. ieee_is_finite(x))) return
        if (present(distance)) then
            if (distance < 1) return
        end if
        if (present(wlen)) then
            if (wlen < 2) return
        end if
        if (present(rel_height)) then
            if (.not. ieee_is_finite(rel_height)) return
            if (rel_height < 0.0_dp) return
        end if
        if (present(plateau_size_min)) then
            if (plateau_size_min < 0) return
        end if
        if (present(plateau_size_max)) then
            if (plateau_size_max < 0) return
        end if
        if (present(height_min)) then
            if (.not. ieee_is_finite(height_min)) return
        end if
        if (present(height_max)) then
            if (.not. ieee_is_finite(height_max)) return
        end if
        if (present(threshold_min)) then
            if (.not. ieee_is_finite(threshold_min)) return
        end if
        if (present(threshold_max)) then
            if (.not. ieee_is_finite(threshold_max)) return
        end if
        if (present(prominence_min)) then
            if (.not. ieee_is_finite(prominence_min)) return
        end if
        if (present(prominence_max)) then
            if (.not. ieee_is_finite(prominence_max)) return
        end if
        if (present(width_min)) then
            if (.not. ieee_is_finite(width_min)) return
        end if
        if (present(width_max)) then
            if (.not. ieee_is_finite(width_max)) return
        end if

        n = size(x)
        ! Dynamic bounds are indexed at the *center* of each peak, just as
        ! scipy.signal.find_peaks indexes position-dependent intervals.
        if (.not. valid_sample_bound(height_min_samples, n)) return
        if (.not. valid_sample_bound(height_max_samples, n)) return
        if (.not. valid_sample_bound(threshold_min_samples, n)) return
        if (.not. valid_sample_bound(threshold_max_samples, n)) return
        if (.not. valid_sample_bound(prominence_min_samples, n)) return
        if (.not. valid_sample_bound(prominence_max_samples, n)) return
        if (.not. valid_sample_bound(width_min_samples, n)) return
        if (.not. valid_sample_bound(width_max_samples, n)) return
        if (.not. valid_sample_bound(plateau_size_min_samples, n)) return
        if (.not. valid_sample_bound(plateau_size_max_samples, n)) return
        if (present(plateau_size_min_samples)) then
            if (any(plateau_size_min_samples < 0.0_dp)) return
        end if
        if (present(plateau_size_max_samples)) then
            if (any(plateau_size_max_samples < 0.0_dp)) return
        end if
        allocate(pks(n), left(n), right(n), plevel(n), active(n), handled(n))
        np = 0
        i = 2
        do while (i < n)
            if (x(i - 1) < x(i)) then
                first = i
                j = i + 1
                do while (j <= n)
                    if (x(j) > x(i) .or. x(j) < x(i)) exit
                    j = j + 1
                end do
                if (j <= n) then
                    if (x(j) < x(i)) then
                        last = j - 1
                        mid = (first + last) / 2
                        np = np + 1
                        pks(np) = mid
                        left(np) = first
                        right(np) = last
                        plevel(np) = last - first + 1
                    end if
                end if
                i = j
            else
                i = i + 1
            end if
        end do
        active = .false.
        active(1:np) = .true.
        want_plateau = present(plateau_size_min) .or. present(plateau_size_max) .or. &
            present(plateau_size_min_samples) .or. present(plateau_size_max_samples)
        want_height = present(height_min) .or. present(height_max) .or. &
            present(height_min_samples) .or. present(height_max_samples)
        want_threshold = present(threshold_min) .or. present(threshold_max) .or. &
            present(threshold_min_samples) .or. present(threshold_max_samples)
        want_prom = present(prominence_min) .or. present(prominence_max) .or. &
            present(width_min) .or. present(width_max) .or. &
            present(prominence_min_samples) .or. present(prominence_max_samples) .or. &
            present(width_min_samples) .or. present(width_max_samples)
        want_width = present(width_min) .or. present(width_max) .or. &
            present(width_min_samples) .or. present(width_max_samples)
        if (present(height_open)) then
            if (height_open) want_height = .true.
        end if
        if (present(threshold_open)) then
            if (threshold_open) want_threshold = .true.
        end if
        if (present(plateau_size_open)) then
            if (plateau_size_open) want_plateau = .true.
        end if
        if (present(prominence_open)) then
            if (prominence_open) want_prom = .true.
        end if
        if (present(width_open)) then
            if (width_open) then
                want_width = .true.
                want_prom = .true.
            end if
        end if

        do k = 1, np
            if (present(plateau_size_min)) then
                if (plevel(k) < plateau_size_min) active(k) = .false.
            end if
            if (present(plateau_size_max)) then
                if (plevel(k) > plateau_size_max) active(k) = .false.
            end if
            if (present(plateau_size_min_samples)) then
                if (real(plevel(k), dp) < plateau_size_min_samples(pks(k))) active(k) = .false.
            end if
            if (present(plateau_size_max_samples)) then
                if (real(plevel(k), dp) > plateau_size_max_samples(pks(k))) active(k) = .false.
            end if
            if (present(height_min)) then
                if (x(pks(k)) < height_min) active(k) = .false.
            end if
            if (present(height_max)) then
                if (x(pks(k)) > height_max) active(k) = .false.
            end if
            if (present(height_min_samples)) then
                if (x(pks(k)) < height_min_samples(pks(k))) active(k) = .false.
            end if
            if (present(height_max_samples)) then
                if (x(pks(k)) > height_max_samples(pks(k))) active(k) = .false.
            end if
            if (present(threshold_min)) then
                if (min(x(pks(k)) - x(pks(k) - 1), &
                        x(pks(k)) - x(pks(k) + 1)) < threshold_min) active(k) = .false.
            end if
            if (present(threshold_max)) then
                if (max(x(pks(k)) - x(pks(k) - 1), &
                        x(pks(k)) - x(pks(k) + 1)) > threshold_max) active(k) = .false.
            end if
            if (present(threshold_min_samples)) then
                if (min(x(pks(k)) - x(pks(k) - 1), &
                        x(pks(k)) - x(pks(k) + 1)) < threshold_min_samples(pks(k))) active(k) = .false.
            end if
            if (present(threshold_max_samples)) then
                if (max(x(pks(k)) - x(pks(k) - 1), &
                        x(pks(k)) - x(pks(k) + 1)) > threshold_max_samples(pks(k))) active(k) = .false.
            end if
        end do

        ! Greedily retain taller peaks before shorter neighbors. SciPy also
        ! prioritizes peak amplitude, not prominence, for distance exclusion.
        spacing = 1
        if (present(distance)) spacing = distance
        if (spacing > 1) then
            handled = .false.
            do i = 1, np
                chosen = 0
                best = -huge(0.0_dp)
                do j = 1, np
                    if (.not. active(j) .or. handled(j)) cycle
                    v = x(pks(j))
                    if (v >= best) then
                        best = v
                        chosen = j
                    end if
                end do
                if (chosen == 0) exit
                handled(chosen) = .true.
                do j = 1, np
                    if (j == chosen .or. .not. active(j)) cycle
                    if (abs(pks(j) - pks(chosen)) < spacing) active(j) = .false.
                end do
            end do
        end if

        ! Evaluate geometric properties only on peaks that survive cheap
        ! filters and distance. Keep the original ascending index order.
        pks = pack(pks(1:np), active(1:np))
        left = pack(left(1:np), active(1:np))
        right = pack(right(1:np), active(1:np))
        plevel = pack(plevel(1:np), active(1:np))
        np = size(pks)
        deallocate(active)
        allocate(active(np), heights(np), lt(np), rt(np))
        active = .true.
        do k = 1, np
            heights(k) = x(pks(k))
            lt(k) = heights(k) - x(pks(k) - 1)
            rt(k) = heights(k) - x(pks(k) + 1)
        end do
        if (want_prom) then
            call peak_prominences(x, pks(1:np), prom, lb, rb, ierr, wlen)
            if (ierr /= peak_status_success) return
            do k = 1, np
                if (present(prominence_min)) then
                    if (prom(k) < prominence_min) active(k) = .false.
                end if
                if (present(prominence_max)) then
                    if (prom(k) > prominence_max) active(k) = .false.
                end if
                if (present(prominence_min_samples)) then
                    if (prom(k) < prominence_min_samples(pks(k))) active(k) = .false.
                end if
                if (present(prominence_max_samples)) then
                    if (prom(k) > prominence_max_samples(pks(k))) active(k) = .false.
                end if
            end do
            if (want_width) then
                call peak_widths(x, pks(1:np), wid, wh, li, ri, ierr, rel_height, &
                    wlen, prom, lb, rb)
                if (ierr /= peak_status_success) return
                do k = 1, np
                    if (present(width_min)) then
                        if (wid(k) < width_min) active(k) = .false.
                    end if
                    if (present(width_max)) then
                        if (wid(k) > width_max) active(k) = .false.
                    end if
                    if (present(width_min_samples)) then
                        if (wid(k) < width_min_samples(pks(k))) active(k) = .false.
                    end if
                    if (present(width_max_samples)) then
                        if (wid(k) > width_max_samples(pks(k))) active(k) = .false.
                    end if
                end do
            end if
        end if
        res%peaks = pack(pks(1:np), active)
        if (want_plateau) then
            res%plateau_sizes = pack(plevel(1:np), active)
            res%left_edges = pack(left(1:np), active)
            res%right_edges = pack(right(1:np), active)
        end if
        if (want_height) res%peak_heights = pack(heights, active)
        if (want_threshold) then
            res%left_thresholds = pack(lt, active)
            res%right_thresholds = pack(rt, active)
        end if
        if (want_prom) then
            res%prominences = pack(prom, active)
            res%left_bases = pack(lb, active)
            res%right_bases = pack(rb, active)
        end if
        if (want_width) then
            res%widths = pack(wid, active)
            res%width_heights = pack(wh, active)
            res%left_ips = pack(li, active)
            res%right_ips = pack(ri, active)
        end if
        res%status = peak_status_success
    end function find_peaks

end module scifort_signal_peaks
