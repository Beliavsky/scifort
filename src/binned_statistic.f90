! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_binned_statistic
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    integer, parameter, public :: binned_statistic_status_success = 0
    integer, parameter, public :: binned_statistic_status_invalid = 1

    abstract interface
        pure function binned_statistic_callback(values) result(value)
            import dp
            real(dp), intent(in) :: values(:) !! values assigned to one bin
            real(dp) :: value
        end function binned_statistic_callback
    end interface

    type, public :: binned_statistic_result
        real(dp), allocatable :: statistic(:) !! statistic for each interior bin
        real(dp), allocatable :: bin_edges(:) !! bin edges, including both endpoints
        integer, allocatable :: binnumber(:) !! SciPy-style bin number for each observation
        integer :: status = binned_statistic_status_success !! zero on success
    end type binned_statistic_result

    type, public :: binned_statistic_2d_result
        real(dp), allocatable :: statistic(:, :) !! statistic indexed by x-bin then y-bin
        real(dp), allocatable :: x_edge(:) !! x-dimension bin edges
        real(dp), allocatable :: y_edge(:) !! y-dimension bin edges
        integer, allocatable :: binnumber(:) !! SciPy C-order flattened bin numbers
        integer, allocatable :: binnumber_expanded(:, :) !! per-dimension bin numbers, shape (2,n)
        integer :: status = binned_statistic_status_success !! zero on success
    end type binned_statistic_2d_result

    type, public :: binned_statistic_dd_result
        real(dp), allocatable :: statistic(:) !! C-order flattened interior-bin statistic
        integer, allocatable :: statistic_shape(:) !! number of interior bins per dimension
        real(dp), allocatable :: bin_edges(:, :) !! padded edge columns, one column per dimension
        integer, allocatable :: edge_counts(:) !! valid edge count in each edge column
        integer, allocatable :: binnumber(:) !! SciPy C-order flattened bin numbers
        integer, allocatable :: binnumber_expanded(:, :) !! per-dimension bin numbers, shape (d,n)
        integer :: status = binned_statistic_status_success !! zero on success
    end type binned_statistic_dd_result

    public :: binned_statistic
    public :: binned_statistic_2d
    public :: binned_statistic_dd
    public :: binned_statistic_callback

contains

    pure function binned_statistic(x, values, statistic, bins, range, bin_edges, stat_func) result(res)
        real(dp), intent(in) :: x(:) !! coordinates to bin
        real(dp), intent(in), optional :: values(:) !! values summarized within bins
        character(len=*), intent(in), optional :: statistic !! mean/std/median/count/sum/min/max
        integer, intent(in), optional :: bins !! equal-width bin count, default 10
        real(dp), intent(in), optional :: range(2) !! lower and upper limits for equal-width bins
        real(dp), intent(in), optional :: bin_edges(:) !! explicit strictly increasing bin edges
        procedure(binned_statistic_callback), optional :: stat_func !! custom statistic callback
        type(binned_statistic_result) :: res

        integer :: i
        integer :: nbins
        integer :: status
        integer, allocatable :: core_bin(:)
        real(dp), allocatable :: edges(:)

        call prepare_edges_1d(x, bins, range, bin_edges, edges, status)
        if (status /= binned_statistic_status_success) then
            call invalidate_1d(res)
            return
        end if

        nbins = size(edges) - 1
        allocate(res%bin_edges(size(edges)), res%binnumber(size(x)), core_bin(size(x)))
        res%bin_edges = edges
        core_bin = 0
        do i = 1, size(x)
            res%binnumber(i) = locate_bin(x(i), edges)
            if (res%binnumber(i) >= 1 .and. res%binnumber(i) <= nbins) then
                core_bin(i) = res%binnumber(i)
            end if
        end do

        call calculate_statistic(core_bin, nbins, values, statistic, stat_func, &
            res%statistic, status)
        res%status = status
        if (status /= binned_statistic_status_success) then
            if (allocated(res%statistic)) deallocate(res%statistic)
            allocate(res%statistic(0))
        end if
    end function binned_statistic

    pure function binned_statistic_2d(x, y, values, statistic, bins, range, x_edges, y_edges, &
        stat_func) result(res)
        real(dp), intent(in) :: x(:) !! first coordinates to bin
        real(dp), intent(in) :: y(:) !! second coordinates to bin
        real(dp), intent(in), optional :: values(:) !! values summarized within bins
        character(len=*), intent(in), optional :: statistic !! mean/std/median/count/sum/min/max
        integer, intent(in), optional :: bins(:) !! one shared or two per-axis equal-width counts
        real(dp), intent(in), optional :: range(2, 2) !! lower/upper limits by dimension
        real(dp), intent(in), optional :: x_edges(:) !! explicit x-axis bin edges
        real(dp), intent(in), optional :: y_edges(:) !! explicit y-axis bin edges
        procedure(binned_statistic_callback), optional :: stat_func !! custom statistic callback
        type(binned_statistic_2d_result) :: res

        integer :: bx
        integer :: by
        integer :: i
        integer :: nx
        integer :: ny
        integer :: status
        integer :: bx_count
        integer :: by_count
        integer, allocatable :: core_bin(:)
        integer, allocatable :: bins_value(:)
        real(dp), allocatable :: flat_stat(:)
        real(dp), allocatable :: xe(:)
        real(dp), allocatable :: ye(:)

        if (size(x) /= size(y)) then
            call invalidate_2d(res)
            return
        end if

        allocate(bins_value(2))
        bins_value = 10
        if (present(bins)) then
            if (size(bins) == 1) then
                bins_value = bins(1)
            else if (size(bins) == 2) then
                bins_value = bins
            else
                call invalidate_2d(res)
                return
            end if
        end if

        if (present(range)) then
            call prepare_edges_1d(x, bins_value(1), range(:, 1), x_edges, xe, status)
        else
            call prepare_edges_1d(x, bins_value(1), explicit_edges=x_edges, edges=xe, status=status)
        end if
        if (status /= binned_statistic_status_success) then
            call invalidate_2d(res)
            return
        end if
        if (present(range)) then
            call prepare_edges_1d(y, bins_value(2), range(:, 2), y_edges, ye, status)
        else
            call prepare_edges_1d(y, bins_value(2), explicit_edges=y_edges, edges=ye, status=status)
        end if
        if (status /= binned_statistic_status_success) then
            call invalidate_2d(res)
            return
        end if

        nx = size(xe) - 1
        ny = size(ye) - 1
        bx_count = nx + 2
        by_count = ny + 2
        allocate(res%x_edge(size(xe)), res%y_edge(size(ye)))
        allocate(res%binnumber(size(x)), res%binnumber_expanded(2, size(x)))
        allocate(core_bin(size(x)))
        res%x_edge = xe
        res%y_edge = ye
        core_bin = 0

        do i = 1, size(x)
            bx = locate_bin(x(i), xe)
            by = locate_bin(y(i), ye)
            res%binnumber_expanded(1, i) = bx
            res%binnumber_expanded(2, i) = by
            res%binnumber(i) = bx * by_count + by
            if (bx >= 1 .and. bx <= nx .and. by >= 1 .and. by <= ny) then
                core_bin(i) = (bx - 1) * ny + by
            end if
        end do

        call calculate_statistic(core_bin, nx * ny, values, statistic, stat_func, flat_stat, status)
        res%status = status
        if (status /= binned_statistic_status_success) then
            allocate(res%statistic(0, 0))
            return
        end if

        allocate(res%statistic(nx, ny))
        do bx = 1, nx
            do by = 1, ny
                res%statistic(bx, by) = flat_stat((bx - 1) * ny + by)
            end do
        end do
    end function binned_statistic_2d

    pure function binned_statistic_dd(sample, values, statistic, bins, range, bin_edges, &
        edge_counts, stat_func) result(res)
        real(dp), intent(in) :: sample(:, :) !! observations by dimension, shape (n,d)
        real(dp), intent(in), optional :: values(:) !! values summarized within bins
        character(len=*), intent(in), optional :: statistic !! mean/std/median/count/sum/min/max
        integer, intent(in), optional :: bins(:) !! one shared or d per-axis equal-width counts
        real(dp), intent(in), optional :: range(:, :) !! shape (2,d), lower/upper limits
        real(dp), intent(in), optional :: bin_edges(:, :) !! padded explicit edge columns
        integer, intent(in), optional :: edge_counts(:) !! valid edge count per explicit column
        procedure(binned_statistic_callback), optional :: stat_func !! custom statistic callback
        type(binned_statistic_dd_result) :: res

        integer :: d
        integer :: dim
        integer :: i
        integer :: idx0
        integer :: max_edges
        integer :: n
        integer :: ncore
        integer :: status
        integer, allocatable :: bins_value(:)
        integer, allocatable :: core_bin(:)
        integer, allocatable :: counts(:)
        integer, allocatable :: expanded(:, :)
        real(dp), allocatable :: edges(:, :)
        real(dp), allocatable :: one_edges(:)

        n = size(sample, 1)
        d = size(sample, 2)
        if (d < 1) then
            call invalidate_dd(res)
            return
        end if
        if (present(range)) then
            if (size(range, 1) /= 2 .or. size(range, 2) /= d) then
                call invalidate_dd(res)
                return
            end if
        end if

        allocate(bins_value(d))
        bins_value = 10
        if (present(bins)) then
            if (size(bins) == 1) then
                bins_value = bins(1)
            else if (size(bins) == d) then
                bins_value = bins
            else
                call invalidate_dd(res)
                return
            end if
        end if
        if (any(bins_value < 1)) then
            call invalidate_dd(res)
            return
        end if

        allocate(counts(d))
        if (present(bin_edges)) then
            if (size(bin_edges, 2) /= d) then
                call invalidate_dd(res)
                return
            end if
            if (present(edge_counts)) then
                if (size(edge_counts) /= d) then
                    call invalidate_dd(res)
                    return
                end if
                counts = edge_counts
            else
                counts = size(bin_edges, 1)
            end if
            if (any(counts < 2) .or. any(counts > size(bin_edges, 1))) then
                call invalidate_dd(res)
                return
            end if
            max_edges = maxval(counts)
        else
            counts = bins_value + 1
            max_edges = maxval(counts)
        end if

        allocate(edges(max_edges, d))
        edges = quiet_nan(0.0_dp)
        do dim = 1, d
            if (present(bin_edges)) then
                call prepare_edges_1d(sample(:, dim), explicit_edges=bin_edges(1:counts(dim), dim), &
                    edges=one_edges, status=status)
            else if (present(range)) then
                call prepare_edges_1d(sample(:, dim), bins_value(dim), range(:, dim), &
                    edges=one_edges, status=status)
            else
                call prepare_edges_1d(sample(:, dim), bins_value(dim), edges=one_edges, status=status)
            end if
            if (status /= binned_statistic_status_success) then
                call invalidate_dd(res)
                return
            end if
            counts(dim) = size(one_edges)
            edges(1:counts(dim), dim) = one_edges
        end do

        allocate(res%statistic_shape(d), res%edge_counts(d))
        res%statistic_shape = counts - 1
        res%edge_counts = counts
        ncore = 1
        do dim = 1, d
            if (res%statistic_shape(dim) > huge(ncore) / ncore) then
                call invalidate_dd(res)
                return
            end if
            ncore = ncore * res%statistic_shape(dim)
        end do

        allocate(res%bin_edges(max_edges, d))
        allocate(res%binnumber(n), res%binnumber_expanded(d, n))
        allocate(expanded(d, n), core_bin(n))
        res%bin_edges = edges
        core_bin = 0

        do i = 1, n
            do dim = 1, d
                expanded(dim, i) = locate_bin(sample(i, dim), edges(1:counts(dim), dim))
            end do
            res%binnumber_expanded(:, i) = expanded(:, i)
            idx0 = expanded(1, i)
            do dim = 2, d
                idx0 = idx0 * (res%statistic_shape(dim) + 2) + expanded(dim, i)
            end do
            res%binnumber(i) = idx0

            if (all(expanded(:, i) >= 1 .and. expanded(:, i) <= res%statistic_shape)) then
                idx0 = expanded(1, i) - 1
                do dim = 2, d
                    idx0 = idx0 * res%statistic_shape(dim) + expanded(dim, i) - 1
                end do
                core_bin(i) = idx0 + 1
            end if
        end do

        call calculate_statistic(core_bin, ncore, values, statistic, stat_func, &
            res%statistic, status)
        res%status = status
        if (status /= binned_statistic_status_success) then
            if (allocated(res%statistic)) deallocate(res%statistic)
            allocate(res%statistic(0))
        end if
    end function binned_statistic_dd

    pure subroutine prepare_edges_1d(x, bins, range, explicit_edges, edges, status)
        real(dp), intent(in) :: x(:) !! coordinates used to infer automatic limits
        integer, intent(in), optional :: bins !! equal-width bin count
        real(dp), intent(in), optional :: range(2) !! requested lower and upper limits
        real(dp), intent(in), optional :: explicit_edges(:) !! requested bin edges
        real(dp), allocatable, intent(out) :: edges(:) !! validated/generated bin edges
        integer, intent(out) :: status !! zero on success

        integer :: i
        integer :: nbins
        real(dp) :: lower
        real(dp) :: upper

        status = binned_statistic_status_success
        if (present(explicit_edges)) then
            if (size(explicit_edges) < 2) then
                status = binned_statistic_status_invalid
                allocate(edges(0))
                return
            end if
            if (any(.not. ieee_is_finite(explicit_edges))) then
                status = binned_statistic_status_invalid
                allocate(edges(0))
                return
            end if
            do i = 2, size(explicit_edges)
                if (explicit_edges(i) <= explicit_edges(i - 1)) then
                    status = binned_statistic_status_invalid
                    allocate(edges(0))
                    return
                end if
            end do
            allocate(edges(size(explicit_edges)))
            edges = explicit_edges
            return
        end if

        nbins = 10
        if (present(bins)) nbins = bins
        if (nbins < 1) then
            status = binned_statistic_status_invalid
            allocate(edges(0))
            return
        end if

        if (present(range)) then
            lower = range(1)
            upper = range(2)
            if (.not. ieee_is_finite(lower) .or. .not. ieee_is_finite(upper) .or. upper < lower) then
                status = binned_statistic_status_invalid
                allocate(edges(0))
                return
            end if
        else
            if (size(x) == 0 .or. any(.not. ieee_is_finite(x))) then
                status = binned_statistic_status_invalid
                allocate(edges(0))
                return
            end if
            lower = minval(x)
            upper = maxval(x)
        end if

        if (lower == upper) then
            lower = lower - 0.5_dp
            upper = upper + 0.5_dp
        end if

        allocate(edges(nbins + 1))
        if (nbins == 1) then
            edges = [lower, upper]
        else
            do i = 1, nbins + 1
                edges(i) = lower + (upper - lower) * real(i - 1, dp) / real(nbins, dp)
            end do
            edges(nbins + 1) = upper
        end if
    end subroutine prepare_edges_1d

    pure integer function locate_bin(value, edges) result(bin)
        real(dp), intent(in) :: value !! coordinate to locate
        real(dp), intent(in) :: edges(:) !! strictly increasing bin edges

        integer :: i
        integer :: nbins
        real(dp) :: scale
        real(dp) :: tolerance

        nbins = size(edges) - 1
        if (ieee_is_nan(value)) then
            bin = nbins + 1
            return
        end if
        if (value < edges(1)) then
            bin = 0
            return
        end if
        if (value > edges(nbins + 1)) then
            scale = max(1.0_dp, abs(edges(nbins + 1)))
            tolerance = 32.0_dp * epsilon(1.0_dp) * scale
            if (value - edges(nbins + 1) <= tolerance) then
                bin = nbins
            else
                bin = nbins + 1
            end if
            return
        end if
        if (value == edges(nbins + 1)) then
            bin = nbins
            return
        end if

        do i = 2, nbins + 1
            if (value < edges(i)) then
                bin = i - 1
                return
            end if
        end do
        bin = nbins
    end function locate_bin

    pure subroutine calculate_statistic(bin_ids, n_bins, values, statistic, stat_func, output, status)
        integer, intent(in) :: bin_ids(:) !! interior flat bin id, or zero for outliers
        integer, intent(in) :: n_bins !! total number of interior bins
        real(dp), intent(in), optional :: values(:) !! observation values
        character(len=*), intent(in), optional :: statistic !! built-in statistic name
        procedure(binned_statistic_callback), optional :: stat_func !! custom statistic callback
        real(dp), allocatable, intent(out) :: output(:) !! statistic by flat interior bin
        integer, intent(out) :: status !! zero on success

        character(len=16) :: name
        integer :: b
        integer :: i
        integer :: k
        integer, allocatable :: counts(:)
        logical, allocatable :: has_nan(:)
        real(dp), allocatable :: scratch(:)
        real(dp), allocatable :: sums(:)

        status = binned_statistic_status_success
        allocate(output(n_bins))
        if (present(stat_func)) then
            if (.not. present(values)) then
                status = binned_statistic_status_invalid
                return
            end if
            if (size(values) /= size(bin_ids)) then
                status = binned_statistic_status_invalid
                return
            end if
            do b = 1, n_bins
                k = count(bin_ids == b)
                allocate(scratch(k))
                if (k > 0) scratch = pack(values, bin_ids == b)
                output(b) = stat_func(scratch)
                deallocate(scratch)
            end do
            return
        end if

        name = 'mean'
        if (present(statistic)) name = lowercase(trim(adjustl(statistic)))
        if (name /= 'count') then
            if (.not. present(values)) then
                status = binned_statistic_status_invalid
                return
            end if
            if (size(values) /= size(bin_ids)) then
                status = binned_statistic_status_invalid
                return
            end if
        end if

        select case (trim(name))
        case ('count')
            output = 0.0_dp
            do i = 1, size(bin_ids)
                if (bin_ids(i) > 0) output(bin_ids(i)) = output(bin_ids(i)) + 1.0_dp
            end do

        case ('sum')
            output = 0.0_dp
            do i = 1, size(bin_ids)
                if (bin_ids(i) > 0) output(bin_ids(i)) = output(bin_ids(i)) + values(i)
            end do

        case ('mean')
            allocate(counts(n_bins), sums(n_bins))
            counts = 0
            sums = 0.0_dp
            do i = 1, size(bin_ids)
                if (bin_ids(i) > 0) then
                    counts(bin_ids(i)) = counts(bin_ids(i)) + 1
                    sums(bin_ids(i)) = sums(bin_ids(i)) + values(i)
                end if
            end do
            output = quiet_nan(0.0_dp)
            do b = 1, n_bins
                if (counts(b) > 0) output(b) = sums(b) / real(counts(b), dp)
            end do

        case ('std')
            allocate(counts(n_bins), sums(n_bins))
            counts = 0
            sums = 0.0_dp
            do i = 1, size(bin_ids)
                if (bin_ids(i) > 0) then
                    counts(bin_ids(i)) = counts(bin_ids(i)) + 1
                    sums(bin_ids(i)) = sums(bin_ids(i)) + values(i)
                end if
            end do
            output = 0.0_dp
            do i = 1, size(bin_ids)
                if (bin_ids(i) > 0) then
                    b = bin_ids(i)
                    output(b) = output(b) + (values(i) - sums(b) / real(counts(b), dp)) ** 2
                end if
            end do
            do b = 1, n_bins
                if (counts(b) > 0) then
                    output(b) = sqrt(output(b) / real(counts(b), dp))
                else
                    output(b) = quiet_nan(0.0_dp)
                end if
            end do

        case ('median')
            output = quiet_nan(0.0_dp)
            do b = 1, n_bins
                k = count(bin_ids == b)
                if (k > 0) then
                    allocate(scratch(k))
                    scratch = pack(values, bin_ids == b)
                    if (any(ieee_is_nan(scratch))) then
                        output(b) = quiet_nan(0.0_dp)
                    else
                        call sort_real(scratch)
                        if (mod(k, 2) == 1) then
                            output(b) = scratch((k + 1) / 2)
                        else
                            output(b) = 0.5_dp * (scratch(k / 2) + scratch(k / 2 + 1))
                        end if
                    end if
                    deallocate(scratch)
                end if
            end do

        case ('min')
            output = quiet_nan(0.0_dp)
            do b = 1, n_bins
                k = count(bin_ids == b .and. .not. ieee_is_nan(values))
                if (k > 0) then
                    allocate(scratch(k))
                    scratch = pack(values, bin_ids == b .and. .not. ieee_is_nan(values))
                    output(b) = minval(scratch)
                    deallocate(scratch)
                end if
            end do

        case ('max')
            allocate(has_nan(n_bins))
            has_nan = .false.
            output = quiet_nan(0.0_dp)
            do b = 1, n_bins
                k = count(bin_ids == b)
                if (k > 0) then
                    has_nan(b) = any(ieee_is_nan(values) .and. bin_ids == b)
                    if (.not. has_nan(b)) then
                        allocate(scratch(k))
                        scratch = pack(values, bin_ids == b)
                        output(b) = maxval(scratch)
                        deallocate(scratch)
                    end if
                end if
            end do

        case default
            status = binned_statistic_status_invalid
        end select
    end subroutine calculate_statistic

    pure subroutine sort_real(values)
        real(dp), intent(inout) :: values(:) !! values sorted ascending in place

        integer :: i
        integer :: j
        real(dp) :: key

        do i = 2, size(values)
            key = values(i)
            j = i - 1
            do while (j >= 1)
                if (values(j) <= key) exit
                values(j + 1) = values(j)
                j = j - 1
            end do
            values(j + 1) = key
        end do
    end subroutine sort_real

    pure function lowercase(text) result(lower)
        character(len=*), intent(in) :: text !! text converted to lowercase
        character(len=len(text)) :: lower

        integer :: code
        integer :: i

        lower = text
        do i = 1, len(text)
            code = iachar(text(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                lower(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
        end do
    end function lowercase

    pure subroutine invalidate_1d(res)
        type(binned_statistic_result), intent(out) :: res !! result marked invalid

        allocate(res%statistic(0), res%bin_edges(0), res%binnumber(0))
        res%status = binned_statistic_status_invalid
    end subroutine invalidate_1d

    pure subroutine invalidate_2d(res)
        type(binned_statistic_2d_result), intent(out) :: res !! result marked invalid

        allocate(res%statistic(0, 0), res%x_edge(0), res%y_edge(0))
        allocate(res%binnumber(0), res%binnumber_expanded(0, 0))
        res%status = binned_statistic_status_invalid
    end subroutine invalidate_2d

    pure subroutine invalidate_dd(res)
        type(binned_statistic_dd_result), intent(out) :: res !! result marked invalid

        allocate(res%statistic(0), res%statistic_shape(0), res%bin_edges(0, 0))
        allocate(res%edge_counts(0), res%binnumber(0), res%binnumber_expanded(0, 0))
        res%status = binned_statistic_status_invalid
    end subroutine invalidate_dd

end module scifort_binned_statistic
