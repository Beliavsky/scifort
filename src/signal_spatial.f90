! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! SciPy-style zero-padded order-statistic and separable two-dimensional FIR filters.
! Independent Fortran implementation; SciPy is used only for verification.
module scifort_signal_spatial
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: spatial_status_success = 0
    integer, parameter, public :: spatial_status_invalid_input = 1
    public :: order_filter, sepfir2d

    interface order_filter
        module procedure order_real_int_1d
        module procedure order_real_int_2d
        module procedure order_real_int_3d
        module procedure order_real_real_1d
        module procedure order_real_real_2d
        module procedure order_real_real_3d
        module procedure order_real_logical_1d
        module procedure order_real_logical_2d
        module procedure order_real_logical_3d
        module procedure order_int_int_1d
        module procedure order_int_int_2d
        module procedure order_int_int_3d
        module procedure order_int_real_1d
        module procedure order_int_real_2d
        module procedure order_int_real_3d
        module procedure order_int_logical_1d
        module procedure order_int_logical_2d
        module procedure order_int_logical_3d
    end interface

    interface sepfir2d
        module procedure sep_real_real_real
        module procedure sep_real_real_complex
        module procedure sep_real_complex_real
        module procedure sep_real_complex_complex
        module procedure sep_complex_real_real
        module procedure sep_complex_real_complex
        module procedure sep_complex_complex_real
        module procedure sep_complex_complex_complex
        module procedure sep_int_real_real
        module procedure sep_int_real_complex
        module procedure sep_int_complex_real
        module procedure sep_int_complex_complex
    end interface

contains

    subroutine shape_strides(dims, stride)
        integer, intent(in) :: dims(:) !! original extents
        integer, intent(out) :: stride(:) !! Fortran column-major strides
        integer :: d

        if (size(dims) == 0) return
        stride(1) = 1
        do d = 2, size(dims)
            stride(d) = stride(d-1) * dims(d-1)
        end do
    end subroutine shape_strides

    subroutine centered_mask_offsets(dims, mask, offsets)
        integer, intent(in) :: dims(:) !! odd footprint dimensions
        logical, intent(in) :: mask(:) !! flat footprint mask
        integer, intent(out) :: offsets(:,:) !! active centered offsets
        integer :: i, q, d, k

        k = 0
        do i = 1, size(mask)
            if (.not. mask(i)) cycle
            k = k + 1
            q = i - 1
            do d = 1, size(dims)
                offsets(d,k) = mod(q, dims(d)) - dims(d)/2
                q = q / dims(d)
            end do
        end do
    end subroutine centered_mask_offsets

    subroutine order_index(i, dims, stride, offsets, idx)
        integer, intent(in) :: i !! one-based flat output element
        integer, intent(in) :: dims(:) !! input dimensions
        integer, intent(in) :: stride(:) !! flat strides
        integer, intent(in) :: offsets(:) !! neighbor displacements
        integer, intent(out) :: idx !! one-based neighbor flat index, zero for padding
        integer :: q, d, c

        idx = 1
        q = i - 1
        do d = 1, size(dims)
            c = mod(q, dims(d)) + 1 + offsets(d)
            q = q / dims(d)
            if (c < 1 .or. c > dims(d)) then
                idx = 0
                return
            end if
            idx = idx + (c - 1) * stride(d)
        end do
    end subroutine order_index

    real(dp) function kth_real(work, rank) result(value)
        real(dp), intent(inout) :: work(:) !! scratch array modified in place
        integer, intent(in) :: rank !! zero-based rank
        real(dp) :: pivot, temp
        integer :: left, right, low, high, i, target

        target = rank + 1
        left = 1
        right = size(work)
        do while (left < right)
            pivot = work(left + (right - left)/2)
            low = left
            high = right
            i = left
            do while (i <= high)
                if (work(i) < pivot) then
                    temp = work(low)
                    work(low) = work(i)
                    work(i) = temp
                    low = low + 1
                    i = i + 1
                else if (work(i) > pivot) then
                    temp = work(high)
                    work(high) = work(i)
                    work(i) = temp
                    high = high - 1
                else
                    i = i + 1
                end if
            end do
            if (target < low) then
                right = low - 1
            else if (target > high) then
                left = high + 1
            else
                exit
            end if
        end do
        value = work(target)
    end function kth_real

    integer function kth_int(work, rank) result(value)
        integer, intent(inout) :: work(:) !! scratch array modified in place
        integer, intent(in) :: rank !! zero-based rank
        integer :: pivot, temp
        integer :: left, right, low, high, i, target

        target = rank + 1
        left = 1
        right = size(work)
        do while (left < right)
            pivot = work(left + (right - left)/2)
            low = left
            high = right
            i = left
            do while (i <= high)
                if (work(i) < pivot) then
                    temp = work(low)
                    work(low) = work(i)
                    work(i) = temp
                    low = low + 1
                    i = i + 1
                else if (work(i) > pivot) then
                    temp = work(high)
                    work(high) = work(i)
                    work(i) = temp
                    high = high - 1
                else
                    i = i + 1
                end if
            end do
            if (target < low) then
                right = low - 1
            else if (target > high) then
                left = high + 1
            else
                exit
            end if
        end do
        value = work(target)
    end function kth_int

    subroutine prepare_mask(dims, domain_dims, mask, rank, offsets, stride, status)
        integer, intent(in) :: dims(:) !! signal extents
        integer, intent(in) :: domain_dims(:) !! mask dimensions
        logical, intent(in) :: mask(:) !! flattened mask
        integer, intent(in) :: rank !! zero-based order rank
        integer, allocatable, intent(out) :: offsets(:,:) !! active mask offsets
        integer, intent(out) :: stride(:) !! Fortran column-major strides
        integer, intent(out) :: status !! success or invalid input
        integer :: count_active

        status = spatial_status_invalid_input
        if (size(dims) /= size(domain_dims)) return
        if (any(domain_dims < 1)) return
        if (any(mod(domain_dims, 2) /= 1)) return
        count_active = count(mask)
        if (count_active < 1) return
        if (rank < 0 .or. rank >= count_active) return
        call shape_strides(dims, stride)
        allocate(offsets(size(dims), count_active))
        call centered_mask_offsets(domain_dims, mask, offsets)
        status = spatial_status_success
    end subroutine prepare_mask

    subroutine order_real_flat(x, dims, domain_dims, mask, rank, y, status)
        real(dp), intent(in) :: x(:) !! flattened real signal
        integer, intent(in) :: dims(:) !! signal extents
        integer, intent(in) :: domain_dims(:) !! footprint extents
        logical, intent(in) :: mask(:) !! flattened footprint
        integer, intent(in) :: rank !! zero-based selection rank
        real(dp), intent(out) :: y(:) !! flattened result
        integer, intent(out) :: status !! success or invalid input
        integer, allocatable :: offsets(:,:), stride(:)
        real(dp), allocatable :: window(:)
        integer :: i, k, idx

        y = 0.0_dp
        allocate(stride(size(dims)))
        call prepare_mask(dims, domain_dims, mask, rank, offsets, stride, status)
        if (status /= spatial_status_success) return
        allocate(window(size(offsets, 2)))
        do i = 1, size(x)
            do k = 1, size(window)
                call order_index(i, dims, stride, offsets(:,k), idx)
                window(k) = 0.0_dp
                if (idx > 0) window(k) = x(idx)
            end do
            y(i) = kth_real(window, rank)
        end do
    end subroutine order_real_flat

    subroutine order_int_flat(x, dims, domain_dims, mask, rank, y, status)
        integer, intent(in) :: x(:) !! flattened integer signal
        integer, intent(in) :: dims(:) !! signal extents
        integer, intent(in) :: domain_dims(:) !! footprint extents
        logical, intent(in) :: mask(:) !! flattened footprint
        integer, intent(in) :: rank !! zero-based selection rank
        integer, intent(out) :: y(:) !! flattened result
        integer, intent(out) :: status !! success or invalid input
        integer, allocatable :: offsets(:,:), stride(:)
        integer, allocatable :: window(:)
        integer :: i, k, idx

        y = 0
        allocate(stride(size(dims)))
        call prepare_mask(dims, domain_dims, mask, rank, offsets, stride, status)
        if (status /= spatial_status_success) return
        allocate(window(size(offsets, 2)))
        do i = 1, size(x)
            do k = 1, size(window)
                call order_index(i, dims, stride, offsets(:,k), idx)
                window(k) = 0
                if (idx > 0) window(k) = x(idx)
            end do
            y(i) = kth_int(window, rank)
        end do
    end subroutine order_int_flat

    ! SciPy sepfir2d uses edge-repeating (half-sample symmetric) reflection.
    ! The modulo also handles filters longer than the input dimension.
    integer function mirror_index(i, n) result(j)
        integer, intent(in) :: i !! possibly outside one-based index
        integer, intent(in) :: n !! nonempty dimension length
        integer :: q

        q = modulo(i - 1, 2*n)
        if (q < n) then
            j = q + 1
        else
            j = 2*n - q
        end if
    end function mirror_index

    function sepfir2d_real_core(x, hrow, hcol, status) result(y)
        real(dp), intent(in) :: x(:,:) !! real input image
        real(dp), intent(in) :: hrow(:) !! odd horizontal FIR taps
        real(dp), intent(in) :: hcol(:) !! odd vertical FIR taps
        integer, intent(out) :: status !! zero on success
        real(dp), allocatable :: y(:,:) !! filtered image
        real(dp), allocatable :: temp(:,:)
        integer :: i, j, k, src, nr, nc

        nr = size(x,1)
        nc = size(x,2)
        allocate(y(nr,nc))
        y = 0.0_dp
        status = spatial_status_invalid_input
        if (nr < 1 .or. nc < 1) return
        if (size(hrow) < 1 .or. size(hcol) < 1) return
        if (mod(size(hrow),2) /= 1 .or. mod(size(hcol),2) /= 1) return
        allocate(temp(nr,nc))
        temp = 0.0_dp
        do j = 1, nc
            do i = 1, nr
                do k = 1, size(hrow)
                    src = mirror_index(j + size(hrow)/2 + 1 - k, nc)
                    temp(i,j) = temp(i,j) + x(i,src)*hrow(k)
                end do
            end do
        end do
        do j = 1, nc
            do i = 1, nr
                do k = 1, size(hcol)
                    src = mirror_index(i + size(hcol)/2 + 1 - k, nr)
                    y(i,j) = y(i,j) + temp(src,j)*hcol(k)
                end do
            end do
        end do
        status = spatial_status_success
    end function sepfir2d_real_core

    function sepfir2d_complex_core(x, hrow, hcol, status) result(y)
        complex(dp), intent(in) :: x(:,:) !! complex input image
        complex(dp), intent(in) :: hrow(:) !! odd horizontal FIR taps
        complex(dp), intent(in) :: hcol(:) !! odd vertical FIR taps
        integer, intent(out) :: status !! zero on success
        complex(dp), allocatable :: y(:,:) !! filtered image
        complex(dp), allocatable :: temp(:,:)
        integer :: i, j, k, src, nr, nc

        nr = size(x,1)
        nc = size(x,2)
        allocate(y(nr,nc))
        y = cmplx(0.0_dp,0.0_dp,dp)
        status = spatial_status_invalid_input
        if (nr < 1 .or. nc < 1) return
        if (size(hrow) < 1 .or. size(hcol) < 1) return
        if (mod(size(hrow),2) /= 1 .or. mod(size(hcol),2) /= 1) return
        allocate(temp(nr,nc))
        temp = cmplx(0.0_dp,0.0_dp,dp)
        do j = 1, nc
            do i = 1, nr
                do k = 1, size(hrow)
                    src = mirror_index(j + size(hrow)/2 + 1 - k, nc)
                    temp(i,j) = temp(i,j) + x(i,src)*hrow(k)
                end do
            end do
        end do
        do j = 1, nc
            do i = 1, nr
                do k = 1, size(hcol)
                    src = mirror_index(i + size(hcol)/2 + 1 - k, nr)
                    y(i,j) = y(i,j) + temp(src,j)*hcol(k)
                end do
            end do
        end do
        status = spatial_status_success
    end function sepfir2d_complex_core

    function order_real_int_1d(x, domain, rank, status) result(y)
        real(dp), intent(in) :: x(:) !! input signal
        integer, intent(in) :: domain(:) !! nonzero/true domain footprint
        integer, intent(in) :: rank !! zero-based order statistic
        integer, intent(out), optional :: status !! success or invalid input
        real(dp), allocatable :: y(:) !! filtered output
        real(dp), allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call order_real_flat(reshape(x, [size(x)]), shape(x), shape(domain), &
            reshape(domain /= 0, [size(domain)]), rank, flat, code)
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function order_real_int_1d

    function order_real_int_2d(x, domain, rank, status) result(y)
        real(dp), intent(in) :: x(:,:) !! input signal
        integer, intent(in) :: domain(:,:) !! nonzero/true domain footprint
        integer, intent(in) :: rank !! zero-based order statistic
        integer, intent(out), optional :: status !! success or invalid input
        real(dp), allocatable :: y(:,:) !! filtered output
        real(dp), allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call order_real_flat(reshape(x, [size(x)]), shape(x), shape(domain), &
            reshape(domain /= 0, [size(domain)]), rank, flat, code)
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function order_real_int_2d

    function order_real_int_3d(x, domain, rank, status) result(y)
        real(dp), intent(in) :: x(:,:,:) !! input signal
        integer, intent(in) :: domain(:,:,:) !! nonzero/true domain footprint
        integer, intent(in) :: rank !! zero-based order statistic
        integer, intent(out), optional :: status !! success or invalid input
        real(dp), allocatable :: y(:,:,:) !! filtered output
        real(dp), allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call order_real_flat(reshape(x, [size(x)]), shape(x), shape(domain), &
            reshape(domain /= 0, [size(domain)]), rank, flat, code)
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function order_real_int_3d

    function order_real_real_1d(x, domain, rank, status) result(y)
        real(dp), intent(in) :: x(:) !! input signal
        real(dp), intent(in) :: domain(:) !! nonzero/true domain footprint
        integer, intent(in) :: rank !! zero-based order statistic
        integer, intent(out), optional :: status !! success or invalid input
        real(dp), allocatable :: y(:) !! filtered output
        real(dp), allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call order_real_flat(reshape(x, [size(x)]), shape(x), shape(domain), &
            reshape(domain /= 0.0_dp, [size(domain)]), rank, flat, code)
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function order_real_real_1d

    function order_real_real_2d(x, domain, rank, status) result(y)
        real(dp), intent(in) :: x(:,:) !! input signal
        real(dp), intent(in) :: domain(:,:) !! nonzero/true domain footprint
        integer, intent(in) :: rank !! zero-based order statistic
        integer, intent(out), optional :: status !! success or invalid input
        real(dp), allocatable :: y(:,:) !! filtered output
        real(dp), allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call order_real_flat(reshape(x, [size(x)]), shape(x), shape(domain), &
            reshape(domain /= 0.0_dp, [size(domain)]), rank, flat, code)
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function order_real_real_2d

    function order_real_real_3d(x, domain, rank, status) result(y)
        real(dp), intent(in) :: x(:,:,:) !! input signal
        real(dp), intent(in) :: domain(:,:,:) !! nonzero/true domain footprint
        integer, intent(in) :: rank !! zero-based order statistic
        integer, intent(out), optional :: status !! success or invalid input
        real(dp), allocatable :: y(:,:,:) !! filtered output
        real(dp), allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call order_real_flat(reshape(x, [size(x)]), shape(x), shape(domain), &
            reshape(domain /= 0.0_dp, [size(domain)]), rank, flat, code)
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function order_real_real_3d

    function order_real_logical_1d(x, domain, rank, status) result(y)
        real(dp), intent(in) :: x(:) !! input signal
        logical, intent(in) :: domain(:) !! nonzero/true domain footprint
        integer, intent(in) :: rank !! zero-based order statistic
        integer, intent(out), optional :: status !! success or invalid input
        real(dp), allocatable :: y(:) !! filtered output
        real(dp), allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call order_real_flat(reshape(x, [size(x)]), shape(x), shape(domain), &
            reshape(domain, [size(domain)]), rank, flat, code)
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function order_real_logical_1d

    function order_real_logical_2d(x, domain, rank, status) result(y)
        real(dp), intent(in) :: x(:,:) !! input signal
        logical, intent(in) :: domain(:,:) !! nonzero/true domain footprint
        integer, intent(in) :: rank !! zero-based order statistic
        integer, intent(out), optional :: status !! success or invalid input
        real(dp), allocatable :: y(:,:) !! filtered output
        real(dp), allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call order_real_flat(reshape(x, [size(x)]), shape(x), shape(domain), &
            reshape(domain, [size(domain)]), rank, flat, code)
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function order_real_logical_2d

    function order_real_logical_3d(x, domain, rank, status) result(y)
        real(dp), intent(in) :: x(:,:,:) !! input signal
        logical, intent(in) :: domain(:,:,:) !! nonzero/true domain footprint
        integer, intent(in) :: rank !! zero-based order statistic
        integer, intent(out), optional :: status !! success or invalid input
        real(dp), allocatable :: y(:,:,:) !! filtered output
        real(dp), allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call order_real_flat(reshape(x, [size(x)]), shape(x), shape(domain), &
            reshape(domain, [size(domain)]), rank, flat, code)
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function order_real_logical_3d

    function order_int_int_1d(x, domain, rank, status) result(y)
        integer, intent(in) :: x(:) !! input signal
        integer, intent(in) :: domain(:) !! nonzero/true domain footprint
        integer, intent(in) :: rank !! zero-based order statistic
        integer, intent(out), optional :: status !! success or invalid input
        integer, allocatable :: y(:) !! filtered output
        integer, allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call order_int_flat(reshape(x, [size(x)]), shape(x), shape(domain), &
            reshape(domain /= 0, [size(domain)]), rank, flat, code)
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function order_int_int_1d

    function order_int_int_2d(x, domain, rank, status) result(y)
        integer, intent(in) :: x(:,:) !! input signal
        integer, intent(in) :: domain(:,:) !! nonzero/true domain footprint
        integer, intent(in) :: rank !! zero-based order statistic
        integer, intent(out), optional :: status !! success or invalid input
        integer, allocatable :: y(:,:) !! filtered output
        integer, allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call order_int_flat(reshape(x, [size(x)]), shape(x), shape(domain), &
            reshape(domain /= 0, [size(domain)]), rank, flat, code)
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function order_int_int_2d

    function order_int_int_3d(x, domain, rank, status) result(y)
        integer, intent(in) :: x(:,:,:) !! input signal
        integer, intent(in) :: domain(:,:,:) !! nonzero/true domain footprint
        integer, intent(in) :: rank !! zero-based order statistic
        integer, intent(out), optional :: status !! success or invalid input
        integer, allocatable :: y(:,:,:) !! filtered output
        integer, allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call order_int_flat(reshape(x, [size(x)]), shape(x), shape(domain), &
            reshape(domain /= 0, [size(domain)]), rank, flat, code)
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function order_int_int_3d

    function order_int_real_1d(x, domain, rank, status) result(y)
        integer, intent(in) :: x(:) !! input signal
        real(dp), intent(in) :: domain(:) !! nonzero/true domain footprint
        integer, intent(in) :: rank !! zero-based order statistic
        integer, intent(out), optional :: status !! success or invalid input
        integer, allocatable :: y(:) !! filtered output
        integer, allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call order_int_flat(reshape(x, [size(x)]), shape(x), shape(domain), &
            reshape(domain /= 0.0_dp, [size(domain)]), rank, flat, code)
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function order_int_real_1d

    function order_int_real_2d(x, domain, rank, status) result(y)
        integer, intent(in) :: x(:,:) !! input signal
        real(dp), intent(in) :: domain(:,:) !! nonzero/true domain footprint
        integer, intent(in) :: rank !! zero-based order statistic
        integer, intent(out), optional :: status !! success or invalid input
        integer, allocatable :: y(:,:) !! filtered output
        integer, allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call order_int_flat(reshape(x, [size(x)]), shape(x), shape(domain), &
            reshape(domain /= 0.0_dp, [size(domain)]), rank, flat, code)
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function order_int_real_2d

    function order_int_real_3d(x, domain, rank, status) result(y)
        integer, intent(in) :: x(:,:,:) !! input signal
        real(dp), intent(in) :: domain(:,:,:) !! nonzero/true domain footprint
        integer, intent(in) :: rank !! zero-based order statistic
        integer, intent(out), optional :: status !! success or invalid input
        integer, allocatable :: y(:,:,:) !! filtered output
        integer, allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call order_int_flat(reshape(x, [size(x)]), shape(x), shape(domain), &
            reshape(domain /= 0.0_dp, [size(domain)]), rank, flat, code)
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function order_int_real_3d

    function order_int_logical_1d(x, domain, rank, status) result(y)
        integer, intent(in) :: x(:) !! input signal
        logical, intent(in) :: domain(:) !! nonzero/true domain footprint
        integer, intent(in) :: rank !! zero-based order statistic
        integer, intent(out), optional :: status !! success or invalid input
        integer, allocatable :: y(:) !! filtered output
        integer, allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call order_int_flat(reshape(x, [size(x)]), shape(x), shape(domain), &
            reshape(domain, [size(domain)]), rank, flat, code)
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function order_int_logical_1d

    function order_int_logical_2d(x, domain, rank, status) result(y)
        integer, intent(in) :: x(:,:) !! input signal
        logical, intent(in) :: domain(:,:) !! nonzero/true domain footprint
        integer, intent(in) :: rank !! zero-based order statistic
        integer, intent(out), optional :: status !! success or invalid input
        integer, allocatable :: y(:,:) !! filtered output
        integer, allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call order_int_flat(reshape(x, [size(x)]), shape(x), shape(domain), &
            reshape(domain, [size(domain)]), rank, flat, code)
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function order_int_logical_2d

    function order_int_logical_3d(x, domain, rank, status) result(y)
        integer, intent(in) :: x(:,:,:) !! input signal
        logical, intent(in) :: domain(:,:,:) !! nonzero/true domain footprint
        integer, intent(in) :: rank !! zero-based order statistic
        integer, intent(out), optional :: status !! success or invalid input
        integer, allocatable :: y(:,:,:) !! filtered output
        integer, allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call order_int_flat(reshape(x, [size(x)]), shape(x), shape(domain), &
            reshape(domain, [size(domain)]), rank, flat, code)
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function order_int_logical_3d

    function sep_real_real_real(x, hrow, hcol, status) result(y)
        real(dp), intent(in) :: x(:,:) !! input image
        real(dp), intent(in) :: hrow(:) !! horizontal FIR coefficients
        real(dp), intent(in) :: hcol(:) !! vertical FIR coefficients
        integer, intent(out), optional :: status !! success or invalid input
        real(dp), allocatable :: y(:,:) !! filtered output
        integer :: code

        y = sepfir2d_real_core(x, &
            hrow, hcol, code)
        if (present(status)) status = code
    end function sep_real_real_real

    function sep_real_real_complex(x, hrow, hcol, status) result(y)
        real(dp), intent(in) :: x(:,:) !! input image
        real(dp), intent(in) :: hrow(:) !! horizontal FIR coefficients
        complex(dp), intent(in) :: hcol(:) !! vertical FIR coefficients
        integer, intent(out), optional :: status !! success or invalid input
        complex(dp), allocatable :: y(:,:) !! filtered output
        integer :: code

        y = sepfir2d_complex_core(cmplx(x,kind=dp), &
            cmplx(hrow,kind=dp), hcol, code)
        if (present(status)) status = code
    end function sep_real_real_complex

    function sep_real_complex_real(x, hrow, hcol, status) result(y)
        real(dp), intent(in) :: x(:,:) !! input image
        complex(dp), intent(in) :: hrow(:) !! horizontal FIR coefficients
        real(dp), intent(in) :: hcol(:) !! vertical FIR coefficients
        integer, intent(out), optional :: status !! success or invalid input
        complex(dp), allocatable :: y(:,:) !! filtered output
        integer :: code

        y = sepfir2d_complex_core(cmplx(x,kind=dp), &
            hrow, cmplx(hcol,kind=dp), code)
        if (present(status)) status = code
    end function sep_real_complex_real

    function sep_real_complex_complex(x, hrow, hcol, status) result(y)
        real(dp), intent(in) :: x(:,:) !! input image
        complex(dp), intent(in) :: hrow(:) !! horizontal FIR coefficients
        complex(dp), intent(in) :: hcol(:) !! vertical FIR coefficients
        integer, intent(out), optional :: status !! success or invalid input
        complex(dp), allocatable :: y(:,:) !! filtered output
        integer :: code

        y = sepfir2d_complex_core(cmplx(x,kind=dp), &
            hrow, hcol, code)
        if (present(status)) status = code
    end function sep_real_complex_complex

    function sep_complex_real_real(x, hrow, hcol, status) result(y)
        complex(dp), intent(in) :: x(:,:) !! input image
        real(dp), intent(in) :: hrow(:) !! horizontal FIR coefficients
        real(dp), intent(in) :: hcol(:) !! vertical FIR coefficients
        integer, intent(out), optional :: status !! success or invalid input
        complex(dp), allocatable :: y(:,:) !! filtered output
        integer :: code

        y = sepfir2d_complex_core(x, &
            cmplx(hrow,kind=dp), cmplx(hcol,kind=dp), code)
        if (present(status)) status = code
    end function sep_complex_real_real

    function sep_complex_real_complex(x, hrow, hcol, status) result(y)
        complex(dp), intent(in) :: x(:,:) !! input image
        real(dp), intent(in) :: hrow(:) !! horizontal FIR coefficients
        complex(dp), intent(in) :: hcol(:) !! vertical FIR coefficients
        integer, intent(out), optional :: status !! success or invalid input
        complex(dp), allocatable :: y(:,:) !! filtered output
        integer :: code

        y = sepfir2d_complex_core(x, &
            cmplx(hrow,kind=dp), hcol, code)
        if (present(status)) status = code
    end function sep_complex_real_complex

    function sep_complex_complex_real(x, hrow, hcol, status) result(y)
        complex(dp), intent(in) :: x(:,:) !! input image
        complex(dp), intent(in) :: hrow(:) !! horizontal FIR coefficients
        real(dp), intent(in) :: hcol(:) !! vertical FIR coefficients
        integer, intent(out), optional :: status !! success or invalid input
        complex(dp), allocatable :: y(:,:) !! filtered output
        integer :: code

        y = sepfir2d_complex_core(x, &
            hrow, cmplx(hcol,kind=dp), code)
        if (present(status)) status = code
    end function sep_complex_complex_real

    function sep_complex_complex_complex(x, hrow, hcol, status) result(y)
        complex(dp), intent(in) :: x(:,:) !! input image
        complex(dp), intent(in) :: hrow(:) !! horizontal FIR coefficients
        complex(dp), intent(in) :: hcol(:) !! vertical FIR coefficients
        integer, intent(out), optional :: status !! success or invalid input
        complex(dp), allocatable :: y(:,:) !! filtered output
        integer :: code

        y = sepfir2d_complex_core(x, &
            hrow, hcol, code)
        if (present(status)) status = code
    end function sep_complex_complex_complex

    function sep_int_real_real(x, hrow, hcol, status) result(y)
        integer, intent(in) :: x(:,:) !! input image
        real(dp), intent(in) :: hrow(:) !! horizontal FIR coefficients
        real(dp), intent(in) :: hcol(:) !! vertical FIR coefficients
        integer, intent(out), optional :: status !! success or invalid input
        real(dp), allocatable :: y(:,:) !! filtered output
        integer :: code

        y = sepfir2d_real_core(real(x,dp), &
            hrow, hcol, code)
        if (present(status)) status = code
    end function sep_int_real_real

    function sep_int_real_complex(x, hrow, hcol, status) result(y)
        integer, intent(in) :: x(:,:) !! input image
        real(dp), intent(in) :: hrow(:) !! horizontal FIR coefficients
        complex(dp), intent(in) :: hcol(:) !! vertical FIR coefficients
        integer, intent(out), optional :: status !! success or invalid input
        complex(dp), allocatable :: y(:,:) !! filtered output
        integer :: code

        y = sepfir2d_complex_core(cmplx(real(x,dp),kind=dp), &
            cmplx(hrow,kind=dp), hcol, code)
        if (present(status)) status = code
    end function sep_int_real_complex

    function sep_int_complex_real(x, hrow, hcol, status) result(y)
        integer, intent(in) :: x(:,:) !! input image
        complex(dp), intent(in) :: hrow(:) !! horizontal FIR coefficients
        real(dp), intent(in) :: hcol(:) !! vertical FIR coefficients
        integer, intent(out), optional :: status !! success or invalid input
        complex(dp), allocatable :: y(:,:) !! filtered output
        integer :: code

        y = sepfir2d_complex_core(cmplx(real(x,dp),kind=dp), &
            hrow, cmplx(hcol,kind=dp), code)
        if (present(status)) status = code
    end function sep_int_complex_real

    function sep_int_complex_complex(x, hrow, hcol, status) result(y)
        integer, intent(in) :: x(:,:) !! input image
        complex(dp), intent(in) :: hrow(:) !! horizontal FIR coefficients
        complex(dp), intent(in) :: hcol(:) !! vertical FIR coefficients
        integer, intent(out), optional :: status !! success or invalid input
        complex(dp), allocatable :: y(:,:) !! filtered output
        integer :: code

        y = sepfir2d_complex_core(cmplx(real(x,dp),kind=dp), &
            hrow, hcol, code)
        if (present(status)) status = code
    end function sep_int_complex_complex

end module scifort_signal_spatial
