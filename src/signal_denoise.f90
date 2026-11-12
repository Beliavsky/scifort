! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Zero-padded median and local Wiener filters for ranks one through three.
! The implementation is independent of SciPy, with SciPy used as a reference.
module scifort_signal_denoise
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: denoise_status_success = 0
    integer, parameter, public :: denoise_status_invalid_input = 1

    public :: medfilt, medfilt2d, wiener

    interface medfilt
        module procedure medfilt_real_1d, medfilt_real_2d, medfilt_real_3d
        module procedure medfilt_int_1d, medfilt_int_2d, medfilt_int_3d
    end interface
    interface medfilt2d
        module procedure medfilt_real_2d, medfilt_int_2d
    end interface
    interface wiener
        module procedure wiener_real_1d, wiener_real_2d, wiener_real_3d
        module procedure wiener_int_1d, wiener_int_2d, wiener_int_3d
        module procedure wiener_complex_1d, wiener_complex_2d, wiener_complex_3d
        module procedure wiener_complex_noise_1d, wiener_complex_noise_2d, wiener_complex_noise_3d
    end interface

contains

    subroutine kernel_dimensions(ndim, kernel_size, kernel_shape, odd, ks, nwin, status)
        integer, intent(in) :: ndim !! number of signal dimensions
        integer, intent(in), optional :: kernel_size !! same window size on every axis
        integer, intent(in), optional :: kernel_shape(:) !! window size for each axis
        logical, intent(in) :: odd !! whether all window lengths must be odd
        integer, intent(out) :: ks(ndim) !! validated kernel dimensions
        integer, intent(out) :: nwin !! number of elements in a window
        integer, intent(out) :: status !! zero on success

        integer :: d
        integer(int64) :: count

        status = denoise_status_invalid_input
        nwin = 0
        ks = 0
        if (present(kernel_shape) .and. present(kernel_size)) return
        if (present(kernel_shape)) then
            if (size(kernel_shape) /= ndim) return
            ks = kernel_shape
        else if (present(kernel_size)) then
            ks = kernel_size
        else
            ks = 3
        end if
        if (any(ks <= 0)) return
        if (odd .and. any(mod(ks, 2) == 0)) return
        count = 1_int64
        do d = 1, ndim
            if (count > int(huge(nwin), int64) / int(ks(d), int64)) return
            count = count * int(ks(d), int64)
        end do
        nwin = int(count)
        status = denoise_status_success
    end subroutine kernel_dimensions

    subroutine window_offsets(ks, offsets)
        integer, intent(in) :: ks(:) !! kernel extents
        integer, intent(out) :: offsets(:,:) !! offset of every sample in each window

        integer :: i, j, q

        do i = 1, size(offsets, 2)
            q = i - 1
            do j = 1, size(ks)
                offsets(j, i) = mod(q, ks(j)) - ks(j) / 2
                q = q / ks(j)
            end do
        end do
    end subroutine window_offsets

    subroutine linear_coordinates(i, dims, coords)
        integer, intent(in) :: i !! one-based flat position
        integer, intent(in) :: dims(:) !! extents
        integer, intent(out) :: coords(:) !! one-based multi-index
        integer :: q, d

        q = i - 1
        do d = 1, size(dims)
            coords(d) = mod(q, dims(d)) + 1
            q = q / dims(d)
        end do
    end subroutine linear_coordinates

    integer function neighbor_index(coords, dims, stride, offsets) result(idx)
        integer, intent(in) :: coords(:) !! center coordinates
        integer, intent(in) :: dims(:) !! source extents
        integer, intent(in) :: stride(:) !! Fortran linear strides
        integer, intent(in) :: offsets(:) !! displacement relative to center
        integer :: d, q

        idx = 1
        do d = 1, size(dims)
            q = coords(d) + offsets(d)
            if (q < 1 .or. q > dims(d)) then
                idx = 0
                return
            end if
            idx = idx + (q - 1) * stride(d)
        end do
    end function neighbor_index

    subroutine dimension_strides(dims, stride)
        integer, intent(in) :: dims(:) !! extents
        integer, intent(out) :: stride(:) !! column-major strides
        integer :: d

        stride(1) = 1
        do d = 2, size(dims)
            stride(d) = stride(d - 1) * dims(d - 1)
        end do
    end subroutine dimension_strides

    real(dp) function select_median(work) result(value)
        real(dp), intent(inout) :: work(:) !! modified in place
        real(dp) :: pivot, tmp
        integer :: lo, hi, lower, upper, i, target

        target = (size(work) + 1) / 2
        lo = 1
        hi = size(work)
        do while (lo < hi)
            pivot = work(lo + (hi - lo) / 2)
            lower = lo
            upper = hi
            i = lo
            do while (i <= upper)
                if (work(i) < pivot) then
                    tmp = work(lower)
                    work(lower) = work(i)
                    work(i) = tmp
                    lower = lower + 1
                    i = i + 1
                else if (work(i) > pivot) then
                    tmp = work(upper)
                    work(upper) = work(i)
                    work(i) = tmp
                    upper = upper - 1
                else
                    i = i + 1
                end if
            end do
            if (target < lower) then
                hi = lower - 1
            else if (target > upper) then
                lo = upper + 1
            else
                exit
            end if
        end do
        value = work(target)
    end function select_median

    subroutine medfilt_flat(x, dims, kernel_size, kernel_shape, y, status)
        real(dp), intent(in) :: x(:) !! flattened input
        integer, intent(in) :: dims(:) !! original Fortran array shape
        integer, intent(in), optional :: kernel_size !! scalar odd window length
        integer, intent(in), optional :: kernel_shape(:) !! per-axis odd window lengths
        real(dp), intent(out) :: y(:) !! flattened result
        integer, intent(out) :: status !! zero on success

        integer :: nd, nwin, i, j, idx
        integer, allocatable :: ks(:), offset(:,:), coord(:), stride(:)
        real(dp), allocatable :: window(:)

        y = 0.0_dp
        nd = size(dims)
        allocate(ks(nd), coord(nd), stride(nd))
        call kernel_dimensions(nd, kernel_size, kernel_shape, .true., ks, nwin, status)
        if (status /= denoise_status_success) return
        allocate(offset(nd, nwin), window(nwin))
        call dimension_strides(dims, stride)
        call window_offsets(ks, offset)
        do i = 1, size(x)
            call linear_coordinates(i, dims, coord)
            do j = 1, nwin
                idx = neighbor_index(coord, dims, stride, offset(:, j))
                window(j) = 0.0_dp
                if (idx > 0) window(j) = x(idx)
            end do
            y(i) = select_median(window)
        end do
    end subroutine medfilt_flat

    subroutine wiener_flat(x, dims, kernel_size, kernel_shape, noise, y, status)
        real(dp), intent(in) :: x(:) !! flattened input
        integer, intent(in) :: dims(:) !! original Fortran array shape
        integer, intent(in), optional :: kernel_size !! scalar window length
        integer, intent(in), optional :: kernel_shape(:) !! per-axis window lengths
        real(dp), intent(in), optional :: noise !! supplied noise power
        real(dp), intent(out) :: y(:) !! flattened filtered output
        integer, intent(out) :: status !! zero on success

        integer :: nd, nwin, i, j, idx
        integer, allocatable :: ks(:), offset(:,:), coord(:), stride(:)
        real(dp), allocatable :: means(:), variances(:)
        real(dp) :: mean_value, second_moment, sample, noise_power

        y = 0.0_dp
        nd = size(dims)
        allocate(ks(nd), coord(nd), stride(nd))
        call kernel_dimensions(nd, kernel_size, kernel_shape, .false., ks, nwin, status)
        if (status /= denoise_status_success) return
        if (size(x) == 0) return
        allocate(offset(nd, nwin), means(size(x)), variances(size(x)))
        call dimension_strides(dims, stride)
        call window_offsets(ks, offset)
        do i = 1, size(x)
            call linear_coordinates(i, dims, coord)
            mean_value = 0.0_dp
            second_moment = 0.0_dp
            do j = 1, nwin
                idx = neighbor_index(coord, dims, stride, offset(:, j))
                sample = 0.0_dp
                if (idx > 0) sample = x(idx)
                mean_value = mean_value + sample
                second_moment = second_moment + sample * sample
            end do
            means(i) = mean_value / real(nwin, dp)
            variances(i) = second_moment / real(nwin, dp) - means(i) * means(i)
        end do
        if (present(noise)) then
            noise_power = noise
        else
            noise_power = sum(variances) / real(size(x), dp)
        end if
        do i = 1, size(x)
            if (variances(i) < noise_power) then
                y(i) = means(i)
            else if (variances(i) == 0.0_dp) then
                ! SciPy's 0/0 gain is NaN; assign without tripping FPE traps.
                y(i) = ieee_value(0.0_dp, ieee_quiet_nan)
            else
                y(i) = means(i) + (x(i) - means(i)) * (1.0_dp - noise_power / variances(i))
            end if
        end do
    end subroutine wiener_flat

    ! SciPy uses the *pseudo*-variance E[x*x]-E[x]**2 for complex Wiener
    ! filtering, rather than E[abs(x)**2]-abs(E[x])**2. Both the local
    ! variance and the automatic noise estimate are therefore complex.
    ! NumPy compares complex values lexicographically (real part first).
    pure logical function complex_less(a, b) result(less)
        complex(dp), intent(in) :: a !! local pseudo-variance
        complex(dp), intent(in) :: b !! estimated or supplied noise

        less = real(a, dp) < real(b, dp)
        if (real(a, dp) == real(b, dp)) less = aimag(a) < aimag(b)
    end function complex_less

    subroutine wiener_complex_flat(x, dims, kernel_size, kernel_shape, noise, y, status)
        complex(dp), intent(in) :: x(:) !! flattened complex input
        integer, intent(in) :: dims(:) !! original signal dimensions
        integer, intent(in), optional :: kernel_size !! scalar window length, default three
        integer, intent(in), optional :: kernel_shape(:) !! per-axis window lengths
        complex(dp), intent(in), optional :: noise !! optional complex noise pseudo-variance
        complex(dp), intent(out) :: y(:) !! flattened filtered output
        integer, intent(out) :: status !! success or invalid input

        integer :: nd, nwin, i, j, idx
        integer, allocatable :: ks(:), offset(:,:), coord(:), stride(:)
        complex(dp), allocatable :: means(:), variances(:)
        complex(dp) :: mean_value, second_moment, sample, noise_power
        real(dp) :: nan

        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        y = cmplx(nan, nan, dp)
        status = denoise_status_invalid_input
        if (.not. all(ieee_is_finite(real(x, dp)))) return
        if (.not. all(ieee_is_finite(aimag(x)))) return
        if (present(noise)) then
            if (.not. ieee_is_finite(real(noise, dp))) return
            if (.not. ieee_is_finite(aimag(noise))) return
        end if
        nd = size(dims)
        allocate(ks(nd), coord(nd), stride(nd))
        call kernel_dimensions(nd, kernel_size, kernel_shape, .false., ks, nwin, status)
        if (status /= denoise_status_success) return
        if (size(x) == 0) return
        allocate(offset(nd, nwin), means(size(x)), variances(size(x)))
        call dimension_strides(dims, stride)
        call window_offsets(ks, offset)
        do i = 1, size(x)
            call linear_coordinates(i, dims, coord)
            mean_value = (0.0_dp, 0.0_dp)
            second_moment = (0.0_dp, 0.0_dp)
            do j = 1, nwin
                idx = neighbor_index(coord, dims, stride, offset(:, j))
                sample = (0.0_dp, 0.0_dp)
                if (idx > 0) sample = x(idx)
                mean_value = mean_value + sample
                second_moment = second_moment + sample * sample
            end do
            means(i) = mean_value / real(nwin, dp)
            variances(i) = second_moment / real(nwin, dp) - means(i)*means(i)
        end do
        if (present(noise)) then
            noise_power = noise
        else
            noise_power = sum(variances) / real(size(x), dp)
        end if
        do i = 1, size(x)
            if (complex_less(variances(i), noise_power)) then
                y(i) = means(i)
            else if (abs(variances(i)) == 0.0_dp) then
                y(i) = cmplx(nan, nan, dp)
            else
                y(i) = means(i) + (x(i)-means(i)) * (1.0_dp-noise_power/variances(i))
            end if
        end do
    end subroutine wiener_complex_flat

    function medfilt_real_1d(x, kernel_size, kernel_shape, status) result(y)
        real(dp), intent(in) :: x(:) !! input signal
        integer, intent(in), optional :: kernel_size !! uniform window length (default 3)
        integer, intent(in), optional :: kernel_shape(:) !! per-axis lengths (overrides default)
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: y(:)
        real(dp), allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call medfilt_flat(reshape(x, [size(x)]), &
            shape(x), kernel_size, kernel_shape, flat, code)
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function medfilt_real_1d

    function medfilt_real_2d(x, kernel_size, kernel_shape, status) result(y)
        real(dp), intent(in) :: x(:,:) !! input signal
        integer, intent(in), optional :: kernel_size !! uniform window length (default 3)
        integer, intent(in), optional :: kernel_shape(:) !! per-axis lengths (overrides default)
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: y(:,:)
        real(dp), allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call medfilt_flat(reshape(x, [size(x)]), &
            shape(x), kernel_size, kernel_shape, flat, code)
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function medfilt_real_2d

    function medfilt_real_3d(x, kernel_size, kernel_shape, status) result(y)
        real(dp), intent(in) :: x(:,:,:) !! input signal
        integer, intent(in), optional :: kernel_size !! uniform window length (default 3)
        integer, intent(in), optional :: kernel_shape(:) !! per-axis lengths (overrides default)
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: y(:,:,:)
        real(dp), allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call medfilt_flat(reshape(x, [size(x)]), &
            shape(x), kernel_size, kernel_shape, flat, code)
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function medfilt_real_3d

    function medfilt_int_1d(x, kernel_size, kernel_shape, status) result(y)
        integer, intent(in) :: x(:) !! input signal
        integer, intent(in), optional :: kernel_size !! uniform window length (default 3)
        integer, intent(in), optional :: kernel_shape(:) !! per-axis lengths (overrides default)
        integer, intent(out), optional :: status !! zero on success
        integer, allocatable :: y(:)
        real(dp), allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call medfilt_flat(real(reshape(x, [size(x)]), dp), &
            shape(x), kernel_size, kernel_shape, flat, code)
        y = reshape(nint(flat), shape(x))
        if (present(status)) status = code
    end function medfilt_int_1d

    function medfilt_int_2d(x, kernel_size, kernel_shape, status) result(y)
        integer, intent(in) :: x(:,:) !! input signal
        integer, intent(in), optional :: kernel_size !! uniform window length (default 3)
        integer, intent(in), optional :: kernel_shape(:) !! per-axis lengths (overrides default)
        integer, intent(out), optional :: status !! zero on success
        integer, allocatable :: y(:,:)
        real(dp), allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call medfilt_flat(real(reshape(x, [size(x)]), dp), &
            shape(x), kernel_size, kernel_shape, flat, code)
        y = reshape(nint(flat), shape(x))
        if (present(status)) status = code
    end function medfilt_int_2d

    function medfilt_int_3d(x, kernel_size, kernel_shape, status) result(y)
        integer, intent(in) :: x(:,:,:) !! input signal
        integer, intent(in), optional :: kernel_size !! uniform window length (default 3)
        integer, intent(in), optional :: kernel_shape(:) !! per-axis lengths (overrides default)
        integer, intent(out), optional :: status !! zero on success
        integer, allocatable :: y(:,:,:)
        real(dp), allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call medfilt_flat(real(reshape(x, [size(x)]), dp), &
            shape(x), kernel_size, kernel_shape, flat, code)
        y = reshape(nint(flat), shape(x))
        if (present(status)) status = code
    end function medfilt_int_3d

    function wiener_real_1d(x, kernel_size, kernel_shape, noise, status) result(y)
        real(dp), intent(in) :: x(:) !! input signal
        integer, intent(in), optional :: kernel_size !! uniform window length (default 3)
        integer, intent(in), optional :: kernel_shape(:) !! per-axis lengths (overrides default)
        real(dp), intent(in), optional :: noise !! noise power; estimated if absent
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: y(:)
        real(dp), allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call wiener_flat(reshape(x, [size(x)]), &
            shape(x), kernel_size, kernel_shape, noise, flat, code)
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function wiener_real_1d

    function wiener_real_2d(x, kernel_size, kernel_shape, noise, status) result(y)
        real(dp), intent(in) :: x(:,:) !! input signal
        integer, intent(in), optional :: kernel_size !! uniform window length (default 3)
        integer, intent(in), optional :: kernel_shape(:) !! per-axis lengths (overrides default)
        real(dp), intent(in), optional :: noise !! noise power; estimated if absent
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: y(:,:)
        real(dp), allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call wiener_flat(reshape(x, [size(x)]), &
            shape(x), kernel_size, kernel_shape, noise, flat, code)
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function wiener_real_2d

    function wiener_real_3d(x, kernel_size, kernel_shape, noise, status) result(y)
        real(dp), intent(in) :: x(:,:,:) !! input signal
        integer, intent(in), optional :: kernel_size !! uniform window length (default 3)
        integer, intent(in), optional :: kernel_shape(:) !! per-axis lengths (overrides default)
        real(dp), intent(in), optional :: noise !! noise power; estimated if absent
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: y(:,:,:)
        real(dp), allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call wiener_flat(reshape(x, [size(x)]), &
            shape(x), kernel_size, kernel_shape, noise, flat, code)
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function wiener_real_3d

    function wiener_int_1d(x, kernel_size, kernel_shape, noise, status) result(y)
        integer, intent(in) :: x(:) !! input signal
        integer, intent(in), optional :: kernel_size !! uniform window length (default 3)
        integer, intent(in), optional :: kernel_shape(:) !! per-axis lengths (overrides default)
        real(dp), intent(in), optional :: noise !! noise power; estimated if absent
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: y(:)
        real(dp), allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call wiener_flat(real(reshape(x, [size(x)]), dp), &
            shape(x), kernel_size, kernel_shape, noise, flat, code)
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function wiener_int_1d

    function wiener_int_2d(x, kernel_size, kernel_shape, noise, status) result(y)
        integer, intent(in) :: x(:,:) !! input signal
        integer, intent(in), optional :: kernel_size !! uniform window length (default 3)
        integer, intent(in), optional :: kernel_shape(:) !! per-axis lengths (overrides default)
        real(dp), intent(in), optional :: noise !! noise power; estimated if absent
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: y(:,:)
        real(dp), allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call wiener_flat(real(reshape(x, [size(x)]), dp), &
            shape(x), kernel_size, kernel_shape, noise, flat, code)
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function wiener_int_2d

    function wiener_int_3d(x, kernel_size, kernel_shape, noise, status) result(y)
        integer, intent(in) :: x(:,:,:) !! input signal
        integer, intent(in), optional :: kernel_size !! uniform window length (default 3)
        integer, intent(in), optional :: kernel_shape(:) !! per-axis lengths (overrides default)
        real(dp), intent(in), optional :: noise !! noise power; estimated if absent
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: y(:,:,:)
        real(dp), allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call wiener_flat(real(reshape(x, [size(x)]), dp), &
            shape(x), kernel_size, kernel_shape, noise, flat, code)
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function wiener_int_3d

    function wiener_complex_1d(x, kernel_size, kernel_shape, noise, status) result(y)
        complex(dp), intent(in) :: x(:) !! complex input signal
        integer, intent(in), optional :: kernel_size !! uniform window length, default three
        integer, intent(in), optional :: kernel_shape(:) !! per-axis window lengths
        real(dp), intent(in), optional :: noise !! user-supplied noise pseudo-variance
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:)
        complex(dp), allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        if (present(noise)) then
            call wiener_complex_flat(reshape(x, [size(x)]), &
                shape(x), kernel_size, kernel_shape, cmplx(noise, 0.0_dp, dp), flat, code)
        else
            call wiener_complex_flat(reshape(x, [size(x)]), &
                shape(x), kernel_size, kernel_shape, y=flat, status=code)
        end if
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function wiener_complex_1d

    function wiener_complex_noise_1d(x, kernel_size, kernel_shape, noise, status) result(y)
        complex(dp), intent(in) :: x(:) !! complex input signal
        integer, intent(in), optional :: kernel_size !! uniform window length, default three
        integer, intent(in), optional :: kernel_shape(:) !! per-axis window lengths
        complex(dp), intent(in) :: noise !! user-supplied noise pseudo-variance
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:)
        complex(dp), allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call wiener_complex_flat(reshape(x, [size(x)]), &
            shape(x), kernel_size, kernel_shape, noise, flat, code)
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function wiener_complex_noise_1d

    function wiener_complex_2d(x, kernel_size, kernel_shape, noise, status) result(y)
        complex(dp), intent(in) :: x(:,:) !! complex input signal
        integer, intent(in), optional :: kernel_size !! uniform window length, default three
        integer, intent(in), optional :: kernel_shape(:) !! per-axis window lengths
        real(dp), intent(in), optional :: noise !! user-supplied noise pseudo-variance
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:,:)
        complex(dp), allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        if (present(noise)) then
            call wiener_complex_flat(reshape(x, [size(x)]), &
                shape(x), kernel_size, kernel_shape, cmplx(noise, 0.0_dp, dp), flat, code)
        else
            call wiener_complex_flat(reshape(x, [size(x)]), &
                shape(x), kernel_size, kernel_shape, y=flat, status=code)
        end if
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function wiener_complex_2d

    function wiener_complex_noise_2d(x, kernel_size, kernel_shape, noise, status) result(y)
        complex(dp), intent(in) :: x(:,:) !! complex input signal
        integer, intent(in), optional :: kernel_size !! uniform window length, default three
        integer, intent(in), optional :: kernel_shape(:) !! per-axis window lengths
        complex(dp), intent(in) :: noise !! user-supplied noise pseudo-variance
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:,:)
        complex(dp), allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call wiener_complex_flat(reshape(x, [size(x)]), &
            shape(x), kernel_size, kernel_shape, noise, flat, code)
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function wiener_complex_noise_2d

    function wiener_complex_3d(x, kernel_size, kernel_shape, noise, status) result(y)
        complex(dp), intent(in) :: x(:,:,:) !! complex input signal
        integer, intent(in), optional :: kernel_size !! uniform window length, default three
        integer, intent(in), optional :: kernel_shape(:) !! per-axis window lengths
        real(dp), intent(in), optional :: noise !! user-supplied noise pseudo-variance
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:,:,:)
        complex(dp), allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        if (present(noise)) then
            call wiener_complex_flat(reshape(x, [size(x)]), &
                shape(x), kernel_size, kernel_shape, cmplx(noise, 0.0_dp, dp), flat, code)
        else
            call wiener_complex_flat(reshape(x, [size(x)]), &
                shape(x), kernel_size, kernel_shape, y=flat, status=code)
        end if
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function wiener_complex_3d

    function wiener_complex_noise_3d(x, kernel_size, kernel_shape, noise, status) result(y)
        complex(dp), intent(in) :: x(:,:,:) !! complex input signal
        integer, intent(in), optional :: kernel_size !! uniform window length, default three
        integer, intent(in), optional :: kernel_shape(:) !! per-axis window lengths
        complex(dp), intent(in) :: noise !! user-supplied noise pseudo-variance
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:,:,:)
        complex(dp), allocatable :: flat(:)
        integer :: code

        allocate(flat(size(x)))
        call wiener_complex_flat(reshape(x, [size(x)]), &
            shape(x), kernel_size, kernel_shape, noise, flat, code)
        y = reshape(flat, shape(x))
        if (present(status)) status = code
    end function wiener_complex_noise_3d


end module scifort_signal_denoise
