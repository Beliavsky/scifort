! SPDX-License-Identifier: MIT
! Driver for independent comparison against SciPy 1.17 find_peaks.
program peak_dynamic_bounds_driver
    use scifort_signal, only : dp, find_peaks, peak_result
    implicit none

    type :: sample_bound
        real(dp), allocatable :: values(:)
    end type sample_bound
    type(sample_bound) :: bounds(10)
    real(dp), allocatable :: x(:), raw(:,:)
    type(peak_result) :: res
    integer :: n, mask, distance, window, i, ios
    real(dp) :: rel_height

    do
        read(*, *, iostat=ios) n, mask, distance, window, rel_height
        if (ios /= 0) exit
        allocate(x(n), raw(n, 10))
        read(*, *) x
        do i = 1, 10
            read(*, *) raw(:, i)
            if (btest(mask, i - 1)) then
                bounds(i)%values = raw(:, i)
            end if
        end do
        if (window <= 0) then
            res = find_peaks(x, distance=distance, rel_height=rel_height, &
                height_min_samples=bounds(1)%values, &
                height_max_samples=bounds(2)%values, &
                threshold_min_samples=bounds(3)%values, &
                threshold_max_samples=bounds(4)%values, &
                prominence_min_samples=bounds(5)%values, &
                prominence_max_samples=bounds(6)%values, &
                width_min_samples=bounds(7)%values, &
                width_max_samples=bounds(8)%values, &
                plateau_size_min_samples=bounds(9)%values, &
                plateau_size_max_samples=bounds(10)%values)
        else
            res = find_peaks(x, distance=distance, wlen=window, rel_height=rel_height, &
                height_min_samples=bounds(1)%values, &
                height_max_samples=bounds(2)%values, &
                threshold_min_samples=bounds(3)%values, &
                threshold_max_samples=bounds(4)%values, &
                prominence_min_samples=bounds(5)%values, &
                prominence_max_samples=bounds(6)%values, &
                width_min_samples=bounds(7)%values, &
                width_max_samples=bounds(8)%values, &
                plateau_size_min_samples=bounds(9)%values, &
                plateau_size_max_samples=bounds(10)%values)
        end if
        print '(i0,1x,i0)', res%status, size(res%peaks)
        print '(*(i0,1x))', res%peaks
        print '(*(es25.16e3,1x))', res%peak_heights
        print '(*(es25.16e3,1x))', res%left_thresholds
        print '(*(es25.16e3,1x))', res%right_thresholds
        print '(*(es25.16e3,1x))', res%prominences
        print '(*(i0,1x))', res%left_bases
        print '(*(i0,1x))', res%right_bases
        print '(*(es25.16e3,1x))', res%widths
        print '(*(es25.16e3,1x))', res%width_heights
        print '(*(es25.16e3,1x))', res%left_ips
        print '(*(es25.16e3,1x))', res%right_ips
        print '(*(i0,1x))', res%plateau_sizes
        print '(*(i0,1x))', res%left_edges
        print '(*(i0,1x))', res%right_edges
        do i = 1, 10
            if (allocated(bounds(i)%values)) deallocate(bounds(i)%values)
        end do
        deallocate(x, raw)
    end do
end program peak_dynamic_bounds_driver
