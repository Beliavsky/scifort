program compare_driver
    use scifort_signal
    implicit none
    integer :: n, dist, window, ios
    real(dp) :: hm, tm, pm, wm, rh
    real(dp), allocatable :: x(:)
    type(peak_result) :: r
    do
        read(*,*,iostat=ios) n, dist, window, hm, tm, pm, wm, rh
        if (ios /= 0) exit
        allocate(x(n))
        read(*,*) x
        if (window <= 0) then
            r = find_peaks(x, height_min=hm, threshold_min=tm, distance=dist, &
                prominence_min=pm, width_min=wm, rel_height=rh, plateau_size_min=1)
        else
            r = find_peaks(x, height_min=hm, threshold_min=tm, distance=dist, &
                prominence_min=pm, width_min=wm, rel_height=rh, plateau_size_min=1, &
                wlen=window)
        end if
        print '(i0,1x,i0)', r%status, size(r%peaks)
        print '(*(i0,1x))', r%peaks
        print '(*(es23.15e3,1x))', r%prominences
        print '(*(i0,1x))', r%left_bases
        print '(*(i0,1x))', r%right_bases
        print '(*(es23.15e3,1x))', r%widths
        print '(*(es23.15e3,1x))', r%width_heights
        print '(*(es23.15e3,1x))', r%left_ips
        print '(*(es23.15e3,1x))', r%right_ips
        print '(*(i0,1x))', r%left_edges
        print '(*(i0,1x))', r%right_edges
        print '(*(es23.15e3,1x))', r%left_thresholds
        print '(*(es23.15e3,1x))', r%right_thresholds
        deallocate(x)
    end do
end program
