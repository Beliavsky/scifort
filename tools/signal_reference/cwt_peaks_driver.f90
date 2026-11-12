program cwt_driver
    use scifort_kinds, only : dp
    use scifort_signal_cwt_peaks, only : find_peaks_cwt, cwt_ricker
    implicit none
    real(dp), allocatable :: x(:), widths(:), distances(:)
    integer, allocatable :: p(:)
    integer :: n, m, i, status, win, minlen, use_dist
    real(dp) :: gap, snr, perc
    do
        read(*,*,iostat=status) n, m, minlen, win, use_dist
        if (status /= 0) exit
        allocate(x(n),widths(m))
        read(*,*) x
        read(*,*) widths
        read(*,*) gap, snr, perc
        if (use_dist == 1) then
            allocate(distances(m))
            read(*,*) distances
            p = find_peaks_cwt(x,widths,max_distances=distances, &
                gap_thresh=gap,min_length=minlen,min_snr=snr,noise_perc=perc,window_size=win,status=status)
        else
            p = find_peaks_cwt(x,widths,gap_thresh=gap, &
                min_length=minlen,min_snr=snr,noise_perc=perc,window_size=win,status=status)
        end if
        write(*,'(i0,1x,i0)',advance='no') status,size(p)
        do i = 1,size(p)
            write(*,'(1x,i0)',advance='no') p(i)
        end do
        write(*,*)
        deallocate(x,widths)
        if (allocated(distances)) deallocate(distances)
    end do
end program cwt_driver
