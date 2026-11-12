! SPDX-License-Identifier: MIT
! Stream real/complex Morlet/custom ridge-detection cases to a Python checker.
program cwt_complex_peaks_driver
    use scifort_kinds, only : dp
    use scifort_signal_cwt_peaks, only : find_peaks_cwt_complex
    implicit none
    complex(dp), allocatable :: x(:)
    real(dp), allocatable :: widths(:), distances(:)
    integer, allocatable :: peaks(:)
    integer :: n, nw, mode, minlen, win, use_dist, status, io, j
    real(dp) :: gap, snr, perc, omega0

    do
        read(*,*,iostat=io) n, nw, mode, minlen, win, use_dist
        if (io /= 0) exit
        allocate(x(n),widths(nw),distances(nw))
        read(*,*) x
        read(*,*) widths
        read(*,*) gap, snr, perc, omega0
        if (use_dist /= 0) then
            read(*,*) distances
        else
            distances = widths/4.0_dp
        end if
        select case(mode)
        case(0)
            peaks = find_peaks_cwt_complex(real(x,dp),widths, &
                max_distances=distances,gap_thresh=gap,min_length=minlen, &
                min_snr=snr,noise_perc=perc,window_size=win,omega0=omega0,status=status)
        case(1)
            peaks = find_peaks_cwt_complex(x,widths, &
                max_distances=distances,gap_thresh=gap,min_length=minlen, &
                min_snr=snr,noise_perc=perc,window_size=win,omega0=omega0,status=status)
        case(2)
            peaks = find_peaks_cwt_complex(real(x,dp),widths, &
                max_distances=distances,gap_thresh=gap,min_length=minlen, &
                min_snr=snr,noise_perc=perc,window_size=win,wavelet=custom_wavelet,status=status)
        case(3)
            peaks = find_peaks_cwt_complex(x,widths, &
                max_distances=distances,gap_thresh=gap,min_length=minlen, &
                min_snr=snr,noise_perc=perc,window_size=win,wavelet=custom_wavelet,status=status)
        case default
            error stop 'unrecognized complex ridge test mode'
        end select
        write(*,'(i0,1x,i0)',advance='no') status,size(peaks)
        do j=1,size(peaks)
            write(*,'(1x,i0)',advance='no') peaks(j)
        end do
        write(*,*)
        deallocate(x,widths,distances)
    end do
contains
    function custom_wavelet(points,width) result(wavelet)
        use scifort_signal_cwt_peaks, only : morlet2_wavelet
        real(dp), intent(in) :: points !! unrounded requested length
        real(dp), intent(in) :: width !! positive wavelet scale
        complex(dp), allocatable :: wavelet(:) !! complex asymmetric wavelet
        real(dp) :: r
        integer :: k

        wavelet = morlet2_wavelet(points,width)
        do k=1,size(wavelet)
            r=(real(k-1,dp)-0.5_dp*(points-1.0_dp))/width
            wavelet(k) = wavelet(k)*cmplx(1.0_dp+0.08_dp*r,0.2_dp,dp)
        end do
    end function custom_wavelet
end program cwt_complex_peaks_driver
