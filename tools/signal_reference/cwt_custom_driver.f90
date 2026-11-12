! SPDX-License-Identifier: MIT
! Fortran half of the SciPy custom-wavelet reference comparison.
module cwt_custom_reference_wavelets
    use scifort_kinds, only : dp
    implicit none
    private
    public :: gaussian_wavelet, shifted_wavelet, modulated_wavelet, derivative_wavelet
contains
    function gaussian_wavelet(points, width) result(w)
        real(dp), intent(in) :: points !! requested sample count
        real(dp), intent(in) :: width !! positive wavelet scale
        real(dp), allocatable :: w(:) !! Gaussian kernel
        real(dp) :: u
        integer :: k
        allocate(w(ceiling(points)))
        do k = 1,size(w)
            u = (real(k-1,dp) - (points-1.0_dp)/2.0_dp)/width
            w(k) = exp(-0.5_dp*u*u)/(sqrt(2.0_dp*acos(-1.0_dp))*width)
        end do
    end function gaussian_wavelet

    function shifted_wavelet(points, width) result(w)
        real(dp), intent(in) :: points !! requested sample count
        real(dp), intent(in) :: width !! positive wavelet scale
        real(dp), allocatable :: w(:) !! shifted Gaussian kernel
        real(dp) :: u
        integer :: k
        allocate(w(ceiling(points)))
        do k = 1,size(w)
            u = (real(k-1,dp) - (points-1.0_dp)/2.0_dp)/width - 0.35_dp
            w(k) = exp(-0.5_dp*u*u)/(sqrt(2.0_dp*acos(-1.0_dp))*width)
        end do
    end function shifted_wavelet

    function modulated_wavelet(points, width) result(w)
        real(dp), intent(in) :: points !! requested sample count
        real(dp), intent(in) :: width !! positive wavelet scale
        real(dp), allocatable :: w(:) !! real oscillating Gaussian kernel
        real(dp) :: u
        integer :: k
        allocate(w(ceiling(points)))
        do k = 1,size(w)
            u = (real(k-1,dp) - (points-1.0_dp)/2.0_dp)/width
            w(k) = cos(3.3_dp*u)*exp(-0.5_dp*u*u)/sqrt(width)
        end do
    end function modulated_wavelet

    function derivative_wavelet(points, width) result(w)
        real(dp), intent(in) :: points !! requested sample count
        real(dp), intent(in) :: width !! positive wavelet scale
        real(dp), allocatable :: w(:) !! derivative Gaussian kernel
        real(dp) :: u
        integer :: k
        allocate(w(ceiling(points)))
        do k = 1,size(w)
            u = (real(k-1,dp) - (points-1.0_dp)/2.0_dp)/width
            w(k) = u*exp(-0.5_dp*u*u)/sqrt(width)
        end do
    end function derivative_wavelet
end module cwt_custom_reference_wavelets

program cwt_custom_driver
    use scifort_kinds, only : dp
    use scifort_signal_cwt_peaks, only : find_peaks_cwt, cwt_wavelet, cwt_peak_success
    use cwt_custom_reference_wavelets
    implicit none
    real(dp), allocatable :: x(:), widths(:), cwt(:,:)
    integer, allocatable :: peaks(:)
    integer :: n, nw, wave_type, minlen, win, i, status, transform_status, io
    real(dp) :: gap, snr, perc
    do
        read(*,*,iostat=io) n, nw, wave_type, minlen, win
        if (io /= 0) exit
        allocate(x(n),widths(nw))
        read(*,*) x
        read(*,*) widths
        read(*,*) gap,snr,perc
        select case (wave_type)
        case (1)
            cwt = cwt_wavelet(x,widths,gaussian_wavelet,transform_status)
            peaks = find_peaks_cwt(x,widths,gap_thresh=gap,min_length=minlen, &
                min_snr=snr,noise_perc=perc,window_size=win,status=status,wavelet=gaussian_wavelet)
        case (2)
            cwt = cwt_wavelet(x,widths,shifted_wavelet,transform_status)
            peaks = find_peaks_cwt(x,widths,gap_thresh=gap,min_length=minlen, &
                min_snr=snr,noise_perc=perc,window_size=win,status=status,wavelet=shifted_wavelet)
        case (3)
            cwt = cwt_wavelet(x,widths,modulated_wavelet,transform_status)
            peaks = find_peaks_cwt(x,widths,gap_thresh=gap,min_length=minlen, &
                min_snr=snr,noise_perc=perc,window_size=win,status=status,wavelet=modulated_wavelet)
        case (4)
            cwt = cwt_wavelet(x,widths,derivative_wavelet,transform_status)
            peaks = find_peaks_cwt(x,widths,gap_thresh=gap,min_length=minlen, &
                min_snr=snr,noise_perc=perc,window_size=win,status=status,wavelet=derivative_wavelet)
        case default
            error stop 'unsupported reference wavelet'
        end select
        write(*,'(i0,1x,i0,1x,i0)') status,size(peaks),transform_status
        write(*,'(*(1x,i0))') peaks
        if (status == cwt_peak_success .and. transform_status == cwt_peak_success) then
            do i=1,nw
                write(*,'(*(1x,es25.16e3))') cwt(i,:)
            end do
        end if
        deallocate(x,widths)
    end do
end program cwt_custom_driver
