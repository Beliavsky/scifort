program test_signal_cwt_peaks
    use scifort_kinds, only : dp
    use scifort_signal, only : find_peaks_cwt, cwt_ricker, ricker_wavelet, &
        cwt_peak_success, cwt_peak_invalid_input
    implicit none
    real(dp) :: x(101), widths(7), t, w(1)
    real(dp), allocatable :: cwt(:,:)
    integer, allocatable :: peaks(:)
    integer :: i, status

    widths = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp]
    do i=1,101
        t = real(i,dp)
        x(i) = 4.0_dp * exp(-0.5_dp*((t-29.0_dp)/3.0_dp)**2) + &
               3.0_dp * exp(-0.5_dp*((t-76.0_dp)/5.0_dp)**2)
    end do
    peaks = find_peaks_cwt(x,widths,status=status)
    if (status /= cwt_peak_success) error stop 'cwt status'
    ! Symmetric nearly tied maxima can move by one bin under -O3 roundoff.
    if (.not. any(abs(peaks - 29) <= 1) .or. &
        .not. any(abs(peaks - 76) <= 1)) error stop 'cwt peaks'
    if (any(peaks < 1) .or. any(peaks > size(x))) error stop 'peak range'
    cwt = cwt_ricker(x,widths,status)
    if (status /= cwt_peak_success) error stop 'cwt transform'
    if (size(cwt,1) /= 7 .or. size(cwt,2) /= 101) error stop 'cwt shape'
    if (cwt(4,29) < 1.0_dp) error stop 'wavelet response'
    w = ricker_wavelet(1.0_dp, 1.0_dp)
    if (abs(w(1) - 2.0_dp/(sqrt(3.0_dp)*acos(-1.0_dp)**0.25_dp)) > 1.0e-12_dp) &
        error stop 'ricker normalisation'
    peaks = find_peaks_cwt(x, [-1.0_dp], status=status)
    if (status /= cwt_peak_invalid_input .or. size(peaks) /= 0) error stop 'invalid width'
    peaks = find_peaks_cwt(x, widths, noise_perc=101.0_dp, status=status)
    if (status /= cwt_peak_invalid_input) error stop 'invalid percentile'
    peaks = find_peaks_cwt(x, widths, max_distances=[0.1_dp], status=status)
    if (status /= cwt_peak_invalid_input) error stop 'invalid distance shape'
    print *, 'signal cwt peaks passed'
end program test_signal_cwt_peaks
