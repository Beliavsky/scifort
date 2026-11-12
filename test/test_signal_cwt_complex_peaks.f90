! SPDX-License-Identifier: MIT
! Complex Morlet / complex user wavelet peak detection, using CWT modulus.
program test_signal_cwt_complex_peaks
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_cwt_peaks, only : find_peaks_cwt_complex, &
        find_peaks_cwt_magnitude, find_peaks_cwt, cwt_morlet2, &
        cwt_complex_wavelet, morlet2_wavelet, cwt_peak_success, &
        cwt_peak_invalid_input
    implicit none
    integer, parameter :: n = 101
    complex(dp) :: z(n), constant(4,64), artificial(5,64)
    complex(dp), allocatable :: coefficients(:,:), bad(:,:)
    real(dp) :: widths(5), x(n), t, amplitude, phase, nanval
    integer, allocatable :: peaks(:), expected(:)
    integer :: j, row, status

    nanval = ieee_value(0.0_dp, ieee_quiet_nan)
    widths = [2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
    do j = 1, n
        t = real(j,dp)
        x(j) = 3.0_dp*exp(-0.5_dp*((t-26.0_dp)/5.0_dp)**2) + &
               2.4_dp*exp(-0.5_dp*((t-77.0_dp)/7.0_dp)**2)
        z(j) = x(j)*cmplx(cos(0.88_dp*t), sin(0.88_dp*t), dp)
    end do

    ! A synthetic coefficient matrix has unambiguous modulus peaks despite
    ! strongly varying phases and negative real-valued coefficients.
    do row = 1, size(artificial,1)
        do j = 1, size(artificial,2)
            t = real(j,dp)
            amplitude = 0.05_dp + 3.0_dp*exp(-0.5_dp*((t-18.0_dp)/3.0_dp)**2) + &
                        2.0_dp*exp(-0.5_dp*((t-43.0_dp)/4.0_dp)**2)
            phase = 0.29_dp*t + real(row,dp)*0.8_dp
            artificial(row,j) = amplitude*cmplx(cos(phase),sin(phase),dp)
        end do
    end do
    peaks = find_peaks_cwt_magnitude(artificial,widths,min_length=4, &
        min_snr=1.1_dp,window_size=7,status=status)
    call require(status == cwt_peak_success, 'synthetic coefficient status')
    call require(size(peaks) == 2, 'synthetic coefficient peak count')
    call require(all(peaks == [18,43]), 'synthetic coefficient peak positions')
    coefficients = artificial*cmplx(cos(0.75_dp),sin(0.75_dp),dp)
    peaks = find_peaks_cwt_magnitude(coefficients,widths,min_length=4, &
        min_snr=1.1_dp,window_size=7,status=status)
    call require(status == cwt_peak_success .and. size(peaks)==2, 'phase-rotated status/count')
    call require(all(peaks == [18,43]), 'invariance to complex phase rotation')

    ! Noncontiguous, sliced CWT inputs must match an explicit contiguous copy.
    coefficients = artificial(:,1:64:2)
    expected = find_peaks_cwt_magnitude(coefficients,widths,min_length=3,status=status)
    peaks = find_peaks_cwt_magnitude(artificial(:,1:64:2),widths,min_length=3,status=status)
    call require(status == cwt_peak_success, 'sliced coefficient status')
    call require(size(peaks) == size(expected), 'sliced coefficient count')
    call require(all(peaks == expected), 'sliced coefficient positions')

    coefficients = cwt_morlet2(z,widths,omega0=4.5_dp,status=status)
    call require(status == cwt_peak_success, 'Morlet transform status')
    expected = find_peaks_cwt_magnitude(coefficients,widths,min_length=3, &
        min_snr=0.8_dp,status=status)
    call require(status == cwt_peak_success, 'computed magnitude status')
    peaks = find_peaks_cwt_complex(z,widths,omega0=4.5_dp, &
        min_length=3,min_snr=0.8_dp,status=status)
    call require(status == cwt_peak_success,'complex input peak status')
    call require(all(peaks == expected) .and. size(peaks) == size(expected), &
        'Morlet and precomputed magnitude equivalence')
    call require(any(abs(peaks - 26) <= 2) .and. any(abs(peaks - 77) <= 2), &
        'Morlet complex peak localization')

    ! Explicit custom wavelets share their modulus-ridge path with the
    ! precomputed-coefficient interface, and tolerate strided data.
    coefficients = cwt_complex_wavelet(z(1:n:2),widths,custom_wavelet,status)
    expected = find_peaks_cwt_magnitude(coefficients,widths,min_length=2,status=status)
    peaks = find_peaks_cwt_complex(z(1:n:2),widths,wavelet=custom_wavelet, &
        min_length=2,status=status)
    call require(status == cwt_peak_success, 'strided custom callback status')
    call require(size(peaks) == size(expected), 'custom peak count')
    call require(all(peaks == expected), 'custom peak positions')

    coefficients = cwt_morlet2(x,widths,omega0=2.8_dp)
    expected = find_peaks_cwt_magnitude(coefficients,widths,status=status)
    peaks = find_peaks_cwt_complex(x,widths,omega0=2.8_dp,status=status)
    call require(status == cwt_peak_success, 'real input peak status')
    call require(size(peaks) == size(expected), 'real input peak count')
    call require(all(peaks == expected), 'real input peak positions')

    ! Constant modulus: no strict interior extrema even for rotating phase.
    ! Exact constant complex coefficients provide a deterministic no-peak
    ! case without modulus rounding differences from unit phasors.
    constant = cmplx(-2.0_dp,3.0_dp,dp)
    peaks = find_peaks_cwt_magnitude(constant,[1.0_dp,2.0_dp,3.0_dp,4.0_dp],status=status)
    call require(status == cwt_peak_success .and. size(peaks)==0, 'flat modulus')

    peaks = find_peaks_cwt_complex(z,[-1.0_dp],status=status)
    call invalid(peaks,status,'negative scale')
    peaks = find_peaks_cwt_complex(z,widths,omega0=nanval,status=status)
    call invalid(peaks,status,'nonfinite frequency')
    peaks = find_peaks_cwt_complex(z,widths,wavelet=bad_length_wavelet,status=status)
    call invalid(peaks,status,'wrong complex wavelet length')
    peaks = find_peaks_cwt_complex(z,widths,wavelet=nan_wavelet,status=status)
    call invalid(peaks,status,'nonfinite complex wavelet')
    peaks = find_peaks_cwt_complex(z,widths,wavelet=custom_wavelet, &
        omega0=4.0_dp,status=status)
    call invalid(peaks,status,'conflicting Morlet and callback')
    peaks = find_peaks_cwt_complex(z,widths,noise_perc=101.0_dp,status=status)
    call invalid(peaks,status,'bad noise percentile')
    peaks = find_peaks_cwt_complex(z,widths,min_length=0,status=status)
    call invalid(peaks,status,'zero ridge length')
    peaks = find_peaks_cwt_complex(z,widths,max_distances=[1.0_dp],status=status)
    call invalid(peaks,status,'short max_distances')
    peaks = find_peaks_cwt_complex(z,widths,min_snr=nanval,status=status)
    call invalid(peaks,status,'nonfinite SNR')
    peaks = find_peaks_cwt_complex(z,widths,window_size=0,status=status)
    call invalid(peaks,status,'invalid noise window')
    peaks = find_peaks_cwt_magnitude(artificial,widths(:4),status=status)
    call invalid(peaks,status,'mismatched coefficient rows')
    bad = artificial
    bad(1,1) = cmplx(0.0_dp,nanval,dp)
    peaks = find_peaks_cwt_magnitude(bad,widths,status=status)
    call invalid(peaks,status,'nonfinite precomputed coefficient')
    peaks = find_peaks_cwt_magnitude(artificial(:,1:0),widths,status=status)
    call invalid(peaks,status,'empty coefficient matrix')
    peaks = find_peaks_cwt_magnitude(artificial,widths,gap_thresh=-1.0_dp,status=status)
    call invalid(peaks,status,'negative gap threshold')

    ! The old real-valued Ricker API is unchanged after ridge-code reuse.
    peaks = find_peaks_cwt(x,widths,status=status)
    call require(status==cwt_peak_success,'legacy real-wavelet peak status')
    call require(any(abs(peaks-26)<=2) .and. any(abs(peaks-77)<=2), &
        'legacy real-wavelet peak locations')

    print '(a)', 'PASS complex magnitude-based CWT peak tests'
contains
    subroutine require(ok,message)
        logical, intent(in) :: ok !! true when assertion holds
        character(*), intent(in) :: message !! failure diagnostic
        if (.not. ok) then
            print '(a)', 'FAIL: '//message
            error stop 1
        end if
    end subroutine require

    subroutine invalid(p,code,message)
        integer, intent(in) :: p(:) !! expected empty result
        integer, intent(in) :: code !! returned status
        character(*), intent(in) :: message !! invalid-case diagnostic
        call require(code==cwt_peak_invalid_input .and. size(p)==0,message)
    end subroutine invalid

    function custom_wavelet(points,width) result(wavelet)
        real(dp), intent(in) :: points !! unrounded wavelet length
        real(dp), intent(in) :: width !! positive scale
        complex(dp), allocatable :: wavelet(:) !! asymmetric complex samples
        integer :: k
        real(dp) :: r
        wavelet = morlet2_wavelet(points,width)
        do k=1,size(wavelet)
            r=(real(k-1,dp)-0.5_dp*(points-1.0_dp))/width
            wavelet(k) = wavelet(k)*cmplx(1.0_dp+0.08_dp*r,0.2_dp,dp)
        end do
    end function custom_wavelet

    function bad_length_wavelet(points,width) result(wavelet)
        real(dp), intent(in) :: points !! desired wavelet length
        real(dp), intent(in) :: width !! scale (unused)
        complex(dp), allocatable :: wavelet(:) !! incorrectly sized samples
        allocate(wavelet(max(0,ceiling(points)-1)))
        wavelet = cmplx(width,0.0_dp,dp)
    end function bad_length_wavelet

    function nan_wavelet(points,width) result(wavelet)
        real(dp), intent(in) :: points !! unrounded wavelet length
        real(dp), intent(in) :: width !! wavelet scale
        complex(dp), allocatable :: wavelet(:) !! nonfinite samples
        wavelet = custom_wavelet(points,width)
        wavelet(1) = cmplx(nanval,0.0_dp,dp)
    end function nan_wavelet
end program test_signal_cwt_complex_peaks
