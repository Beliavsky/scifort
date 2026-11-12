! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Independent SciPy reference driver for per-segment detrending callbacks.
program signal_detrend_callback_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : periodogram, welch, csd, coherence, stft, &
        spectrogram, spectral_result, cross_spectral_result, stft_result, spectrogram_result, &
        spectral_status_success
    use scifort_signal_spectral_complex, only : periodogram_complex, &
        welch_complex, csd_complex, coherence_complex, stft_complex, spectrogram_complex
    implicit none
    interface emit
        procedure emit_vector
        procedure emit_matrix
    end interface
    integer, parameter :: n = 64
    integer :: j, config, ns, no, nf, status, cnt
    real(dp) :: x(n), y(n), t
    complex(dp) :: z(n), w(n)
    logical :: onesided
    type(spectral_result) :: p
    type(cross_spectral_result) :: c
    type(stft_result) :: f
    type(spectrogram_result) :: s
    character(len=12) :: label

    do j = 1, n
        t = real(j, dp)
        x(j) = sin(0.19_dp*t) + 0.25_dp*cos(0.72_dp*t) + &
            0.003_dp*t*t + 0.05_dp*t
        y(j) = cos(0.25_dp*t) + 0.14_dp*sin(0.34_dp*t) - &
            0.0004_dp*t*t + 0.025_dp*t
    end do
    z = cmplx(x, 0.45_dp*y, dp)
    w = cmplx(y, -0.2_dp*x + 0.1_dp*y, dp)
    cnt = 0
    do config = 1, 4
        select case(config)
        case (1)
            ns = 12
            no = 6
        case (2)
            ns = 16
            no = 8
        case (3)
            ns = 23
            no = 9
        case (4)
            ns = 32
            no = 13
        end select
        nf = ns + 3
        onesided = mod(config, 2) == 1
        ! An 8th-order periodogram FFT also tests input truncation.
        p = periodogram(x, fs=9.0_dp, window='hann', nfft=nf, &
            detrend_callback=quadratic_real)
        call check(p%status)
        cnt = cnt + 1
        call emit(cnt, p%frequencies, cmplx(p%power, 0.0_dp, dp))
        p = welch(x, fs=9.0_dp, window='hann', nperseg=ns, &
            noverlap=no, nfft=nf, average='median', return_onesided=onesided, &
            detrend_callback=quadratic_real)
        call check(p%status)
        cnt = cnt + 1
        call emit(cnt, p%frequencies, cmplx(p%power, 0.0_dp, dp))
        c = csd(x, y, fs=9.0_dp, window='hann', nperseg=ns, &
            noverlap=no, nfft=nf, return_onesided=onesided, &
            detrend_callback=quadratic_real)
        call check(c%status)
        cnt = cnt + 1
        call emit(cnt, c%frequencies, c%power)
        p = coherence(x, y, fs=9.0_dp, window='hann', nperseg=ns, &
            noverlap=no, nfft=nf, detrend_callback=quadratic_real)
        call check(p%status)
        cnt = cnt + 1
        call emit(cnt, p%frequencies, cmplx(p%power, 0.0_dp, dp))
        f = stft(x, fs=9.0_dp, window='hann', nperseg=ns, &
            noverlap=no, nfft=nf, boundary='none', padded=.false., &
            return_onesided=onesided, detrend_callback=quadratic_real)
        call check(f%status)
        cnt = cnt + 1
        call emit(cnt, f%frequencies, f%values, f%times)
        s = spectrogram(x, fs=9.0_dp, window='hann', nperseg=ns, &
            noverlap=no, nfft=nf, mode='complex', return_onesided=onesided, &
            detrend_callback=quadratic_real)
        call check(s%status)
        cnt = cnt + 1
        call emit(cnt, s%frequencies, s%values, s%times)
        p = periodogram_complex(z, fs=9.0_dp, window='hann', nfft=nf, &
            detrend_callback=quadratic_complex)
        call check(p%status)
        cnt = cnt + 1
        call emit(cnt, p%frequencies, cmplx(p%power, 0.0_dp, dp))
        p = welch_complex(z, fs=9.0_dp, window='hann', nperseg=ns, &
            noverlap=no, nfft=nf, average='median', detrend_callback=quadratic_complex)
        call check(p%status)
        cnt = cnt + 1
        call emit(cnt, p%frequencies, cmplx(p%power, 0.0_dp, dp))
        c = csd_complex(z, w, fs=9.0_dp, window='hann', nperseg=ns, &
            noverlap=no, nfft=nf, detrend_callback=quadratic_complex)
        call check(c%status)
        cnt = cnt + 1
        call emit(cnt, c%frequencies, c%power)
        p = coherence_complex(z, w, fs=9.0_dp, window='hann', nperseg=ns, &
            noverlap=no, nfft=nf, detrend_callback=quadratic_complex)
        call check(p%status)
        cnt = cnt + 1
        call emit(cnt, p%frequencies, cmplx(p%power, 0.0_dp, dp))
        f = stft_complex(z, fs=9.0_dp, window='hann', nperseg=ns, &
            noverlap=no, nfft=nf, boundary='none', padded=.false., &
            detrend_callback=quadratic_complex)
        call check(f%status)
        cnt = cnt + 1
        call emit(cnt, f%frequencies, f%values, f%times)
        s = spectrogram_complex(z, fs=9.0_dp, window='hann', nperseg=ns, &
            noverlap=no, nfft=nf, mode='complex', detrend_callback=quadratic_complex)
        call check(s%status)
        cnt = cnt + 1
        call emit(cnt, s%frequencies, s%values, s%times)
    end do
contains
    subroutine check(code)
        integer, intent(in) :: code !! spectral result status
        if (code /= spectral_status_success) error stop 'reference spectral status failed'
    end subroutine check

    subroutine quadratic_real(samples)
        real(dp), intent(inout) :: samples(:) !! unwindowed segment
        real(dp) :: coordinate(size(samples)), square(size(samples)), av, slope, curve
        integer :: i, size_line
        size_line = size(samples)
        do i = 1, size_line
            coordinate(i) = (real(i,dp)-0.5_dp*real(size_line+1,dp))/real(size_line,dp)
        end do
        square = coordinate*coordinate
        square = square - sum(square)/real(size_line,dp)
        av = sum(samples)/real(size_line,dp)
        slope = sum(samples*coordinate)/sum(coordinate*coordinate)
        curve = sum(samples*square)/sum(square*square)
        samples = samples - av - slope*coordinate - curve*square
    end subroutine quadratic_real

    subroutine quadratic_complex(samples)
        complex(dp), intent(inout) :: samples(:) !! complex unwindowed segment
        real(dp) :: coordinate(size(samples)), square(size(samples))
        complex(dp) :: av, slope, curve
        integer :: i, size_line
        size_line = size(samples)
        do i = 1, size_line
            coordinate(i) = (real(i,dp)-0.5_dp*real(size_line+1,dp))/real(size_line,dp)
        end do
        square = coordinate*coordinate
        square = square - sum(square)/real(size_line,dp)
        av = sum(samples)/real(size_line,dp)
        slope = sum(samples*coordinate)/sum(coordinate*coordinate)
        curve = sum(samples*square)/sum(square*square)
        samples = samples - av - slope*coordinate - curve*square
    end subroutine quadratic_complex

    subroutine emit_vector(num, frequency, spectrum, times)
        integer, intent(in) :: num !! comparison case number
        real(dp), intent(in) :: frequency(:) !! spectral frequency grid
        complex(dp), intent(in) :: spectrum(:) !! spectral estimates
        real(dp), intent(in), optional :: times(:) !! optional time grid
        complex(dp) :: mat(size(spectrum),1)
        mat(:,1) = spectrum
        call emit_matrix(num,frequency,mat,times)
    end subroutine emit_vector

    subroutine emit_matrix(num, frequency, spectrum, times)
        integer, intent(in) :: num !! case number
        real(dp), intent(in) :: frequency(:) !! spectral frequency grid
        complex(dp), intent(in) :: spectrum(:,:) !! spectral estimates
        real(dp), intent(in), optional :: times(:) !! optional time grid
        integer :: i, j
        real(dp) :: timestamp
        write(*,'(a,1x,i0,1x,i0,1x,i0)') 'CASE', num, size(spectrum,1), size(spectrum,2)
        do j = 1,size(spectrum,2)
            timestamp = 0.0_dp
            if (present(times)) timestamp = times(j)
            do i = 1,size(spectrum,1)
                write(*,'(4(es25.16,1x))') frequency(i), timestamp, &
                    real(spectrum(i,j),dp), aimag(spectrum(i,j))
            end do
        end do
    end subroutine emit_matrix
end program signal_detrend_callback_reference_driver
