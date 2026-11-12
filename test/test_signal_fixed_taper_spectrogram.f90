! SPDX-License-Identifier: MIT
program test_signal_fixed_taper_spectrogram
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : spectral_status_success, spectral_status_invalid_input
    use scifort_signal_multitaper, only : multitaper_spectral_matrix_spectrogram_result
    use scifort_signal_fixed_taper_spectrogram, only : fixed_taper_spectral_matrix_spectrogram
    implicit none
    type(multitaper_spectral_matrix_spectrogram_result) :: result, single, real_result, bad
    real(dp) :: x(24,3), taper(8,3), weights(3), bands(2,3), freqs(3)
    complex(dp) :: z(24,3)
    integer :: i, j, k

    do j = 1, 3
        do i = 1, 24
            x(i,j) = sin(real(i*j,dp)*0.17_dp)+cos(real(i+3*j,dp)*0.13_dp)
            z(i,j) = cmplx(x(i,j),0.1_dp*cos(real(i+j,dp)*0.4_dp),dp)
        end do
        do i = 1, 8
            taper(i,j) = 0.3_dp+real(j,dp)*sin(real(i+j,dp)*0.15_dp)
        end do
    end do
    bands(:,1) = [-0.5_dp,0.5_dp]
    bands(:,2) = [0.02_dp,0.021_dp]
    bands(:,3) = [-0.26_dp,0.13_dp]
    freqs = [-0.17_dp,0.0_dp,0.26_dp]
    weights = [1.0_dp,0.0_dp,3.0_dp]
    result = fixed_taper_spectral_matrix_spectrogram(z,taper,freqs, &
        bands=bands,taper_weights=weights)
    call check(result%status == spectral_status_success,'complex status')
    call check(all(shape(result%power) == [6,3,3,3]),'complex dimensions')
    call check(abs(result%effective_tapers-1.6_dp) < 1.0e-13_dp,'effective K')
    call check(all(result%times == [4.0_dp,12.0_dp,20.0_dp]),'time centers')
    do k = 1, size(result%power,2)
        do j = 1, size(result%power,1)
            call check(maxval(abs(result%power(j,k,:,:)- &
                conjg(transpose(result%power(j,k,:,:))))) < 1.0e-12_dp, &
                'complex Hermitian spectral matrix')
        end do
    end do
    ! An unused zero-energy taper is permitted, and weights are scale-invariant.
    taper(:,2) = 0.0_dp
    single = fixed_taper_spectral_matrix_spectrogram(z,taper,freqs, &
        bands=bands,taper_weights=weights*100000.0_dp)
    call check(single%status == spectral_status_success,'unused zero-energy taper')
    call check(maxval(abs(result%power-single%power)) < 1.0e-11_dp, &
        'relative weights independent of scaling')
    real_result = fixed_taper_spectral_matrix_spectrogram(x,taper,freqs, &
        bands=bands,taper_weights=weights)
    single = fixed_taper_spectral_matrix_spectrogram(cmplx(x,0.0_dp,dp), &
        taper,freqs,bands=bands,taper_weights=weights)
    call check(maxval(abs(real_result%power-single%power)) < 1.0e-12_dp, &
        'real/complex overload agreement')
    call check(all(real_result%frequencies(1:3) == freqs),'point frequencies')
    call check(abs(real_result%frequencies(4)) < 1.0e-12_dp,'band center')
    single = fixed_taper_spectral_matrix_spectrogram(x,taper,freqs, &
        bands=bands,taper_weights=weights,noverlap=4)
    call check(size(single%power,2) == 5,'overlapping segments')
    bad = fixed_taper_spectral_matrix_spectrogram(x,taper,freqs, &
        taper_weights=[0.0_dp,0.0_dp,0.0_dp])
    call check(bad%status == spectral_status_invalid_input,'zero taper weights invalid')
    bad = fixed_taper_spectral_matrix_spectrogram(x,taper,freqs, &
        taper_weights=[1.0_dp,-0.3_dp,0.0_dp])
    call check(bad%status == spectral_status_invalid_input,'negative taper weights invalid')
    bad = fixed_taper_spectral_matrix_spectrogram(x,taper,freqs,noverlap=8)
    call check(bad%status == spectral_status_invalid_input,'invalid overlap')
    bad = fixed_taper_spectral_matrix_spectrogram(x,taper,freqs, &
        bands=reshape([0.2_dp,-0.2_dp],[2,1]))
    call check(bad%status == spectral_status_invalid_input,'reversed band')
    taper(1,1) = ieee_value(0.0_dp,ieee_quiet_nan)
    bad = fixed_taper_spectral_matrix_spectrogram(x,taper,freqs)
    call check(bad%status == spectral_status_invalid_input,'NaN taper')
    print '(a)', 'PASS fixed-taper exact spectral matrices, contracts and overloads'
contains
    subroutine check(ok,message)
        logical, intent(in) :: ok !! condition that must hold
        character(*), intent(in) :: message !! description on failure
        if (ok) return
        print '(a)', 'FAIL: '//message
        error stop 1
    end subroutine check
end program test_signal_fixed_taper_spectrogram
