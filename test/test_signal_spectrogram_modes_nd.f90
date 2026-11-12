! SPDX-License-Identifier: MIT
! Verify that rank-independent spectrogram modes use the 1-D spectral kernel.
program test_signal_spectrogram_modes_nd
    use scifort_kinds, only : dp
    use scifort_signal_timefreq_nd, only : spectrogram_into
    use scifort_signal_spectral, only : spectrogram_real => spectrogram, &
        spectrogram_result, spectral_status_success, spectral_status_invalid_input
    implicit none

    real(dp) :: x(33,4), real_output(7,2,5), failure(7,2,5)
    complex(dp) :: complex_output(12,2,5), z(33,2)
    real(dp), allocatable :: f(:), t(:)
    type(spectrogram_result) :: one
    integer :: i,j,status

    do j=1,4
        do i=1,33
            x(i,j) = cos(0.71_dp*real(i*j,dp)) - 0.3_dp*sin(0.15_dp*real(i,dp))
        end do
    end do
    real_output=-999.0_dp
    ! Noncontiguous inputs and outputs; negative axis selects the time axis.
    call spectrogram_into(x(:,1:4:2),real_output,f,t,nperseg=9,noverlap=3,nfft=12, &
        mode='phase',axis=-2,status=status,return_onesided=.true.)
    if (status /= spectral_status_success) error stop 'rank-independent phase status'
    do j=1,2
        one=spectrogram_real(x(:,2*j-1),nperseg=9,noverlap=3,nfft=12, &
            mode='phase',return_onesided=.true.)
        if (maxval(abs(real_output(:,j,:)-one%power))>2.0e-12_dp) &
            error stop 'rank-independent phase mismatch'
    end do

    complex_output=(-999.0_dp,111.0_dp)
    call spectrogram_into(x(:,1:4:2),complex_output,f,t,nperseg=9,noverlap=3,nfft=12, &
        mode='complex',axis=1,status=status,return_onesided=.false.)
    if (status /= spectral_status_success) error stop 'rank-independent complex status'
    do j=1,2
        one=spectrogram_real(x(:,2*j-1),nperseg=9,noverlap=3,nfft=12, &
            mode='complex',return_onesided=.false.)
        if (maxval(abs(complex_output(:,j,:)-one%values))>2.0e-12_dp) &
            error stop 'rank-independent complex mismatch'
    end do

    ! Fortran rank-2 complex inputs, the frequency axis is full and ordered.
    z=cmplx(x(:,1:4:2),0.29_dp*x(:,2:4:2),dp)
    call spectrogram_into(z,complex_output,f,t,nperseg=9,noverlap=3,nfft=12, &
        mode='complex',axis=1,status=status)
    if (status /= spectral_status_success) error stop 'complex source status'
    if (size(f) /= 12) error stop 'complex source frequency count'

    ! Preserve the destination after a bad mode or incorrect output shape.
    failure=-999.0_dp
    call spectrogram_into(x(:,1:4:2),failure,f,t,nperseg=9,noverlap=3,nfft=12, &
        mode='invalid',axis=1,status=status)
    if (status /= spectral_status_invalid_input) error stop 'invalid mode status'
    if (any(failure/=-999.0_dp)) error stop 'invalid mode corrupted destination'
    call spectrogram_into(x(:,1:4:2),failure,f,t,nperseg=9,noverlap=3,nfft=12, &
        mode='complex',axis=1,status=status)
    if (status /= spectral_status_invalid_input) error stop 'wrong result type status'
    if (any(failure/=-999.0_dp)) error stop 'wrong result type corrupted destination'
    print *, 'PASS: rank-independent spectrogram modes'
end program test_signal_spectrogram_modes_nd
