! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Confirm that FFT plans are reused across convolution lines, multitaper
! channels, and stream resets/reinitialization, including prime FFT lengths.
program test_signal_fft_plan_reuse_consumers
    use scifort_kinds, only : dp
    use scifort_signal, only : fftconvolve_nd, oaconvolve_nd, &
        multitaper_stream, multitaper_spectral_matrix_spectrogram, &
        multitaper_spectral_matrix_spectrogram_result, &
        spectral_status_success, spectral_status_invalid_input
    implicit none
    integer, parameter :: ns=225, nc=2, nseg=80
    real(dp) :: samples(ns,nc)
    complex(dp) :: complex_samples(ns,nc)
    real(dp) :: a(129,7), b(5,3), y1(133,9), y2(133,9)
    real(dp) :: a15(129,1,1,1,1,1,1,1,1,1,1,1,1,1,7)
    real(dp) :: b15(5,1,1,1,1,1,1,1,1,1,1,1,1,1,3)
    real(dp) :: y15(133,1,1,1,1,1,1,1,1,1,1,1,1,1,9)
    type(multitaper_stream) :: stream
    type(multitaper_spectral_matrix_spectrogram_result) :: batch, part1, part2, full
    complex(dp), allocatable :: observed(:,:,:,:)
    integer :: i, ch, code, nfft, j, nfirst
    integer, parameter :: lengths(2) = [263,315]
    real(dp), parameter :: tol=6.0e-11_dp

    do ch=1,nc
        do i=1,ns
            samples(i,ch) = sin(0.057_dp*real(i*ch,dp))+0.3_dp*cos(0.17_dp*real(i+ch,dp))
            complex_samples(i,ch) = cmplx(samples(i,ch),0.2_dp*sin(0.093_dp*real(i+ch,dp)),dp)
        end do
    end do
    do j=1,size(a,2)
        do i=1,size(a,1)
            a(i,j) = sin(0.027_dp*real(i*j,dp))
        end do
    end do
    do j=1,size(b,2)
        do i=1,size(b,1)
            b(i,j) = cos(0.13_dp*real(i+j,dp))
        end do
    end do
    call fftconvolve_nd(a,b,y1,axes=[1,2],status=code)
    if (code /= 0) error stop 'planned 2d convolution failed'
    call oaconvolve_nd(a,b,y2,axes=[1,2],status=code)
    if (code /= 0 .or. maxval(abs(y1-y2)) > 3.0e-11_dp) &
        error stop 'planned overlap-add mismatch'
    a15 = reshape(a,shape(a15))
    b15 = reshape(b,shape(b15))
    call fftconvolve_nd(a15,b15,y15,axes=[1,15],status=code)
    if (code /= 0) error stop 'planned rank15 convolution failed'
    if (maxval(abs(reshape(y15,[size(y1)])-reshape(y1,[size(y1)]))) > 3.0e-11_dp) &
        error stop 'planned rank15 values differ'

    do j=1,size(lengths)
        nfft = lengths(j)
        batch = multitaper_spectral_matrix_spectrogram(samples,3.5_dp,nseg, &
            noverlap=45,kmax=5,nfft=nfft,weighting='adaptive')
        if (batch%status /= spectral_status_success) error stop 'batch multitaper failed'
        call stream%initialize(3.5_dp,nseg,nc,code,noverlap=45,kmax=5,nfft=nfft,weighting='adaptive')
        if (code /= spectral_status_success) error stop 'stream init failed'
        part1 = stream%push(samples(1:137,:))
        part2 = stream%push(samples(138:ns,:))
        if (part1%status /= spectral_status_success .or. part2%status /= spectral_status_success) &
            error stop 'stream chunk failed'
        nfirst = size(part1%times)
        if (nfirst+size(part2%times) /= size(batch%times)) error stop 'stream frame count'
        allocate(observed(size(batch%power,1),size(batch%power,2),nc,nc))
        observed(:,1:nfirst,:,:) = part1%power
        observed(:,nfirst+1:,:,:) = part2%power
        if (maxval(abs(observed-batch%power)) > tol) error stop 'stream versus batch mismatch'
        deallocate(observed)
        call stream%reset()
        full = stream%push(samples)
        if (full%status /= spectral_status_success) error stop 'stream after reset failed'
        if (maxval(abs(full%power-batch%power)) > tol) error stop 'reset changed FFT plan'
    end do

    ! Reinitialize the same instance as a complex, two-sided prime-length stream.
    batch = multitaper_spectral_matrix_spectrogram(complex_samples,3.5_dp,nseg, &
        noverlap=40,kmax=5,nfft=263,weighting='adaptive')
    call stream%initialize(3.5_dp,nseg,nc,code,noverlap=40,kmax=5,nfft=263, &
        weighting='adaptive',complex_input=.true.)
    if (code /= spectral_status_success) error stop 'complex stream init'
    full = stream%push(complex_samples)
    if (full%status /= spectral_status_success .or. maxval(abs(full%power-batch%power)) > tol) &
        error stop 'complex plan reuse mismatch'

    ! Reject invalid FFT lengths without retaining an accidentally usable plan.
    call stream%initialize(3.5_dp,nseg,nc,code,nfft=20)
    if (code /= spectral_status_invalid_input) error stop 'invalid reinit accepted'
    full = stream%push(samples)
    if (full%status /= spectral_status_invalid_input) error stop 'invalid stream remained active'
    call stream%initialize(3.5_dp,nseg,nc,code,nfft=263)
    if (code /= spectral_status_success) error stop 'recover after failed init'
    full = stream%push(samples)
    if (full%status /= spectral_status_success) error stop 'recovered stream push failed'
    print '(a)', 'test_signal_fft_plan_reuse_consumers passed'
end program test_signal_fft_plan_reuse_consumers
