! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_multitaper_stream
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : spectral_status_success, spectral_status_invalid_input
    use scifort_signal_multitaper, only : multitaper_stream, multitaper_spectral_matrix_spectrogram, &
        multitaper_spectral_matrix_spectrogram_result, multitaper_matrix_result
    implicit none
    integer, parameter :: n=103, nc=3, nseg=24, nk=16
    real(dp) :: data(n,nc), x(n)
    complex(dp) :: complex_data(n,nc), xc(n)
    type(multitaper_stream) :: stream
    type(multitaper_spectral_matrix_spectrogram_result) :: batch, partial
    type(multitaper_matrix_result) :: avg
    complex(dp), allocatable :: collected(:,:,:,:), complex_collected(:,:,:,:)
    real(dp), allocatable :: time_collected(:)
    integer :: i, c, j, nf, first, last, block, cursor, status, nwindows, weight_index
    integer, parameter :: chunks(9) = [7,1,0,18,23,3,15,9,27]
    character(8), parameter :: weightings(3) = [character(8) :: 'equal', 'eigen', 'adaptive']
    real(dp), parameter :: tol=1.0e-10_dp
    real(dp) :: pi
    pi = acos(-1.0_dp)
    do c = 1, nc
        do i = 1, n
            data(i,c) = sin(0.13_dp*real(c*i,dp))+0.4_dp*cos(0.032_dp*real(i*i+c,dp))
            complex_data(i,c) = cmplx(data(i,c),sin(0.23_dp*real(i+c,dp)),dp)
        end do
    end do
    x = data(:,1)
    xc = complex_data(:,1)

    ! Each of three weighting schemes must agree with batch even when chunks
    ! terminate in the middle of windows, include empty calls, or emit many windows.
    do weight_index = 1, 3
        batch = multitaper_spectral_matrix_spectrogram(data,2.5_dp,nseg,noverlap=15, &
            kmax=4,fs=8.0_dp,nfft=31,weighting=trim(weightings(weight_index)))
        if (batch%status /= spectral_status_success) error stop 'real batch failed'
        call stream%initialize(2.5_dp,nseg,nc,status,noverlap=15,kmax=4,fs=8.0_dp,nfft=31, &
            weighting=trim(weightings(weight_index)))
        if (status /= spectral_status_success) error stop 'real stream init failed'
        nwindows = size(batch%times)
        allocate(collected(size(batch%power,1),nwindows,nc,nc),time_collected(nwindows))
        cursor = 0
        first = 1
        do j = 1, size(chunks)
            block = chunks(j)
            last = first+block-1
            partial = stream%push(data(first:last,:))
            if (partial%status /= spectral_status_success) error stop 'real push failed'
            nf = size(partial%times)
            if (nf > 0) then
                collected(:,cursor+1:cursor+nf,:,:) = partial%power
                time_collected(cursor+1:cursor+nf) = partial%times
            end if
            cursor = cursor+nf
            first = last+1
        end do
        if (first /= n+1 .or. cursor /= nwindows) error stop 'real frame count mismatch'
        if (stream%sample_count() /= int(n,int64)) error stop 'real sample count'
        if (stream%window_count() /= int(nwindows,int64)) error stop 'real window count'
        if (maxval(abs(collected-batch%power)) > tol) error stop 'real streamed spectrum mismatch'
        if (maxval(abs(time_collected-batch%times)) > tol) error stop 'real streamed times mismatch'
        avg = stream%mean_matrix()
        if (avg%status /= spectral_status_success) error stop 'real mean status'
        if (maxval(abs(avg%power-sum(batch%power,dim=2)/real(nwindows,dp))) > tol) &
            error stop 'real mean mismatches batch frames'
        deallocate(collected,time_collected)
    end do

    ! Rank-two push accepts axis=2, matching the batch API convention.
    batch = multitaper_spectral_matrix_spectrogram(data,2.5_dp,nseg,noverlap=13, &
        kmax=4,weighting='adaptive',nfft=28)
    call stream%initialize(2.5_dp,nseg,nc,status,noverlap=13,kmax=4,weighting='adaptive',nfft=28)
    nwindows = size(batch%times)
    allocate(collected(size(batch%power,1),nwindows,nc,nc))
    partial = stream%push(transpose(data(1:47,:)),axis=2)
    cursor = size(partial%times)
    if (cursor > 0) collected(:,1:cursor,:,:) = partial%power
    partial = stream%push(transpose(data(48:n,:)),axis=2)
    nf = size(partial%times)
    if (nf > 0) collected(:,cursor+1:cursor+nf,:,:) = partial%power
    if (cursor+nf /= nwindows) error stop 'axis=2 window count'
    if (maxval(abs(collected-batch%power)) > tol) error stop 'axis=2 mismatch'
    j = int(stream%sample_count())
    partial = stream%push(transpose(data(1:4,:)),axis=3)
    if (partial%status /= spectral_status_invalid_input) error stop 'invalid axis accepted'
    if (stream%sample_count() /= int(j,int64)) error stop 'invalid axis consumed samples'
    deallocate(collected)

    ! A complex stream is two-sided, including negative-frequency FFT bins.
    batch = multitaper_spectral_matrix_spectrogram(complex_data,2.5_dp,nseg, &
        noverlap=20,kmax=4,fs=16.0_dp,nfft=32,weighting='adaptive',detrend=.false.)
    call stream%initialize(2.5_dp,nseg,nc,status,noverlap=20,kmax=4,fs=16.0_dp, &
        nfft=32,weighting='adaptive',detrend=.false.,complex_input=.true.)
    if (status /= spectral_status_success) error stop 'complex stream init failed'
    nwindows = size(batch%times)
    allocate(complex_collected(size(batch%power,1),nwindows,nc,nc))
    first = 1
    cursor = 0
    do while (first <= n)
        last = min(n,first+mod(first,11)+2)
        partial = stream%push(complex_data(first:last,:))
        if (partial%status /= spectral_status_success) error stop 'complex push status'
        nf = size(partial%times)
        if (nf > 0) complex_collected(:,cursor+1:cursor+nf,:,:) = partial%power
        cursor = cursor+nf
        first = last+1
    end do
    if (cursor /= nwindows) error stop 'complex window count'
    if (maxval(abs(complex_collected-batch%power)) > tol) error stop 'complex streamed spectrum mismatch'
    if (maxval(abs(partial%frequencies-batch%frequencies)) > tol) error stop 'complex stream frequency mismatch'
    avg = stream%mean_matrix()
    if (maxval(abs(avg%power-sum(batch%power,dim=2)/real(nwindows,dp))) > tol) &
        error stop 'complex running average mismatch'
    deallocate(complex_collected)

    ! The single-channel generics accept rank-one real and complex arrays.
    batch = multitaper_spectral_matrix_spectrogram(reshape(x,[n,1]),2.5_dp,nseg, &
        noverlap=0,kmax=4,weighting='eigen')
    call stream%initialize(2.5_dp,nseg,1,status,noverlap=0,kmax=4,weighting='eigen')
    partial = stream%push(x(1:12))
    if (size(partial%times) /= 0) error stop 'partial window should not emit'
    avg = stream%mean_matrix()
    if (avg%status /= spectral_status_invalid_input) error stop 'empty mean not rejected'
    partial = stream%push(x(13:n))
    if (maxval(abs(partial%power(:,:,1,1)-batch%power(:,:,1,1))) > tol) &
        error stop 'real vector push mismatch'
    call stream%reset()
    if (stream%sample_count() /= 0_int64 .or. stream%window_count() /= 0_int64) &
        error stop 'reset did not clear counters'
    partial = stream%push(x)
    if (maxval(abs(partial%power-batch%power)) > tol) error stop 'reset mismatch'

    batch = multitaper_spectral_matrix_spectrogram(reshape(xc,[n,1]),2.5_dp,nseg, &
        noverlap=23,kmax=4,nfft=27,weighting='adaptive')
    call stream%initialize(2.5_dp,nseg,1,status,noverlap=23,kmax=4,nfft=27, &
        weighting='adaptive',complex_input=.true.)
    partial = stream%push(xc)
    if (size(partial%times) /= size(batch%times)) error stop 'complex vector count'
    if (maxval(abs(partial%power-batch%power)) > tol) error stop 'complex vector mismatch'

    ! Invalid chunks never consume samples or modify the running mean.
    first = int(stream%sample_count())
    partial = stream%push(x(1:9))
    if (partial%status /= spectral_status_invalid_input) error stop 'accepted real samples on complex stream'
    partial = stream%push(reshape(xc(1:9),[9,1]))
    if (partial%status /= spectral_status_success) error stop 'complex rank-two push rejected'
    if (stream%sample_count() /= int(first+9,int64)) error stop 'invalid push changed count'
    partial = stream%push(reshape(xc(1:9),[3,3]))
    if (partial%status /= spectral_status_invalid_input) error stop 'accepted wrong channel count'
    if (stream%sample_count() /= int(first+9,int64)) error stop 'invalid channel push changed count'

    call stream%initialize(2.5_dp,nseg,0,status)
    if (status /= spectral_status_invalid_input) error stop 'accepted zero channels'
    partial = stream%push(x)
    if (partial%status /= spectral_status_invalid_input) error stop 'push on uninitialized stream'
    call stream%initialize(2.5_dp,nseg,1,status,noverlap=nseg)
    if (status /= spectral_status_invalid_input) error stop 'accepted overlap equal window size'
    call stream%initialize(2.5_dp,nseg,1,status,nfft=nseg-1)
    if (status /= spectral_status_invalid_input) error stop 'accepted too-small FFT length'
    call stream%initialize(2.5_dp,nseg,1,status,weighting='unknown')
    if (status /= spectral_status_invalid_input) error stop 'accepted unknown weighting'

    print '(a)', 'Streaming multitaper spectra agree with batch spectral matrices'
end program test_signal_multitaper_stream
