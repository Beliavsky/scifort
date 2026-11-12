! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Compare the end-to-end cost of repeated multitaper stream frames for
! composite and prime FFT lengths. Run this same driver against the pre-plan
! and plan-reuse versions for a fair A/B comparison.
program bench_signal_multitaper_fft_plan_reuse
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : spectral_status_success
    use scifort_signal_multitaper, only : multitaper_stream, &
        multitaper_spectral_matrix_spectrogram_result
    implicit none
    integer, parameter :: nsample = 2480, nchannels = 2, segment = 80, nrepeat = 10
    real(dp) :: signal(nsample,nchannels), t0, t1, elapsed, checksum
    integer :: i, channel, nfft, k, status
    character(20) :: arg
    type(multitaper_stream) :: estimator
    type(multitaper_spectral_matrix_spectrogram_result) :: result

    nfft = 263
    if (command_argument_count() > 0) then
        call get_command_argument(1,arg)
        read(arg,*) nfft
    end if
    do channel = 1,nchannels
        do i = 1,nsample
            signal(i,channel) = sin(0.028_dp*real(i*channel,dp)) + &
                0.3_dp*cos(0.093_dp*real(i+channel,dp))
        end do
    end do
    call estimator%initialize(3.5_dp,segment,nchannels,status,noverlap=60,kmax=5, &
        nfft=nfft,weighting='adaptive')
    if (status /= spectral_status_success) error stop 'benchmark init failed'
    checksum = 0.0_dp
    call cpu_time(t0)
    do k = 1,nrepeat
        call estimator%reset()
        result = estimator%push(signal)
        if (result%status /= spectral_status_success) error stop 'benchmark push failed'
        checksum = checksum + sum(real(result%power(:,1,1,1),dp))
    end do
    call cpu_time(t1)
    elapsed = t1-t0
    print '(a,i0,a,f11.6,a,es16.8)', 'nfft=',nfft,' seconds=',elapsed,' checksum=',checksum
end program bench_signal_multitaper_fft_plan_reuse
