! SPDX-License-Identifier: MIT
! Illustrative single-process throughput comparison; not a test assertion.
program benchmark_signal_harmonic_plan
    use scifort_kinds, only : dp
    use scifort_signal_irregular, only : harmonic_fit_result, harmonic_fit_plan, &
        lombscargle_harmonic_prepare, lombscargle_harmonic_fit_batch, lombscargle_harmonic_fit
    implicit none

    integer, parameter :: n = 160, nfreq = 64, nchannels = 20, repetitions = 10
    integer :: i, j, iter
    real(dp) :: x(n), y(n,nchannels), weights(n), freq(nfreq), before, after
    real(dp) :: direct_seconds, planned_seconds
    type(harmonic_fit_plan) :: plan
    type(harmonic_fit_result) :: one
    type(harmonic_fit_result), allocatable :: fitted(:)

    do i = 1, n
        x(i) = 0.137_dp*real(i,dp)+0.05_dp*cos(real(i*i,dp))
        weights(i) = 0.1_dp+real(mod(i,11),dp)/10.0_dp
    end do
    do j = 1, nchannels
        y(:,j) = real(j,dp)+cos(0.1_dp*real(j,dp)*x)+ &
            0.08_dp*sin(real(j,dp)*x)
    end do
    do j = 1, nfreq
        freq(j) = 0.15_dp+0.019_dp*real(j,dp)
    end do

    call cpu_time(before)
    do iter = 1, repetitions
        do j = 1, nchannels
            one = lombscargle_harmonic_fit(x,y(:,j),freq,3, &
                weights=weights,floating_mean=.true.)
        end do
    end do
    call cpu_time(after)
    direct_seconds = after-before

    call cpu_time(before)
    do iter = 1, repetitions
        plan = lombscargle_harmonic_prepare(x,freq,3, &
            weights=weights,floating_mean=.true.)
        fitted = lombscargle_harmonic_fit_batch(plan,y)
    end do
    call cpu_time(after)
    planned_seconds = after-before

    if (one%status /= 0 .or. fitted(1)%status /= 0) error stop 'benchmark fit failed'
    print '(a,es12.4)', 'direct CPU seconds:  ',direct_seconds
    print '(a,es12.4)', 'planned CPU seconds: ',planned_seconds
    print '(a,f9.2)', 'speedup factor:     ',direct_seconds/max(planned_seconds,tiny(1.0_dp))
end program benchmark_signal_harmonic_plan
