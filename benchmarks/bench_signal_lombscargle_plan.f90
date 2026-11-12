! SPDX-License-Identifier: MIT
! Repeated observation-grid benchmark; compare cached and uncached paths.
program bench_signal_lombscargle_plan
    use scifort_kinds, only : dp
    use scifort_signal_irregular_nd, only : lombscargle_plan, &
        lombscargle_into, lombscargle_plan_into
    implicit none
    integer, parameter :: nt=384, nf=160, nch=16, repeats=12
    real(dp) :: times(nt), frequencies(nf), weights(nt), signals(nt,nch), original(nt,nch)
    real(dp) :: uncached(nf,nch), cached(nf,nch)
    type(lombscargle_plan) :: plan
    real(dp) :: before, after, baseline_sec, planned_sec, initialize_sec
    real(dp) :: checksum_uncached, checksum_cached
    integer :: i, k, r, status

    do i=1,nt
        times(i) = 0.12_dp*(i-1)+0.001_dp*i*i
        weights(i) = 0.1_dp+0.03_dp*mod(i*7,29)
    end do
    do i=1,nf
        frequencies(i) = 0.07_dp+0.044_dp*i
    end do
    do k=1,nch
        do i=1,nt
            signals(i,k) = sin(1.1_dp*times(i)+0.02_dp*k) &
                +0.11_dp*cos(0.7_dp*times(i)+0.4_dp*k)+0.25_dp*k
        end do
    end do
    original = signals
    call cpu_time(before)
    do r=1,repeats
        signals = original+1.0e-6_dp*real(r,dp)
        call lombscargle_into(times,signals,frequencies,uncached, &
            axis=1,weights=weights,floating_mean=.true., &
            normalize='normalize',status=status)
        if (status/=0) error stop 'baseline failed'
    end do
    call cpu_time(after)
    baseline_sec = after-before
    checksum_uncached = sum(uncached)
    call cpu_time(before)
    call plan%initialize(times,frequencies,weights=weights, &
        floating_mean=.true.,status=status)
    if (status/=0) error stop 'plan init failed'
    call cpu_time(after)
    initialize_sec = after-before
    call cpu_time(before)
    do r=1,repeats
        signals = original+1.0e-6_dp*real(r,dp)
        call lombscargle_plan_into(plan,signals,cached,axis=1, &
            normalize='normalize',status=status)
        if (status/=0) error stop 'planned evaluation failed'
    end do
    call cpu_time(after)
    planned_sec = after-before
    checksum_cached = sum(cached)
    if (maxval(abs(cached-uncached)) > 1.0e-10_dp) error stop 'benchmark mismatch'
    print '(a,i0,a,i0,a,i0,a,i0)', 'SAMPLES=',nt,' FREQUENCIES=',nf, &
        ' CHANNELS=',nch,' EVALUATIONS=',repeats
    print '(a,f11.6)', 'BASELINE_SECONDS=',baseline_sec
    print '(a,f11.6)', 'PLAN_SETUP_SECONDS=',initialize_sec
    print '(a,f11.6)', 'PLANNED_SECONDS=',planned_sec
    print '(a,f11.6)', 'SPEEDUP_REPEATED=',baseline_sec/planned_sec
    print '(a,f11.6)', 'SPEEDUP_INCLUDING_SETUP=', &
        baseline_sec/(planned_sec+initialize_sec)
    print '(a,es19.11)', 'UNCACHED_CHECKSUM=',checksum_uncached
    print '(a,es19.11)', 'CACHED_CHECKSUM=',checksum_cached
end program bench_signal_lombscargle_plan
