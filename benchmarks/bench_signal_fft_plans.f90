! SPDX-License-Identifier: MIT
! Amortized cost of reusable Bluestein plans compared to the unplanned FFT.
program bench_signal_fft_plans
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : fft_plan, forward_fft, forward_fft_planned, spectral_status_success
    implicit none
    integer, parameter :: lengths(*) = [257,1009,4093,8191]
    type(fft_plan) :: plan
    complex(dp), allocatable :: x(:), y(:), z(:)
    real(dp) :: t0, t1, unplanned_seconds, planned_seconds, tol
    integer :: i, j, k, n, status, repetitions

    print '(a)', '    length    unplanned sec       planned sec     speedup'
    do i = 1, size(lengths)
        n = lengths(i)
        repetitions = max(20, 200000/n)
        x = [(cmplx(sin(real(j,dp)*0.13_dp),cos(real(j,dp)*0.19_dp),dp),j=1,n)]
        call plan%initialize(n,status)
        if (status /= spectral_status_success) error stop 'plan initialization failure'
        call cpu_time(t0)
        do k = 1,repetitions
            y = x
            call forward_fft(y)
        end do
        call cpu_time(t1)
        unplanned_seconds = (t1-t0)/real(repetitions,dp)
        call cpu_time(t0)
        do k = 1,repetitions
            z = x
            call forward_fft_planned(z,plan,status)
            if (status /= spectral_status_success) error stop 'planned FFT status'
        end do
        call cpu_time(t1)
        planned_seconds = (t1-t0)/real(repetitions,dp)
        tol = 2.0e-12_dp*max(1.0_dp,maxval(abs(y)))
        if (maxval(abs(y-z)) > tol) error stop 'benchmark spectra differ'
        print '(i10,2(4x,es13.6),4x,f9.2)',n,unplanned_seconds,planned_seconds, &
            unplanned_seconds/max(planned_seconds,tiny(1.0_dp))
    end do
end program bench_signal_fft_plans
