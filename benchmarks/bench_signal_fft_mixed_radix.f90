! SPDX-License-Identifier: MIT
! Measure reusable mixed-radix FFT versus one-shot FFT on smooth lengths.
program bench_signal_fft_mixed_radix
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : fft_plan, forward_fft, spectral_status_success
    implicit none
    integer, parameter :: lengths(*) = [225,343,625,1000,1728,2401,5040,8192]
    type(fft_plan) :: plan
    complex(dp), allocatable :: x(:), y(:), z(:)
    real(dp) :: t0,t1,plain_sec,planned_sec,tol
    integer :: n, i, j, k, repeats, status

    print '(a)', ' length   unplanned seconds   planned seconds     unplanned / planned'
    do i = 1, size(lengths)
        n = lengths(i)
        repeats = max(12,60000/n)
        x = [(cmplx(sin(0.13_dp*real(j,dp)),cos(0.19_dp*real(j,dp)),dp),j=1,n)]
        call plan%initialize(n,status)
        if (status /= spectral_status_success) error stop 'invalid plan'
        call cpu_time(t0)
        do k=1,repeats
            y=x
            call forward_fft(y)
        end do
        call cpu_time(t1)
        plain_sec=(t1-t0)/real(repeats,dp)
        call cpu_time(t0)
        do k=1,repeats
            z=x
            call plan%execute(z,status)
            if (status /= spectral_status_success) error stop 'planned failure'
        end do
        call cpu_time(t1)
        planned_sec=(t1-t0)/real(repeats,dp)
        tol = 4.0e-12_dp*max(1.0_dp,maxval(abs(y)))
        if (maxval(abs(y-z)) > tol) error stop 'planned and unplanned differ'
        print '(i7,2(4x,es15.6),6x,f9.2)',n,plain_sec,planned_sec, &
             plain_sec/max(planned_sec,tiny(1.0_dp))
    end do
end program bench_signal_fft_mixed_radix
