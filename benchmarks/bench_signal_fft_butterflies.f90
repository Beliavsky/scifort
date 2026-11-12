! SPDX-License-Identifier: MIT
! Planned mixed-radix FFT timings; run against old/new signal_spectral
! builds with exactly the same flags. Report minimum of three CPU-time runs.
program bench_signal_fft_butterflies
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : fft_plan, spectral_status_success
    implicit none
    integer, parameter :: lengths(*) = [243,343,625,729,1000,1225,1728,2401,3125,5040]
    type(fft_plan) :: plan
    complex(dp), allocatable :: x(:), z(:)
    real(dp) :: first, last, best, checksum
    integer :: i, j, k, n, nrepeat, pass, status

    print '(a)', ' length,seconds_per_transform'
    do i = 1, size(lengths)
        n = lengths(i)
        allocate(x(n),z(n))
        x = [(cmplx(sin(0.031_dp*real(j,dp)),cos(0.021_dp*real(j,dp)),dp),j=1,n)]
        call plan%initialize(n,status)
        if (status /= spectral_status_success) error stop 'FFT benchmark initialization'
        nrepeat = max(320,1000000/n)
        z = x
        call plan%execute(z,status)
        best = huge(1.0_dp)
        checksum = 0.0_dp
        do pass = 1, 3
            call cpu_time(first)
            do k = 1, nrepeat
                z = x
                call plan%execute(z,status)
                if (status /= spectral_status_success) error stop 'FFT benchmark execute'
                checksum = checksum + real(z(1),dp)
            end do
            call cpu_time(last)
            best = min(best,(last-first)/real(nrepeat,dp))
        end do
        if (.not.(checksum > -huge(1.0_dp))) error stop 'FFT benchmark checksum'
        print '(i0,a,es16.8)', n, ',', best
        call plan%clear()
        deallocate(x,z)
    end do
end program bench_signal_fft_butterflies
