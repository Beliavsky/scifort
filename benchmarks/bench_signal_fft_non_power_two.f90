! SPDX-License-Identifier: MIT
! Compare the improved shared FFT to the previous direct-DFT fallback.
program bench_signal_fft_non_power_two
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : forward_fft
    implicit none
    integer, parameter :: lengths(*) = [257,1009,4093]
    complex(dp), allocatable :: x(:), y(:), direct(:)
    real(dp) :: t0, t1, fast_seconds, slow_seconds, checksum, pi
    integer :: i, j, k, n, repeats_fast

    pi = acos(-1.0_dp)
    print '(a)', '  length    FFT seconds   DFT seconds   speedup'
    do i = 1, size(lengths)
        n = lengths(i)
        allocate(x(n),y(n),direct(n))
        x = [(cmplx(sin(real(j,dp)*0.13_dp),cos(real(j,dp)*0.19_dp),dp),j=1,n)]
        repeats_fast = max(20, 500000/n)
        call cpu_time(t0)
        do k = 1, repeats_fast
            y = x
            call forward_fft(y)
        end do
        call cpu_time(t1)
        fast_seconds = (t1-t0)/real(repeats_fast,dp)
        call cpu_time(t0)
        direct = x
        call old_direct_dft(direct,pi)
        call cpu_time(t1)
        slow_seconds = t1-t0
        checksum = maxval(abs(y-direct))
        if (checksum > 5.0e-10_dp*max(1.0_dp,maxval(abs(direct)))) error stop 'benchmark mismatch'
        print '(i8,2(3x,es12.4),3x,f8.1)', n,fast_seconds,slow_seconds,slow_seconds/fast_seconds
        deallocate(x,y,direct)
    end do
contains
    ! Preserved O(n**2) algorithm from the previous SciFort forward_fft.
    subroutine old_direct_dft(values, pi_local)
        complex(dp), intent(inout) :: values(:) !! input and output DFT
        real(dp), intent(in) :: pi_local !! circle constant
        complex(dp) :: temp(size(values))
        real(dp) :: angle
        integer :: n,j,k
        n = size(values)
        do k = 0, n-1
            temp(k+1) = (0.0_dp,0.0_dp)
            do j = 0, n-1
                angle = -2.0_dp*pi_local*real(k,dp)*real(j,dp)/real(n,dp)
                temp(k+1) = temp(k+1) + values(j+1)*cmplx(cos(angle),sin(angle),dp)
            end do
        end do
        values = temp
    end subroutine old_direct_dft
end program bench_signal_fft_non_power_two
