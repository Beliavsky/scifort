! SPDX-License-Identifier: MIT
! CPU benchmark for many Savitzky-Golay channels; compare original vs shared weights.
program benchmark_savgol_multichannel
    use scifort_kinds, only : dp
    use scifort_signal_savgol, only : savgol_filter
    implicit none
    integer, parameter :: nrep = 30
    real(dp) :: x(24,32,160), t0, t1, checksum
    real(dp), allocatable :: y(:,:,:)
    integer :: i, j, k, rep

    do k = 1, size(x,3)
        do j = 1, size(x,2)
            do i = 1, size(x,1)
                x(i,j,k) = sin(0.07_dp * real(i+j+k,dp)) + 0.01_dp * real(i*j,dp)
            end do
        end do
    end do
    call cpu_time(t0)
    do rep = 1, nrep
        y = savgol_filter(x,21,3,axis=3,mode='interp')
    end do
    call cpu_time(t1)
    checksum = sum(y)
    print '(a,f12.6,a,es24.15)', 'CPU seconds: ', t1-t0, ' checksum: ', checksum
end program benchmark_savgol_multichannel
