! SPDX-License-Identifier: MIT
! Text-based independent FFT comparison driver, input length then complex samples.
program driver_signal_fft_bluestein
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : forward_fft
    implicit none
    complex(dp), allocatable :: x(:)
    real(dp) :: xr, xi
    integer :: n, i, status
    read(*,*,iostat=status) n
    if (status /= 0 .or. n < 1) error stop 'invalid input FFT length'
    allocate(x(n))
    do i = 1, n
        read(*,*,iostat=status) xr,xi
        if (status /= 0) error stop 'invalid complex sample'
        x(i) = cmplx(xr,xi,dp)
    end do
    call forward_fft(x)
    do i = 1, n
        print '(es26.17e3,1x,es26.17e3)',real(x(i),dp),aimag(x(i))
    end do
end program driver_signal_fft_bluestein
