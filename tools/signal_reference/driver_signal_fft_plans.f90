! SPDX-License-Identifier: MIT
! Text protocol for direct SciPy comparison of repeated planned FFTs.
program driver_signal_fft_plans
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : fft_plan, forward_fft, forward_fft_planned, spectral_status_success
    implicit none
    type(fft_plan) :: plan
    complex(dp), allocatable :: x(:), original(:)
    real(dp) :: re, im
    integer :: n, i, k, repeats, status, ios

    read(*,*,iostat=ios) n,repeats
    if (ios /= 0 .or. n < 1 .or. n > 32768 .or. repeats < 1 .or. repeats > 20) error stop 'invalid input'
    allocate(x(n),original(n))
    do i=1,n
        read(*,*,iostat=ios) re,im
        if (ios /= 0) error stop 'invalid coefficient'
        original(i)=cmplx(re,im,dp)
    end do
    call plan%initialize(n,status)
    if (status /= spectral_status_success) error stop 'invalid plan'
    do k=1,repeats
        x=original*cmplx(real(k,dp),0.3_dp*real(k,dp),dp)
        call forward_fft_planned(x,plan,status)
        if (status /= spectral_status_success) error stop 'fft status'
    end do
    do i=1,n
        write(*,'(es26.17,1x,es26.17)') real(x(i),dp),aimag(x(i))
    end do
end program driver_signal_fft_plans
