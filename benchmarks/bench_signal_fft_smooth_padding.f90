! SPDX-License-Identifier: MIT
! Compare convolution runtimes with the current and prior FFT padding selection.
! Compile and run against each version of signal_fftconvolve.f90 separately.
program bench_signal_fft_smooth_padding
    use scifort_kinds, only : dp
    use scifort_signal_fftconvolve, only : fftconvolve, fftconvolve_nd, oaconvolve
    implicit none
    real(dp), allocatable :: x(:), h(:), y(:)
    real(dp), allocatable :: a(:,:), b(:,:), c(:,:)
    real(dp) :: t0, t1, checksum
    integer :: i, j, k, rep, status

    checksum = 0.0_dp
    do k = 1,4
        select case(k)
        case(1)
            allocate(x(6000),h(71))
        case(2)
            allocate(x(20000),h(101))
        case(3)
            allocate(x(20000),h(257))
        case(4)
            allocate(x(8000),h(503))
        end select
        do i = 1,size(x)
            x(i) = cos(0.013_dp*real(i,dp)) + sin(0.067_dp*real(i,dp))
        end do
        do i = 1,size(h)
            h(i) = sin(0.025_dp*real(i,dp))
        end do
        call cpu_time(t0)
        do rep=1,3
            if (k==1) then
                y = fftconvolve(x,h,status=status)
            else
                y = oaconvolve(x,h,status=status)
            end if
            if (status/=0) error stop 'fft convolution failed'
            checksum = checksum + sum(abs(y))
        end do
        call cpu_time(t1)
        print '(a,i0,a,f10.5,a,es15.7)', '1D case ',k,' seconds=',t1-t0,' checksum=',checksum
        deallocate(x,h,y)
    end do

    do k=1,4
        select case(k)
        case(1)
            allocate(a(127,65),b(43,21),c(169,85))
        case(2)
            allocate(a(270,191),b(53,39),c(322,229))
        case(3)
            allocate(a(384,320),b(73,57),c(456,376))
        case(4)
            allocate(a(101,251),b(29,59),c(129,309))
        end select
        do j=1,size(a,2)
            do i=1,size(a,1)
                a(i,j) = cos(0.017_dp*real(i,dp))+sin(0.011_dp*real(i*j,dp))
            end do
        end do
        do j=1,size(b,2)
            do i=1,size(b,1)
                b(i,j) = sin(0.025_dp*real(i+j,dp))
            end do
        end do
        call cpu_time(t0)
        do rep=1,2
            call fftconvolve_nd(a,b,c,status=status)
            if (status/=0) error stop '2D fft convolution failed'
            checksum = checksum + sum(abs(c))
        end do
        call cpu_time(t1)
        print '(a,i0,a,f10.5,a,es15.7)', '2D case ',k,' seconds=',t1-t0,' checksum=',checksum
        deallocate(a,b,c)
    end do
end program bench_signal_fft_smooth_padding
