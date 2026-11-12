! SPDX-License-Identifier: MIT
! SciPy reference driver for FFT and overlap-add convolution.
program fftconvolve_driver
    use scifort_signal
    implicit none
    character(16) :: method, mode
    integer :: n,m,as_complex,i,status
    real(dp) :: re,im
    complex(dp), allocatable :: x(:),y(:),z(:)
    real(dp), allocatable :: zr(:)

    read(*,'(a)') method
    read(*,'(a)') mode
    read(*,*) n,m,as_complex
    allocate(x(n),y(m))
    do i=1,n
        read(*,*) re,im
        x(i) = cmplx(re,im,dp)
    end do
    do i=1,m
        read(*,*) re,im
        y(i) = cmplx(re,im,dp)
    end do
    status = convolution_status_invalid_input
    if (as_complex/=0) then
        if (trim(method)=='fftconvolve') then
            z = fftconvolve(x,y,mode=trim(mode),status=status)
        else
            z = oaconvolve(x,y,mode=trim(mode),status=status)
        end if
    else
        if (trim(method)=='fftconvolve') then
            zr = fftconvolve(real(x,dp),real(y,dp),mode=trim(mode),status=status)
        else
            zr = oaconvolve(real(x,dp),real(y,dp),mode=trim(mode),status=status)
        end if
        z = cmplx(zr,0.0_dp,dp)
    end if
    write(*,*) status, size(z)
    do i=1,size(z)
        write(*,'(2(es26.17e3,1x))') real(z(i),dp),aimag(z(i))
    end do
end program fftconvolve_driver
