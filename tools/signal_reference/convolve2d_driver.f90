! SPDX-License-Identifier: MIT
! Machine-readable two-dimensional convolution/correlation driver for SciPy checks.
program convolve2d_driver
    use scifort_kinds, only : dp
    use scifort_signal_convolve2d, only : convolve2d, correlate2d
    implicit none
    integer :: method, n1, n2, m1, m2, md, bd, status, io, i, j
    real(dp) :: re, im, fre, fim
    complex(dp), allocatable :: a(:,:), b(:,:), y(:,:)
    character(len=5) :: modes(3) = [character(len=5) :: 'full', 'same', 'valid']
    character(len=4) :: bounds(3) = [character(len=4) :: 'fill', 'wrap', 'symm']

    do
        read(*,*,iostat=io) method, n1, n2, m1, m2, md, bd
        if (io /= 0) exit
        read(*,*) fre, fim
        allocate(a(n1,n2),b(m1,m2))
        do j = 1, n2
            do i = 1, n1
                read(*,*) re, im
                a(i,j) = cmplx(re,im,dp)
            end do
        end do
        do j = 1, m2
            do i = 1, m1
                read(*,*) re, im
                b(i,j) = cmplx(re,im,dp)
            end do
        end do
        if (method == 1) then
            y = convolve2d(a,b,trim(modes(md)),trim(bounds(bd)),cmplx(fre,fim,dp),status)
        else
            y = correlate2d(a,b,trim(modes(md)),trim(bounds(bd)),cmplx(fre,fim,dp),status)
        end if
        write(*,'(3(i0,1x))') status,size(y,1),size(y,2)
        do j = 1,size(y,2)
            do i = 1,size(y,1)
                write(*,'(2(es26.17e3,1x))') real(y(i,j),dp),aimag(y(i,j))
            end do
        end do
        deallocate(a,b,y)
    end do
end program convolve2d_driver
