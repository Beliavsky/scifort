! SPDX-License-Identifier: MIT
! Streaming cross-language probe for rank-two and rank-three FFT convolution.
program fftconvolve_nd_driver
    use scifort_kinds, only : dp
    use scifort_signal_fftconvolve, only : fftconvolve
    implicit none
    integer :: rank, code, mid, mask, n(3), m(3), axis(3), naxis, st, io
    integer :: i,j,k, q
    real(dp) :: re, im
    complex(dp), allocatable :: a2(:,:), b2(:,:), y2(:,:)
    complex(dp), allocatable :: a3(:,:,:), b3(:,:,:), y3(:,:,:)
    character(5), parameter :: modes(3) = ['full ', 'same ', 'valid']

    do
        read(*,*,iostat=io) rank,code,mid,mask,n,m
        if (io /= 0) exit
        naxis = 0
        do i = 1, rank
            if (.not. btest(mask,i-1)) cycle
            naxis = naxis + 1
            axis(naxis) = i
        end do
        if (rank == 2) then
            allocate(a2(n(1),n(2)),b2(m(1),m(2)))
            do j = 1, n(2)
                do i = 1, n(1)
                    read(*,*) re, im
                    a2(i,j) = cmplx(re,im,dp)
                end do
            end do
            do j = 1, m(2)
                do i = 1, m(1)
                    read(*,*) re, im
                    b2(i,j) = cmplx(re,im,dp)
                end do
            end do
            select case (code)
            case (0)
                y2 = fftconvolve(real(a2,dp),real(b2,dp),modes(mid),st,axis(:naxis))
            case (1)
                y2 = fftconvolve(a2,b2,modes(mid),st,axis(:naxis))
            case (2)
                y2 = fftconvolve(real(a2,dp),b2,modes(mid),st,axis(:naxis))
            case default
                y2 = fftconvolve(a2,real(b2,dp),modes(mid),st,axis(:naxis))
            end select
            print '(3(i0,1x))', st,size(y2,1),size(y2,2)
            do j = 1, size(y2,2)
                do i = 1, size(y2,1)
                    print '(es25.17,1x,es25.17)', real(y2(i,j),dp),aimag(y2(i,j))
                end do
            end do
            deallocate(a2,b2,y2)
        else
            allocate(a3(n(1),n(2),n(3)),b3(m(1),m(2),m(3)))
            do k = 1, n(3)
                do j = 1, n(2)
                    do i = 1, n(1)
                        read(*,*) re, im
                        a3(i,j,k) = cmplx(re,im,dp)
                    end do
                end do
            end do
            do k = 1, m(3)
                do j = 1, m(2)
                    do i = 1, m(1)
                        read(*,*) re, im
                        b3(i,j,k) = cmplx(re,im,dp)
                    end do
                end do
            end do
            select case (code)
            case (0)
                y3 = fftconvolve(real(a3,dp),real(b3,dp),modes(mid),st,axis(:naxis))
            case (1)
                y3 = fftconvolve(a3,b3,modes(mid),st,axis(:naxis))
            case (2)
                y3 = fftconvolve(real(a3,dp),b3,modes(mid),st,axis(:naxis))
            case default
                y3 = fftconvolve(a3,real(b3,dp),modes(mid),st,axis(:naxis))
            end select
            print '(4(i0,1x))', st,size(y3,1),size(y3,2),size(y3,3)
            do k = 1, size(y3,3)
                do j = 1, size(y3,2)
                    do i = 1, size(y3,1)
                        print '(es25.17,1x,es25.17)', real(y3(i,j,k),dp),aimag(y3(i,j,k))
                    end do
                end do
            end do
            deallocate(a3,b3,y3)
        end if
    end do
end program fftconvolve_nd_driver
