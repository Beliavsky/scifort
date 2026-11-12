! SPDX-License-Identifier: MIT
! Cross-language rank-15 FFT convolution probe; singleton dimensions permit
! arbitrary broadcast patterns without allocating large dense arrays.
program fftconvolve_assumed_rank_driver
    use scifort_kinds, only : dp
    use scifort_signal_fftconvolve, only : fftconvolve_nd
    implicit none
    complex(dp), allocatable :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
    complex(dp), allocatable :: b(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
    complex(dp), allocatable :: y(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
    real(dp), allocatable :: yr(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
    complex(dp), allocatable :: aa(:), bb(:), flat(:)
    integer :: n(15), m(15), o(15), axis(15), naxis, mask, mid, kind
    integer :: code, ios, i, k
    real(dp) :: re, im
    character(5), parameter :: modes(3) = ['full ', 'same ', 'valid']

    do
        read(*,*,iostat=ios) kind, mid, mask
        if (ios /= 0) exit
        read(*,*) n
        read(*,*) m
        allocate(aa(product(n)),bb(product(m)))
        do i = 1,size(aa)
            read(*,*) re, im
            aa(i) = cmplx(re,im,dp)
        end do
        do i = 1,size(bb)
            read(*,*) re, im
            bb(i) = cmplx(re,im,dp)
        end do
        naxis = 0
        o = max(n,m)
        do i = 1,15
            if (.not. btest(mask,i-1)) cycle
            naxis = naxis+1
            axis(naxis) = i
            select case (mid)
            case (1)
                o(i) = n(i)+m(i)-1
            case (2)
                o(i) = n(i)
            case default
                o(i) = abs(n(i)-m(i))+1
            end select
        end do
        if (mid == 2) then
            do i = 1,15
                if (.not. btest(mask,i-1)) o(i) = n(i)
            end do
        end if
        allocate(a(n(1),n(2),n(3),n(4),n(5),n(6),n(7),n(8),n(9),n(10), &
            n(11),n(12),n(13),n(14),n(15)))
        allocate(b(m(1),m(2),m(3),m(4),m(5),m(6),m(7),m(8),m(9),m(10), &
            m(11),m(12),m(13),m(14),m(15)))
        a = reshape(aa,n)
        b = reshape(bb,m)
        deallocate(aa,bb)
        allocate(y(o(1),o(2),o(3),o(4),o(5),o(6),o(7),o(8),o(9),o(10), &
            o(11),o(12),o(13),o(14),o(15)))
        select case (kind)
        case (0)
            allocate(yr(o(1),o(2),o(3),o(4),o(5),o(6),o(7),o(8),o(9),o(10), &
                o(11),o(12),o(13),o(14),o(15)))
            call fftconvolve_nd(real(a,dp),real(b,dp),yr,modes(mid),axis(:naxis),code)
            y = cmplx(yr,0.0_dp,dp)
            deallocate(yr)
        case (1)
            call fftconvolve_nd(a,b,y,modes(mid),axis(:naxis),code)
        case (2)
            call fftconvolve_nd(real(a,dp),b,y,modes(mid),axis(:naxis),code)
        case (3)
            call fftconvolve_nd(a,real(b,dp),y,modes(mid),axis(:naxis),code)
        end select
        print '(*(i0,1x))', code,o
        flat = reshape(y,[product(o)])
        do k = 1,size(flat)
            print '(es25.17,1x,es25.17)', real(flat(k),dp),aimag(flat(k))
        end do
        deallocate(a,b,y,flat)
    end do
end program fftconvolve_assumed_rank_driver
