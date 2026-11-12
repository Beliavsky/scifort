! SPDX-License-Identifier: MIT
! Batch probe for continuous/discrete state-space frequency responses.
program lti_resolvent_driver
    use scifort_signal, only : dp, freqresp, dfreqresp, frequency_response
    implicit none
    complex(dp), allocatable :: a(:,:), b(:,:), c(:,:), d(:,:)
    real(dp), allocatable :: ar(:,:), br(:,:), cr(:,:), dr(:,:), w(:)
    type(frequency_response) :: response
    integer :: mode, n, nw, i, j, ios
    real(dp) :: re, im
    do
        read(*,*,iostat=ios) mode, n, nw
        if (ios < 0) exit
        if (ios /= 0 .or. n < 0 .or. nw < 0) error stop 'invalid case header'
        allocate(a(n,n),b(n,1),c(1,n),d(1,1),w(nw))
        do j = 1,n
            do i = 1,n
                read(*,*) re, im
                a(i,j) = cmplx(re,im,dp)
            end do
        end do
        do i = 1,n
            read(*,*) re, im
            b(i,1) = cmplx(re,im,dp)
        end do
        do i = 1,n
            read(*,*) re, im
            c(1,i) = cmplx(re,im,dp)
        end do
        read(*,*) re, im
        d(1,1) = cmplx(re,im,dp)
        do i = 1,nw
            read(*,*) w(i)
        end do
        select case(mode)
        case(1)
            ar = real(a,dp)
            br = real(b,dp)
            cr = real(c,dp)
            dr = real(d,dp)
            response = freqresp(ar,br,cr,dr,w=w)
        case(2)
            ar = real(a,dp)
            br = real(b,dp)
            cr = real(c,dp)
            dr = real(d,dp)
            response = dfreqresp(ar,br,cr,dr,w=w)
        case(3)
            response = freqresp(a,b,c,d,w=w)
        case(4)
            response = dfreqresp(a,b,c,d,w=w)
        case default
            error stop 'invalid mode'
        end select
        write(*,'(i0,1x,i0)') response%status,size(response%response)
        do i = 1,size(response%response)
            write(*,'(3(es25.16,1x))') response%frequencies(i), &
                real(response%response(i),dp),aimag(response%response(i))
        end do
        deallocate(a,b,c,d,w)
    end do
end program lti_resolvent_driver
