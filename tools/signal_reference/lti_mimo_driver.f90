! SPDX-License-Identifier: MIT
! Streamed MIMO state-space frequency and Bode comparison probe.
program lti_mimo_driver
    use scifort_signal, only : dp, freqresp_mimo, dfreqresp_mimo, &
        bode_mimo, dbode_mimo, mimo_frequency_response, mimo_bode_result
    implicit none
    complex(dp), allocatable :: a(:,:), b(:,:), c(:,:), d(:,:)
    real(dp), allocatable :: ar(:,:), br(:,:), cr(:,:), dr(:,:), w(:)
    type(mimo_frequency_response) :: response
    type(mimo_bode_result) :: plot
    integer :: mode, ns, ni, no, nw, i, j, k, ios
    real(dp) :: re, im

    do
        read(*,*,iostat=ios) mode, ns, ni, no, nw
        if (ios < 0) exit
        if (ios /= 0 .or. ns < 0 .or. ni < 1 .or. no < 1 .or. nw < 0) &
            error stop 'invalid case header'
        allocate(a(ns,ns),b(ns,ni),c(no,ns),d(no,ni),w(nw))
        do j = 1,ns
            do i = 1,ns
                read(*,*) re, im
                a(i,j) = cmplx(re,im,dp)
            end do
        end do
        do j = 1,ni
            do i = 1,ns
                read(*,*) re, im
                b(i,j) = cmplx(re,im,dp)
            end do
        end do
        do j = 1,ns
            do i = 1,no
                read(*,*) re, im
                c(i,j) = cmplx(re,im,dp)
            end do
        end do
        do j = 1,ni
            do i = 1,no
                read(*,*) re, im
                d(i,j) = cmplx(re,im,dp)
            end do
        end do
        do i = 1,nw
            read(*,*) w(i)
        end do
        select case(mode)
        case(1,2)
            ar = real(a,dp)
            br = real(b,dp)
            cr = real(c,dp)
            dr = real(d,dp)
            if (mode == 1) then
                response = freqresp_mimo(ar,br,cr,dr,w=w)
                plot = bode_mimo(ar,br,cr,dr,w=w)
            else
                response = dfreqresp_mimo(ar,br,cr,dr,w=w)
                plot = dbode_mimo(ar,br,cr,dr,w=w,dt=0.2_dp)
            end if
        case(3,4)
            if (mode == 3) then
                response = freqresp_mimo(a,b,c,d,w=w)
                plot = bode_mimo(a,b,c,d,w=w)
            else
                response = dfreqresp_mimo(a,b,c,d,w=w)
                plot = dbode_mimo(a,b,c,d,w=w,dt=0.2_dp)
            end if
        case default
            error stop 'invalid mode'
        end select
        write(*,'(3(i0,1x))') response%status, plot%status, size(response%frequencies)
        do k = 1,nw
            do j = 1,ni
                do i = 1,no
                    write(*,'(5(es25.16,1x))') plot%frequencies(k), &
                        real(response%response(i,j,k),dp), aimag(response%response(i,j,k)), &
                        plot%magnitude(i,j,k), plot%phase(i,j,k)
                end do
            end do
        end do
        deallocate(a,b,c,d,w)
    end do
end program lti_mimo_driver
