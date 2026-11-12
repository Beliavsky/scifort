! SPDX-License-Identifier: MIT
! Input: n, output length, channels, case mode followed by channel sample pairs.
! Mode 0: real/time/real callback, 1: complex/time/real callback,
! 2: complex/time/complex callback, 3: complex/freq/real callback,
! 4: complex/freq/complex callback.
program resample_callable_window_driver
    use scifort_kinds, only : dp
    use scifort_signal, only : resample_into
    implicit none
    real(dp), allocatable :: xr(:,:), yr(:,:)
    complex(dp), allocatable :: xc(:,:), yc(:,:)
    integer :: n, m, nc, mode, ios, i, j, st
    real(dp) :: ar, ai

    do
        read(*,*,iostat=ios) n,m,nc,mode
        if (ios < 0) exit
        if (ios /= 0) error stop 'invalid reference-case header'
        if (n < 1 .or. m < 1 .or. nc < 1) error stop 'invalid reference-case dimensions'
        allocate(xr(n,nc),xc(n,nc),yr(m,nc),yc(m,nc))
        do j=1,nc
            do i=1,n
                read(*,*) ar,ai
                xr(i,j) = ar
                xc(i,j) = cmplx(ar,ai,dp)
            end do
        end do
        select case (mode)
        case (0)
            call resample_into(xr,yr,m,axis=1,window_function=real_response,status=st)
            yc = cmplx(yr,0.0_dp,dp)
        case (1)
            call resample_into(xc,yc,m,axis=1,window_function=real_response,status=st)
        case (2)
            call resample_into(xc,yc,m,axis=1,complex_window_function=complex_response,status=st)
        case (3)
            call resample_into(xc,yc,m,axis=1,domain='freq',window_function=real_response,status=st)
        case (4)
            call resample_into(xc,yc,m,axis=1,domain='freq',complex_window_function=complex_response,status=st)
        case default
            error stop 'invalid reference mode'
        end select
        print '(3(i0,1x))', st,m,nc
        do j=1,nc
            do i=1,m
                write(*,'(2(es26.17,1x))') real(yc(i,j),dp),aimag(yc(i,j))
            end do
        end do
        deallocate(xr,xc,yr,yc)
    end do
contains
    function real_response(f) result(w)
        real(dp), intent(in) :: f(:) !! FFT-bin frequencies, cycles/sample
        real(dp) :: w(size(f))
        w = 0.75_dp+0.25_dp*sin(11.0_dp*f)+0.3_dp*f*f
    end function real_response

    function complex_response(f) result(w)
        real(dp), intent(in) :: f(:) !! FFT-bin frequencies, cycles/sample
        complex(dp) :: w(size(f))
        w = cmplx(1.0_dp+0.2_dp*cos(8.0_dp*f),0.4_dp*sin(7.0_dp*f)-0.1_dp*f,dp)
    end function complex_response
end program resample_callable_window_driver
