! SPDX-License-Identifier: MIT
! Text protocol for independent SciPy waveform comparison. See verify_waveforms.py.
program waveforms_driver
    use scifort_kinds, only : dp
    use scifort_signal_waveforms, only : chirp, chirp_complex, sweep_poly, &
        gausspulse, gausspulse_cutoff
    implicit none

    real(dp), allocatable :: t(:), poly(:), y(:), q(:), envelope(:)
    complex(dp), allocatable :: z(:)
    real(dp) :: f0, f1, t1, phi, fc, bw, bwr, tpr, cutoff
    integer :: mode, n, degree, vertex, method_id, i, status, ios
    character(len=16) :: method

    do
        read(*,*,iostat=ios) mode
        if (ios /= 0) exit
        select case(mode)
        case (1,2)
            read(*,*) n, f0, t1, f1, phi, vertex, method_id
            allocate(t(n))
            if (n > 0) read(*,*) t
            select case(method_id)
            case (1)
                method = 'linear'
            case (2)
                method = 'quadratic'
            case (3)
                method = 'logarithmic'
            case (4)
                method = 'hyperbolic'
            case default
                error stop 'waveforms_driver invalid method'
            end select
            if (mode == 1) then
                y = chirp(t, f0, t1, f1, method=trim(method), phi=phi, &
                    vertex_zero=vertex /= 0, status=status)
                print '(i0,1x,i0)', status, size(y)
                do i = 1, size(y)
                    print '(es26.17e3)', y(i)
                end do
                deallocate(y)
            else
                z = chirp_complex(t, f0, t1, f1, method=trim(method), phi=phi, &
                    vertex_zero=vertex /= 0, status=status)
                print '(i0,1x,i0)', status, size(z)
                do i = 1, size(z)
                    print '(es26.17e3,1x,es26.17e3)', real(z(i),dp), aimag(z(i))
                end do
                deallocate(z)
            end if
            deallocate(t)
        case (3)
            read(*,*) n, degree, phi
            allocate(poly(degree),t(n))
            if (degree > 0) read(*,*) poly
            if (n > 0) read(*,*) t
            y = sweep_poly(t,poly,phi=phi,status=status)
            print '(i0,1x,i0)', status, size(y)
            do i = 1, size(y)
                print '(es26.17e3)', y(i)
            end do
            deallocate(poly,t,y)
        case (4)
            read(*,*) n, fc, bw, bwr
            allocate(t(n))
            if (n > 0) read(*,*) t
            y = gausspulse(t, fc=fc, bw=bw, bwr=bwr, quadrature=q, &
                envelope=envelope, status=status)
            print '(i0,1x,i0)', status, size(y)
            do i = 1, size(y)
                print '(es26.17e3,1x,es26.17e3,1x,es26.17e3)', y(i),q(i),envelope(i)
            end do
            deallocate(t,y,q,envelope)
        case (5)
            read(*,*) fc, bw, bwr, tpr
            cutoff = gausspulse_cutoff(fc=fc,bw=bw,bwr=bwr,tpr=tpr,status=status)
            print '(i0,1x,es26.17e3)', status,cutoff
        case default
            error stop 'waveforms_driver invalid mode'
        end select
    end do
end program waveforms_driver
