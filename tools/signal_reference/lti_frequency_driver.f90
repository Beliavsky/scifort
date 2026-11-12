! SPDX-License-Identifier: MIT
! Batch probe for comparing SciFort LTI frequency responses to SciPy.
program lti_frequency_driver
    use scifort_signal, only : dp, freqresp, dfreqresp, bode, dbode, &
        freqresp_zpk, dfreqresp_zpk, bode_result, frequency_response
    implicit none
    integer :: mode, nb, na, nw, j, ios
    real(dp) :: dt, x, y
    real(dp), allocatable :: w(:)
    complex(dp), allocatable :: b(:), a(:)
    complex(dp) :: gain
    type(frequency_response) :: resp
    type(bode_result) :: diagram
    do
        read(*,*,iostat=ios) mode, nb, na, nw, dt
        if (ios /= 0) exit
        allocate(b(nb),a(na),w(nw))
        do j=1,nb
            read(*,*) x, y
            b(j) = cmplx(x,y,dp)
        end do
        do j=1,na
            read(*,*) x, y
            a(j) = cmplx(x,y,dp)
        end do
        do j=1,nw
            read(*,*) w(j)
        end do
        select case(mode)
        case(1)
            resp = freqresp(b,a,w=w)
        case(2)
            resp = dfreqresp(b,a,w=w,dt=dt)
        case(3)
            gain = cmplx(dt,0.0_dp,dp)
            resp = freqresp_zpk(b,a,gain,w=w)
        case(4)
            gain = cmplx(dt,0.0_dp,dp)
            resp = dfreqresp_zpk(b,a,gain,w=w)
        case(5)
            diagram = bode(b,a,w=w)
        case(6)
            diagram = dbode(b,a,w=w,dt=dt)
        end select
        if (mode < 5) then
            write(*,'(i0,1x,i0)') resp%status, size(resp%response)
            if (allocated(resp%response)) then
                do j=1,size(resp%response)
                    write(*,'(3(es24.16,1x))') resp%frequencies(j), &
                        real(resp%response(j),dp), aimag(resp%response(j))
                end do
            end if
        else
            write(*,'(i0,1x,i0)') diagram%status, size(diagram%phase)
            if (allocated(diagram%phase)) then
                do j=1,size(diagram%phase)
                    write(*,'(3(es24.16,1x))') diagram%frequencies(j), &
                        diagram%magnitude(j), diagram%phase(j)
                end do
            end if
        end if
        deallocate(b,a,w)
    end do
end program lti_frequency_driver
