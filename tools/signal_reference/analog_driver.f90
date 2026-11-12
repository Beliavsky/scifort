! SPDX-License-Identifier: MIT
! Stdin/stdout comparison driver for analog and digital ZPK frequency responses.
program analog_driver
    use scifort_kinds, only : dp
    use scifort_signal_analog, only : freqs, freqs_zpk
    use scifort_signal_design, only : freqz_zpk, frequency_response
    implicit none
    character(8) :: command
    integer :: nb, na, nw, grid_mode, whole_flag, ios, i
    real(dp) :: fs, gain_real, gain_im
    real(dp), allocatable :: data(:,:), freqs_grid(:)
    complex(dp), allocatable :: numerator(:), denominator(:)
    type(frequency_response) :: res

    do
        read(*,*,iostat=ios) command, nb, na, nw, grid_mode, &
            gain_real, gain_im, fs, whole_flag
        if (ios /= 0) exit
        allocate(numerator(nb),denominator(na),freqs_grid(0))
        if (nb > 0) then
            allocate(data(2,nb))
            read(*,*) data
            numerator = cmplx(data(1,:),data(2,:),dp)
            deallocate(data)
        end if
        if (na > 0) then
            allocate(data(2,na))
            read(*,*) data
            denominator = cmplx(data(1,:),data(2,:),dp)
            deallocate(data)
        end if
        if (grid_mode == 0) then
            deallocate(freqs_grid)
            allocate(freqs_grid(nw))
            if (nw > 0) read(*,*) freqs_grid
        end if
        select case (command)
        case ('ATF')
            if (grid_mode == 0) then
                res = freqs(real(numerator,dp),real(denominator,dp),w=freqs_grid)
            else
                res = freqs(real(numerator,dp),real(denominator,dp),worn=nw)
            end if
        case ('CTF')
            if (grid_mode == 0) then
                res = freqs(numerator,denominator,w=freqs_grid)
            else
                res = freqs(numerator,denominator,worn=nw)
            end if
        case ('AZPK')
            if (grid_mode == 0) then
                res = freqs_zpk(numerator,denominator,gain_real,w=freqs_grid)
            else
                res = freqs_zpk(numerator,denominator,gain_real,worn=nw)
            end if
        case ('CZPK')
            if (grid_mode == 0) then
                res = freqs_zpk(numerator,denominator,cmplx(gain_real,gain_im,dp), &
                    w=freqs_grid)
            else
                res = freqs_zpk(numerator,denominator,cmplx(gain_real,gain_im,dp),worn=nw)
            end if
        case ('DZPK')
            if (grid_mode == 0) then
                res = freqz_zpk(numerator,denominator,gain_real,fs=fs,w=freqs_grid)
            else
                res = freqz_zpk(numerator,denominator,gain_real,fs=fs, &
                    whole=whole_flag==1,worn=nw)
            end if
        end select
        print '(i0,1x,i0)',res%status,size(res%response)
        do i = 1, size(res%response)
            print '(3(es27.17e3,1x))',res%frequencies(i), &
                real(res%response(i),dp),aimag(res%response(i))
        end do
        deallocate(numerator,denominator,freqs_grid)
    end do
end program analog_driver
