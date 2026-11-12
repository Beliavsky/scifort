! SPDX-License-Identifier: MIT
! Evaluates SciFort's scalar band_stop_obj on line-oriented cases from stdin.
program band_stop_obj_driver
    use scifort_kinds, only : dp
    use scifort_signal_band_stop_obj, only : band_stop_obj
    implicit none
    real(dp) :: wp, passb(2), stopb(2), gpass, gstop, order
    integer :: ind, status, io
    character(16) :: family

    do
        read(*,*,iostat=io) wp, ind, passb, stopb, gpass, gstop, family
        if (io < 0) exit
        if (io /= 0) error stop 'invalid reference input'
        order = band_stop_obj(wp,ind,passb,stopb,gpass,gstop,trim(family),status)
        write(*,'(ES25.16E3,1X,I0)') order, status
    end do
end program band_stop_obj_driver
