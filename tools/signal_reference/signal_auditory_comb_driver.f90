! SPDX-License-Identifier: MIT
! Compare gammatone and comb filter coefficients with SciPy.
program signal_auditory_comb_driver
    use scifort_kinds, only : dp
    use scifort_signal_iir_extra, only : gammatone, iircomb
    use scifort_signal_design, only : butter_result
    implicit none
    character(1) :: command
    character(16) :: mode
    real(dp) :: frequency, factor, fs
    integer :: iostat, order, numtaps, passflag, n_b, n_a
    type(butter_result) :: filt

    do
        read(*,*,iostat=iostat) command, mode, frequency, factor, fs, order, numtaps, passflag
        if (iostat < 0) exit
        if (iostat /= 0) error stop 'reference input read failed'
        select case(command)
        case ('c')
            filt = iircomb(frequency, factor, trim(mode), fs=fs, pass_zero=(passflag /= 0))
        case ('g')
            filt = gammatone(frequency, trim(mode), order=order, numtaps=numtaps, fs=fs)
        case default
            error stop 'unknown reference command'
        end select
        n_b = 0
        n_a = 0
        if (allocated(filt%b)) n_b = size(filt%b)
        if (allocated(filt%a)) n_a = size(filt%a)
        print '(3(i0,1x))', filt%status, n_b, n_a
        if (n_b > 0) then
            print '(*(es27.18,1x))', filt%b
        else
            print '(a)', '-'
        end if
        if (n_a > 0) then
            print '(*(es27.18,1x))', filt%a
        else
            print '(a)', '-'
        end if
    end do
end program signal_auditory_comb_driver
