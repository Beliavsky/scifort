! SPDX-License-Identifier: MIT
! Test-only stdin/stdout driver for tools/signal_reference/verify_design.py.
program design_driver
    use scifort_signal_design, only : firwin, butter, freqz, sosfreqz, &
        butter_result, frequency_response
    use scifort_kinds, only : dp
    implicit none
    type(butter_result) :: design
    type(frequency_response) :: response
    character(8) :: command, band, window
    real(dp), allocatable :: cutoff(:), taps(:), frequencies(:), pairs(:,:)
    complex(dp), allocatable :: b(:), a(:), sos(:,:)
    integer :: n, na, nf, ios, status, scale_flag, j
    real(dp) :: fs

    do
        read(*,*,iostat=ios) command, n, na, nf, fs
        if (ios/=0) exit
        select case (command)
        case ('FIR')
            allocate(cutoff(na))
            read(*,*) band, window, scale_flag
            read(*,*) cutoff
            taps = firwin(n,cutoff,pass_zero=trim(band),window=trim(window), &
                scale=scale_flag==1,fs=fs,status=status)
            print '(i0,1x,i0)',status,size(taps)
            if (size(taps)>0) print '(*(es26.17e3,1x))',taps
            deallocate(cutoff)
        case ('BUT')
            allocate(cutoff(na))
            read(*,*) band
            read(*,*) cutoff
            design = butter(n,cutoff,btype=trim(band),fs=fs)
            print '(i0,1x,i0)',design%status,size(design%b)
            if (design%status==0) then
                print '(*(es26.17e3,1x))',design%b
                print '(*(es26.17e3,1x))',design%a
                do j=1,size(design%sos,1)
                    print '(*(es26.17e3,1x))',design%sos(j,:)
                end do
            end if
            deallocate(cutoff)
        case ('FREQ', 'CFREQ')
            allocate(frequencies(nf),pairs(2,n),b(n))
            read(*,*) frequencies
            read(*,*) pairs
            b=cmplx(pairs(1,:),pairs(2,:),dp)
            deallocate(pairs)
            allocate(pairs(2,na),a(na))
            read(*,*) pairs
            a=cmplx(pairs(1,:),pairs(2,:),dp)
            if (command=='FREQ') then
                response=freqz(real(b,dp),real(a,dp),w=frequencies,fs=fs)
            else
                response=freqz(b,a,w=frequencies,fs=fs)
            end if
            call print_response(response)
            deallocate(pairs,frequencies,b,a)
        case ('SOS', 'CSOS')
            allocate(frequencies(nf),pairs(2,6*n),sos(n,6))
            read(*,*) frequencies
            read(*,*) pairs
            sos=transpose(reshape(cmplx(pairs(1,:),pairs(2,:),dp),[6,n]))
            if (command=='SOS') then
                response=sosfreqz(real(sos,dp),w=frequencies,fs=fs)
            else
                response=sosfreqz(sos,w=frequencies,fs=fs)
            end if
            call print_response(response)
            deallocate(frequencies,pairs,sos)
        end select
    end do
contains
    subroutine print_response(r)
        type(frequency_response), intent(in) :: r !! computed frequency response
        integer :: i
        print '(i0,1x,i0)',r%status,size(r%response)
        do i=1,size(r%response)
            print '(2(es26.17e3,1x))',real(r%response(i),dp),aimag(r%response(i))
        end do
    end subroutine print_response
end program design_driver
