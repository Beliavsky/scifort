! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Batch numerical comparison driver for FIR design and phase conversion.
program fir_advanced_driver
    use scifort_signal_fir_advanced, only : firwin2, firls, minimum_phase
    use scifort_kinds, only : dp
    implicit none
    character(16) :: mode, wname
    integer :: ios, n, nf, nedge, st, anti, half_output, i
    real(dp) :: fs
    real(dp), allocatable :: freq(:), gain(:), weights(:), coeffs(:), taps(:)
    do
        read(*,*,iostat=ios) mode
        if (ios /= 0) exit
        select case (trim(mode))
        case ('WIN2')
            read(*,*) n, nedge, nf, anti, fs
            read(*,*) wname
            allocate(freq(nedge),gain(nedge))
            read(*,*) freq
            read(*,*) gain
            taps = firwin2(n,freq,gain,nfreqs=nf,window=trim(wname), &
                antisymmetric=anti/=0,fs=fs,status=st)
            deallocate(freq,gain)
        case ('LS')
            read(*,*) n, nedge, fs
            allocate(freq(nedge),gain(nedge),weights(nedge/2))
            read(*,*) freq
            read(*,*) gain
            read(*,*) weights
            taps = firls(n,freq,gain,weight=weights,fs=fs,status=st)
            deallocate(freq,gain,weights)
        case ('MIN')
            read(*,*) n, nf, half_output
            read(*,*) wname
            allocate(coeffs(n))
            read(*,*) coeffs
            taps = minimum_phase(coeffs,method=trim(wname),n_fft=nf, &
                half=half_output/=0,status=st)
            deallocate(coeffs)
        case default
            error stop 'invalid comparison command'
        end select
        write(*,'(i0,1x,i0)') st, size(taps)
        do i = 1, size(taps)
            write(*,'(es26.17e3)') taps(i)
        end do
        if (allocated(taps)) deallocate(taps)
    end do
end program fir_advanced_driver
