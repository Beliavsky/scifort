! SPDX-License-Identifier: MIT
! Reference driver for test/compare_signal_multitaper.py
program signal_multitaper_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_psd
    use scifort_signal_spectral, only : spectral_result
    implicit none
    type(spectral_result) :: result
    real(dp), allocatable :: x(:)
    real(dp) :: nw,fs
    integer :: n,k,nfft,idetrend,iweight,ios
    logical :: centered
    character(16) :: weighting
    do
        read(*,*,iostat=ios) n,nw,k,fs,nfft,idetrend,iweight
        if (ios /= 0) exit
        allocate(x(n))
        read(*,*,iostat=ios) x
        if (ios /= 0) exit
        centered = idetrend /= 0
        weighting = 'equal'
        if (iweight /= 0) weighting = 'eigen'
        result = multitaper_psd(x,nw,kmax=k,fs=fs,nfft=nfft, &
            detrend=centered,weighting=trim(weighting))
        write(*,'(a,1x,i0,1x,i0)') 'CASE',size(result%power),result%status
        if (result%status == 0) then
            write(*,'(*(es25.16,1x))') result%frequencies
            write(*,'(*(es25.16,1x))') result%power
        end if
        deallocate(x)
    end do
end program signal_multitaper_reference_driver
