! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Deterministic numeric output for independent comparisons against SciPy.
! Build and run with test/compare_signal_stft_boundaries_scipy.py.
program driver_signal_stft_boundaries
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : stft, stft_result, spectral_status_success
    use scifort_signal_spectral_complex, only : stft_complex
    implicit none
    real(dp), parameter :: fs = 5.5_dp
    integer, parameter :: lengths(7) = [1,2,3,5,8,17,31]
    integer, parameter :: segments(6) = [1,2,3,4,7,9]
    character(8), parameter :: modes(5) = [character(8) :: &
        'zeros', 'even', 'odd', 'constant', 'none']
    character(8) :: mode, scaling, detrend, window
    real(dp), allocatable :: x(:)
    complex(dp), allocatable :: z(:)
    type(stft_result) :: res
    integer :: il, iw, it, im, ip, iq, n, nw, no, nf, j, k, i
    logical :: padded

    do il=1,size(lengths)
        n = lengths(il)
        allocate(x(n),z(n))
        do i=1,n
            x(i) = 0.8_dp*cos(0.43_dp*real(i,dp)) + 0.07_dp*real(i,dp) &
                - 0.02_dp*sin(real(i*i,dp)*0.31_dp)
            z(i) = cmplx(x(i), 0.23_dp*sin(0.21_dp*real(i,dp)) &
                - 0.03_dp*real(i,dp), dp)
        end do
        do iw=1,size(segments)
            nw=segments(iw)
            if (nw > n) cycle
            do ip=0,1
                padded=(ip==1)
                do iq=0,1
                    no = nw/2
                    if (iq==1) no=max(0,nw-1)
                    nf = nw + 2*iq + 1
                    detrend = 'none'
                    if (mod(iw+il+iq,3)==1) detrend = 'constant'
                    if (mod(iw+il+iq,3)==2) detrend = 'linear'
                    window = 'boxcar'
                    if (mod(iw+il,2)==0 .and. nw >= 2) window = 'hann'
                    scaling = 'spectrum'
                    if (iq==1) scaling = 'psd'
                    do im=1,size(modes)
                        mode=modes(im)
                        do it=1,2
                            if (it==1) then
                                res=stft(x,fs=fs,window=trim(window),nperseg=nw, &
                                    noverlap=no,nfft=nf,detrend=trim(detrend), &
                                    boundary=trim(mode),padded=padded,scaling=trim(scaling))
                            else
                                res=stft_complex(z,fs=fs,window=trim(window),nperseg=nw, &
                                    noverlap=no,nfft=nf,detrend=trim(detrend), &
                                    boundary=trim(mode),padded=padded,scaling=trim(scaling))
                            end if
                            if (res%status /= spectral_status_success) &
                                error stop 'unexpected boundary comparison status'
                            write(*,'(a,1x,8(i0,1x),5(a,1x))') 'CASE', &
                                it,n,nw,no,nf,ip,size(res%frequencies),size(res%times), &
                                trim(mode),trim(scaling),trim(detrend),trim(window),'END'
                            write(*,'(*(es26.17e3,1x))') res%frequencies
                            write(*,'(*(es26.17e3,1x))') res%times
                            do j=1,size(res%times)
                                do k=1,size(res%frequencies)
                                    write(*,'(2(es26.17e3,1x))',advance='no') &
                                        real(res%values(k,j),dp),aimag(res%values(k,j))
                                end do
                                write(*,*)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        deallocate(x,z)
    end do
end program driver_signal_stft_boundaries
