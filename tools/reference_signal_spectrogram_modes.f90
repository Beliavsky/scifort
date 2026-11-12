! SPDX-License-Identifier: MIT
! Stream seeded SciPy comparison cases through the independent Fortran spectra.
program reference_signal_spectrogram_modes
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : spectrogram_real => spectrogram, spectrogram_result
    use scifort_signal_spectral_complex, only : spectrogram_complex
    implicit none

    integer :: count, case_id, kind_id, mode_id, window_id, n, nw, no, nf, onesided_id, trend_id, scale_id
    integer :: i, j
    integer, allocatable :: nperseg, noverlap, nfft
    real(dp) :: fs
    real(dp), allocatable :: xr(:), xi(:), win(:)
    complex(dp), allocatable :: xc(:)
    character(len=9), parameter :: modes(5) = [character(len=9) :: &
        'psd', 'complex', 'magnitude', 'angle', 'phase']
    character(len=8), parameter :: trends(3) = [character(len=8) :: &
        'constant', 'linear', 'none']
    character(len=8), parameter :: scales(2) = [character(len=8) :: &
        'density', 'spectrum']
    type(spectrogram_result) :: result

    read(*,*) count
    do case_id = 1, count
        read(*,*) kind_id, mode_id, window_id, n, nw, no, nf, onesided_id, trend_id, scale_id, fs
        allocate(xr(n),xi(n),xc(n))
        read(*,*) xr
        read(*,*) xi
        if (nw > 0) nperseg = nw
        if (no >= 0) noverlap = no
        if (nf > 0) nfft = nf
        xc = cmplx(xr,xi,dp)
        if (window_id == 3) then
            allocate(win(nw))
            read(*,*) win
        end if
        if (kind_id == 0) then
            select case (window_id)
            case (0)
                result = spectrogram_real(xr,fs=fs,nperseg=nperseg,noverlap=noverlap,nfft=nfft, &
                    detrend=trim(trends(trend_id)),scaling=trim(scales(scale_id)), &
                    return_onesided=onesided_id==1,mode=trim(modes(mode_id)))
            case (1)
                result = spectrogram_real(xr,fs=fs,nperseg=nperseg,noverlap=noverlap,nfft=nfft, &
                    detrend=trim(trends(trend_id)),scaling=trim(scales(scale_id)), &
                    return_onesided=onesided_id==1,mode=trim(modes(mode_id)),window='hann')
            case (2)
                result = spectrogram_real(xr,fs=fs,nperseg=nperseg,noverlap=noverlap,nfft=nfft, &
                    detrend=trim(trends(trend_id)),scaling=trim(scales(scale_id)), &
                    return_onesided=onesided_id==1,mode=trim(modes(mode_id)),window='gaussian', &
                    window_parameter=4.3_dp)
            case (3)
                result = spectrogram_real(xr,fs=fs,nperseg=nperseg,noverlap=noverlap,nfft=nfft, &
                    detrend=trim(trends(trend_id)),scaling=trim(scales(scale_id)), &
                    return_onesided=onesided_id==1,mode=trim(modes(mode_id)),window_values=win)
            end select
        else
            select case (window_id)
            case (0)
                result = spectrogram_complex(xc,fs=fs,nperseg=nperseg,noverlap=noverlap,nfft=nfft, &
                    detrend=trim(trends(trend_id)),scaling=trim(scales(scale_id)), &
                    mode=trim(modes(mode_id)))
            case (1)
                result = spectrogram_complex(xc,fs=fs,nperseg=nperseg,noverlap=noverlap,nfft=nfft, &
                    detrend=trim(trends(trend_id)),scaling=trim(scales(scale_id)), &
                    mode=trim(modes(mode_id)),window='hann')
            case (2)
                result = spectrogram_complex(xc,fs=fs,nperseg=nperseg,noverlap=noverlap,nfft=nfft, &
                    detrend=trim(trends(trend_id)),scaling=trim(scales(scale_id)), &
                    mode=trim(modes(mode_id)),window='gaussian',window_parameter=4.3_dp)
            case (3)
                result = spectrogram_complex(xc,fs=fs,nperseg=nperseg,noverlap=noverlap,nfft=nfft, &
                    detrend=trim(trends(trend_id)),scaling=trim(scales(scale_id)), &
                    mode=trim(modes(mode_id)),window_values=win)
            end select
        end if
        write(*,'(i0,1x,i0,1x,i0)') result%status, size(result%frequencies), size(result%times)
        do i = 1, size(result%frequencies)
            write(*,'(es26.17e3)') result%frequencies(i)
        end do
        do j = 1, size(result%times)
            write(*,'(es26.17e3)') result%times(j)
        end do
        do j = 1, size(result%times)
            do i = 1, size(result%frequencies)
                if (mode_id == 2) then
                    write(*,'(es26.17e3,1x,es26.17e3)') real(result%values(i,j),dp), aimag(result%values(i,j))
                else
                    write(*,'(es26.17e3,1x,es26.17e3)') result%power(i,j), 0.0_dp
                end if
            end do
        end do
        deallocate(xr,xi,xc)
        if (allocated(nperseg)) deallocate(nperseg)
        if (allocated(noverlap)) deallocate(noverlap)
        if (allocated(nfft)) deallocate(nfft)
        if (allocated(win)) deallocate(win)
    end do
end program reference_signal_spectrogram_modes
