! SPDX-License-Identifier: MIT
! SciPy comparison driver for complex one-dimensional spectral functions.
program spectral_complex_driver
    use scifort_signal
    implicit none
    character(32) :: method, window, trend, scale, average, boundary
    integer :: n, ny, nw, no, nf, pad, i, j
    real(dp) :: fs, re, im
    complex(dp), allocatable :: x(:), y(:)
    type(spectral_result) :: p
    type(cross_spectral_result) :: c
    type(stft_result) :: z
    type(spectrogram_result) :: s

    read(*,'(a)') method
    read(*,*) n, ny, nw, no, nf, fs, pad
    read(*,'(a)') window
    read(*,'(a)') trend
    read(*,'(a)') scale
    read(*,'(a)') average
    read(*,'(a)') boundary
    allocate(x(n),y(ny))
    do i = 1, n
        read(*,*) re, im
        x(i) = cmplx(re,im,dp)
    end do
    do i = 1, ny
        read(*,*) re, im
        y(i) = cmplx(re,im,dp)
    end do
    select case (trim(method))
    case ('periodogram')
        p = periodogram(x,fs=fs,window=trim(window),nfft=nf, &
            detrend=trim(trend),scaling=trim(scale))
        call print_real(p)
    case ('welch')
        p = welch(x,fs=fs,window=trim(window),nperseg=nw,noverlap=no, &
            nfft=nf,detrend=trim(trend),scaling=trim(scale),average=trim(average))
        call print_real(p)
    case ('csd')
        c = csd(x,y,fs=fs,window=trim(window),nperseg=nw,noverlap=no, &
            nfft=nf,detrend=trim(trend),scaling=trim(scale),average=trim(average))
        write(*,*) c%status, size(c%frequencies), 1
        do i = 1, size(c%frequencies)
            write(*,'(3(es26.17e3,1x))') c%frequencies(i), real(c%power(i),dp), aimag(c%power(i))
        end do
    case ('coherence')
        p = coherence(x,y,fs=fs,window=trim(window),nperseg=nw,noverlap=no, &
            nfft=nf,detrend=trim(trend))
        call print_real(p)
    case ('stft')
        z = stft(x,fs=fs,window=trim(window),nperseg=nw,noverlap=no, &
            nfft=nf,detrend=trim(trend),boundary=trim(boundary), &
            padded=(pad/=0),scaling=trim(scale))
        write(*,*) z%status,size(z%frequencies),size(z%times)
        do j = 1, size(z%times)
            do i = 1, size(z%frequencies)
                write(*,'(4(es26.17e3,1x))') z%frequencies(i),z%times(j), &
                    real(z%values(i,j),dp), aimag(z%values(i,j))
            end do
        end do
    case ('spectrogram')
        s = spectrogram(x,fs=fs,window=trim(window),nperseg=nw, &
            noverlap=no,nfft=nf,detrend=trim(trend),scaling=trim(scale))
        write(*,*) s%status,size(s%frequencies),size(s%times)
        do j = 1, size(s%times)
            do i = 1, size(s%frequencies)
                write(*,'(3(es26.17e3,1x))') s%frequencies(i),s%times(j),s%power(i,j)
            end do
        end do
    end select
contains
    subroutine print_real(r)
        type(spectral_result), intent(in) :: r !! result to serialize
        integer :: k
        write(*,*) r%status,size(r%frequencies),1
        do k = 1, size(r%power)
            write(*,'(2(es26.17e3,1x))') r%frequencies(k),r%power(k)
        end do
    end subroutine print_real
end program spectral_complex_driver
