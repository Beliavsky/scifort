! SPDX-License-Identifier: MIT
! Text-based SciPy comparison driver for one-dimensional spectral functions.
program spectral_driver
    use scifort_signal
    implicit none
    character(32) :: method, window, trend, scale, average, boundary
    integer :: n, ny, nw, no, nf, flag, i, j, code
    real(dp) :: fs
    real(dp), allocatable :: x(:), y(:), w(:)
    type(spectral_result) :: p
    type(cross_spectral_result) :: c
    type(stft_result) :: t
    type(spectrogram_result) :: s

    read(*,'(a)') method
    read(*,*) n, ny, nw, no, nf, fs, flag
    read(*,'(a)') window
    read(*,'(a)') trend
    read(*,'(a)') scale
    read(*,'(a)') average
    read(*,'(a)') boundary
    allocate(x(n), y(ny))
    do i = 1, n
        read(*,*) x(i)
    end do
    do i = 1, ny
        read(*,*) y(i)
    end do
    select case(trim(method))
    case('get_window')
        w = get_window(trim(window), nw, fftbins=(flag/=0), status=code)
        write(*,*) code, size(w), 1
        do i = 1, size(w)
            write(*,'(2(es26.17e3,1x))') real(i-1,dp),w(i)
        end do
    case('periodogram')
        p = periodogram(x, fs=fs, window=trim(window), nfft=nf, detrend=trim(trend), scaling=trim(scale))
        call print_real(p)
    case('welch')
        p = welch(x, fs=fs, window=trim(window), nperseg=nw, noverlap=no, nfft=nf, &
            detrend=trim(trend), scaling=trim(scale), average=trim(average))
        call print_real(p)
    case('csd')
        c = csd(x,y,fs=fs,window=trim(window),nperseg=nw,noverlap=no,nfft=nf, &
            detrend=trim(trend),scaling=trim(scale),average=trim(average))
        write(*,*) c%status, size(c%frequencies), 1
        do i = 1,size(c%frequencies)
            write(*,'(3(es26.17e3,1x))') c%frequencies(i), real(c%power(i)), aimag(c%power(i))
        end do
    case('coherence')
        p = coherence(x,y,fs=fs,window=trim(window),nperseg=nw,noverlap=no,nfft=nf,detrend=trim(trend))
        call print_real(p)
    case('stft')
        t = stft(x,fs=fs,window=trim(window),nperseg=nw,noverlap=no,nfft=nf, &
            detrend=trim(trend),boundary=trim(boundary),padded=(flag/=0),scaling=trim(scale))
        write(*,*) t%status,size(t%frequencies),size(t%times)
        do j = 1,size(t%times)
            do i = 1,size(t%frequencies)
                write(*,'(4(es26.17e3,1x))') t%frequencies(i),t%times(j), &
                    real(t%values(i,j)),aimag(t%values(i,j))
            end do
        end do
    case('spectrogram')
        s = spectrogram(x,fs=fs,window=trim(window),nperseg=nw,noverlap=no,nfft=nf, &
            detrend=trim(trend),scaling=trim(scale))
        write(*,*) s%status,size(s%frequencies),size(s%times)
        do j = 1,size(s%times)
            do i = 1,size(s%frequencies)
                write(*,'(3(es26.17e3,1x))') s%frequencies(i),s%times(j),s%power(i,j)
            end do
        end do
    end select
contains
    subroutine print_real(a)
        type(spectral_result), intent(in) :: a !! spectral result to serialize
        integer :: k
        write(*,*) a%status,size(a%frequencies),1
        do k = 1,size(a%frequencies)
            write(*,'(2(es26.17e3,1x))') a%frequencies(k),a%power(k)
        end do
    end subroutine print_real
end program spectral_driver
