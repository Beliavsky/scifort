! SPDX-License-Identifier: MIT
! Accept independently generated SciPy legacy STFT coefficients on standard input.
! Emit inverse sample values and times for comparison against scipy.signal.istft.
program reference_signal_istft_parameter_windows
    use scifort_kinds, only : dp
    use scifort_signal_istft, only : istft
    implicit none
    complex(dp), allocatable :: z(:,:), ycomplex(:)
    real(dp), allocatable :: yreal(:), times(:), coeff(:), weights(:)
    real(dp) :: rate, a, b, zr, zi
    integer :: nc, c, nf, nt, nw, no, nfft, kind, mode, ncw, code, k, j, i
    character(len=50) :: window
    character(len=16) :: scaling

    read(*,*) nc
    do c=1,nc
        read(*,*) nf,nt,nw,no,nfft,kind,mode,rate,a,b,ncw
        read(*,*) window,scaling
        allocate(coeff(ncw),z(nf,nt))
        do i=1,ncw
            read(*,*) coeff(i)
        end do
        do j=1,nt
            do k=1,nf
                read(*,*) zr,zi
                z(k,j)=cmplx(zr,zi,dp)
            end do
        end do
        ! kind=1 -> real, kind=2 -> complex; mode=1 -> first param,
        ! 2 -> first and second params, 3 -> general cosine coefficients,
        ! 4 -> direct caller-supplied real window, 0 -> no parameters.
        if (kind==1) then
            select case(mode)
            case(0)
                call istft(z,yreal,code,fs=rate,window=trim(window), &
                    nperseg=nw,noverlap=no,nfft=nfft,scaling=trim(scaling),times=times)
            case(1)
                call istft(z,yreal,code,fs=rate,window=trim(window), &
                    window_parameter=a,nperseg=nw,noverlap=no,nfft=nfft, &
                    scaling=trim(scaling),times=times)
            case(2)
                call istft(z,yreal,code,fs=rate,window=trim(window), &
                    window_parameter=a,window_parameter2=b, &
                    nperseg=nw,noverlap=no,nfft=nfft,scaling=trim(scaling),times=times)
            case(3)
                call istft(z,yreal,code,fs=rate,window=trim(window), &
                    window_coefficients=coeff,nperseg=nw,noverlap=no,nfft=nfft, &
                    scaling=trim(scaling),times=times)
            case(4)
                call istft(z,yreal,code,fs=rate,window_values=coeff, &
                    nperseg=nw,noverlap=no,nfft=nfft,scaling=trim(scaling),times=times)
            end select
        else
            select case(mode)
            case(0)
                call istft(z,ycomplex,code,fs=rate,window=trim(window), &
                    nperseg=nw,noverlap=no,nfft=nfft,scaling=trim(scaling),times=times)
            case(1)
                call istft(z,ycomplex,code,fs=rate,window=trim(window), &
                    window_parameter=a,nperseg=nw,noverlap=no,nfft=nfft, &
                    scaling=trim(scaling),times=times)
            case(2)
                call istft(z,ycomplex,code,fs=rate,window=trim(window), &
                    window_parameter=a,window_parameter2=b, &
                    nperseg=nw,noverlap=no,nfft=nfft,scaling=trim(scaling),times=times)
            case(3)
                call istft(z,ycomplex,code,fs=rate,window=trim(window), &
                    window_coefficients=coeff,nperseg=nw,noverlap=no,nfft=nfft, &
                    scaling=trim(scaling),times=times)
            case(4)
                call istft(z,ycomplex,code,fs=rate,window_values=coeff, &
                    nperseg=nw,noverlap=no,nfft=nfft,scaling=trim(scaling),times=times)
            end select
        end if
        if (kind == 1) then
            print '(3(i0,1x))', c,code,size(yreal)
            do i=1,size(yreal)
                print '(3(es26.17e3,1x))', times(i),yreal(i),0.0_dp
            end do
            deallocate(yreal)
        else
            print '(3(i0,1x))', c,code,size(ycomplex)
            do i=1,size(ycomplex)
                print '(3(es26.17e3,1x))', times(i),real(ycomplex(i),dp),aimag(ycomplex(i))
            end do
            deallocate(ycomplex)
        end if
        deallocate(z,coeff,times)
    end do
end program reference_signal_istft_parameter_windows
