! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Input/output driver for SciPy multitaper windowed noise-floor comparisons.
program signal_multitaper_noise_floor_window_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_noise_floor_spectrogram, &
        multitaper_noise_floor_spectrogram_channels_result
    implicit none
    integer :: cases, q, n, nc, seg, overlap, k, nf, complex_flag
    integer :: trend_flag, weight_flag, width, axis, i, c, t, f
    real(dp) :: nw, rate, re, im
    real(dp), allocatable :: x(:,:)
    complex(dp), allocatable :: z(:,:)
    character(8) :: trend, weighting
    type(multitaper_noise_floor_spectrogram_channels_result) :: out

    read(*,*) cases
    do q = 1, cases
        read(*,*) n,nc,seg,overlap,k,nf,complex_flag,trend_flag,weight_flag,width,axis,nw,rate
        allocate(x(n,nc),z(n,nc))
        do c = 1, nc
            do i = 1, n
                read(*,*) re, im
                x(i,c) = re
                z(i,c) = cmplx(re,im,dp)
            end do
        end do
        select case(trend_flag)
        case (0)
            trend = 'none'
        case (1)
            trend = 'constant'
        case default
            trend = 'linear'
        end select
        select case(weight_flag)
        case (0)
            weighting = 'equal'
        case (1)
            weighting = 'eigen'
        case default
            weighting = 'adaptive'
        end select
        if (complex_flag == 1) then
            if (axis == 1) then
                out = multitaper_noise_floor_spectrogram(z,nw,width,seg,noverlap=overlap, &
                    kmax=k,fs=rate,nfft=nf,detrend_mode=trim(trend),weighting=trim(weighting))
            else
                out = multitaper_noise_floor_spectrogram(transpose(z),nw,width,seg,noverlap=overlap, &
                    kmax=k,fs=rate,nfft=nf,detrend_mode=trim(trend),weighting=trim(weighting),axis=2)
            end if
        else
            if (axis == 1) then
                out = multitaper_noise_floor_spectrogram(x,nw,width,seg,noverlap=overlap, &
                    kmax=k,fs=rate,nfft=nf,detrend_mode=trim(trend),weighting=trim(weighting))
            else
                out = multitaper_noise_floor_spectrogram(transpose(x),nw,width,seg,noverlap=overlap, &
                    kmax=k,fs=rate,nfft=nf,detrend_mode=trim(trend),weighting=trim(weighting),axis=2)
            end if
        end if
        write(*,'(4(i0,1x))') out%status,size(out%power,1),size(out%power,2),size(out%power,3)
        do f = 1, size(out%frequencies)
            write(*,'(es25.16e3)') out%frequencies(f)
        end do
        do t = 1, size(out%times)
            write(*,'(es25.16e3)') out%times(t)
        end do
        do c = 1, size(out%power,3)
            do t = 1, size(out%power,2)
                do f = 1, size(out%power,1)
                    write(*,'(3(es25.16e3,1x))') out%power(f,t,c), &
                        out%background(f,t,c),out%whitened(f,t,c)
                end do
            end do
        end do
        deallocate(x,z)
    end do
end program signal_multitaper_noise_floor_window_reference_driver
