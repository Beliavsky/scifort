! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program signal_multitaper_noise_floor_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_noise_floor, multitaper_noise_floor_result
    implicit none
    integer :: cases, q, n, k, nf, complex_flag, trend_flag, weight_flag, width, i
    real(dp) :: nw, rate
    real(dp), allocatable :: re(:), im(:)
    character(10) :: trend, weighting
    type(multitaper_noise_floor_result) :: out

    read(*,*) cases
    do q = 1, cases
        read(*,*) n,k,nf,complex_flag,trend_flag,weight_flag,width,nw,rate
        allocate(re(n),im(n))
        do i = 1, n
            read(*,*) re(i),im(i)
        end do
        select case (trend_flag)
        case (0)
            trend = 'none'
        case (1)
            trend = 'constant'
        case (2)
            trend = 'linear'
        end select
        select case (weight_flag)
        case (0)
            weighting = 'equal'
        case (1)
            weighting = 'eigen'
        end select
        if (complex_flag == 1) then
            out = multitaper_noise_floor(cmplx(re,im,dp),nw,width,kmax=k,nfft=nf,fs=rate, &
                detrend_mode=trim(trend),weighting=trim(weighting))
        else
            out = multitaper_noise_floor(re,nw,width,kmax=k,nfft=nf,fs=rate, &
                detrend_mode=trim(trend),weighting=trim(weighting))
        end if
        write(*,'(i0,1x,i0)') out%status,size(out%power)
        do i = 1, size(out%power)
            write(*,'(4(es24.16e3,1x))') out%frequencies(i),out%power(i), &
                out%background(i),out%whitened(i)
        end do
        deallocate(re,im)
    end do
end program signal_multitaper_noise_floor_reference_driver
