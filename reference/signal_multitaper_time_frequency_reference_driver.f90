! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Driver for independent SciPy DPSS and NumPy FFT comparisons.
program signal_multitaper_time_frequency_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_spectrogram, &
        multitaper_spectrogram_channels_result, multitaper_harmonic_ftest, multitaper_harmonic_result
    implicit none
    type(multitaper_spectrogram_channels_result) :: spectrogram
    type(multitaper_harmonic_result) :: ftest
    real(dp), allocatable :: xr(:,:), realx(:)
    complex(dp), allocatable :: xc(:,:), x(:)
    real(dp) :: nw, rate, re, im
    integer :: nspec, ntest, trial, n, nc, seg, overlap, k, nf, kind_code, centre, axis, i, j, f, t, c
    character(16) :: weighting

    read(*,*) nspec
    do trial = 1, nspec
        read(*,*) n,nc,seg,overlap,k,nf,kind_code,centre,axis,nw,rate,weighting
        allocate(xc(n,nc),xr(n,nc))
        do c = 1, nc
            do i = 1, n
                read(*,*) re, im
                xr(i,c) = re
                xc(i,c) = cmplx(re,im,dp)
            end do
        end do
        if (kind_code == 0) then
            if (axis == 1) then
                spectrogram = multitaper_spectrogram(xr,nw,seg,overlap,k,rate,nf,centre == 1,trim(weighting),axis=1)
            else
                spectrogram = multitaper_spectrogram(transpose(xr),nw,seg,overlap,k,rate,nf,centre == 1, &
                    trim(weighting),axis=2)
            end if
        else
            if (axis == 1) then
                spectrogram = multitaper_spectrogram(xc,nw,seg,overlap,k,rate,nf,centre == 1,trim(weighting),axis=1)
            else
                spectrogram = multitaper_spectrogram(transpose(xc),nw,seg,overlap,k,rate,nf,centre == 1, &
                    trim(weighting),axis=2)
            end if
        end if
        write(*,'(4(i0,1x))') spectrogram%status,size(spectrogram%frequencies),size(spectrogram%times), &
            size(spectrogram%power,3)
        do f = 1, size(spectrogram%frequencies)
            write(*,'(es26.17)') spectrogram%frequencies(f)
        end do
        do t = 1, size(spectrogram%times)
            write(*,'(es26.17)') spectrogram%times(t)
        end do
        do c = 1, size(spectrogram%power,3)
            do t = 1, size(spectrogram%times)
                do f = 1, size(spectrogram%frequencies)
                    write(*,'(es26.17)') spectrogram%power(f,t,c)
                end do
            end do
        end do
        deallocate(xc,xr)
    end do
    read(*,*) ntest
    do trial = 1, ntest
        read(*,*) n,k,nf,kind_code,centre,nw,rate
        allocate(x(n),realx(n))
        do i = 1, n
            read(*,*) re, im
            x(i) = cmplx(re,im,dp)
            realx(i) = re
        end do
        if (kind_code == 0) then
            ftest = multitaper_harmonic_ftest(realx,nw,k,rate,nf,centre == 1)
        else
            ftest = multitaper_harmonic_ftest(x,nw,k,rate,nf,centre == 1)
        end if
        write(*,'(2(i0,1x))') ftest%status,size(ftest%frequencies)
        do j = 1, size(ftest%frequencies)
            write(*,'(5(es26.17,1x))') ftest%frequencies(j),ftest%f_statistic(j),ftest%pvalue(j), &
                real(ftest%amplitude(j),dp),aimag(ftest%amplitude(j))
        end do
        deallocate(x,realx)
    end do
end program signal_multitaper_time_frequency_reference_driver
