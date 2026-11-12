! SPDX-License-Identifier: MIT
! Read randomized multichannel test cases for the independent SciPy oracle.
program signal_multitaper_channels_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_spectral_matrix, multitaper_matrix_result
    implicit none
    integer :: n, nc, nf, nt, kind_case, ax, center_case, eigen_case, &
        t, c, a, b, bin, stat
    real(dp) :: nw, fs, re, im
    real(dp), allocatable :: samples_real(:,:), output_values(:)
    complex(dp), allocatable :: samples_complex(:,:)
    type(multitaper_matrix_result) :: r
    character(len=5) :: weighting

    do
        read(*,*,iostat=stat) n,nc,nf,nt,nw,fs,kind_case,ax,center_case,eigen_case
        if (stat /= 0) exit
        allocate(samples_complex(n,nc),samples_real(n,nc))
        do c = 1, nc
            do t = 1, n
                read(*,*) re,im
                samples_complex(t,c) = cmplx(re,im,dp)
                samples_real(t,c) = re
            end do
        end do
        weighting = 'equal'
        if (eigen_case /= 0) weighting = 'eigen'
        if (kind_case == 0) then
            if (ax == 1) then
                r = multitaper_spectral_matrix(samples_real,nw,kmax=nt,fs=fs, &
                    nfft=nf,detrend=center_case/=0,weighting=trim(weighting),axis=ax)
            else
                r = multitaper_spectral_matrix(transpose(samples_real),nw,kmax=nt, &
                    fs=fs,nfft=nf,detrend=center_case/=0,weighting=trim(weighting),axis=ax)
            end if
        else
            if (ax == 1) then
                r = multitaper_spectral_matrix(samples_complex,nw,kmax=nt,fs=fs, &
                    nfft=nf,detrend=center_case/=0,weighting=trim(weighting),axis=ax)
            else
                r = multitaper_spectral_matrix(transpose(samples_complex),nw,kmax=nt, &
                    fs=fs,nfft=nf,detrend=center_case/=0,weighting=trim(weighting),axis=ax)
            end if
        end if
        print '(a,1x,i0,1x,i0,1x,i0)', 'CASE',size(r%frequencies),nc,r%status
        write(*,'(*(es25.16,1x))') r%frequencies
        allocate(output_values(2*size(r%frequencies)))
        do b = 1, nc
            do a = 1, nc
                do bin = 1, size(r%frequencies)
                    output_values(2*bin-1) = real(r%power(bin,a,b),dp)
                    output_values(2*bin) = aimag(r%power(bin,a,b))
                end do
                write(*,'(*(es25.16,1x))') output_values
            end do
        end do
        deallocate(samples_complex,samples_real,output_values)
    end do
end program signal_multitaper_channels_reference_driver
