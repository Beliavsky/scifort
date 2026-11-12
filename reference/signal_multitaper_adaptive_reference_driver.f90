! SPDX-License-Identifier: MIT
! Read independent randomized real/complex adaptive spectral test cases.
program signal_multitaper_adaptive_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_psd_ci, multitaper_spectral_matrix, &
        multitaper_confidence_result, multitaper_matrix_result
    implicit none
    integer :: n, nc, nf, nt, kind_case, ax, center_case, t, c, a, b, j, status
    real(dp) :: nw, fs, conf, re, im
    real(dp), allocatable :: xr(:,:), vec(:)
    complex(dp), allocatable :: xc(:,:)
    type(multitaper_matrix_result) :: sm
    type(multitaper_confidence_result) :: ci

    do
        read(*,*,iostat=status) n,nc,nf,nt,nw,fs,kind_case,ax,center_case,conf
        if (status /= 0) exit
        allocate(xr(n,nc),xc(n,nc))
        do c = 1, nc
            do t = 1, n
                read(*,*) re,im
                xr(t,c) = re
                xc(t,c) = cmplx(re,im,dp)
            end do
        end do
        if (kind_case == 0) then
            if (ax == 1) then
                sm = multitaper_spectral_matrix(xr,nw,kmax=nt,fs=fs,nfft=nf, &
                    detrend=center_case/=0,weighting='adaptive',axis=1)
            else
                sm = multitaper_spectral_matrix(transpose(xr),nw,kmax=nt,fs=fs, &
                    nfft=nf,detrend=center_case/=0,weighting='adaptive',axis=2)
            end if
            ci = multitaper_psd_ci(xr(:,1),nw,kmax=nt,fs=fs,nfft=nf, &
                detrend=center_case/=0,weighting='adaptive',confidence=conf)
        else
            if (ax == 1) then
                sm = multitaper_spectral_matrix(xc,nw,kmax=nt,fs=fs,nfft=nf, &
                    detrend=center_case/=0,weighting='adaptive',axis=1)
            else
                sm = multitaper_spectral_matrix(transpose(xc),nw,kmax=nt,fs=fs, &
                    nfft=nf,detrend=center_case/=0,weighting='adaptive',axis=2)
            end if
            ci = multitaper_psd_ci(xc(:,1),nw,kmax=nt,fs=fs,nfft=nf, &
                detrend=center_case/=0,weighting='adaptive',confidence=conf)
        end if
        print '(a,1x,i0,1x,i0,1x,i0,1x,i0)', &
            'CASE',size(sm%frequencies),nc,sm%status,ci%status
        write(*,'(*(es26.16e3,1x))') sm%frequencies
        allocate(vec(2*size(sm%frequencies)))
        do b = 1, nc
            do a = 1, nc
                do j = 1, size(sm%frequencies)
                    vec(2*j-1) = real(sm%power(j,a,b),dp)
                    vec(2*j) = aimag(sm%power(j,a,b))
                end do
                write(*,'(*(es26.16e3,1x))') vec
            end do
        end do
        write(*,'(*(es26.16e3,1x))') ci%power
        write(*,'(*(es26.16e3,1x))') ci%lower
        write(*,'(*(es26.16e3,1x))') ci%upper
        write(*,'(*(es26.16e3,1x))') ci%stderr_log
        write(*,'(*(es26.16e3,1x))') ci%equivalent_dof
        deallocate(xr,xc,vec)
    end do
end program signal_multitaper_adaptive_reference_driver
