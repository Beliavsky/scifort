! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program signal_multitaper_harmonic_scan_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_harmonic_scan, multitaper_harmonic_scan_result
    implicit none
    integer :: cases, q, n, k, nf, isc, center, method_id, edges, nsel, i
    integer, allocatable :: bins(:)
    real(dp) :: nw, rate, alpha
    real(dp), allocatable :: re(:), im(:)
    type(multitaper_harmonic_scan_result) :: result
    character(12) :: method

    read(*,*) cases
    do q=1,cases
        read(*,*) n,k,nf,isc,center,method_id,edges,nsel,nw,rate,alpha
        allocate(re(n),im(n))
        do i=1,n
            read(*,*) re(i),im(i)
        end do
        select case(method_id)
        case(0)
            method='none'
        case(1)
            method='bonferroni'
        case(2)
            method='bh'
        case(3)
            method='by'
        case(4)
            method='holm'
        end select
        if (nsel >= 0) then
            allocate(bins(nsel))
            if (nsel > 0) read(*,*) bins
            if (isc == 0) then
                result = multitaper_harmonic_scan(re,nw,kmax=k,nfft=nf,fs=rate, &
                    detrend=center /= 0,method=trim(method),alpha=alpha,test_bins=bins)
            else
                result = multitaper_harmonic_scan(cmplx(re,im,dp),nw,kmax=k,nfft=nf,fs=rate, &
                    detrend=center /= 0,method=trim(method),alpha=alpha,test_bins=bins)
            end if
            deallocate(bins)
        else
            if (isc == 0) then
                result = multitaper_harmonic_scan(re,nw,kmax=k,nfft=nf,fs=rate, &
                    detrend=center /= 0,method=trim(method),alpha=alpha,include_edges=edges /= 0)
            else
                result = multitaper_harmonic_scan(cmplx(re,im,dp),nw,kmax=k,nfft=nf,fs=rate, &
                    detrend=center /= 0,method=trim(method),alpha=alpha,include_edges=edges /= 0)
            end if
        end if
        write(*,'(i0,1x,i0,1x,i0)') result%status,size(result%pvalue),result%num_tests
        do i=1,size(result%pvalue)
            write(*,'(6(es24.16e3,1x),i0)') result%frequencies(i),result%pvalue(i), &
                result%adjusted_pvalue(i), result%f_statistic(i), &
                real(result%amplitude(i),dp),aimag(result%amplitude(i)), &
                merge(1,0,result%significant(i))
        end do
        deallocate(re,im)
    end do
end program signal_multitaper_harmonic_scan_reference_driver
