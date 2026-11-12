! SPDX-License-Identifier: MIT
! Two-sided legacy scipy.signal spectral estimators and assumed-rank wrappers.
program test_signal_twosided
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : periodogram, welch, csd, stft, &
        spectrogram, spectral_result, cross_spectral_result, stft_result, &
        spectrogram_result, spectral_status_success, spectral_status_invalid_input
    use scifort_signal_spectral_nd, only : periodogram_into, welch_into, csd_into
    use scifort_signal_timefreq_nd, only : stft_into, spectrogram_into
    implicit none
    real(dp), parameter :: tol=3.0e-12_dp
    real(dp) :: x(2,32), p(2,8), invalid(2,7)
    real(dp) :: powers(2,8,3)
    real(dp), allocatable :: f(:), t(:)
    complex(dp) :: z(2,16,3), zshort(2,8,3)
    complex(dp) :: cross_arr(2,8)
    type(spectral_result) :: narrow, full, reference
    type(cross_spectral_result) :: c
    type(stft_result) :: dft
    type(spectrogram_result) :: psd
    integer :: i, st, k

    do i=1,32
        x(1,i)=0.02_dp*real(i,dp)+sin(real(i,dp)*0.38_dp)+0.2_dp*cos(real(i*i,dp)*0.09_dp)
        x(2,i)=0.2_dp*sin(real(i,dp)*0.71_dp)+cos(real(i,dp)*0.15_dp)
    end do

    narrow = periodogram(x(1,1:8), nfft=8, detrend='none')
    full = periodogram(x(1,1:8), nfft=8, detrend='none',return_onesided=.false.)
    call check(narrow%status==spectral_status_success .and. full%status==spectral_status_success, &
        'periodogram succeeds')
    call check(size(narrow%frequencies)==5 .and. size(full%frequencies)==8, 'one-/two-sided frequency sizes')
    call check(abs(full%frequencies(5)+0.5_dp)<tol, 'negative Nyquist frequency (fs default 1)')
    call check(abs(full%frequencies(8)+0.125_dp)<tol, 'last negative frequency')
    call check(abs(narrow%power(2)-2.0_dp*full%power(2))<tol, 'interior one-sided PSD doubled')
    call check(abs(narrow%power(5)-full%power(5))<tol, 'even Nyquist bin is not doubled')
    call check(abs(narrow%power(1)-full%power(1))<tol, 'DC bin is not doubled')
    call check(abs(full%power(8)-full%power(2))<tol, 'real PSD positive/negative symmetry')

    full=welch(x(1,:),window='gaussian',window_parameter=2.0_dp,nperseg=8, &
               noverlap=3,nfft=9,average='median',return_onesided=.false.)
    call check(full%status==spectral_status_success, 'two-sided median Welch')
    call check(size(full%power)==9 .and. abs(full%frequencies(6)+4.0_dp/9.0_dp)<tol, &
        'odd-length Welch frequencies')

    c=csd(x(1,:),x(2,:),nperseg=8,noverlap=3,nfft=8,return_onesided=.false.)
    call check(c%status==spectral_status_success, 'real two-sided CSD')
    call check(size(c%power)==8,'CSD bin count')
    call check(maxval(abs(c%power(2:4)-conjg(c%power(8:6:-1))))<tol,'CSD Hermitian symmetry')

    ! Assume rank 2 with frequency on the last axis; exact output shape is the
    ! responsibility of the caller and is checked before any writes.
    p=-999.0_dp
    call periodogram_into(x,p,f,nfft=8,detrend='none',axis=2,return_onesided=.false.,status=st)
    call check(st==spectral_status_success,'rank-two periodogram into')
    do k=1,2
        reference=periodogram(x(k,:),nfft=8,detrend='none',return_onesided=.false.)
        call check(maxval(abs(reference%power-p(k,:)))<tol, 'rank-two periodogram lane')
    end do
    call check(maxval(abs(f- reference%frequencies))<tol,'periodogram into frequency grid')

    call welch_into(x,p,f,nperseg=8,noverlap=3,nfft=8,axis=2,return_onesided=.false.,status=st)
    call check(st==spectral_status_success,'rank-two Welch into')
    do k=1,2
        reference=welch(x(k,:),nperseg=8,noverlap=3,nfft=8,return_onesided=.false.)
        call check(maxval(abs(reference%power-p(k,:)))<tol,'rank-two Welch lane')
    end do
    call csd_into(x,x, cross_arr,f,nperseg=8,noverlap=3,nfft=8, &
                  axis=2,return_onesided=.false.,status=st)
    call check(st==spectral_status_success,'rank-two CSD into')
    do k=1,2
        c=csd(x(k,:),x(k,:),nperseg=8,noverlap=3,nfft=8,return_onesided=.false.)
        call check(maxval(abs(c%power-cross_arr(k,:)))<tol,'rank-two CSD lane')
    end do

    invalid=-123.0_dp
    call periodogram_into(x,invalid,f,nfft=8,axis=2,return_onesided=.false.,status=st)
    call check(st==spectral_status_invalid_input,'invalid two-sided output shape detected')
    call check(all(invalid==-123.0_dp),'invalid output shape leaves caller values untouched')
    call check(size(f)==0,'invalid output clears frequency coordinates')

    ! Noncontiguous output view: odd frequency slots are untouched by the call.
    z=cmplx(-17.0_dp,9.0_dp,dp)
    call stft_into(x(:,1:16),z(:,1:16:2,:),f,t,window='hann',nperseg=8, &
                  noverlap=4,nfft=8,axis=2,boundary='none',padded=.false., &
                  return_onesided=.false.,status=st)
    call check(st==spectral_status_success,'strided two-sided STFT into')
    call check(all(z(:,2:16:2,:)==cmplx(-17.0_dp,9.0_dp,dp)), &
               'unchanged elements of strided STFT output')
    do k=1,2
        dft=stft(x(k,1:16),window='hann',nperseg=8,noverlap=4,nfft=8, &
                 boundary='none',padded=.false.,return_onesided=.false.)
        call check(maxval(abs(dft%values-z(k,1:16:2,:)))<tol,'strided STFT lane')
    end do
    call check(abs(f(5)+0.5_dp)<tol, 'strided STFT negative Nyquist')

    powers=-111.0_dp
    call spectrogram_into(x(:,1:16),powers,f,t,window='hann',nperseg=8, &
                          noverlap=4,nfft=8,axis=2,mode='psd', &
                          return_onesided=.false.,status=st)
    call check(st==spectral_status_success,'two-sided rank-two PSD spectrogram')
    do k=1,2
        psd=spectrogram(x(k,1:16),window='hann',nperseg=8,noverlap=4,nfft=8, &
            return_onesided=.false.)
        call check(maxval(abs(psd%power-powers(k,:,:)))<tol,'rank-two PSD spectrogram lane')
    end do
    zshort=cmplx(-10.0_dp,10.0_dp,dp)
    call spectrogram_into(x(:,1:16),zshort,f,t,window='hann',nperseg=8, &
                          noverlap=4,nfft=8,axis=2,mode='complex', &
                          return_onesided=.false.,status=st)
    call check(st==spectral_status_success,'complex-mode two-sided spectrogram')
    call check(size(f)==8 .and. abs(f(5)+0.5_dp)<tol,'complex spectrogram negative Nyquist')
    call check(all(abs(zshort)<huge(1.0_dp)),'complex spectrogram finite')

    call high_rank_check()
    print '(a)', 'PASS: two-sided legacy spectral regressions; 2D, strided, rank-14 and rank-15 paths'
contains
    subroutine check(success,message)
        logical, intent(in) :: success !! condition to test
        character(*), intent(in) :: message !! assertion description
        if (.not. success) then
            print '(2a)', 'FAIL: ',message
            error stop 1
        end if
    end subroutine check

    subroutine high_rank_check()
        ! Rank-14 input STFT becomes rank-15; spectral ranks remain unchanged.
        real(dp) :: high(1,1,1,1,1,1,1,1,1,1,1,1,1,16)
        complex(dp) :: zout(1,1,1,1,1,1,1,1,1,1,1,1,1,8,3)
        real(dp) :: high15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,16)
        real(dp) :: pout15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,8)
        real(dp), allocatable :: fh(:), th(:)
        type(stft_result) :: comparison
        type(spectral_result) :: p_ref
        integer :: ierr
        high(1,1,1,1,1,1,1,1,1,1,1,1,1,:)=x(1,1:16)
        zout=cmplx(-33.0_dp,0.0_dp,dp)
        call stft_into(high,zout,fh,th,nperseg=8,noverlap=4,nfft=8, &
            boundary='none',padded=.false.,axis=14,return_onesided=.false.,status=ierr)
        call check(ierr==spectral_status_success,'rank-14 input, rank-15 output STFT')
        comparison=stft(x(1,1:16),nperseg=8,noverlap=4,nfft=8, &
            boundary='none',padded=.false.,return_onesided=.false.)
        call check(maxval(abs(zout(1,1,1,1,1,1,1,1,1,1,1,1,1,:,:)- &
            comparison%values))<tol,'rank-14 STFT numerical values')
        high15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,:)=x(1,1:16)
        pout15=-99.0_dp
        call periodogram_into(high15,pout15,fh,nfft=8,axis=15, &
            return_onesided=.false.,status=ierr)
        call check(ierr==spectral_status_success,'rank-15 periodogram')
        p_ref=periodogram(x(1,1:16),nfft=8,return_onesided=.false.)
        call check(maxval(abs(pout15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,:)- &
            p_ref%power))<tol,'rank-15 periodogram values')
    end subroutine high_rank_check
end program test_signal_twosided
