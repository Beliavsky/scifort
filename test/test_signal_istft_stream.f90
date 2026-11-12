! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_istft_stream
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_signal, only : dp, istft, istft_stream, istft_stream_init, &
        istft_stream_push, istft_stream_finish, istft_status_success, &
        istft_status_invalid_input, istft_status_noninvertible
    implicit none
    integer :: nw, nf, no, frames, widx, bidx, sidx, sides, fidx, k, cases
    integer :: count, nsum
    character(8) :: win, scaling
    logical :: boundary, onesided
    complex(dp) :: bad(5)
    type(istft_stream) :: state
    real(dp), allocatable :: real_chunk(:), tm(:)
    complex(dp), allocatable :: complex_chunk(:)
    integer :: status

    cases = 0
    do nw = 5,8
        do fidx = 0,1
            nf = nw+2*fidx
            do no = 0,nw-1,max(1,nw/2)
                do frames = 1,5,2
                    do sides = 0,1
                        onesided = sides==1
                        do bidx = 0,1
                            boundary = bidx==1
                            do sidx = 0,1
                                scaling = 'spectrum'
                                if (sidx==1) scaling = 'psd'
                                do widx = 0,1
                                    win = 'hann'
                                    if (widx==1) win = 'boxcar'
                                    call check_case(nw,nf,no,frames,onesided,boundary,trim(scaling),trim(win))
                                    cases = cases+1
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
    call istft_stream_init(state,0,status)
    if (status /= istft_status_invalid_input) error stop 'zero frequencies'
    call istft_stream_init(state,5,status,nperseg=8,noverlap=8)
    if (status /= istft_status_invalid_input) error stop 'invalid overlap'
    call istft_stream_init(state,5,status,nperseg=8,nfft=7)
    if (status /= istft_status_invalid_input) error stop 'invalid fft'
    call istft_stream_init(state,5,status,fs=0.0_dp,nperseg=8)
    if (status /= istft_status_invalid_input) error stop 'invalid rate'
    call istft_stream_init(state,5,status,window='hann',window_values=[(1.0_dp,k=1,8)],nperseg=8)
    if (status /= istft_status_invalid_input) error stop 'two windows'

    call istft_stream_init(state,5,status,nperseg=8,noverlap=4,nfft=8)
    if (status/=istft_status_success) error stop 'valid initialization'
    call istft_stream_init(state,5,status,nperseg=8,noverlap=8,nfft=8)
    if (status/=istft_status_invalid_input) error stop 'invalid reinit must fail'
    call istft_stream_finish(state,real_chunk,status)
    if (status/=istft_status_invalid_input .or. size(real_chunk)/=0) error stop 'empty stream finish'
    call istft_stream_push(state,[(cmplx(0.5_dp,0.0_dp,dp),k=1,4)],real_chunk,status)
    if (status/=istft_status_invalid_input .or. size(real_chunk)/=0) error stop 'wrong frame length'
    bad=cmplx(0.5_dp,0.0_dp,dp)
    bad(3)=cmplx(ieee_value(1.0_dp,ieee_quiet_nan),0.0_dp,dp)
    call istft_stream_push(state,bad,real_chunk,status)
    if (status/=istft_status_invalid_input .or. size(real_chunk)/=0) error stop 'nonfinite frame'
    bad=cmplx(0.5_dp,0.0_dp,dp)
    call istft_stream_push(state,bad,real_chunk,status,tm)
    if (status/=istft_status_success) error stop 'stream not usable after rejected frame'
    call istft_stream_finish(state,real_chunk,status)
    if (status/=istft_status_success) error stop 'finish valid single frame'
    call istft_stream_push(state,bad,real_chunk,status)
    if (status/=istft_status_invalid_input) error stop 'push after finish'
    call istft_stream_finish(state,real_chunk,status)
    if (status/=istft_status_invalid_input) error stop 'double finish'

    call istft_stream_init(state,8,status,nperseg=8,noverlap=4,input_onesided=.false.)
    if (status/=istft_status_success) error stop 'two-sided initialization'
    call istft_stream_push(state,[(cmplx(1.0_dp,0.0_dp,dp),k=1,8)],real_chunk,status)
    if (status/=istft_status_invalid_input) error stop 'real output with two-sided input'
    call istft_stream_push(state,[(cmplx(1.0_dp,0.0_dp,dp),k=1,8)],complex_chunk,status)
    if (status/=istft_status_success) error stop 'complex output after rejected real push'
    call istft_stream_finish(state,complex_chunk,status)
    if (status/=istft_status_success) error stop 'two-sided finish'

    ! Exercise explicit synthesis windows, including odd-length windows and
    ! an odd one-sided FFT which cannot be inferred from its bin count alone.
    block
        real(dp) :: w(7), accum(13)
        real(dp), allocatable :: xr(:)
        complex(dp) :: one(4,4)
        integer :: st, r, j
        do j=1,7
            w(j)=0.2_dp+0.6_dp*sin(real(j,dp))**2
        end do
        one=cmplx(0.0_dp,0.0_dp,dp)
        do j=1,4
            one(:,j)=cmplx(real(j,dp)*0.4_dp,0.09_dp*j,dp)
        end do
        call istft(one,xr,status,nperseg=7,noverlap=3,nfft=7, &
            boundary=.true.,window_values=w)
        if (status/=istft_status_success) error stop 'custom window batch'
        call istft_stream_init(state,4,status,nperseg=7,noverlap=3,nfft=7, &
            boundary=.true.,window_values=w)
        if (status/=istft_status_success) error stop 'custom window init'
        accum=0.0_dp
        r=0
        do j=1,4
            call istft_stream_push(state,one(:,j),real_chunk,st)
            if (st/=istft_status_success) error stop 'custom window push'
            accum(r+1:r+size(real_chunk))=real_chunk
            r=r+size(real_chunk)
        end do
        call istft_stream_finish(state,real_chunk,st)
        if (st/=istft_status_success) error stop 'custom window finish'
        accum(r+1:r+size(real_chunk))=real_chunk
        r=r+size(real_chunk)
        if (r/=13 .or. maxval(abs(xr-accum))>2e-12_dp) error stop 'custom window values'
    end block

    ! Long-running input must not need an array whose size grows with frames.
    call istft_stream_init(state,4,status,nperseg=7,noverlap=3,nfft=7, &
        window='boxcar',boundary=.true.)
    if (status/=istft_status_success) error stop 'stress init'
    nsum=0
    do count=1,20000
        call istft_stream_push(state,[(cmplx(0.7_dp,0.1_dp*k,dp),k=1,4)],real_chunk,status)
        if (status/=istft_status_success) error stop 'stress push'
        nsum=nsum+size(real_chunk)
    end do
    call istft_stream_finish(state,real_chunk,status)
    if (status/=istft_status_success) error stop 'stress finish'
    nsum=nsum+size(real_chunk)
    if (nsum/=7+19999*4-2*3) error stop 'stress sample count'

    print '(a,i0,a)', 'test_signal_istft_stream passed: ',cases,' batch/stream comparisons'
contains
    subroutine check_case(nw,nf,no,frames,onesided,boundary,scaling,win)
        integer, intent(in) :: nw,nf,no,frames
        logical, intent(in) :: onesided,boundary
        character(*), intent(in) :: scaling,win
        type(istft_stream) :: state
        complex(dp), allocatable :: z(:,:), yc(:), expected_c(:), chunk_c(:)
        real(dp), allocatable :: yr(:), expected_r(:), chunk_r(:), times(:), chunk_times(:)
        real(dp), allocatable :: reference_times(:), all_times(:)
        integer :: nfreq, i,j,st,batch_status
        real(dp) :: phase

        nfreq=nf
        if (onesided) nfreq=nf/2+1
        allocate(z(nfreq,frames), expected_c(0), expected_r(0), all_times(0))
        do j=1,frames
            do i=1,nfreq
                phase = real(13*i+7*j+3*nw,dp)*0.119_dp
                z(i,j) = cmplx(cos(phase)+0.3_dp*sin(0.41_dp*phase), &
                    sin(phase)+0.2_dp*cos(0.77_dp*phase), dp)
            end do
        end do
        call istft_stream_init(state,nfreq,st,fs=12.5_dp,window=win, &
            nperseg=nw,noverlap=no,nfft=nf,boundary=boundary,scaling=scaling, &
            input_onesided=onesided)
        if (st /= istft_status_success) error stop 'init failed'
        if (onesided) then
            call istft(z,yr,batch_status,fs=12.5_dp,window=win,nperseg=nw, &
                noverlap=no,nfft=nf,boundary=boundary,scaling=scaling, &
                times=reference_times,input_onesided=.true.)
        else
            call istft(z,yc,batch_status,fs=12.5_dp,window=win,nperseg=nw, &
                noverlap=no,nfft=nf,boundary=boundary,scaling=scaling, &
                times=reference_times,input_onesided=.false.)
        end if
        if (batch_status==istft_status_invalid_input) error stop 'batch invalid'
        do j=1,frames
            if (onesided) then
                call istft_stream_push(state,z(:,j),chunk_r,st,chunk_times)
                expected_r=[expected_r,chunk_r]
            else
                call istft_stream_push(state,z(:,j),chunk_c,st,chunk_times)
                expected_c=[expected_c,chunk_c]
            end if
            if (st==istft_status_invalid_input) error stop 'push invalid'
            all_times=[all_times,chunk_times]
        end do
        if (onesided) then
            call istft_stream_finish(state,chunk_r,st,chunk_times)
            expected_r=[expected_r,chunk_r]
            if (size(expected_r)/=size(yr)) error stop 'real length mismatch'
            if (size(yr)>0) then
                if (maxval(abs(expected_r-yr))>3.0e-12_dp) error stop 'real sample mismatch'
            end if
        else
            call istft_stream_finish(state,chunk_c,st,chunk_times)
            expected_c=[expected_c,chunk_c]
            if (size(expected_c)/=size(yc)) error stop 'complex length mismatch'
            if (size(yc)>0) then
                if (maxval(abs(expected_c-yc))>3.0e-12_dp) error stop 'complex sample mismatch'
            end if
        end if
        if (st /= batch_status) error stop 'final status mismatch'
        all_times=[all_times,chunk_times]
        if (size(all_times)/=size(reference_times)) error stop 'time count mismatch'
        if (size(all_times)>0) then
            if (maxval(abs(all_times-reference_times))>2.0e-13_dp) error stop 'time values mismatch'
        end if
    end subroutine check_case
end program test_signal_istft_stream
