! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_stream_multi
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_signal, only : dp, stft_stream, stft_stream_init, stft_stream_push, &
        stft_stream_finish, stft_stream_multi, stft_stream_multi_init, &
        stft_stream_multi_push, stft_stream_multi_finish, stft_stream_multi_frequencies, &
        stft_stream_frequencies, spectrogram_stream_init, spectrogram_stream_push, &
        spectrogram_stream_finish, spectrogram_stream_multi_init, &
        spectrogram_stream_multi_push, spectrogram_stream_multi_finish, &
        istft_stream, istft_stream_init, istft_stream_push, istft_stream_finish, &
        istft_stream_multi, istft_stream_multi_init, istft_stream_multi_push, &
        istft_stream_multi_finish, spectral_status_success, spectral_status_invalid_input, &
        istft_status_success, istft_status_invalid_input, istft_status_noninvertible
    implicit none
    integer :: nch, nw, no, parity, mode, signal_type, tests

    tests = 0
    do nch=1,3,2
        do nw=7,8
            do parity=0,1
                do no=0,1
                    do mode=0,1
                        do signal_type=0,1
                            call check_forward(nch,nw,nw+parity,no==1,mode==1,signal_type==1)
                            tests=tests+1
                        end do
                    end do
                end do
            end do
        end do
    end do
    call check_rejected_forward()
    call check_rejected_inverse()
    do nch=1,3,2
        call check_inverse(nch,.true.)
        call check_inverse(nch,.false.)
        tests=tests+2
    end do
    call check_round_trip(.false.)
    call check_round_trip(.true.)
    call check_noninvertible()
    call check_long_stream()
    print '(a,i0,a)', 'PASS multichannel streaming ',tests+6,' cases'

contains

    subroutine assert_close(a,b,name)
        complex(dp), intent(in) :: a(:,:) !! computed values
        complex(dp), intent(in) :: b(:,:) !! reference values
        character(*), intent(in) :: name !! label
        if (any(shape(a)/=shape(b))) then
            print *, 'wrong shapes ', name, shape(a), shape(b)
            error stop 'shape mismatch'
        end if
        if (size(a)>0) then
            if (maxval(abs(a-b))>2.0e-11_dp) then
                print *, 'difference ', name,maxval(abs(a-b))
                error stop 'value mismatch'
            end if
        end if
    end subroutine assert_close

    subroutine compare_times(a,b)
        real(dp), intent(in) :: a(:) !! computed frame/sample times
        real(dp), intent(in) :: b(:) !! reference times
        if (size(a)/=size(b)) error stop 'time vector length mismatch'
        if (size(a)>0) then
            if (maxval(abs(a-b))>1.0e-13_dp) error stop 'time coordinates mismatch'
        end if
    end subroutine compare_times

    subroutine check_forward(nch,nw,nf,overlap,power,is_complex)
        integer, intent(in) :: nch, nw, nf !! number of channels and FFT sizes
        logical, intent(in) :: overlap, power, is_complex !! case choices
        integer, parameter :: n = 71
        real(dp) :: x(n,nch), xstride(2*n,nch), t, c
        complex(dp) :: z(n,nch), zstride(2*n,nch)
        type(stft_stream_multi) :: state
        type(stft_stream) :: single(nch)
        complex(dp), allocatable :: frames(:,:,:), f1(:,:), cf1(:,:)
        real(dp), allocatable :: powers(:,:,:), p1(:,:), times(:), one_times(:)
        real(dp), allocatable :: frequencies(:), one_freq(:)
        integer :: no, j, i, start, stop, step, st, k

        no=0
        if (overlap) no=nw/2
        do j=1,nch
            do i=1,n
                t=real(i,dp)/13.0_dp
                c=real(j,dp)
                x(i,j)=sin(t*c)+cos(t*0.6_dp)+c/9.0_dp
                z(i,j)=cmplx(x(i,j),0.8_dp*cos(t*c*0.4_dp)-0.2_dp*c,dp)
            end do
        end do
        xstride=0.0_dp
        zstride=(0.0_dp,0.0_dp)
        xstride(1:2*n:2,:)=x
        zstride(1:2*n:2,:)=z
        if (power) then
            call spectrogram_stream_multi_init(state,nch,nw,st,fs=13.0_dp, &
                window='hann',noverlap=no,nfft=nf,detrend='linear', &
                scaling='density',input_complex=is_complex)
            if (st/=spectral_status_success) error stop 'multichannel spectrogram init'
            do j=1,nch
                call spectrogram_stream_init(single(j),nw,st,fs=13.0_dp, &
                    window='hann',noverlap=no,nfft=nf,detrend='linear', &
                    scaling='density',input_complex=is_complex)
                if (st/=spectral_status_success) error stop 'single spectrogram init'
            end do
        else
            call stft_stream_multi_init(state,nch,nw,st,fs=13.0_dp, &
                window='hann',noverlap=no,nfft=nf,detrend='constant', &
                boundary='zeros',padded=.true.,scaling='psd',input_complex=is_complex)
            if (st/=spectral_status_success) error stop 'multichannel STFT init'
            do j=1,nch
                call stft_stream_init(single(j),nw,st,fs=13.0_dp, &
                    window='hann',noverlap=no,nfft=nf,detrend='constant', &
                    boundary='zeros',padded=.true.,scaling='psd',input_complex=is_complex)
                if (st/=spectral_status_success) error stop 'single STFT init'
            end do
        end if
        call stft_stream_multi_frequencies(state,frequencies,st)
        if (st/=spectral_status_success) error stop 'multichannel frequencies'
        call stft_stream_frequencies(single(1),one_freq,st)
        call compare_times(frequencies,one_freq)

        ! Irregular chunks including a zero-length chunk and noncontiguous views.
        start=1
        do k=0,17
            step=mod(k*7+3,11)
            if (k==0) step=0
            stop=min(n,start+step-1)
            if (power) then
                if (is_complex) then
                    if (mod(k,2)==0) then
                        call spectrogram_stream_multi_push(state,zstride(2*start-1:2*stop-1:2,:),powers,st,times)
                    else
                        call spectrogram_stream_multi_push(state,z(start:stop,:),powers,st,times)
                    end if
                else
                    if (mod(k,2)==0) then
                        call spectrogram_stream_multi_push(state,xstride(2*start-1:2*stop-1:2,:),powers,st,times)
                    else
                        call spectrogram_stream_multi_push(state,x(start:stop,:),powers,st,times)
                    end if
                end if
                if (st/=spectral_status_success) error stop 'spectrogram push'
                if (size(powers,2)/=nch .or. size(powers,3)/=size(times)) &
                    error stop 'spectrogram shape'
                do j=1,nch
                    if (is_complex) then
                        call spectrogram_stream_push(single(j),z(start:stop,j),p1,st,one_times)
                    else
                        call spectrogram_stream_push(single(j),x(start:stop,j),p1,st,one_times)
                    end if
                    if (st/=spectral_status_success) error stop 'single spectrogram push'
                    call compare_times(times,one_times)
                    call assert_close(cmplx(powers(:,j,:),0.0_dp,dp), &
                        cmplx(p1,0.0_dp,dp),'spectrogram channel')
                end do
            else
                if (is_complex) then
                    if (mod(k,2)==0) then
                        call stft_stream_multi_push(state,zstride(2*start-1:2*stop-1:2,:),frames,st,times)
                    else
                        call stft_stream_multi_push(state,z(start:stop,:),frames,st,times)
                    end if
                else
                    if (mod(k,2)==0) then
                        call stft_stream_multi_push(state,xstride(2*start-1:2*stop-1:2,:),frames,st,times)
                    else
                        call stft_stream_multi_push(state,x(start:stop,:),frames,st,times)
                    end if
                end if
                if (st/=spectral_status_success) error stop 'STFT push'
                if (size(frames,2)/=nch .or. size(frames,3)/=size(times)) &
                    error stop 'STFT shape'
                do j=1,nch
                    if (is_complex) then
                        call stft_stream_push(single(j),z(start:stop,j),f1,st,one_times)
                    else
                        call stft_stream_push(single(j),x(start:stop,j),f1,st,one_times)
                    end if
                    if (st/=spectral_status_success) error stop 'single STFT push'
                    call compare_times(times,one_times)
                    call assert_close(frames(:,j,:),f1,'STFT channel')
                end do
            end if
            start=stop+1
            if (start>n) exit
        end do
        if (start<=n) then
            if (power) then
                if (is_complex) then
                    call spectrogram_stream_multi_push(state,zstride(2*start-1:2*n-1:2,:),powers,st,times)
                else
                    call spectrogram_stream_multi_push(state,xstride(2*start-1:2*n-1:2,:),powers,st,times)
                end if
            else
                if (is_complex) then
                    call stft_stream_multi_push(state,zstride(2*start-1:2*n-1:2,:),frames,st,times)
                else
                    call stft_stream_multi_push(state,xstride(2*start-1:2*n-1:2,:),frames,st,times)
                end if
            end if
            if (st/=spectral_status_success) error stop 'last data push'
            do j=1,nch
                if (power) then
                    if (is_complex) then
                        call spectrogram_stream_push(single(j),z(start:n,j),p1,st,one_times)
                    else
                        call spectrogram_stream_push(single(j),x(start:n,j),p1,st,one_times)
                    end if
                    call assert_close(cmplx(powers(:,j,:),0.0_dp,dp), &
                        cmplx(p1,0.0_dp,dp),'last spectrogram')
                else
                    if (is_complex) then
                        call stft_stream_push(single(j),z(start:n,j),f1,st,one_times)
                    else
                        call stft_stream_push(single(j),x(start:n,j),f1,st,one_times)
                    end if
                    call assert_close(frames(:,j,:),f1,'last STFT')
                end if
                if (st/=spectral_status_success) error stop 'last single data push'
                call compare_times(times,one_times)
            end do
        end if

        if (power) then
            call spectrogram_stream_multi_finish(state,powers,st,times)
            if (st/=spectral_status_success) error stop 'spectrogram finish'
        else
            call stft_stream_multi_finish(state,frames,st,times)
            if (st/=spectral_status_success) error stop 'STFT finish'
        end if
        do j=1,nch
            if (power) then
                call spectrogram_stream_finish(single(j),p1,st,one_times)
                call assert_close(cmplx(powers(:,j,:),0.0_dp,dp), &
                    cmplx(p1,0.0_dp,dp),'spectrogram finish')
            else
                call stft_stream_finish(single(j),f1,st,one_times)
                call assert_close(frames(:,j,:),f1,'STFT finish')
            end if
            if (st/=spectral_status_success) error stop 'single finish'
            call compare_times(times,one_times)
        end do
        if (power) then
            call spectrogram_stream_multi_finish(state,powers,st)
        else
            call stft_stream_multi_finish(state,frames,st)
        end if
        if (st/=spectral_status_invalid_input) error stop 'double finish accepted'
    end subroutine check_forward

    subroutine check_inverse(nch,onesided)
        integer, intent(in) :: nch !! number of independent signals
        logical, intent(in) :: onesided !! real or complex output mode
        type(istft_stream_multi) :: state
        type(istft_stream) :: single(nch)
        complex(dp) :: frames(8,nch), frame_stride(16,2*nch)
        complex(dp), allocatable :: cm(:,:), cz(:)
        real(dp), allocatable :: rm(:,:), rz(:), times(:), one_times(:)
        integer :: nf, nw, no, st, j, k, i

        nf=8
        nw=8
        no=4
        if (onesided) then
            nf=4
            nw=7
            no=3
        end if
        call istft_stream_multi_init(state,nch,nf,st,window='boxcar', &
            nperseg=nw,noverlap=no,nfft=nw,input_onesided=onesided,fs=11.0_dp)
        if (st/=istft_status_success) error stop 'multichannel inverse init'
        do j=1,nch
            call istft_stream_init(single(j),nf,st,window='boxcar', &
                nperseg=nw,noverlap=no,nfft=nw,input_onesided=onesided,fs=11.0_dp)
            if (st/=istft_status_success) error stop 'single inverse init'
        end do
        do k=1,16
            frames=(0.0_dp,0.0_dp)
            do j=1,nch
                do i=1,nf
                    frames(i,j)=cmplx(sin(real(i*k*j,dp)*0.11_dp), &
                        cos(real(i+k+j,dp)*0.15_dp),dp)
                end do
            end do
            frame_stride=(0.0_dp,0.0_dp)
            frame_stride(1:2*nf:2,1:2*nch:2)=frames(1:nf,:)
            if (onesided) then
                if (mod(k,2)==0) then
                    call istft_stream_multi_push(state, &
                        frame_stride(1:2*nf:2,1:2*nch:2),rm,st,times)
                else
                    call istft_stream_multi_push(state,frames(1:nf,:),rm,st,times)
                end if
            else
                if (mod(k,2)==0) then
                    call istft_stream_multi_push(state, &
                        frame_stride(1:2*nf:2,1:2*nch:2),cm,st,times)
                else
                    call istft_stream_multi_push(state,frames(1:nf,:),cm,st,times)
                end if
            end if
            if (st/=istft_status_success) error stop 'multichannel inverse push'
            do j=1,nch
                if (onesided) then
                    call istft_stream_push(single(j),frames(1:nf,j),rz,st,one_times)
                    call assert_close(cmplx(rm(:,j:j),0.0_dp,dp), &
                        cmplx(reshape(rz,[size(rz),1]),0.0_dp,dp),'real inverse')
                else
                    call istft_stream_push(single(j),frames(1:nf,j),cz,st,one_times)
                    call assert_close(cm(:,j:j),reshape(cz,[size(cz),1]),'complex inverse')
                end if
                if (st/=istft_status_success) error stop 'single inverse push'
                call compare_times(times,one_times)
            end do
        end do
        if (onesided) then
            call istft_stream_multi_finish(state,rm,st,times)
        else
            call istft_stream_multi_finish(state,cm,st,times)
        end if
        if (st/=istft_status_success) error stop 'multichannel inverse finish'
        do j=1,nch
            if (onesided) then
                call istft_stream_finish(single(j),rz,st,one_times)
                call assert_close(cmplx(rm(:,j:j),0.0_dp,dp), &
                    cmplx(reshape(rz,[size(rz),1]),0.0_dp,dp),'real inverse finish')
            else
                call istft_stream_finish(single(j),cz,st,one_times)
                call assert_close(cm(:,j:j),reshape(cz,[size(cz),1]),'complex inverse finish')
            end if
            if (st/=istft_status_success) error stop 'single inverse finish'
            call compare_times(times,one_times)
        end do
        if (onesided) then
            call istft_stream_multi_finish(state,rm,st)
        else
            call istft_stream_multi_finish(state,cm,st)
        end if
        if (st/=istft_status_invalid_input) error stop 'double inverse finish accepted'
    end subroutine check_inverse

    subroutine check_rejected_forward()
        type(stft_stream_multi) :: state
        real(dp) :: x(22,3), bad(2,3)
        complex(dp) :: z(2,3)
        complex(dp), allocatable :: frames(:,:,:)
        real(dp), allocatable :: power(:,:,:)
        integer :: st, j

        call stft_stream_multi_init(state,0,8,st)
        if (st/=spectral_status_invalid_input) error stop 'zero channel init'
        call stft_stream_multi_init(state,3,8,st,window='boxcar',noverlap=4)
        if (st/=spectral_status_success) error stop 'valid init'
        call stft_stream_multi_init(state,3,8,st,noverlap=8)
        if (st/=spectral_status_invalid_input) error stop 'bad overlap accepted'
        call stft_stream_multi_push(state,x(1:2,1:2),frames,st)
        if (st/=spectral_status_invalid_input) error stop 'wrong channel count'
        call stft_stream_multi_push(state,z,frames,st)
        if (st/=spectral_status_invalid_input) error stop 'wrong sample type'
        bad=0.0_dp
        bad(2,3)=ieee_value(0.0_dp,ieee_quiet_nan)
        call stft_stream_multi_push(state,bad,frames,st)
        if (st/=spectral_status_invalid_input) error stop 'nonfinite channel accepted'
        if (size(frames)/=0) error stop 'error returns frames'
        x=0.0_dp
        do j=1,3
            x(:,j)=real(j,dp)
        end do
        call stft_stream_multi_push(state,x,frames,st)
        if (st/=spectral_status_success) error stop 'not usable after rejected input'
        if (size(frames,2)/=3) error stop 'frame channels'
        call stft_stream_multi_finish(state,frames,st)
        if (st/=spectral_status_success) error stop 'finish after rejection'
        call stft_stream_multi_push(state,x(1:1,:),frames,st)
        if (st/=spectral_status_invalid_input) error stop 'push after finish accepted'
        call spectrogram_stream_multi_init(state,3,8,st,window='boxcar')
        if (st/=spectral_status_success) error stop 'spectrogram reinit'
        call stft_stream_multi_push(state,x(1:2,:),frames,st)
        if (st/=spectral_status_invalid_input) error stop 'wrong spectral mode'
        call spectrogram_stream_multi_push(state,x(1:2,:),power,st)
        if (st/=spectral_status_success) error stop 'spectrogram recovered'
    end subroutine check_rejected_forward

    subroutine check_rejected_inverse()
        type(istft_stream_multi) :: state
        complex(dp) :: frames(5,3)
        complex(dp), allocatable :: cm(:,:)
        real(dp), allocatable :: rm(:,:)
        integer :: st, j

        call istft_stream_multi_init(state,0,5,st,nperseg=8,nfft=8)
        if (st/=istft_status_invalid_input) error stop 'inverse zero channel accepted'
        call istft_stream_multi_init(state,3,5,st,nperseg=8,nfft=8,window='boxcar')
        if (st/=istft_status_success) error stop 'inverse valid init'
        call istft_stream_multi_init(state,3,5,st,nperseg=9,nfft=8)
        if (st/=istft_status_invalid_input) error stop 'inverse bad config accepted'
        frames=(1.0_dp,0.0_dp)
        call istft_stream_multi_push(state,frames(:,1:2),rm,st)
        if (st/=istft_status_invalid_input) error stop 'inverse channel mismatch'
        frames(5,3)=cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
        call istft_stream_multi_push(state,frames,rm,st)
        if (st/=istft_status_invalid_input) error stop 'inverse nonfinite accepted'
        if (size(rm)/=0) error stop 'inverse error returns samples'
        frames=(1.0_dp,0.0_dp)
        call istft_stream_multi_push(state,frames,rm,st)
        if (st/=istft_status_success) error stop 'inverse not usable after reject'
        call istft_stream_multi_finish(state,rm,st)
        if (st/=istft_status_success) error stop 'inverse finished after rejection'
        call istft_stream_multi_init(state,3,8,st,nperseg=8,nfft=8, &
            input_onesided=.false.,window='boxcar')
        if (st/=istft_status_success) error stop 'two-sided inverse init'
        call istft_stream_multi_push(state,spread(frames(1:5,1),2,3),rm,st)
        if (st/=istft_status_invalid_input) error stop 'bad two-sided frame accepted'
        call istft_stream_multi_push(state,spread(cmplx([(real(j,dp),j=1,8)], &
            0.0_dp,dp),2,3),rm,st)
        if (st/=istft_status_invalid_input) error stop 'real output two-sided accepted'
    end subroutine check_rejected_inverse

    subroutine check_round_trip(is_complex)
        logical, intent(in) :: is_complex !! real or complex reconstruction
        integer, parameter :: n = 128, nch = 3
        type(stft_stream_multi) :: analysis
        type(istft_stream_multi) :: synthesis
        real(dp) :: x(n,nch), recovered(n+16,nch), yr(nch)
        complex(dp) :: z(n,nch), reconstruct(n+16,nch)
        complex(dp), allocatable :: frames(:,:,:), cm(:,:)
        real(dp), allocatable :: rm(:,:)
        integer :: i, j, k, st, first, last, total

        do j=1,nch
            do i=1,n
                x(i,j)=cos(real(i*j,dp)*0.054_dp)+0.2_dp*sin(real(i+3*j,dp)*0.2_dp)
                z(i,j)=cmplx(x(i,j),0.4_dp*sin(real(i*j,dp)*0.063_dp),dp)
            end do
        end do
        call stft_stream_multi_init(analysis,nch,8,st,window='hann',noverlap=4, &
            input_complex=is_complex)
        if (st/=spectral_status_success) error stop 'roundtrip forward init'
        if (is_complex) then
            call istft_stream_multi_init(synthesis,nch,8,st,nperseg=8,noverlap=4, &
                input_onesided=.false.,window='hann')
        else
            call istft_stream_multi_init(synthesis,nch,5,st,nperseg=8,noverlap=4, &
                window='hann')
        end if
        if (st/=istft_status_success) error stop 'roundtrip inverse init'
        total=0
        first=1
        do while (first<=n)
            last=min(n,first+mod(first*7,17))
            if (is_complex) then
                call stft_stream_multi_push(analysis,z(first:last,:),frames,st)
            else
                call stft_stream_multi_push(analysis,x(first:last,:),frames,st)
            end if
            if (st/=spectral_status_success) error stop 'roundtrip forward push'
            call consume_frames(synthesis,frames,is_complex,recovered,reconstruct,total,st)
            first=last+1
        end do
        call stft_stream_multi_finish(analysis,frames,st)
        if (st/=spectral_status_success) error stop 'roundtrip forward finish'
        call consume_frames(synthesis,frames,is_complex,recovered,reconstruct,total,st)
        if (is_complex) then
            call istft_stream_multi_finish(synthesis,cm,st)
            if (st/=istft_status_success) error stop 'roundtrip inverse finish complex'
            reconstruct(total+1:total+size(cm,1),:)=cm
            total=total+size(cm,1)
            if (total<n) error stop 'short complex reconstruction'
            if (maxval(abs(reconstruct(1:n,:)-z))>1e-12_dp) &
                error stop 'complex roundtrip mismatch'
        else
            call istft_stream_multi_finish(synthesis,rm,st)
            if (st/=istft_status_success) error stop 'roundtrip inverse finish real'
            recovered(total+1:total+size(rm,1),:)=rm
            total=total+size(rm,1)
            if (total<n) error stop 'short real reconstruction'
            if (maxval(abs(recovered(1:n,:)-x))>1e-12_dp) &
                error stop 'real roundtrip mismatch'
        end if
    end subroutine check_round_trip

    subroutine consume_frames(synthesis,frames,is_complex,recovered,reconstruct,total,st)
        type(istft_stream_multi), intent(inout) :: synthesis !! inverse transform stream
        complex(dp), intent(in) :: frames(:,:,:) !! spectral frames to reconstruct
        logical, intent(in) :: is_complex !! complex input signal
        real(dp), intent(inout) :: recovered(:,:) !! accumulated real output
        complex(dp), intent(inout) :: reconstruct(:,:) !! accumulated complex output
        integer, intent(inout) :: total !! samples already reconstructed
        integer, intent(out) :: st !! output status
        complex(dp), allocatable :: cm(:,:)
        real(dp), allocatable :: rm(:,:)
        integer :: f

        do f=1,size(frames,3)
            if (is_complex) then
                call istft_stream_multi_push(synthesis,frames(:,:,f),cm,st)
                if (st/=istft_status_success) error stop 'roundtrip inverse push complex'
                reconstruct(total+1:total+size(cm,1),:)=cm
                total=total+size(cm,1)
            else
                call istft_stream_multi_push(synthesis,frames(:,:,f),rm,st)
                if (st/=istft_status_success) error stop 'roundtrip inverse push real'
                recovered(total+1:total+size(rm,1),:)=rm
                total=total+size(rm,1)
            end if
        end do
    end subroutine consume_frames

    subroutine check_noninvertible()
        type(istft_stream_multi) :: state
        real(dp) :: w(4)
        real(dp), allocatable :: samples(:,:)
        complex(dp) :: frame(3,3)
        integer :: st

        w=[0.0_dp,0.0_dp,0.0_dp,1.0_dp]
        call istft_stream_multi_init(state,3,3,st,nperseg=4,nfft=4, &
            noverlap=0,window_values=w,boundary=.false.)
        if (st/=istft_status_success) error stop 'noninvertible init'
        frame=(1.0_dp,0.0_dp)
        call istft_stream_multi_push(state,frame,samples,st)
        if (st/=istft_status_noninvertible) error stop 'noninvertible warning lost'
        if (size(samples,1)/=4 .or. size(samples,2)/=3) error stop 'noninvertible dimensions'
        call istft_stream_multi_finish(state,samples,st)
        if (st/=istft_status_noninvertible) error stop 'noninvertible finish status lost'
    end subroutine check_noninvertible

    subroutine check_long_stream()
        type(stft_stream_multi) :: analysis
        type(istft_stream_multi) :: synthesis
        real(dp), allocatable :: x(:,:), time(:), samples(:,:)
        complex(dp), allocatable :: frames(:,:,:)
        integer :: st, j, i, frames_seen, reconstructed

        allocate(x(5,4))
        call stft_stream_multi_init(analysis,4,7,st,window='boxcar',noverlap=5, &
            nfft=9,boundary='none',padded=.false.)
        if (st/=spectral_status_success) error stop 'long forward init'
        call istft_stream_multi_init(synthesis,4,5,st,window='boxcar',nperseg=7, &
            noverlap=5,nfft=9,boundary=.false.)
        if (st/=istft_status_success) error stop 'long inverse init'
        frames_seen=0
        reconstructed=0
        do i=1,2000
            do j=1,4
                x(:,j)=sin(real(i*j,dp)*0.17_dp)
            end do
            call stft_stream_multi_push(analysis,x,frames,st,time)
            if (st/=spectral_status_success) error stop 'long forward push'
            frames_seen=frames_seen+size(frames,3)
            do j=1,size(frames,3)
                call istft_stream_multi_push(synthesis,frames(:,:,j),samples,st)
                if (st/=istft_status_success) error stop 'long inverse push'
                reconstructed=reconstructed+size(samples,1)
            end do
        end do
        call stft_stream_multi_finish(analysis,frames,st)
        if (st/=spectral_status_success) error stop 'long forward finish'
        do j=1,size(frames,3)
            call istft_stream_multi_push(synthesis,frames(:,:,j),samples,st)
            if (st/=istft_status_success) error stop 'long inverse final push'
            reconstructed=reconstructed+size(samples,1)
        end do
        call istft_stream_multi_finish(synthesis,samples,st)
        if (st/=istft_status_success) error stop 'long inverse finish'
        reconstructed=reconstructed+size(samples,1)
        if (frames_seen<4000 .or. reconstructed<9900) &
            error stop 'long stream did not advance'
    end subroutine check_long_stream

end program test_signal_stream_multi
