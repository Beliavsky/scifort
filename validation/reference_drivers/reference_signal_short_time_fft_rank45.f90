! SPDX-License-Identifier: MIT
! SciPy reference producer for rank-four and rank-five ShortTimeFFT.
program reference_signal_short_time_fft_rank45
    use scifort_kinds, only : dp
    use scifort_signal_short_time_fft, only : short_time_fft, &
        short_time_fft_result_5d, short_time_fft_result_6d, short_fft_success
    implicit none
    type(short_time_fft) :: plan
    type(short_time_fft_result_5d) :: s4
    type(short_time_fft_result_6d) :: s5
    real(dp), allocatable :: x4(:,:,:,:), x5(:,:,:,:,:)
    complex(dp), allocatable :: z4(:,:,:,:), z5(:,:,:,:,:)
    character(24) :: mode, scale, pad, kind, trend
    integer :: ic, ncases, rank, axis, n, nw, hop, nf, phase, p0, p1, offset
    integer :: dims4(4), dims5(5), ix(5), a,b,c,d,e,f,fslice
    real(dp) :: fs, value, arg1, arg2, arg3
    real(dp), parameter :: weights1(5) = [0.13_dp,0.07_dp,0.11_dp,0.05_dp,0.03_dp]
    real(dp), parameter :: weights2(5) = [0.17_dp,-0.1_dp,0.03_dp,0.08_dp,0.07_dp]
    real(dp), parameter :: weights3(5) = [0.08_dp,0.12_dp,-0.09_dp,0.13_dp,0.05_dp]
    complex(dp) :: zvalue

    open(unit=11,file='short_fft_rank45_cases.txt',status='old',action='read')
    open(unit=12,file='short_fft_rank45_results.txt',status='replace',action='write')
    read(11,*) ncases
    do ic=1,ncases
        read(11,*) rank,axis,n,nw,hop,nf,fs,mode,scale,phase,pad,kind,p0,p1,offset,trend
        call plan%initialize_window('hann',nw,hop,fs,mfft=nf,fft_mode=trim(mode), &
            scale_to=trim(scale),phase_shift=phase)
        if (plan%status /= short_fft_success) error stop 'invalid reference plan'
        if (rank == 4) then
            dims4 = 2
            dims4(axis) = n
            allocate(x4(dims4(1),dims4(2),dims4(3),dims4(4)))
            allocate(z4(dims4(1),dims4(2),dims4(3),dims4(4)))
            do d=1,dims4(4)
                do c=1,dims4(3)
                    do b=1,dims4(2)
                        do a=1,dims4(1)
                            ix=[a,b,c,d,0]
                            arg1 = sum(weights1(:4)*ix(:4))
                            arg2 = sum(weights2(:4)*ix(:4))
                            arg3 = sum(weights3(:4)*ix(:4))
                            value = sin(arg1)+0.2_dp*cos(arg2)
                            zvalue = cmplx(value,0.17_dp*sin(arg3),dp)
                            x4(a,b,c,d) = value
                            z4(a,b,c,d) = zvalue
                        end do
                    end do
                end do
            end do
            if (kind == 'real') then
                if (trend == 'none') then
                    call plan%stft(x4,s4,p0,p1,offset,trim(pad),axis)
                else
                    call plan%stft_detrend(x4,s4,trim(trend),p0,p1,offset,trim(pad),axis)
                end if
            else
                if (trend == 'none') then
                    call plan%stft(z4,s4,p0,p1,offset,trim(pad),axis)
                else
                    call plan%stft_detrend(z4,s4,trim(trend),p0,p1,offset,trim(pad),axis)
                end if
            end if
            if (s4%status /= short_fft_success) error stop 'rank4 reference invalid'
            write(12,*) ic,rank,shape(s4%values),s4%p0
            do fslice=1,size(s4%values,5)
                do d=1,size(s4%values,4)
                    do c=1,size(s4%values,3)
                        do b=1,size(s4%values,2)
                            do a=1,size(s4%values,1)
                                write(12,'(2(es26.17,1x))') real(s4%values(a,b,c,d,fslice),dp), &
                                    aimag(s4%values(a,b,c,d,fslice))
                            end do
                        end do
                    end do
                end do
            end do
            deallocate(x4,z4)
        else
            dims5 = 2
            dims5(axis) = n
            allocate(x5(dims5(1),dims5(2),dims5(3),dims5(4),dims5(5)))
            allocate(z5(dims5(1),dims5(2),dims5(3),dims5(4),dims5(5)))
            do e=1,dims5(5)
                do d=1,dims5(4)
                    do c=1,dims5(3)
                        do b=1,dims5(2)
                            do a=1,dims5(1)
                                ix=[a,b,c,d,e]
                                arg1 = sum(weights1*ix)
                                arg2 = sum(weights2*ix)
                                arg3 = sum(weights3*ix)
                                value = sin(arg1)+0.2_dp*cos(arg2)
                                zvalue = cmplx(value,0.17_dp*sin(arg3),dp)
                                x5(a,b,c,d,e) = value
                                z5(a,b,c,d,e) = zvalue
                            end do
                        end do
                    end do
                end do
            end do
            if (kind == 'real') then
                if (trend == 'none') then
                    call plan%stft(x5,s5,p0,p1,offset,trim(pad),axis)
                else
                    call plan%stft_detrend(x5,s5,trim(trend),p0,p1,offset,trim(pad),axis)
                end if
            else
                if (trend == 'none') then
                    call plan%stft(z5,s5,p0,p1,offset,trim(pad),axis)
                else
                    call plan%stft_detrend(z5,s5,trim(trend),p0,p1,offset,trim(pad),axis)
                end if
            end if
            if (s5%status /= short_fft_success) error stop 'rank5 reference invalid'
            write(12,*) ic,rank,shape(s5%values),s5%p0
            do fslice=1,size(s5%values,6)
                do e=1,size(s5%values,5)
                    do d=1,size(s5%values,4)
                        do c=1,size(s5%values,3)
                            do b=1,size(s5%values,2)
                                do a=1,size(s5%values,1)
                                    write(12,'(2(es26.17,1x))') real(s5%values(a,b,c,d,e,fslice),dp), &
                                        aimag(s5%values(a,b,c,d,e,fslice))
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            deallocate(x5,z5)
        end if
    end do
    close(11)
    close(12)
    print *, 'rank-four/five reference cases:',ncases
end program reference_signal_short_time_fft_rank45
