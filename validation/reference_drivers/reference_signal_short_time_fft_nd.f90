! SPDX-License-Identifier: MIT
! Text-format reference producer for rank-two/three ShortTimeFFT against SciPy.
program reference_signal_short_time_fft_nd
    use scifort_kinds, only : dp
    use scifort_signal_short_time_fft, only : short_time_fft, short_time_fft_result_3d, &
        short_time_fft_result_4d, short_fft_success
    implicit none
    type(short_time_fft) :: plan
    type(short_time_fft_result_3d) :: s2
    type(short_time_fft_result_4d) :: s3
    real(dp), allocatable :: x2(:,:), x3(:,:,:)
    complex(dp), allocatable :: z2(:,:), z3(:,:,:)
    character(24) :: mode, scale, pad, kind, trend
    integer :: ic, ncases, rank, axis, n, nw, hop, nf, phase, p0, p1, offset
    integer :: dims(3), i, j, k, a, b, c, d
    real(dp) :: fs, value
    complex(dp) :: zvalue

    open(unit=11,file='short_fft_nd_cases.txt',status='old',action='read')
    open(unit=12,file='short_fft_nd_results.txt',status='replace',action='write')
    read(11,*) ncases
    do ic=1,ncases
        read(11,*) rank,axis,n,nw,hop,nf,fs,mode,scale,phase,pad,kind,p0,p1,offset,trend
        call plan%initialize_window('hann',nw,hop,fs,mfft=nf,fft_mode=trim(mode), &
            scale_to=trim(scale),phase_shift=phase)
        if (plan%status /= short_fft_success) error stop 'invalid plan'
        dims = [2,3,n]
        if (rank == 2) then
            if (axis == 1) then
                allocate(x2(n,3),z2(n,3))
            else
                allocate(x2(3,n),z2(3,n))
            end if
            do j=1,size(x2,2)
                do i=1,size(x2,1)
                    value = sin(0.13_dp*i+0.07_dp*j)+0.2_dp*cos(0.19_dp*i-0.1_dp*j)
                    zvalue = cmplx(value,0.17_dp*sin(0.11_dp*i+0.03_dp*j),dp)
                    x2(i,j) = value
                    z2(i,j) = zvalue
                end do
            end do
            if (kind == 'real') then
                if (trend == 'none') then
                    call plan%stft(x2,s2,p0,p1,offset,trim(pad),axis)
                else
                    call plan%stft_detrend(x2,s2,trim(trend),p0,p1,offset,trim(pad),axis)
                end if
            else
                if (trend == 'none') then
                    call plan%stft(z2,s2,p0,p1,offset,trim(pad),axis)
                else
                    call plan%stft_detrend(z2,s2,trim(trend),p0,p1,offset,trim(pad),axis)
                end if
            end if
            if (s2%status /= short_fft_success) error stop 'invalid rank-2 STFT'
            write(12,*) ic,rank,shape(s2%values),1,s2%p0
            do c=1,size(s2%values,3)
                do b=1,size(s2%values,2)
                    do a=1,size(s2%values,1)
                        write(12,'(2(es26.17,1x))') real(s2%values(a,b,c),dp), &
                            aimag(s2%values(a,b,c))
                    end do
                end do
            end do
            deallocate(x2,z2)
        else
            select case (axis)
            case (1)
                dims = [n,2,3]
            case (2)
                dims = [2,n,3]
            case (3)
                dims = [2,3,n]
            end select
            allocate(x3(dims(1),dims(2),dims(3)),z3(dims(1),dims(2),dims(3)))
            do k=1,dims(3)
                do j=1,dims(2)
                    do i=1,dims(1)
                        value = sin(0.13_dp*i+0.07_dp*j+0.09_dp*k) + &
                            0.2_dp*cos(0.19_dp*i-0.1_dp*j+0.05_dp*k)
                        zvalue = cmplx(value,0.17_dp*sin(0.11_dp*i+0.03_dp*j-0.02_dp*k),dp)
                        x3(i,j,k) = value
                        z3(i,j,k) = zvalue
                    end do
                end do
            end do
            if (kind == 'real') then
                if (trend == 'none') then
                    call plan%stft(x3,s3,p0,p1,offset,trim(pad),axis)
                else
                    call plan%stft_detrend(x3,s3,trim(trend),p0,p1,offset,trim(pad),axis)
                end if
            else
                if (trend == 'none') then
                    call plan%stft(z3,s3,p0,p1,offset,trim(pad),axis)
                else
                    call plan%stft_detrend(z3,s3,trim(trend),p0,p1,offset,trim(pad),axis)
                end if
            end if
            if (s3%status /= short_fft_success) error stop 'invalid rank-3 STFT'
            write(12,*) ic,rank,shape(s3%values),s3%p0
            do d=1,size(s3%values,4)
                do c=1,size(s3%values,3)
                    do b=1,size(s3%values,2)
                        do a=1,size(s3%values,1)
                            write(12,'(2(es26.17,1x))') real(s3%values(a,b,c,d),dp), &
                                aimag(s3%values(a,b,c,d))
                        end do
                    end do
                end do
            end do
            deallocate(x3,z3)
        end if
    end do
    close(12)
    close(11)
    print *, 'rank-2 and rank-3 ShortTimeFFT cases:',ncases
end program reference_signal_short_time_fft_nd
