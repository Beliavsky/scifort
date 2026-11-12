! SPDX-License-Identifier: MIT
! SciPy comparison driver for input ranks 6--14 and spectrum ranks 7--15.
program reference_signal_short_time_fft_high_rank
    use scifort_kinds, only : dp
    use scifort_signal_short_time_fft
    implicit none
    type(short_time_fft) :: sft
    character(24) :: mode, scale, pad, kind, trend
    integer :: ic, ncases, rank, axis, n, nw, hop, nf, phase, p0, p1, offset
    integer :: uin, uout
    real(dp) :: fs

    open(newunit=uin,file='stft_high_rank_cases.txt',status='old',action='read')
    open(newunit=uout,file='stft_high_rank_results.txt',status='replace',action='write')
    read(uin,*) ncases
    do ic=1,ncases
        read(uin,*) rank,axis,n,nw,hop,nf,fs,mode,scale,phase,pad,kind,p0,p1,offset,trend
        call sft%initialize_window('hann',nw,hop,fs,mfft=nf,fft_mode=trim(mode), &
            scale_to=trim(scale),phase_shift=phase)
        if (sft%status /= short_fft_success) error stop 'invalid high-rank plan'
        select case(rank)
        case(6)
            call run_rank6(sft,axis,n,p0,p1,offset,pad,kind,trend,uout,ic)
        case(7)
            call run_rank7(sft,axis,n,p0,p1,offset,pad,kind,trend,uout,ic)
        case(8)
            call run_rank8(sft,axis,n,p0,p1,offset,pad,kind,trend,uout,ic)
        case(9)
            call run_rank9(sft,axis,n,p0,p1,offset,pad,kind,trend,uout,ic)
        case(10)
            call run_rank10(sft,axis,n,p0,p1,offset,pad,kind,trend,uout,ic)
        case(11)
            call run_rank11(sft,axis,n,p0,p1,offset,pad,kind,trend,uout,ic)
        case(12)
            call run_rank12(sft,axis,n,p0,p1,offset,pad,kind,trend,uout,ic)
        case(13)
            call run_rank13(sft,axis,n,p0,p1,offset,pad,kind,trend,uout,ic)
        case(14)
            call run_rank14(sft,axis,n,p0,p1,offset,pad,kind,trend,uout,ic)
        case default
            error stop 'invalid rank'
        end select
    end do
    close(uin)
    close(uout)
    print '(a,i0)', 'High-rank STFT reference cases: ',ncases
contains
    subroutine run_rank6(sft,axis,n,p0,p1,offset,pad,kind,trend,uout,ic)
        type(short_time_fft), intent(in) :: sft !! configured STFT
        integer, intent(in) :: axis,n,p0,p1,offset,uout,ic !! case parameters
        character(*), intent(in) :: pad,kind,trend !! case settings
        type(short_time_fft_result_7d) :: frame
        integer :: dims(6), j, k
        real(dp), allocatable :: x(:,:,:,:,:,:)
        complex(dp), allocatable :: z(:,:,:,:,:,:), flat(:)

        dims = 1
        dims(axis) = n
        j = 1
        if (j == axis) j = 2
        dims(j) = 2
        j = 6
        if (j == axis) j = 5
        dims(j) = 2
        allocate(x(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)), z(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)))
        x = reshape([(sin(0.071_dp*real(k,dp)) + &
            0.2_dp*cos(0.13_dp*real(k,dp)),k=1,product(dims))],dims)
        z = cmplx(x,0.15_dp*reshape([(sin(0.11_dp*real(k,dp)), &
            k=1,product(dims))],dims),dp)
        if (kind == 'real') then
            if (trend == 'none') then
                call sft%stft(x,frame,p0=p0,p1=p1,k_offset=offset,padding=trim(pad),axis=axis)
            else
                call sft%stft_detrend(x,frame,trim(trend),p0=p0,p1=p1, &
                    k_offset=offset,padding=trim(pad),axis=axis)
            end if
        else
            if (trend == 'none') then
                call sft%stft(z,frame,p0=p0,p1=p1,k_offset=offset,padding=trim(pad),axis=axis)
            else
                call sft%stft_detrend(z,frame,trim(trend),p0=p0,p1=p1, &
                    k_offset=offset,padding=trim(pad),axis=axis)
            end if
        end if
        if (frame%status /= short_fft_success) error stop 'high-rank reference STFT'
        write(uout,'(*(i0,1x))') ic,6,shape(frame%values),frame%p0
        flat = reshape(frame%values,[size(frame%values)])
        do k=1,size(flat)
            write(uout,'(2(es26.17,1x))') real(flat(k),dp), aimag(flat(k))
        end do
    end subroutine run_rank6

    subroutine run_rank7(sft,axis,n,p0,p1,offset,pad,kind,trend,uout,ic)
        type(short_time_fft), intent(in) :: sft !! configured STFT
        integer, intent(in) :: axis,n,p0,p1,offset,uout,ic !! case parameters
        character(*), intent(in) :: pad,kind,trend !! case settings
        type(short_time_fft_result_8d) :: frame
        integer :: dims(7), j, k
        real(dp), allocatable :: x(:,:,:,:,:,:,:)
        complex(dp), allocatable :: z(:,:,:,:,:,:,:), flat(:)

        dims = 1
        dims(axis) = n
        j = 1
        if (j == axis) j = 2
        dims(j) = 2
        j = 7
        if (j == axis) j = 6
        dims(j) = 2
        allocate(x(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7)), z(dims(1),dims(2),dims(3), &
            dims(4),dims(5),dims(6),dims(7)))
        x = reshape([(sin(0.071_dp*real(k,dp)) + &
            0.2_dp*cos(0.13_dp*real(k,dp)),k=1,product(dims))],dims)
        z = cmplx(x,0.15_dp*reshape([(sin(0.11_dp*real(k,dp)), &
            k=1,product(dims))],dims),dp)
        if (kind == 'real') then
            if (trend == 'none') then
                call sft%stft(x,frame,p0=p0,p1=p1,k_offset=offset,padding=trim(pad),axis=axis)
            else
                call sft%stft_detrend(x,frame,trim(trend),p0=p0,p1=p1, &
                    k_offset=offset,padding=trim(pad),axis=axis)
            end if
        else
            if (trend == 'none') then
                call sft%stft(z,frame,p0=p0,p1=p1,k_offset=offset,padding=trim(pad),axis=axis)
            else
                call sft%stft_detrend(z,frame,trim(trend),p0=p0,p1=p1, &
                    k_offset=offset,padding=trim(pad),axis=axis)
            end if
        end if
        if (frame%status /= short_fft_success) error stop 'high-rank reference STFT'
        write(uout,'(*(i0,1x))') ic,7,shape(frame%values),frame%p0
        flat = reshape(frame%values,[size(frame%values)])
        do k=1,size(flat)
            write(uout,'(2(es26.17,1x))') real(flat(k),dp), aimag(flat(k))
        end do
    end subroutine run_rank7

    subroutine run_rank8(sft,axis,n,p0,p1,offset,pad,kind,trend,uout,ic)
        type(short_time_fft), intent(in) :: sft !! configured STFT
        integer, intent(in) :: axis,n,p0,p1,offset,uout,ic !! case parameters
        character(*), intent(in) :: pad,kind,trend !! case settings
        type(short_time_fft_result_9d) :: frame
        integer :: dims(8), j, k
        real(dp), allocatable :: x(:,:,:,:,:,:,:,:)
        complex(dp), allocatable :: z(:,:,:,:,:,:,:,:), flat(:)

        dims = 1
        dims(axis) = n
        j = 1
        if (j == axis) j = 2
        dims(j) = 2
        j = 8
        if (j == axis) j = 7
        dims(j) = 2
        allocate(x(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8)), z(dims(1),dims(2), &
            dims(3),dims(4),dims(5),dims(6),dims(7),dims(8)))
        x = reshape([(sin(0.071_dp*real(k,dp)) + &
            0.2_dp*cos(0.13_dp*real(k,dp)),k=1,product(dims))],dims)
        z = cmplx(x,0.15_dp*reshape([(sin(0.11_dp*real(k,dp)), &
            k=1,product(dims))],dims),dp)
        if (kind == 'real') then
            if (trend == 'none') then
                call sft%stft(x,frame,p0=p0,p1=p1,k_offset=offset,padding=trim(pad),axis=axis)
            else
                call sft%stft_detrend(x,frame,trim(trend),p0=p0,p1=p1, &
                    k_offset=offset,padding=trim(pad),axis=axis)
            end if
        else
            if (trend == 'none') then
                call sft%stft(z,frame,p0=p0,p1=p1,k_offset=offset,padding=trim(pad),axis=axis)
            else
                call sft%stft_detrend(z,frame,trim(trend),p0=p0,p1=p1, &
                    k_offset=offset,padding=trim(pad),axis=axis)
            end if
        end if
        if (frame%status /= short_fft_success) error stop 'high-rank reference STFT'
        write(uout,'(*(i0,1x))') ic,8,shape(frame%values),frame%p0
        flat = reshape(frame%values,[size(frame%values)])
        do k=1,size(flat)
            write(uout,'(2(es26.17,1x))') real(flat(k),dp), aimag(flat(k))
        end do
    end subroutine run_rank8

    subroutine run_rank9(sft,axis,n,p0,p1,offset,pad,kind,trend,uout,ic)
        type(short_time_fft), intent(in) :: sft !! configured STFT
        integer, intent(in) :: axis,n,p0,p1,offset,uout,ic !! case parameters
        character(*), intent(in) :: pad,kind,trend !! case settings
        type(short_time_fft_result_10d) :: frame
        integer :: dims(9), j, k
        real(dp), allocatable :: x(:,:,:,:,:,:,:,:,:)
        complex(dp), allocatable :: z(:,:,:,:,:,:,:,:,:), flat(:)

        dims = 1
        dims(axis) = n
        j = 1
        if (j == axis) j = 2
        dims(j) = 2
        j = 9
        if (j == axis) j = 8
        dims(j) = 2
        allocate(x(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9)), z(dims(1), &
            dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9)))
        x = reshape([(sin(0.071_dp*real(k,dp)) + &
            0.2_dp*cos(0.13_dp*real(k,dp)),k=1,product(dims))],dims)
        z = cmplx(x,0.15_dp*reshape([(sin(0.11_dp*real(k,dp)), &
            k=1,product(dims))],dims),dp)
        if (kind == 'real') then
            if (trend == 'none') then
                call sft%stft(x,frame,p0=p0,p1=p1,k_offset=offset,padding=trim(pad),axis=axis)
            else
                call sft%stft_detrend(x,frame,trim(trend),p0=p0,p1=p1, &
                    k_offset=offset,padding=trim(pad),axis=axis)
            end if
        else
            if (trend == 'none') then
                call sft%stft(z,frame,p0=p0,p1=p1,k_offset=offset,padding=trim(pad),axis=axis)
            else
                call sft%stft_detrend(z,frame,trim(trend),p0=p0,p1=p1, &
                    k_offset=offset,padding=trim(pad),axis=axis)
            end if
        end if
        if (frame%status /= short_fft_success) error stop 'high-rank reference STFT'
        write(uout,'(*(i0,1x))') ic,9,shape(frame%values),frame%p0
        flat = reshape(frame%values,[size(frame%values)])
        do k=1,size(flat)
            write(uout,'(2(es26.17,1x))') real(flat(k),dp), aimag(flat(k))
        end do
    end subroutine run_rank9

    subroutine run_rank10(sft,axis,n,p0,p1,offset,pad,kind,trend,uout,ic)
        type(short_time_fft), intent(in) :: sft !! configured STFT
        integer, intent(in) :: axis,n,p0,p1,offset,uout,ic !! case parameters
        character(*), intent(in) :: pad,kind,trend !! case settings
        type(short_time_fft_result_11d) :: frame
        integer :: dims(10), j, k
        real(dp), allocatable :: x(:,:,:,:,:,:,:,:,:,:)
        complex(dp), allocatable :: z(:,:,:,:,:,:,:,:,:,:), flat(:)

        dims = 1
        dims(axis) = n
        j = 1
        if (j == axis) j = 2
        dims(j) = 2
        j = 10
        if (j == axis) j = 9
        dims(j) = 2
        allocate(x(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9),dims(10)), &
            z(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9),dims(10)))
        x = reshape([(sin(0.071_dp*real(k,dp)) + &
            0.2_dp*cos(0.13_dp*real(k,dp)),k=1,product(dims))],dims)
        z = cmplx(x,0.15_dp*reshape([(sin(0.11_dp*real(k,dp)), &
            k=1,product(dims))],dims),dp)
        if (kind == 'real') then
            if (trend == 'none') then
                call sft%stft(x,frame,p0=p0,p1=p1,k_offset=offset,padding=trim(pad),axis=axis)
            else
                call sft%stft_detrend(x,frame,trim(trend),p0=p0,p1=p1, &
                    k_offset=offset,padding=trim(pad),axis=axis)
            end if
        else
            if (trend == 'none') then
                call sft%stft(z,frame,p0=p0,p1=p1,k_offset=offset,padding=trim(pad),axis=axis)
            else
                call sft%stft_detrend(z,frame,trim(trend),p0=p0,p1=p1, &
                    k_offset=offset,padding=trim(pad),axis=axis)
            end if
        end if
        if (frame%status /= short_fft_success) error stop 'high-rank reference STFT'
        write(uout,'(*(i0,1x))') ic,10,shape(frame%values),frame%p0
        flat = reshape(frame%values,[size(frame%values)])
        do k=1,size(flat)
            write(uout,'(2(es26.17,1x))') real(flat(k),dp), aimag(flat(k))
        end do
    end subroutine run_rank10

    subroutine run_rank11(sft,axis,n,p0,p1,offset,pad,kind,trend,uout,ic)
        type(short_time_fft), intent(in) :: sft !! configured STFT
        integer, intent(in) :: axis,n,p0,p1,offset,uout,ic !! case parameters
        character(*), intent(in) :: pad,kind,trend !! case settings
        type(short_time_fft_result_12d) :: frame
        integer :: dims(11), j, k
        real(dp), allocatable :: x(:,:,:,:,:,:,:,:,:,:,:)
        complex(dp), allocatable :: z(:,:,:,:,:,:,:,:,:,:,:), flat(:)

        dims = 1
        dims(axis) = n
        j = 1
        if (j == axis) j = 2
        dims(j) = 2
        j = 11
        if (j == axis) j = 10
        dims(j) = 2
        allocate(x(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9),dims(10), &
            dims(11)), z(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9),dims(10),dims(11)))
        x = reshape([(sin(0.071_dp*real(k,dp)) + &
            0.2_dp*cos(0.13_dp*real(k,dp)),k=1,product(dims))],dims)
        z = cmplx(x,0.15_dp*reshape([(sin(0.11_dp*real(k,dp)), &
            k=1,product(dims))],dims),dp)
        if (kind == 'real') then
            if (trend == 'none') then
                call sft%stft(x,frame,p0=p0,p1=p1,k_offset=offset,padding=trim(pad),axis=axis)
            else
                call sft%stft_detrend(x,frame,trim(trend),p0=p0,p1=p1, &
                    k_offset=offset,padding=trim(pad),axis=axis)
            end if
        else
            if (trend == 'none') then
                call sft%stft(z,frame,p0=p0,p1=p1,k_offset=offset,padding=trim(pad),axis=axis)
            else
                call sft%stft_detrend(z,frame,trim(trend),p0=p0,p1=p1, &
                    k_offset=offset,padding=trim(pad),axis=axis)
            end if
        end if
        if (frame%status /= short_fft_success) error stop 'high-rank reference STFT'
        write(uout,'(*(i0,1x))') ic,11,shape(frame%values),frame%p0
        flat = reshape(frame%values,[size(frame%values)])
        do k=1,size(flat)
            write(uout,'(2(es26.17,1x))') real(flat(k),dp), aimag(flat(k))
        end do
    end subroutine run_rank11

    subroutine run_rank12(sft,axis,n,p0,p1,offset,pad,kind,trend,uout,ic)
        type(short_time_fft), intent(in) :: sft !! configured STFT
        integer, intent(in) :: axis,n,p0,p1,offset,uout,ic !! case parameters
        character(*), intent(in) :: pad,kind,trend !! case settings
        type(short_time_fft_result_13d) :: frame
        integer :: dims(12), j, k
        real(dp), allocatable :: x(:,:,:,:,:,:,:,:,:,:,:,:)
        complex(dp), allocatable :: z(:,:,:,:,:,:,:,:,:,:,:,:), flat(:)

        dims = 1
        dims(axis) = n
        j = 1
        if (j == axis) j = 2
        dims(j) = 2
        j = 12
        if (j == axis) j = 11
        dims(j) = 2
        allocate(x(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9),dims(10),dims(11), &
            dims(12)), z(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9),dims(10),dims(11),dims(12)))
        x = reshape([(sin(0.071_dp*real(k,dp)) + &
            0.2_dp*cos(0.13_dp*real(k,dp)),k=1,product(dims))],dims)
        z = cmplx(x,0.15_dp*reshape([(sin(0.11_dp*real(k,dp)), &
            k=1,product(dims))],dims),dp)
        if (kind == 'real') then
            if (trend == 'none') then
                call sft%stft(x,frame,p0=p0,p1=p1,k_offset=offset,padding=trim(pad),axis=axis)
            else
                call sft%stft_detrend(x,frame,trim(trend),p0=p0,p1=p1, &
                    k_offset=offset,padding=trim(pad),axis=axis)
            end if
        else
            if (trend == 'none') then
                call sft%stft(z,frame,p0=p0,p1=p1,k_offset=offset,padding=trim(pad),axis=axis)
            else
                call sft%stft_detrend(z,frame,trim(trend),p0=p0,p1=p1, &
                    k_offset=offset,padding=trim(pad),axis=axis)
            end if
        end if
        if (frame%status /= short_fft_success) error stop 'high-rank reference STFT'
        write(uout,'(*(i0,1x))') ic,12,shape(frame%values),frame%p0
        flat = reshape(frame%values,[size(frame%values)])
        do k=1,size(flat)
            write(uout,'(2(es26.17,1x))') real(flat(k),dp), aimag(flat(k))
        end do
    end subroutine run_rank12

    subroutine run_rank13(sft,axis,n,p0,p1,offset,pad,kind,trend,uout,ic)
        type(short_time_fft), intent(in) :: sft !! configured STFT
        integer, intent(in) :: axis,n,p0,p1,offset,uout,ic !! case parameters
        character(*), intent(in) :: pad,kind,trend !! case settings
        type(short_time_fft_result_14d) :: frame
        integer :: dims(13), j, k
        real(dp), allocatable :: x(:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(dp), allocatable :: z(:,:,:,:,:,:,:,:,:,:,:,:,:), flat(:)

        dims = 1
        dims(axis) = n
        j = 1
        if (j == axis) j = 2
        dims(j) = 2
        j = 13
        if (j == axis) j = 12
        dims(j) = 2
        allocate(x(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9),dims(10),dims(11), &
            dims(12),dims(13)), z(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9), &
            dims(10),dims(11),dims(12),dims(13)))
        x = reshape([(sin(0.071_dp*real(k,dp)) + &
            0.2_dp*cos(0.13_dp*real(k,dp)),k=1,product(dims))],dims)
        z = cmplx(x,0.15_dp*reshape([(sin(0.11_dp*real(k,dp)), &
            k=1,product(dims))],dims),dp)
        if (kind == 'real') then
            if (trend == 'none') then
                call sft%stft(x,frame,p0=p0,p1=p1,k_offset=offset,padding=trim(pad),axis=axis)
            else
                call sft%stft_detrend(x,frame,trim(trend),p0=p0,p1=p1, &
                    k_offset=offset,padding=trim(pad),axis=axis)
            end if
        else
            if (trend == 'none') then
                call sft%stft(z,frame,p0=p0,p1=p1,k_offset=offset,padding=trim(pad),axis=axis)
            else
                call sft%stft_detrend(z,frame,trim(trend),p0=p0,p1=p1, &
                    k_offset=offset,padding=trim(pad),axis=axis)
            end if
        end if
        if (frame%status /= short_fft_success) error stop 'high-rank reference STFT'
        write(uout,'(*(i0,1x))') ic,13,shape(frame%values),frame%p0
        flat = reshape(frame%values,[size(frame%values)])
        do k=1,size(flat)
            write(uout,'(2(es26.17,1x))') real(flat(k),dp), aimag(flat(k))
        end do
    end subroutine run_rank13

    subroutine run_rank14(sft,axis,n,p0,p1,offset,pad,kind,trend,uout,ic)
        type(short_time_fft), intent(in) :: sft !! configured STFT
        integer, intent(in) :: axis,n,p0,p1,offset,uout,ic !! case parameters
        character(*), intent(in) :: pad,kind,trend !! case settings
        type(short_time_fft_result_15d) :: frame
        integer :: dims(14), j, k
        real(dp), allocatable :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(dp), allocatable :: z(:,:,:,:,:,:,:,:,:,:,:,:,:,:), flat(:)

        dims = 1
        dims(axis) = n
        j = 1
        if (j == axis) j = 2
        dims(j) = 2
        j = 14
        if (j == axis) j = 13
        dims(j) = 2
        allocate(x(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9),dims(10),dims(11), &
            dims(12),dims(13),dims(14)), z(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8), &
            dims(9),dims(10),dims(11),dims(12),dims(13),dims(14)))
        x = reshape([(sin(0.071_dp*real(k,dp)) + &
            0.2_dp*cos(0.13_dp*real(k,dp)),k=1,product(dims))],dims)
        z = cmplx(x,0.15_dp*reshape([(sin(0.11_dp*real(k,dp)), &
            k=1,product(dims))],dims),dp)
        if (kind == 'real') then
            if (trend == 'none') then
                call sft%stft(x,frame,p0=p0,p1=p1,k_offset=offset,padding=trim(pad),axis=axis)
            else
                call sft%stft_detrend(x,frame,trim(trend),p0=p0,p1=p1, &
                    k_offset=offset,padding=trim(pad),axis=axis)
            end if
        else
            if (trend == 'none') then
                call sft%stft(z,frame,p0=p0,p1=p1,k_offset=offset,padding=trim(pad),axis=axis)
            else
                call sft%stft_detrend(z,frame,trim(trend),p0=p0,p1=p1, &
                    k_offset=offset,padding=trim(pad),axis=axis)
            end if
        end if
        if (frame%status /= short_fft_success) error stop 'high-rank reference STFT'
        write(uout,'(*(i0,1x))') ic,14,shape(frame%values),frame%p0
        flat = reshape(frame%values,[size(frame%values)])
        do k=1,size(flat)
            write(uout,'(2(es26.17,1x))') real(flat(k),dp), aimag(flat(k))
        end do
    end subroutine run_rank14

end program reference_signal_short_time_fft_high_rank
