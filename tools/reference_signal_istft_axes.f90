! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Reference protocol for independent SciPy inverse-STFT axis permutations.
program reference_signal_istft_axes
    use scifort_signal, only : dp, istft_into
    implicit none
    complex(dp), allocatable :: z(:,:,:), yc(:,:)
    real(dp), allocatable :: yr(:,:), times(:)
    integer :: n1,n2,n3,nout,fa,ta,nw,no,nfft,real_flag,boundary_flag,onesided_flag,mode
    integer :: dims(3), dout(2), a, i, j, k, st
    real(dp) :: re,im,fs
    character(32) :: win,scale

    read(*,*) n1,n2,n3,nout,fa,ta,nw,no,nfft,real_flag,boundary_flag,onesided_flag,fs,mode
    read(*,*) win,scale
    dims = [n1,n2,n3]
    a = fa
    if (fa > ta) a = fa-1
    dout = pack(dims, [ta /= 1,ta /= 2,ta /= 3])
    dout(a) = nout
    allocate(z(n1,n2,n3))
    do k=1,n3
        do j=1,n2
            do i=1,n1
                read(*,*) re,im
                z(i,j,k)=cmplx(re,im,dp)
            end do
        end do
    end do
    if (real_flag == 1) then
        allocate(yr(dout(1),dout(2)))
        yr=-999.0_dp
        select case(mode)
        case(0)
            call istft_into(z,yr,times,fs=fs,window=trim(win),nperseg=nw, &
                noverlap=no,nfft=nfft,boundary=boundary_flag==1,scaling=trim(scale), &
                freq_axis=fa,time_axis=ta,input_onesided=onesided_flag==1,status=st)
        case(1)
            call istft_into(z,yr,times,fs=fs,window=trim(win),nperseg=nw, &
                noverlap=no,nfft=nfft,boundary=boundary_flag==1,scaling=trim(scale), &
                axis=a,time_axis=ta,input_onesided=onesided_flag==1,status=st)
        case default
            call istft_into(z,yr,times,fs=fs,window=trim(win),nperseg=nw, &
                noverlap=no,nfft=nfft,boundary=boundary_flag==1,scaling=trim(scale), &
                axis=a,freq_axis=fa,time_axis=ta,input_onesided=onesided_flag==1,status=st)
        end select
        print '(i0)',st
        if (st/=0) stop
        do i=1,size(times)
            print '(es26.17)',times(i)
        end do
        do j=1,size(yr,2)
            do i=1,size(yr,1)
                print '(es26.17,1x,es26.17)',yr(i,j),0.0_dp
            end do
        end do
    else
        allocate(yc(dout(1),dout(2)))
        yc=(-999.0_dp,0.0_dp)
        select case(mode)
        case(0)
            call istft_into(z,yc,times,fs=fs,window=trim(win),nperseg=nw, &
                noverlap=no,nfft=nfft,boundary=boundary_flag==1,scaling=trim(scale), &
                freq_axis=fa,time_axis=ta,input_onesided=onesided_flag==1,status=st)
        case(1)
            call istft_into(z,yc,times,fs=fs,window=trim(win),nperseg=nw, &
                noverlap=no,nfft=nfft,boundary=boundary_flag==1,scaling=trim(scale), &
                axis=a,time_axis=ta,input_onesided=onesided_flag==1,status=st)
        case default
            call istft_into(z,yc,times,fs=fs,window=trim(win),nperseg=nw, &
                noverlap=no,nfft=nfft,boundary=boundary_flag==1,scaling=trim(scale), &
                axis=a,freq_axis=fa,time_axis=ta,input_onesided=onesided_flag==1,status=st)
        end select
        print '(i0)',st
        if (st/=0) stop
        do i=1,size(times)
            print '(es26.17)',times(i)
        end do
        do j=1,size(yc,2)
            do i=1,size(yc,1)
                print '(es26.17,1x,es26.17)',real(yc(i,j),dp),aimag(yc(i,j))
            end do
        end do
    end if
end program reference_signal_istft_axes
