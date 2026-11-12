! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Simple stdin/stdout numerical-reference protocol for independent SciPy checks.
program reference_signal_istft_nd
    use scifort_signal, only : dp, istft_into
    implicit none
    complex(dp), allocatable :: z(:,:,:), yc(:,:)
    real(dp), allocatable :: yr(:,:), times(:)
    integer :: n1,n2,n3,nout,a,nw,no,nfft,real_flag,boundary_flag,onesided_flag
    integer :: i,j,k,st
    real(dp) :: re,im,fs
    character(32) :: win,scale

    read(*,*) n1,n2,n3,nout,a,nw,no,nfft,real_flag,boundary_flag,onesided_flag,fs
    read(*,*) win,scale
    allocate(z(n1,n2,n3))
    do k=1,n3
        do j=1,n2
            do i=1,n1
                read(*,*) re,im
                z(i,j,k)=cmplx(re,im,dp)
            end do
        end do
    end do
    if (a == 1) then
        allocate(yr(nout,n2),yc(nout,n2))
    else
        allocate(yr(n1,nout),yc(n1,nout))
    end if
    if (real_flag == 1) then
        call istft_into(z,yr,times,fs=fs,window=trim(win),nperseg=nw, &
            noverlap=no,nfft=nfft,boundary=boundary_flag==1,scaling=trim(scale), &
            axis=a,input_onesided=onesided_flag==1,status=st)
        print '(i0)',st
        if (st/=0) stop
        do j=1,size(yr,2)
            do i=1,size(yr,1)
                print '(es26.17,1x,es26.17)',yr(i,j),0.0_dp
            end do
        end do
    else
        call istft_into(z,yc,times,fs=fs,window=trim(win),nperseg=nw, &
            noverlap=no,nfft=nfft,boundary=boundary_flag==1,scaling=trim(scale), &
            axis=a,input_onesided=onesided_flag==1,status=st)
        print '(i0)',st
        if (st/=0) stop
        do j=1,size(yc,2)
            do i=1,size(yc,1)
                print '(es26.17,1x,es26.17)',real(yc(i,j),dp),aimag(yc(i,j))
            end do
        end do
    end if
end program reference_signal_istft_nd
