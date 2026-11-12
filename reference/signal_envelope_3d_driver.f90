! Reference driver, not an FPM test program.
program signal_envelope_3d_driver
    use scifort_kinds, only : dp
    use scifort_signal_envelope, only : envelope
    implicit none
    integer :: nc, k, n1,n2,n3,i,j,h,ax,nout,low,high,cx,square,res,status
    real(dp) :: re,im
    real(dp), allocatable :: x(:,:,:), y(:,:,:), rr(:,:,:)
    complex(dp), allocatable :: z(:,:,:), cr(:,:,:)
    character(7) :: mode
    read(*,*) nc
    do k=1,nc
        read(*,*) n1,n2,n3,ax,nout,low,high,cx,square,res
        allocate(x(n1,n2,n3),z(n1,n2,n3))
        do h=1,n3
            do j=1,n2
                do i=1,n1
                    read(*,*) re,im
                    x(i,j,h)=re
                    z(i,j,h)=cmplx(re,im,dp)
                end do
            end do
        end do
        select case(res)
        case(0)
            mode='none'
        case(1)
            mode='all'
        case default
            mode='lowpass'
        end select
        if(cx==0) then
            y=envelope(x,bp_low=low,bp_high=high,n_out=nout,squared=square==1, &
                residual=trim(mode),axis=ax,remainder=rr,status=status)
            write(*,'(7(i0,1x))') status,shape(y),shape(rr)
            do h=1,size(y,3)
                do j=1,size(y,2)
                    do i=1,size(y,1)
                        re=0.0_dp
                        if(size(rr)>0) re=rr(i,j,h)
                        write(*,'(3(es26.17,1x))') y(i,j,h),re,0.0_dp
                    end do
                end do
            end do
        else
            y=envelope(z,bp_low=low,bp_high=high,n_out=nout,squared=square==1, &
                residual=trim(mode),axis=ax,remainder=cr,status=status)
            write(*,'(7(i0,1x))') status,shape(y),shape(cr)
            do h=1,size(y,3)
                do j=1,size(y,2)
                    do i=1,size(y,1)
                        re=0.0_dp
                        im=0.0_dp
                        if(size(cr)>0) then
                            re=real(cr(i,j,h),dp)
                            im=aimag(cr(i,j,h))
                        end if
                        write(*,'(3(es26.17,1x))') y(i,j,h),re,im
                    end do
                end do
            end do
        end if
        deallocate(x,z,y)
        if(allocated(rr)) deallocate(rr)
        if(allocated(cr)) deallocate(cr)
    end do
end program signal_envelope_3d_driver
