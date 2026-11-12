! SPDX-License-Identifier: MIT
! Streaming probe for rank-three real and complex filtering.
program lfilter_rank3_driver
    use scifort_kinds, only : dp
    use scifort_signal_lfilter, only : lfilter
    implicit none
    integer :: kind,n(3),nb,na,axis,with_state,st,io,sh(3),i,j,k
    real(dp) :: re,im
    complex(dp), allocatable :: a(:),b(:),x(:,:,:),zi(:,:,:),y(:,:,:),zf(:,:,:)
    real(dp), allocatable :: ry(:,:,:),rzf(:,:,:)

    do
        read(*,*,iostat=io) kind,n,nb,na,axis,with_state
        if (io /= 0) exit
        allocate(a(na),b(nb),x(n(1),n(2),n(3)))
        do i=1,nb
            read(*,*) re,im
            b(i) = cmplx(re,im,dp)
        end do
        do i=1,na
            read(*,*) re,im
            a(i) = cmplx(re,im,dp)
        end do
        do k=1,n(3)
            do j=1,n(2)
                do i=1,n(1)
                    read(*,*) re,im
                    x(i,j,k) = cmplx(re,im,dp)
                end do
            end do
        end do
        sh=n
        sh(axis)=max(na,nb)-1
        if (with_state==1) then
            allocate(zi(sh(1),sh(2),sh(3)))
            do k=1,sh(3)
                do j=1,sh(2)
                    do i=1,sh(1)
                        read(*,*) re,im
                        zi(i,j,k) = cmplx(re,im,dp)
                    end do
                end do
            end do
        end if
        if (kind == 0) then
            if (with_state == 1) then
                ry = lfilter(real(b,dp),real(a,dp),real(x,dp), &
                    zi=real(zi,dp),zf=rzf,status=st,axis=axis)
            else
                ry = lfilter(real(b,dp),real(a,dp),real(x,dp), &
                    zf=rzf,status=st,axis=axis)
            end if
            y = cmplx(ry,0.0_dp,dp)
            zf = cmplx(rzf,0.0_dp,dp)
            deallocate(ry,rzf)
        else
            if (with_state == 1) then
                y = lfilter(b,a,x,zi=zi,zf=zf,status=st,axis=axis)
            else
                y = lfilter(b,a,x,zf=zf,status=st,axis=axis)
            end if
        end if
        print '(7(i0,1x))',st,shape(y),shape(zf)
        do k=1,size(y,3)
            do j=1,size(y,2)
                do i=1,size(y,1)
                    print '(es25.17,1x,es25.17)',real(y(i,j,k),dp),aimag(y(i,j,k))
                end do
            end do
        end do
        do k=1,size(zf,3)
            do j=1,size(zf,2)
                do i=1,size(zf,1)
                    print '(es25.17,1x,es25.17)',real(zf(i,j,k),dp),aimag(zf(i,j,k))
                end do
            end do
        end do
        deallocate(a,b,x,y,zf)
        if (allocated(zi)) deallocate(zi)
    end do
end program lfilter_rank3_driver
