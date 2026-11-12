! SPDX-License-Identifier: MIT
! Streaming probe for SciFort rank-two lfilter and final DF-II states.
program lfilter_axes_driver
    use scifort_kinds, only : dp
    use scifort_signal_lfilter, only : lfilter
    implicit none
    integer :: kind, nx, ny, nb, na, axis, with_state, st, io
    integer :: i,j
    real(dp) :: re, im
    complex(dp), allocatable :: b(:),a(:),x(:,:),zi(:,:),y(:,:),zf(:,:)

    do
        read(*,*,iostat=io) kind,nx,ny,nb,na,axis,with_state
        if (io /= 0) exit
        allocate(b(nb),a(na),x(nx,ny))
        do i = 1, nb
            read(*,*) re,im
            b(i) = cmplx(re,im,dp)
        end do
        do i = 1, na
            read(*,*) re,im
            a(i) = cmplx(re,im,dp)
        end do
        do j = 1, ny
            do i = 1, nx
                read(*,*) re,im
                x(i,j) = cmplx(re,im,dp)
            end do
        end do
        if (with_state == 1) then
            if (axis == 1) then
                allocate(zi(max(na,nb)-1,ny))
            else
                allocate(zi(nx,max(na,nb)-1))
            end if
            do j = 1, size(zi,2)
                do i = 1, size(zi,1)
                    read(*,*) re,im
                    zi(i,j) = cmplx(re,im,dp)
                end do
            end do
        end if
        if (kind == 0) then
            if (with_state == 1) then
                call run_real(b,a,x,zi,axis,st,y,zf)
            else
                call run_real(b,a,x,axis=axis,st=st,y=y,zf=zf)
            end if
        else
            if (with_state == 1) then
                y = lfilter(b,a,x,zi=zi,zf=zf,status=st,axis=axis)
            else
                y = lfilter(b,a,x,zf=zf,status=st,axis=axis)
            end if
        end if
        print '(5(i0,1x))', st,size(y,1),size(y,2),size(zf,1),size(zf,2)
        do j = 1, size(y,2)
            do i = 1, size(y,1)
                print '(es25.17,1x,es25.17)',real(y(i,j),dp),aimag(y(i,j))
            end do
        end do
        do j = 1, size(zf,2)
            do i = 1, size(zf,1)
                print '(es25.17,1x,es25.17)',real(zf(i,j),dp),aimag(zf(i,j))
            end do
        end do
        deallocate(a,b,x,y,zf)
        if (allocated(zi)) deallocate(zi)
    end do
contains
    subroutine run_real(b,a,x,zi,axis,st,y,zf)
        complex(dp), intent(in) :: b(:) !! numerator
        complex(dp), intent(in) :: a(:) !! denominator
        complex(dp), intent(in) :: x(:,:) !! input samples
        complex(dp), intent(in), optional :: zi(:,:) !! initial states
        integer, intent(in) :: axis !! dimension to filter
        integer, intent(out) :: st !! status code
        complex(dp), allocatable, intent(out) :: y(:,:) !! filtered signals
        complex(dp), allocatable, intent(out) :: zf(:,:) !! final states
        real(dp), allocatable :: yr(:,:),zr(:,:)

        if (present(zi)) then
            yr = lfilter(real(b,dp),real(a,dp),real(x,dp),zi=real(zi,dp), &
                zf=zr,status=st,axis=axis)
        else
            yr = lfilter(real(b,dp),real(a,dp),real(x,dp),zf=zr,status=st,axis=axis)
        end if
        y = cmplx(yr,0.0_dp,dp)
        zf = cmplx(zr,0.0_dp,dp)
    end subroutine run_real
end program lfilter_axes_driver
