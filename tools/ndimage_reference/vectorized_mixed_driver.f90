! SPDX-License-Identifier: MIT
program vectorized_mixed_driver
    use scifort_kinds, only : dp
    use scifort_ndimage, only : vectorized_filter_nd, vectorized_filter_multi_nd
    implicit none
    real(dp) :: x(5,4,3),r(5,4,3),rm(5,4,3,2),rv(4,2,2),rmv(4,2,2,2)
    complex(dp) :: z(5,4,3),c(5,4,3),cm(5,4,3,2),cv(4,2,2),cmv(4,2,2,2)
    character(len=16), parameter :: modes(5) = &
        [character(len=16) :: 'reflect','nearest','mirror','wrap','constant']
    integer :: i,j,k,n,unit,code

    do k=1,3
        do j=1,4
            do i=1,5
                x(i,j,k)=sin(0.31_dp*i+0.23_dp*j)+0.17_dp*k+0.09_dp*i*j
            end do
        end do
    end do
    z=cmplx(x,0.37_dp*x+0.21_dp,dp)
    open(newunit=unit,file='vectorized_mixed_reference.bin', &
        form='unformatted',access='stream',status='replace')
    do n=1,size(modes)
        call vectorized_filter_nd(x,rc_scalar,c,status=code,sizes=[2,3,2], &
            mode=trim(modes(n)),cval=-1.75_dp,batch_size=7)
        if (code/=0) error stop 'real to complex scalar'
        write(unit) c
        call vectorized_filter_nd(z,cr_scalar,r,status=code,sizes=[2,3,2], &
            mode=trim(modes(n)),cval=(-1.75_dp,0.6_dp),batch_size=9)
        if (code/=0) error stop 'complex to real scalar'
        write(unit) r
        call vectorized_filter_multi_nd(x,rc_multi,cm,status=code,sizes=[2,3,2], &
            mode=trim(modes(n)),cval=-1.75_dp,batch_size=8)
        if (code/=0) error stop 'real to complex multi'
        write(unit) cm
        call vectorized_filter_multi_nd(z,cr_multi,rm,status=code,sizes=[2,3,2], &
            mode=trim(modes(n)),cval=(-1.75_dp,0.6_dp),batch_size=5)
        if (code/=0) error stop 'complex to real multi'
        write(unit) rm
    end do
    call vectorized_filter_nd(x,rc_scalar,cv,status=code,sizes=[2,3,2],mode='valid')
    if (code/=0) error stop 'real to complex valid'
    write(unit) cv
    call vectorized_filter_nd(z,cr_scalar,rv,status=code,sizes=[2,3,2],mode='valid')
    if (code/=0) error stop 'complex to real valid'
    write(unit) rv
    call vectorized_filter_multi_nd(x,rc_multi,cmv,status=code,sizes=[2,3,2],mode='valid')
    if (code/=0) error stop 'real to complex multi valid'
    write(unit) cmv
    call vectorized_filter_multi_nd(z,cr_multi,rmv,status=code,sizes=[2,3,2],mode='valid')
    if (code/=0) error stop 'complex to real multi valid'
    write(unit) rmv
    call vectorized_filter_nd(x,rc_scalar,c,status=code,sizes=[2,3], &
        axes=[3,1],mode='nearest')
    if (code/=0) error stop 'real to complex axes'
    write(unit) c
    call vectorized_filter_nd(z,cr_scalar,r,status=code,sizes=[2,3], &
        axes=[3,1],mode='nearest')
    if (code/=0) error stop 'complex to real axes'
    write(unit) r
    close(unit)
contains
    subroutine rc_scalar(w,v)
        real(dp), intent(in) :: w(:,:) !! real windows
        complex(dp), intent(out) :: v(:) !! complex mean and mean square
        v=cmplx(sum(w,dim=1),sum(w*w,dim=1),dp)/real(size(w,1),dp)
    end subroutine rc_scalar
    subroutine cr_scalar(w,v)
        complex(dp), intent(in) :: w(:,:) !! complex windows
        real(dp), intent(out) :: v(:) !! mean squared modulus
        v=sum(abs(w)**2,dim=1)/real(size(w,1),dp)
    end subroutine cr_scalar
    subroutine rc_multi(w,v)
        real(dp), intent(in) :: w(:,:) !! real windows
        complex(dp), intent(out) :: v(:,:) !! complex output channels
        v(1,:)=cmplx(sum(w,dim=1),sum(w*w,dim=1),dp)/real(size(w,1),dp)
        v(2,:)=cmplx(sum(w*w,dim=1),sum(w,dim=1),dp)/real(size(w,1),dp)
    end subroutine rc_multi
    subroutine cr_multi(w,v)
        complex(dp), intent(in) :: w(:,:) !! complex windows
        real(dp), intent(out) :: v(:,:) !! real output channels
        v(1,:)=sum(abs(w)**2,dim=1)/real(size(w,1),dp)
        v(2,:)=sum(real(w,dp),dim=1)/real(size(w,1),dp)
    end subroutine cr_multi
end program vectorized_mixed_driver
