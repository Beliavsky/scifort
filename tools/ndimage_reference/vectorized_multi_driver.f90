! SPDX-License-Identifier: MIT
! Emit reference arrays for real/complex multi-result vectorized filters.
program vectorized_multi_driver
    use scifort_kinds, only : dp
    use scifort_ndimage, only : vectorized_filter_multi_nd, vectorized_filter_nd
    implicit none
    real(dp) :: x(5,4,3),y(5,4,3,2),y_flat(60,2),valid_y(4,2,2,2),select_valid(3,4,2,2)
    real(dp) :: scalar(5,4,3),strided(10,8,3,2)
    complex(dp) :: z(5,4,3),zy(5,4,3,2),valid_z(4,2,2,2)
    logical :: footprint(2,3)
    character(len=16), parameter :: modes(5)=[character(len=16) :: &
        'reflect','nearest','mirror','wrap','constant']
    integer :: i,j,k,u,code

    do k=1,3
        do j=1,4
            do i=1,5
                x(i,j,k)=sin(0.31_dp*i+0.23_dp*j)+0.17_dp*k+0.09_dp*i*j
            end do
        end do
    end do
    z=cmplx(x,0.37_dp*x+0.21_dp,dp)
    footprint=.false.
    footprint(1,1)=.true.
    footprint(1,3)=.true.
    footprint(2,2)=.true.
    footprint(2,3)=.true.
    open(newunit=u,file='vectorized_multi_reference.bin',access='stream',form='unformatted', &
         status='replace',action='write')
    do k=1,size(modes)
        call vectorized_filter_multi_nd(x,real_features,y,status=code,sizes=[2,3,2], &
            mode=trim(modes(k)),cval=-1.75_dp,batch_size=7)
        call check(code==0,'real all modes')
        write(u) y
        call vectorized_filter_multi_nd(z,complex_features,zy,status=code,sizes=[2,3,2], &
            mode=trim(modes(k)),cval=cmplx(-1.75_dp,0.6_dp,dp),batch_size=11)
        call check(code==0,'complex all modes')
        write(u) zy
    end do
    call vectorized_filter_multi_nd(x,real_features,y,status=code,footprint=footprint, &
        axes=[3,1],origin=[0,1],mode='nearest',batch_size=3)
    call check(code==0,'sparse reordered axes')
    write(u) y
    call vectorized_filter_multi_nd(x,real_features,valid_y,status=code,sizes=[2,3,2], &
        mode='valid',batch_size=5)
    call check(code==0,'valid real')
    write(u) valid_y
    call vectorized_filter_multi_nd(z,complex_features,valid_z,status=code,sizes=[2,3,2], &
        mode='valid',batch_size=5)
    call check(code==0,'valid complex')
    write(u) valid_z
    call vectorized_filter_multi_nd(x,real_features,select_valid,status=code, &
        footprint=footprint,axes=[3,1],mode='valid',batch_size=2)
    call check(code==0,'selected valid')
    write(u) select_valid
    call vectorized_filter_multi_nd(x,real_features,y,status=code,sizes=[2,3,2], &
        origin=[-1,1,0],mode='mirror',batch_size=2)
    call check(code==0,'shifted origins')
    write(u) y
    call vectorized_filter_multi_nd(x,real_features,y_flat,status=code,sizes=[2,3,2], &
        mode='nearest',batch_size=4)
    call check(code==0,'flattened output')
    write(u) y_flat
    strided=-99.0_dp
    call vectorized_filter_multi_nd(x,real_features,strided(1:10:2,2:8:2,:,:),status=code, &
        sizes=[2,3,2],mode='nearest',batch_size=4)
    call check(code==0,'noncontiguous output')
    write(u) strided(1:10:2,2:8:2,:,:)
    call vectorized_filter_nd(x,scalar_mean,scalar,status=code,sizes=[2,3,2],mode='nearest')
    call check(code==0,'original scalar batch fast path')
    write(u) scalar
    close(u)
    print '(a)', 'PASS: vectorized multi-output SciPy reference driver'
contains
    subroutine check(ok,msg)
        logical, intent(in) :: ok !! assertion
        character(len=*), intent(in) :: msg !! message
        if (.not.ok) then
            print '(2a)', 'FAIL: ',msg
            error stop 1
        end if
    end subroutine check
    subroutine real_features(w,r)
        real(dp), intent(in) :: w(:,:) !! windows (neighbors,batch)
        real(dp), intent(out) :: r(:,:) !! channels by batch
        r(1,:)=sum(w,dim=1)/real(size(w,1),dp)
        r(2,:)=sum(w*w,dim=1)/real(size(w,1),dp)
    end subroutine real_features
    subroutine complex_features(w,r)
        complex(dp), intent(in) :: w(:,:) !! complex neighborhoods
        complex(dp), intent(out) :: r(:,:) !! complex channels x batch
        r(1,:)=sum(w,dim=1)/real(size(w,1),dp)
        r(2,:)=sum(w*w,dim=1)/real(size(w,1),dp)
    end subroutine complex_features
    subroutine scalar_mean(w,r)
        real(dp), intent(in) :: w(:,:) !! scalar neighborhoods
        real(dp), intent(out) :: r(:) !! one result per pixel
        r=sum(w,dim=1)/real(size(w,1),dp)
    end subroutine scalar_mean
end program vectorized_multi_driver
