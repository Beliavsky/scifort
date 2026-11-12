! SPDX-License-Identifier: MIT
! Generate binary reference cases for scipy.ndimage.vectorized_filter.
program vectorized_filter_driver
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_ndimage, only : vectorized_filter_nd
    implicit none
    real(dp) :: x(5,4,3),y(5,4,3),nan_x(5,4,3),valid_y(4,2,2),selected_valid(3,4,2)
    complex(dp) :: z(5,4,3),zy(5,4,3),valid_z(4,2,2)
    logical :: footprint(2,3)
    character(len=16), parameter :: modes(5)=[character(len=16) :: &
        'reflect','nearest','mirror','wrap','constant']
    integer :: i,j,k,m,code

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
    open(newunit=m,file='vectorized_reference.bin',access='stream',form='unformatted', &
         status='replace',action='write')
    do k=1,size(modes)
        call vectorized_filter_nd(x,mean_real,y,status=code,sizes=[2,3,2], &
            mode=trim(modes(k)),cval=-1.75_dp,batch_size=7)
        call check(code==0,'real all modes')
        write(m) y
        call vectorized_filter_nd(z,mean_complex,zy,status=code,sizes=[2,3,2], &
            mode=trim(modes(k)),cval=cmplx(-1.75_dp,0.6_dp,dp),batch_size=11)
        call check(code==0,'complex all modes')
        write(m) zy
    end do
    call vectorized_filter_nd(x,weighted_real,y,status=code,footprint=footprint, &
        axes=[3,1],origin=[0,1],mode='nearest',batch_size=3)
    call check(code==0,'sparse reordered footprint')
    write(m) y
    call vectorized_filter_nd(x,mean_real,valid_y,status=code,sizes=[2,3,2], &
        mode='valid',batch_size=5)
    call check(code==0,'valid real')
    write(m) valid_y
    call vectorized_filter_nd(z,mean_complex,valid_z,status=code,sizes=[2,3,2], &
        mode='valid',batch_size=5)
    call check(code==0,'valid complex')
    write(m) valid_z
    call vectorized_filter_nd(x,weighted_real,selected_valid,status=code, &
        footprint=footprint,axes=[3,1],mode='valid',batch_size=2)
    call check(code==0,'selected-axis sparse valid')
    write(m) selected_valid
    call vectorized_filter_nd(x,mean_real,y,status=code,sizes=[2,3,2], &
        origin=[-1,1,0],mode='mirror',batch_size=2)
    call check(code==0,'shifted origins')
    write(m) y
    nan_x=x
    nan_x(2,3,1)=ieee_value(0.0_dp,ieee_quiet_nan)
    call vectorized_filter_nd(nan_x,nanmean_real,y,status=code,sizes=[2,2,3], &
        mode='constant',cval=ieee_value(0.0_dp,ieee_quiet_nan),batch_size=4)
    call check(code==0,'nanmean and NaN boundary value')
    write(m) y
    call vectorized_filter_nd(x,mean_real,y,status=code,sizes=[integer::],axes=[integer::])
    call check(code==0,'empty axis selection')
    write(m) y
    close(m)
    print '(a)', 'PASS: vectorized SciPy reference driver'
contains
    subroutine check(ok,msg)
        logical, intent(in) :: ok !! assertion
        character(len=*), intent(in) :: msg !! failure description
        if (.not.ok) then
            print '(2a)', 'FAIL: ',msg
            error stop 1
        end if
    end subroutine
    subroutine mean_real(w,r)
        real(dp), intent(in) :: w(:,:) !! windows
        real(dp), intent(out) :: r(:) !! mean of each window
        r=sum(w,dim=1)/real(size(w,1),dp)
    end subroutine
    subroutine mean_complex(w,r)
        complex(dp), intent(in) :: w(:,:) !! windows
        complex(dp), intent(out) :: r(:) !! mean of each window
        r=sum(w,dim=1)/real(size(w,1),dp)
    end subroutine
    subroutine weighted_real(w,r)
        real(dp), intent(in) :: w(:,:) !! sparse window in C-order
        real(dp), intent(out) :: r(:) !! order-sensitive weighted response
        integer :: c
        do c=1,size(r)
            r(c)=dot_product([1.0_dp,2.0_dp,3.0_dp,4.0_dp],w(:,c))
        end do
    end subroutine
    subroutine nanmean_real(w,r)
        real(dp), intent(in) :: w(:,:) !! input may contain NaN
        real(dp), intent(out) :: r(:) !! means ignoring NaNs
        integer :: c,k,n
        real(dp) :: total
        do c=1,size(r)
            total=0.0_dp
            n=0
            do k=1,size(w,1)
                if (ieee_is_nan(w(k,c))) cycle
                total=total+w(k,c)
                n=n+1
            end do
            r(c)=ieee_value(0.0_dp,ieee_quiet_nan)
            if (n>0) r(c)=total/real(n,dp)
        end do
    end subroutine
end program vectorized_filter_driver
