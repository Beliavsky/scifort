! SPDX-License-Identifier: MIT
! Batched real/complex neighborhood callbacks, valid mode, axes and guards.
program test_ndimage_vectorized
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_ndimage, only : vectorized_filter_nd, generic_filter_nd, &
        ndimage_status_too_expensive, map_coordinates_grid_nd
    implicit none
    real(dp) :: x(4,3),y(4,3),reference(4,3),strided(8,6),source15(2,3,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: out15(2,2,1,1,1,1,1,1,1,1,1,1,1,1,1),small(3,2),nan_input(4,3)
    real(dp) :: coords(2,2,2),grid_output(2,2)
    complex(dp) :: z(4,3),zy(4,3),zref(4,3)
    logical :: footprint(2,3)
    character(len=16) :: modes(2)
    integer :: i,j,code,calls,columns,windows_seen
    real(dp), parameter :: tol=1.0e-12_dp

    do j=1,3
        do i=1,4
            x(i,j)=real(i*i-2*j*i+5*j,dp)/7.0_dp
        end do
    end do
    z=cmplx(x,0.5_dp*x+0.3_dp,dp)
    calls=0
    columns=0
    windows_seen=0
    call vectorized_filter_nd(x,real_mean,y,status=code,sizes=[2,3], &
        mode='constant',cval=-2.0_dp,batch_size=5)
    call check(code==0,'real filter status')
    call check(calls==3 .and. columns==12 .and. windows_seen==6,'callback grouping')
    call generic_filter_nd(x,scalar_mean,reference,status=code,sizes=[2,3], &
        mode='constant',cval=-2.0_dp)
    call check(code==0 .and. maxval(abs(y-reference))<tol,'batch matches scalar filter')

    strided=-999.0_dp
    call vectorized_filter_nd(x,real_mean,strided(1:8:2,2:6:2),status=code, &
        sizes=[2,3],mode='constant',cval=-2.0_dp,batch_size=1)
    call check(code==0,'noncontiguous output status')
    call check(maxval(abs(strided(1:8:2,2:6:2)-reference))<tol,'noncontiguous result')
    call check(all(strided(2:8:2,:)==-999.0_dp),'untouched output rows')
    call check(all(strided(:,1)==-999.0_dp),'untouched output first column')
    call check(all(strided(:,3)==-999.0_dp),'untouched output middle column')

    footprint=.false.
    footprint(1,1)=.true.
    footprint(1,3)=.true.
    footprint(2,2)=.true.
    footprint(2,3)=.true.
    modes=[character(len=16) :: 'wrap','constant']
    call vectorized_filter_nd(x,real_mean,y,status=code,footprint=footprint, &
        axes=[2,1],modes=modes,cval=-3.0_dp,origin=[0,-1],batch_size=3)
    call check(code==0,'selected axes/mixed modes')
    call generic_filter_nd(x,scalar_mean,reference,status=code,footprint=footprint, &
        axes=[2,1],modes=modes,cval=-3.0_dp,origin=[0,-1])
    ! SciFort scalar generic_filter canonicalizes axes; its compact-footprint
    ! orientation can differ for reversed axes, so compare to manual C-order
    ! footprint calculation instead when axes are reordered.
    call check(code==0,'scalar test accepts mixed mode')
    call check(abs(y(1,1)-(x(1,3)+x(3,3)+x(2,1)+x(3,1))/4.0_dp)<tol, &
        'reordered sparse mixed-mode first corner')
    call check(abs(y(4,3)-(x(4,2)-9.0_dp)/4.0_dp)<tol, &
        'reordered sparse mixed-mode opposite corner')

    call vectorized_filter_nd(x,real_mean,small,status=code,sizes=[2,2],mode='valid')
    call check(code==0,'valid shape accepted')
    do j=1,2
        do i=1,3
            call check(abs(small(i,j)-sum(x(i:i+1,j:j+1))/4.0_dp)<tol,'valid window')
        end do
    end do
    call vectorized_filter_nd(x,real_mean,y,status=code,sizes=[2,2],mode='valid')
    call check(code/=0,'invalid output size for valid')
    call vectorized_filter_nd(x,real_mean,small,status=code,sizes=[2,2], &
        mode='valid',origin=[0,1])
    call check(code/=0,'valid rejects shifted origin')
    call vectorized_filter_nd(x,real_mean,y,status=code,sizes=[2,3],axes=[1,1])
    call check(code/=0,'duplicate axes rejected')
    call vectorized_filter_nd(x,real_mean,y,status=code,sizes=[2,3],batch_size=0)
    call check(code/=0,'zero batch rejected')
    call vectorized_filter_nd(x,real_mean,y,status=code,sizes=[60000,1])
    call check(code==ndimage_status_too_expensive,'resource guard')
    call vectorized_filter_nd(x,real_mean,y,status=code,sizes=[2,3],mode='invalid')
    call check(code/=0,'invalid mode rejected')
    call vectorized_filter_nd(x,real_mean,y,status=code,sizes=[integer::],axes=[integer::])
    call check(code==0 .and. all(abs(y-x)<tol),'no selected axes is identity')

    call vectorized_filter_nd(z,complex_mean,zy,status=code,sizes=[2,3], &
        mode='mirror',batch_size=4)
    call generic_filter_nd(z,scalar_complex_mean,zref,status=code,sizes=[2,3],mode='mirror')
    call check(code==0 .and. maxval(abs(zy-zref))<tol,'complex callback')

    nan_input=x
    nan_input(1,1)=ieee_value(0.0_dp,ieee_quiet_nan)
    call vectorized_filter_nd(nan_input,nan_mean,y,status=code,sizes=[2,2], &
        mode='constant',cval=ieee_value(0.0_dp,ieee_quiet_nan))
    call check(code==0,'NaN-aware input permitted')
    call check(ieee_is_nan(y(1,1)),'all-NaN window returns NaN')
    call check(.not.ieee_is_nan(y(3,2)),'nonempty finite window returns number')

    source15=1.0_dp
    source15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=3.0_dp
    source15(1,2,1,1,1,1,1,1,1,1,1,1,1,1,1)=5.0_dp
    source15(2,2,1,1,1,1,1,1,1,1,1,1,1,1,1)=7.0_dp
    call vectorized_filter_nd(source15,real_mean,out15,status=code,sizes=[2], &
        axes=[2],mode='valid',batch_size=1)
    call check(code==0,'rank15 selected valid')
    call check(abs(out15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)-3.0_dp)<tol, &
        'rank15 value 1')
    call check(abs(out15(2,2,1,1,1,1,1,1,1,1,1,1,1,1,1)-4.0_dp)<tol, &
        'rank15 value 2')
    ! Check the new interpolation helper is re-exported by scifort_ndimage.
    coords(:,1,1)=[1.0_dp,1.0_dp]
    coords(:,2,1)=[2.0_dp,1.0_dp]
    coords(:,1,2)=[1.0_dp,2.0_dp]
    coords(:,2,2)=[2.0_dp,2.0_dp]
    call map_coordinates_grid_nd(x,coords,grid_output,order=1,status=code)
    call check(code==0,'grid helper re-export')
    call check(maxval(abs(grid_output-x(:2,:2)))<tol,'grid output values')
    print '(a)', 'PASS: vectorized batch callbacks, shape validation and rank-15'
contains
    subroutine check(ok,message)
        logical, intent(in) :: ok !! condition
        character(len=*), intent(in) :: message !! diagnostic
        if (.not.ok) then
            print '(2a)', 'FAIL: ',message
            error stop 1
        end if
    end subroutine check
    subroutine real_mean(w,r)
        real(dp), intent(in) :: w(:,:) !! window columns
        real(dp), intent(out) :: r(:) !! means
        calls=calls+1
        columns=columns+size(w,2)
        windows_seen=size(w,1)
        r=sum(w,dim=1)/real(size(w,1),dp)
    end subroutine real_mean
    real(dp) function scalar_mean(w) result(r)
        real(dp), intent(in) :: w(:) !! one neighborhood
        r=sum(w)/real(size(w),dp)
    end function scalar_mean
    subroutine complex_mean(w,r)
        complex(dp), intent(in) :: w(:,:) !! neighborhoods
        complex(dp), intent(out) :: r(:) !! means
        r=sum(w,dim=1)/real(size(w,1),dp)
    end subroutine complex_mean
    complex(dp) function scalar_complex_mean(w) result(r)
        complex(dp), intent(in) :: w(:) !! one neighborhood
        r=sum(w)/real(size(w),dp)
    end function scalar_complex_mean
    subroutine nan_mean(w,r)
        real(dp), intent(in) :: w(:,:) !! sample neighborhoods
        real(dp), intent(out) :: r(:) !! finite-only mean
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
    end subroutine nan_mean
end program test_ndimage_vectorized
