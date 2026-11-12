! SPDX-License-Identifier: MIT
! Reusable spline plans for repeated coordinate and grid queries.
program test_signal_spline_plan
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_spline, only : spline_status_success, spline_status_invalid_input
    use scifort_signal_spline_nd, only : spline_interpolation_plan, &
        prepare_spline_plan_nd, evaluate_spline_plan_nd, &
        evaluate_spline_plan_grid_nd, map_coordinates_nd, map_coordinates_grid_nd
    implicit none
    type(spline_interpolation_plan) :: plan
    real(dp) :: x(8,7), cx(2,6), a(6), b(6), old(8,7), &
        grid(2,3,2), y(3,2), strided(6,4), &
        big_coord(2,12), rank15(2,2,1,1,1,1,1,1,1,1,1,1,1,1,1), v15(2)
    complex(dp) :: z(8,7), za(6), zb(6), zy(3,2), zyr(3,2)
    character(len=16) :: modes(5)
    integer :: i,j,k,order,code
    real(dp), parameter :: tol=1.e-11_dp

    do j=1,7
        do i=1,8
            x(i,j)=sin(0.23_dp*i+0.57_dp*j)+0.1_dp*i*j
        end do
    end do
    z=cmplx(x,0.4_dp*x+0.2_dp,dp)
    cx(:,1)=[1.0_dp,1.0_dp]
    cx(:,2)=[3.25_dp,2.5_dp]
    cx(:,3)=[6.75_dp,5.2_dp]
    cx(:,4)=[8.8_dp,7.3_dp]
    cx(:,5)=[-0.8_dp,5.0_dp]
    cx(:,6)=[3.9_dp,-1.2_dp]
    grid=reshape(cx,[2,3,2])
    modes=[character(len=16) :: 'mirror','reflect','nearest','grid-wrap','grid-constant']
    old=x
    do k=1,size(modes)
        do order=0,5
            call prepare_spline_plan_nd(x,plan,order=order,mode=trim(modes(k)), &
                cval=-2.5_dp,status=code)
            call check(code==spline_status_success,'real plan preparation')
            call evaluate_spline_plan_nd(plan,cx,a,status=code)
            call check(code==spline_status_success,'real flat coordinates')
            call map_coordinates_nd(x,cx,b,order=order,mode=trim(modes(k)), &
                cval=-2.5_dp,prefilter=.true.,status=code)
            call check(code==spline_status_success,'direct coordinate reference')
            call check(maxval(abs(a-b))<tol,'real plan matches direct')
            call evaluate_spline_plan_grid_nd(plan,grid,y,status=code)
            call check(code==spline_status_success,'real grid coordinates')
            call check(maxval(abs(y-reshape(a,[3,2])))<tol,'grid equals flat')
        end do
    end do
    ! The plan must own the samples; changes to x cannot change queries.
    call prepare_spline_plan_nd(x,plan,order=3,mode='grid-wrap',status=code)
    call evaluate_spline_plan_nd(plan,cx,a,status=code)
    x=-100.0_dp
    call evaluate_spline_plan_nd(plan,cx,b,status=code)
    call check(maxval(abs(a-b))<tol,'prepared plan owns its data')
    x=old

    call prepare_spline_plan_nd(z,plan,order=5,mode='grid-constant', &
        cval=(-3.0_dp,1.25_dp),status=code)
    call check(code==spline_status_success,'complex plan preparation')
    call evaluate_spline_plan_nd(plan,cx,za,status=code)
    call map_coordinates_nd(z,cx,zb,order=5,mode='grid-constant', &
        cval=(-3.0_dp,1.25_dp),prefilter=.true.,status=code)
    call check(maxval(abs(za-zb))<tol,'complex planar reference')
    call evaluate_spline_plan_grid_nd(plan,grid,zy,status=code)
    call check(code==spline_status_success,'complex grid status')
    call map_coordinates_grid_nd(z,grid,zyr,order=5,mode='grid-constant', &
        cval=(-3.0_dp,1.25_dp),status=code)
    call check(maxval(abs(zy-zyr))<tol,'complex grid reference')
    call evaluate_spline_plan_nd(plan,cx,b,status=code)
    call check(code==spline_status_invalid_input,'real output rejects complex cache')

    call prepare_spline_plan_nd(x,plan,order=3,mode='nearest',status=code)
    call evaluate_spline_plan_grid_nd(plan,grid,strided(1:6:2,1:4:2),status=code)
    call check(code==spline_status_success,'strided plan destination')
    call evaluate_spline_plan_grid_nd(plan,grid,y,status=code)
    call check(maxval(abs(y-strided(1:6:2,1:4:2)))<tol,'strided values')
    big_coord(:,1:12:2)=cx
    big_coord(:,2:12:2)=cx
    call evaluate_spline_plan_nd(plan,big_coord(:,1:12:2),a,status=code)
    call check(code==spline_status_success,'noncontiguous coordinates')
    call evaluate_spline_plan_nd(plan,cx,b,status=code)
    call check(maxval(abs(a-b))<tol,'noncontiguous matches contiguous')

    rank15=1.0_dp
    rank15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=2.0_dp
    rank15(1,2,1,1,1,1,1,1,1,1,1,1,1,1,1)=3.0_dp
    call prepare_spline_plan_nd(rank15,plan,order=1,mode='mirror',status=code)
    block
        real(dp) :: coords15(15,2)
        coords15=1.0_dp
        coords15(1,2)=2.0_dp
        coords15(2,2)=2.0_dp
        call evaluate_spline_plan_nd(plan,coords15,v15,status=code)
    end block
    call check(code==spline_status_success,'rank15 plan status')
    call check(abs(v15(1)-1.0_dp)<tol,'rank15 interpolation')

    ! Many singleton axes must not explode the padded-array allocation.
    block
        real(dp) :: coords15(15,2), original15(2)
        coords15=1.0_dp
        coords15(1,:)=[1.0_dp,2.5_dp]
        coords15(2,:)=[1.0_dp,1.7_dp]
        call prepare_spline_plan_nd(rank15,plan,order=3, &
            mode='grid-constant',cval=-1.25_dp,status=code)
        call check(code==spline_status_success,'rank15 grid-constant plan status')
        call evaluate_spline_plan_nd(plan,coords15,v15,status=code)
        call map_coordinates_nd(rank15,coords15,original15,order=3, &
            mode='grid-constant',prefilter=.true.,cval=-1.25_dp,status=code)
        call check(code==spline_status_success,'rank15 direct grid-constant')
        call check(maxval(abs(v15-original15))<tol,'rank15 padded plan equals direct')
    end block

    call prepare_spline_plan_nd(x,plan,order=7,status=code)
    call check(code==spline_status_invalid_input,'invalid degree')
    call evaluate_spline_plan_nd(plan,cx,a,status=code)
    call check(code==spline_status_invalid_input,'invalid build clears old plan')
    call prepare_spline_plan_nd(x,plan,mode='not-a-mode',status=code)
    call check(code==spline_status_invalid_input,'invalid boundary')
    call prepare_spline_plan_nd(x,plan,status=code)
    call evaluate_spline_plan_nd(plan,cx(1:1,:),a,status=code)
    call check(code==spline_status_invalid_input,'invalid coordinate rank')
    call evaluate_spline_plan_nd(plan,cx,a(1:5),status=code)
    call check(code==spline_status_invalid_input,'invalid output length')
    call evaluate_spline_plan_grid_nd(plan,grid,strided,status=code)
    call check(code==spline_status_invalid_input,'invalid output shape')
    cx(1,1)=ieee_value(0.0_dp,ieee_quiet_nan)
    call evaluate_spline_plan_nd(plan,cx,a,status=code)
    call check(code==spline_status_invalid_input,'nonfinite coordinate')
    print '(a)', 'PASS: reusable real/complex spline plans and grid sampling'
contains
    subroutine check(ok,message)
        logical, intent(in) :: ok !! asserted condition
        character(len=*), intent(in) :: message !! diagnostic text
        if (.not.ok) then
            print '(2a)', 'FAIL: ',message
            error stop 1
        end if
    end subroutine check
end program test_signal_spline_plan
