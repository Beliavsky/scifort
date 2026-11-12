! SPDX-License-Identifier: MIT
program spline_plan_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_spline_nd, only : spline_interpolation_plan, &
        prepare_spline_plan_nd, evaluate_spline_plan_nd
    implicit none
    type(spline_interpolation_plan) :: plan
    real(dp) :: x(8,7), coordinates(2,96), values(96)
    integer :: i,j,k,point,code
    integer, parameter :: orders(15)=[0,1,3,5,0,1,3,5,3,5,3,5,1,3,5]
    character(len=16), parameter :: modes(15)= &
        [character(len=16) :: 'mirror','mirror','mirror','mirror', &
        'reflect','reflect','nearest','nearest','grid-wrap','grid-wrap', &
        'grid-constant','grid-constant','constant','constant','constant']

    do j=1,7
        do i=1,8
            x(i,j)=sin(0.23_dp*i+0.57_dp*j)+0.1_dp*i*j
        end do
    end do
    coordinates(:,1)=[1.0_dp,1.0_dp]
    coordinates(:,2)=[3.25_dp,2.5_dp]
    coordinates(:,3)=[6.75_dp,5.2_dp]
    coordinates(:,4)=[8.8_dp,7.3_dp]
    coordinates(:,5)=[-0.8_dp,5.0_dp]
    coordinates(:,6)=[3.9_dp,-1.2_dp]
    do point=7,96
        coordinates(1,point)=-2.0_dp+real(mod(point*17,117),dp)/10.0_dp
        coordinates(2,point)=-1.75_dp+real(mod(point*31,108),dp)/10.0_dp
    end do
    do k=1,size(orders)
        call prepare_spline_plan_nd(x,plan,order=orders(k),mode=trim(modes(k)), &
            cval=-2.5_dp,status=code)
        if (code/=0) error stop 'cannot prepare plan'
        call evaluate_spline_plan_nd(plan,coordinates,values,status=code)
        if (code/=0) error stop 'cannot evaluate plan'
        do point=1,size(values)
            print '(i0,1x,i0,1x,es24.16)',k,point,values(point)
        end do
    end do
end program spline_plan_reference_driver
