! SPDX-License-Identifier: MIT
! Values for independent SciPy comparisons of a batched nonlinear warp.
program spline_plan_batch_mapping_driver
    use scifort_kinds, only : dp
    use scifort_signal_spline_nd, only : spline_interpolation_plan, &
        prepare_spline_plan_nd, evaluate_spline_plan_geometric_batch_nd
    implicit none
    character(len=16), parameter :: modes(8)=[character(len=16) :: &
        'mirror','constant','grid-constant','nearest','reflect', &
        'grid-mirror','wrap','grid-wrap']
    type(spline_interpolation_plan) :: plan
    real(dp) :: x(19,17),y(12,11)
    complex(dp) :: z(19,17),w(12,11)
    integer :: i,j,order,mode_index,kind,code

    do j=1,17
        do i=1,19
            x(i,j)=sin(0.17_dp*i+0.31_dp*j)+0.004_dp*i*j
        end do
    end do
    z=cmplx(x,-0.37_dp*x+0.29_dp,dp)
    do mode_index=1,8
        do order=0,5
            do kind=0,1
                if (kind==0) then
                    call prepare_spline_plan_nd(x,plan,order=order, &
                        mode=trim(modes(mode_index)),cval=1.125_dp,status=code)
                    if (code/=0) error stop 'prepare real'
                    call evaluate_spline_plan_geometric_batch_nd(plan,map,y,code,batch_size=13)
                    if (code/=0) error stop 'evaluate real'
                else
                    call prepare_spline_plan_nd(z,plan,order=order, &
                        mode=trim(modes(mode_index)), &
                        cval=cmplx(1.125_dp,-0.35_dp,dp),status=code)
                    if (code/=0) error stop 'prepare complex'
                    call evaluate_spline_plan_geometric_batch_nd(plan,map,w,code,batch_size=13)
                    if (code/=0) error stop 'evaluate complex'
                end if
                do j=1,11
                    do i=1,12
                        if (kind==0) then
                            print '(5(i0,1x),2(es25.16e3,1x))', &
                                mode_index,order,kind,i,j,y(i,j),0.0_dp
                        else
                            print '(5(i0,1x),2(es25.16e3,1x))', &
                                mode_index,order,kind,i,j,real(w(i,j)),aimag(w(i,j))
                        end if
                    end do
                end do
            end do
        end do
    end do
contains
    subroutine map(output_coords,input_coords)
        integer, intent(in) :: output_coords(:,:) !! one-based destination index columns
        real(dp), intent(out) :: input_coords(:,:) !! one-based source locations
        input_coords(1,:)=0.52_dp+0.86_dp*output_coords(1,:) &
            -0.13_dp*output_coords(2,:)+0.11_dp*sin(0.13_dp*output_coords(2,:))
        input_coords(2,:)=0.72_dp+0.22_dp*output_coords(1,:) &
            +0.79_dp*output_coords(2,:)-0.17_dp*cos(0.21_dp*output_coords(1,:))
    end subroutine map
end program spline_plan_batch_mapping_driver
