! SPDX-License-Identifier: MIT
! Print real/complex spline-plan displacement results for SciPy comparison.
program spline_plan_displacements_driver
    use scifort_kinds, only : dp
    use scifort_signal_spline_nd, only : spline_interpolation_plan, &
        prepare_spline_plan_nd, evaluate_spline_plan_displacement_nd, &
        evaluate_spline_plan_displacement_grid_nd
    implicit none
    character(len=16), parameter :: modes(8)=[character(len=16) :: &
        'mirror','constant','grid-constant','nearest','reflect', &
        'grid-mirror','wrap','grid-wrap']
    type(spline_interpolation_plan) :: plan
    real(dp) :: x(16,13),y(14,12),delta(2,14,12),line(27),line_offsets(2,27)
    real(dp) :: matrix(2,1),offset(2)
    complex(dp) :: z(16,13),zy(14,12)
    integer :: i,j,m,degree,code

    do j=1,13
        do i=1,16
            x(i,j)=sin(0.17_dp*i+0.29_dp*j)+0.002_dp*i*j
        end do
    end do
    z=cmplx(x,0.32_dp*x+0.1_dp,dp)
    do j=1,12
        do i=1,14
            delta(1,i,j)=0.21_dp*sin(0.4_dp*(i+j))
            delta(2,i,j)=-0.31_dp*cos(0.27_dp*(2*i-j))
        end do
    end do
    do m=1,8
        do degree=0,5
            call prepare_spline_plan_nd(x,plan,order=degree,mode=trim(modes(m)), &
                cval=-0.4_dp,status=code)
            if (code/=0) error stop 'real plan failure'
            call evaluate_spline_plan_displacement_grid_nd(plan,delta,y,code,batch_size=11)
            if (code/=0) error stop 'real warp failure'
            call prepare_spline_plan_nd(z,plan,order=degree,mode=trim(modes(m)), &
                cval=cmplx(-0.4_dp,0.25_dp,dp),status=code)
            if (code/=0) error stop 'complex plan failure'
            call evaluate_spline_plan_displacement_grid_nd(plan,delta,zy,code,batch_size=17)
            if (code/=0) error stop 'complex warp failure'
            do j=1,12
                do i=1,14
                    write (*,'(5(i0,1x),2(es25.16e3,1x))') m,degree,0,i,j,y(i,j),0.0_dp
                    write (*,'(5(i0,1x),2(es25.16e3,1x))') m,degree,1,i,j,real(zy(i,j),dp),aimag(zy(i,j))
                end do
            end do
        end do
    end do

    matrix(:,1)=[0.45_dp,0.27_dp]
    offset=[0.1_dp,-0.3_dp]
    do i=1,27
        line_offsets(:,i)=[0.1_dp*sin(0.31_dp*i),0.15_dp*cos(0.23_dp*i)]
    end do
    call prepare_spline_plan_nd(x,plan,order=3,mode='grid-wrap',status=code)
    if (code/=0) error stop 'rank-change plan failure'
    call evaluate_spline_plan_displacement_nd(plan,line_offsets,line,code, &
        matrix=matrix,offset=offset,batch_size=5)
    if (code/=0) error stop 'rank-change warp failure'
    do i=1,27
        write (*,'(5(i0,1x),2(es25.16e3,1x))') 0,3,0,i,1,line(i),0.0_dp
    end do
end program spline_plan_displacements_driver
