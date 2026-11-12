! SPDX-License-Identifier: MIT
! Emit tensor-grid results for an independent SciPy comparison.
program spline_tensor_driver
    use scifort_kinds, only : dp
    use scifort_signal_spline_nd, only : spline_interpolation_plan, spline_axis_grid, &
        prepare_spline_plan_nd, evaluate_spline_plan_tensor_nd
    implicit none
    character(len=16), parameter :: modes(8)=[character(len=16) :: &
        'mirror','constant','grid-constant','nearest','reflect','grid-mirror','wrap','grid-wrap']
    real(dp) :: x(11,9)
    complex(dp) :: cx(11,9), result(13,9)
    type(spline_axis_grid) :: grid(2)
    type(spline_interpolation_plan) :: plan
    integer :: i,j,degree,modeid,filtered,code

    do j=1,9
        do i=1,11
            x(i,j)=sin(0.13_dp*i+0.27_dp*j)+0.001_dp*i*j
        end do
    end do
    cx=cmplx(x,-0.37_dp*x+0.2_dp,dp)
    allocate(grid(1)%coordinates(13),grid(2)%coordinates(9))
    do i=1,13
        grid(1)%coordinates(i)=real(i-2,dp)+0.31_dp
    end do
    do j=1,9
        grid(2)%coordinates(j)=real(j-2,dp)+0.27_dp
    end do
    do modeid=1,size(modes)
        do degree=0,5
            do filtered=0,1
                call prepare_spline_plan_nd(cx,plan,degree,code,trim(modes(modeid)), &
                    filtered==1,cmplx(1.3_dp,-0.5_dp,dp))
                if (code/=0) error stop 'preparation failed'
                call evaluate_spline_plan_tensor_nd(plan,grid,result,code)
                if (code/=0) error stop 'tensor interpolation failed'
                do j=1,9
                    do i=1,13
                        write(*,'(5(i0,1x),2(es25.16e3,1x))') modeid,degree,filtered,i,j, &
                            real(result(i,j),dp),aimag(result(i,j))
                    end do
                end do
            end do
        end do
    end do
end program spline_tensor_driver
