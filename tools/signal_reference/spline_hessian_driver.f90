! SPDX-License-Identifier: MIT
! Numerical records for independent SciPy Hessian comparisons.
program spline_hessian_driver
    use scifort_kinds, only : dp
    use scifort_signal_spline_nd, only : spline_interpolation_plan, &
        prepare_spline_plan_nd, evaluate_spline_plan_with_hessian_nd
    implicit none
    character(len=16), parameter :: modes(8)=[character(len=16) :: &
        'mirror','constant','grid-constant','nearest', &
        'reflect','grid-mirror','wrap','grid-wrap']
    real(dp) :: x(17,16),c(2,8),v(8),g(2,8),h(2,2,8)
    complex(dp) :: z(17,16),zv(8),zg(2,8),zh(2,2,8)
    complex(dp) :: result(6)
    type(spline_interpolation_plan) :: plan
    integer :: i,j,k,m,d,kind,code
    do j=1,16
        do i=1,17
            x(i,j)=sin(0.19_dp*i+0.23_dp*j)+0.003_dp*i*j
        end do
    end do
    z=cmplx(x,0.37_dp*x-0.12_dp,dp)
    c(:,1)=[2.35_dp,3.42_dp]
    c(:,2)=[7.28_dp,5.67_dp]
    c(:,3)=[14.17_dp,12.39_dp]
    c(:,4)=[-1.23_dp,4.51_dp]
    c(:,5)=[18.22_dp,10.73_dp]
    c(:,6)=[5.33_dp,-0.74_dp]
    c(:,7)=[8.41_dp,16.32_dp]
    c(:,8)=[1.23_dp,1.47_dp]
    do m=1,8
        do d=0,5
            do kind=0,1
                if (kind==0) then
                    call prepare_spline_plan_nd(x,plan,order=d,mode=trim(modes(m)), &
                        cval=0.43_dp,status=code)
                    if (code/=0) error stop 'real plan'
                    call evaluate_spline_plan_with_hessian_nd(plan,c,v,g,h,code)
                    if (code/=0) error stop 'real Hessian'
                else
                    call prepare_spline_plan_nd(z,plan,order=d,mode=trim(modes(m)), &
                        cval=cmplx(0.43_dp,-0.29_dp,dp),status=code)
                    if (code/=0) error stop 'complex plan'
                    call evaluate_spline_plan_with_hessian_nd(plan,c,zv,zg,zh,code)
                    if (code/=0) error stop 'complex Hessian'
                end if
                do k=1,8
                    if (kind==0) then
                        result=[cmplx(v(k),0.0_dp,dp),cmplx(g(1,k),0.0_dp,dp), &
                            cmplx(g(2,k),0.0_dp,dp),cmplx(h(1,1,k),0.0_dp,dp), &
                            cmplx(h(1,2,k),0.0_dp,dp),cmplx(h(2,2,k),0.0_dp,dp)]
                    else
                        result=[zv(k),zg(1,k),zg(2,k),zh(1,1,k),zh(1,2,k),zh(2,2,k)]
                    end if
                    write(*,'(4(i0,1x),12(es24.15e3,1x))') m,d,kind,k, &
                        (real(result(i),dp),aimag(result(i)),i=1,6)
                end do
            end do
        end do
    end do
end program spline_hessian_driver
