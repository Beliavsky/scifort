! SPDX-License-Identifier: MIT
! Emit cached values and analytic derivatives for SciPy reference checks.
program spline_plan_gradient_driver
    use scifort_kinds, only : dp
    use scifort_signal_spline_nd, only : spline_interpolation_plan, &
        prepare_spline_plan_nd, evaluate_spline_plan_with_gradient_nd
    implicit none
    character(len=16), parameter :: modes(8)=[character(len=16) :: &
        'mirror','constant','grid-constant','nearest', &
        'reflect','grid-mirror','wrap','grid-wrap']
    type(spline_interpolation_plan) :: plan
    real(dp) :: x(7,5),xs(1,5),coords(2,10),v(10),g(2,10)
    complex(dp) :: z(7,5),cv(10),cg(2,10),fill
    integer :: i,j,k,kind,degree,imode,filtered,code

    do j=1,5
        do i=1,7
            x(i,j)=sin(0.21_dp*i+0.37_dp*j)+0.01_dp*i*j
        end do
    end do
    xs=x(1:1,:)
    z=cmplx(x,0.3_dp*x+0.04_dp,dp)
    coords(:,1)=[1.21_dp,1.37_dp]
    coords(:,2)=[3.42_dp,2.71_dp]
    coords(:,3)=[6.68_dp,4.21_dp]
    coords(:,4)=[-0.31_dp,1.83_dp]
    coords(:,5)=[8.34_dp,2.34_dp]
    coords(:,6)=[2.14_dp,-1.23_dp]
    coords(:,7)=[3.29_dp,7.41_dp]
    coords(:,8)=[-3.27_dp,8.31_dp]
    coords(:,9)=[5.32_dp,4.89_dp]
    coords(:,10)=[4.0_dp,3.0_dp]
    fill=cmplx(0.25_dp,-0.13_dp,dp)
    do filtered=0,1
        do imode=1,8
            do degree=0,5
                do kind=0,2
                    select case(kind)
                    case(0)
                        call prepare_spline_plan_nd(x,plan,order=degree,mode=trim(modes(imode)), &
                            prefilter=filtered==1,cval=real(fill,dp),status=code)
                    case(1)
                        call prepare_spline_plan_nd(z,plan,order=degree,mode=trim(modes(imode)), &
                            prefilter=filtered==1,cval=fill,status=code)
                    case(2)
                        call prepare_spline_plan_nd(xs,plan,order=degree,mode=trim(modes(imode)), &
                            prefilter=filtered==1,cval=real(fill,dp),status=code)
                    end select
                    if (code/=0) error stop 'prepare'
                    if (kind==1) then
                        call evaluate_spline_plan_with_gradient_nd(plan,coords,cv,cg,code)
                    else
                        call evaluate_spline_plan_with_gradient_nd(plan,coords,v,g,code)
                        cv=cmplx(v,0.0_dp,dp)
                        cg=cmplx(g,0.0_dp,dp)
                    end if
                    if (code/=0) error stop 'gradient'
                    do k=1,size(v)
                        write(*,'(5(i0,1x),6(es24.16,1x))') imode,degree,kind,filtered,k, &
                            real(cv(k),dp),aimag(cv(k)),real(cg(1,k),dp),aimag(cg(1,k)), &
                            real(cg(2,k),dp),aimag(cg(2,k))
                    end do
                end do
            end do
        end do
    end do
end program spline_plan_gradient_driver
