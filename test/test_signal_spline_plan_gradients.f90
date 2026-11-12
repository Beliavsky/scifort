! SPDX-License-Identifier: MIT
! Analytic value/coordinate-gradient checks for cached spline interpolation.
program test_signal_spline_plan_gradients
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_spline_nd, only : spline_interpolation_plan, &
        prepare_spline_plan_nd, evaluate_spline_plan_nd, &
        evaluate_spline_plan_with_gradient_nd
    implicit none
    character(len=16), parameter :: modes(8)=[character(len=16) :: &
        'mirror','constant','grid-constant','nearest', &
        'reflect','grid-mirror','wrap','grid-wrap']
    type(spline_interpolation_plan) :: plan,unset
    real(dp) :: x(7,5),xs(1,5),coords(2,10),p(2,1),v(10),g(2,10)
    real(dp) :: left(1),right(1),vp(10),eps,err,wrong(3,10)
    real(dp) :: rank15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: highpos(15,2),highv(2),highgrad(15,2)
    complex(dp) :: z(7,5),zc(10),zg(2,10),zl(1),zr(1)
    integer :: i,j,k,imode,degree,code,checks

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
    checks=0
    eps=1e-5_dp
    do imode=1,8
        do degree=0,5
            call prepare_spline_plan_nd(x,plan,order=degree,mode=trim(modes(imode)), &
                cval=0.25_dp,status=code)
            call check(code==0,'prepare real plan')
            call evaluate_spline_plan_with_gradient_nd(plan,coords,v,g,code)
            call check(code==0,'evaluate real gradients')
            call evaluate_spline_plan_nd(plan,coords,vp,code)
            call check(code==0,'evaluate original values')
            call check(maxval(abs(v-vp))<2e-12_dp,'real value identity')
            do k=1,size(v)
                do j=1,2
                    p(:,1)=coords(:,k)
                    p(j,1)=coords(j,k)+eps
                    call evaluate_spline_plan_nd(plan,p,right,code)
                    call check(code==0,'right interpolation')
                    p(j,1)=coords(j,k)-eps
                    call evaluate_spline_plan_nd(plan,p,left,code)
                    call check(code==0,'left interpolation')
                    err=abs(g(j,k)-(right(1)-left(1))/(2*eps))
                    call check(err<3e-4_dp,'real finite-difference derivative')
                    checks=checks+1
                end do
            end do
            call prepare_spline_plan_nd(z,plan,order=degree,mode=trim(modes(imode)), &
                cval=cmplx(0.25_dp,-0.13_dp,dp),status=code)
            call check(code==0,'prepare complex plan')
            call evaluate_spline_plan_with_gradient_nd(plan,coords,zc,zg,code)
            call check(code==0,'evaluate complex gradients')
            do k=1,size(v)
                do j=1,2
                    p(:,1)=coords(:,k)
                    p(j,1)=coords(j,k)+eps
                    call evaluate_spline_plan_nd(plan,p,zr,code)
                    call check(code==0,'right complex interpolation')
                    p(j,1)=coords(j,k)-eps
                    call evaluate_spline_plan_nd(plan,p,zl,code)
                    call check(code==0,'left complex interpolation')
                    err=abs(zg(j,k)-(zr(1)-zl(1))/(2*eps))
                    call check(err<3e-4_dp,'complex finite-difference derivative')
                    checks=checks+1
                end do
            end do
        end do
        do degree=0,5
            call prepare_spline_plan_nd(xs,plan,order=degree,mode=trim(modes(imode)), &
                cval=0.25_dp,status=code)
            call check(code==0,'prepare singleton plan')
            call evaluate_spline_plan_with_gradient_nd(plan,coords,v,g,code)
            call check(code==0,'evaluate singleton gradients')
            do k=1,size(v)
                do j=1,2
                    p(:,1)=coords(:,k)
                    p(j,1)=coords(j,k)+eps
                    call evaluate_spline_plan_nd(plan,p,right,code)
                    call check(code==0,'right singleton interpolation')
                    p(j,1)=coords(j,k)-eps
                    call evaluate_spline_plan_nd(plan,p,left,code)
                    call check(code==0,'left singleton interpolation')
                    err=abs(g(j,k)-(right(1)-left(1))/(2*eps))
                    call check(err<3e-4_dp,'singleton finite-difference derivative')
                    checks=checks+1
                end do
            end do
        end do
    end do
    call evaluate_spline_plan_with_gradient_nd(unset,coords,v,g,code)
    call check(code/=0,'reject uninitialized plan')
    call evaluate_spline_plan_with_gradient_nd(plan,coords,v,wrong,code)
    call check(code/=0,'reject wrong gradient shape')
    call evaluate_spline_plan_with_gradient_nd(plan,coords,zc,zg,code)
    call check(code/=0,'reject wrong output kind')
    coords(1,1)=ieee_value(1.0_dp,ieee_quiet_nan)
    call evaluate_spline_plan_with_gradient_nd(plan,coords,v,g,code)
    call check(code/=0,'reject nonfinite query coordinates')
    rank15=3.5_dp
    rank15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=6.5_dp
    call prepare_spline_plan_nd(rank15,plan,order=1,mode='nearest',status=code)
    call check(code==0,'prepare high-rank')
    highpos=1.21_dp
    highpos(1,:)=[1.2_dp,1.7_dp]
    call evaluate_spline_plan_with_gradient_nd(plan,highpos,highv,highgrad,code)
    call check(code==0,'rank-15 derivatives')
    call check(maxval(abs(highgrad(2:,:)))<1e-12_dp,'singleton axis derivatives')
    call check(maxval(abs(highgrad(1,:)-3.0_dp))<1e-12_dp,'nontrivial derivative')
    print '(a,i0)', 'PASS: spline plan analytic gradient comparisons ',checks
contains
    subroutine check(ok,message)
        logical, intent(in) :: ok !! result of comparison
        character(len=*), intent(in) :: message !! assertion description
        if (.not.ok) then
            print '(2a)', 'FAIL: ',message
            error stop 1
        end if
    end subroutine check
end program test_signal_spline_plan_gradients
