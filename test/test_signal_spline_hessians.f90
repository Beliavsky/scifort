! SPDX-License-Identifier: MIT
! Analytic spline Hessians and grid displacement sensitivities.
program test_signal_spline_hessians
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_spline_nd, only : spline_interpolation_plan, &
        prepare_spline_plan_nd, evaluate_spline_plan_with_gradient_nd, &
        evaluate_spline_plan_with_hessian_nd, evaluate_spline_plan_grid_with_hessian_nd, &
        evaluate_spline_plan_displacement_grid_with_gradient_nd, &
        evaluate_spline_plan_displacement_with_gradient_nd
    implicit none
    character(len=16), parameter :: modes(8)=[character(len=16) :: &
        'mirror','constant','grid-constant','nearest', &
        'reflect','grid-mirror','wrap','grid-wrap']
    real(dp), parameter :: step=1.0e-5_dp
    type(spline_interpolation_plan) :: plan,unset
    real(dp) :: x(8,7),coords(2,6),xp(2,6),v(6),g(2,6),h(2,2,6), &
        vp(6),gm(2,6),gp(2,6),value_old(6),grad_old(2,6)
    real(dp) :: gc(2,3,2),vg(3,2),gg(2,3,2),hg(2,2,3,2)
    real(dp) :: bigv(6,4),bigg(4,6,4),bigh(4,4,6,4),bigc(4,6,4)
    real(dp) :: flow(2,3,2),w(3,2),dw(2,3,2),wflat(3,2),df(2,6), &
        mat(2,2),off(2),bad(3,3,2),singleton(1,6),csmall(2,1), &
        bigflow(4,6,4),bigw(6,4),bigdw(4,6,4), &
        svalue(1),sg(2,1),sh(2,2,1),sp(2,1),sm(2,1), &
        hp(2,1),hm(2,1), v15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1), &
        cc15(15,1),vv15(1),gg15(15,1),hh15(15,15,1)
    complex(dp) :: z(8,7),zv(6),zg(2,6),zh(2,2,6), &
        zold(6),zgold(2,6),zvg(3,2),zgg(2,3,2),zhg(2,2,3,2), &
        zw(3,2),zdw(2,3,2),zflat(3,2),zdf(2,6)
    integer :: i,j,k,d,m,code,checks

    checks=0
    do j=1,7
        do i=1,8
            x(i,j)=sin(0.23_dp*i+0.17_dp*j)+0.015_dp*i*j
        end do
    end do
    z=cmplx(x,0.37_dp*x,dp)
    coords(:,1)=[2.23_dp,3.34_dp]
    coords(:,2)=[5.19_dp,4.76_dp]
    coords(:,3)=[-0.23_dp,2.28_dp]
    coords(:,4)=[8.35_dp,7.26_dp]
    coords(:,5)=[4.27_dp,-0.45_dp]
    coords(:,6)=[3.38_dp,6.72_dp]
    gc=reshape(coords,shape(gc))
    mat(1,:)=[0.23_dp,0.04_dp]
    mat(2,:)=[0.07_dp,0.27_dp]
    off=[2.13_dp,1.67_dp]
    flow=0.0_dp
    flow(1,:,:)=0.09_dp
    flow(2,:,:)=-0.04_dp
    do m=1,8
        do d=0,5
            call prepare_spline_plan_nd(x,plan,order=d,mode=trim(modes(m)), &
                cval=0.48_dp,status=code)
            call check(code==0,'real plan')
            call evaluate_spline_plan_with_hessian_nd(plan,coords,v,g,h,code)
            call check(code==0,'real Hessian')
            call evaluate_spline_plan_with_gradient_nd(plan,coords,value_old,grad_old,code)
            call check(code==0,'existing gradient')
            call check(maxval(abs(v-value_old))<2e-12_dp,'values unchanged')
            call check(maxval(abs(g-grad_old))<2e-12_dp,'gradients unchanged')
            call check(maxval(abs(h(1,2,:)-h(2,1,:)))<1e-13_dp,'symmetric Hessian')
            if (d==0) call check(maxval(abs(h))<1e-12_dp,'zero-order Hessian')
            if (d==1) then
                call check(maxval(abs(h(1,1,:)))<1e-12_dp,'linear spline diagonal 1')
                call check(maxval(abs(h(2,2,:)))<1e-12_dp,'linear spline diagonal 2')
            end if
            if (d>=1) then
                do i=1,2
                    xp=coords
                    xp(i,:)=coords(i,:)+step
                    call evaluate_spline_plan_with_gradient_nd(plan,xp,vp,gp,code)
                    call check(code==0,'gradient plus')
                    xp(i,:)=coords(i,:)-step
                    call evaluate_spline_plan_with_gradient_nd(plan,xp,vp,gm,code)
                    call check(code==0,'gradient minus')
                    call check(maxval(abs(h(i,:,:)-(gp-gm)/(2.0_dp*step)))<1e-4_dp, &
                        'analytic Hessian agrees with central derivative')
                end do
            end if
            call evaluate_spline_plan_grid_with_hessian_nd(plan,gc,vg,gg,hg,code)
            call check(code==0,'real grid Hessian')
            call check(maxval(abs(reshape(vg,[6])-v))<2e-12_dp,'grid values')
            call check(maxval(abs(reshape(gg,[12])-reshape(g,[12])))<2e-12_dp,'grid gradient')
            call check(maxval(abs(reshape(hg,[24])-reshape(h,[24])))<2e-12_dp,'grid Hessian')
            bigc=-99.0_dp
            bigc(1:4:2,1:6:2,1:4:2)=gc
            bigv=-77.0_dp
            bigg=-77.0_dp
            bigh=-77.0_dp
            call evaluate_spline_plan_grid_with_hessian_nd(plan, &
                bigc(1:4:2,1:6:2,1:4:2),bigv(1:6:2,1:4:2), &
                bigg(1:4:2,1:6:2,1:4:2),bigh(1:4:2,1:4:2,1:6:2,1:4:2),code)
            call check(code==0,'noncontiguous Hessian output')
            call check(maxval(abs(bigh(1:4:2,1:4:2,1:6:2,1:4:2)-hg))<2e-12_dp, &
                'noncontiguous Hessian values')
            checks=checks+6*7
            call prepare_spline_plan_nd(z,plan,order=d,mode=trim(modes(m)), &
                cval=cmplx(0.48_dp,0.37_dp*0.48_dp,dp),status=code)
            call check(code==0,'complex plan')
            call evaluate_spline_plan_with_hessian_nd(plan,coords,zv,zg,zh,code)
            call check(code==0,'complex Hessian')
            call evaluate_spline_plan_with_gradient_nd(plan,coords,zold,zgold,code)
            call check(code==0,'complex existing gradient')
            call check(maxval(abs(zv-zold))<2e-12_dp,'complex value')
            call check(maxval(abs(zg-zgold))<2e-12_dp,'complex gradient')
            call check(maxval(abs(aimag(zh)-0.37_dp*real(zh,dp)))<2e-10_dp, &
                'complex Hessian linearity')
            call evaluate_spline_plan_grid_with_hessian_nd(plan,gc,zvg,zgg,zhg,code)
            call check(code==0,'complex grid Hessian')
            call check(maxval(abs(reshape(zhg,[24])-reshape(zh,[24])))<2e-12_dp, &
                'complex grid matches flat Hessian')
            checks=checks+6*7
        end do
    end do
    ! Unfiltered coefficients use the same analytic derivative basis.
    call prepare_spline_plan_nd(x,plan,order=5,mode='mirror', &
        prefilter=.false.,status=code)
    call check(code==0,'unfiltered plan')
    call evaluate_spline_plan_with_hessian_nd(plan,coords,v,g,h,code)
    call check(code==0,'unfiltered Hessian')
    xp=coords
    xp(1,:)=coords(1,:)+step
    call evaluate_spline_plan_with_gradient_nd(plan,xp,vp,gp,code)
    call check(code==0,'unfiltered plus')
    xp(1,:)=coords(1,:)-step
    call evaluate_spline_plan_with_gradient_nd(plan,xp,vp,gm,code)
    call check(code==0,'unfiltered minus')
    call check(maxval(abs(h(1,:,:)-(gp-gm)/(2*step)))<1e-4_dp, &
        'unfiltered Hessian finite difference')
    ! Dense displacements with grid-shaped sensitivities; compare the flat API.
    call prepare_spline_plan_nd(x,plan,order=3,mode='nearest',status=code)
    call check(code==0,'displacement real plan')
    call evaluate_spline_plan_displacement_grid_with_gradient_nd(plan,flow,w,dw,code, &
        matrix=mat,offset=off,batch_size=2)
    call check(code==0,'grid displacement gradients')
    call evaluate_spline_plan_displacement_with_gradient_nd(plan,reshape(flow,[2,6]), &
        wflat,df,code,matrix=mat,offset=off,batch_size=3)
    call check(code==0,'flat displacement gradients')
    call check(maxval(abs(w-wflat))<1e-12_dp,'displacement image')
    call check(maxval(abs(reshape(dw,[12])-reshape(df,[12])))<1e-12_dp, &
        'displacement derivatives')
    call evaluate_spline_plan_displacement_grid_with_gradient_nd(plan,bad,w,dw,code)
    call check(code/=0,'reject wrong displacement-component count')
    call evaluate_spline_plan_displacement_grid_with_gradient_nd(plan,flow,w,bad,code)
    call check(code/=0,'reject wrong gradient component count')
    bigflow=-99.0_dp
    bigflow(1:4:2,1:6:2,1:4:2)=flow
    bigw=-77.0_dp
    bigdw=-77.0_dp
    call evaluate_spline_plan_displacement_grid_with_gradient_nd(plan, &
        bigflow(1:4:2,1:6:2,1:4:2),bigw(1:6:2,1:4:2), &
        bigdw(1:4:2,1:6:2,1:4:2),code,matrix=mat,offset=off)
    call check(code==0,'strided displacement input and output')
    call check(maxval(abs(bigw(1:6:2,1:4:2)-w))<1e-12_dp,'strided warp values')
    call check(maxval(abs(bigdw(1:4:2,1:6:2,1:4:2)-dw))<1e-12_dp, &
        'strided warp gradient')
    call check(all(bigw(2:6:2,:)==-77.0_dp),'strided destination preserves padding')
    call prepare_spline_plan_nd(z,plan,order=3,mode='nearest',status=code)
    call check(code==0,'displacement complex plan')
    call evaluate_spline_plan_displacement_grid_with_gradient_nd(plan,flow,zw,zdw,code, &
        matrix=mat,offset=off)
    call check(code==0,'complex grid displacement')
    call evaluate_spline_plan_displacement_with_gradient_nd(plan,reshape(flow,[2,6]), &
        zflat,zdf,code,matrix=mat,offset=off)
    call check(code==0,'complex flat displacement')
    call check(maxval(abs(zw-zflat))<1e-12_dp,'complex grid image')
    call check(maxval(abs(reshape(zdw,[12])-reshape(zdf,[12])))<1e-12_dp, &
        'complex grid derivatives')
    call evaluate_spline_plan_with_hessian_nd(unset,coords,v,g,h,code)
    call check(code/=0,'reject unset plan')
    call evaluate_spline_plan_grid_with_hessian_nd(plan,gc,vg,gg,gg,code)
    call check(code/=0,'reject wrong Hessian shape')
    cc15=1.0_dp
    v15=4.0_dp
    v15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=5.0_dp
    call prepare_spline_plan_nd(v15,plan,order=1,mode='nearest',status=code)
    call check(code==0,'rank-15 plan')
    cc15(1,1)=1.23_dp
    call evaluate_spline_plan_with_hessian_nd(plan,cc15,vv15,gg15,hh15,code)
    call check(code==0,'rank-15 Hessian')
    call check(maxval(abs(hh15))<1e-12_dp,'rank-15 first order Hessian zero')
    ! Grid-constant singleton factors contribute second derivatives even when
    ! the other dimensions have been collapsed to reduce stencil costs.
    singleton(1,:)=[0.15_dp,0.7_dp,1.4_dp,2.1_dp,1.1_dp,0.2_dp]
    call prepare_spline_plan_nd(singleton,plan,order=3,mode='grid-constant', &
        cval=0.33_dp,status=code)
    call check(code==0,'singleton grid-constant plan')
    csmall(:,1)=[1.28_dp,3.12_dp]
    call evaluate_spline_plan_with_hessian_nd(plan,csmall,svalue,sg,sh,code)
    call check(code==0,'singleton Hessian')
    do i=1,2
        sp=csmall; sp(i,1)=sp(i,1)+step
        call evaluate_spline_plan_with_gradient_nd(plan,sp,vp(:1),hp,code)
        call check(code==0,'singleton gradient plus')
        sm=csmall; sm(i,1)=sm(i,1)-step
        call evaluate_spline_plan_with_gradient_nd(plan,sm,vp(:1),hm,code)
        call check(code==0,'singleton gradient minus')
        call check(maxval(abs(sh(i,:,1)-(hp(:,1)-hm(:,1))/(2*step)))<3e-4_dp, &
            'singleton Hessian finite differences')
    end do
    print '(a,i0)', 'PASS: analytic spline Hessian and displacement-grid checks ',checks
contains
    subroutine check(ok,what)
        logical, intent(in) :: ok !! assertion result
        character(len=*), intent(in) :: what !! failure context
        if (.not.ok) then
            print '(2a)', 'FAIL: ',what
            error stop 1
        end if
    end subroutine check
end program test_signal_spline_hessians
