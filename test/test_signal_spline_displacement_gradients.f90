! SPDX-License-Identifier: MIT
! Gradients of cached displacement warps with respect to each offset.
program test_signal_spline_displacement_gradients
    use scifort_kinds, only : dp
    use scifort_signal_spline_nd, only : spline_interpolation_plan, &
        prepare_spline_plan_nd, evaluate_spline_plan_with_gradient_nd, &
        evaluate_spline_plan_displacement_with_gradient_nd
    implicit none
    character(len=16), parameter :: modes(8)=[character(len=16) :: &
        'mirror','constant','grid-constant','nearest', &
        'reflect','grid-mirror','wrap','grid-wrap']
    type(spline_interpolation_plan) :: plan,unused
    real(dp) :: x(11,9),dense(2,48),coords(2,48),y(8,6),v(48),g(2,48),gg(2,48)
    real(dp) :: large_g(2,96),noncontiguous(16,12)
    real(dp) :: highx(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: highy(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1),highdisp(15,2),highgrad(15,2)
    real(dp) :: a(2,1),offset(2),line(23),line_grad(2,23),ld(2,23),line_coords(2,23)
    complex(dp) :: z(11,9),zy(8,6),zg(2,48),zv(48),zgr(2,48)
    integer :: i,j,k,mode,degree,code,checks

    do j=1,9
        do i=1,11
            x(i,j)=sin(0.18_dp*i+0.31_dp*j)+0.01_dp*i*j
        end do
    end do
    z=cmplx(x,0.24_dp*x+0.07_dp,dp)
    k=0
    do j=1,6
        do i=1,8
            k=k+1
            dense(:,k)=[0.27_dp*sin(real(i+j,dp)*0.41_dp), &
                        -0.32_dp*cos(real(2*i-j,dp)*0.18_dp)]
            coords(:,k)=[real(i,dp),real(j,dp)]+dense(:,k)
        end do
    end do
    checks=0
    do mode=1,8
        do degree=0,5
            call prepare_spline_plan_nd(x,plan,order=degree, &
                mode=trim(modes(mode)),cval=0.37_dp,status=code)
            call check(code==0,'prepare real plan')
            call evaluate_spline_plan_with_gradient_nd(plan,coords,v,g,code)
            call check(code==0,'sample reference values and gradients')
            call evaluate_spline_plan_displacement_with_gradient_nd(plan,dense,y,gg,code,batch_size=5)
            call check(code==0,'real displacement value and gradients')
            call check(maxval(abs(reshape(y,[48])-v))<1e-12_dp,'real displacement samples')
            call check(maxval(abs(gg-g))<1e-12_dp,'real displacement gradients')
            noncontiguous=-999.0_dp
            large_g=-777.0_dp
            call evaluate_spline_plan_displacement_with_gradient_nd(plan,dense, &
                noncontiguous(1:16:2,1:12:2),large_g(:,1:96:2),code,batch_size=1)
            call check(code==0,'noncontiguous real image and gradients')
            call check(maxval(abs(noncontiguous(1:16:2,1:12:2)-y))<1e-12_dp, &
                'noncontiguous real image values')
            call check(maxval(abs(large_g(:,1:96:2)-g))<1e-12_dp, &
                'noncontiguous real gradients')
            call check(all(large_g(:,2:96:2)==-777.0_dp),'gradient stride untouched')
            call prepare_spline_plan_nd(z,plan,order=degree, &
                mode=trim(modes(mode)),cval=cmplx(0.37_dp,-0.21_dp,dp),status=code)
            call check(code==0,'prepare complex plan')
            call evaluate_spline_plan_with_gradient_nd(plan,coords,zv,zgr,code)
            call check(code==0,'sample complex reference')
            call evaluate_spline_plan_displacement_with_gradient_nd(plan,dense,zy,zg,code,batch_size=7)
            call check(code==0,'complex displacement and gradients')
            call check(maxval(abs(reshape(zy,[48])-zv))<1e-12_dp,'complex output')
            call check(maxval(abs(zg-zgr))<1e-12_dp,'complex displacement gradients')
            checks=checks+size(g)*2+size(y)*2
        end do
    end do
    call evaluate_spline_plan_displacement_with_gradient_nd(unused,dense,y,gg,code)
    call check(code/=0,'reject unprepared plan')
    call evaluate_spline_plan_displacement_with_gradient_nd(plan,dense,y,gg,code)
    call check(code/=0,'reject real output for complex source')
    call evaluate_spline_plan_displacement_with_gradient_nd(plan,dense(:,1:47),zy,zg,code)
    call check(code/=0,'reject wrong displacement length')
    call evaluate_spline_plan_displacement_with_gradient_nd(plan,dense,zy,zg(:,1:47),code)
    call check(code/=0,'reject wrong derivative length')
    call evaluate_spline_plan_displacement_with_gradient_nd(plan,dense,zy,zg,code,batch_size=0)
    call check(code/=0,'reject zero batch size')
    call prepare_spline_plan_nd(x,plan,order=3,mode='nearest',status=code)
    call check(code==0,'prepare rank-change source')
    a(:,1)=[0.37_dp,0.23_dp]
    offset=[0.2_dp,-0.1_dp]
    do i=1,size(line)
        ld(:,i)=[0.03_dp*sin(real(i,dp)),0.05_dp*cos(real(i,dp))]
        line_coords(:,i)=1.0_dp+offset+a(:,1)*real(i-1,dp)+ld(:,i)
    end do
    call evaluate_spline_plan_with_gradient_nd(plan,line_coords,v(:23),g(:,:23),code)
    call check(code==0,'rank-change reference')
    call evaluate_spline_plan_displacement_with_gradient_nd(plan,ld,line,line_grad,code, &
        matrix=a,offset=offset,batch_size=3)
    call check(code==0,'rank-changing affine and displacement')
    call check(maxval(abs(line-v(:23)))<1e-12_dp,'rank-change output')
    call check(maxval(abs(line_grad-g(:,:23)))<1e-12_dp,'rank-change gradient')
    highx=2.5_dp
    highx(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=5.5_dp
    call prepare_spline_plan_nd(highx,plan,order=1,mode='nearest',status=code)
    call check(code==0,'prepare rank-15 plan')
    highdisp=0.2_dp
    highdisp(1,2)=-0.2_dp
    call evaluate_spline_plan_displacement_with_gradient_nd(plan,highdisp, &
        highy,highgrad,code)
    call check(code==0,'rank-15 displacement gradient')
    call check(maxval(abs(highgrad(2:,:)))<1e-12_dp,'singleton derivatives')
    call check(maxval(abs(highgrad(1,:)-3.0_dp))<1e-12_dp,'rank-15 nonzero derivative')
    print '(a,i0)', 'PASS: spline displacement gradient comparisons ',checks
contains
    subroutine check(ok,message)
        logical, intent(in) :: ok !! test outcome
        character(len=*), intent(in) :: message !! diagnostic
        if (.not.ok) then
            print '(2a)', 'FAIL: ',message
            error stop 1
        end if
    end subroutine check
end program test_signal_spline_displacement_gradients
