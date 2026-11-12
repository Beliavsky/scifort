! SPDX-License-Identifier: MIT
! Grid-shaped spline gradients and affine-parameter Jacobians.
program test_signal_spline_grid_affine_gradients
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_spline_nd, only : spline_interpolation_plan, &
        prepare_spline_plan_nd, evaluate_spline_plan_nd, &
        evaluate_spline_plan_with_gradient_nd, &
        evaluate_spline_plan_grid_with_gradient_nd, &
        evaluate_spline_plan_affine_nd, evaluate_spline_plan_affine_with_gradient_nd
    implicit none
    character(len=16), parameter :: modes(8)=[character(len=16) :: &
        'mirror','constant','grid-constant','nearest', &
        'reflect','grid-mirror','wrap','grid-wrap']
    type(spline_interpolation_plan) :: plan,unset
    real(dp) :: x(8,7),coords(2,5,4),cflat(2,20),y(5,4),g(2,5,4), &
        flat(20),gf(2,20),a(2,2),off(2),dm(2,2,20),db(2,20)
    real(dp) :: big(10,8),bigg(4,10,8),bigc(4,10,8),vplus(5,4),vminus(5,4)
    real(dp) :: wrong(3,2,20),eps,hi_a(15,15),hi_off(15),hi_dm(15,15,2),hi_db(15,2)
    real(dp) :: rank15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: rank15y(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: grid14(15,2,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: grad14(15,2,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: missing14(2,2,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: rankchange(2,1),rankoffset(2),ranky(12),rankdm(2,1,12),rankdb(2,12)
    real(dp) :: rankcoords(2,12),rankvalues(12),rankgrads(2,12)
    real(dp) :: img14(2,1,1,1,1,1,1,1,1,1,1,1,1,1)
    complex(dp) :: z(8,7),zy(5,4),zg(2,5,4),zdm(2,2,20),zdb(2,20), &
        zflat(20),zflatg(2,20),zplus(5,4),zminus(5,4)
    integer :: i,j,k,m,d,code,checks

    do j=1,7
        do i=1,8
            x(i,j)=sin(0.16_dp*i+0.39_dp*j)+0.04_dp*i*j
        end do
    end do
    z=cmplx(x,0.31_dp*x-0.15_dp,dp)
    a(1,:)=[0.23_dp,0.12_dp]
    a(2,:)=[-0.06_dp,0.28_dp]
    off=[1.47_dp,1.28_dp]
    do j=1,4
        do i=1,5
            coords(:,i,j)=1.0_dp+off+matmul(a,real([i-1,j-1],dp))
        end do
    end do
    cflat=reshape(coords,[2,20])
    bigc=-999.0_dp
    bigc(1:4:2,1:10:2,1:8:2)=coords
    eps=2.0e-5_dp
    checks=0
    do m=1,8
        do d=0,5
            call prepare_spline_plan_nd(x,plan,order=d,mode=trim(modes(m)), &
                cval=0.35_dp,status=code)
            call check(code==0,'prepare real')
            call evaluate_spline_plan_grid_with_gradient_nd(plan,coords,y,g,code)
            call check(code==0,'real grid evaluation')
            call evaluate_spline_plan_with_gradient_nd(plan,cflat,flat,gf,code)
            call check(code==0,'flat reference')
            call check(maxval(abs(reshape(y,[20])-flat))<1e-12_dp,'grid values')
            call check(maxval(abs(reshape(g,[40])-reshape(gf,[40])))<1e-12_dp,'grid derivatives')
            big=-123.0_dp
            bigg=-456.0_dp
            call evaluate_spline_plan_grid_with_gradient_nd(plan, &
                bigc(1:4:2,1:10:2,1:8:2),big(1:10:2,1:8:2), &
                bigg(1:4:2,1:10:2,1:8:2),code)
            call check(code==0,'noncontiguous coordinate and output')
            call check(maxval(abs(big(1:10:2,1:8:2)-y))<1e-12_dp,'strided grid values')
            call check(maxval(abs(bigg(1:4:2,1:10:2,1:8:2)-g))<1e-12_dp,'strided grid gradient')
            call check(all(big(2:10:2,:)==-123.0_dp),'unselected output untouched')
            call evaluate_spline_plan_affine_with_gradient_nd(plan,a,off,y,dm,db,code,batch_size=3)
            call check(code==0,'affine value and parameter derivatives')
            call check(maxval(abs(reshape(y,[20])-flat))<1e-12_dp,'affine values')
            call check(maxval(abs(db-gf))<1e-12_dp,'offset Jacobian')
            k=0
            do j=1,4
                do i=1,5
                    k=k+1
                    call check(maxval(abs(dm(:,1,k)-gf(:,k)*real(i-1,dp)))<1e-12_dp, &
                        'matrix Jacobian first axis')
                    call check(maxval(abs(dm(:,2,k)-gf(:,k)*real(j-1,dp)))<1e-12_dp, &
                        'matrix Jacobian second axis')
                end do
            end do
            checks=checks+180
            call prepare_spline_plan_nd(z,plan,order=d,mode=trim(modes(m)), &
                cval=cmplx(0.35_dp,-0.24_dp,dp),status=code)
            call check(code==0,'prepare complex')
            call evaluate_spline_plan_grid_with_gradient_nd(plan,coords,zy,zg,code)
            call check(code==0,'complex grid')
            call evaluate_spline_plan_with_gradient_nd(plan,cflat,zflat,zflatg,code)
            call check(code==0,'complex flat')
            call check(maxval(abs(reshape(zy,[20])-zflat))<1e-12_dp,'complex grid values')
            call check(maxval(abs(reshape(zg,[40])-reshape(zflatg,[40])))<1e-12_dp, &
                'complex grid gradient')
            call evaluate_spline_plan_affine_with_gradient_nd(plan,a,off,zy,zdm,zdb,code,batch_size=7)
            call check(code==0,'complex affine')
            call check(maxval(abs(reshape(zy,[20])-zflat))<1e-12_dp,'complex affine values')
            call check(maxval(abs(zdb-zflatg))<1e-12_dp,'complex offset Jacobian')
            k=0
            do j=1,4
                do i=1,5
                    k=k+1
                    call check(maxval(abs(zdm(:,1,k)-zflatg(:,k)*real(i-1,dp)))<1e-12_dp, &
                        'complex matrix Jacobian axis 1')
                    call check(maxval(abs(zdm(:,2,k)-zflatg(:,k)*real(j-1,dp)))<1e-12_dp, &
                        'complex matrix Jacobian axis 2')
                end do
            end do
            checks=checks+180
        end do
    end do
    ! Finite differences independently verify several affine parameters.
    call prepare_spline_plan_nd(x,plan,order=3,mode='mirror',status=code)
    call check(code==0,'finite-difference plan')
    call evaluate_spline_plan_affine_with_gradient_nd(plan,a,off,y,dm,db,code)
    call check(code==0,'finite-difference reference')
    do i=1,2
        off(i)=off(i)+eps
        call evaluate_spline_plan_affine_nd(plan,a,off,vplus,code)
        call check(code==0,'offset finite diff plus')
        off(i)=off(i)-2*eps
        call evaluate_spline_plan_affine_nd(plan,a,off,vminus,code)
        call check(code==0,'offset finite diff minus')
        off(i)=off(i)+eps
        call check(maxval(abs((reshape(vplus,[20])-reshape(vminus,[20]))/(2*eps)-db(i,:))) &
            <5e-6_dp,'offset finite difference')
        do j=1,2
            a(i,j)=a(i,j)+eps
            call evaluate_spline_plan_affine_nd(plan,a,off,vplus,code)
            call check(code==0,'matrix finite diff plus')
            a(i,j)=a(i,j)-2*eps
            call evaluate_spline_plan_affine_nd(plan,a,off,vminus,code)
            call check(code==0,'matrix finite diff minus')
            a(i,j)=a(i,j)+eps
            call check(maxval(abs((reshape(vplus,[20])-reshape(vminus,[20]))/(2*eps)-dm(i,j,:))) &
                <5e-6_dp,'matrix finite difference')
        end do
    end do
    ! Error paths and rank-15 affine sensitivity.
    call evaluate_spline_plan_grid_with_gradient_nd(unset,coords,y,g,code)
    call check(code/=0,'uninitialized plan')
    call evaluate_spline_plan_grid_with_gradient_nd(plan,coords,y,wrong,code)
    call check(code/=0,'reject malformed grid gradient')
    call evaluate_spline_plan_affine_with_gradient_nd(plan,a,off,y,wrong,db,code)
    call check(code/=0,'reject malformed affine derivative')
    call evaluate_spline_plan_affine_with_gradient_nd(plan,a,off,y,dm,db,code,batch_size=0)
    call check(code/=0,'reject bad batch size')
    coords(1,1,1)=ieee_value(1.0_dp,ieee_quiet_nan)
    call evaluate_spline_plan_grid_with_gradient_nd(plan,coords,y,g,code)
    call check(code/=0,'nonfinite grid coordinate')
    rank15=3.0_dp
    rank15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=6.0_dp
    call prepare_spline_plan_nd(rank15,plan,order=1,mode='nearest',status=code)
    call check(code==0,'rank-15 plan')
    hi_a=0.0_dp
    do i=1,15
        hi_a(i,i)=1.0_dp
    end do
    hi_a(1,1)=0.3_dp
    hi_off=0.0_dp
    hi_off(1)=0.25_dp
    call evaluate_spline_plan_affine_with_gradient_nd(plan,hi_a,hi_off, &
        rank15y,hi_dm,hi_db,code,batch_size=1)
    call check(code==0,'rank-15 affine gradients')
    call check(maxval(abs(hi_db(1,:)-3.0_dp))<1e-12_dp,'rank-15 derivative first dimension')
    call check(maxval(abs(hi_db(2:,:)))<1e-12_dp,'rank-15 singleton derivatives')
    ! Fourteen destination dimensions permit a leading source coordinate dimension.
    grid14=1.0_dp
    grid14(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=1.25_dp
    grid14(1,2,1,1,1,1,1,1,1,1,1,1,1,1,1)=1.75_dp
    call evaluate_spline_plan_grid_with_gradient_nd(plan,grid14,img14,grad14,code)
    call check(code==0,'rank-14 output with rank-15 source')
    call check(abs(img14(1,1,1,1,1,1,1,1,1,1,1,1,1,1)-3.75_dp)<1e-12_dp, &
        'high rank grid first value')
    call check(abs(grad14(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)-3.0_dp)<1e-12_dp, &
        'high rank grid derivative')
    missing14=1.0_dp
    call evaluate_spline_plan_grid_with_gradient_nd(plan,missing14,img14,grad14,code)
    call check(code/=0,'reject grid with wrong input rank')
    ! A 1-D output from a 2-D input also has two independent affine offsets.
    call prepare_spline_plan_nd(x,plan,order=3,mode='mirror',status=code)
    call check(code==0,'prepare rank-change')
    rankchange(:,1)=[0.31_dp,0.22_dp]
    rankoffset=[1.17_dp,0.92_dp]
    do i=1,12
        rankcoords(:,i)=1.0_dp+rankoffset+rankchange(:,1)*real(i-1,dp)
    end do
    call evaluate_spline_plan_with_gradient_nd(plan,rankcoords,rankvalues,rankgrads,code)
    call check(code==0,'rank-change reference')
    call evaluate_spline_plan_affine_with_gradient_nd(plan,rankchange,rankoffset, &
        ranky,rankdm,rankdb,code,batch_size=2)
    call check(code==0,'rank-changing affine parameter derivatives')
    call check(maxval(abs(ranky-rankvalues))<1e-12_dp,'rank-changing values')
    call check(maxval(abs(rankdb-rankgrads))<1e-12_dp,'rank-changing offsets')
    do i=1,12
        call check(maxval(abs(rankdm(:,1,i)-real(i-1,dp)*rankgrads(:,i))) &
            <1e-12_dp,'rank-changing matrix gradient')
    end do
    print '(a,i0)', 'PASS: grid and affine parameter gradient comparisons ',checks
contains
    subroutine check(ok,message)
        logical, intent(in) :: ok !! assertion
        character(len=*), intent(in) :: message !! diagnostic
        if (.not.ok) then
            print '(2a)', 'FAIL: ',message
            error stop 1
        end if
    end subroutine check
end program test_signal_spline_grid_affine_gradients
