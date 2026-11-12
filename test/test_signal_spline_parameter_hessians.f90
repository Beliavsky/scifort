! SPDX-License-Identifier: MIT
! Affine-parameter and displacement-field analytic spline Hessians.
program test_signal_spline_parameter_hessians
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_spline_nd, only : spline_interpolation_plan, &
        prepare_spline_plan_nd, evaluate_spline_plan_with_hessian_nd, &
        evaluate_spline_plan_affine_with_gradient_nd, &
        evaluate_spline_plan_affine_with_hessian_nd, &
        evaluate_spline_plan_displacement_with_gradient_nd, &
        evaluate_spline_plan_displacement_with_hessian_nd, &
        evaluate_spline_plan_displacement_grid_with_hessian_nd
    implicit none
    character(len=16), parameter :: modes(8)=[character(len=16) :: &
        'mirror','constant','grid-constant','nearest', &
        'reflect','grid-mirror','wrap','grid-wrap']
    real(dp), parameter :: step=1.0e-5_dp
    type(spline_interpolation_plan) :: plan,unset
    real(dp) :: image(8,7),matrix(2,2),offset(2),coords(2,20),y(5,4), &
        dm(2,2,20),db(2,20),hmm(2,2,2,2,20),hmo(2,2,2,20),hoo(2,2,20), &
        v(20),g(2,20),h(2,2,20),mp(2,2),op(2),yp(5,4), &
        dmp(2,2,20),dbp(2,20),dmn(2,2,20),dbn(2,20), &
        flow(2,5,4),flat_flow(2,20),yg(5,4),gd(2,20),hd(2,2,20), &
        ygrid(5,4),gg(2,5,4),hg(2,2,5,4), &
        gp(2,20),gm(2,20),bad_h(3,2,20),wrong_m(2,2,2,2,19), &
        big_y(10,8),big_field(4,10,8),big_g(4,10,8),big_h(4,4,10,8), &
        rank15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1), &
        high_matrix(15,15),high_offset(15),high_y(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), &
        high_dm(15,15,1),high_db(15,1),high_hm(15,15,15,15,1), &
        high_hc(15,15,15,1),high_hb(15,15,1), &
        low_a(2,1),low_y(6),low_dm(2,1,6),low_db(2,6), &
        low_hm(2,1,2,1,6),low_hc(2,1,2,6),low_hb(2,2,6), &
        low_coords(2,6),low_v(6),low_g(2,6),low_h(2,2,6), &
        high_flow(15,1),high_grad(15,1),high_hess(15,15,1)
    complex(dp) :: z(8,7),zy(5,4),zdm(2,2,20),zdb(2,20), &
        zhmm(2,2,2,2,20),zhmo(2,2,2,20),zhoo(2,2,20), &
        zyg(5,4),zgd(2,20),zhd(2,2,20),zgrid(5,4), &
        zgg(2,5,4),zhg(2,2,5,4)
    integer :: i,j,k,a,b,c,d,degree,mode,code,checks

    checks=0
    do j=1,7
        do i=1,8
            image(i,j)=sin(0.19_dp*i+0.29_dp*j)+0.015_dp*i*j
        end do
    end do
    z=cmplx(image,0.43_dp*image+0.1_dp,dp)
    matrix(1,:)=[0.32_dp,0.12_dp]
    matrix(2,:)=[0.09_dp,0.26_dp]
    offset=[1.19_dp,1.31_dp]
    k=0
    do j=1,4
        do i=1,5
            k=k+1
            coords(:,k)=1.0_dp+offset+matmul(matrix,real([i-1,j-1],dp))
        end do
    end do
    flow=0.0_dp
    do j=1,4
        do i=1,5
            flow(:,i,j)=[0.03_dp*i,-0.02_dp*j]
        end do
    end do
    flat_flow=reshape(flow,[2,20])
    do mode=1,8
        do degree=0,5
            call prepare_spline_plan_nd(image,plan,order=degree,mode=trim(modes(mode)), &
                cval=0.37_dp,status=code)
            call check(code==0,'real preparation')
            call evaluate_spline_plan_affine_with_hessian_nd(plan,matrix,offset,y,dm,db, &
                hmm,hmo,hoo,code,batch_size=3)
            call check(code==0,'affine Hessian real')
            call evaluate_spline_plan_with_hessian_nd(plan,coords,v,g,h,code)
            call check(code==0,'coordinate Hessian')
            call check(maxval(abs(reshape(y,[20])-v))<5e-12_dp,'real values')
            call check(maxval(abs(db-g))<5e-12_dp,'offset gradient')
            call check(maxval(abs(hoo-h))<5e-12_dp,'offset Hessian')
            k=0
            do j=1,4
                do i=1,5
                    k=k+1
                    do a=1,2
                        do b=1,2
                            call check(abs(hmo(a,1,b,k)-real(i-1,dp)*h(a,b,k))<5e-12_dp, &
                                'matrix-offset first output axis')
                            call check(abs(hmo(a,2,b,k)-real(j-1,dp)*h(a,b,k))<5e-12_dp, &
                                'matrix-offset second output axis')
                            do c=1,2
                                call check(abs(hmm(a,1,b,c,k)-hmo(a,1,b,k)*real( &
                                    merge(i-1,j-1,c==1),dp))<5e-12_dp,'matrix-matrix Hessian')
                            end do
                        end do
                    end do
                end do
            end do
            call evaluate_spline_plan_displacement_with_hessian_nd(plan,flat_flow,yg,gd,hd, &
                code,matrix=matrix,offset=offset,batch_size=7)
            call check(code==0,'displacement Hessian flat')
            call evaluate_spline_plan_displacement_grid_with_hessian_nd(plan,flow,ygrid,gg,hg, &
                code,matrix=matrix,offset=offset,batch_size=4)
            call check(code==0,'displacement Hessian grid')
            call check(maxval(abs(yg-ygrid))<5e-12_dp,'displacement grid values')
            call check(maxval(abs(reshape(gg,[40])-reshape(gd,[40])))<5e-12_dp,'grid gradients')
            call check(maxval(abs(reshape(hg,[80])-reshape(hd,[80])))<5e-12_dp,'grid Hessians')
            call check(maxval(abs(hd(1,2,:)-hd(2,1,:)))<5e-12_dp,'displacement symmetry')
            checks=checks+20*17
            call prepare_spline_plan_nd(z,plan,order=degree,mode=trim(modes(mode)), &
                cval=cmplx(0.37_dp,0.43_dp*0.37_dp+0.1_dp,dp),status=code)
            call check(code==0,'complex preparation')
            call evaluate_spline_plan_affine_with_hessian_nd(plan,matrix,offset,zy,zdm,zdb, &
                zhmm,zhmo,zhoo,code,batch_size=5)
            call check(code==0,'complex affine Hessian')
            call check(maxval(abs(real(zhmm,dp)-hmm))<5e-11_dp,'complex affine real part')
            call check(maxval(abs(aimag(zhmm)-0.43_dp*hmm))<5e-11_dp,'complex affine imaginary')
            call check(maxval(abs(aimag(zhoo)-0.43_dp*hoo))<5e-11_dp,'complex offset Hessian')
            call evaluate_spline_plan_displacement_with_hessian_nd(plan,flat_flow,zyg,zgd,zhd, &
                code,matrix=matrix,offset=offset)
            call check(code==0,'complex displacement Hessian')
            call evaluate_spline_plan_displacement_grid_with_hessian_nd(plan,flow,zgrid,zgg,zhg, &
                code,matrix=matrix,offset=offset)
            call check(code==0,'complex grid displacement Hessian')
            call check(maxval(abs(reshape(zhg,[80])-reshape(zhd,[80])))<5e-12_dp, &
                'complex grid Hessian')
            checks=checks+20*12
        end do
    end do
    ! Independent derivatives: finite differences of the earlier gradient API.
    call prepare_spline_plan_nd(image,plan,order=3,mode='mirror',status=code)
    call check(code==0,'finite difference prepare')
    call evaluate_spline_plan_affine_with_hessian_nd(plan,matrix,offset,y,dm,db, &
        hmm,hmo,hoo,code)
    call check(code==0,'finite difference affine')
    do a=1,2
        op=offset; op(a)=op(a)+step
        call evaluate_spline_plan_affine_with_gradient_nd(plan,matrix,op,yp,dmp,dbp,code)
        call check(code==0,'positive offset gradient')
        op(a)=offset(a)-step
        call evaluate_spline_plan_affine_with_gradient_nd(plan,matrix,op,yp,dmn,dbn,code)
        call check(code==0,'negative offset gradient')
        call check(maxval(abs((dbp-dbn)/(2*step)-hoo(:,a,:)))<2e-4_dp, &
            'offset-offset finite difference')
        call check(maxval(abs((dmp-dmn)/(2*step)-hmo(:,:,a,:)))<2e-4_dp, &
            'matrix-offset finite difference')
        do b=1,2
            mp=matrix; mp(a,b)=matrix(a,b)+step
            call evaluate_spline_plan_affine_with_gradient_nd(plan,mp,offset,yp,dmp,dbp,code)
            call check(code==0,'positive matrix gradient')
            mp(a,b)=matrix(a,b)-step
            call evaluate_spline_plan_affine_with_gradient_nd(plan,mp,offset,yp,dmn,dbn,code)
            call check(code==0,'negative matrix gradient')
            call check(maxval(abs((dmp-dmn)/(2*step)-hmm(:,:,a,b,:)))<2e-4_dp, &
                'matrix-matrix finite difference')
            call check(maxval(abs((dbp-dbn)/(2*step)-hmo(a,b,:,:)))<2e-4_dp, &
                'offset-matrix finite difference')
        end do
    end do
    call evaluate_spline_plan_displacement_with_hessian_nd(plan,flat_flow,yg,gd,hd, &
        code,matrix=matrix,offset=offset)
    call check(code==0,'finite difference displacement')
    do a=1,2
        flat_flow(a,:)=flat_flow(a,:)+step
        call evaluate_spline_plan_displacement_with_gradient_nd(plan,flat_flow,yp,gp,code, &
            matrix=matrix,offset=offset)
        call check(code==0,'positive displacement gradient')
        flat_flow(a,:)=flat_flow(a,:)-2*step
        call evaluate_spline_plan_displacement_with_gradient_nd(plan,flat_flow,yp,gm,code, &
            matrix=matrix,offset=offset)
        call check(code==0,'negative displacement gradient')
        flat_flow(a,:)=flat_flow(a,:)+step
        call check(maxval(abs((gp-gm)/(2*step)-hd(:,a,:)))<2e-4_dp, &
            'displacement Hessian finite difference')
    end do
    ! Different input/output ranks: two input axes, one output axis.
    call prepare_spline_plan_nd(image,plan,order=3,mode='mirror',status=code)
    call check(code==0,'rank-changing preparation')
    low_a(:,1)=[0.24_dp,0.14_dp]
    call evaluate_spline_plan_affine_with_hessian_nd(plan,low_a,offset,low_y, &
        low_dm,low_db,low_hm,low_hc,low_hb,code,batch_size=2)
    call check(code==0,'rank-changing affine Hessian')
    do i=1,6
        low_coords(:,i)=1.0_dp+offset+low_a(:,1)*real(i-1,dp)
    end do
    call evaluate_spline_plan_with_hessian_nd(plan,low_coords,low_v,low_g,low_h,code)
    call check(code==0,'rank-changing coordinate reference')
    call check(maxval(abs(low_y-low_v))<1e-12_dp,'rank-changing values')
    call check(maxval(abs(low_hb-low_h))<1e-12_dp,'rank-changing offset Hessian')
    do i=1,6
        call check(maxval(abs(low_hm(:,1,:,1,i)-low_h(:,:,i)*real((i-1)**2,dp)))<1e-12_dp, &
            'rank-changing matrix Hessian')
    end do
    ! Noncontiguous grid outputs.
    big_field=-9.0_dp
    big_field(1:4:2,1:10:2,1:8:2)=flow
    big_y=-7.0_dp
    big_g=-7.0_dp
    big_h=-7.0_dp
    call evaluate_spline_plan_displacement_grid_with_hessian_nd(plan, &
        big_field(1:4:2,1:10:2,1:8:2),big_y(1:10:2,1:8:2), &
        big_g(1:4:2,1:10:2,1:8:2),big_h(1:4:2,1:4:2,1:10:2,1:8:2),code, &
        matrix=matrix,offset=offset)
    call check(code==0,'noncontiguous displacement grid')
    call check(maxval(abs(big_y(1:10:2,1:8:2)-yg))<5e-12_dp, &
        'strided displacement values')
    call check(maxval(abs(big_h(1:4:2,1:4:2,1:10:2,1:8:2)- &
        reshape(hd,[2,2,5,4])))<5e-12_dp,'strided displacement Hessian')
    call check(all(big_y(2:10:2,:)==-7.0_dp),'untouched output padding')
    ! Validation and error paths.
    call evaluate_spline_plan_affine_with_hessian_nd(unset,matrix,offset,y,dm,db, &
        hmm,hmo,hoo,code)
    call check(code/=0,'unset affine plan')
    call evaluate_spline_plan_affine_with_hessian_nd(plan,matrix,offset,y,dm,db, &
        wrong_m,hmo,hoo,code)
    call check(code/=0,'reject wrong matrix Hessian shape')
    call evaluate_spline_plan_affine_with_hessian_nd(plan,matrix,offset,y,dm,db, &
        hmm,hmo,hoo,code,batch_size=0)
    call check(code/=0,'reject bad batch size')
    call evaluate_spline_plan_displacement_with_hessian_nd(plan,flat_flow,yg,gd,bad_h,code)
    call check(code/=0,'reject wrong displacement Hessian dimensions')
    call evaluate_spline_plan_displacement_grid_with_hessian_nd(plan,flow,ygrid,gg,gg,code)
    call check(code/=0,'reject wrong grid Hessian shape')
    flat_flow(1,1)=ieee_value(1.0_dp,ieee_quiet_nan)
    call evaluate_spline_plan_displacement_with_hessian_nd(plan,flat_flow,yg,gd,hd,code)
    call check(code/=0,'reject NaN displacement')
    ! Rank-15 data with a rank-15 output must use a flattened Hessian layout.
    rank15=2.0_dp
    rank15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=3.0_dp
    call prepare_spline_plan_nd(rank15,plan,order=1,mode='nearest',status=code)
    call check(code==0,'rank-15 preparation')
    high_matrix=0.0_dp
    do i=1,15
        high_matrix(i,i)=1.0_dp
    end do
    high_matrix(1,1)=0.3_dp
    high_offset=0.0_dp
    high_offset(1)=0.23_dp
    call evaluate_spline_plan_affine_with_hessian_nd(plan,high_matrix,high_offset,high_y, &
        high_dm,high_db,high_hm,high_hc,high_hb,code,batch_size=1)
    call check(code==0,'rank-15 affine Hessian')
    call check(maxval(abs(high_hm))<1e-12_dp,'rank-15 first-order diagonal Hessian')
    high_flow=0.0_dp
    high_flow(1,1)=0.12_dp
    call evaluate_spline_plan_displacement_with_hessian_nd(plan,high_flow,high_y, &
        high_grad,high_hess,code,batch_size=1)
    call check(code==0,'rank-15 displacement Hessian')
    call check(maxval(abs(high_hess))<1e-12_dp,'rank-15 first-order displacement Hessian')
    print '(a,i0)', 'PASS: analytic affine and displacement Hessian values ',checks
contains
    subroutine check(ok,what)
        logical, intent(in) :: ok !! assertion result
        character(len=*), intent(in) :: what !! failing check description
        if (.not.ok) then
            print '(2a)', 'FAIL: ',what
            error stop 1
        end if
    end subroutine check
end program test_signal_spline_parameter_hessians
