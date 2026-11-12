! SPDX-License-Identifier: MIT
! Public ndimage re-exports of the spline Hessian and warp-gradient APIs.
program test_ndimage_spline_hessian_exports
    use scifort_kinds, only : dp
    use scifort_ndimage, only : spline_interpolation_plan, prepare_spline_plan_nd, &
        evaluate_spline_plan_with_hessian_nd, evaluate_spline_plan_grid_with_hessian_nd, &
        evaluate_spline_plan_displacement_grid_with_gradient_nd, &
        evaluate_spline_plan_affine_with_hessian_nd, &
        evaluate_spline_plan_displacement_grid_with_hessian_nd
    implicit none
    type(spline_interpolation_plan) :: plan
    real(dp) :: source(5,4), coords(2,2,2),v(4),g(2,4),h(2,2,4)
    real(dp) :: flatcoords(2,4),y(2,2),gg(2,2,2),hg(2,2,2,2)
    real(dp) :: displacements(2,2,2),warp(2,2),sens(2,2,2)
    real(dp) :: m(2,2),off(2),dm(2,2,4),db(2,4),hm(2,2,2,2,4), &
        hc(2,2,2,4),hb(2,2,4),hs(2,2,2,2)
    integer :: code,i,j

    do j=1,4
        do i=1,5
            source(i,j)=sin(real(i+2*j,dp)*0.2_dp)
        end do
    end do
    coords(:,1,1)=[2.1_dp,1.2_dp]
    coords(:,2,1)=[3.2_dp,1.7_dp]
    coords(:,1,2)=[2.3_dp,2.4_dp]
    coords(:,2,2)=[3.4_dp,2.1_dp]
    flatcoords=reshape(coords,[2,4])
    call prepare_spline_plan_nd(source,plan,order=3,mode='mirror',status=code)
    if (code/=0) error stop 'public plan construction'
    call evaluate_spline_plan_with_hessian_nd(plan,flatcoords,v,g,h,code)
    if (code/=0) error stop 'public Hessian'
    call evaluate_spline_plan_grid_with_hessian_nd(plan,coords,y,gg,hg,code)
    if (code/=0) error stop 'public shaped Hessian'
    if (maxval(abs(y-reshape(v,[2,2])))>1e-12_dp) error stop 'public values'
    if (maxval(abs(hg-reshape(h,[2,2,2,2])))>1e-12_dp) error stop 'public Hessian mismatch'
    displacements=0.0_dp
    call evaluate_spline_plan_displacement_grid_with_gradient_nd(plan,displacements,warp,sens,code)
    if (code/=0) error stop 'public displacement Jacobian'
    if (maxval(abs(sens))>1.0e3_dp) error stop 'unexpected warp gradient'
    m=0.0_dp
    m(1,1)=0.4_dp
    m(2,2)=0.4_dp
    off=[1.1_dp,0.3_dp]
    call evaluate_spline_plan_affine_with_hessian_nd(plan,m,off,warp,dm,db, &
        hm,hc,hb,code)
    if (code/=0) error stop 'public affine Hessians'
    if (maxval(abs(hm))>1.0e3_dp) error stop 'unexpected affine Hessian'
    call evaluate_spline_plan_displacement_grid_with_hessian_nd(plan,displacements, &
        warp,sens,hs,code)
    if (code/=0) error stop 'public grid displacement Hessians'
    if (maxval(abs(hs))>1.0e3_dp) error stop 'unexpected displacement Hessian'
    print '(a)', 'PASS: public ndimage spline Hessians and displacement gradients'
end program test_ndimage_spline_hessian_exports
