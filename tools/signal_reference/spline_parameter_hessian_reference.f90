! SPDX-License-Identifier: MIT
! Machine-readable analytic affine/displacement derivatives for SciPy checks.
program spline_parameter_hessian_reference
    use scifort_kinds, only : dp
    use scifort_signal_spline_nd, only : spline_interpolation_plan, &
        prepare_spline_plan_nd, evaluate_spline_plan_affine_with_hessian_nd, &
        evaluate_spline_plan_displacement_with_hessian_nd
    implicit none
    character(len=16), parameter :: modes(4)=[character(len=16) :: &
        'mirror','grid-wrap','nearest','grid-constant']
    integer, parameter :: degrees(3)=[2,3,5]
    real(dp) :: source(8,7),a(2,2),off(2),f(2,20),y(5,4), &
        dm(2,2,20),db(2,20),hm(2,2,2,2,20),hc(2,2,2,20),hb(2,2,20), &
        wy(5,4),wg(2,20),wh(2,2,20),coords(2),in_warp(2)
    type(spline_interpolation_plan) :: plan
    integer :: i,j,k,m,d,code

    do j=1,7
        do i=1,8
            source(i,j)=sin(0.19_dp*i+0.29_dp*j)+0.015_dp*i*j
        end do
    end do
    a(1,:)=[0.32_dp,0.12_dp]
    a(2,:)=[0.09_dp,0.26_dp]
    off=[1.19_dp,1.31_dp]
    k=0
    do j=1,4
        do i=1,5
            k=k+1
            f(:,k)=[0.03_dp*i,-0.02_dp*j]
        end do
    end do
    do m=1,4
        do d=1,3
            call prepare_spline_plan_nd(source,plan,order=degrees(d), &
                mode=trim(modes(m)),cval=0.37_dp,status=code)
            if (code/=0) error stop 'failed to prepare plan'
            call evaluate_spline_plan_affine_with_hessian_nd(plan,a,off,y,dm,db,hm,hc,hb,code,batch_size=7)
            if (code/=0) error stop 'affine reference failed'
            call evaluate_spline_plan_displacement_with_hessian_nd(plan,f,wy,wg,wh,code, &
                matrix=a,offset=off,batch_size=3)
            if (code/=0) error stop 'displacement reference failed'
            k=0
            do j=1,4
                do i=1,5
                    k=k+1
                    coords=1.0_dp+off+matmul(a,real([i-1,j-1],dp))
                    in_warp=coords+f(:,k)
                    write(*,'(3(i0,1x),*(es24.16,1x))') m,degrees(d),k, &
                        coords(1),coords(2),in_warp(1),in_warp(2), &
                        y(i,j),dm(1,1,k),db(1,k),hb(1,1,k),hb(1,2,k), &
                        hc(1,1,2,k),hm(1,1,2,2,k),wy(i,j),wg(1,k),wh(1,1,k),wh(1,2,k)
                end do
            end do
        end do
    end do
end program spline_parameter_hessian_reference
