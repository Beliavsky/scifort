! SPDX-License-Identifier: MIT
! Regression tests for separable tensor-product spline interpolation.
program test_signal_spline_tensor
    use scifort_kinds, only : dp
    use scifort_signal_spline_nd, only : spline_interpolation_plan, spline_axis_grid, &
        prepare_spline_plan_nd, evaluate_spline_plan_nd, evaluate_spline_plan_tensor_nd
    implicit none
    character(len=16), parameter :: modes(8) = [character(len=16) :: &
        'mirror', 'constant', 'grid-constant', 'nearest', 'reflect', &
        'grid-mirror', 'wrap', 'grid-wrap']
    type(spline_interpolation_plan) :: plan
    type(spline_axis_grid) :: grid(2), grid15(15)
    real(dp) :: x(5,4), y(7,6), baseline(7,6), large(14,12)
    real(dp) :: coords(2,42), flat(42), xsingle(1,4), ys(7,6)
    real(dp) :: x15(2,2,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: y15(2,2,1,1,1,1,1,1,1,1,1,1,1,1,1)
    complex(dp) :: cx(5,4), cy(7,6), cbase(7,6)
    integer :: i, j, k, degree, code, modeid, filtered, cases
    real(dp) :: tol

    do j=1,4
        do i=1,5
            x(i,j)=sin(real(2*i+3*j,dp))+0.07_dp*i*j
        end do
    end do
    cx=cmplx(x,-0.37_dp*x+0.2_dp,dp)
    grid(1)%coordinates=[-2.1_dp,0.15_dp,1.0_dp,2.7_dp,4.99_dp,5.4_dp,8.0_dp]
    grid(2)%coordinates=[-1.8_dp,0.3_dp,1.0_dp,2.66_dp,4.15_dp,6.1_dp]
    do j=1,6
        do i=1,7
            k=i+7*(j-1)
            coords(:,k)=[grid(1)%coordinates(i),grid(2)%coordinates(j)]
        end do
    end do
    cases=0
    do modeid=1,size(modes)
        do degree=0,5
            do filtered=0,1
                call prepare_spline_plan_nd(x,plan,order=degree, &
                    mode=trim(modes(modeid)),prefilter=filtered==1, &
                    cval=2.625_dp,status=code)
                if (code/=0) error stop 'real prepare failed'
                call evaluate_spline_plan_nd(plan,coords,flat,code)
                if (code/=0) error stop 'real pointwise failed'
                baseline=reshape(flat,shape(baseline))
                call evaluate_spline_plan_tensor_nd(plan,grid,y,code)
                if (code/=0) error stop 'real tensor failed'
                tol=3.0e-10_dp
                if (maxval(abs(y-baseline))>tol) then
                    print *, trim(modes(modeid)),degree,filtered,maxval(abs(y-baseline))
                    error stop 'real tensor mismatch'
                end if
                call evaluate_spline_plan_tensor_nd(plan,grid,large(1:14:2,1:12:2),code)
                if (code/=0) error stop 'noncontiguous tensor failed'
                if (maxval(abs(large(1:14:2,1:12:2)-baseline))>tol) &
                    error stop 'noncontiguous mismatch'
                cases=cases+42
                call prepare_spline_plan_nd(cx,plan,order=degree, &
                    mode=trim(modes(modeid)),prefilter=filtered==1, &
                    cval=cmplx(2.625_dp,-0.4_dp,dp),status=code)
                if (code/=0) error stop 'complex prepare failed'
                block
                    complex(dp) :: cflat(42)
                    call evaluate_spline_plan_nd(plan,coords,cflat,code)
                    if (code/=0) error stop 'complex pointwise failed'
                    cbase=reshape(cflat,shape(cbase))
                end block
                call evaluate_spline_plan_tensor_nd(plan,grid,cy,code)
                if (code/=0) error stop 'complex tensor failed'
                if (maxval(abs(cy-cbase))>tol) then
                    print *, trim(modes(modeid)),degree,filtered,maxval(abs(cy-cbase))
                    error stop 'complex tensor mismatch'
                end if
                cases=cases+42
            end do
        end do
    end do
    ! A grid-constant singleton axis is evaluated using its filtered impulse.
    xsingle(1,:)=x(1,:)
    call prepare_spline_plan_nd(xsingle,plan,order=5,mode='grid-constant', &
        prefilter=.true.,cval=1.5_dp,status=code)
    if (code/=0) error stop 'singleton prepare failed'
    call evaluate_spline_plan_nd(plan,coords,flat,code)
    if (code/=0) error stop 'singleton pointwise failed'
    call evaluate_spline_plan_tensor_nd(plan,grid,ys,code)
    if (code/=0) error stop 'singleton tensor failed'
    if (maxval(abs(ys-reshape(flat,shape(ys))))>1.e-10_dp) &
        error stop 'singleton tensor mismatch'
    ! The rank-15 evaluator works without constructing rank-specific wrappers.
    x15=0.0_dp
    x15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=2.0_dp
    x15(2,2,1,1,1,1,1,1,1,1,1,1,1,1,1)=3.0_dp
    do i=1,15
        grid15(i)%coordinates=[1.0_dp]
    end do
    grid15(1)%coordinates=[1.0_dp,2.0_dp]
    grid15(2)%coordinates=[1.0_dp,2.0_dp]
    call prepare_spline_plan_nd(x15,plan,order=0,status=code)
    if (code/=0) error stop 'rank-15 prepare failed'
    call evaluate_spline_plan_tensor_nd(plan,grid15,y15,code)
    if (code/=0 .or. any(y15/=x15)) error stop 'rank-15 mismatch'
    deallocate(grid15(1)%coordinates)
    call evaluate_spline_plan_tensor_nd(plan,grid15,y15,code)
    if (code==0) error stop 'missing coordinates should fail'
    print '(a,i0)', 'tensor reference values passed: ',cases+42+size(x15)
end program test_signal_spline_tensor
