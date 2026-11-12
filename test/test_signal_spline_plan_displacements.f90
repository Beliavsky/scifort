! SPDX-License-Identifier: MIT
! Test cached spline warps with dense/flat displacement fields.
program test_signal_spline_plan_displacements
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_ndimage, only : spline_interpolation_plan, prepare_spline_plan_nd, &
        evaluate_spline_plan_nd, evaluate_spline_plan_displacement_nd, &
        evaluate_spline_plan_displacement_grid_nd
    implicit none
    character(len=16), parameter :: modes(8)=[character(len=16) :: &
        'mirror','constant','grid-constant','nearest','reflect', &
        'grid-mirror','wrap','grid-wrap']
    type(spline_interpolation_plan) :: plan,unused
    real(dp) :: x(16,13), y(14,12),dense(2,14,12),flat(2,14*12)
    real(dp) :: output_flat(14,12),reference(14*12),coords(2,14*12)
    real(dp) :: noncontiguous(28,24),ax(2,1),off(2)
    real(dp) :: line(27),shift_line(2,27),target_line(27),coord_line(2,27)
    real(dp) :: highrank(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: highout(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1),highfield(15,2)
    real(dp) :: highout14(2,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: highgrid(15,2,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: highmatrix(15,14),highoffset(15)
    complex(dp) :: z(16,13),zy(14,12),zreference(14*12)
    integer :: i,j,k,code,degree,imode,checks
    real(dp) :: nan_value

    do j=1,13
        do i=1,16
            x(i,j)=sin(0.17_dp*i+0.29_dp*j)+0.002_dp*i*j
        end do
    end do
    z=cmplx(x,0.32_dp*x+0.1_dp,dp)
    k=0
    do j=1,12
        do i=1,14
            k=k+1
            dense(1,i,j)=0.21_dp*sin(real(i+j,dp)*0.4_dp)
            dense(2,i,j)=-0.31_dp*cos(real(2*i-j,dp)*0.27_dp)
            coords(:,k)=[real(i,dp),real(j,dp)]+dense(:,i,j)
        end do
    end do
    flat=reshape(dense,shape(flat))
    checks=0
    do imode=1,8
        do degree=0,5
            call prepare_spline_plan_nd(x,plan,order=degree,mode=trim(modes(imode)), &
                cval=-0.4_dp,status=code)
            call check(code==0,'prepare real plan')
            call evaluate_spline_plan_nd(plan,coords,reference,code)
            call check(code==0,'reference sample')
            call evaluate_spline_plan_displacement_nd(plan,flat,y,code,batch_size=7)
            call check(code==0,'flat displacement')
            call check(maxval(abs(reshape(y,[size(y)])-reference))<3e-13_dp, &
                'flat matches existing interpolation')
            call evaluate_spline_plan_displacement_grid_nd(plan,dense,output_flat,code,batch_size=3)
            call check(code==0,'dense displacement')
            call check(maxval(abs(y-output_flat))<3e-13_dp,'grid and flat results')
            noncontiguous=-999.0_dp
            call evaluate_spline_plan_displacement_nd(plan,flat, &
                noncontiguous(1:28:2,1:24:2),code,batch_size=1)
            call check(code==0,'strided output')
            call check(maxval(abs(noncontiguous(1:28:2,1:24:2)-y))<3e-13_dp, &
                'strided numerical results')
            call check(all(noncontiguous(2:28:2,:)==-999.0_dp),'strided untouched rows')
            checks=checks+3*size(y)
        end do
    end do

    call prepare_spline_plan_nd(z,plan,order=5,mode='nearest', &
        cval=cmplx(-0.4_dp,0.25_dp,dp),status=code)
    call check(code==0,'prepare complex plan')
    call evaluate_spline_plan_nd(plan,coords,zreference,code)
    call check(code==0,'complex reference')
    call evaluate_spline_plan_displacement_grid_nd(plan,dense,zy,code,batch_size=5)
    call check(code==0,'complex grid displacements')
    call check(maxval(abs(reshape(zy,[size(zy)])-zreference))<3e-13_dp,'complex grid values')
    checks=checks+2*size(zy)
    call evaluate_spline_plan_displacement_nd(plan,flat,y,code)
    call check(code/=0,'reject wrong output type')

    ! Affine base can map a rank-1 output into a rank-2 input.
    call prepare_spline_plan_nd(x,plan,order=3,mode='grid-wrap',status=code)
    call check(code==0,'prepare rank-change plan')
    ax(:,1)=[0.45_dp,0.27_dp]
    off=[0.1_dp,-0.3_dp]
    do i=1,27
        shift_line(:,i)=[0.1_dp*sin(0.31_dp*i),0.15_dp*cos(0.23_dp*i)]
        coord_line(:,i)=1.0_dp+off+ax(:,1)*real(i-1,dp)+shift_line(:,i)
    end do
    call evaluate_spline_plan_nd(plan,coord_line,target_line,code)
    call check(code==0,'rank-change interpolation reference')
    call evaluate_spline_plan_displacement_nd(plan,shift_line,line,code, &
        matrix=ax,offset=off,batch_size=4)
    call check(code==0,'rank-changing displacement')
    call check(maxval(abs(line-target_line))<3e-13_dp,'rank-changing numerical values')
    checks=checks+size(line)
    call evaluate_spline_plan_displacement_nd(plan,shift_line,line,code)
    call check(code/=0,'reject identity mapping with mismatched ranks')
    call evaluate_spline_plan_displacement_nd(plan,shift_line,line,code,matrix=ax)
    call check(code/=0,'reject matrix without offset')
    call evaluate_spline_plan_displacement_nd(plan,shift_line,line,code,offset=off)
    call check(code/=0,'reject offset without matrix')
    call evaluate_spline_plan_displacement_nd(plan,shift_line(:,1:26),line,code, &
        matrix=ax,offset=off)
    call check(code/=0,'reject incorrect number of displacement pixels')
    call evaluate_spline_plan_displacement_grid_nd(plan,dense(:,1:13,:),y,code)
    call check(code/=0,'reject wrong dense displacement shape')
    call evaluate_spline_plan_displacement_nd(unused,flat,y,code)
    call check(code/=0,'reject unprepared plan')
    call evaluate_spline_plan_displacement_nd(plan,flat,y,code,batch_size=0)
    call check(code/=0,'reject zero batch size')
    nan_value=ieee_value(1.0_dp,ieee_quiet_nan)
    flat(2,9)=nan_value
    call evaluate_spline_plan_displacement_nd(plan,flat,y,code)
    call check(code/=0,'reject nonfinite displacement')

    highrank=2.0_dp
    highrank(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=5.0_dp
    call prepare_spline_plan_nd(highrank,plan,order=1,mode='nearest',status=code)
    call check(code==0,'prepare rank-15 plan')
    highfield=0.0_dp
    call evaluate_spline_plan_displacement_nd(plan,highfield,highout,code)
    call check(code==0,'rank-15 output and displacement')
    call check(maxval(abs(highout-highrank))<3e-13_dp,'rank-15 interpolation')
    highgrid=0.0_dp
    highmatrix=0.0_dp
    highoffset=0.0_dp
    do i=1,14
        highmatrix(i,i)=1.0_dp
    end do
    call evaluate_spline_plan_displacement_grid_nd(plan,highgrid,highout14,code, &
        matrix=highmatrix,offset=highoffset)
    call check(code==0,'rank-15 field and rank-14 output')
    call check(maxval(abs(reshape(highout14,[2])-reshape(highrank,[2])))<3e-13_dp, &
        'rank-changing high-rank grid values')
    print '(a,i0)', 'PASS: spline displacement comparisons ',checks
contains
    subroutine check(ok,message)
        logical, intent(in) :: ok !! test condition
        character(len=*), intent(in) :: message !! failure detail
        if (.not.ok) then
            print '(2a)', 'FAIL: ',message
            error stop 1
        end if
    end subroutine check
end program test_signal_spline_plan_displacements
