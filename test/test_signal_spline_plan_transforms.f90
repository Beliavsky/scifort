! SPDX-License-Identifier: MIT
! Regression checks for streaming affine and callback spline-plan evaluation.
program test_signal_spline_plan_transforms
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_spline_nd, only : spline_interpolation_plan, &
        prepare_spline_plan_nd, evaluate_spline_plan_affine_nd, &
        evaluate_spline_plan_geometric_nd, evaluate_spline_plan_nd, &
        affine_transform_nd, spline_nd_status_too_expensive
    implicit none
    character(len=16), parameter :: modes(8)=[character(len=16) :: &
        'mirror','constant','grid-constant','nearest','reflect', &
        'grid-mirror','wrap','grid-wrap']
    type(spline_interpolation_plan) :: plan
    real(dp) :: source(8,7), ref(10,9), a(10,9), mapped(10,9), tall(20,18)
    real(dp) :: matrix(2,2),offset(2),rect(2,1),line(17),points(2,17),checkline(17)
    real(dp) :: invalid_mat(3,2),nonfinite_offset(2), nan_value
    complex(dp) :: z(8,7),zr(10,9),za(10,9),zc(10,9)
    real(dp) :: x15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: y15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: eye15(15,15), zeros15(15)
    integer :: i,j,k,modeid,degree,code,called,checked

    do j=1,7
        do i=1,8
            source(i,j)=sin(0.19_dp*i+0.31_dp*j)+0.02_dp*i*j
        end do
    end do
    z=cmplx(source,0.17_dp-0.43_dp*source,dp)
    matrix=reshape([1.08_dp,0.23_dp,-0.16_dp,0.93_dp],[2,2])
    offset=[-0.72_dp,0.42_dp]
    checked=0
    do modeid=1,8
        do degree=0,5
            call prepare_spline_plan_nd(source,plan,order=degree, &
                mode=trim(modes(modeid)),cval=1.23_dp,status=code)
            call check(code==0,'prepare real plan')
            call affine_transform_nd(source,matrix,offset,ref,order=degree, &
                mode=trim(modes(modeid)),prefilter=.true.,cval=1.23_dp,status=code)
            call check(code==0,'direct affine transformation')
            call evaluate_spline_plan_affine_nd(plan,matrix,offset,a,code,batch_size=1)
            call check(code==0,'single-coordinate batches')
            call check(maxval(abs(a-ref))<2e-10_dp,'affine matches direct')
            called=0
            call evaluate_spline_plan_geometric_nd(plan,affine_mapping,mapped,code,batch_size=7)
            call check(code==0,'geometric batches')
            call check(called==size(mapped),'one callback per output pixel')
            call check(maxval(abs(mapped-ref))<2e-10_dp,'mapping matches direct')
            tall=-888.0_dp
            call evaluate_spline_plan_affine_nd(plan,matrix,offset, &
                tall(1:20:2,1:18:2),code,batch_size=17)
            call check(code==0,'strided output')
            call check(maxval(abs(tall(1:20:2,1:18:2)-ref))<2e-10_dp, &
                'strided output numerical values')
            call check(all(tall(2:20:2,:)==-888.0_dp),'strided destination untouched')
            checked=checked+size(ref)*3
        end do
    end do
    call prepare_spline_plan_nd(z,plan,order=5,mode='nearest', &
        cval=cmplx(0.3_dp,1.2_dp,dp),status=code)
    call check(code==0,'complex plan')
    call affine_transform_nd(z,matrix,offset,zr,order=5,mode='nearest', &
        prefilter=.true.,cval=cmplx(0.3_dp,1.2_dp,dp),status=code)
    call check(code==0,'complex direct affine')
    call evaluate_spline_plan_affine_nd(plan,matrix,offset,za,code,batch_size=19)
    call check(code==0,'complex affine plan')
    call evaluate_spline_plan_geometric_nd(plan,affine_mapping,zc,code,batch_size=5)
    call check(code==0,'complex geometric plan')
    call check(maxval(abs(za-zr))<2e-10_dp,'complex affine agreement')
    call check(maxval(abs(zc-zr))<2e-10_dp,'complex callback agreement')
    checked=checked+4*size(zr)
    call evaluate_spline_plan_affine_nd(plan,matrix,offset,a,code)
    call check(code/=0,'real output rejects complex plan')
    call prepare_spline_plan_nd(source,plan,order=3,mode='reflect',status=code)
    call check(code==0,'real plan for rank change')
    rect(:,1)=[0.37_dp,0.12_dp]
    do i=1,size(line)
        points(:,i)=1.0_dp+offset+rect(:,1)*real(i-1,dp)
    end do
    call evaluate_spline_plan_nd(plan,points,checkline,code)
    call check(code==0,'flat nonseparable mapping reference')
    call evaluate_spline_plan_affine_nd(plan,rect,offset,line,code,batch_size=3)
    call check(code==0,'affine reduces 2-D source to 1-D output')
    call check(maxval(abs(line-checkline))<2e-10_dp,'rectangular affine matrix')
    call evaluate_spline_plan_geometric_nd(plan,rectangular_mapping,line,code)
    call check(code==0,'callback changes output rank')
    call check(maxval(abs(line-checkline))<2e-10_dp,'rank-changing callback')

    invalid_mat=0.0_dp
    call evaluate_spline_plan_affine_nd(plan,invalid_mat,offset,a,code)
    call check(code/=0,'reject wrong matrix shape')
    call evaluate_spline_plan_affine_nd(plan,matrix,offset,a,code,batch_size=0)
    call check(code/=0,'reject zero batch size')
    call evaluate_spline_plan_geometric_nd(plan,affine_mapping,a,code,batch_size=9000)
    call check(code/=0,'reject excessive batch size')
    nan_value=ieee_value(0.0_dp,ieee_quiet_nan)
    nonfinite_offset=[nan_value,0.0_dp]
    call evaluate_spline_plan_affine_nd(plan,matrix,nonfinite_offset,a,code)
    call check(code/=0,'reject nonfinite offset')
    call evaluate_spline_plan_geometric_nd(plan,nan_mapping,a,code)
    call check(code/=0,'reject callback NaN')

    x15=0.0_dp
    x15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=4.0_dp
    x15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=7.0_dp
    eye15=0.0_dp
    do i=1,15
        eye15(i,i)=1.0_dp
    end do
    zeros15=0.0_dp
    call prepare_spline_plan_nd(x15,plan,order=1,mode='nearest',status=code)
    call check(code==0,'rank-15 preparation')
    call evaluate_spline_plan_affine_nd(plan,eye15,zeros15,y15,code)
    call check(code==0,'rank-15 evaluation')
    call check(maxval(abs(y15-x15))<2e-10_dp,'rank-15 values')
    print '(a,i0)', 'PASS: streaming spline-plan transforms, scalar values checked ',checked
contains
    subroutine check(ok,message)
        logical, intent(in) :: ok !! assertion
        character(len=*), intent(in) :: message !! failure description
        if (.not.ok) then
            print '(2a)', 'FAIL: ',message
            error stop 1
        end if
    end subroutine check

    subroutine affine_mapping(output_coords,input_coords)
        integer, intent(in) :: output_coords(:) !! one-based destination coordinate
        real(dp), intent(out) :: input_coords(:) !! one-based source coordinate
        called=called+1
        input_coords=1.0_dp+offset+matmul(matrix,real(output_coords-1,dp))
    end subroutine affine_mapping

    subroutine rectangular_mapping(output_coords,input_coords)
        integer, intent(in) :: output_coords(:) !! one-based destination coordinate
        real(dp), intent(out) :: input_coords(:) !! one-based source coordinate
        input_coords=1.0_dp+offset+matmul(rect,real(output_coords-1,dp))
    end subroutine rectangular_mapping

    subroutine nan_mapping(output_coords,input_coords)
        integer, intent(in) :: output_coords(:) !! one-based destination coordinate
        real(dp), intent(out) :: input_coords(:) !! deliberately invalid coordinate
        input_coords=nan_value
    end subroutine nan_mapping
end program test_signal_spline_plan_transforms
