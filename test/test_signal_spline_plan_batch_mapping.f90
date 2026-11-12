! SPDX-License-Identifier: MIT
! Checked native regressions for vectorized spline-coordinate callbacks.
program test_signal_spline_plan_batch_mapping
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_spline_nd, only : spline_interpolation_plan, &
        prepare_spline_plan_nd, evaluate_spline_plan_geometric_nd, evaluate_spline_plan_nd
    use scifort_ndimage, only : evaluate_spline_plan_geometric_batch_nd
    implicit none
    character(len=16), parameter :: modes(8)=[character(len=16) :: &
        'mirror','constant','grid-constant','nearest','reflect', &
        'grid-mirror','wrap','grid-wrap']
    type(spline_interpolation_plan) :: plan
    real(dp) :: x(16,13),a(19,17),b(19,17),strided(38,34)
    real(dp) :: line(27),reference(27),points(2,27),nan_value
    real(dp) :: highrank(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: highout(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    complex(dp) :: z(16,13),za(19,17),zb(19,17)
    integer :: i,j,d,degree,code,calls,points_seen,mode_id,checked
    integer :: batch_calls,largest_block

    do j=1,size(x,2)
        do i=1,size(x,1)
            x(i,j)=sin(0.19_dp*i+0.23_dp*j)+0.005_dp*i*j
        end do
    end do
    z=cmplx(x,0.43_dp*x-0.27_dp,dp)
    checked=0
    do mode_id=1,8
        do degree=0,5
            call prepare_spline_plan_nd(x,plan,order=degree,mode=trim(modes(mode_id)), &
                cval=-0.73_dp,status=code)
            call check(code==0,'real prepare')
            calls=0
            call evaluate_spline_plan_geometric_nd(plan,scalar_map,a,code,batch_size=23)
            call check(code==0,'scalar evaluation')
            call check(calls==size(a),'scalar callback count')
            batch_calls=0
            largest_block=0
            points_seen=0
            call evaluate_spline_plan_geometric_batch_nd(plan,batched_map,b,code,batch_size=7)
            call check(code==0,'batched real evaluation')
            call check(batch_calls==(size(b)+6)/7,'callback called once per block')
            call check(points_seen==size(b),'each point visited')
            call check(largest_block==7,'block size respected')
            call check(maxval(abs(a-b))<2e-13_dp,'batched vs scalar real values')
            strided=-999.0_dp
            call evaluate_spline_plan_geometric_batch_nd(plan,batched_map, &
                strided(1:38:2,1:34:2),code,batch_size=1)
            call check(code==0,'noncontiguous real destination')
            call check(maxval(abs(a-strided(1:38:2,1:34:2)))<2e-13_dp, &
                'noncontiguous values')
            call check(all(strided(2:38:2,:)==-999.0_dp),'untouched output cells')
            checked=checked+3*size(a)
        end do
    end do
    call prepare_spline_plan_nd(z,plan,order=5,mode='nearest', &
        cval=cmplx(-0.73_dp,0.18_dp,dp),status=code)
    call check(code==0,'complex prepare')
    call evaluate_spline_plan_geometric_nd(plan,scalar_map,za,code)
    call check(code==0,'complex scalar mapping')
    call evaluate_spline_plan_geometric_batch_nd(plan,batched_map,zb,code,batch_size=5)
    call check(code==0,'complex batch mapping')
    call check(maxval(abs(za-zb))<2e-13_dp,'complex callback values')
    checked=checked+2*size(zb)

    call evaluate_spline_plan_geometric_batch_nd(plan,batched_map,a,code)
    call check(code/=0,'reject output with wrong numeric type')
    call prepare_spline_plan_nd(x,plan,order=3,mode='reflect',status=code)
    call check(code==0,'rank-change preparation')
    do i=1,size(line)
        points(:,i)=[1.0_dp+0.43_dp*i,1.0_dp+0.17_dp*i]
    end do
    call evaluate_spline_plan_nd(plan,points,reference,code)
    call check(code==0,'rank-change flat-coordinate reference')
    call evaluate_spline_plan_geometric_batch_nd(plan,line_map,line,code,batch_size=4)
    call check(code==0,'rank-changing batched transform')
    call check(maxval(abs(line-reference))<1e-13_dp,'rank change numerical result')

    call evaluate_spline_plan_geometric_batch_nd(plan,batched_map,a,code,batch_size=0)
    call check(code/=0,'reject zero batch size')
    call evaluate_spline_plan_geometric_batch_nd(plan,batched_map,a,code,batch_size=8193)
    call check(code/=0,'reject excessive batch size')
    nan_value=ieee_value(1.0_dp,ieee_quiet_nan)
    call evaluate_spline_plan_geometric_batch_nd(plan,nan_map,a,code)
    call check(code/=0,'reject nonfinite callback coordinate')

    highrank=3.0_dp
    highrank(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=5.0_dp
    call prepare_spline_plan_nd(highrank,plan,order=1,mode='nearest',status=code)
    call check(code==0,'rank-15 source')
    call evaluate_spline_plan_geometric_batch_nd(plan,identity_15,highout,code)
    call check(code==0,'rank-15 batched transform')
    call check(maxval(abs(highrank-highout))<2e-13_dp,'rank-15 numerical values')
    print '(a,i0)', 'PASS: batched spline mapping native comparisons ',checked
contains
    subroutine check(ok,message)
        logical, intent(in) :: ok !! condition to verify
        character(len=*), intent(in) :: message !! failure description
        if (.not.ok) then
            print '(2a)', 'FAIL: ',message
            error stop 1
        end if
    end subroutine check

    subroutine scalar_map(out_coords,in_coords)
        integer, intent(in) :: out_coords(:) !! destination one-based indices
        real(dp), intent(out) :: in_coords(:) !! source one-based coordinates
        calls=calls+1
        in_coords(1)=0.52_dp+0.86_dp*out_coords(1)-0.13_dp*out_coords(2) &
            +0.11_dp*sin(0.13_dp*out_coords(2))
        in_coords(2)=0.72_dp+0.22_dp*out_coords(1)+0.79_dp*out_coords(2) &
            -0.17_dp*cos(0.21_dp*out_coords(1))
    end subroutine scalar_map

    subroutine batched_map(out_coords,in_coords)
        integer, intent(in) :: out_coords(:,:) !! destination indices by column
        real(dp), intent(out) :: in_coords(:,:) !! one-based mapped source positions
        batch_calls=batch_calls+1
        points_seen=points_seen+size(out_coords,2)
        largest_block=max(largest_block,size(out_coords,2))
        in_coords(1,:)=0.52_dp+0.86_dp*out_coords(1,:)-0.13_dp*out_coords(2,:) &
            +0.11_dp*sin(0.13_dp*out_coords(2,:))
        in_coords(2,:)=0.72_dp+0.22_dp*out_coords(1,:)+0.79_dp*out_coords(2,:) &
            -0.17_dp*cos(0.21_dp*out_coords(1,:))
    end subroutine batched_map

    subroutine line_map(out_coords,in_coords)
        integer, intent(in) :: out_coords(:,:) !! one-dimensional output index
        real(dp), intent(out) :: in_coords(:,:) !! two-dimensional image location
        in_coords(1,:)=1.0_dp+0.43_dp*out_coords(1,:)
        in_coords(2,:)=1.0_dp+0.17_dp*out_coords(1,:)
    end subroutine line_map

    subroutine nan_map(out_coords,in_coords)
        integer, intent(in) :: out_coords(:,:) !! destination positions
        real(dp), intent(out) :: in_coords(:,:) !! intentionally invalid coordinates
        in_coords=nan_value
    end subroutine nan_map

    subroutine identity_15(out_coords,in_coords)
        integer, intent(in) :: out_coords(:,:) !! one-based source coordinates
        real(dp), intent(out) :: in_coords(:,:) !! one-based copy
        in_coords=real(out_coords,dp)
    end subroutine identity_15
end program test_signal_spline_plan_batch_mapping
