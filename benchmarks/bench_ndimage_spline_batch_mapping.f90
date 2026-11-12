! SPDX-License-Identifier: MIT
! Batched versus scalar callback overhead with the same cached spline plan.
program bench_ndimage_spline_batch_mapping
    use scifort_kinds, only : dp
    use scifort_signal_spline_nd, only : spline_interpolation_plan, &
        prepare_spline_plan_nd, evaluate_spline_plan_geometric_nd, &
        evaluate_spline_plan_geometric_batch_nd
    implicit none
    type(spline_interpolation_plan) :: plan
    real(dp) :: source(160,144),scalar(170,150),batched(170,150),t0,t1,t2
    integer :: i,j,k,code,scalar_calls,batch_calls
    integer, parameter :: repeats=12

    do j=1,size(source,2)
        do i=1,size(source,1)
            source(i,j)=0.1_dp*i+sin(0.13_dp*i+0.21_dp*j)
        end do
    end do
    call prepare_spline_plan_nd(source,plan,order=1,mode='nearest',status=code)
    if (code/=0) error stop 'prepare plan'
    scalar_calls=0
    batch_calls=0
    call cpu_time(t0)
    do k=1,repeats
        call evaluate_spline_plan_geometric_nd(plan,point_map,scalar,code,batch_size=256)
        if (code/=0) error stop 'scalar mapping'
    end do
    call cpu_time(t1)
    do k=1,repeats
        call evaluate_spline_plan_geometric_batch_nd(plan,block_map,batched,code,batch_size=256)
        if (code/=0) error stop 'batch mapping'
    end do
    call cpu_time(t2)
    if (maxval(abs(scalar-batched))>2e-12_dp) error stop 'mismatch'
    print '(a,es12.5)', 'Scalar callback CPU seconds: ',t1-t0
    print '(a,es12.5)', 'Batched callback CPU seconds: ',t2-t1
    if (t2>t1 .and. t1>t0) print '(a,f9.3)', 'Scalar/batch CPU ratio: ',(t1-t0)/(t2-t1)
    print '(a,i0)', 'Scalar callback calls: ',scalar_calls
    print '(a,i0)', 'Batched callback calls: ',batch_calls
    print '(a,es12.5)', 'Maximum output difference: ',maxval(abs(scalar-batched))
contains
    subroutine point_map(out_coords,in_coords)
        integer, intent(in) :: out_coords(:) !! one-based output index
        real(dp), intent(out) :: in_coords(:) !! one-based image position
        scalar_calls=scalar_calls+1
        in_coords(1)=1.0_dp+0.88_dp*out_coords(1)-0.07_dp*out_coords(2) &
            +0.2_dp*sin(0.04_dp*out_coords(2))
        in_coords(2)=1.0_dp+0.08_dp*out_coords(1)+0.81_dp*out_coords(2) &
            -0.1_dp*cos(0.05_dp*out_coords(1))
    end subroutine point_map

    subroutine block_map(out_coords,in_coords)
        integer, intent(in) :: out_coords(:,:) !! destination indices by column
        real(dp), intent(out) :: in_coords(:,:) !! source coordinates by column
        batch_calls=batch_calls+1
        in_coords(1,:)=1.0_dp+0.88_dp*out_coords(1,:)-0.07_dp*out_coords(2,:) &
            +0.2_dp*sin(0.04_dp*out_coords(2,:))
        in_coords(2,:)=1.0_dp+0.08_dp*out_coords(1,:)+0.81_dp*out_coords(2,:) &
            -0.1_dp*cos(0.05_dp*out_coords(1,:))
    end subroutine block_map
end program bench_ndimage_spline_batch_mapping
