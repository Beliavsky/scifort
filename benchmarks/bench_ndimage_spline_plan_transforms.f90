! SPDX-License-Identifier: MIT
! CPU benchmark for repeated nonseparable transforms on a fixed source.
program bench_ndimage_spline_plan_transforms
    use scifort_kinds, only : dp
    use scifort_signal_spline_nd, only : spline_interpolation_plan, &
        prepare_spline_plan_nd, evaluate_spline_plan_affine_nd, affine_transform_nd
    implicit none
    integer, parameter :: repeat_count=12
    real(dp) :: x(240,180),y(72,64),reference(72,64)
    real(dp) :: matrix(2,2),offset(2),t0,t1,t2,t3
    type(spline_interpolation_plan) :: plan
    integer :: i,j,n,status

    do j=1,180
        do i=1,240
            x(i,j)=sin(real(i,dp)*0.02_dp)+cos(real(j,dp)*0.04_dp)+0.00002_dp*i*j
        end do
    end do
    matrix=reshape([1.03_dp,0.12_dp,-0.23_dp,0.98_dp],[2,2])
    offset=[55.13_dp,35.27_dp]
    call cpu_time(t0)
    do n=1,repeat_count
        call affine_transform_nd(x,matrix,offset,reference,order=5, &
            prefilter=.true.,mode='mirror',status=status)
        if (status/=0) error stop 'direct affine failed'
    end do
    call cpu_time(t1)
    call prepare_spline_plan_nd(x,plan,order=5,mode='mirror',status=status)
    if (status/=0) error stop 'plan preparation failed'
    call cpu_time(t2)
    do n=1,repeat_count
        call evaluate_spline_plan_affine_nd(plan,matrix,offset,y,status)
        if (status/=0) error stop 'plan affine failed'
    end do
    call cpu_time(t3)
    if (maxval(abs(y-reference))>1e-10_dp) error stop 'benchmark disagreement'
    print '(a,i0)', 'iterations: ',repeat_count
    print '(a,es15.7)', 'direct seconds: ',t1-t0
    print '(a,es15.7)', 'planned total seconds: ',t3-t1
    print '(a,es15.7)', 'plan setup seconds: ',t2-t1
    print '(a,es15.7)', 'maximum absolute error: ',maxval(abs(y-reference))
    print '(a,f8.3)', 'direct / planned: ',(t1-t0)/max(t3-t1,tiny(1.0_dp))
end program bench_ndimage_spline_plan_transforms
