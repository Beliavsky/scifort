! SPDX-License-Identifier: MIT
! Compare direct tensor-product queries with successive 1-D contractions.
program bench_ndimage_spline_tensor
    use scifort_kinds, only : dp
    use scifort_signal_spline_nd, only : spline_interpolation_plan, spline_axis_grid, &
        prepare_spline_plan_nd, evaluate_spline_plan_nd, evaluate_spline_plan_tensor_nd
    implicit none
    integer, parameter :: nx=180,ny=160,mx=210,my=180,repeats=4
    real(dp) :: source(nx,ny), coords(2,mx*my), direct(mx*my), tensor(mx,my)
    type(spline_interpolation_plan) :: plan
    type(spline_axis_grid) :: grid(2)
    real(dp) :: t0,t1,t2,err
    integer :: i,j,k,code

    do j=1,ny
        do i=1,nx
            source(i,j)=sin(0.014_dp*i+0.029_dp*j)+0.0001_dp*i*j
        end do
    end do
    allocate(grid(1)%coordinates(mx),grid(2)%coordinates(my))
    do i=1,mx
        grid(1)%coordinates(i)=real(i,dp)*0.84_dp-1.3_dp
    end do
    do j=1,my
        grid(2)%coordinates(j)=real(j,dp)*0.88_dp-0.7_dp
    end do
    do j=1,my
        do i=1,mx
            k=i+mx*(j-1)
            coords(:,k)=[grid(1)%coordinates(i),grid(2)%coordinates(j)]
        end do
    end do
    call prepare_spline_plan_nd(source,plan,order=5,mode='mirror',status=code)
    if (code/=0) error stop 'cannot prepare plan'
    call cpu_time(t0)
    do k=1,repeats
        call evaluate_spline_plan_nd(plan,coords,direct,code)
        if (code/=0) error stop 'pointwise interpolation failed'
    end do
    call cpu_time(t1)
    do k=1,repeats
        call evaluate_spline_plan_tensor_nd(plan,grid,tensor,code)
        if (code/=0) error stop 'tensor interpolation failed'
    end do
    call cpu_time(t2)
    err=maxval(abs(reshape(direct,shape(tensor))-tensor))
    if (err>1.0e-10_dp) error stop 'tensor results disagree'
    print '(a,es12.4)', 'Pointwise plan seconds: ',t1-t0
    print '(a,es12.4)', 'Tensor plan seconds:    ',t2-t1
    if (t2>t1) print '(a,f8.2,a)', 'Measured speedup: ', &
        (t1-t0)/(t2-t1),'x'
    print '(a,es12.4)', 'Maximum difference:    ',err
end program bench_ndimage_spline_tensor
