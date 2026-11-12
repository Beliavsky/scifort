! SPDX-License-Identifier: MIT
! Benchmark of repeated cubic spline interpolation on a fixed source array.
program bench_ndimage_spline_plan
    use scifort_kinds, only : dp
    use scifort_signal_spline_nd, only : spline_interpolation_plan, &
        prepare_spline_plan_nd, evaluate_spline_plan_nd, map_coordinates_nd
    implicit none
    integer, parameter :: nx=80,ny=88,npoints=64,repeats=32
    real(dp) :: source(nx,ny), coords(2,npoints), cached(npoints), direct(npoints)
    real(dp) :: t0,t1,t2,t3,prepared_seconds,unprepared_seconds,error
    type(spline_interpolation_plan) :: plan
    integer :: i,j,k,code

    do j=1,ny
        do i=1,nx
            source(i,j)=sin(0.013_dp*i+0.043_dp*j)+0.0004_dp*i*j
        end do
    end do
    do k=1,npoints
        coords(1,k)=1.0_dp+real(mod(17*k,nx),dp)+0.12_dp
        coords(2,k)=1.0_dp+real(mod(13*k,ny),dp)+0.27_dp
    end do
    call cpu_time(t0)
    call prepare_spline_plan_nd(source,plan,order=5,mode='grid-wrap',status=code)
    if (code/=0) error stop 'cannot prepare plan'
    do k=1,repeats
        call evaluate_spline_plan_nd(plan,coords,cached,status=code)
        if (code/=0) error stop 'cannot evaluate plan'
    end do
    call cpu_time(t1)
    do k=1,repeats
        call map_coordinates_nd(source,coords,direct,order=5,mode='grid-wrap', &
            prefilter=.true.,status=code)
        if (code/=0) error stop 'cannot interpolate'
    end do
    call cpu_time(t2)
    error=maxval(abs(cached-direct))
    if (error>1.0e-11_dp) error stop 'cached and direct results disagree'
    prepared_seconds=t1-t0
    unprepared_seconds=t2-t1
    write(*,'(a,es13.5)') 'Cached preparation + queries CPU seconds: ', prepared_seconds
    write(*,'(a,es13.5)') 'Repeated direct-map CPU seconds:         ', unprepared_seconds
    if (prepared_seconds>0.0_dp) then
        write(*,'(a,f9.2)') 'Observed speedup: ', unprepared_seconds/prepared_seconds
    end if
    write(*,'(a,es13.5)') 'Maximum absolute numerical difference: ',error
end program bench_ndimage_spline_plan
