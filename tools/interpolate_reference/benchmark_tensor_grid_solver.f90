! SPDX-License-Identifier: MIT
! Benchmark direct 1D B-spline collocation solver (banded vs historical dense).
program benchmark_tensor_grid_solver
    use scifort_kinds, only : dp
    use scifort_interpolate_regular_grid, only : grid_axis, regular_grid_interpolator, &
        regular_grid_interpolator_type, grid_status_success
    implicit none
    type(grid_axis) :: axes(1)
    type(regular_grid_interpolator_type) :: grid
    real(dp), allocatable :: values(:)
    real(dp) :: elapsed, start, finish, max_error, x, y
    character(len=24) :: arg, method
    integer :: n, i, status

    n = 256
    method = 'quintic'
    call get_command_argument(1,arg)
    if (len_trim(arg)>0) read(arg,*) n
    call get_command_argument(2,arg)
    if (len_trim(arg)>0) method = arg
    if (n < 6) error stop 'Need at least six knots'
    allocate(axes(1)%points(n),values(n))
    do i=1,n
        axes(1)%points(i) = real(i-1,dp)/real(n-1,dp) + &
            0.0001_dp*sin(real(i,dp)) / real(n,dp)
    end do
    values = sin(8.0_dp*axes(1)%points)+axes(1)%points**2
    call cpu_time(start)
    grid = regular_grid_interpolator(axes, values, method=trim(method), status=status)
    call cpu_time(finish)
    if (status /= grid_status_success) error stop 'grid fit failed'
    elapsed = finish-start
    max_error = 0.0_dp
    do i=1,n
        x = axes(1)%points(i)
        y = grid%evaluate([x],status=status)
        if (status /= grid_status_success) error stop 'grid evaluation failed'
        max_error = max(max_error,abs(y-values(i)))
    end do
    print '(a,1x,a,1x,i0,1x,a,1x,es16.8e3,1x,a,1x,es16.8e3)', &
        'METHOD',trim(method),n,'CPU_SECONDS',elapsed,'MAX_NODE_ERROR',max_error
end program benchmark_tensor_grid_solver
