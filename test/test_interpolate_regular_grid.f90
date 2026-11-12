! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_interpolate_regular_grid
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_interpolate, only : grid_axis, regular_grid_interpolator_type, &
        complex_regular_grid_interpolator_type, regular_grid_interpolator, complex_regular_grid_interpolator, &
        grid_status_success, grid_status_invalid_input, grid_status_out_of_bounds, interpn, complex_interpn
    implicit none
    type(grid_axis) :: axes(2), six(6), invalid_axis(1)
    type(regular_grid_interpolator_type) :: grid, invalid, high
    type(complex_regular_grid_interpolator_type) :: cgrid
    real(dp) :: y(3,4), yy(3,4,2), x(2), q(3,2), data6(2,1,2,2,2,2)
    real(dp), allocatable :: vals(:,:), channels(:), values(:)
    type(grid_axis) :: one_axis(1)
    complex(dp) :: z(3,4), zz(3,4,2)
    complex(dp), allocatable :: cvals(:,:)
    integer :: i, j, code

    axes(1) = grid_axis([0.0_dp,1.0_dp,3.0_dp])
    axes(2) = grid_axis([7.0_dp,4.0_dp,1.0_dp,-2.0_dp]) ! descending irregular grids
    do j = 1, 4
        do i = 1, 3
            y(i,j) = 5.0_dp+2.0_dp*axes(1)%points(i)-3.0_dp*axes(2)%points(j)
        end do
    end do
    grid = regular_grid_interpolator(axes,y,status=code)
    if (code /= grid_status_success .or. .not. grid%valid() .or. grid%channels() /= 1) &
        error stop 'grid construction failed'
    x = [1.5_dp,2.75_dp]
    if (abs(grid%evaluate(x)-(5.0_dp+2.0_dp*x(1)-3.0_dp*x(2))) > 5.0e-14_dp) &
        error stop 'multilinear affine exactness'
    q(1,:) = x
    q(2,:) = [0.5_dp,3.0_dp]
    q(3,:) = [2.0_dp,-0.5_dp]
    values = grid%evaluate_values(q,status=code)
    if (code /= 0 .or. maxval(abs(values - (5.0_dp+2.0_dp*q(:,1)-3.0_dp*q(:,2)))) > 4.0e-14_dp) &
        error stop 'batch linear exactness'
    vals = interpn(axes,y,q,status=code)
    if (code /= 0 .or. maxval(abs(vals(:,1)-values)) > 1.0e-15_dp) &
        error stop 'interpn helper disagrees'

    if (abs(grid%evaluate([1.5_dp,5.5_dp],method='nearest')-y(2,2)) > 0.0_dp) &
        error stop 'descending nearest tie chooses smaller physical coordinate'
    if (.not. ieee_is_nan(grid%evaluate([-1.0_dp,0.0_dp],status=code))) &
        error stop 'invalid outside evaluation should be nan'
    if (code /= grid_status_out_of_bounds) error stop 'bounds error status'

    grid = regular_grid_interpolator(axes,y,bounds_error=.false.,fill_value=-432.0_dp,status=code)
    if (code /= 0 .or. abs(grid%evaluate([-1.0_dp,0.0_dp])+432.0_dp) > 0.0_dp) &
        error stop 'explicit outside fill'
    grid = regular_grid_interpolator(axes,y,bounds_error=.false.,extrapolate=.true.,status=code)
    if (abs(grid%evaluate([-1.0_dp,8.0_dp])-(5.0_dp-2.0_dp-24.0_dp)) > 6.0e-14_dp) &
        error stop 'linear extrapolation'

    yy(:,:,1) = y
    yy(:,:,2) = -2.0_dp*y+0.5_dp
    grid = regular_grid_interpolator(axes,yy,status=code)
    if (code /= 0 .or. grid%channels() /= 2) error stop 'multichannel grid construction'
    channels = grid%evaluate_channels(x,status=code)
    if (code /= 0 .or. maxval(abs(channels-[5.0_dp+2.0_dp*x(1)-3.0_dp*x(2), &
        -2.0_dp*(5.0_dp+2.0_dp*x(1)-3.0_dp*x(2))+0.5_dp])) > 1.0e-13_dp) &
        error stop 'multichannel evaluation'
    vals = grid%evaluate(q)
    if (any(shape(vals) /= [3,2])) error stop 'multichannel batch shape'
    if (maxval(abs(vals(:,2)-(-2.0_dp*vals(:,1)+0.5_dp))) > 1.0e-13_dp) &
        error stop 'multichannel values'
    channels = grid%evaluate_channels([ieee_value(0.0_dp,ieee_quiet_nan),1.0_dp],status=code)
    if (code /= 0 .or. .not. all(ieee_is_nan(channels))) error stop 'nan query'
    invalid = regular_grid_interpolator(axes,y,method='invalid',status=code)
    if (code /= grid_status_invalid_input .or. invalid%valid()) error stop 'unknown method validation'
    invalid_axis(1) = grid_axis([1.0_dp,1.0_dp])
    invalid = regular_grid_interpolator(invalid_axis,[3.0_dp,4.0_dp],status=code)
    if (code /= grid_status_invalid_input .or. invalid%valid()) error stop 'duplicate knots validation'
    ! Match SciPy's linear NaN propagation even on a zero-weight neighbor.
    one_axis(1) = grid_axis([0.0_dp,1.0_dp])
    grid = regular_grid_interpolator(one_axis, &
        [1.0_dp,ieee_value(0.0_dp,ieee_quiet_nan)],status=code)
    if (code /= 0 .or. .not. ieee_is_nan(grid%evaluate([0.0_dp]))) &
        error stop 'zero-weight NaN neighbor propagation'

    z = cmplx(y,-0.75_dp*y+1.0_dp,kind=dp)
    zz(:,:,1) = z
    zz(:,:,2) = -0.25_dp*z
    cgrid = complex_regular_grid_interpolator(axes,zz,status=code)
    if (code /= 0 .or. .not. cgrid%valid() .or. cgrid%channels() /= 2) &
        error stop 'complex grid construction'
    cvals = cgrid%evaluate(q,status=code)
    if (code /= 0 .or. any(shape(cvals) /= [3,2])) error stop 'complex matrix shape'
    if (maxval(abs(cvals(:,2)+0.25_dp*cvals(:,1))) > 7.0e-14_dp) &
        error stop 'complex matrix interpolation'
    cvals = complex_interpn(axes,zz,q,status=code)
    if (code /= grid_status_success .or. &
        maxval(abs(cvals(:,2)+0.25_dp*cvals(:,1))) > 7.0e-14_dp) &
        error stop 'complex interpn helper'
    cgrid = complex_regular_grid_interpolator(axes,z,bounds_error=.false.,status=code)
    if (code /= 0) error stop 'single-channel complex grid constructor'
    if (.not. ieee_is_nan(real(cgrid%evaluate([-5.0_dp,0.0_dp]),dp))) &
        error stop 'complex default fill nan'
    if (abs(aimag(cgrid%evaluate([-5.0_dp,0.0_dp]))) > 0.0_dp) &
        error stop 'complex default imaginary fill zero'

    do i = 1, 6
        if (i == 2) then
            six(i) = grid_axis([5.0_dp])
        else
            six(i) = grid_axis([0.0_dp,2.0_dp])
        end if
    end do
    data6 = 2.5_dp
    high = regular_grid_interpolator(six,data6,bounds_error=.false.,extrapolate=.true.,status=code)
    if (code /= 0 .or. high%dimensions() /= 6 .or. &
        abs(high%evaluate([0.3_dp,12.0_dp,1.0_dp,-2.0_dp,4.0_dp,3.0_dp])-2.5_dp) > 1.0e-14_dp) &
        error stop 'six-dimensional singleton extrapolation'
    print '(a)', 'regular_grid_interpolator native regression passed'
end program test_interpolate_regular_grid
