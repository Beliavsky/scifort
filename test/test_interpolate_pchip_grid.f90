! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_interpolate_pchip_grid
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_interpolate, only : grid_axis, regular_grid_interpolator, &
        regular_grid_interpolator_type, complex_regular_grid_interpolator, &
        complex_regular_grid_interpolator_type, interpn, &
        grid_status_success, grid_status_invalid_input, grid_status_out_of_bounds
    implicit none
    type(grid_axis) :: axes(2), short(1), high_axes(6)
    type(regular_grid_interpolator_type) :: grid, invalid, override
    type(complex_regular_grid_interpolator_type) :: cgrid
    real(dp) :: a(5), b(4), data(5,4), multi(5,4,2)
    real(dp) :: q(2), truth, got, query(3,2)
    real(dp), allocatable :: output(:,:), channels(:)
    integer :: i, j, k, code

    a = [-2.0_dp,-0.5_dp,0.7_dp,1.9_dp,3.5_dp]
    b = [9.0_dp,5.0_dp,1.0_dp,-2.0_dp]
    axes(1) = grid_axis(a)
    axes(2) = grid_axis(b)
    do j=1,size(b)
        do i=1,size(a)
            data(i,j) = 3.0_dp+2.0_dp*a(i)-4.0_dp*b(j)
        end do
    end do
    multi(:,:,1) = data
    multi(:,:,2) = 2.0_dp*data+0.5_dp
    grid = regular_grid_interpolator(axes,data,method='pchip',status=code)
    if (code /= grid_status_success .or. .not. grid%valid()) error stop 'PCHIP construction'
    q = [0.34_dp,2.1_dp]
    truth = 3.0_dp+2.0_dp*q(1)-4.0_dp*q(2)
    got = grid%evaluate(q,status=code)
    if (code /= grid_status_success .or. abs(got-truth)>2.0e-13_dp) &
        error stop 'PCHIP affine reproduction with descending axis'
    query(1,:) = q
    query(2,:) = [-0.9_dp,4.3_dp]
    query(3,:) = [2.0_dp,0.1_dp]
    output = interpn(axes,data,query,method='pchip',status=code)
    if (code /= 0) error stop 'PCHIP interpn status'
    if (maxval(abs(output(:,1)-(3.0_dp+2.0_dp*query(:,1)-4.0_dp*query(:,2)))) > 3.0e-13_dp) &
        error stop 'PCHIP interpn affine reproduction'

    grid = regular_grid_interpolator(axes,multi,method='pchip',status=code)
    if (code /= 0 .or. grid%channels() /= 2) error stop 'PCHIP multichannel setup'
    channels = grid%evaluate_channels(q,status=code)
    if (code /= 0 .or. maxval(abs(channels-[truth,2.0_dp*truth+0.5_dp])) > 3.0e-13_dp) &
        error stop 'PCHIP multichannel agreement'
    output = grid%evaluate(query,status=code)
    if (code /= 0 .or. any(shape(output) /= [3,2])) error stop 'PCHIP batch output shape'
    if (maxval(abs(output(:,2)-2.0_dp*output(:,1)-0.5_dp)) > 4.0e-13_dp) &
        error stop 'PCHIP batch multichannel agreement'
    channels = grid%evaluate_channels(q,nu=[0,0],status=code)
    if (code /= grid_status_invalid_input) error stop 'PCHIP derivative request must fail'
    channels = grid%evaluate_channels([ieee_value(0.0_dp,ieee_quiet_nan),2.0_dp],status=code)
    if (code /= grid_status_success .or. .not. all(ieee_is_nan(channels))) &
        error stop 'PCHIP NaN query'

    grid = regular_grid_interpolator(axes,data,method='pchip',bounds_error=.false., &
        fill_value=-17.0_dp,status=code)
    if (abs(grid%evaluate([-5.0_dp,1.2_dp])+17.0_dp)>0.0_dp) &
        error stop 'PCHIP out-of-range fill'
    grid = regular_grid_interpolator(axes,data,method='pchip',bounds_error=.false., &
        extrapolate=.true.,status=code)
    if (abs(grid%evaluate([-5.0_dp,1.2_dp])-(3.0_dp-10.0_dp-4.8_dp))>1.0e-12_dp) &
        error stop 'PCHIP affine extrapolation'
    grid = regular_grid_interpolator(axes,data,method='pchip',status=code)
    got = grid%evaluate([-5.0_dp,1.2_dp],status=code)
    if (code /= grid_status_out_of_bounds .or. .not. ieee_is_nan(got)) &
        error stop 'PCHIP bounds error'

    ! Evaluation-time method changes from linear to PCHIP and vice versa.
    override = regular_grid_interpolator(axes,data,method='linear',status=code)
    got = override%evaluate(q,method='pchip',status=code)
    if (code /= 0 .or. abs(got-truth)>2.0e-13_dp) error stop 'PCHIP method override'
    got = grid%evaluate(q,method='linear',status=code)
    if (code /= 0 .or. abs(got-truth)>2.0e-13_dp) error stop 'linear method override'
    short(1) = grid_axis([0.0_dp,1.0_dp,2.0_dp])
    invalid = regular_grid_interpolator(short,[1.0_dp,2.0_dp,3.0_dp],method='pchip',status=code)
    if (code /= grid_status_invalid_input .or. invalid%valid()) error stop 'PCHIP needs four knots'
    override = regular_grid_interpolator(short,[1.0_dp,2.0_dp,3.0_dp],method='linear',status=code)
    got = override%evaluate([0.5_dp],method='pchip',status=code)
    if (code /= grid_status_invalid_input .or. .not. ieee_is_nan(got)) &
        error stop 'PCHIP override requires four knots'
    cgrid = complex_regular_grid_interpolator(axes,cmplx(data,0.0_dp,kind=dp), &
        method='pchip',status=code)
    if (code /= grid_status_invalid_input .or. cgrid%valid()) error stop 'complex PCHIP must be rejected'
    cgrid = complex_regular_grid_interpolator(axes,cmplx(data,0.0_dp,kind=dp), &
        method='linear',status=code)
    block
        complex(dp), allocatable :: complex_out(:)
        complex_out = cgrid%evaluate_channels(q,method='pchip',status=code)
        if (code /= grid_status_invalid_input .or. .not. all(ieee_is_nan(real(complex_out,dp)))) &
            error stop 'complex PCHIP method override must be rejected'
    end block

    do k=1,6
        if (mod(k,2)==0) then
            high_axes(k) = grid_axis([3.0_dp,2.0_dp,0.3_dp,-2.0_dp])
        else
            high_axes(k) = grid_axis([-2.0_dp,0.3_dp,2.0_dp,3.0_dp])
        end if
    end do
    block
        real(dp) :: data6(4,4,4,4,4,4)
        data6 = 4.25_dp
        grid = regular_grid_interpolator(high_axes,data6,method='pchip',status=code)
        if (code /= 0) error stop '6D PCHIP setup'
        if (abs(grid%evaluate([1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp])-4.25_dp) > 2.0e-14_dp) &
            error stop '6D PCHIP constant reproduction'
    end block
    print '(a)', 'PCHIP tensor grid native regression passed'
end program test_interpolate_pchip_grid
