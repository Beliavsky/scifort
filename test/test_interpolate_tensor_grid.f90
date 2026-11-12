! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_interpolate_tensor_grid
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use scifort_kinds, only: dp
    use scifort_interpolate, only: grid_axis, regular_grid_interpolator, &
        complex_regular_grid_interpolator, regular_grid_interpolator_type, &
        complex_regular_grid_interpolator_type, interpn, grid_status_success, &
        grid_status_invalid_input
    implicit none
    type(grid_axis) :: a(2), b(2), one(1), six(6)
    type(regular_grid_interpolator_type) :: grid, bad
    type(complex_regular_grid_interpolator_type) :: complex_grid
    real(dp) :: sample(7,6), affine(4,3), xv(7), yv(6), val, x(2), queries(2,2)
    real(dp), allocatable :: batch(:,:)
    complex(dp) :: z(7,6,2), expected
    integer :: i, j, code

    xv = [-2.0_dp, -0.6_dp, 0.1_dp, 0.8_dp, 1.5_dp, 2.3_dp, 3.6_dp]
    yv = [4.0_dp, 3.2_dp, 2.1_dp, 0.7_dp, -0.2_dp, -1.0_dp]
    b(1) = grid_axis(xv)
    b(2) = grid_axis(yv)
    do j = 1, 6
        do i = 1, 7
            sample(i,j) = xv(i)**4 + 0.1_dp*xv(i)*yv(j)**3 + yv(j)**2
        end do
    end do
    grid = regular_grid_interpolator(b, sample, method='quintic', status=code)
    if (.not. grid%valid() .or. code /= grid_status_success) error stop 'quintic construction'
    x = [0.34_dp, 1.14_dp]
    val = grid%evaluate(x, status=code)
    if (code /= 0) error stop 'quintic evaluation status'
    if (abs(val - (x(1)**4+0.1_dp*x(1)*x(2)**3+x(2)**2)) > 2.0e-11_dp) &
        error stop 'quintic polynomial reproduction'
    queries(1,:) = x
    queries(2,:) = [1.3_dp, -0.4_dp]
    batch = interpn(b, sample, queries, method='quintic', status=code)
    if (code /= 0) error stop 'quintic interpn status'
    if (abs(batch(1,1) - val) > 1.0e-12_dp) error stop 'quintic interpn agreement'
    val = grid%evaluate(x, nu=[1, 0])
    if (abs(val-(4.0_dp*x(1)**3+0.1_dp*x(2)**3)) > 2.0e-11_dp) &
        error stop 'quintic first derivative'
    val = grid%evaluate(x, nu=[0, 2])
    if (abs(val - (0.6_dp*x(1)*x(2)+2.0_dp)) > 2.0e-10_dp) &
        error stop 'quintic second partial'
    val = grid%evaluate(x, nu=[1, 3])
    if (abs(val-0.6_dp) > 3.0e-11_dp) error stop 'quintic mixed partial'
    val = grid%evaluate(x, nu=[6, 0])
    if (abs(val) > 1.0e-14_dp) error stop 'higher than degree derivative must vanish'
    val = grid%evaluate(x, nu=[-1, 0], status=code)
    if (code /= grid_status_invalid_input .or. .not. ieee_is_nan(val)) &
        error stop 'invalid negative derivative order'
    val = grid%evaluate(x, nu=[1], status=code)
    if (code /= grid_status_invalid_input) error stop 'invalid derivative rank'

    ! The original values are retained for fast per-call method overrides.
    a(1) = grid_axis([-1.0_dp, 0.0_dp, 2.0_dp, 4.0_dp])
    a(2) = grid_axis([5.0_dp, 2.0_dp, -0.5_dp])
    do j = 1,3
        do i = 1,4
            affine(i,j) = 1.5_dp+2.0_dp*a(1)%points(i)-3.0_dp*a(2)%points(j)
        end do
    end do
    grid = regular_grid_interpolator(a, affine, method='linear', status=code)
    if (code /= grid_status_success) error stop 'multilinear construction'
    x = [1.1_dp, 0.75_dp]
    val = grid%evaluate(x,method='slinear',nu=[0,1],status=code)
    if (code /= 0 .or. abs(val+3.0_dp) > 5.0e-14_dp) &
        error stop 'slinear method override with derivative'
    queries(1,:) = x
    queries(2,:) = [-0.7_dp, 3.4_dp]
    batch = grid%evaluate(queries,method='slinear',nu=[1,0],status=code)
    if (code /= 0 .or. maxval(abs(batch(:,1)-2.0_dp)) > 1.0e-13_dp) &
        error stop 'slinear batch method override derivative'
    val = grid%evaluate(x,method='slinear',status=code)
    if (code /= 0 .or. abs(val-(1.5_dp+2*x(1)-3*x(2))) > 5.0e-14_dp) &
        error stop 'slinear method override function'
    grid = regular_grid_interpolator(a, affine, method='cubic', status=code)
    if (code /= grid_status_invalid_input .or. grid%valid()) &
        error stop 'cubic requires four nodes on every axis'
    grid = regular_grid_interpolator(a, affine, method='slinear', status=code)
    if (code /= 0 .or. .not. grid%valid()) error stop 'slinear construction'
    val = grid%evaluate(x,method='linear',status=code)
    if (code /= 0 .or. abs(val-(1.5_dp+2*x(1)-3*x(2))) > 1.0e-14_dp) &
        error stop 'linear method override from spline'
    val = grid%evaluate(x, method='nearest', nu=[0,0], status=code)
    if (code /= grid_status_invalid_input) error stop 'nearest derivative rejection'
    val = grid%evaluate(x, method='quintic', status=code)
    if (code /= grid_status_invalid_input) error stop 'invalid quintic override size'

    ! Complex data and explicit fill values must survive the two-component fit.
    do j=1,6
        do i=1,7
            z(i,j,1) = cmplx(sample(i,j), -2.0_dp*sample(i,j), kind=dp)
            z(i,j,2) = cmplx(0.5_dp*sample(i,j), sample(i,j), kind=dp)
        end do
    end do
    complex_grid = complex_regular_grid_interpolator(b,z,method='quintic', &
        bounds_error=.false.,fill_value=cmplx(8.0_dp,3.0_dp,kind=dp),status=code)
    if (code /= 0 .or. complex_grid%channels() /= 2) error stop 'complex quintic constructor'
    expected = cmplx(4.0_dp*0.34_dp**3+0.1_dp*1.14_dp**3, &
        -2.0_dp*(4.0_dp*0.34_dp**3+0.1_dp*1.14_dp**3),kind=dp)
    block
        complex(dp), allocatable :: zv(:)
        zv = complex_grid%evaluate_channels([0.34_dp,1.14_dp],nu=[1,0])
        if (abs(zv(1)-expected) > 1.0e-10_dp) error stop 'complex quintic derivative'
    end block
    block
        complex(dp) :: bad_scalar
        bad_scalar = complex_grid%evaluate([10.0_dp,2.0_dp],status=code)
    end block
    if (code /= grid_status_invalid_input) error stop 'complex scalar invalid status'
    if (maxval(abs(complex_grid%evaluate_channels([10.0_dp,2.0_dp],status=code) - &
        cmplx(8.0_dp,3.0_dp,kind=dp))) > 0.0_dp) error stop 'complex fill value'
    if (code /= 0) error stop 'complex out of range fill status'

    ! Six-dimensional smallest possible cubic grid: 4**6 points.
    do i=1,6
        six(i) = grid_axis([0.0_dp,0.5_dp,1.3_dp,2.0_dp])
    end do
    block
        real(dp) :: c6(4,4,4,4,4,4)
        c6 = 1.25_dp
        grid = regular_grid_interpolator(six,c6,method='cubic',status=code)
        if (code /= 0 .or. .not. grid%valid()) error stop 'six-dimensional cubic construction'
        if (abs(grid%evaluate([1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp])-1.25_dp) > 1.0e-13_dp) &
            error stop 'six-dimensional cubic constant reproduction'
    end block
    do i=1,6
        if (mod(i,2) == 0) then
            six(i) = grid_axis([2.6_dp,1.7_dp,1.0_dp,0.2_dp,-0.7_dp,-1.5_dp])
        else
            six(i) = grid_axis([-1.5_dp,-0.7_dp,0.2_dp,1.0_dp,1.7_dp,2.6_dp])
        end if
    end do
    block
        real(dp), allocatable :: c6(:,:,:,:,:,:)
        allocate(c6(6,6,6,6,6,6))
        c6 = 1.25_dp
        grid = regular_grid_interpolator(six,c6,method='quintic',status=code)
        if (code /= 0 .or. .not. grid%valid()) error stop 'six-dimensional quintic construction'
        if (abs(grid%evaluate([1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp])-1.25_dp) > 1.0e-11_dp) &
            error stop 'six-dimensional quintic constant reproduction'
        if (abs(grid%evaluate([1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp], &
            nu=[0,0,1,0,0,0])) > 1.0e-11_dp) error stop 'six-dimensional quintic derivative'
    end block
    one(1) = grid_axis([0.0_dp,1.0_dp])
    bad = regular_grid_interpolator(one,[0.0_dp,1.0_dp],method='quintic',status=code)
    if (code /= grid_status_invalid_input .or. bad%valid()) error stop 'too-short quintic validation'
    print '(a)', 'tensor grid cubic/quintic/partial derivatives native regression passed'
end program test_interpolate_tensor_grid
