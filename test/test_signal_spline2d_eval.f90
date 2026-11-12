! SPDX-License-Identifier: MIT
! Tests public 2-D quadratic and cubic spline coefficient evaluators.
program test_signal_spline2d_eval
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_signal, only : dp, cspline2d_eval, qspline2d_eval, &
        cspline2d_grid, qspline2d_grid, spline_status_success, spline_status_invalid_input
    implicit none
    real(dp) :: coeff(13,17), rows(4,5), cols(4,5), nanval
    real(dp), allocatable :: values(:,:), grid(:,:), vector(:), sub(:,:)
    complex(dp), allocatable :: cv(:,:), cgrid(:,:), cvector(:)
    complex(dp) :: z
    integer :: i, j, status

    do j = 1, 17
        do i = 1, 13
            coeff(i,j) = sin(0.21_dp*real(i,dp)) + cos(0.13_dp*real(j,dp))
        end do
    end do
    do j = 1, 5
        do i = 1, 4
            rows(i,j) = 0.31_dp*real(i,dp) + 7.0_dp*real(j-2,dp)
            cols(i,j) = 1.27_dp*real(j,dp) - 4.25_dp*real(i,dp)
        end do
    end do
    values = cspline2d_eval(coeff, rows, cols, status=status)
    if (status /= spline_status_success .or. any(shape(values) /= [4,5])) &
        error stop 'cspline2d_eval array shape/status'
    do j = 1, 5
        do i = 1, 4
            if (abs(values(i,j) - cspline2d_eval(coeff, rows(i,j), cols(i,j))) > 2.e-14_dp) &
                error stop 'cspline2d_eval scalar/array mismatch'
        end do
    end do
    vector = cspline2d_eval(coeff, rows(:,2), cols(:,2), status=status)
    if (status /= spline_status_success .or. size(vector) /= 4) &
        error stop 'cspline2d_eval vector shape'
    if (maxval(abs(vector-values(:,2))) > 2.e-14_dp) error stop 'vector/array mismatch'

    grid = cspline2d_grid(coeff, rows(:,2), cols(3,:), status=status)
    if (status /= spline_status_success .or. any(shape(grid) /= [4,5])) &
        error stop 'cspline2d_grid shape/status'
    do j = 1, 5
        do i = 1, 4
            if (abs(grid(i,j)-cspline2d_eval(coeff,rows(i,2),cols(3,j))) > 2.e-14_dp) &
                error stop 'cspline2d_grid values'
        end do
    end do

    sub = qspline2d_eval(coeff(1:13:2,1:17:2), rows, cols, status=status)
    if (status /= spline_status_success) error stop 'noncontiguous input'
    values = coeff(1:13:2,1:17:2)
    values = qspline2d_eval(values, rows, cols)
    if (maxval(abs(sub-values)) > 3.e-14_dp) error stop 'noncontiguous coeff mismatch'

    cv = cmplx(coeff, -0.7_dp*coeff, dp)
    cgrid = qspline2d_grid(cv, rows(:,1), cols(1,:), status=status)
    if (status /= spline_status_success) error stop 'complex quadratic grid'
    grid = qspline2d_grid(coeff, rows(:,1), cols(1,:))
    if (maxval(abs(cgrid - cmplx(grid,-0.7_dp*grid,dp))) > 1.e-13_dp) &
        error stop 'complex quadratic linearity'
    cvector = cspline2d_eval(cv, rows(:,1), cols(:,1))
    vector = cspline2d_eval(coeff, rows(:,1), cols(:,1))
    if (maxval(abs(cvector - cmplx(vector,-0.7_dp*vector,dp))) > 1.e-13_dp) &
        error stop 'complex cubic linearity'
    z = cspline2d_eval(cv, -8.5_dp, 19.125_dp, status=status)
    if (status /= spline_status_success) error stop 'complex scalar'
    if (abs(z-cmplx(cspline2d_eval(coeff,-8.5_dp,19.125_dp), &
        -0.7_dp*cspline2d_eval(coeff,-8.5_dp,19.125_dp),dp)) > 1.e-13_dp) &
        error stop 'complex scalar linearity'

    ! Constant-coefficient images must be constant at every mirrored coordinate.
    values = cspline2d_eval(reshape([3.0_dp],[1,1]), rows, cols, status=status)
    if (status /= spline_status_success .or. &
        maxval(abs(values-3.0_dp)) > 1.e-13_dp) error stop 'one-pixel cubic'
    values = qspline2d_eval(reshape([3.0_dp, 3.0_dp],[1,2]), rows, cols)
    if (maxval(abs(values-3.0_dp)) > 1.e-13_dp) error stop 'one-row quadratic'
    values = qspline2d_eval(reshape([3.0_dp],[1,1]), rows, cols)
    if (maxval(abs(values-3.0_dp)) > 1.e-13_dp) error stop 'one-pixel quadratic'

    nanval = ieee_value(0.0_dp, ieee_quiet_nan)
    vector = cspline2d_eval(coeff, [1.0_dp,nanval], [1.0_dp,2.0_dp], status=status)
    if (status /= spline_status_invalid_input .or. size(vector) /= 0) &
        error stop 'NaN coordinate validation'
    values = qspline2d_grid(coeff, [1.0_dp,2.0_dp], [1.0_dp,nanval], status=status)
    if (status /= spline_status_invalid_input .or. size(values) /= 0) &
        error stop 'NaN grid validation'
    grid = cspline2d_eval(coeff(1:0,:), rows, cols, status=status)
    if (status /= spline_status_invalid_input .or. size(grid) /= 0) &
        error stop 'empty coeff validation'
    grid = cspline2d_eval(coeff, rows, cols(1:2,:), status=status)
    if (status /= spline_status_invalid_input .or. size(grid) /= 0) &
        error stop 'mismatched coord shape validation'
    grid = qspline2d_grid(coeff, [real(dp)::], [1.0_dp])
    if (size(grid,1) /= 0 .or. size(grid,2) /= 1) &
        error stop 'empty grid axis'
    print '(a)', 'PASS: 2-D spline evaluation scalar/vector/matrix/grid regression'
end program test_signal_spline2d_eval
