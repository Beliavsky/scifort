! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_savgol
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_signal
    implicit none

    real(dp), allocatable :: coeff(:), y(:)
    real(dp) :: x(21), cubic(21)
    integer :: i, status

    coeff = savgol_coeffs(5, 2, status=status)
    if (status /= savgol_status_success) error stop 'coefficients status'
    call assert_close(coeff, [-3.0_dp, 12.0_dp, 17.0_dp, 12.0_dp, -3.0_dp] / 35.0_dp)
    coeff = savgol_coeffs(5, 2, deriv=1, use='conv')
    call assert_close(coeff, [2.0_dp, 1.0_dp, 0.0_dp, -1.0_dp, -2.0_dp] / 10.0_dp)
    coeff = savgol_coeffs(5, 2, deriv=1, use='dot')
    call assert_close(coeff, [-2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp] / 10.0_dp)
    coeff = savgol_coeffs(4, 1)
    call assert_close(coeff, [0.25_dp, 0.25_dp, 0.25_dp, 0.25_dp])
    coeff = savgol_coeffs(5, 2, deriv=3)
    call assert_close(coeff, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp])
    coeff = savgol_coeffs(4, 3, pos=1, use='dot')
    call assert_close(coeff, [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp])

    do i = 1, size(x)
        x(i) = real(i - 1, dp)
        cubic(i) = 5.0_dp + 2.0_dp * x(i) - 0.25_dp * x(i)**2 + 0.125_dp * x(i)**3
    end do
    y = savgol_filter(cubic, 7, 3)
    call assert_close(y, cubic, 2.0e-10_dp)
    y = savgol_filter(cubic, 7, 3, deriv=1, delta=1.0_dp)
    call assert_close(y, 2.0_dp - 0.5_dp * x + 0.375_dp * x**2, 2.0e-10_dp)
    y = savgol_filter(cubic, 7, 3, deriv=2)
    call assert_close(y, -0.5_dp + 0.75_dp * x, 2.0e-10_dp)
    y = savgol_filter(cubic, 7, 3, deriv=3)
    call assert_close(y, spread(0.75_dp, 1, 21), 2.0e-10_dp)

    ! Even-length windows evaluate the interior halfway between data samples.
    y = savgol_filter(x, 4, 1, mode='nearest')
    call assert_close(y(3:19), x(3:19) + 0.5_dp)
    y = savgol_filter(x, 4, 1, mode='interp')
    call assert_close(y(1:2), x(1:2))
    call assert_close(y(20:21), x(20:21))
    call assert_close(y(3:19), x(3:19) + 0.5_dp)

    ! All boundary policies must produce finite, properly sized output.
    do i = 1, 4
        select case (i)
        case (1)
            y = savgol_filter(x, 7, 2, mode='mirror', status=status)
        case (2)
            y = savgol_filter(x, 7, 2, mode='wrap', status=status)
        case (3)
            y = savgol_filter(x, 7, 2, mode='nearest', status=status)
        case (4)
            y = savgol_filter(x, 7, 2, mode='constant', cval=11.0_dp, status=status)
        end select
        if (status /= savgol_status_success) error stop 'boundary status'
        if (size(y) /= size(x)) error stop 'boundary length'
        if (any(ieee_is_nan(y))) error stop 'boundary NaN'
    end do

    coeff = savgol_coeffs(0, 0, status=status)
    if (status /= savgol_status_invalid_input) error stop 'zero window status'
    coeff = savgol_coeffs(5, 5, status=status)
    if (status /= savgol_status_invalid_input) error stop 'invalid order status'
    if (.not. all(ieee_is_nan(coeff))) error stop 'invalid order coefficient'
    coeff = savgol_coeffs(5, 2, pos=0, status=status)
    if (status /= savgol_status_invalid_input) error stop 'invalid position status'
    coeff = savgol_coeffs(5, 2, use='bad', status=status)
    if (status /= savgol_status_invalid_input) error stop 'invalid coefficient ordering'
    y = savgol_filter(x, 5, 2, mode='unknown', status=status)
    if (status /= savgol_status_invalid_input) error stop 'invalid boundary status'
    if (.not. all(ieee_is_nan(y))) error stop 'invalid output values'
    y = savgol_filter(x(1:3), 5, 2, status=status)
    if (status /= savgol_status_invalid_input) error stop 'short interpolation window'
    y = savgol_filter(x(1:3), 5, 2, mode='nearest', status=status)
    if (status /= savgol_status_success) error stop 'short nearest window'
    if (size(y) /= 3) error stop 'short nearest length'
    print '(a)', 'test_signal_savgol passed'
contains
    subroutine assert_close(a, b, tolerance)
        real(dp), intent(in) :: a(:) !! actual values
        real(dp), intent(in) :: b(:) !! expected values
        real(dp), intent(in), optional :: tolerance !! allowable absolute difference
        real(dp) :: tol
        tol = 2.0e-12_dp
        if (present(tolerance)) tol = tolerance
        if (size(a) /= size(b)) error stop 'savgol length mismatch'
        if (any(ieee_is_nan(a)) .or. any(abs(a - b) > tol)) then
            print *, a
            print *, b
            error stop 'savgol numerical mismatch'
        end if
    end subroutine assert_close
end program test_signal_savgol
