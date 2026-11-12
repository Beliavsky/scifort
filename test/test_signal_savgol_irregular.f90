! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_savgol_irregular
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_savgol, only : savgol_irregular, savgol_filter, &
        savgol_status_success, savgol_status_invalid_input, savgol_status_singular
    implicit none
    integer, parameter :: n = 31
    real(dp) :: x(n), y(n), w(n), q(17), bad(n), uniform(n)
    real(dp), allocatable :: result(:), benchmark(:)
    complex(dp) :: yc(n)
    complex(dp), allocatable :: complex_result(:)
    integer :: i, stat

    do i = 1, n
        x(i) = 0.12_dp * real(i,dp) + 0.001_dp * real(i,dp)**2
        uniform(i) = real(i-1,dp)
    end do
    y = 2.0_dp - 0.5_dp*x + 1.25_dp*x**2 - 0.2_dp*x**3
    yc = cmplx(y, -2.0_dp*y + 0.75_dp*x, dp)
    w = 1.0_dp + 0.03_dp * real([(i,i=1,n)],dp)
    q = [(x(2) + real(i-2,dp)*0.19_dp, i=1,17)]

    result = savgol_irregular(x,y,7,3,query=q,weights=w,status=stat)
    if (stat /= savgol_status_success) error stop 'real polynomial status'
    call assert_real(result,2.0_dp-0.5_dp*q+1.25_dp*q**2-0.2_dp*q**3,2.0e-11_dp)
    result = savgol_irregular(x,y,7,3,deriv=1,query=q,weights=w,status=stat)
    if (stat /= savgol_status_success) error stop 'real first derivative status'
    call assert_real(result,-0.5_dp+2.5_dp*q-0.6_dp*q**2,2.0e-11_dp)
    result = savgol_irregular(x,y,7,3,deriv=2,query=q,status=stat)
    if (stat /= savgol_status_success) error stop 'real second derivative status'
    call assert_real(result,2.5_dp-1.2_dp*q,1.0e-10_dp)
    result = savgol_irregular(x,y,7,3,deriv=4,status=stat)
    if (stat /= savgol_status_success) error stop 'above polynomial degree'
    if (maxval(abs(result)) > 1.0e-12_dp) error stop 'degree+1 derivative must vanish'
    complex_result = savgol_irregular(x,yc,7,3,query=q,weights=w,status=stat)
    if (stat /= savgol_status_success) error stop 'complex polynomial status'
    if (maxval(abs(real(complex_result,dp) - (2.0_dp-0.5_dp*q+1.25_dp*q**2-0.2_dp*q**3))) > &
        3.0e-11_dp) error stop 'complex real part'
    if (maxval(abs(aimag(complex_result) - (-2.0_dp*(2.0_dp-0.5_dp*q+1.25_dp*q**2-0.2_dp*q**3) &
        + 0.75_dp*q))) > 3.0e-11_dp) error stop 'complex imaginary part'

    ! Selecting nearest neighbors from strided inputs and returning derivatives.
    result = savgol_irregular(x(1:n:2),y(1:n:2),7,3,deriv=1,query=q,status=stat)
    if (stat /= savgol_status_success) error stop 'strided real status'
    call assert_real(result,-0.5_dp+2.5_dp*q-0.6_dp*q**2,5.0e-11_dp)
    complex_result = savgol_irregular(x(1:n:2),yc(1:n:2),7,3,deriv=2,query=q,status=stat)
    if (stat /= savgol_status_success) error stop 'strided complex status'
    if (maxval(abs(complex_result - cmplx(2.5_dp-1.2_dp*q, -2.0_dp*(2.5_dp-1.2_dp*q), dp))) > &
        2.0e-10_dp) error stop 'strided complex second derivative'

    ! Uniform-grid odd-window equivalence to scipy.signal.savgol_filter interp mode.
    benchmark = savgol_filter(sin(uniform/6.0_dp),9,3,mode='interp')
    result = savgol_irregular(uniform,sin(uniform/6.0_dp),9,3,status=stat)
    if (stat /= savgol_status_success) error stop 'uniform grid status'
    call assert_real(result,benchmark,1.0e-12_dp)

    bad = x
    bad(3) = bad(2)
    result = savgol_irregular(bad,y,7,3,status=stat)
    if (stat /= savgol_status_invalid_input .or. .not. all(ieee_is_nan(result))) &
        error stop 'reject duplicate timestamps'
    result = savgol_irregular(x,y,7,3,weights=-w,status=stat)
    if (stat /= savgol_status_invalid_input .or. .not. all(ieee_is_nan(result))) &
        error stop 'reject negative weights'
    result = savgol_irregular(x,y,7,3,weights=w(:n-1),status=stat)
    if (stat /= savgol_status_invalid_input) error stop 'reject wrong weights shape'
    result = savgol_irregular(x,y,7,3,weights=0.0_dp*w,status=stat)
    if (stat /= savgol_status_singular .or. .not. all(ieee_is_nan(result))) &
        error stop 'reject singular local polynomial'
    result = savgol_irregular(x,y,7,7,status=stat)
    if (stat /= savgol_status_invalid_input) error stop 'reject invalid polynomial degree'
    bad = y
    bad(6) = ieee_value(0.0_dp,ieee_quiet_nan)
    result = savgol_irregular(x,bad,7,3,status=stat)
    if (stat /= savgol_status_invalid_input) error stop 'reject nonfinite samples'
    print *, 'Irregular Savitzky-Golay native tests passed.'
contains
    subroutine assert_real(a,b,tolerance)
        real(dp), intent(in) :: a(:) !! numerical result
        real(dp), intent(in) :: b(:) !! expected values
        real(dp), intent(in) :: tolerance !! maximum absolute difference
        if (size(a) /= size(b)) error stop 'shape mismatch'
        if (maxval(abs(a-b)) > tolerance) error stop 'numeric mismatch'
    end subroutine assert_real
end program test_signal_savgol_irregular
