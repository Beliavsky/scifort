! SPDX-License-Identifier: MIT
! Tests the public 2-D B-spline interface and SciPy reference values.
program test_signal_spline2d
    use scifort_signal, only : dp, cspline2d, qspline2d, spline_filter, &
        spline_status_success, spline_status_invalid_input
    implicit none
    integer, parameter :: m = 77, n = 83
    integer, parameter :: indices(2,6) = reshape([ &
        1,1, 11,15, 35,40, 77,83, 1,43, 51,1], [2,6])
    real(dp), parameter :: refc(6) = [ &
        1.0862596456869642_dp, 1.6004822435930883_dp, &
        -0.19493069115583572_dp, 8.125342813496976_dp, &
        -0.8714969750878969_dp, 0.4143681022264074_dp]
    real(dp), parameter :: refq(6) = [ &
        1.0930934093883493_dp, 1.5999077056655495_dp, &
        -0.19440927301505084_dp, 8.116864830418635_dp, &
        -0.8614152368303225_dp, 0.41762928185360126_dp]
    real(dp), parameter :: refsmooth(6) = [ &
        1.1736853629651343_dp, 1.6008825603235324_dp, &
        -0.1943405984513996_dp, 7.995075559053377_dp, &
        -0.7505935746208111_dp, 0.45208319054676965_dp]
    real(dp), parameter :: reffilter(6) = [ &
        1.1809152594431185_dp, 1.5985400220136856_dp, &
        -0.19225845101273115_dp, 7.985577658180979_dp, &
        -0.7385307373521709_dp, 0.45621494231855136_dp]
    real(dp) :: x(m,n)
    real(dp), allocatable :: a(:,:), b(:,:), q(:,:), s(:,:), f(:,:)
    complex(dp), allocatable :: z(:,:), zf(:,:)
    integer :: i, j, status

    do j = 1, n
        do i = 1, m
            x(i,j) = sin(real(i,dp)*0.11_dp) + cos(real(j,dp)*0.07_dp) + &
                0.001_dp*real(i,dp)*real(j,dp)
        end do
    end do

    a = cspline2d(x, status=status)
    if (status /= spline_status_success) error stop 'cspline2d: status'
    call check_values(a, refc, 3.0e-6_dp)
    q = qspline2d(x, status=status)
    if (status /= spline_status_success) error stop 'qspline2d: status'
    call check_values(q, refq, 3.0e-6_dp)
    s = cspline2d(x, lamb=5.0_dp, status=status)
    if (status /= spline_status_success) error stop 'cspline2d smooth: status'
    call check_values(s, refsmooth, 3.0e-6_dp)
    f = spline_filter(x, lmbda=5.0_dp, status=status)
    if (status /= spline_status_success) error stop 'spline_filter: status'
    call check_values(f, reffilter, 3.0e-5_dp)

    ! Linear complex extension, including the smoothed case that SciPy does
    ! not currently support for complex double-precision images.
    z = cspline2d(cmplx(x, 0.2_dp*x, dp), lamb=5.0_dp, status=status)
    if (status /= spline_status_success) error stop 'complex smoothing: status'
    if (maxval(abs(z-cmplx(s,0.2_dp*s,dp))) > 2.0e-12_dp) &
        error stop 'complex smoothing: linearity'
    zf = spline_filter(cmplx(x,0.2_dp*x,dp), lmbda=5.0_dp, status=status)
    if (status /= spline_status_success) error stop 'complex spline_filter: status'
    if (maxval(abs(zf-cmplx(f,0.2_dp*f,dp))) > 3.0e-12_dp) &
        error stop 'complex spline_filter: linearity'

    ! Noncontiguous row/column sections have the same result as a copy.
    a = cspline2d(x(1:m:2,2:n:2))
    b = x(1:m:2,2:n:2)
    b = cspline2d(b)
    if (maxval(abs(a-b)) > 1.0e-12_dp) error stop 'strided image'

    ! Invalid parameters must not return a shape-conforming false success.
    a = qspline2d(x, lamb=0.1_dp, status=status)
    if (status /= spline_status_invalid_input .or. size(a) /= 0) &
        error stop 'invalid quadratic smoothing'
    a = cspline2d(x, lamb=-1.0_dp, status=status)
    if (status /= spline_status_invalid_input .or. size(a) /= 0) &
        error stop 'invalid negative smoothing'
    a = spline_filter(x, lmbda=-1.0_dp, status=status)
    if (status /= spline_status_invalid_input .or. size(a) /= 0) &
        error stop 'invalid filter smoothing'
    a = cspline2d(x(1:0,:), status=status)
    if (status /= spline_status_invalid_input .or. size(a) /= 0) &
        error stop 'invalid empty dimension'

    a = cspline2d(reshape([2.0_dp],[1,1]), status=status)
    if (status /= spline_status_success .or. size(a) /= 1) &
        error stop 'single pixel image'
    if (abs(a(1,1)-2.0_dp) > 1.0e-14_dp) error stop 'single pixel preservation'

    print '(a)', 'PASS: spline2d coefficient, smoothing, and spatial regression'

contains
    subroutine check_values(values, expected, tol)
        real(dp), intent(in) :: values(:,:) !! computed spline output
        real(dp), intent(in) :: expected(:) !! reference samples
        real(dp), intent(in) :: tol !! absolute tolerance
        integer :: k
        do k = 1, size(expected)
            if (abs(values(indices(1,k),indices(2,k)) - expected(k)) > tol) then
                print *, 'Mismatch', k, values(indices(1,k),indices(2,k)), expected(k)
                error stop 'SciPy reference mismatch'
            end if
        end do
    end subroutine check_values
end program test_signal_spline2d
