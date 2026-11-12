! SPDX-License-Identifier: MIT
! Native regression checks for weighted multi-harmonic Lomb-Scargle fitting.
program test_signal_multiharmonic
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_signal_irregular, only : lombscargle, lombscargle_multiharmonic, &
        irregular_status_success, irregular_status_invalid_input, irregular_status_singular
    implicit none
    integer, parameter :: n = 57
    real(dp) :: x(n), y(n), weights(n), xs(2*n), ys(2*n)
    real(dp), parameter :: f(5) = [0.27_dp, 0.84_dp, 1.37_dp, 1.93_dp, 3.42_dp]
    real(dp), allocatable :: p(:), p1(:), q(:)
    integer :: i, st

    do i = 1, n
        x(i) = 0.31_dp*real(i, dp) + 0.017_dp*sin(real(i*i,dp))
        weights(i) = 0.2_dp + 0.8_dp*(1.0_dp+sin(real(i,dp)))
    end do
    y = 2.0_dp + 0.8_dp*cos(1.37_dp*x) - 0.45_dp*sin(1.37_dp*x) &
        + 0.35_dp*cos(2.74_dp*x) + 0.2_dp*sin(2.74_dp*x) &
        - 0.18_dp*cos(4.11_dp*x)

    p = lombscargle_multiharmonic(x, y, f, 3, weights=weights, &
        floating_mean=.true., normalize='normalize', status=st)
    if (st /= irregular_status_success) error stop 'multi-harmonic status'
    if (abs(p(3)-1.0_dp) > 2.0e-12_dp) error stop 'exact multi-harmonic reconstruction'
    if (any(p < 0.0_dp) .or. any(p > 1.0_dp)) error stop 'normalized bounds'
    if (p(3) <= maxval(p([1,2,4,5]))) error stop 'wrong multiharmonic peak'

    ! Frequency-specific weighted offsets should absorb arbitrary DC shifts.
    q = lombscargle_multiharmonic(x, y+700.0_dp, f, 3, weights=weights, &
        floating_mean=.true., normalize='normalize', status=st)
    call check(q, p, 2.0e-12_dp)
    q = lombscargle_multiharmonic(x, y, f, 3, weights=10.0_dp*weights, &
        floating_mean=.true., normalize='normalize', status=st)
    call check(q, p, 1.0e-12_dp)

    ! The one-harmonic specialization agrees with the legacy SciPy model.
    p = lombscargle(x, y, f, nterms=3, weights=weights, &
        floating_mean=.true., normalize='normalize', status=st)
    q = lombscargle_multiharmonic(x, y, f, 3, weights=weights, &
        floating_mean=.true., normalize='normalize', status=st)
    call check(p, q, 1.0e-12_dp)

    p = lombscargle_multiharmonic(x, y, f, 1, weights=weights, &
        floating_mean=.true., status=st)
    if (st /= irregular_status_success) error stop 'one harmonic status'
    p1 = lombscargle(x, y, f, weights=weights, floating_mean=.true.)
    call check(p, p1, 3.0e-12_dp)
    p = lombscargle_multiharmonic(x, y, f, 1, weights=weights, precenter=.true., status=st)
    p1 = lombscargle(x, y, f, weights=weights, precenter=.true.)
    call check(p, p1, 3.0e-12_dp)

    ! Fortran array sections may be strided on the sample dimension.
    xs = -99.0_dp
    ys = -99.0_dp
    xs(1:2*n:2) = x
    ys(1:2*n:2) = y
    p = lombscargle_multiharmonic(xs(1:2*n:2), ys(1:2*n:2), f, 2, status=st)
    q = lombscargle_multiharmonic(x, y, f, 2, status=st)
    call check(p, q, 1.0e-12_dp)

    p = lombscargle_multiharmonic(x, y, f, 0, status=st)
    call reject(p, st)
    p = lombscargle_multiharmonic(x, y(:10), f, 2, status=st)
    call reject(p, st)
    p = lombscargle_multiharmonic(x, y, f, 2, weights=-weights, status=st)
    call reject(p, st)
    p = lombscargle_multiharmonic(x, y, f, 2, weights=0.0_dp*weights, status=st)
    call reject(p, st)
    p = lombscargle_multiharmonic(x, y, f, 2, normalize='invalid', status=st)
    call reject(p, st)
    p = lombscargle_multiharmonic(x, y, f, n, status=st)
    call reject(p, st)
    p = lombscargle_multiharmonic(x, 3.0_dp+0.0_dp*y, f, 2, &
        floating_mean=.true., normalize='normalize', status=st)
    if (st /= irregular_status_singular) error stop 'constant data status'
    if (.not. all(ieee_is_nan(p))) error stop 'constant data NaNs'
    p = lombscargle_multiharmonic(x, y, [0.0_dp, 1.37_dp], 2, &
        floating_mean=.true., normalize='normalize', status=st)
    if (st /= irregular_status_singular) error stop 'rank-deficient frequency status'
    if (.not. ieee_is_nan(p(1))) error stop 'rank-deficient frequency NaN'
    if (.not. (p(2) > 0.0_dp .and. p(2) <= 1.0_dp)) error stop 'valid frequency preserved'
    print '(a)', 'test_signal_multiharmonic passed'
contains
    subroutine check(a, b, tolerance)
        real(dp), intent(in) :: a(:) !! observed spectral powers
        real(dp), intent(in) :: b(:) !! reference spectral powers
        real(dp), intent(in) :: tolerance !! maximum absolute error
        if (size(a) /= size(b)) error stop 'spectral shape mismatch'
        if (.not. all(abs(a-b) <= tolerance)) error stop 'spectral values mismatch'
    end subroutine check

    subroutine reject(a, code)
        real(dp), intent(in) :: a(:) !! invalid-result spectral array
        integer, intent(in) :: code !! expected invalid-input status
        if (code /= irregular_status_invalid_input) error stop 'invalid input status'
        if (.not. all(ieee_is_nan(a))) error stop 'invalid input quiet NaNs'
    end subroutine reject
end program test_signal_multiharmonic
