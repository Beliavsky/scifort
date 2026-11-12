! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_lfilter
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_signal
    implicit none

    real(dp) :: b(3), a(3), x(15)
    real(dp), allocatable :: y(:), y2(:), zi(:), zf(:), intermediate(:)
    complex(dp) :: bc(2), ac(2), xc(10)
    complex(dp), allocatable :: yc(:), zic(:), zfc(:)
    integer :: status, i

    b = [0.2_dp, 0.1_dp, 0.05_dp]
    a = [2.0_dp, -0.7_dp, 0.14_dp]
    x = [(real(i, dp) / 10.0_dp, i = 1, 15)]
    y = lfilter(b, a, x, status=status)
    if (status /= lfilter_status_success) error stop 'real filter status'
    if (size(y) /= size(x)) error stop 'real filter shape'
    y2 = lfilter(b, a, x(1:6), zf=zf, status=status)
    if (status /= lfilter_status_success) error stop 'first block status'
    intermediate = lfilter(b, a, x(7:), zi=zf, zf=zi, status=status)
    if (status /= lfilter_status_success) error stop 'second block status'
    call assert_real(y, [y2, intermediate])
    if (size(zi) /= 2) error stop 'final filter state size'

    zi = lfilter_zi(b, a, status=status)
    if (status /= lfilter_status_success) error stop 'real step state status'
    y = lfilter(b, a, spread(1.0_dp, 1, 15), zi=zi)
    call assert_real(y, spread(sum(b) / sum(a), 1, 15))
    if (maxval(abs(lfilter([1.0_dp], [2.0_dp], x) - x / 2.0_dp)) > 1.0e-13_dp) &
        error stop 'zero-state scalar filter'
    zi = lfilter_zi([1.0_dp], [2.0_dp], status=status)
    if (size(zi) /= 0 .or. status /= lfilter_status_success) error stop 'scalar step state'

    bc = [cmplx(1.0_dp, -0.25_dp, dp), cmplx(0.5_dp, 0.4_dp, dp)]
    ac = [cmplx(2.0_dp, 0.0_dp, dp), cmplx(-0.4_dp, -0.3_dp, dp)]
    do i = 1, size(xc)
        xc(i) = cmplx(real(i, dp), real(i * i, dp) / 5.0_dp, dp)
    end do
    yc = lfilter(bc, ac, xc, status=status)
    if (status /= lfilter_status_success) error stop 'complex filter status'
    zic = lfilter_zi(bc, ac, status=status)
    if (status /= lfilter_status_success) error stop 'complex state status'
    if (size(zic) /= 1) error stop 'complex state length'
    yc = lfilter(bc, ac, xc(1:5), zf=zfc)
    zic = lfilter(bc, ac, xc(6:), zi=zfc, zf=zic)
    yc = [yc, zic]
    call assert_complex(yc, lfilter(bc, ac, xc))

    y = lfilter([1.0_dp], [0.0_dp], x, status=status)
    if (status /= lfilter_status_invalid_input) error stop 'zero denominator status'
    if (.not. all(ieee_is_nan(y))) error stop 'zero denominator output'
    y = lfilter(b, a, x, zi=[0.0_dp], status=status)
    if (status /= lfilter_status_invalid_input) error stop 'wrong state size'
    zi = lfilter_zi([1.0_dp], [1.0_dp, -1.0_dp], status=status)
    if (status /= lfilter_status_no_steady_state) error stop 'unit pole state'
    if (.not. all(ieee_is_nan(zi))) error stop 'unit pole output'
    print '(a)', 'test_signal_lfilter passed'
contains
    subroutine assert_real(actual, expected)
        real(dp), intent(in) :: actual(:) !! actual output
        real(dp), intent(in) :: expected(:) !! expected output
        if (size(actual) /= size(expected)) error stop 'lfilter real length'
        if (any(abs(actual - expected) > 1.0e-12_dp)) error stop 'lfilter real mismatch'
    end subroutine assert_real
    subroutine assert_complex(actual, expected)
        complex(dp), intent(in) :: actual(:) !! actual complex output
        complex(dp), intent(in) :: expected(:) !! expected complex output
        if (size(actual) /= size(expected)) error stop 'lfilter complex length'
        if (any(abs(actual - expected) > 1.0e-12_dp)) error stop 'lfilter complex mismatch'
    end subroutine assert_complex
end program test_signal_lfilter
