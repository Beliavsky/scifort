! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_detrend
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_signal
    implicit none

    integer :: ierr, i
    real(dp), allocatable :: y(:)
    real(dp) :: x(8)

    y = detrend([1.0_dp, 3.0_dp, 4.0_dp, 8.0_dp], status=ierr)
    if (ierr /= detrend_status_success) error stop 'default detrend status'
    call assert_close(y, [0.3_dp, 0.1_dp, -1.1_dp, 0.7_dp])
    y = detrend([1.0_dp, 3.0_dp, 4.0_dp, 8.0_dp], kind='constant')
    call assert_close(y, [-3.0_dp, -1.0_dp, 0.0_dp, 4.0_dp])
    x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 4.0_dp, 5.0_dp, 7.0_dp, 9.0_dp]
    y = detrend(x, breaks=[5])
    call assert_close(y, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
        0.3_dp, -0.4_dp, -0.1_dp, 0.2_dp])
    y = detrend(x, breaks=[5, 5, 1, 9, 5])
    call assert_close(y, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
        0.3_dp, -0.4_dp, -0.1_dp, 0.2_dp])
    y = detrend([real(dp) :: ], status=ierr)
    if (ierr /= detrend_status_success) error stop 'empty input'
    if (size(y) /= 0) error stop 'empty output'
    y = detrend([5.0_dp])
    call assert_close(y, [0.0_dp])
    y = detrend([5.0_dp, 5.0_dp, 5.0_dp])
    call assert_close(y, [0.0_dp, 0.0_dp, 0.0_dp])
    y = detrend(x, kind='invalid', status=ierr)
    if (ierr /= detrend_status_invalid_input) error stop 'invalid kind status'
    if (.not. all(ieee_is_nan(y))) error stop 'invalid kind output'
    y = detrend(x, breaks=[0], status=ierr)
    if (ierr /= detrend_status_invalid_input) error stop 'invalid breakpoint status'
    do i = 2, 100
        block
            real(dp), allocatable :: a(:), b(:)
            integer :: j
            allocate(a(i))
            do j = 1, i
                a(j) = 7.0_dp + 2.5_dp * real(j, dp)
            end do
            b = detrend(a)
            if (maxval(abs(b)) > 1.0e-12_dp) error stop 'affine not removed'
        end block
    end do
    print '(a)', 'test_signal_detrend passed'
contains
    subroutine assert_close(a, b)
        real(dp), intent(in) :: a(:), b(:)
        if (size(a) /= size(b)) error stop 'detrend size mismatch'
        if (any(abs(a - b) > 1.0e-12_dp)) then
            print *, a
            print *, b
            error stop 'detrend mismatch'
        end if
    end subroutine
end program test_signal_detrend
