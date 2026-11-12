! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! One-based reference cases against scipy.signal in SciPy 1.17.0.
program test_signal_peaks
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_signal
    implicit none

    integer :: i, ierr
    integer, allocatable :: lb(:), rb(:)
    real(dp), allocatable :: p(:), w(:), wh(:), li(:), ri(:)
    type(peak_result) :: res
    real(dp) :: x(9), nanval

    x = [0.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, &
        0.0_dp, 2.0_dp, 0.0_dp]
    res = find_peaks(x)
    call assert_status(res%status)
    call assert_int(res%peaks, [2, 5, 8])
    if (size(res%prominences) /= 0) error stop 'unrequested property returned'

    res = find_peaks(x, height_min=-100.0_dp, threshold_min=-100.0_dp, &
        prominence_min=0.0_dp, width_min=0.0_dp, plateau_size_min=1)
    call assert_status(res%status)
    call assert_int(res%peaks, [2, 5, 8])
    call assert_int(res%plateau_sizes, [1, 3, 1])
    call assert_int(res%left_edges, [2, 4, 8])
    call assert_int(res%right_edges, [2, 6, 8])
    call assert_real(res%peak_heights, [1.0_dp, 2.0_dp, 2.0_dp])
    call assert_real(res%left_thresholds, [1.0_dp, 0.0_dp, 2.0_dp])
    call assert_real(res%right_thresholds, [1.0_dp, 0.0_dp, 2.0_dp])
    call assert_real(res%prominences, [1.0_dp, 2.0_dp, 2.0_dp])
    call assert_int(res%left_bases, [1, 3, 7])
    call assert_int(res%right_bases, [3, 7, 9])
    call assert_real(res%widths, [1.0_dp, 3.0_dp, 1.0_dp])
    call assert_real(res%left_ips, [1.5_dp, 3.5_dp, 7.5_dp])
    call assert_real(res%right_ips, [2.5_dp, 6.5_dp, 8.5_dp])

    res = find_peaks(x, threshold_min=0.1_dp)
    call assert_int(res%peaks, [2, 8])
    res = find_peaks(x, plateau_size_min=2)
    call assert_int(res%peaks, [5])
    res = find_peaks(x, plateau_size_max=1)
    call assert_int(res%peaks, [2, 8])
    res = find_peaks(x, prominence_min=1.5_dp)
    call assert_int(res%peaks, [5, 8])
    res = find_peaks(x, width_min=2.0_dp)
    call assert_int(res%peaks, [5])
    res = find_peaks(x, height_max=1.0_dp)
    call assert_int(res%peaks, [2])

    x(1:7) = [0.0_dp, 4.0_dp, 1.0_dp, 3.0_dp, 1.0_dp, 5.0_dp, 0.0_dp]
    res = find_peaks(x(1:7), distance=3)
    call assert_int(res%peaks, [2, 6])
    res = find_peaks(x(1:7), prominence_min=1.0_dp, width_min=0.0_dp)
    call assert_int(res%peaks, [2, 4, 6])
    call assert_real(res%prominences, [3.0_dp, 2.0_dp, 5.0_dp])
    call assert_real(res%widths, [0.875_dp, 1.0_dp, 1.125_dp])
    call assert_real(res%left_ips, [1.625_dp, 3.5_dp, 5.375_dp])
    call assert_real(res%right_ips, [2.5_dp, 4.5_dp, 6.5_dp])
    call peak_prominences(x(1:7), [2, 4, 6], p, lb, rb, ierr, wlen=3)
    call assert_status(ierr)
    call assert_real(p, [3.0_dp, 2.0_dp, 4.0_dp])
    call assert_int(lb, [1, 3, 5])
    call assert_int(rb, [3, 5, 7])

    call peak_widths(x(1:7), [2, 4, 6], w, wh, li, ri, ierr, &
        rel_height=1.0_dp, prominences=p, left_bases=lb, right_bases=rb)
    call assert_status(ierr)
    call assert_real(wh, [1.0_dp, 1.0_dp, 1.0_dp])
    call assert_real(w, [1.75_dp, 2.0_dp, 1.8_dp])

    res = find_peaks([real(dp) :: ])
    call assert_status(res%status)
    if (size(res%peaks) /= 0) error stop 'empty signal'
    res = find_peaks([1.0_dp, 2.0_dp])
    if (size(res%peaks) /= 0) error stop 'two sample signal'
    res = find_peaks(x, distance=0)
    if (res%status /= peak_status_invalid_input) error stop 'invalid distance'
    res = find_peaks(x, rel_height=-1.0_dp)
    if (res%status /= peak_status_invalid_input) error stop 'invalid relative height'
    res = find_peaks(x, wlen=1)
    if (res%status /= peak_status_invalid_input) error stop 'invalid window'
    nanval = ieee_value(0.0_dp, ieee_quiet_nan)
    res = find_peaks([0.0_dp, nanval, 1.0_dp])
    if (res%status /= peak_status_invalid_input) error stop 'NaN not rejected'
    call peak_prominences(x, [0], p, lb, rb, ierr)
    if (ierr /= peak_status_invalid_input) error stop 'invalid peak index'
    call peak_widths(x, [2], w, wh, li, ri, ierr, rel_height=-0.5_dp)
    if (ierr /= peak_status_invalid_input) error stop 'invalid width height'

    do i = 1, 12
        x = 0.0_dp
        x(2:8:2) = real(i, dp)
        res = find_peaks(x, prominence_min=0.1_dp, width_min=0.0_dp)
        call assert_status(res%status)
        if (size(res%peaks) /= 4) error stop 'peak count'
    end do
    print '(a)', 'test_signal_peaks passed'

contains
    subroutine assert_status(status)
        integer, intent(in) :: status
        if (status /= peak_status_success) error stop 'peak status failure'
    end subroutine

    subroutine assert_int(a, b)
        integer, intent(in) :: a(:), b(:)
        if (size(a) /= size(b)) error stop 'peak integer size mismatch'
        if (any(a /= b)) then
            print *, 'got: ', a
            print *, 'want:', b
            error stop 'peak integer mismatch'
        end if
    end subroutine

    subroutine assert_real(a, b)
        real(dp), intent(in) :: a(:), b(:)
        if (size(a) /= size(b)) error stop 'peak real size mismatch'
        if (any(abs(a - b) > 1.0e-12_dp)) then
            print *, 'got: ', a
            print *, 'want:', b
            error stop 'peak real mismatch'
        end if
    end subroutine
end program test_signal_peaks
