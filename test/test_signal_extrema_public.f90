! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_extrema_public
    use scifort_kinds, only : dp
    use scifort_signal, only : argrelmax, argrelmin, argrelextrema, peak_detect, &
        extrema_status_success, excursion_status_success
    implicit none

    integer, allocatable :: idx(:), imx(:), imn(:)
    real(dp), allocatable :: mx(:), mn(:)
    integer :: status

    idx = argrelmax([0.0_dp, 2.0_dp, 0.0_dp], status=status)
    if (status /= extrema_status_success) error stop 'argrelmax failed from public signal module'
    if (size(idx) /= 1) error stop 'argrelmax count failure'
    if (idx(1) /= 2) error stop 'argrelmax location failure'
    idx = argrelmin([0.0_dp, -2.0_dp, 0.0_dp])
    if (size(idx) /= 1) error stop 'argrelmin count failure'
    if (idx(1) /= 2) error stop 'argrelmin location failure'
    idx = argrelextrema([1.0_dp, 1.0_dp], greater_or_equal)
    if (size(idx) /= 2) error stop 'argrelextrema custom comparator failure'
    call peak_detect([0.0_dp, 2.0_dp, 0.0_dp], 1.0_dp, imx, mx, imn, mn, status)
    if (status /= excursion_status_success) error stop 'peak_detect failed from public signal module'
    if (size(imx) /= 1) error stop 'peak_detect count failure'
    if (imx(1) /= 2) error stop 'peak_detect location failure'
    print '(a)', 'test_signal_extrema_public passed'

contains

    pure logical function greater_or_equal(a, b)
        real(dp), intent(in) :: a !! candidate
        real(dp), intent(in) :: b !! neighbor
        greater_or_equal = a >= b
    end function greater_or_equal
end program test_signal_extrema_public
