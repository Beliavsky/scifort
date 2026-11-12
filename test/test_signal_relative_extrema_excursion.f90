! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_relative_extrema_excursion
    use scifort_kinds, only : dp
    use scifort_signal_relative_extrema
    use scifort_signal_peak_excursion
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none

    real(dp) :: x(9), qnan
    real(dp), allocatable :: mx(:), mn(:)
    integer, allocatable :: imx(:), imn(:), ind(:)
    integer :: status

    x = [0.0_dp, 2.0_dp, 0.0_dp, 3.0_dp, -1.0_dp, 4.0_dp, 0.0_dp, 2.0_dp, 0.0_dp]

    ind = argrelmax(x, status=status)
    call check_status(status, extrema_status_success)
    call check_indices(ind, [2, 4, 6, 8])
    ind = argrelmin(x)
    call check_indices(ind, [3, 5, 7])
    ind = argrelmax(x, order=2)
    call check_indices(ind, [6])
    ind = argrelmin(x, order=2)
    call check_indices(ind, [5])
    ind = argrelmax([1.0_dp, 0.0_dp, 0.0_dp], mode='wrap')
    call check_indices(ind, [1])
    ind = argrelmax([1.0_dp, 0.0_dp, 0.0_dp], mode='clip')
    call check_indices(ind, [integer :: ])
    ind = argrelextrema([1.0_dp, 1.0_dp, 1.0_dp], greater_or_equal)
    call check_indices(ind, [1, 2, 3])
    ind = argrelmax(x(2:8:2))
    call check_indices(ind, [3])
    ind = argrelmax(x, order=0, status=status)
    call check_status(status, extrema_status_invalid_input)
    call check_indices(ind, [integer :: ])
    ind = argrelmin(x, mode='no-such-mode', status=status)
    call check_status(status, extrema_status_invalid_input)
    call check_indices(ind, [integer :: ])
    ind = argrelmin([real(dp) :: ], mode='wrap', status=status)
    call check_status(status, extrema_status_success)
    call check_indices(ind, [integer :: ])

    call peak_detect(x, 1.0_dp, imx, mx, imn, mn, status)
    call check_status(status, excursion_status_success)
    call check_indices(imx, [2, 4, 6, 8])
    call check_indices(imn, [3, 5, 7])
    call check_values(mx, [2.0_dp, 3.0_dp, 4.0_dp, 2.0_dp])
    call check_values(mn, [0.0_dp, -1.0_dp, 0.0_dp])

    call peak_detect([1.0_dp, 2.0_dp, 3.0_dp], 0.2_dp, imx, mx, imn, mn)
    if (size(imx) /= 0 .or. size(imn) /= 0) error stop 'monotonic signal should return no peaks'
    call peak_detect([real(dp) :: ], 1.0_dp, imx, mx, imn, mn, status)
    call check_status(status, excursion_status_success)
    call peak_detect(x, 0.0_dp, imx, mx, imn, mn, status)
    call check_status(status, excursion_status_invalid_input)
    if (size(imx) /= 0) error stop 'invalid delta should return empty outputs'
    qnan = ieee_value(0.0_dp, ieee_quiet_nan)
    call peak_detect([0.0_dp, qnan, 1.0_dp], 1.0_dp, imx, mx, imn, mn, status)
    call check_status(status, excursion_status_invalid_input)
    print '(a)', 'test_signal_relative_extrema_excursion passed'

contains

    pure logical function greater_or_equal(a, b)
        real(dp), intent(in) :: a !! candidate
        real(dp), intent(in) :: b !! neighbor
        greater_or_equal = a >= b
    end function greater_or_equal

    subroutine check_indices(a, b)
        integer, intent(in) :: a(:) !! actual indices
        integer, intent(in) :: b(:) !! expected indices
        if (size(a) /= size(b)) then
            print *, 'received: ', a
            print *, 'expected: ', b
            error stop 'index count mismatch'
        end if
        if (any(a /= b)) error stop 'incorrect extrema indices'
    end subroutine check_indices

    subroutine check_values(a, b)
        real(dp), intent(in) :: a(:) !! actual values
        real(dp), intent(in) :: b(:) !! expected values
        if (size(a) /= size(b)) error stop 'value count mismatch'
        if (any(abs(a - b) > 1.0e-12_dp)) error stop 'incorrect extrema values'
    end subroutine check_values

    subroutine check_status(actual, expected)
        integer, intent(in) :: actual !! returned status
        integer, intent(in) :: expected !! expected status
        if (actual /= expected) error stop 'unexpected extrema status'
    end subroutine check_status
end program test_signal_relative_extrema_excursion
