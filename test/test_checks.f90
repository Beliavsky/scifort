! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module test_checks
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    implicit none
    private

    public :: check_close
    public :: check_true

contains

    subroutine check_close(name, actual, expected, atol, rtol, failures)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual
        real(dp), intent(in) :: expected
        real(dp), intent(in) :: atol
        real(dp), intent(in) :: rtol
        integer, intent(inout) :: failures

        real(dp) :: tolerance

        if (actual == expected) return
        tolerance = atol + rtol * abs(expected)
        if (ieee_is_nan(actual) .or. .not. (abs(actual - expected) <= tolerance)) then
            failures = failures + 1
            print '(a,1x,a)', 'FAIL:', name
            print '(a,1x,es26.17e3)', '  actual   =', actual
            print '(a,1x,es26.17e3)', '  expected =', expected
            print '(a,1x,es26.17e3)', '  tolerance=', tolerance
        end if
    end subroutine check_close

    subroutine check_true(name, condition, failures)
        character(len=*), intent(in) :: name
        logical, intent(in) :: condition
        integer, intent(inout) :: failures

        if (.not. condition) then
            failures = failures + 1
            print '(a,1x,a)', 'FAIL:', name
        end if
    end subroutine check_true

end module test_checks
