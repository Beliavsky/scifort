! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Expected error termination is checked in child processes so the parent FPM
! test succeeds.

program test_power_transform_errors
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_stats, only : boxcox_fit, boxcox_llf, boxcox_normmax, &
        yeojohnson_llf, yeojohnson_normmax
    implicit none

    character(len=1024) :: executable
    character(len=64) :: mode
    character(len=64) :: error_case
    character(len=32), parameter :: cases(*) = [character(len=32) :: &
        'boxcox_nan_raise', 'yeojohnson_nan_raise', 'boxcox_bad_method', &
        'boxcox_nonpositive_fit', 'yeojohnson_nonfinite_normmax']
    integer :: argc
    integer :: i

    argc = command_argument_count()
    if (argc >= 1) then
        call get_command_argument(1, mode)
        if (trim(mode) == '--error-case') then
            if (argc < 2) error stop 2
            call get_command_argument(2, error_case)
            call run_error_case(trim(error_case))
            stop 0
        end if
    end if

    call get_command_argument(0, executable)
    do i = 1, size(cases)
        call expect_error_stop(trim(executable), trim(cases(i)))
    end do

    print *, 'power-transform error tests passed'

contains

    subroutine run_error_case(name)
        character(len=*), intent(in) :: name
        real(dp) :: qnan
        real(dp) :: scalar
        real(dp) :: lmbda
        real(dp), allocatable :: values(:)

        qnan = ieee_value(0.0_dp, ieee_quiet_nan)
        select case (name)
        case ('boxcox_nan_raise')
            scalar = boxcox_llf(0.5_dp, [1.0_dp, qnan, 3.0_dp], nan_policy='raise')
        case ('yeojohnson_nan_raise')
            scalar = yeojohnson_llf(0.5_dp, [-1.0_dp, qnan, 3.0_dp], nan_policy='raise')
        case ('boxcox_bad_method')
            scalar = boxcox_normmax([1.0_dp, 2.0_dp, 3.0_dp], method='unsupported')
        case ('boxcox_nonpositive_fit')
            call boxcox_fit([0.0_dp, 1.0_dp, 2.0_dp], values, lmbda)
        case ('yeojohnson_nonfinite_normmax')
            scalar = yeojohnson_normmax([1.0_dp, qnan, 2.0_dp])
        case default
            error stop 4
        end select
    end subroutine run_error_case

    subroutine expect_error_stop(program_path, name)
        character(len=*), intent(in) :: program_path
        character(len=*), intent(in) :: name
        character(len=2200) :: command
        integer :: cmdstat
        integer :: envstat
        integer :: exitstat
        character(len=64) :: os_name

        call get_environment_variable('OS', os_name, status=envstat)
        command = '"' // trim(program_path) // '" --error-case ' // trim(name)
        if (envstat == 0 .and. index(os_name, 'Windows_NT') > 0) then
            command = trim(command) // ' >NUL 2>&1'
        else
            command = trim(command) // ' >/dev/null 2>&1'
        end if
        call execute_command_line(trim(command), wait=.true., exitstat=exitstat, &
            cmdstat=cmdstat)
        if (cmdstat /= 0 .or. exitstat == 0) then
            print *, 'FAIL expected power-transform error: ', trim(name)
            print *, ' cmdstat = ', cmdstat
            print *, ' exitstat= ', exitstat
            error stop 1
        end if
    end subroutine expect_error_stop

end program test_power_transform_errors
