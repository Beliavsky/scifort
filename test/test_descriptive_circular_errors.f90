! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Error-path regression for circular descriptive statistics. Expected error
! termination is checked in child processes so the parent FPM test succeeds.

program test_descriptive_circular_errors
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_stats, only : circmean, circstd, circvar
    implicit none

    character(len=1024) :: executable
    character(len=64) :: mode
    character(len=64) :: error_case
    character(len=24), parameter :: cases(*) = [character(len=24) :: &
        'circmean_nan_raise', 'circvar_nan_raise', 'circstd_nan_raise', &
        'circmean_zero_period', 'circstd_zero_period']
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

    print *, 'descriptive circular error tests passed'

contains

    subroutine run_error_case(name)
        character(len=*), intent(in) :: name !! child test case

        real(dp) :: qnan
        real(dp) :: scalar

        qnan = ieee_value(0.0_dp, ieee_quiet_nan)
        select case (name)
        case ('circmean_nan_raise')
            scalar = circmean([0.0_dp, qnan], nan_policy='raise')
        case ('circvar_nan_raise')
            scalar = circvar([0.0_dp, qnan], nan_policy='raise')
        case ('circstd_nan_raise')
            scalar = circstd([0.0_dp, qnan], nan_policy='raise')
        case ('circmean_zero_period')
            scalar = circmean([0.0_dp, 1.0_dp], high=1.0_dp, low=1.0_dp)
        case ('circstd_zero_period')
            scalar = circstd([0.0_dp, 1.0_dp], high=2.0_dp, low=2.0_dp)
        case default
            error stop 4
        end select
    end subroutine run_error_case

    subroutine expect_error_stop(program_path, name)
        character(len=*), intent(in) :: program_path !! current executable path
        character(len=*), intent(in) :: name !! child case expected to terminate

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
            print *, 'FAIL expected circular error: ', trim(name)
            print *, ' cmdstat = ', cmdstat
            print *, ' exitstat= ', exitstat
            error stop 1
        end if
    end subroutine expect_error_stop

end program test_descriptive_circular_errors
