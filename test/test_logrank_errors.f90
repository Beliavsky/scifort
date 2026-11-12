! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_logrank_errors
    use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_stats, only : censored_data, logrank, logrank_result
    implicit none

    character(len=1024) :: executable
    character(len=64) :: mode
    character(len=64) :: error_case
    character(len=32), parameter :: cases(*) = [character(len=32) :: &
        'bad_alternative', 'nan_array', 'infinite_event', 'nan_censor']
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

    print *, 'logrank error tests passed'

contains

    subroutine run_error_case(name)
        character(len=*), intent(in) :: name
        type(logrank_result) :: result
        type(censored_data) :: sample

        select case (name)
        case ('bad_alternative')
            result = logrank([1.0_dp], [2.0_dp], alternative='bad')
        case ('nan_array')
            result = logrank([1.0_dp, ieee_value(0.0_dp, ieee_quiet_nan)], [2.0_dp])
        case ('infinite_event')
            sample = censored_data(uncensored=[ieee_value(0.0_dp, ieee_positive_inf)], &
                right=[real(dp) ::])
            result = logrank(sample, [2.0_dp])
        case ('nan_censor')
            sample = censored_data(uncensored=[1.0_dp], &
                right=[ieee_value(0.0_dp, ieee_quiet_nan)])
            result = logrank(sample, [2.0_dp])
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
        call execute_command_line(trim(command), wait=.true., exitstat=exitstat, cmdstat=cmdstat)
        if (cmdstat /= 0 .or. exitstat == 0) then
            print *, 'FAIL expected logrank error: ', trim(name)
            error stop 1
        end if
    end subroutine expect_error_stop

end program test_logrank_errors
