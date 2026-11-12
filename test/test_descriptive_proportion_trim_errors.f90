! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Error-path regression for SciPy-style proportion trimming. Expected error
! termination is checked in child processes so the parent FPM test succeeds.

program test_descriptive_proportion_trim_errors
    use scifort_kinds, only : dp
    use scifort_stats, only : trim1, trim_mean, trimboth
    implicit none

    character(len=1024) :: executable
    character(len=64) :: mode
    character(len=64) :: error_case
    character(len=24), parameter :: cases(*) = [character(len=24) :: &
        'trim_mean_large', 'trimboth_half', 'trim1_bad_tail', 'negative_fraction']
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

    print *, 'descriptive proportion trimming error tests passed'

contains

    subroutine run_error_case(name)
        character(len=*), intent(in) :: name !! child test case

        real(dp), allocatable :: values(:)
        real(dp) :: scalar

        select case (name)
        case ('trim_mean_large')
            scalar = trim_mean([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], 0.8_dp)
        case ('trimboth_half')
            values = trimboth([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], 0.5_dp)
        case ('trim1_bad_tail')
            values = trim1([1.0_dp, 2.0_dp, 3.0_dp], 0.1_dp, tail='middle')
        case ('negative_fraction')
            scalar = trim_mean([1.0_dp, 2.0_dp, 3.0_dp], -0.1_dp)
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
            print *, 'FAIL expected trimming error: ', trim(name)
            print *, ' cmdstat = ', cmdstat
            print *, ' exitstat= ', exitstat
            error stop 1
        end if
    end subroutine expect_error_stop

end program test_descriptive_proportion_trim_errors
