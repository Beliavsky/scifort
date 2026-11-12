! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Expected error termination is checked in child processes so the parent FPM
! test succeeds.

program test_descriptive_lmoment_errors
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_stats, only : lmoment
    implicit none

    character(len=1024) :: executable
    character(len=64) :: mode
    character(len=64) :: error_case
    character(len=24), parameter :: cases(*) = [character(len=24) :: &
        'scalar_order_zero', 'vector_order_zero', 'empty_order_vector', &
        'default_nan_raise', 'scalar_nan_raise']
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

    print *, 'descriptive lmoment error tests passed'

contains

    subroutine run_error_case(name)
        character(len=*), intent(in) :: name !! child test case

        real(dp) :: qnan
        real(dp) :: scalar
        real(dp), allocatable :: values(:)

        qnan = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(values(1))
        values = 0.0_dp
        select case (name)
        case ('scalar_order_zero')
            scalar = lmoment([1.0_dp, 2.0_dp, 3.0_dp], 0)
        case ('vector_order_zero')
            values = lmoment([1.0_dp, 2.0_dp, 3.0_dp], [1, 0, 2])
        case ('empty_order_vector')
            values = lmoment([1.0_dp, 2.0_dp, 3.0_dp], [integer ::])
        case ('default_nan_raise')
            values = lmoment([1.0_dp, qnan, 3.0_dp], nan_policy='raise')
        case ('scalar_nan_raise')
            scalar = lmoment([1.0_dp, qnan, 3.0_dp], 2, nan_policy='raise')
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
            print *, 'FAIL expected lmoment error: ', trim(name)
            print *, ' cmdstat = ', cmdstat
            print *, ' exitstat= ', exitstat
            error stop 1
        end if
    end subroutine expect_error_stop

end program test_descriptive_lmoment_errors
