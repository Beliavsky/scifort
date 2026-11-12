! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Expected error termination is checked in child processes so the parent test succeeds.

program test_sample_distance_nd_errors
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_stats, only : wasserstein_distance_nd
    implicit none

    character(len=1024) :: executable
    character(len=64) :: mode
    character(len=64) :: error_case
    character(len=32), parameter :: cases(*) = [character(len=32) :: &
        'empty_distribution', 'dimension_mismatch', 'nan_value', 'negative_weight']
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

    print *, 'multidimensional sample-distance error tests passed'

contains

    subroutine run_error_case(name)
        character(len=*), intent(in) :: name
        real(dp), allocatable :: empty(:, :)
        real(dp) :: a(2, 2)
        real(dp) :: b(2, 3)
        real(dp) :: value

        a(1, :) = [0.0_dp, 1.0_dp]
        a(2, :) = [2.0_dp, 3.0_dp]
        b(1, :) = [0.0_dp, 1.0_dp, 2.0_dp]
        b(2, :) = [3.0_dp, 4.0_dp, 5.0_dp]

        select case (name)
        case ('empty_distribution')
            allocate(empty(0, 2))
            value = wasserstein_distance_nd(empty, a)
        case ('dimension_mismatch')
            value = wasserstein_distance_nd(a, b)
        case ('nan_value')
            a(2, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
            value = wasserstein_distance_nd(a, a)
        case ('negative_weight')
            value = wasserstein_distance_nd(a, a, [1.0_dp, -1.0_dp], [1.0_dp, 1.0_dp])
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
            print *, 'FAIL expected multidimensional sample-distance error: ', trim(name)
            print *, ' cmdstat = ', cmdstat
            print *, ' exitstat= ', exitstat
            error stop 1
        end if
    end subroutine expect_error_stop

end program test_sample_distance_nd_errors
