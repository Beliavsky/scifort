! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Expected error termination is checked in child processes so the parent test succeeds.

program test_sample_distance_errors
    use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_stats, only : energy_distance, wasserstein_distance
    implicit none

    character(len=1024) :: executable
    character(len=64) :: mode
    character(len=64) :: error_case
    character(len=32), parameter :: cases(*) = [character(len=32) :: &
        'empty_distribution', 'weight_size', 'negative_weight', 'zero_weight_sum', &
        'infinite_weight', 'nan_weight', 'overflow_weight_sum']
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

    print *, 'sample-distance error tests passed'

contains

    subroutine run_error_case(name)
        character(len=*), intent(in) :: name
        real(dp), allocatable :: empty(:)
        real(dp) :: value

        select case (name)
        case ('empty_distribution')
            allocate(empty(0))
            value = wasserstein_distance(empty, [1.0_dp])
        case ('weight_size')
            value = energy_distance([1.0_dp, 2.0_dp], [3.0_dp], [1.0_dp], [1.0_dp])
        case ('negative_weight')
            value = wasserstein_distance([1.0_dp, 2.0_dp], [3.0_dp], &
                [1.0_dp, -1.0_dp], [1.0_dp])
        case ('zero_weight_sum')
            value = wasserstein_distance([1.0_dp, 2.0_dp], [3.0_dp], &
                [0.0_dp, 0.0_dp], [1.0_dp])
        case ('infinite_weight')
            value = wasserstein_distance([1.0_dp, 2.0_dp], [3.0_dp], &
                [1.0_dp, ieee_value(0.0_dp, ieee_positive_inf)], [1.0_dp])
        case ('nan_weight')
            value = wasserstein_distance([1.0_dp, 2.0_dp], [3.0_dp], &
                [1.0_dp, ieee_value(0.0_dp, ieee_quiet_nan)], [1.0_dp])
        case ('overflow_weight_sum')
            value = wasserstein_distance([1.0_dp, 2.0_dp], [3.0_dp], &
                [huge(0.0_dp), huge(0.0_dp)], [1.0_dp])
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
            print *, 'FAIL expected sample-distance error: ', trim(name)
            print *, ' cmdstat = ', cmdstat
            print *, ' exitstat= ', exitstat
            error stop 1
        end if
    end subroutine expect_error_stop

end program test_sample_distance_errors
