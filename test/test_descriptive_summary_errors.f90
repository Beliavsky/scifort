! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Error-path regression for mode and describe.

program test_descriptive_summary_errors
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_stats, only : describe, describe_all, mode, mode_all
    implicit none

    character(len=1024) :: executable
    character(len=64) :: launch_mode
    character(len=64) :: error_case
    character(len=40), parameter :: cases(*) = [character(len=40) :: &
        'mode_nan_raise', 'mode_bad_policy', 'mode_bad_axis', &
        'describe_nan_raise', 'describe_bad_policy', 'describe_bad_axis', &
        'describe_empty']
    integer :: argc
    integer :: i

    argc = command_argument_count()
    if (argc >= 1) then
        call get_command_argument(1, launch_mode)
        if (trim(launch_mode) == '--error-case') then
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

    print *, 'descriptive summary error tests passed'

contains

    subroutine run_error_case(name)
        character(len=*), intent(in) :: name
        real(dp) :: qnan
        real(dp) :: mode_scalar
        integer :: count_scalar
        real(dp), allocatable :: modes(:)
        integer, allocatable :: counts(:)
        integer, allocatable :: nobs(:)
        real(dp), allocatable :: minimum_value(:)
        real(dp), allocatable :: maximum_value(:)
        real(dp), allocatable :: mean_value(:)
        real(dp), allocatable :: variance_value(:)
        real(dp), allocatable :: skewness_value(:)
        real(dp), allocatable :: kurtosis_value(:)
        integer :: nobs_scalar
        real(dp) :: minimum_scalar
        real(dp) :: maximum_scalar
        real(dp) :: mean_scalar
        real(dp) :: variance_scalar
        real(dp) :: skewness_scalar
        real(dp) :: kurtosis_scalar
        real(dp) :: matrix(2, 2)
        real(dp) :: empty(0)

        qnan = ieee_value(0.0_dp, ieee_quiet_nan)
        matrix = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], shape(matrix))
        select case (name)
        case ('mode_nan_raise')
            call mode_all([1.0_dp, qnan], mode_scalar, count_scalar, nan_policy='raise')
        case ('mode_bad_policy')
            call mode_all([1.0_dp, 2.0_dp], mode_scalar, count_scalar, nan_policy='bad')
        case ('mode_bad_axis')
            call mode(matrix, modes, counts, axis=3)
        case ('describe_nan_raise')
            call describe_all([1.0_dp, qnan], nobs_scalar, minimum_scalar, maximum_scalar, &
                mean_scalar, variance_scalar, skewness_scalar, kurtosis_scalar, &
                nan_policy='raise')
        case ('describe_bad_policy')
            call describe_all([1.0_dp, 2.0_dp], nobs_scalar, minimum_scalar, maximum_scalar, &
                mean_scalar, variance_scalar, skewness_scalar, kurtosis_scalar, &
                nan_policy='bad')
        case ('describe_bad_axis')
            call describe(matrix, nobs, minimum_value, maximum_value, mean_value, &
                variance_value, skewness_value, kurtosis_value, axis=3)
        case ('describe_empty')
            call describe_all(empty, nobs_scalar, minimum_scalar, maximum_scalar, &
                mean_scalar, variance_scalar, skewness_scalar, kurtosis_scalar)
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
            print *, 'FAIL expected summary error: ', trim(name)
            error stop 1
        end if
    end subroutine expect_error_stop

end program test_descriptive_summary_errors
