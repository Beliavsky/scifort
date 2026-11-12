! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Expected error termination is checked in child processes so the parent test succeeds.

program test_probability_plot_errors
    use scifort_kinds, only : dp
    use scifort_stats, only : boxcox_normplot, ppcc_max, ppcc_plot, probplot
    implicit none

    character(len=1024) :: executable
    character(len=64) :: mode
    character(len=64) :: error_case
    character(len=32), parameter :: cases(*) = [character(len=32) :: &
        'ppcc_bad_interval', 'ppcc_bad_dist', 'ppcc_empty_max', &
        'boxcox_nonpositive', 'probplot_bad_params']
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

    print *, 'probability-plot error tests passed'

contains

    subroutine run_error_case(name)
        character(len=*), intent(in) :: name
        real(dp), allocatable :: a(:)
        real(dp), allocatable :: b(:)
        real(dp), allocatable :: empty(:)
        real(dp) :: value

        select case (name)
        case ('ppcc_bad_interval')
            call ppcc_plot([1.0_dp, 2.0_dp], 1.0_dp, -1.0_dp, a, b)
        case ('ppcc_bad_dist')
            value = ppcc_max([1.0_dp, 2.0_dp, 3.0_dp], dist='norm')
        case ('ppcc_empty_max')
            allocate(empty(0))
            value = ppcc_max(empty)
        case ('boxcox_nonpositive')
            call boxcox_normplot([0.0_dp, 1.0_dp, 2.0_dp], -1.0_dp, 1.0_dp, a, b)
        case ('probplot_bad_params')
            call probplot([1.0_dp, 2.0_dp], a, b, dist='t', sparams=[real(dp) ::])
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
            print *, 'FAIL expected probability-plot error: ', trim(name)
            print *, ' cmdstat = ', cmdstat
            print *, ' exitstat= ', exitstat
            error stop 1
        end if
    end subroutine expect_error_stop

end program test_probability_plot_errors
