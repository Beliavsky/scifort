! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Error-path regression for information-theory descriptive statistics.

program test_descriptive_entropy_errors
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_stats, only : differential_entropy, entropy
    implicit none

    character(len=1024) :: executable
    character(len=64) :: mode
    character(len=64) :: error_case
    character(len=40), parameter :: cases(*) = [character(len=40) :: &
        'entropy_nan_raise', 'entropy_q_nan_raise', 'entropy_bad_policy', &
        'entropy_bad_base', 'entropy_bad_broadcast', 'diff_nan_raise', &
        'diff_bad_policy', 'diff_bad_base', 'diff_bad_method']
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

    print *, 'descriptive entropy error tests passed'

contains

    subroutine run_error_case(name)
        character(len=*), intent(in) :: name
        real(dp) :: qnan
        real(dp) :: scalar
        real(dp), allocatable :: values(:)
        real(dp) :: matrix(2, 3)
        real(dp) :: bad_q(2, 2)

        qnan = ieee_value(0.0_dp, ieee_quiet_nan)
        matrix = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp], shape(matrix))
        bad_q = 1.0_dp
        select case (name)
        case ('entropy_nan_raise')
            scalar = entropy([1.0_dp, qnan], nan_policy='raise')
        case ('entropy_q_nan_raise')
            scalar = entropy([1.0_dp, 2.0_dp], qk=[1.0_dp, qnan], nan_policy='raise')
        case ('entropy_bad_policy')
            scalar = entropy([1.0_dp, 2.0_dp], nan_policy='bad')
        case ('entropy_bad_base')
            scalar = entropy([1.0_dp, 2.0_dp], base=0.0_dp)
        case ('entropy_bad_broadcast')
            values = entropy(matrix, qk=bad_q, axis=2)
        case ('diff_nan_raise')
            scalar = differential_entropy([1.0_dp, qnan, 2.0_dp, 3.0_dp, 4.0_dp], &
                nan_policy='raise')
        case ('diff_bad_policy')
            scalar = differential_entropy([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], &
                nan_policy='bad')
        case ('diff_bad_base')
            scalar = differential_entropy([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], &
                base=-2.0_dp)
        case ('diff_bad_method')
            scalar = differential_entropy([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], &
                method='bad')
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
            print *, 'FAIL expected entropy error: ', trim(name)
            error stop 1
        end if
    end subroutine expect_error_stop

end program test_descriptive_entropy_errors
