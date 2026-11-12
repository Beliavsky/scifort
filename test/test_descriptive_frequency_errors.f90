! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Error-path regression for frequency descriptive statistics. Expected error
! termination is checked in child processes so the parent FPM test succeeds.

program test_descriptive_frequency_errors
    use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_stats, only : cumfreq, cumfreq_result, percentileofscore, relfreq, &
        relfreq_result, scoreatpercentile
    implicit none

    character(len=1024) :: executable
    character(len=64) :: mode
    character(len=64) :: error_case
    character(len=32), parameter :: cases(*) = [character(len=32) :: &
        'percentile_nan_raise', 'percentile_bad_kind', 'percentile_bad_policy', &
        'score_bad_percentile', 'score_bad_method', 'score_bad_limit', &
        'cumfreq_bad_bins', 'cumfreq_auto_nan', 'relfreq_auto_inf', &
        'cumfreq_weight_shape']
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

    print *, 'descriptive frequency error tests passed'

contains

    subroutine run_error_case(name)
        character(len=*), intent(in) :: name !! child test case

        type(cumfreq_result) :: cf
        real(dp) :: qnan
        type(relfreq_result) :: rf
        real(dp) :: scalar

        qnan = ieee_value(0.0_dp, ieee_quiet_nan)
        select case (name)
        case ('percentile_nan_raise')
            scalar = percentileofscore([1.0_dp, qnan], 2.0_dp, nan_policy='raise')
        case ('percentile_bad_kind')
            scalar = percentileofscore([1.0_dp, 2.0_dp], 2.0_dp, kind='bad')
        case ('percentile_bad_policy')
            scalar = percentileofscore([1.0_dp, 2.0_dp], 2.0_dp, nan_policy='bad')
        case ('score_bad_percentile')
            scalar = scoreatpercentile([1.0_dp, 2.0_dp], -1.0_dp)
        case ('score_bad_method')
            scalar = scoreatpercentile([1.0_dp, 2.0_dp], 50.0_dp, interpolation_method='bad')
        case ('score_bad_limit')
            scalar = scoreatpercentile([1.0_dp, 2.0_dp], 50.0_dp, limit=[2.0_dp, 1.0_dp])
        case ('cumfreq_bad_bins')
            cf = cumfreq([1.0_dp, 2.0_dp], numbins=0)
        case ('cumfreq_auto_nan')
            cf = cumfreq([1.0_dp, qnan], numbins=4)
        case ('relfreq_auto_inf')
            rf = relfreq([1.0_dp, ieee_value(0.0_dp, ieee_positive_inf)], numbins=4)
        case ('cumfreq_weight_shape')
            cf = cumfreq(reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [2, 2]), &
                weights=[1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp])
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
            print *, 'FAIL expected frequency error: ', trim(name)
            print *, ' cmdstat = ', cmdstat
            print *, ' exitstat= ', exitstat
            error stop 1
        end if
    end subroutine expect_error_stop

end program test_descriptive_frequency_errors
