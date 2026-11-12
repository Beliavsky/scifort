! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Expected error termination is checked in child processes so the parent test succeeds.

program test_empirical_distribution_errors
    use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_stats, only : confidence_interval, ecdf, ecdf_confidence_interval, ecdf_result, &
        normal_cdf, normal_pdf, normal_ppf, order_statistic, &
        order_statistic_distribution, quantile_test, quantile_test_result
    implicit none

    character(len=1024) :: executable
    character(len=64) :: mode
    character(len=64) :: error_case
    character(len=32), parameter :: cases(*) = [character(len=32) :: &
        'ecdf_nan', 'ecdf_infinite', 'ecdf_bad_method', 'ecdf_nested_ci', &
        'quantile_bad_p', 'quantile_bad_alt', 'quantile_bad_ci', 'quantile_nan_ci', &
        'order_negative']
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

    print *, 'empirical-distribution error tests passed'

contains

    subroutine run_error_case(name)
        character(len=*), intent(in) :: name
        type(ecdf_confidence_interval) :: ci
        type(ecdf_result) :: er
        type(order_statistic_distribution) :: os
        type(quantile_test_result) :: qr
        type(confidence_interval) :: qci

        select case (name)
        case ('ecdf_nan')
            er = ecdf([1.0_dp, ieee_value(0.0_dp, ieee_quiet_nan)])
        case ('ecdf_infinite')
            er = ecdf([1.0_dp, ieee_value(0.0_dp, ieee_positive_inf)])
        case ('ecdf_bad_method')
            er = ecdf([1.0_dp, 2.0_dp])
            ci = er%cdf%confidence_interval(method='bad')
        case ('ecdf_nested_ci')
            er = ecdf([1.0_dp, 2.0_dp])
            ci = er%cdf%confidence_interval()
            ci = ci%low%confidence_interval()
        case ('quantile_bad_p')
            qr = quantile_test([1.0_dp, 2.0_dp], p=1.0_dp)
        case ('quantile_bad_alt')
            qr = quantile_test([1.0_dp, 2.0_dp], alternative='bad')
        case ('quantile_bad_ci')
            qr = quantile_test([1.0_dp, 2.0_dp])
            qci = qr%confidence_interval(1.0_dp)
        case ('quantile_nan_ci')
            qr = quantile_test([1.0_dp, 2.0_dp], p=ieee_value(0.0_dp, ieee_quiet_nan))
            qci = qr%confidence_interval()
        case ('order_negative')
            os = order_statistic(base_cdf, base_pdf, base_ppf, -1, 5)
        case default
            error stop 4
        end select
    end subroutine run_error_case

    pure function base_cdf(value, parameters) result(probability)
        real(dp), intent(in) :: value
        real(dp), intent(in) :: parameters(:)
        real(dp) :: probability

        if (size(parameters) /= 0) error stop 'unexpected callback parameters'
        probability = normal_cdf(value)
    end function base_cdf

    pure function base_pdf(value, parameters) result(probability)
        real(dp), intent(in) :: value
        real(dp), intent(in) :: parameters(:)
        real(dp) :: probability

        if (size(parameters) /= 0) error stop 'unexpected callback parameters'
        probability = normal_pdf(value)
    end function base_pdf

    pure function base_ppf(probability, parameters) result(value)
        real(dp), intent(in) :: probability
        real(dp), intent(in) :: parameters(:)
        real(dp) :: value

        if (size(parameters) /= 0) error stop 'unexpected callback parameters'
        value = normal_ppf(probability)
    end function base_ppf

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
            print *, 'FAIL expected empirical-distribution error: ', trim(name)
            print *, ' cmdstat = ', cmdstat
            print *, ' exitstat= ', exitstat
            error stop 1
        end if
    end subroutine expect_error_stop

end program test_empirical_distribution_errors
