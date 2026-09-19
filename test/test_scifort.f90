! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_scifort
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_negative_inf, ieee_positive_inf, ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_c_binding, only : c_double, c_int, c_size_t
    use scifort_c_api, only : scifort_normal_cdf_f64, &
        scifort_normal_cdf_vec_f64, scifort_version_major, &
        scifort_version_minor, scifort_version_patch
    use scifort_stats
    implicit none

    integer :: failures

    failures = 0
    call test_normal(failures)
    call test_uniform(failures)
    call test_exponential(failures)
    call test_laplace(failures)
    call test_logistic(failures)
    call test_cauchy(failures)
    call test_invalid_inputs(failures)
    call test_ieee_inputs(failures)
    call test_elemental_arrays(failures)
    call test_distribution_identities(failures)
    call test_c_api(failures)

    if (failures == 0) then
        print '(a)', 'test_scifort: PASS'
    else
        print '(a,1x,i0)', 'test_scifort: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_normal(failures)
        integer, intent(inout) :: failures

        integer :: i
        real(dp) :: p(7)
        real(dp) :: x

        call check_close('normal pdf 0', normal_pdf(0.0_dp), &
            0.39894228040143267794_dp, 2.0e-16_dp, 2.0e-15_dp, failures)
        call check_close('normal cdf 0', normal_cdf(0.0_dp), &
            0.5_dp, 0.0_dp, 0.0_dp, failures)
        call check_close('normal sf 0', normal_sf(0.0_dp), &
            0.5_dp, 0.0_dp, 0.0_dp, failures)
        call check_close('normal ppf 0.975', normal_ppf(0.975_dp), &
            1.9599639845400542355_dp, 5.0e-15_dp, 5.0e-15_dp, failures)
        call check_close('normal logcdf -10', normal_logcdf(-10.0_dp), &
            -53.231285150512470578_dp, 2.0e-13_dp, 5.0e-15_dp, failures)
        call check_close('normal logsf 10', normal_logsf(10.0_dp), &
            -53.231285150512470578_dp, 2.0e-13_dp, 5.0e-15_dp, failures)
        call check_close('normal loc scale', normal_cdf(3.0_dp, 1.0_dp, 2.0_dp), &
            normal_cdf(1.0_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('normal extreme log round trip', &
            normal_logcdf(normal_ppf(1.0e-300_dp)), log(1.0e-300_dp), &
            3.0e-13_dp, 2.0e-15_dp, failures)

        p = [1.0e-12_dp, 1.0e-6_dp, 0.1_dp, 0.5_dp, 0.9_dp, &
            1.0_dp - 1.0e-6_dp, 1.0_dp - 1.0e-12_dp]
        do i = 1, size(p)
            x = normal_ppf(p(i))
            call check_close('normal cdf ppf round trip', normal_cdf(x), p(i), &
                3.0e-16_dp, 2.0e-13_dp, failures)
            call check_close('normal sf isf round trip', normal_sf(normal_isf(p(i))), &
                p(i), 3.0e-16_dp, 2.0e-13_dp, failures)
        end do

        call check_true('normal ppf zero is negative infinity', &
            .not. ieee_is_finite(normal_ppf(0.0_dp)) .and. &
            normal_ppf(0.0_dp) < 0.0_dp, failures)
        call check_true('normal ppf one is positive infinity', &
            .not. ieee_is_finite(normal_ppf(1.0_dp)) .and. &
            normal_ppf(1.0_dp) > 0.0_dp, failures)
    end subroutine test_normal

    subroutine test_uniform(failures)
        integer, intent(inout) :: failures

        call check_close('uniform pdf center', uniform_pdf(0.5_dp), &
            1.0_dp, 0.0_dp, 0.0_dp, failures)
        call check_close('uniform cdf center', uniform_cdf(0.5_dp), &
            0.5_dp, 0.0_dp, 0.0_dp, failures)
        call check_close('uniform sf center', uniform_sf(0.5_dp), &
            0.5_dp, 0.0_dp, 0.0_dp, failures)
        call check_close('uniform ppf', uniform_ppf(0.25_dp, 2.0_dp, 4.0_dp), &
            3.0_dp, 0.0_dp, 0.0_dp, failures)
        call check_close('uniform isf', uniform_isf(0.25_dp, 2.0_dp, 4.0_dp), &
            5.0_dp, 0.0_dp, 0.0_dp, failures)
        call check_close('uniform outside pdf', uniform_pdf(-1.0_dp), &
            0.0_dp, 0.0_dp, 0.0_dp, failures)
    end subroutine test_uniform

    subroutine test_exponential(failures)
        integer, intent(inout) :: failures

        call check_close('exponential pdf zero', exponential_pdf(0.0_dp), &
            1.0_dp, 0.0_dp, 0.0_dp, failures)
        call check_close('exponential cdf one', exponential_cdf(1.0_dp), &
            0.63212055882855767840_dp, 2.0e-16_dp, 2.0e-15_dp, failures)
        call check_close('exponential sf one', exponential_sf(1.0_dp), &
            exp(-1.0_dp), 2.0e-16_dp, 2.0e-15_dp, failures)
        call check_close('exponential ppf cdf', &
            exponential_cdf(exponential_ppf(0.9_dp)), 0.9_dp, &
            2.0e-16_dp, 2.0e-15_dp, failures)
        call check_close('exponential isf sf', &
            exponential_sf(exponential_isf(1.0e-12_dp)), 1.0e-12_dp, &
            1.0e-27_dp, 2.0e-14_dp, failures)
        call check_close('exponential below support pdf', &
            exponential_pdf(-1.0_dp), 0.0_dp, 0.0_dp, 0.0_dp, failures)
    end subroutine test_exponential

    subroutine test_laplace(failures)
        integer, intent(inout) :: failures

        call check_close('laplace pdf zero', laplace_pdf(0.0_dp), &
            0.5_dp, 0.0_dp, 0.0_dp, failures)
        call check_close('laplace cdf minus one', laplace_cdf(-1.0_dp), &
            0.18393972058572116080_dp, 2.0e-16_dp, 2.0e-15_dp, failures)
        call check_close('laplace cdf one', laplace_cdf(1.0_dp), &
            0.81606027941427883920_dp, 2.0e-16_dp, 2.0e-15_dp, failures)
        call check_close('laplace ppf round trip', &
            laplace_cdf(laplace_ppf(0.01_dp)), 0.01_dp, &
            2.0e-17_dp, 2.0e-15_dp, failures)
        call check_close('laplace isf round trip', &
            laplace_sf(laplace_isf(0.01_dp)), 0.01_dp, &
            2.0e-17_dp, 2.0e-15_dp, failures)
    end subroutine test_laplace

    subroutine test_logistic(failures)
        integer, intent(inout) :: failures

        call check_close('logistic pdf zero', logistic_pdf(0.0_dp), &
            0.25_dp, 2.0e-16_dp, 2.0e-15_dp, failures)
        call check_close('logistic cdf zero', logistic_cdf(0.0_dp), &
            0.5_dp, 0.0_dp, 0.0_dp, failures)
        call check_close('logistic cdf ppf', &
            logistic_cdf(logistic_ppf(1.0e-10_dp)), 1.0e-10_dp, &
            1.0e-24_dp, 2.0e-14_dp, failures)
        call check_close('logistic sf isf', &
            logistic_sf(logistic_isf(1.0e-10_dp)), 1.0e-10_dp, &
            1.0e-24_dp, 2.0e-14_dp, failures)
        call check_close('logistic logcdf tail', logistic_logcdf(-100.0_dp), &
            -100.0_dp, 1.0e-14_dp, 1.0e-15_dp, failures)
        call check_close('logistic logsf tail', logistic_logsf(100.0_dp), &
            -100.0_dp, 1.0e-14_dp, 1.0e-15_dp, failures)
    end subroutine test_logistic

    subroutine test_cauchy(failures)
        integer, intent(inout) :: failures

        call check_close('cauchy pdf zero', cauchy_pdf(0.0_dp), &
            0.31830988618379067154_dp, 2.0e-16_dp, 2.0e-15_dp, failures)
        call check_close('cauchy cdf minus one', cauchy_cdf(-1.0_dp), &
            0.25_dp, 2.0e-16_dp, 2.0e-15_dp, failures)
        call check_close('cauchy cdf one', cauchy_cdf(1.0_dp), &
            0.75_dp, 2.0e-16_dp, 2.0e-15_dp, failures)
        call check_close('cauchy ppf quarter', cauchy_ppf(0.25_dp), &
            -1.0_dp, 3.0e-16_dp, 2.0e-15_dp, failures)
        call check_close('cauchy isf quarter', cauchy_isf(0.25_dp), &
            1.0_dp, 3.0e-16_dp, 2.0e-15_dp, failures)
        call check_close('cauchy tail round trip', &
            cauchy_cdf(cauchy_ppf(1.0e-12_dp)), 1.0e-12_dp, &
            1.0e-27_dp, 2.0e-14_dp, failures)
        call check_close('cauchy right logcdf tail', cauchy_logcdf(1.0e20_dp), &
            -3.1830988618379067154e-21_dp, 2.0e-36_dp, 2.0e-15_dp, failures)
        call check_close('cauchy left logsf tail', cauchy_logsf(-1.0e20_dp), &
            -3.1830988618379067154e-21_dp, 2.0e-36_dp, 2.0e-15_dp, failures)
    end subroutine test_cauchy

    subroutine test_invalid_inputs(failures)
        integer, intent(inout) :: failures

        call check_true('normal invalid scale', &
            ieee_is_nan(normal_pdf(0.0_dp, scale=0.0_dp)), failures)
        call check_true('uniform invalid scale', &
            ieee_is_nan(uniform_cdf(0.0_dp, scale=-1.0_dp)), failures)
        call check_true('exponential invalid p', &
            ieee_is_nan(exponential_ppf(1.1_dp)), failures)
        call check_true('laplace invalid p', &
            ieee_is_nan(laplace_ppf(-0.1_dp)), failures)
        call check_true('logistic invalid scale', &
            ieee_is_nan(logistic_pdf(0.0_dp, scale=0.0_dp)), failures)
        call check_true('cauchy invalid p', &
            ieee_is_nan(cauchy_ppf(2.0_dp)), failures)
    end subroutine test_invalid_inputs

    subroutine test_ieee_inputs(failures)
        integer, intent(inout) :: failures

        real(dp) :: nan_value
        real(dp) :: negative_inf
        real(dp) :: positive_inf

        nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
        negative_inf = ieee_value(0.0_dp, ieee_negative_inf)
        positive_inf = ieee_value(0.0_dp, ieee_positive_inf)

        call check_true('normal NaN propagation', &
            ieee_is_nan(normal_cdf(nan_value)), failures)
        call check_true('uniform NaN propagation', &
            ieee_is_nan(uniform_pdf(nan_value)), failures)
        call check_true('exponential NaN propagation', &
            ieee_is_nan(exponential_cdf(nan_value)), failures)
        call check_true('laplace NaN propagation', &
            ieee_is_nan(laplace_cdf(nan_value)), failures)
        call check_true('logistic NaN propagation', &
            ieee_is_nan(logistic_cdf(nan_value)), failures)
        call check_true('cauchy NaN propagation', &
            ieee_is_nan(cauchy_cdf(nan_value)), failures)

        call check_close('normal cdf negative infinity', normal_cdf(negative_inf), &
            0.0_dp, 0.0_dp, 0.0_dp, failures)
        call check_close('normal cdf positive infinity', normal_cdf(positive_inf), &
            1.0_dp, 0.0_dp, 0.0_dp, failures)
        call check_close('exponential cdf negative infinity', &
            exponential_cdf(negative_inf), 0.0_dp, 0.0_dp, 0.0_dp, failures)
        call check_close('exponential cdf positive infinity', &
            exponential_cdf(positive_inf), 1.0_dp, 0.0_dp, 0.0_dp, failures)
        call check_close('laplace cdf negative infinity', laplace_cdf(negative_inf), &
            0.0_dp, 0.0_dp, 0.0_dp, failures)
        call check_close('laplace cdf positive infinity', laplace_cdf(positive_inf), &
            1.0_dp, 0.0_dp, 0.0_dp, failures)
        call check_close('logistic cdf negative infinity', &
            logistic_cdf(negative_inf), 0.0_dp, 0.0_dp, 0.0_dp, failures)
        call check_close('logistic cdf positive infinity', &
            logistic_cdf(positive_inf), 1.0_dp, 0.0_dp, 0.0_dp, failures)
        call check_close('cauchy cdf negative infinity', cauchy_cdf(negative_inf), &
            0.0_dp, 0.0_dp, 0.0_dp, failures)
        call check_close('cauchy cdf positive infinity', cauchy_cdf(positive_inf), &
            1.0_dp, 0.0_dp, 0.0_dp, failures)
    end subroutine test_ieee_inputs

    subroutine test_elemental_arrays(failures)
        integer, intent(inout) :: failures

        real(dp) :: actual(3)
        real(dp) :: expected(3)
        real(dp) :: x(3)

        x = [-1.0_dp, 0.0_dp, 1.0_dp]
        actual = normal_cdf(x)
        expected = [normal_cdf(-1.0_dp), 0.5_dp, normal_cdf(1.0_dp)]
        call check_array_close('normal elemental array', actual, expected, &
            0.0_dp, 0.0_dp, failures)

        actual = logistic_cdf(x, scale=2.0_dp)
        expected = [logistic_cdf(-1.0_dp, scale=2.0_dp), 0.5_dp, &
            logistic_cdf(1.0_dp, scale=2.0_dp)]
        call check_array_close('logistic elemental array', actual, expected, &
            0.0_dp, 0.0_dp, failures)
    end subroutine test_elemental_arrays

    subroutine test_distribution_identities(failures)
        integer, intent(inout) :: failures

        integer :: i
        real(dp) :: p(6)
        real(dp) :: x(5)

        p = [1.0e-8_dp, 0.01_dp, 0.25_dp, 0.5_dp, 0.9_dp, &
            1.0_dp - 1.0e-8_dp]
        do i = 1, size(p)
            call check_close('normal cdf-ppf identity', &
                normal_cdf(normal_ppf(p(i))), p(i), &
                4.0e-16_dp, 3.0e-13_dp, failures)
            call check_close('uniform cdf-ppf identity', &
                uniform_cdf(uniform_ppf(p(i))), p(i), &
                2.0e-16_dp, 2.0e-15_dp, failures)
            call check_close('exponential cdf-ppf identity', &
                exponential_cdf(exponential_ppf(p(i))), p(i), &
                2.0e-16_dp, 3.0e-14_dp, failures)
            call check_close('laplace cdf-ppf identity', &
                laplace_cdf(laplace_ppf(p(i))), p(i), &
                2.0e-16_dp, 3.0e-14_dp, failures)
            call check_close('logistic cdf-ppf identity', &
                logistic_cdf(logistic_ppf(p(i))), p(i), &
                2.0e-16_dp, 3.0e-14_dp, failures)
            call check_close('cauchy cdf-ppf identity', &
                cauchy_cdf(cauchy_ppf(p(i))), p(i), &
                2.0e-16_dp, 3.0e-14_dp, failures)
        end do

        x = [-10.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 10.0_dp]
        do i = 1, size(x)
            call check_close('normal cdf plus sf', &
                normal_cdf(x(i)) + normal_sf(x(i)), 1.0_dp, &
                2.0e-16_dp, 2.0e-15_dp, failures)
            call check_close('uniform cdf plus sf', &
                uniform_cdf(x(i)) + uniform_sf(x(i)), 1.0_dp, &
                0.0_dp, 0.0_dp, failures)
            call check_close('exponential cdf plus sf', &
                exponential_cdf(x(i)) + exponential_sf(x(i)), 1.0_dp, &
                2.0e-16_dp, 2.0e-15_dp, failures)
            call check_close('laplace cdf plus sf', &
                laplace_cdf(x(i)) + laplace_sf(x(i)), 1.0_dp, &
                2.0e-16_dp, 2.0e-15_dp, failures)
            call check_close('logistic cdf plus sf', &
                logistic_cdf(x(i)) + logistic_sf(x(i)), 1.0_dp, &
                2.0e-16_dp, 2.0e-15_dp, failures)
            call check_close('cauchy cdf plus sf', &
                cauchy_cdf(x(i)) + cauchy_sf(x(i)), 1.0_dp, &
                2.0e-16_dp, 2.0e-15_dp, failures)
        end do
    end subroutine test_distribution_identities

    subroutine test_c_api(failures)
        integer, intent(inout) :: failures

        integer(c_int) :: status
        real(c_double) :: x(3)
        real(c_double) :: y(3)

        call check_true('C API version major', scifort_version_major() == 0_c_int, failures)
        call check_true('C API version minor', scifort_version_minor() == 1_c_int, failures)
        call check_true('C API version patch', scifort_version_patch() == 0_c_int, failures)
        call check_close('C API scalar normal cdf', &
            real(scifort_normal_cdf_f64(0.0_c_double, 0.0_c_double, &
            1.0_c_double), dp), 0.5_dp, 0.0_dp, 0.0_dp, failures)

        x = [-1.0_c_double, 0.0_c_double, 1.0_c_double]
        call scifort_normal_cdf_vec_f64(3_c_size_t, x, 0.0_c_double, &
            1.0_c_double, y, status)
        call check_true('C API vector status', status == 0_c_int, failures)
        call check_close('C API vector center', real(y(2), dp), &
            0.5_dp, 0.0_dp, 0.0_dp, failures)

        call scifort_normal_cdf_vec_f64(3_c_size_t, x, 0.0_c_double, &
            0.0_c_double, y, status)
        call check_true('C API invalid status', status == 1_c_int, failures)
    end subroutine test_c_api

    subroutine check_close(name, actual, expected, atol, rtol, failures)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual
        real(dp), intent(in) :: expected
        real(dp), intent(in) :: atol
        real(dp), intent(in) :: rtol
        integer, intent(inout) :: failures

        real(dp) :: tolerance

        tolerance = atol + rtol * abs(expected)
        if (ieee_is_nan(actual) .or. abs(actual - expected) > tolerance) then
            failures = failures + 1
            print '(a,1x,a)', 'FAIL:', name
            print '(a,1x,es24.16)', '  actual  =', actual
            print '(a,1x,es24.16)', '  expected=', expected
            print '(a,1x,es24.16)', '  tolerance=', tolerance
        end if
    end subroutine check_close

    subroutine check_array_close(name, actual, expected, atol, rtol, failures)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual(:)
        real(dp), intent(in) :: expected(:)
        real(dp), intent(in) :: atol
        real(dp), intent(in) :: rtol
        integer, intent(inout) :: failures

        integer :: i

        if (size(actual) /= size(expected)) then
            failures = failures + 1
            print '(a,1x,a)', 'FAIL:', name
            print '(a)', '  array sizes differ'
            return
        end if

        do i = 1, size(actual)
            call check_close(name, actual(i), expected(i), atol, rtol, failures)
        end do
    end subroutine check_array_close

    subroutine check_true(name, condition, failures)
        character(len=*), intent(in) :: name
        logical, intent(in) :: condition
        integer, intent(inout) :: failures

        if (.not. condition) then
            failures = failures + 1
            print '(a,1x,a)', 'FAIL:', name
        end if
    end subroutine check_true

end program test_scifort
