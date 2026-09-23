! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Full-stack regression coverage for Bradford, truncated exponential,
! Fisk/log-logistic, and double-Weibull distributions.

program test_closed_form_continuous
    use scifort_c_api
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state, rng_uniform
    use scifort_stats
    use test_checks, only : check_close, check_true
    use test_closed_form_reference
    implicit none

    integer :: failures

    failures = 0
    call test_reference_values(failures)
    call test_tail_identities(failures)
    call test_likelihood_scores(failures)
    call test_random_fit_c_api(failures)

    if (failures == 0) then
        print '(a)', 'test_closed_form_continuous: PASS'
    else
        print '(a,1x,i0)', 'test_closed_form_continuous: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference_values(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        integer :: i
        character(len=64) :: label

        do i = 1, size(br_x)
            write(label, '(a,i0)') 'bradford pdf ref ', i
            call check_close(trim(label), bradford_pdf(br_x(i), reference_shape), &
                br_pdf(i), 3.0e-14_dp, 3.0e-13_dp, failures)
            write(label, '(a,i0)') 'bradford cdf ref ', i
            call check_close(trim(label), bradford_cdf(br_x(i), reference_shape), &
                br_cdf(i), 3.0e-14_dp, 3.0e-13_dp, failures)
            write(label, '(a,i0)') 'bradford sf ref ', i
            call check_close(trim(label), bradford_sf(br_x(i), reference_shape), &
                br_sf(i), 3.0e-14_dp, 3.0e-13_dp, failures)
            write(label, '(a,i0)') 'bradford logcdf ref ', i
            call check_close(trim(label), bradford_logcdf(br_x(i), reference_shape), &
                br_logcdf(i), 3.0e-14_dp, 3.0e-13_dp, failures)
            write(label, '(a,i0)') 'bradford logsf ref ', i
            call check_close(trim(label), bradford_logsf(br_x(i), reference_shape), &
                br_logsf(i), 3.0e-14_dp, 3.0e-13_dp, failures)

            write(label, '(a,i0)') 'truncexpon pdf ref ', i
            call check_close(trim(label), truncexpon_pdf(te_x(i), reference_shape), &
                te_pdf(i), 3.0e-14_dp, 3.0e-13_dp, failures)
            write(label, '(a,i0)') 'truncexpon cdf ref ', i
            call check_close(trim(label), truncexpon_cdf(te_x(i), reference_shape), &
                te_cdf(i), 3.0e-14_dp, 3.0e-13_dp, failures)
            write(label, '(a,i0)') 'truncexpon sf ref ', i
            call check_close(trim(label), truncexpon_sf(te_x(i), reference_shape), &
                te_sf(i), 3.0e-14_dp, 3.0e-13_dp, failures)
            write(label, '(a,i0)') 'truncexpon logcdf ref ', i
            call check_close(trim(label), truncexpon_logcdf(te_x(i), reference_shape), &
                te_logcdf(i), 3.0e-14_dp, 3.0e-13_dp, failures)
            write(label, '(a,i0)') 'truncexpon logsf ref ', i
            call check_close(trim(label), truncexpon_logsf(te_x(i), reference_shape), &
                te_logsf(i), 3.0e-14_dp, 3.0e-13_dp, failures)

            write(label, '(a,i0)') 'fisk pdf ref ', i
            call check_close(trim(label), fisk_pdf(fi_x(i), reference_shape), &
                fi_pdf(i), 3.0e-14_dp, 3.0e-13_dp, failures)
            write(label, '(a,i0)') 'fisk cdf ref ', i
            call check_close(trim(label), fisk_cdf(fi_x(i), reference_shape), &
                fi_cdf(i), 3.0e-14_dp, 3.0e-13_dp, failures)
            write(label, '(a,i0)') 'fisk sf ref ', i
            call check_close(trim(label), fisk_sf(fi_x(i), reference_shape), &
                fi_sf(i), 3.0e-14_dp, 3.0e-13_dp, failures)
            write(label, '(a,i0)') 'fisk logcdf ref ', i
            call check_close(trim(label), fisk_logcdf(fi_x(i), reference_shape), &
                fi_logcdf(i), 3.0e-14_dp, 3.0e-13_dp, failures)
            write(label, '(a,i0)') 'fisk logsf ref ', i
            call check_close(trim(label), fisk_logsf(fi_x(i), reference_shape), &
                fi_logsf(i), 3.0e-14_dp, 3.0e-13_dp, failures)

            write(label, '(a,i0)') 'dweibull pdf ref ', i
            call check_close(trim(label), dweibull_pdf(dw_x(i), reference_shape), &
                dw_pdf(i), 3.0e-14_dp, 3.0e-13_dp, failures)
            write(label, '(a,i0)') 'dweibull cdf ref ', i
            call check_close(trim(label), dweibull_cdf(dw_x(i), reference_shape), &
                dw_cdf(i), 3.0e-14_dp, 3.0e-13_dp, failures)
            write(label, '(a,i0)') 'dweibull sf ref ', i
            call check_close(trim(label), dweibull_sf(dw_x(i), reference_shape), &
                dw_sf(i), 3.0e-14_dp, 3.0e-13_dp, failures)
            write(label, '(a,i0)') 'dweibull logcdf ref ', i
            call check_close(trim(label), dweibull_logcdf(dw_x(i), reference_shape), &
                dw_logcdf(i), 3.0e-14_dp, 3.0e-13_dp, failures)
            write(label, '(a,i0)') 'dweibull logsf ref ', i
            call check_close(trim(label), dweibull_logsf(dw_x(i), reference_shape), &
                dw_logsf(i), 3.0e-14_dp, 3.0e-13_dp, failures)
        end do

        do i = 1, size(reference_probs)
            write(label, '(a,i0)') 'bradford ppf ref ', i
            call check_close(trim(label), bradford_ppf(reference_probs(i), reference_shape), &
                br_ppf(i), 3.0e-14_dp, 3.0e-13_dp, failures)
            write(label, '(a,i0)') 'bradford isf ref ', i
            call check_close(trim(label), bradford_isf(reference_probs(i), reference_shape), &
                br_isf(i), 3.0e-14_dp, 3.0e-13_dp, failures)
            write(label, '(a,i0)') 'truncexpon ppf ref ', i
            call check_close(trim(label), truncexpon_ppf(reference_probs(i), reference_shape), &
                te_ppf(i), 3.0e-14_dp, 3.0e-13_dp, failures)
            write(label, '(a,i0)') 'truncexpon isf ref ', i
            call check_close(trim(label), truncexpon_isf(reference_probs(i), reference_shape), &
                te_isf(i), 3.0e-14_dp, 3.0e-13_dp, failures)
            write(label, '(a,i0)') 'fisk ppf ref ', i
            call check_close(trim(label), fisk_ppf(reference_probs(i), reference_shape), &
                fi_ppf(i), 3.0e-14_dp, 3.0e-13_dp, failures)
            write(label, '(a,i0)') 'fisk isf ref ', i
            call check_close(trim(label), fisk_isf(reference_probs(i), reference_shape), &
                fi_isf(i), 3.0e-14_dp, 3.0e-13_dp, failures)
            write(label, '(a,i0)') 'dweibull ppf ref ', i
            call check_close(trim(label), dweibull_ppf(reference_probs(i), reference_shape), &
                dw_ppf(i), 3.0e-14_dp, 3.0e-13_dp, failures)
            write(label, '(a,i0)') 'dweibull isf ref ', i
            call check_close(trim(label), dweibull_isf(reference_probs(i), reference_shape), &
                dw_isf(i), 3.0e-14_dp, 3.0e-13_dp, failures)
        end do
    end subroutine test_reference_values

    subroutine test_tail_identities(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: probs(5) = [ &
            1.0e-8_dp, 0.1_dp, 0.5_dp, 0.9_dp, 1.0_dp - 1.0e-8_dp]
        integer :: i
        real(dp) :: x

        do i = 1, size(probs)
            x = bradford_ppf(probs(i), 1.6_dp, -0.3_dp, 1.4_dp)
            call check_close('bradford ppf/cdf', bradford_cdf(x, 1.6_dp, -0.3_dp, 1.4_dp), &
                probs(i), 3.0e-12_dp, 3.0e-9_dp, failures)
            x = truncexpon_isf(probs(i), 1.8_dp, -0.3_dp, 1.4_dp)
            call check_close('truncexpon isf/sf', truncexpon_sf(x, 1.8_dp, -0.3_dp, 1.4_dp), &
                probs(i), 3.0e-12_dp, 3.0e-9_dp, failures)
            x = fisk_ppf(probs(i), 1.6_dp, -0.3_dp, 1.4_dp)
            call check_close('fisk ppf/cdf', fisk_cdf(x, 1.6_dp, -0.3_dp, 1.4_dp), &
                probs(i), 3.0e-12_dp, 3.0e-9_dp, failures)
            x = dweibull_isf(probs(i), 1.8_dp, -0.3_dp, 1.4_dp)
            call check_close('dweibull isf/sf', dweibull_sf(x, 1.8_dp, -0.3_dp, 1.4_dp), &
                probs(i), 3.0e-12_dp, 3.0e-9_dp, failures)
        end do
    end subroutine test_tail_identities

    subroutine test_likelihood_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: xb(5) = [0.2_dp, 0.4_dp, 0.7_dp, 1.0_dp, 1.2_dp]
        real(dp), parameter :: xt(5) = [0.2_dp, 0.5_dp, 0.9_dp, 1.3_dp, 1.8_dp]
        real(dp), parameter :: xf(5) = [0.3_dp, 0.7_dp, 1.1_dp, 2.0_dp, 3.0_dp]
        real(dp), parameter :: xd(5) = [-2.0_dp, -0.8_dp, 0.4_dp, 1.3_dp, 2.5_dp]
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp), parameter :: tol = 5.0e-5_dp
        real(dp) :: s3(3)

        s3 = bradford_score(xb, 1.7_dp, 0.1_dp, 1.2_dp)
        call check_close('bradford score c', s3(1), centered( &
            bradford_loglikelihood(xb, 1.7_dp + h, 0.1_dp, 1.2_dp), &
            bradford_loglikelihood(xb, 1.7_dp - h, 0.1_dp, 1.2_dp), h), tol, tol, failures)
        call check_close('bradford score loc', s3(2), centered( &
            bradford_loglikelihood(xb, 1.7_dp, 0.1_dp + h, 1.2_dp), &
            bradford_loglikelihood(xb, 1.7_dp, 0.1_dp - h, 1.2_dp), h), tol, tol, failures)
        call check_close('bradford score scale', s3(3), centered( &
            bradford_loglikelihood(xb, 1.7_dp, 0.1_dp, 1.2_dp + h), &
            bradford_loglikelihood(xb, 1.7_dp, 0.1_dp, 1.2_dp - h), h), tol, tol, failures)

        s3 = truncexpon_score(xt, 1.7_dp, 0.1_dp, 1.2_dp)
        call check_close('truncexpon score b', s3(1), centered( &
            truncexpon_loglikelihood(xt, 1.7_dp + h, 0.1_dp, 1.2_dp), &
            truncexpon_loglikelihood(xt, 1.7_dp - h, 0.1_dp, 1.2_dp), h), tol, tol, failures)
        call check_close('truncexpon score loc', s3(2), centered( &
            truncexpon_loglikelihood(xt, 1.7_dp, 0.1_dp + h, 1.2_dp), &
            truncexpon_loglikelihood(xt, 1.7_dp, 0.1_dp - h, 1.2_dp), h), tol, tol, failures)
        call check_close('truncexpon score scale', s3(3), centered( &
            truncexpon_loglikelihood(xt, 1.7_dp, 0.1_dp, 1.2_dp + h), &
            truncexpon_loglikelihood(xt, 1.7_dp, 0.1_dp, 1.2_dp - h), h), tol, tol, failures)

        s3 = fisk_score(xf, 1.7_dp, 0.1_dp, 1.2_dp)
        call check_close('fisk score c', s3(1), centered( &
            fisk_loglikelihood(xf, 1.7_dp + h, 0.1_dp, 1.2_dp), &
            fisk_loglikelihood(xf, 1.7_dp - h, 0.1_dp, 1.2_dp), h), tol, tol, failures)
        call check_close('fisk score loc', s3(2), centered( &
            fisk_loglikelihood(xf, 1.7_dp, 0.1_dp + h, 1.2_dp), &
            fisk_loglikelihood(xf, 1.7_dp, 0.1_dp - h, 1.2_dp), h), tol, tol, failures)
        call check_close('fisk score scale', s3(3), centered( &
            fisk_loglikelihood(xf, 1.7_dp, 0.1_dp, 1.2_dp + h), &
            fisk_loglikelihood(xf, 1.7_dp, 0.1_dp, 1.2_dp - h), h), tol, tol, failures)

        s3 = dweibull_score(xd, 1.7_dp, 0.1_dp, 1.2_dp)
        call check_close('dweibull score c', s3(1), centered( &
            dweibull_loglikelihood(xd, 1.7_dp + h, 0.1_dp, 1.2_dp), &
            dweibull_loglikelihood(xd, 1.7_dp - h, 0.1_dp, 1.2_dp), h), tol, tol, failures)
        call check_close('dweibull score loc', s3(2), centered( &
            dweibull_loglikelihood(xd, 1.7_dp, 0.1_dp + h, 1.2_dp), &
            dweibull_loglikelihood(xd, 1.7_dp, 0.1_dp - h, 1.2_dp), h), tol, tol, failures)
        call check_close('dweibull score scale', s3(3), centered( &
            dweibull_loglikelihood(xd, 1.7_dp, 0.1_dp, 1.2_dp + h), &
            dweibull_loglikelihood(xd, 1.7_dp, 0.1_dp, 1.2_dp - h), h), tol, tol, failures)
    end subroutine test_likelihood_scores

    subroutine test_random_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: xb(5) = [0.2_dp, 0.4_dp, 0.7_dp, 1.0_dp, 1.2_dp]
        real(dp), parameter :: xt(5) = [0.2_dp, 0.5_dp, 0.9_dp, 1.3_dp, 1.8_dp]
        real(dp), parameter :: xf(5) = [0.3_dp, 0.7_dp, 1.1_dp, 2.0_dp, 3.0_dp]
        real(dp), parameter :: xd(5) = [-2.0_dp, -0.8_dp, 0.4_dp, 1.3_dp, 2.5_dp]
        real(dp) :: actual
        real(dp) :: expected
        type(fit_result) :: result
        type(rng_state) :: reference
        type(rng_state) :: state

        call rng_seed(state, 97531)
        call rng_seed(reference, 97531)
        actual = bradford_rvs(state, 1.6_dp, -0.2_dp, 1.3_dp)
        expected = bradford_ppf(rng_uniform(reference), 1.6_dp, -0.2_dp, 1.3_dp)
        call check_close('rvs bradford', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = truncexpon_rvs(state, 1.8_dp, -0.2_dp, 1.3_dp)
        expected = truncexpon_ppf(rng_uniform(reference), 1.8_dp, -0.2_dp, 1.3_dp)
        call check_close('rvs truncexpon', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = fisk_rvs(state, 1.6_dp, -0.2_dp, 1.3_dp)
        expected = fisk_ppf(rng_uniform(reference), 1.6_dp, -0.2_dp, 1.3_dp)
        call check_close('rvs fisk', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = dweibull_rvs(state, 1.8_dp, -0.2_dp, 1.3_dp)
        expected = dweibull_ppf(rng_uniform(reference), 1.8_dp, -0.2_dp, 1.3_dp)
        call check_close('rvs dweibull', actual, expected, 0.0_dp, 0.0_dp, failures)

        call bradford_fit(xb, [1.7_dp, 0.1_dp, 1.2_dp], &
            [1.7_dp, 0.1_dp, 1.2_dp], result)
        call check_true('fit bradford fixed success', result%success, failures)
        call truncexpon_fit(xt, [1.7_dp, 0.1_dp, 1.2_dp], &
            [1.7_dp, 0.1_dp, 1.2_dp], result)
        call check_true('fit truncexpon fixed success', result%success, failures)
        call fisk_fit(xf, [1.7_dp, 0.1_dp, 1.2_dp], &
            [1.7_dp, 0.1_dp, 1.2_dp], result)
        call check_true('fit fisk fixed success', result%success, failures)
        call dweibull_fit(xd, [1.7_dp, 0.1_dp, 1.2_dp], &
            [1.7_dp, 0.1_dp, 1.2_dp], result)
        call check_true('fit dweibull fixed success', result%success, failures)

        call check_close('C ABI bradford cdf', &
            scifort_bradford_cdf_f64(0.6_dp, 1.7_dp, -0.2_dp, 1.3_dp), &
            bradford_cdf(0.6_dp, 1.7_dp, -0.2_dp, 1.3_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI truncexpon ppf', &
            scifort_truncexpon_ppf_f64(0.4_dp, 1.7_dp, -0.2_dp, 1.3_dp), &
            truncexpon_ppf(0.4_dp, 1.7_dp, -0.2_dp, 1.3_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI fisk pdf', &
            scifort_fisk_pdf_f64(1.2_dp, 1.7_dp, -0.2_dp, 1.3_dp), &
            fisk_pdf(1.2_dp, 1.7_dp, -0.2_dp, 1.3_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI dweibull cdf', &
            scifort_dweibull_cdf_f64(0.8_dp, 1.7_dp, -0.2_dp, 1.3_dp), &
            dweibull_cdf(0.8_dp, 1.7_dp, -0.2_dp, 1.3_dp), 0.0_dp, 0.0_dp, failures)
    end subroutine test_random_fit_c_api

    pure function centered(plus, minus, h) result(value)
        real(dp), intent(in) :: plus !! objective at theta + h
        real(dp), intent(in) :: minus !! objective at theta - h
        real(dp), intent(in) :: h !! positive centered-difference half step
        real(dp) :: value

        value = (plus - minus) / (2.0_dp * h)
    end function centered

end program test_closed_form_continuous
