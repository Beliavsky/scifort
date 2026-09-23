! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Full-stack regression coverage for Gompertz, inverse Weibull, beta-prime, and Burr XII.
program test_positive_shape_continuous
    use scifort_c_api
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state, rng_uniform
    use scifort_stats
    use test_checks, only : check_close, check_true
    use test_positive_shape_reference
    implicit none

    integer :: failures

    failures = 0
    call test_reference_values(failures)
    call test_tail_identities(failures)
    call test_likelihood_scores(failures)
    call test_random_fit_c_api(failures)

    if (failures == 0) then
        print '(a)', 'test_positive_shape_continuous: PASS'
    else
        print '(a,1x,i0)', 'test_positive_shape_continuous: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference_values(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        integer :: i
        character(len=64) :: label

        do i = 1, 5
            write(label, '(a,i0)') 'gompertz pdf ref ', i
            call check_close(trim(label), gompertz_pdf(go_x(i), reference_c), &
                go_pdf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            write(label, '(a,i0)') 'gompertz cdf ref ', i
            call check_close(trim(label), gompertz_cdf(go_x(i), reference_c), &
                go_cdf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            write(label, '(a,i0)') 'gompertz sf ref ', i
            call check_close(trim(label), gompertz_sf(go_x(i), reference_c), &
                go_sf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            if (i == 1) then
                call check_true('gompertz logcdf endpoint', &
                    gompertz_logcdf(go_x(i), reference_c) < -huge(1.0_dp), failures)
            else
                write(label, '(a,i0)') 'gompertz logcdf ref ', i
                call check_close(trim(label), gompertz_logcdf(go_x(i), reference_c), &
                    go_logcdf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            end if
            write(label, '(a,i0)') 'gompertz logsf ref ', i
            call check_close(trim(label), gompertz_logsf(go_x(i), reference_c), &
                go_logsf(i), 5.0e-12_dp, 5.0e-11_dp, failures)

            write(label, '(a,i0)') 'invweibull pdf ref ', i
            call check_close(trim(label), invweibull_pdf(iw_x(i), reference_c), &
                iw_pdf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            write(label, '(a,i0)') 'invweibull cdf ref ', i
            call check_close(trim(label), invweibull_cdf(iw_x(i), reference_c), &
                iw_cdf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            write(label, '(a,i0)') 'invweibull sf ref ', i
            call check_close(trim(label), invweibull_sf(iw_x(i), reference_c), &
                iw_sf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            write(label, '(a,i0)') 'invweibull logcdf ref ', i
            call check_close(trim(label), invweibull_logcdf(iw_x(i), reference_c), &
                iw_logcdf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            write(label, '(a,i0)') 'invweibull logsf ref ', i
            call check_close(trim(label), invweibull_logsf(iw_x(i), reference_c), &
                iw_logsf(i), 5.0e-12_dp, 5.0e-11_dp, failures)

            write(label, '(a,i0)') 'betaprime pdf ref ', i
            call check_close(trim(label), betaprime_pdf(bp_x(i), reference_c, reference_d), &
                bp_pdf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            write(label, '(a,i0)') 'betaprime cdf ref ', i
            call check_close(trim(label), betaprime_cdf(bp_x(i), reference_c, reference_d), &
                bp_cdf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            write(label, '(a,i0)') 'betaprime sf ref ', i
            call check_close(trim(label), betaprime_sf(bp_x(i), reference_c, reference_d), &
                bp_sf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            write(label, '(a,i0)') 'betaprime logcdf ref ', i
            call check_close(trim(label), betaprime_logcdf(bp_x(i), reference_c, reference_d), &
                bp_logcdf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            write(label, '(a,i0)') 'betaprime logsf ref ', i
            call check_close(trim(label), betaprime_logsf(bp_x(i), reference_c, reference_d), &
                bp_logsf(i), 5.0e-12_dp, 5.0e-11_dp, failures)

            write(label, '(a,i0)') 'burr12 pdf ref ', i
            call check_close(trim(label), burr12_pdf(bu_x(i), reference_c, reference_d), &
                bu_pdf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            write(label, '(a,i0)') 'burr12 cdf ref ', i
            call check_close(trim(label), burr12_cdf(bu_x(i), reference_c, reference_d), &
                bu_cdf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            write(label, '(a,i0)') 'burr12 sf ref ', i
            call check_close(trim(label), burr12_sf(bu_x(i), reference_c, reference_d), &
                bu_sf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            write(label, '(a,i0)') 'burr12 logcdf ref ', i
            call check_close(trim(label), burr12_logcdf(bu_x(i), reference_c, reference_d), &
                bu_logcdf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            write(label, '(a,i0)') 'burr12 logsf ref ', i
            call check_close(trim(label), burr12_logsf(bu_x(i), reference_c, reference_d), &
                bu_logsf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
        end do

        do i = 1, 5
            call check_close('gompertz ppf ref', gompertz_ppf(reference_probs(i), reference_c), &
                go_ppf(i), 2.0e-10_dp, 2.0e-9_dp, failures)
            call check_close('gompertz isf ref', gompertz_isf(reference_probs(i), reference_c), &
                go_isf(i), 2.0e-10_dp, 2.0e-9_dp, failures)
            call check_close('invweibull ppf ref', &
                invweibull_ppf(reference_probs(i), reference_c), iw_ppf(i), &
                2.0e-10_dp, 2.0e-9_dp, failures)
            call check_close('invweibull isf ref', &
                invweibull_isf(reference_probs(i), reference_c), iw_isf(i), &
                2.0e-10_dp, 2.0e-9_dp, failures)
            call check_close('betaprime ppf ref', &
                betaprime_ppf(reference_probs(i), reference_c, reference_d), bp_ppf(i), &
                2.0e-10_dp, 2.0e-9_dp, failures)
            call check_close('betaprime isf ref', &
                betaprime_isf(reference_probs(i), reference_c, reference_d), bp_isf(i), &
                2.0e-10_dp, 2.0e-9_dp, failures)
            call check_close('burr12 ppf ref', &
                burr12_ppf(reference_probs(i), reference_c, reference_d), bu_ppf(i), &
                2.0e-10_dp, 2.0e-9_dp, failures)
            call check_close('burr12 isf ref', &
                burr12_isf(reference_probs(i), reference_c, reference_d), bu_isf(i), &
                2.0e-10_dp, 2.0e-9_dp, failures)
        end do
    end subroutine test_reference_values

    subroutine test_tail_identities(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: probs(5) = [1.0e-8_dp, 0.1_dp, 0.5_dp, 0.9_dp, &
            1.0_dp - 1.0e-8_dp]
        integer :: i
        real(dp) :: x

        do i = 1, size(probs)
            x = gompertz_ppf(probs(i), 1.6_dp, -0.2_dp, 1.3_dp)
            call check_close('gompertz ppf/cdf', gompertz_cdf(x, 1.6_dp, -0.2_dp, 1.3_dp), &
                probs(i), 2.0e-9_dp, 2.0e-7_dp, failures)
            x = invweibull_isf(probs(i), 1.6_dp, -0.2_dp, 1.3_dp)
            call check_close('invweibull isf/sf', invweibull_sf(x, 1.6_dp, -0.2_dp, 1.3_dp), &
                probs(i), 2.0e-9_dp, 2.0e-7_dp, failures)
            x = betaprime_ppf(probs(i), 1.6_dp, 2.2_dp, -0.2_dp, 1.3_dp)
            call check_close('betaprime ppf/cdf', &
                betaprime_cdf(x, 1.6_dp, 2.2_dp, -0.2_dp, 1.3_dp), &
                probs(i), 2.0e-9_dp, 2.0e-7_dp, failures)
            x = burr12_isf(probs(i), 1.6_dp, 2.2_dp, -0.2_dp, 1.3_dp)
            call check_close('burr12 isf/sf', burr12_sf(x, 1.6_dp, 2.2_dp, -0.2_dp, 1.3_dp), &
                probs(i), 2.0e-9_dp, 2.0e-7_dp, failures)
        end do
    end subroutine test_tail_identities

    subroutine test_likelihood_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: x1(5) = [0.3_dp, 0.6_dp, 1.0_dp, 1.7_dp, 2.8_dp]
        real(dp), parameter :: x2(5) = [0.25_dp, 0.5_dp, 1.0_dp, 2.0_dp, 4.0_dp]
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp), parameter :: tol = 2.0e-4_dp
        real(dp) :: s3(3), s4(4)

        s3 = gompertz_score(x1, 1.7_dp, 0.1_dp, 1.2_dp)
        call compare_score3('gompertz', s3, &
            centered(gompertz_loglikelihood(x1, 1.7_dp + h, 0.1_dp, 1.2_dp), &
                gompertz_loglikelihood(x1, 1.7_dp - h, 0.1_dp, 1.2_dp), h), &
            centered(gompertz_loglikelihood(x1, 1.7_dp, 0.1_dp + h, 1.2_dp), &
                gompertz_loglikelihood(x1, 1.7_dp, 0.1_dp - h, 1.2_dp), h), &
            centered(gompertz_loglikelihood(x1, 1.7_dp, 0.1_dp, 1.2_dp + h), &
                gompertz_loglikelihood(x1, 1.7_dp, 0.1_dp, 1.2_dp - h), h), tol, failures)

        s3 = invweibull_score(x2, 1.7_dp, 0.1_dp, 1.2_dp)
        call compare_score3('invweibull', s3, &
            centered(invweibull_loglikelihood(x2, 1.7_dp + h, 0.1_dp, 1.2_dp), &
                invweibull_loglikelihood(x2, 1.7_dp - h, 0.1_dp, 1.2_dp), h), &
            centered(invweibull_loglikelihood(x2, 1.7_dp, 0.1_dp + h, 1.2_dp), &
                invweibull_loglikelihood(x2, 1.7_dp, 0.1_dp - h, 1.2_dp), h), &
            centered(invweibull_loglikelihood(x2, 1.7_dp, 0.1_dp, 1.2_dp + h), &
                invweibull_loglikelihood(x2, 1.7_dp, 0.1_dp, 1.2_dp - h), h), tol, failures)

        s4 = betaprime_score(x2, 1.7_dp, 2.2_dp, 0.1_dp, 1.2_dp)
        call compare_score4('betaprime', s4, &
            centered(betaprime_loglikelihood(x2, 1.7_dp + h, 2.2_dp, 0.1_dp, 1.2_dp), &
                betaprime_loglikelihood(x2, 1.7_dp - h, 2.2_dp, 0.1_dp, 1.2_dp), h), &
            centered(betaprime_loglikelihood(x2, 1.7_dp, 2.2_dp + h, 0.1_dp, 1.2_dp), &
                betaprime_loglikelihood(x2, 1.7_dp, 2.2_dp - h, 0.1_dp, 1.2_dp), h), &
            centered(betaprime_loglikelihood(x2, 1.7_dp, 2.2_dp, 0.1_dp + h, 1.2_dp), &
                betaprime_loglikelihood(x2, 1.7_dp, 2.2_dp, 0.1_dp - h, 1.2_dp), h), &
            centered(betaprime_loglikelihood(x2, 1.7_dp, 2.2_dp, 0.1_dp, 1.2_dp + h), &
                betaprime_loglikelihood(x2, 1.7_dp, 2.2_dp, 0.1_dp, 1.2_dp - h), h), tol, failures)

        s4 = burr12_score(x2, 1.7_dp, 2.2_dp, 0.1_dp, 1.2_dp)
        call compare_score4('burr12', s4, &
            centered(burr12_loglikelihood(x2, 1.7_dp + h, 2.2_dp, 0.1_dp, 1.2_dp), &
                burr12_loglikelihood(x2, 1.7_dp - h, 2.2_dp, 0.1_dp, 1.2_dp), h), &
            centered(burr12_loglikelihood(x2, 1.7_dp, 2.2_dp + h, 0.1_dp, 1.2_dp), &
                burr12_loglikelihood(x2, 1.7_dp, 2.2_dp - h, 0.1_dp, 1.2_dp), h), &
            centered(burr12_loglikelihood(x2, 1.7_dp, 2.2_dp, 0.1_dp + h, 1.2_dp), &
                burr12_loglikelihood(x2, 1.7_dp, 2.2_dp, 0.1_dp - h, 1.2_dp), h), &
            centered(burr12_loglikelihood(x2, 1.7_dp, 2.2_dp, 0.1_dp, 1.2_dp + h), &
                burr12_loglikelihood(x2, 1.7_dp, 2.2_dp, 0.1_dp, 1.2_dp - h), h), tol, failures)
    end subroutine test_likelihood_scores

    subroutine compare_score3(name, score, dshape, dloc, dscale, tol, failures)
        character(len=*), intent(in) :: name !! family name used in check labels
        real(dp), intent(in) :: score(3) !! analytic score vector
        real(dp), intent(in) :: dshape !! finite-difference shape derivative
        real(dp), intent(in) :: dloc !! finite-difference location derivative
        real(dp), intent(in) :: dscale !! finite-difference scale derivative
        real(dp), intent(in) :: tol !! comparison tolerance
        integer, intent(inout) :: failures !! running count of failed checks
        call check_close(name // ' score shape', score(1), dshape, tol, tol, failures)
        call check_close(name // ' score loc', score(2), dloc, tol, tol, failures)
        call check_close(name // ' score scale', score(3), dscale, tol, tol, failures)
    end subroutine compare_score3

    subroutine compare_score4(name, score, d1, d2, dloc, dscale, tol, failures)
        character(len=*), intent(in) :: name !! family name used in check labels
        real(dp), intent(in) :: score(4) !! analytic score vector
        real(dp), intent(in) :: d1 !! first shape finite-difference derivative
        real(dp), intent(in) :: d2 !! second shape finite-difference derivative
        real(dp), intent(in) :: dloc !! finite-difference location derivative
        real(dp), intent(in) :: dscale !! finite-difference scale derivative
        real(dp), intent(in) :: tol !! comparison tolerance
        integer, intent(inout) :: failures !! running count of failed checks
        call check_close(name // ' score shape1', score(1), d1, tol, tol, failures)
        call check_close(name // ' score shape2', score(2), d2, tol, tol, failures)
        call check_close(name // ' score loc', score(3), dloc, tol, tol, failures)
        call check_close(name // ' score scale', score(4), dscale, tol, tol, failures)
    end subroutine compare_score4

    subroutine test_random_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: x1(5) = [0.3_dp, 0.6_dp, 1.0_dp, 1.7_dp, 2.8_dp]
        real(dp), parameter :: x2(5) = [0.25_dp, 0.5_dp, 1.0_dp, 2.0_dp, 4.0_dp]
        real(dp) :: actual, expected
        type(fit_result) :: result
        type(rng_state) :: reference, state

        call rng_seed(state, 97531)
        call rng_seed(reference, 97531)
        actual = gompertz_rvs(state, 1.6_dp, -0.2_dp, 1.3_dp)
        expected = gompertz_ppf(rng_uniform(reference), 1.6_dp, -0.2_dp, 1.3_dp)
        call check_close('rvs gompertz', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = invweibull_rvs(state, 1.6_dp, -0.2_dp, 1.3_dp)
        expected = invweibull_ppf(rng_uniform(reference), 1.6_dp, -0.2_dp, 1.3_dp)
        call check_close('rvs invweibull', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = betaprime_rvs(state, 1.6_dp, 2.2_dp, -0.2_dp, 1.3_dp)
        expected = betaprime_ppf(rng_uniform(reference), 1.6_dp, 2.2_dp, -0.2_dp, 1.3_dp)
        call check_close('rvs betaprime', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = burr12_rvs(state, 1.6_dp, 2.2_dp, -0.2_dp, 1.3_dp)
        expected = burr12_ppf(rng_uniform(reference), 1.6_dp, 2.2_dp, -0.2_dp, 1.3_dp)
        call check_close('rvs burr12', actual, expected, 0.0_dp, 0.0_dp, failures)

        call gompertz_fit(x1, [1.7_dp, 0.1_dp, 1.2_dp], [1.7_dp, 0.1_dp, 1.2_dp], result)
        call check_true('fit gompertz fixed success', result%success, failures)
        call invweibull_fit(x2, [1.7_dp, 0.1_dp, 1.2_dp], [1.7_dp, 0.1_dp, 1.2_dp], result)
        call check_true('fit invweibull fixed success', result%success, failures)
        call betaprime_fit(x2, [1.7_dp, 2.2_dp, 0.1_dp, 1.2_dp], &
            [1.7_dp, 2.2_dp, 0.1_dp, 1.2_dp], result)
        call check_true('fit betaprime fixed success', result%success, failures)
        call burr12_fit(x2, [1.7_dp, 2.2_dp, 0.1_dp, 1.2_dp], &
            [1.7_dp, 2.2_dp, 0.1_dp, 1.2_dp], result)
        call check_true('fit burr12 fixed success', result%success, failures)

        call check_close('C ABI gompertz cdf', &
            scifort_gompertz_cdf_f64(0.8_dp, 1.7_dp, -0.2_dp, 1.3_dp), &
            gompertz_cdf(0.8_dp, 1.7_dp, -0.2_dp, 1.3_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI invweibull ppf', &
            scifort_invweibull_ppf_f64(0.4_dp, 1.7_dp, -0.2_dp, 1.3_dp), &
            invweibull_ppf(0.4_dp, 1.7_dp, -0.2_dp, 1.3_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI betaprime pdf', &
            scifort_betaprime_pdf_f64(0.8_dp, 1.7_dp, 2.2_dp, -0.2_dp, 1.3_dp), &
            betaprime_pdf(0.8_dp, 1.7_dp, 2.2_dp, -0.2_dp, 1.3_dp), &
            0.0_dp, 0.0_dp, failures)
        call check_close('C ABI burr12 cdf', &
            scifort_burr12_cdf_f64(0.8_dp, 1.7_dp, 2.2_dp, -0.2_dp, 1.3_dp), &
            burr12_cdf(0.8_dp, 1.7_dp, 2.2_dp, -0.2_dp, 1.3_dp), 0.0_dp, 0.0_dp, failures)
    end subroutine test_random_fit_c_api

    pure function centered(plus, minus, h) result(value)
        real(dp), intent(in) :: plus !! objective at theta+h
        real(dp), intent(in) :: minus !! objective at theta-h
        real(dp), intent(in) :: h !! positive finite-difference step
        real(dp) :: value
        value = (plus - minus) / (2.0_dp * h)
    end function centered

end program test_positive_shape_continuous
