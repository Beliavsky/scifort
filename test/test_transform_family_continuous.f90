! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Full-stack regression coverage for Nakagami, power-normal, log-gamma, and Wald.
program test_transform_family_continuous
    use scifort_c_api
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state, rng_uniform
    use scifort_stats
    use test_checks, only : check_close, check_true
    use test_transform_family_reference
    implicit none

    integer :: failures

    failures = 0
    call test_reference_values(failures)
    call test_tail_identities(failures)
    call test_likelihood_scores(failures)
    call test_random_fit_c_api(failures)

    if (failures == 0) then
        print '(a)', 'test_transform_family_continuous: PASS'
    else
        print '(a,1x,i0)', 'test_transform_family_continuous: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference_values(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        integer :: i
        character(len=64) :: label

        do i = 1, 5
            write(label, '(a,i0)') 'nakagami pdf ref ', i
            call check_close(trim(label), nakagami_pdf(nk_x(i), reference_shape), &
                nk_pdf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            write(label, '(a,i0)') 'nakagami cdf ref ', i
            call check_close(trim(label), nakagami_cdf(nk_x(i), reference_shape), &
                nk_cdf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            write(label, '(a,i0)') 'nakagami sf ref ', i
            call check_close(trim(label), nakagami_sf(nk_x(i), reference_shape), &
                nk_sf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            write(label, '(a,i0)') 'nakagami logcdf ref ', i
            call check_close(trim(label), nakagami_logcdf(nk_x(i), reference_shape), &
                nk_logcdf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            write(label, '(a,i0)') 'nakagami logsf ref ', i
            call check_close(trim(label), nakagami_logsf(nk_x(i), reference_shape), &
                nk_logsf(i), 5.0e-12_dp, 5.0e-11_dp, failures)

            write(label, '(a,i0)') 'powernorm pdf ref ', i
            call check_close(trim(label), powernorm_pdf(pn_x(i), reference_shape), &
                pn_pdf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            write(label, '(a,i0)') 'powernorm cdf ref ', i
            call check_close(trim(label), powernorm_cdf(pn_x(i), reference_shape), &
                pn_cdf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            write(label, '(a,i0)') 'powernorm sf ref ', i
            call check_close(trim(label), powernorm_sf(pn_x(i), reference_shape), &
                pn_sf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            write(label, '(a,i0)') 'powernorm logcdf ref ', i
            call check_close(trim(label), powernorm_logcdf(pn_x(i), reference_shape), &
                pn_logcdf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            write(label, '(a,i0)') 'powernorm logsf ref ', i
            call check_close(trim(label), powernorm_logsf(pn_x(i), reference_shape), &
                pn_logsf(i), 5.0e-12_dp, 5.0e-11_dp, failures)

            write(label, '(a,i0)') 'loggamma pdf ref ', i
            call check_close(trim(label), loggamma_pdf(lg_x(i), reference_shape), &
                lg_pdf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            write(label, '(a,i0)') 'loggamma cdf ref ', i
            call check_close(trim(label), loggamma_cdf(lg_x(i), reference_shape), &
                lg_cdf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            write(label, '(a,i0)') 'loggamma sf ref ', i
            call check_close(trim(label), loggamma_sf(lg_x(i), reference_shape), &
                lg_sf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            write(label, '(a,i0)') 'loggamma logcdf ref ', i
            call check_close(trim(label), loggamma_logcdf(lg_x(i), reference_shape), &
                lg_logcdf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            write(label, '(a,i0)') 'loggamma logsf ref ', i
            call check_close(trim(label), loggamma_logsf(lg_x(i), reference_shape), &
                lg_logsf(i), 5.0e-12_dp, 5.0e-11_dp, failures)

            write(label, '(a,i0)') 'wald pdf ref ', i
            call check_close(trim(label), wald_pdf(wd_x(i)), wd_pdf(i), &
                5.0e-12_dp, 5.0e-11_dp, failures)
            write(label, '(a,i0)') 'wald cdf ref ', i
            call check_close(trim(label), wald_cdf(wd_x(i)), wd_cdf(i), &
                5.0e-12_dp, 5.0e-11_dp, failures)
            write(label, '(a,i0)') 'wald sf ref ', i
            call check_close(trim(label), wald_sf(wd_x(i)), wd_sf(i), &
                5.0e-12_dp, 5.0e-11_dp, failures)
            write(label, '(a,i0)') 'wald logcdf ref ', i
            call check_close(trim(label), wald_logcdf(wd_x(i)), wd_logcdf(i), &
                5.0e-12_dp, 5.0e-11_dp, failures)
            write(label, '(a,i0)') 'wald logsf ref ', i
            call check_close(trim(label), wald_logsf(wd_x(i)), wd_logsf(i), &
                5.0e-12_dp, 5.0e-11_dp, failures)
        end do

        do i = 1, 5
            call check_close('nakagami ppf ref', &
                nakagami_ppf(reference_probs(i), reference_shape), nk_ppf(i), &
                2.0e-10_dp, 2.0e-9_dp, failures)
            call check_close('nakagami isf ref', &
                nakagami_isf(reference_probs(i), reference_shape), nk_isf(i), &
                2.0e-10_dp, 2.0e-9_dp, failures)
            call check_close('powernorm ppf ref', &
                powernorm_ppf(reference_probs(i), reference_shape), pn_ppf(i), &
                2.0e-10_dp, 2.0e-9_dp, failures)
            call check_close('powernorm isf ref', &
                powernorm_isf(reference_probs(i), reference_shape), pn_isf(i), &
                2.0e-10_dp, 2.0e-9_dp, failures)
            call check_close('loggamma ppf ref', &
                loggamma_ppf(reference_probs(i), reference_shape), lg_ppf(i), &
                2.0e-10_dp, 2.0e-9_dp, failures)
            call check_close('loggamma isf ref', &
                loggamma_isf(reference_probs(i), reference_shape), lg_isf(i), &
                2.0e-10_dp, 2.0e-9_dp, failures)
            call check_close('wald ppf ref', wald_ppf(reference_probs(i)), wd_ppf(i), &
                2.0e-10_dp, 2.0e-9_dp, failures)
            call check_close('wald isf ref', wald_isf(reference_probs(i)), wd_isf(i), &
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
            x = nakagami_ppf(probs(i), 1.6_dp, -0.2_dp, 1.3_dp)
            call check_close('nakagami ppf/cdf', nakagami_cdf(x, 1.6_dp, -0.2_dp, 1.3_dp), &
                probs(i), 2.0e-9_dp, 2.0e-7_dp, failures)
            x = powernorm_isf(probs(i), 1.6_dp, -0.2_dp, 1.3_dp)
            call check_close('powernorm isf/sf', powernorm_sf(x, 1.6_dp, -0.2_dp, 1.3_dp), &
                probs(i), 2.0e-9_dp, 2.0e-7_dp, failures)
            x = loggamma_ppf(probs(i), 1.6_dp, -0.2_dp, 1.3_dp)
            call check_close('loggamma ppf/cdf', loggamma_cdf(x, 1.6_dp, -0.2_dp, 1.3_dp), &
                probs(i), 2.0e-9_dp, 2.0e-7_dp, failures)
            x = wald_isf(probs(i), -0.2_dp, 1.3_dp)
            call check_close('wald isf/sf', wald_sf(x, -0.2_dp, 1.3_dp), &
                probs(i), 2.0e-9_dp, 2.0e-7_dp, failures)
        end do
        call check_close('wald/invgauss identity', wald_cdf(1.4_dp, -0.2_dp, 1.3_dp), &
            invgauss_cdf(1.4_dp, 1.0_dp, -0.2_dp, 1.3_dp), 0.0_dp, 0.0_dp, failures)
    end subroutine test_tail_identities

    subroutine test_likelihood_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: xn(5) = [0.3_dp, 0.6_dp, 1.0_dp, 1.7_dp, 2.8_dp]
        real(dp), parameter :: xp(5) = [-2.0_dp, -0.7_dp, 0.2_dp, 1.1_dp, 2.4_dp]
        real(dp), parameter :: xl(5) = [-2.0_dp, -0.6_dp, 0.2_dp, 1.0_dp, 2.0_dp]
        real(dp), parameter :: xw(5) = [0.3_dp, 0.6_dp, 1.0_dp, 1.8_dp, 3.0_dp]
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp), parameter :: tol = 1.0e-4_dp
        real(dp) :: s2(2)
        real(dp) :: s3(3)

        s3 = nakagami_score(xn, 1.7_dp, 0.1_dp, 1.2_dp)
        call compare_score3('nakagami', s3, &
            centered(nakagami_loglikelihood(xn, 1.7_dp + h, 0.1_dp, 1.2_dp), &
                nakagami_loglikelihood(xn, 1.7_dp - h, 0.1_dp, 1.2_dp), h), &
            centered(nakagami_loglikelihood(xn, 1.7_dp, 0.1_dp + h, 1.2_dp), &
                nakagami_loglikelihood(xn, 1.7_dp, 0.1_dp - h, 1.2_dp), h), &
            centered(nakagami_loglikelihood(xn, 1.7_dp, 0.1_dp, 1.2_dp + h), &
                nakagami_loglikelihood(xn, 1.7_dp, 0.1_dp, 1.2_dp - h), h), tol, failures)

        s3 = powernorm_score(xp, 1.7_dp, 0.1_dp, 1.2_dp)
        call compare_score3('powernorm', s3, &
            centered(powernorm_loglikelihood(xp, 1.7_dp + h, 0.1_dp, 1.2_dp), &
                powernorm_loglikelihood(xp, 1.7_dp - h, 0.1_dp, 1.2_dp), h), &
            centered(powernorm_loglikelihood(xp, 1.7_dp, 0.1_dp + h, 1.2_dp), &
                powernorm_loglikelihood(xp, 1.7_dp, 0.1_dp - h, 1.2_dp), h), &
            centered(powernorm_loglikelihood(xp, 1.7_dp, 0.1_dp, 1.2_dp + h), &
                powernorm_loglikelihood(xp, 1.7_dp, 0.1_dp, 1.2_dp - h), h), tol, failures)

        s3 = loggamma_score(xl, 1.7_dp, 0.1_dp, 1.2_dp)
        call compare_score3('loggamma', s3, &
            centered(loggamma_loglikelihood(xl, 1.7_dp + h, 0.1_dp, 1.2_dp), &
                loggamma_loglikelihood(xl, 1.7_dp - h, 0.1_dp, 1.2_dp), h), &
            centered(loggamma_loglikelihood(xl, 1.7_dp, 0.1_dp + h, 1.2_dp), &
                loggamma_loglikelihood(xl, 1.7_dp, 0.1_dp - h, 1.2_dp), h), &
            centered(loggamma_loglikelihood(xl, 1.7_dp, 0.1_dp, 1.2_dp + h), &
                loggamma_loglikelihood(xl, 1.7_dp, 0.1_dp, 1.2_dp - h), h), tol, failures)

        s2 = wald_score(xw, 0.1_dp, 1.2_dp)
        call check_close('wald score loc', s2(1), &
            centered(wald_loglikelihood(xw, 0.1_dp + h, 1.2_dp), &
                wald_loglikelihood(xw, 0.1_dp - h, 1.2_dp), h), tol, tol, failures)
        call check_close('wald score scale', s2(2), &
            centered(wald_loglikelihood(xw, 0.1_dp, 1.2_dp + h), &
                wald_loglikelihood(xw, 0.1_dp, 1.2_dp - h), h), tol, tol, failures)
    end subroutine test_likelihood_scores

    subroutine compare_score3(name, score, dshape, dloc, dscale, tol, failures)
        character(len=*), intent(in) :: name !! family name used in check labels
        real(dp), intent(in) :: score(3) !! analytic score vector
        real(dp), intent(in) :: dshape !! finite-difference shape derivative
        real(dp), intent(in) :: dloc !! finite-difference location derivative
        real(dp), intent(in) :: dscale !! finite-difference scale derivative
        real(dp), intent(in) :: tol !! absolute and relative comparison tolerance
        integer, intent(inout) :: failures !! running count of failed checks
        call check_close(name // ' score shape', score(1), dshape, tol, tol, failures)
        call check_close(name // ' score loc', score(2), dloc, tol, tol, failures)
        call check_close(name // ' score scale', score(3), dscale, tol, tol, failures)
    end subroutine compare_score3

    subroutine test_random_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: xn(5) = [0.3_dp, 0.6_dp, 1.0_dp, 1.7_dp, 2.8_dp]
        real(dp), parameter :: xp(5) = [-2.0_dp, -0.7_dp, 0.2_dp, 1.1_dp, 2.4_dp]
        real(dp), parameter :: xl(5) = [-2.0_dp, -0.6_dp, 0.2_dp, 1.0_dp, 2.0_dp]
        real(dp), parameter :: xw(5) = [0.3_dp, 0.6_dp, 1.0_dp, 1.8_dp, 3.0_dp]
        real(dp) :: actual
        real(dp) :: expected
        type(fit_result) :: result
        type(rng_state) :: reference
        type(rng_state) :: state

        call rng_seed(state, 24680)
        call rng_seed(reference, 24680)
        actual = nakagami_rvs(state, 1.6_dp, -0.2_dp, 1.3_dp)
        expected = nakagami_ppf(rng_uniform(reference), 1.6_dp, -0.2_dp, 1.3_dp)
        call check_close('rvs nakagami', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = powernorm_rvs(state, 1.6_dp, -0.2_dp, 1.3_dp)
        expected = powernorm_ppf(rng_uniform(reference), 1.6_dp, -0.2_dp, 1.3_dp)
        call check_close('rvs powernorm', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = loggamma_rvs(state, 1.6_dp, -0.2_dp, 1.3_dp)
        expected = loggamma_ppf(rng_uniform(reference), 1.6_dp, -0.2_dp, 1.3_dp)
        call check_close('rvs loggamma', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = wald_rvs(state, -0.2_dp, 1.3_dp)
        expected = wald_ppf(rng_uniform(reference), -0.2_dp, 1.3_dp)
        call check_close('rvs wald', actual, expected, 0.0_dp, 0.0_dp, failures)

        call nakagami_fit(xn, [1.7_dp, 0.1_dp, 1.2_dp], &
            [1.7_dp, 0.1_dp, 1.2_dp], result)
        call check_true('fit nakagami fixed success', result%success, failures)
        call powernorm_fit(xp, [1.7_dp, 0.1_dp, 1.2_dp], &
            [1.7_dp, 0.1_dp, 1.2_dp], result)
        call check_true('fit powernorm fixed success', result%success, failures)
        call loggamma_fit(xl, [1.7_dp, 0.1_dp, 1.2_dp], &
            [1.7_dp, 0.1_dp, 1.2_dp], result)
        call check_true('fit loggamma fixed success', result%success, failures)
        call wald_fit(xw, [0.1_dp, 1.2_dp], [0.1_dp, 1.2_dp], result)
        call check_true('fit wald fixed success', result%success, failures)

        call check_close('C ABI nakagami cdf', &
            scifort_nakagami_cdf_f64(0.8_dp, 1.7_dp, -0.2_dp, 1.3_dp), &
            nakagami_cdf(0.8_dp, 1.7_dp, -0.2_dp, 1.3_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI powernorm ppf', &
            scifort_powernorm_ppf_f64(0.4_dp, 1.7_dp, -0.2_dp, 1.3_dp), &
            powernorm_ppf(0.4_dp, 1.7_dp, -0.2_dp, 1.3_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI loggamma pdf', &
            scifort_loggamma_pdf_f64(0.3_dp, 1.7_dp, -0.2_dp, 1.3_dp), &
            loggamma_pdf(0.3_dp, 1.7_dp, -0.2_dp, 1.3_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI wald cdf', scifort_wald_cdf_f64(0.8_dp, -0.2_dp, 1.3_dp), &
            wald_cdf(0.8_dp, -0.2_dp, 1.3_dp), 0.0_dp, 0.0_dp, failures)
    end subroutine test_random_fit_c_api

    pure function centered(plus, minus, h) result(value)
        real(dp), intent(in) :: plus !! objective at theta+h
        real(dp), intent(in) :: minus !! objective at theta-h
        real(dp), intent(in) :: h !! positive finite-difference step
        real(dp) :: value
        value = (plus - minus) / (2.0_dp * h)
    end function centered

end program test_transform_family_continuous
