! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Regression tests for the Milestone 7 continuous distributions and their
! likelihood, score, fitting, random-variate, and C-ABI integration.

program test_extended_continuous
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_c_api
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state, rng_uniform
    use scifort_stats
    use test_checks, only : check_close, check_true
    use test_extended_continuous_reference
    implicit none

    integer :: failures

    failures = 0
    call test_reference_values(failures)
    call test_round_trips(failures)
    call test_identities_and_endpoints(failures)
    call test_likelihood_and_scores(failures)
    call test_random_variates(failures)
    call test_fit_dispatch(failures)
    call test_c_api(failures)

    if (failures == 0) then
        print '(a)', 'test_extended_continuous: PASS'
    else
        print '(a,1x,i0)', 'test_extended_continuous: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference_values(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        integer :: i
        character(len=64) :: label
        real(dp), parameter :: atol = 8.0e-15_dp
        real(dp), parameter :: rtol = 8.0e-14_dp

        do i = 1, size(gr_x)
            write(label, '(a,i0)') 'gumbel_r pdf ref ', i
            call check_close(trim(label), gumbel_r_pdf(gr_x(i)), gr_pdf(i), atol, rtol, failures)
            write(label, '(a,i0)') 'gumbel_r cdf ref ', i
            call check_close(trim(label), gumbel_r_cdf(gr_x(i)), gr_cdf(i), atol, rtol, failures)
            write(label, '(a,i0)') 'gumbel_r sf ref ', i
            call check_close(trim(label), gumbel_r_sf(gr_x(i)), gr_sf(i), atol, rtol, failures)
            write(label, '(a,i0)') 'gumbel_r logcdf ref ', i
            call check_close(trim(label), gumbel_r_logcdf(gr_x(i)), gr_logcdf(i), atol, rtol, failures)
            write(label, '(a,i0)') 'gumbel_r logsf ref ', i
            call check_close(trim(label), gumbel_r_logsf(gr_x(i)), gr_logsf(i), atol, rtol, failures)
        end do

        do i = 1, size(gl_x)
            write(label, '(a,i0)') 'gumbel_l pdf ref ', i
            call check_close(trim(label), gumbel_l_pdf(gl_x(i)), gl_pdf(i), atol, rtol, failures)
            write(label, '(a,i0)') 'gumbel_l cdf ref ', i
            call check_close(trim(label), gumbel_l_cdf(gl_x(i)), gl_cdf(i), atol, rtol, failures)
            write(label, '(a,i0)') 'gumbel_l sf ref ', i
            call check_close(trim(label), gumbel_l_sf(gl_x(i)), gl_sf(i), atol, rtol, failures)
            write(label, '(a,i0)') 'gumbel_l logcdf ref ', i
            call check_close(trim(label), gumbel_l_logcdf(gl_x(i)), gl_logcdf(i), atol, rtol, failures)
            write(label, '(a,i0)') 'gumbel_l logsf ref ', i
            call check_close(trim(label), gumbel_l_logsf(gl_x(i)), gl_logsf(i), atol, rtol, failures)
        end do

        do i = 1, size(pw_x)
            write(label, '(a,i0)') 'powerlaw pdf ref ', i
            call check_close(trim(label), powerlaw_pdf(pw_x(i), pw_shape(i)), pw_pdf(i), atol, rtol, failures)
            write(label, '(a,i0)') 'powerlaw cdf ref ', i
            call check_close(trim(label), powerlaw_cdf(pw_x(i), pw_shape(i)), pw_cdf(i), atol, rtol, failures)
            write(label, '(a,i0)') 'powerlaw sf ref ', i
            call check_close(trim(label), powerlaw_sf(pw_x(i), pw_shape(i)), pw_sf(i), atol, rtol, failures)
            write(label, '(a,i0)') 'powerlaw logcdf ref ', i
            call check_close(trim(label), powerlaw_logcdf(pw_x(i), pw_shape(i)), pw_logcdf(i), atol, rtol, failures)
            write(label, '(a,i0)') 'powerlaw logsf ref ', i
            call check_close(trim(label), powerlaw_logsf(pw_x(i), pw_shape(i)), pw_logsf(i), atol, rtol, failures)
        end do

        do i = 1, size(tr_x)
            write(label, '(a,i0)') 'triang pdf ref ', i
            call check_close(trim(label), triang_pdf(tr_x(i), tr_shape(i)), tr_pdf(i), atol, rtol, failures)
            write(label, '(a,i0)') 'triang cdf ref ', i
            call check_close(trim(label), triang_cdf(tr_x(i), tr_shape(i)), tr_cdf(i), atol, rtol, failures)
            write(label, '(a,i0)') 'triang sf ref ', i
            call check_close(trim(label), triang_sf(tr_x(i), tr_shape(i)), tr_sf(i), atol, rtol, failures)
            write(label, '(a,i0)') 'triang logcdf ref ', i
            call check_close(trim(label), triang_logcdf(tr_x(i), tr_shape(i)), tr_logcdf(i), atol, rtol, failures)
            write(label, '(a,i0)') 'triang logsf ref ', i
            call check_close(trim(label), triang_logsf(tr_x(i), tr_shape(i)), tr_logsf(i), atol, rtol, failures)
        end do

        do i = 1, size(gp_x)
            write(label, '(a,i0)') 'genpareto pdf ref ', i
            call check_close(trim(label), genpareto_pdf(gp_x(i), gp_shape(i)), gp_pdf(i), atol, rtol, failures)
            write(label, '(a,i0)') 'genpareto cdf ref ', i
            call check_close(trim(label), genpareto_cdf(gp_x(i), gp_shape(i)), gp_cdf(i), atol, rtol, failures)
            write(label, '(a,i0)') 'genpareto sf ref ', i
            call check_close(trim(label), genpareto_sf(gp_x(i), gp_shape(i)), gp_sf(i), atol, rtol, failures)
            write(label, '(a,i0)') 'genpareto logcdf ref ', i
            call check_close(trim(label), genpareto_logcdf(gp_x(i), gp_shape(i)), gp_logcdf(i), atol, rtol, failures)
            write(label, '(a,i0)') 'genpareto logsf ref ', i
            call check_close(trim(label), genpareto_logsf(gp_x(i), gp_shape(i)), gp_logsf(i), atol, rtol, failures)
        end do
    end subroutine test_reference_values

    subroutine test_round_trips(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        integer :: i
        real(dp), parameter :: p(6) = [1.0e-10_dp, 1.0e-5_dp, 0.1_dp, 0.5_dp, 0.9_dp, 1.0_dp - 1.0e-10_dp]
        real(dp) :: x

        do i = 1, size(p)
            x = gumbel_r_ppf(p(i), -0.3_dp, 1.7_dp)
            call check_close('gumbel_r ppf/cdf', gumbel_r_cdf(x, -0.3_dp, 1.7_dp), p(i), &
                2.0e-15_dp, 3.0e-12_dp, failures)
            x = gumbel_l_isf(p(i), 0.4_dp, 1.2_dp)
            call check_close('gumbel_l isf/sf', gumbel_l_sf(x, 0.4_dp, 1.2_dp), p(i), &
                2.0e-15_dp, 3.0e-12_dp, failures)
            x = powerlaw_ppf(p(i), 2.3_dp, -0.2_dp, 2.0_dp)
            call check_close('powerlaw ppf/cdf', powerlaw_cdf(x, 2.3_dp, -0.2_dp, 2.0_dp), p(i), &
                2.0e-15_dp, 3.0e-12_dp, failures)
            x = triang_ppf(p(i), 0.35_dp, -0.2_dp, 2.0_dp)
            call check_close('triang ppf/cdf', triang_cdf(x, 0.35_dp, -0.2_dp, 2.0_dp), p(i), &
                2.0e-15_dp, 3.0e-12_dp, failures)
            x = genpareto_isf(p(i), 0.4_dp, -0.2_dp, 1.6_dp)
            call check_close('genpareto isf/sf', genpareto_sf(x, 0.4_dp, -0.2_dp, 1.6_dp), p(i), &
                2.0e-15_dp, 3.0e-12_dp, failures)
        end do
    end subroutine test_round_trips

    subroutine test_identities_and_endpoints(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: x(5) = [-2.0_dp, -0.5_dp, 0.0_dp, 0.7_dp, 2.5_dp]
        integer :: i

        do i = 1, size(x)
            call check_close('gumbel reflection pdf', gumbel_l_pdf(x(i)), gumbel_r_pdf(-x(i)), &
                2.0e-15_dp, 2.0e-14_dp, failures)
            call check_close('gumbel reflection cdf', gumbel_l_cdf(x(i)), gumbel_r_sf(-x(i)), &
                2.0e-15_dp, 2.0e-14_dp, failures)
            call check_close('genpareto c=0 exponential', genpareto_cdf(max(x(i), 0.0_dp), 0.0_dp), &
                exponential_cdf(max(x(i), 0.0_dp)), 2.0e-15_dp, 2.0e-14_dp, failures)
        end do

        call check_close('powerlaw a=1 uniform cdf', powerlaw_cdf(0.37_dp, 1.0_dp), &
            uniform_cdf(0.37_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('genpareto c=-1 uniform cdf', genpareto_cdf(0.37_dp, -1.0_dp), &
            uniform_cdf(0.37_dp), 2.0e-15_dp, 2.0e-14_dp, failures)
        call check_close('triang c=0 pdf left', triang_pdf(0.0_dp, 0.0_dp), 2.0_dp, &
            0.0_dp, 0.0_dp, failures)
        call check_close('triang c=1 pdf right', triang_pdf(1.0_dp, 1.0_dp), 2.0_dp, &
            0.0_dp, 0.0_dp, failures)
        call check_close('genpareto finite upper ppf', genpareto_ppf(1.0_dp, -0.5_dp), 2.0_dp, &
            0.0_dp, 0.0_dp, failures)
        call check_true('genpareto positive shape ppf 1 inf', &
            .not. ieee_is_finite(genpareto_ppf(1.0_dp, 0.5_dp)), failures)
        call check_true('powerlaw invalid shape nan', ieee_is_nan(powerlaw_pdf(0.5_dp, 0.0_dp)), failures)
        call check_true('triang invalid shape nan', ieee_is_nan(triang_cdf(0.5_dp, 1.1_dp)), failures)
        call check_true('gumbel invalid scale nan', ieee_is_nan(gumbel_r_pdf(0.0_dp, scale=0.0_dp)), failures)
        call check_true('genpareto invalid p nan', ieee_is_nan(genpareto_ppf(-0.1_dp, 0.2_dp)), failures)
    end subroutine test_identities_and_endpoints

    subroutine test_likelihood_and_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: xg(5) = [-1.2_dp, -0.3_dp, 0.4_dp, 1.1_dp, 2.0_dp]
        real(dp), parameter :: xp(5) = [0.1_dp, 0.4_dp, 0.8_dp, 1.2_dp, 1.6_dp]
        real(dp), parameter :: xt(5) = [0.1_dp, 0.35_dp, 0.75_dp, 1.25_dp, 1.65_dp]
        real(dp), parameter :: xgp(5) = [0.2_dp, 0.5_dp, 1.0_dp, 1.8_dp, 3.0_dp]
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp), parameter :: tol = 7.0e-6_dp
        real(dp) :: s2(2)
        real(dp) :: s3(3)

        call check_close('gumbel_r loglike sum', gumbel_r_loglikelihood(xg, 0.1_dp, 1.3_dp), &
            sum(gumbel_r_logpdf(xg, 0.1_dp, 1.3_dp)), 3.0e-15_dp, 3.0e-15_dp, failures)
        call check_close('powerlaw loglike sum', powerlaw_loglikelihood(xp, 2.2_dp, -0.2_dp, 2.0_dp), &
            sum(powerlaw_logpdf(xp, 2.2_dp, -0.2_dp, 2.0_dp)), 3.0e-15_dp, 3.0e-15_dp, failures)
        call check_close('triang loglike sum', triang_loglikelihood(xt, 0.4_dp, -0.2_dp, 2.0_dp), &
            sum(triang_logpdf(xt, 0.4_dp, -0.2_dp, 2.0_dp)), 3.0e-15_dp, 3.0e-15_dp, failures)
        call check_close('genpareto loglike sum', genpareto_loglikelihood(xgp, 0.2_dp, -0.1_dp, 1.5_dp), &
            sum(genpareto_logpdf(xgp, 0.2_dp, -0.1_dp, 1.5_dp)), 3.0e-15_dp, 3.0e-15_dp, failures)

        s2 = gumbel_r_score(xg, 0.1_dp, 1.3_dp)
        call check_close('gumbel_r score loc', s2(1), centered( &
            gumbel_r_loglikelihood(xg, 0.1_dp + h, 1.3_dp), &
            gumbel_r_loglikelihood(xg, 0.1_dp - h, 1.3_dp), h), tol, tol, failures)
        call check_close('gumbel_r score scale', s2(2), centered( &
            gumbel_r_loglikelihood(xg, 0.1_dp, 1.3_dp + h), &
            gumbel_r_loglikelihood(xg, 0.1_dp, 1.3_dp - h), h), tol, tol, failures)

        s2 = gumbel_l_score(xg, 0.1_dp, 1.3_dp)
        call check_close('gumbel_l score loc', s2(1), centered( &
            gumbel_l_loglikelihood(xg, 0.1_dp + h, 1.3_dp), &
            gumbel_l_loglikelihood(xg, 0.1_dp - h, 1.3_dp), h), tol, tol, failures)
        call check_close('gumbel_l score scale', s2(2), centered( &
            gumbel_l_loglikelihood(xg, 0.1_dp, 1.3_dp + h), &
            gumbel_l_loglikelihood(xg, 0.1_dp, 1.3_dp - h), h), tol, tol, failures)

        s3 = powerlaw_score(xp, 2.2_dp, -0.2_dp, 2.0_dp)
        call check_close('powerlaw score a', s3(1), centered( &
            powerlaw_loglikelihood(xp, 2.2_dp + h, -0.2_dp, 2.0_dp), &
            powerlaw_loglikelihood(xp, 2.2_dp - h, -0.2_dp, 2.0_dp), h), tol, tol, failures)
        call check_close('powerlaw score loc', s3(2), centered( &
            powerlaw_loglikelihood(xp, 2.2_dp, -0.2_dp + h, 2.0_dp), &
            powerlaw_loglikelihood(xp, 2.2_dp, -0.2_dp - h, 2.0_dp), h), tol, tol, failures)
        call check_close('powerlaw score scale', s3(3), centered( &
            powerlaw_loglikelihood(xp, 2.2_dp, -0.2_dp, 2.0_dp + h), &
            powerlaw_loglikelihood(xp, 2.2_dp, -0.2_dp, 2.0_dp - h), h), tol, tol, failures)

        s3 = triang_score(xt, 0.4_dp, -0.2_dp, 2.0_dp)
        call check_close('triang score c', s3(1), centered( &
            triang_loglikelihood(xt, 0.4_dp + h, -0.2_dp, 2.0_dp), &
            triang_loglikelihood(xt, 0.4_dp - h, -0.2_dp, 2.0_dp), h), tol, tol, failures)
        call check_close('triang score loc', s3(2), centered( &
            triang_loglikelihood(xt, 0.4_dp, -0.2_dp + h, 2.0_dp), &
            triang_loglikelihood(xt, 0.4_dp, -0.2_dp - h, 2.0_dp), h), tol, tol, failures)
        call check_close('triang score scale', s3(3), centered( &
            triang_loglikelihood(xt, 0.4_dp, -0.2_dp, 2.0_dp + h), &
            triang_loglikelihood(xt, 0.4_dp, -0.2_dp, 2.0_dp - h), h), tol, tol, failures)

        s3 = genpareto_score(xgp, 0.2_dp, -0.1_dp, 1.5_dp)
        call check_close('genpareto score c', s3(1), centered( &
            genpareto_loglikelihood(xgp, 0.2_dp + h, -0.1_dp, 1.5_dp), &
            genpareto_loglikelihood(xgp, 0.2_dp - h, -0.1_dp, 1.5_dp), h), tol, tol, failures)
        call check_close('genpareto score loc', s3(2), centered( &
            genpareto_loglikelihood(xgp, 0.2_dp, -0.1_dp + h, 1.5_dp), &
            genpareto_loglikelihood(xgp, 0.2_dp, -0.1_dp - h, 1.5_dp), h), tol, tol, failures)
        call check_close('genpareto score scale', s3(3), centered( &
            genpareto_loglikelihood(xgp, 0.2_dp, -0.1_dp, 1.5_dp + h), &
            genpareto_loglikelihood(xgp, 0.2_dp, -0.1_dp, 1.5_dp - h), h), tol, tol, failures)
    end subroutine test_likelihood_and_scores

    subroutine test_random_variates(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp) :: actual
        real(dp) :: expected
        type(rng_state) :: reference
        type(rng_state) :: state

        call rng_seed(state, 86420)
        call rng_seed(reference, 86420)
        actual = gumbel_r_rvs(state, -0.2_dp, 1.4_dp)
        expected = gumbel_r_ppf(rng_uniform(reference), -0.2_dp, 1.4_dp)
        call check_close('rvs gumbel_r', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = gumbel_l_rvs(state, 0.3_dp, 1.1_dp)
        expected = gumbel_l_ppf(rng_uniform(reference), 0.3_dp, 1.1_dp)
        call check_close('rvs gumbel_l', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = powerlaw_rvs(state, 2.1_dp, -0.2_dp, 2.0_dp)
        expected = powerlaw_ppf(rng_uniform(reference), 2.1_dp, -0.2_dp, 2.0_dp)
        call check_close('rvs powerlaw', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = triang_rvs(state, 0.35_dp, -0.2_dp, 2.0_dp)
        expected = triang_ppf(rng_uniform(reference), 0.35_dp, -0.2_dp, 2.0_dp)
        call check_close('rvs triang', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = genpareto_rvs(state, 0.25_dp, -0.2_dp, 1.3_dp)
        expected = genpareto_ppf(rng_uniform(reference), 0.25_dp, -0.2_dp, 1.3_dp)
        call check_close('rvs genpareto', actual, expected, 0.0_dp, 0.0_dp, failures)
    end subroutine test_random_variates

    subroutine test_fit_dispatch(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: xg(5) = [-1.2_dp, -0.3_dp, 0.4_dp, 1.1_dp, 2.0_dp]
        real(dp), parameter :: xp(5) = [0.1_dp, 0.4_dp, 0.8_dp, 1.2_dp, 1.6_dp]
        real(dp), parameter :: xt(5) = [0.1_dp, 0.35_dp, 0.75_dp, 1.25_dp, 1.65_dp]
        real(dp), parameter :: xgp(5) = [0.2_dp, 0.5_dp, 1.0_dp, 1.8_dp, 3.0_dp]
        type(fit_result) :: result

        call gumbel_r_fit(xg, [0.1_dp, 1.3_dp], [0.1_dp, 1.3_dp], result)
        call check_true('fit gumbel_r fixed success', result%success, failures)
        call check_true('fit gumbel_r fixed params', all(result%params == [0.1_dp, 1.3_dp]), failures)
        call gumbel_l_fit(xg, [0.1_dp, 1.3_dp], [0.1_dp, 1.3_dp], result)
        call check_true('fit gumbel_l fixed success', result%success, failures)
        call powerlaw_fit(xp, [2.2_dp, -0.2_dp, 2.0_dp], [2.2_dp, -0.2_dp, 2.0_dp], result)
        call check_true('fit powerlaw fixed success', result%success, failures)
        call triang_fit(xt, [0.4_dp, -0.2_dp, 2.0_dp], [0.4_dp, -0.2_dp, 2.0_dp], result)
        call check_true('fit triang fixed success', result%success, failures)
        call genpareto_fit(xgp, [0.2_dp, -0.1_dp, 1.5_dp], [0.2_dp, -0.1_dp, 1.5_dp], result)
        call check_true('fit genpareto fixed success', result%success, failures)
    end subroutine test_fit_dispatch

    subroutine test_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        call check_close('C ABI gumbel_r pdf', scifort_gumbel_r_pdf_f64(0.3_dp, -0.2_dp, 1.4_dp), &
            gumbel_r_pdf(0.3_dp, -0.2_dp, 1.4_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI gumbel_l cdf', scifort_gumbel_l_cdf_f64(0.3_dp, -0.2_dp, 1.4_dp), &
            gumbel_l_cdf(0.3_dp, -0.2_dp, 1.4_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI powerlaw ppf', scifort_powerlaw_ppf_f64(0.4_dp, 2.2_dp, -0.2_dp, 2.0_dp), &
            powerlaw_ppf(0.4_dp, 2.2_dp, -0.2_dp, 2.0_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI triang pdf', scifort_triang_pdf_f64(0.5_dp, 0.4_dp, -0.2_dp, 2.0_dp), &
            triang_pdf(0.5_dp, 0.4_dp, -0.2_dp, 2.0_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI genpareto ppf', scifort_genpareto_ppf_f64(0.4_dp, 0.2_dp, -0.2_dp, 1.4_dp), &
            genpareto_ppf(0.4_dp, 0.2_dp, -0.2_dp, 1.4_dp), 0.0_dp, 0.0_dp, failures)
    end subroutine test_c_api

    pure function centered(plus, minus, h) result(value)
        real(dp), intent(in) :: plus !! objective at theta + h
        real(dp), intent(in) :: minus !! objective at theta - h
        real(dp), intent(in) :: h !! positive centered-difference half step
        real(dp) :: value

        value = (plus - minus) / (2.0_dp * h)
    end function centered

end program test_extended_continuous
