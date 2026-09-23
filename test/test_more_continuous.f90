! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Compact regression coverage for the Milestone 7 arcsine and half-normal
! additions and their likelihood, score, fit, RNG, and scalar C-ABI layers.

program test_more_continuous
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_c_api
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state, rng_uniform
    use scifort_stats
    use test_checks, only : check_close, check_true
    implicit none

    integer :: failures

    failures = 0
    call test_reference_values(failures)
    call test_identities_and_round_trips(failures)
    call test_likelihood_scores(failures)
    call test_random_fit_c_api(failures)

    if (failures == 0) then
        print '(a)', 'test_more_continuous: PASS'
    else
        print '(a,1x,i0)', 'test_more_continuous: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference_values(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: ax(5) = [0.01_dp, 0.2_dp, 0.5_dp, 0.8_dp, 0.99_dp]
        real(dp), parameter :: apdf(5) = [ &
            3.19913472585565417_dp, 0.795774715459476645_dp, &
            0.636619772367581382_dp, 0.795774715459476756_dp, &
            3.19913472585565284_dp]
        real(dp), parameter :: acdf(5) = [ &
            0.0637685608585198543_dp, 0.295167235300866526_dp, &
            0.500000000000000111_dp, 0.704832764699133474_dp, &
            0.936231439141480326_dp]
        real(dp), parameter :: asf(5) = [ &
            0.936231439141480104_dp, 0.704832764699133474_dp, &
            0.499999999999999889_dp, 0.295167235300866526_dp, &
            0.0637685608585196739_dp]
        real(dp), parameter :: hx(5) = [0.1_dp, 0.5_dp, 1.0_dp, 2.0_dp, 5.0_dp]
        real(dp), parameter :: hpdf(5) = [ &
            0.793905094954023616_dp, 0.704130653528598938_dp, &
            0.483941449038286731_dp, 0.107981933026376126_dp, &
            2.97343902946859536e-6_dp]
        real(dp), parameter :: hcdf(5) = [ &
            0.0796556745540579619_dp, 0.382924922548026125_dp, &
            0.682689492137085852_dp, 0.954499736103641583_dp, &
            0.999999426696856264_dp]
        real(dp), parameter :: hsf(5) = [ &
            0.920344325445942024_dp, 0.617075077451973764_dp, &
            0.317310507862914148_dp, 0.0455002638963583894_dp, &
            5.73303143758386553e-7_dp]
        integer :: i
        character(len=48) :: label

        do i = 1, size(ax)
            write(label, '(a,i0)') 'arcsine pdf ref ', i
            call check_close(trim(label), arcsine_pdf(ax(i)), apdf(i), 4.0e-15_dp, 4.0e-14_dp, failures)
            write(label, '(a,i0)') 'arcsine cdf ref ', i
            call check_close(trim(label), arcsine_cdf(ax(i)), acdf(i), 4.0e-15_dp, 4.0e-14_dp, failures)
            write(label, '(a,i0)') 'arcsine sf ref ', i
            call check_close(trim(label), arcsine_sf(ax(i)), asf(i), 4.0e-15_dp, 4.0e-14_dp, failures)
        end do

        do i = 1, size(hx)
            write(label, '(a,i0)') 'halfnorm pdf ref ', i
            call check_close(trim(label), halfnorm_pdf(hx(i)), hpdf(i), 4.0e-15_dp, 4.0e-14_dp, failures)
            write(label, '(a,i0)') 'halfnorm cdf ref ', i
            call check_close(trim(label), halfnorm_cdf(hx(i)), hcdf(i), 4.0e-15_dp, 4.0e-14_dp, failures)
            write(label, '(a,i0)') 'halfnorm sf ref ', i
            call check_close(trim(label), halfnorm_sf(hx(i)), hsf(i), 4.0e-15_dp, 4.0e-14_dp, failures)
        end do

        call check_close('halfnorm boundary pdf', halfnorm_pdf(0.0_dp), &
            0.797884560802865406_dp, 3.0e-15_dp, 3.0e-14_dp, failures)
        call check_true('arcsine boundary pdf inf', &
            .not. ieee_is_finite(arcsine_pdf(0.0_dp)), failures)
    end subroutine test_reference_values

    subroutine test_identities_and_round_trips(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: probs(5) = [1.0e-8_dp, 0.1_dp, 0.5_dp, 0.9_dp, 1.0_dp - 1.0e-8_dp]
        real(dp), parameter :: points(4) = [0.05_dp, 0.2_dp, 0.6_dp, 0.95_dp]
        integer :: i
        real(dp) :: x

        do i = 1, size(points)
            call check_close('arcsine symmetry pdf', arcsine_pdf(points(i)), &
                arcsine_pdf(1.0_dp - points(i)), 3.0e-15_dp, 3.0e-14_dp, failures)
            call check_close('arcsine symmetry tail', arcsine_sf(points(i)), &
                arcsine_cdf(1.0_dp - points(i)), 3.0e-15_dp, 3.0e-14_dp, failures)
            call check_close('halfnorm normal pdf identity', halfnorm_pdf(points(i)), &
                2.0_dp * normal_pdf(points(i)), 0.0_dp, 0.0_dp, failures)
            call check_close('halfnorm normal sf identity', halfnorm_sf(points(i)), &
                2.0_dp * normal_sf(points(i)), 3.0e-15_dp, 3.0e-14_dp, failures)
        end do

        do i = 1, size(probs)
            x = arcsine_ppf(probs(i))
            call check_close('arcsine ppf/cdf', arcsine_cdf(x), &
                probs(i), 4.0e-15_dp, 1.0e-9_dp, failures)
            x = halfnorm_isf(probs(i), -0.2_dp, 1.7_dp)
            call check_close('halfnorm isf/sf', halfnorm_sf(x, -0.2_dp, 1.7_dp), &
                probs(i), 4.0e-15_dp, 5.0e-10_dp, failures)
        end do

        call check_true('arcsine invalid scale nan', &
            ieee_is_nan(arcsine_cdf(0.5_dp, scale=0.0_dp)), failures)
        call check_true('halfnorm invalid p nan', ieee_is_nan(halfnorm_ppf(-0.1_dp)), failures)
    end subroutine test_identities_and_round_trips

    subroutine test_likelihood_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: xa(5) = [0.1_dp, 0.4_dp, 0.8_dp, 1.2_dp, 1.6_dp]
        real(dp), parameter :: xh(5) = [0.1_dp, 0.4_dp, 0.8_dp, 1.2_dp, 2.0_dp]
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp), parameter :: tol = 8.0e-6_dp
        real(dp) :: s(2)

        call check_close('arcsine loglike sum', arcsine_loglikelihood(xa, -0.2_dp, 2.0_dp), &
            sum(arcsine_logpdf(xa, -0.2_dp, 2.0_dp)), 3.0e-15_dp, 3.0e-15_dp, failures)
        call check_close('halfnorm loglike sum', halfnorm_loglikelihood(xh, -0.2_dp, 1.4_dp), &
            sum(halfnorm_logpdf(xh, -0.2_dp, 1.4_dp)), 3.0e-15_dp, 3.0e-15_dp, failures)

        s = arcsine_score(xa, -0.2_dp, 2.0_dp)
        call check_close('arcsine score loc', s(1), centered( &
            arcsine_loglikelihood(xa, -0.2_dp + h, 2.0_dp), &
            arcsine_loglikelihood(xa, -0.2_dp - h, 2.0_dp), h), tol, tol, failures)
        call check_close('arcsine score scale', s(2), centered( &
            arcsine_loglikelihood(xa, -0.2_dp, 2.0_dp + h), &
            arcsine_loglikelihood(xa, -0.2_dp, 2.0_dp - h), h), tol, tol, failures)

        s = halfnorm_score(xh, -0.2_dp, 1.4_dp)
        call check_close('halfnorm score loc', s(1), centered( &
            halfnorm_loglikelihood(xh, -0.2_dp + h, 1.4_dp), &
            halfnorm_loglikelihood(xh, -0.2_dp - h, 1.4_dp), h), tol, tol, failures)
        call check_close('halfnorm score scale', s(2), centered( &
            halfnorm_loglikelihood(xh, -0.2_dp, 1.4_dp + h), &
            halfnorm_loglikelihood(xh, -0.2_dp, 1.4_dp - h), h), tol, tol, failures)
    end subroutine test_likelihood_scores

    subroutine test_random_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: xa(5) = [0.1_dp, 0.4_dp, 0.8_dp, 1.2_dp, 1.6_dp]
        real(dp), parameter :: xh(5) = [0.1_dp, 0.4_dp, 0.8_dp, 1.2_dp, 2.0_dp]
        real(dp) :: actual
        real(dp) :: expected
        type(fit_result) :: result
        type(rng_state) :: reference
        type(rng_state) :: state

        call rng_seed(state, 97531)
        call rng_seed(reference, 97531)
        actual = arcsine_rvs(state, -0.2_dp, 1.8_dp)
        expected = arcsine_ppf(rng_uniform(reference), -0.2_dp, 1.8_dp)
        call check_close('rvs arcsine', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = halfnorm_rvs(state, 0.1_dp, 1.3_dp)
        expected = halfnorm_ppf(rng_uniform(reference), 0.1_dp, 1.3_dp)
        call check_close('rvs halfnorm', actual, expected, 0.0_dp, 0.0_dp, failures)

        call arcsine_fit(xa, [-0.2_dp, 2.0_dp], [-0.2_dp, 2.0_dp], result)
        call check_true('fit arcsine fixed success', result%success, failures)
        call halfnorm_fit(xh, [-0.2_dp, 1.4_dp], [-0.2_dp, 1.4_dp], result)
        call check_true('fit halfnorm fixed success', result%success, failures)

        call check_close('C ABI arcsine cdf', scifort_arcsine_cdf_f64(0.4_dp, -0.2_dp, 1.8_dp), &
            arcsine_cdf(0.4_dp, -0.2_dp, 1.8_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI halfnorm ppf', scifort_halfnorm_ppf_f64(0.4_dp, -0.2_dp, 1.8_dp), &
            halfnorm_ppf(0.4_dp, -0.2_dp, 1.8_dp), 0.0_dp, 0.0_dp, failures)
    end subroutine test_random_fit_c_api

    pure function centered(plus, minus, h) result(value)
        real(dp), intent(in) :: plus !! objective at theta + h
        real(dp), intent(in) :: minus !! objective at theta - h
        real(dp), intent(in) :: h !! positive centered-difference half step
        real(dp) :: value

        value = (plus - minus) / (2.0_dp * h)
    end function centered

end program test_more_continuous
