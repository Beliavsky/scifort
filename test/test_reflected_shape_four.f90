! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_reflected_shape_four
    use scifort_c_api
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state, rng_uniform
    use scifort_stats
    use test_checks, only : check_close, check_true
    implicit none

    integer :: failures
    real(dp), parameter :: probs(5) = [1.0e-6_dp, 0.1_dp, 0.5_dp, 0.9_dp, 1.0_dp - 1.0e-6_dp]

    failures = 0
    call test_reference(failures)
    call test_scores(failures)
    call test_rng_fit_c_api(failures)

    if (failures == 0) then
        print '(a)', 'test_reflected_shape_four: PASS'
    else
        print '(a,1x,i0)', 'test_reflected_shape_four: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: ll_x(5) = [-10.0_dp, -2.0_dp, -1.0_dp, -0.4_dp, -0.1_dp]
        real(dp), parameter :: ll_pdf(5) = [0.01200038948430136_dp, 0.1098478223669306_dp, &
            0.24197072451914337_dp, 0.45180598167045316_dp, 0.08500366602520341_dp]
        real(dp), parameter :: ll_cdf(5) = [0.24817036595415076_dp, 0.5204998778130465_dp, &
            0.6826894921370859_dp, 0.8861537019933419_dp, 0.9984345977419975_dp]
        real(dp), parameter :: ll_ppf(5) = [-6.3661977247199316e11_dp, -63.328117677016621_dp, &
            -2.1981093383177326_dp, -0.3696115094681951_dp, -0.041791821021978919_dp]
        real(dp), parameter :: wm_x(5) = [-2.0_dp, -1.0_dp, -0.7_dp, -0.2_dp, -0.05_dp]
        real(dp), parameter :: wm_pdf(5) = [0.10718721420110489_dp, 0.6253950499914519_dp, &
            0.7676791945248881_dp, 0.5164356682052841_dp, 0.20752043119749888_dp]
        real(dp), parameter :: wm_cdf(5) = [0.0388126293957927_dp, 0.36787944117144233_dp, &
            0.5796451697945955_dp, 0.9372302799973403_dp, 0.9938776781289106_dp]
        real(dp), parameter :: wm_ppf(5) = [-4.6860193954781675_dp, -1.6333078605318136_dp, &
            -0.80606101718592083_dp, -0.26613662169674629_dp, -0.00029552101044325507_dp]
        real(dp), parameter :: rd_x(5) = [-0.9_dp, -0.4_dp, 0.0_dp, 0.4_dp, 0.9_dp]
        real(dp), parameter :: rd_pdf(5) = [0.33477260088407235_dp, 0.5632275580759529_dp, &
            0.5986681395059306_dp, 0.5632275580759529_dp, 0.3347726008840727_dp]
        real(dp), parameter :: rd_cdf(5) = [0.02499040897068806_dp, 0.26515183176298585_dp, &
            0.5_dp, 0.7348481682370143_dp, 0.9750095910293122_dp]
        real(dp), parameter :: rd_ppf(5) = [-0.9999451592548781_dp, -0.7165807860967973_dp, &
            0.0_dp, 0.7165807860967974_dp, 0.9999451592548769_dp]
        real(dp), parameter :: sc_x(5) = [-3.0_dp, -0.5_dp, 0.0_dp, 0.7_dp, 4.0_dp]
        real(dp), parameter :: sc_pdf(5) = [0.01427284976520579_dp, 0.1999790734760618_dp, &
            0.3183098861837907_dp, 0.2508626021924145_dp, 0.03254985370009587_dp]
        real(dp), parameter :: sc_cdf(5) = [0.04414627632904577_dp, 0.18933563989590071_dp, &
            0.325_dp, 0.530556815783638_dp, 0.8601284593446146_dp]
        real(dp), parameter :: sc_ppf(5) = [-134485.92690901543_dp, -1.2384715273469951_dp, &
            0.58233315680943731_dp, 5.6960978923685897_dp, 580119.76745552348_dp]
        integer :: i

        do i = 1, 5
            call check_close('levy_l pdf', levy_l_pdf(ll_x(i)), ll_pdf(i), &
                8.0e-13_dp, 8.0e-12_dp, failures)
            call check_close('levy_l cdf', levy_l_cdf(ll_x(i)), ll_cdf(i), &
                8.0e-13_dp, 8.0e-12_dp, failures)
            call check_close('levy_l ppf', levy_l_ppf(probs(i)), ll_ppf(i), &
                3.0e-6_dp, 5.0e-10_dp, failures)
            call check_close('weibull_max pdf', weibull_max_pdf(wm_x(i), 1.7_dp), wm_pdf(i), &
                8.0e-13_dp, 8.0e-12_dp, failures)
            call check_close('weibull_max cdf', weibull_max_cdf(wm_x(i), 1.7_dp), wm_cdf(i), &
                8.0e-13_dp, 8.0e-12_dp, failures)
            call check_close('weibull_max ppf', weibull_max_ppf(probs(i), 1.7_dp), wm_ppf(i), &
                2.0e-11_dp, 2.0e-10_dp, failures)
            call check_close('rdist pdf', rdist_pdf(rd_x(i), 2.7_dp), rd_pdf(i), &
                2.0e-12_dp, 2.0e-11_dp, failures)
            call check_close('rdist cdf', rdist_cdf(rd_x(i), 2.7_dp), rd_cdf(i), &
                2.0e-12_dp, 2.0e-11_dp, failures)
            call check_close('rdist ppf', rdist_ppf(probs(i), 2.7_dp), rd_ppf(i), &
                3.0e-10_dp, 3.0e-9_dp, failures)
            call check_close('skewcauchy pdf', skewcauchy_pdf(sc_x(i), 0.35_dp), sc_pdf(i), &
                2.0e-12_dp, 2.0e-11_dp, failures)
            call check_close('skewcauchy cdf', skewcauchy_cdf(sc_x(i), 0.35_dp), sc_cdf(i), &
                2.0e-12_dp, 2.0e-11_dp, failures)
            call check_close('skewcauchy ppf', skewcauchy_ppf(probs(i), 0.35_dp), sc_ppf(i), &
                4.0e-6_dp, 4.0e-10_dp, failures)
        end do
        call check_close('levy_l logsf tail', levy_l_logsf(-0.1_dp), &
            -6.459612454150123_dp, 3.0e-12_dp, 3.0e-11_dp, failures)
        call check_close('weibull_max logcdf', weibull_max_logcdf(-2.0_dp, 1.7_dp), &
            -3.249009585424942_dp, 3.0e-12_dp, 3.0e-11_dp, failures)
        call check_close('rdist logsf', rdist_logsf(0.9_dp, 2.7_dp), &
            -3.689263168895514_dp, 3.0e-12_dp, 3.0e-11_dp, failures)
        call check_close('skewcauchy logcdf', skewcauchy_logcdf(-3.0_dp, 0.35_dp), &
            -3.120246696847067_dp, 3.0e-12_dp, 3.0e-11_dp, failures)
    end subroutine test_reference

    subroutine test_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp), parameter :: tol = 5.0e-5_dp
        real(dp), parameter :: left_data(4) = [-3.0_dp, -1.5_dp, -0.8_dp, -0.3_dp]
        real(dp), parameter :: rd_data(4) = [-0.6_dp, -0.2_dp, 0.2_dp, 0.6_dp]
        real(dp), parameter :: sc_data(4) = [-1.2_dp, -0.3_dp, 0.4_dp, 1.5_dp]
        real(dp) :: q2(2), q3(3)

        q2 = levy_l_score(left_data, 0.2_dp, 1.1_dp)
        call check_close('levy_l score loc', q2(1), centered( &
            levy_l_loglikelihood(left_data, 0.2_dp+h, 1.1_dp), &
            levy_l_loglikelihood(left_data, 0.2_dp-h, 1.1_dp), h), tol, tol, failures)
        call check_close('levy_l score scale', q2(2), centered( &
            levy_l_loglikelihood(left_data, 0.2_dp, 1.1_dp+h), &
            levy_l_loglikelihood(left_data, 0.2_dp, 1.1_dp-h), h), tol, tol, failures)

        q3 = weibull_max_score(left_data, 1.7_dp, 0.2_dp, 1.1_dp)
        call check_close('weibull_max score c', q3(1), centered( &
            weibull_max_loglikelihood(left_data, 1.7_dp+h, 0.2_dp, 1.1_dp), &
            weibull_max_loglikelihood(left_data, 1.7_dp-h, 0.2_dp, 1.1_dp), h), tol, tol, failures)
        call check_close('weibull_max score loc', q3(2), centered( &
            weibull_max_loglikelihood(left_data, 1.7_dp, 0.2_dp+h, 1.1_dp), &
            weibull_max_loglikelihood(left_data, 1.7_dp, 0.2_dp-h, 1.1_dp), h), tol, tol, failures)
        call check_close('weibull_max score scale', q3(3), centered( &
            weibull_max_loglikelihood(left_data, 1.7_dp, 0.2_dp, 1.1_dp+h), &
            weibull_max_loglikelihood(left_data, 1.7_dp, 0.2_dp, 1.1_dp-h), h), tol, tol, failures)

        q3 = rdist_score(rd_data, 2.7_dp, 0.05_dp, 1.2_dp)
        call check_close('rdist score c', q3(1), centered( &
            rdist_loglikelihood(rd_data, 2.7_dp+h, 0.05_dp, 1.2_dp), &
            rdist_loglikelihood(rd_data, 2.7_dp-h, 0.05_dp, 1.2_dp), h), tol, tol, failures)
        call check_close('rdist score loc', q3(2), centered( &
            rdist_loglikelihood(rd_data, 2.7_dp, 0.05_dp+h, 1.2_dp), &
            rdist_loglikelihood(rd_data, 2.7_dp, 0.05_dp-h, 1.2_dp), h), tol, tol, failures)
        call check_close('rdist score scale', q3(3), centered( &
            rdist_loglikelihood(rd_data, 2.7_dp, 0.05_dp, 1.2_dp+h), &
            rdist_loglikelihood(rd_data, 2.7_dp, 0.05_dp, 1.2_dp-h), h), tol, tol, failures)

        q3 = skewcauchy_score(sc_data, 0.35_dp, 0.1_dp, 1.2_dp)
        call check_close('skewcauchy score a', q3(1), centered( &
            skewcauchy_loglikelihood(sc_data, 0.35_dp+h, 0.1_dp, 1.2_dp), &
            skewcauchy_loglikelihood(sc_data, 0.35_dp-h, 0.1_dp, 1.2_dp), h), tol, tol, failures)
        call check_close('skewcauchy score loc', q3(2), centered( &
            skewcauchy_loglikelihood(sc_data, 0.35_dp, 0.1_dp+h, 1.2_dp), &
            skewcauchy_loglikelihood(sc_data, 0.35_dp, 0.1_dp-h, 1.2_dp), h), tol, tol, failures)
        call check_close('skewcauchy score scale', q3(3), centered( &
            skewcauchy_loglikelihood(sc_data, 0.35_dp, 0.1_dp, 1.2_dp+h), &
            skewcauchy_loglikelihood(sc_data, 0.35_dp, 0.1_dp, 1.2_dp-h), h), tol, tol, failures)
    end subroutine test_scores

    subroutine test_rng_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: left_data(4) = [-3.0_dp, -1.5_dp, -0.8_dp, -0.3_dp]
        real(dp), parameter :: rd_data(4) = [-0.6_dp, -0.2_dp, 0.2_dp, 0.6_dp]
        real(dp), parameter :: sc_data(4) = [-1.2_dp, -0.3_dp, 0.4_dp, 1.5_dp]
        real(dp) :: actual, expected
        type(rng_state) :: reference, state
        type(fit_result) :: result

        call rng_seed(state, 123456)
        call rng_seed(reference, 123456)
        actual = levy_l_rvs(state, 0.2_dp, 1.1_dp)
        expected = levy_l_ppf(rng_uniform(reference), 0.2_dp, 1.1_dp)
        call check_close('levy_l rvs', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = weibull_max_rvs(state, 1.7_dp, 0.2_dp, 1.1_dp)
        expected = weibull_max_ppf(rng_uniform(reference), 1.7_dp, 0.2_dp, 1.1_dp)
        call check_close('weibull_max rvs', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = rdist_rvs(state, 2.7_dp, 0.05_dp, 1.2_dp)
        expected = rdist_ppf(rng_uniform(reference), 2.7_dp, 0.05_dp, 1.2_dp)
        call check_close('rdist rvs', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = skewcauchy_rvs(state, 0.35_dp, 0.1_dp, 1.2_dp)
        expected = skewcauchy_ppf(rng_uniform(reference), 0.35_dp, 0.1_dp, 1.2_dp)
        call check_close('skewcauchy rvs', actual, expected, 0.0_dp, 0.0_dp, failures)

        call levy_l_fit(left_data, [0.2_dp,1.1_dp], [0.2_dp,1.1_dp], result)
        call check_true('levy_l fixed fit', result%success, failures)
        call weibull_max_fit(left_data, [1.7_dp,0.2_dp,1.1_dp], &
            [1.7_dp,0.2_dp,1.1_dp], result)
        call check_true('weibull_max fixed fit', result%success, failures)
        call rdist_fit(rd_data, [2.7_dp,0.05_dp,1.2_dp], [2.7_dp,0.05_dp,1.2_dp], result)
        call check_true('rdist fixed fit', result%success, failures)
        call skewcauchy_fit(sc_data, [0.35_dp,0.1_dp,1.2_dp], &
            [0.35_dp,0.1_dp,1.2_dp], result)
        call check_true('skewcauchy fixed fit', result%success, failures)

        call check_close('C ABI levy_l cdf', scifort_levy_l_cdf_f64(-1.0_dp,0.0_dp,1.0_dp), &
            levy_l_cdf(-1.0_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI weibull_max pdf', &
            scifort_weibull_max_pdf_f64(-0.7_dp,1.7_dp,0.0_dp,1.0_dp), &
            weibull_max_pdf(-0.7_dp,1.7_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI rdist ppf', scifort_rdist_ppf_f64(0.9_dp,2.7_dp,0.0_dp,1.0_dp), &
            rdist_ppf(0.9_dp,2.7_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI skewcauchy cdf', &
            scifort_skewcauchy_cdf_f64(0.7_dp,0.35_dp,0.0_dp,1.0_dp), &
            skewcauchy_cdf(0.7_dp,0.35_dp), 0.0_dp, 0.0_dp, failures)
    end subroutine test_rng_fit_c_api

    pure function centered(plus, minus, h) result(value)
        real(dp), intent(in) :: plus !! objective at theta+h
        real(dp), intent(in) :: minus !! objective at theta-h
        real(dp), intent(in) :: h !! positive finite-difference step
        real(dp) :: value
        value = (plus - minus) / (2.0_dp * h)
    end function centered

end program test_reflected_shape_four
