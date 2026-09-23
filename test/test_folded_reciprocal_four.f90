! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_folded_reciprocal_four
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
        print '(a)', 'test_folded_reciprocal_four: PASS'
    else
        print '(a,1x,i0)', 'test_folded_reciprocal_four: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: fn_x(5) = [0.0_dp, 0.2_dp, 1.0_dp, 2.0_dp, 4.0_dp]
        real(dp), parameter :: fn_pdf(5) = [0.34273718409561471_dp, 0.34736977269844227_dp, &
            0.40971485320212525_dp, 0.31397650230581498_dp, 0.010421251949344262_dp]
        real(dp), parameter :: fn_cdf(5) = [0.0_dp, 0.068858859677524586_dp, &
            0.37136446778937149_dp, 0.75755292363454318_dp, 0.99653296829561899_dp]
        real(dp), parameter :: fn_ppf(5) = [2.9176875062708102e-6_dp, 0.28906241282139217_dp, &
            1.3113048158704315_dp, 2.5818469693040189_dp, 6.053424328333401_dp]
        real(dp), parameter :: fc_x(5) = [0.0_dp, 0.2_dp, 1.0_dp, 2.0_dp, 8.0_dp]
        real(dp), parameter :: fc_pdf(5) = [0.38818278802901302_dp, 0.39320632999174143_dp, &
            0.38114028461774646_dp, 0.16646277706369889_dp, 0.010082037336456056_dp]
        real(dp), parameter :: fc_cdf(5) = [0.0_dp, 0.077979130377369282_dp, &
            0.40141849097190596_dp, 0.66967002694267497_dp, 0.92005424224524623_dp]
        real(dp), parameter :: fc_ppf(5) = [2.5761059759134373e-6_dp, 0.25580351919673433_dp, &
            1.2806248474865698_dp, 6.411170593547233_dp, 636619.77227825345_dp]
        real(dp), parameter :: ri_x(5) = [0.05_dp, 0.2_dp, 0.7_dp, 2.0_dp, 6.0_dp]
        real(dp), parameter :: ri_pdf(5) = [9.9502890777039083e-9_dp, 0.020492650593672581_dp, &
            0.32635955361677887_dp, 0.25998150722708241_dp, 0.028543818921045842_dp]
        real(dp), parameter :: ri_cdf(5) = [2.2798886245177165e-11_dp, 0.00064673160319552022_dp, &
            0.09655607061239864_dp, 0.52341061046392479_dp, 0.94789792550048813_dp]
        real(dp), parameter :: ri_ppf(5) = [0.096055756061007871_dp, 0.71050545537020804_dp, &
            1.9118862002339738_dp, 4.8150708142323388_dp, 26.6142305007712_dp]
        real(dp), parameter :: tp_x(5) = [1.05_dp, 1.2_dp, 2.0_dp, 3.5_dp, 4.9_dp]
        real(dp), parameter :: tp_pdf(5) = [1.5934759476359115_dp, 1.1111360946921862_dp, &
            0.27975356634565091_dp, 0.061740892084143191_dp, 0.024890132234936461_dp]
        real(dp), parameter :: tp_cdf(5) = [0.085114289535269283_dp, 0.28498866093943603_dp, &
            0.74019817913903707_dp, 0.94220642054891968_dp, 0.99757787604498582_dp]
        real(dp), parameter :: tp_ppf(5) = [1.0000005501026061_dp, 1.0594553440960661_dp, &
            1.4488724079983382_dp, 2.9568082837742513_dp, 4.9999575715255586_dp]
        integer :: i

        do i = 1, 5
            call check_close('foldnorm pdf', foldnorm_pdf(fn_x(i), 1.3_dp), fn_pdf(i), &
                4.0e-12_dp, 4.0e-11_dp, failures)
            call check_close('foldnorm cdf', foldnorm_cdf(fn_x(i), 1.3_dp), fn_cdf(i), &
                4.0e-12_dp, 4.0e-11_dp, failures)
            call check_close('foldnorm ppf', foldnorm_ppf(probs(i), 1.3_dp), fn_ppf(i), &
                2.0e-9_dp, 3.0e-9_dp, failures)
            call check_close('foldcauchy pdf', foldcauchy_pdf(fc_x(i), 0.8_dp), fc_pdf(i), &
                4.0e-12_dp, 4.0e-11_dp, failures)
            call check_close('foldcauchy cdf', foldcauchy_cdf(fc_x(i), 0.8_dp), fc_cdf(i), &
                4.0e-12_dp, 4.0e-11_dp, failures)
            call check_close('foldcauchy ppf', foldcauchy_ppf(probs(i), 0.8_dp), fc_ppf(i), &
                3.0e-8_dp, 3.0e-10_dp, failures)
            call check_close('recipinvgauss pdf', recipinvgauss_pdf(ri_x(i), 0.7_dp), ri_pdf(i), &
                5.0e-12_dp, 5.0e-11_dp, failures)
            call check_close('recipinvgauss cdf', recipinvgauss_cdf(ri_x(i), 0.7_dp), ri_cdf(i), &
                8.0e-12_dp, 8.0e-11_dp, failures)
            call check_close('recipinvgauss ppf', recipinvgauss_ppf(probs(i), 0.7_dp), ri_ppf(i), &
                3.0e-8_dp, 3.0e-9_dp, failures)
            call check_close('truncpareto pdf', truncpareto_pdf(tp_x(i), 1.7_dp, 5.0_dp), tp_pdf(i), &
                5.0e-12_dp, 5.0e-11_dp, failures)
            call check_close('truncpareto cdf', truncpareto_cdf(tp_x(i), 1.7_dp, 5.0_dp), tp_cdf(i), &
                5.0e-12_dp, 5.0e-11_dp, failures)
            call check_close('truncpareto ppf', truncpareto_ppf(probs(i), 1.7_dp, 5.0_dp), tp_ppf(i), &
                6.0e-12_dp, 6.0e-11_dp, failures)
        end do
        call check_close('foldnorm logcdf', foldnorm_logcdf(0.2_dp,1.3_dp), &
            -2.675696381221591_dp, 5.0e-12_dp, 5.0e-11_dp, failures)
        call check_close('foldcauchy logsf', foldcauchy_logsf(2.0_dp,0.8_dp), &
            -1.1076632057514375_dp, 5.0e-12_dp, 5.0e-11_dp, failures)
        call check_close('recipinvgauss logcdf', recipinvgauss_logcdf(0.2_dp,0.7_dp), &
            -7.3435791823038308_dp, 8.0e-12_dp, 8.0e-11_dp, failures)
        call check_close('truncpareto logsf', truncpareto_logsf(3.5_dp,1.7_dp,5.0_dp), &
            -2.8508775916351259_dp, 5.0e-12_dp, 5.0e-11_dp, failures)
        call check_close('truncpareto negative b cdf', truncpareto_cdf(2.0_dp,-1.2_dp,5.0_dp), &
            0.21994814_dp, 2.0e-9_dp, 2.0e-8_dp, failures)
    end subroutine test_reference

    subroutine test_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp), parameter :: tol = 1.2e-4_dp
        real(dp), parameter :: fn_data(4) = [0.4_dp, 1.0_dp, 2.0_dp, 3.2_dp]
        real(dp), parameter :: fc_data(4) = [0.4_dp, 0.9_dp, 2.0_dp, 4.0_dp]
        real(dp), parameter :: ri_data(4) = [0.3_dp, 0.7_dp, 2.0_dp, 5.0_dp]
        real(dp), parameter :: tp_data(4) = [1.6_dp, 2.0_dp, 3.4_dp, 5.0_dp]
        real(dp) :: q3(3), q4(4)

        q3 = foldnorm_score(fn_data,1.3_dp,0.1_dp,1.2_dp)
        call check_close('foldnorm score c',q3(1),centered( &
            foldnorm_loglikelihood(fn_data,1.3_dp+h,0.1_dp,1.2_dp), &
            foldnorm_loglikelihood(fn_data,1.3_dp-h,0.1_dp,1.2_dp),h),tol,tol,failures)
        call check_close('foldnorm score loc',q3(2),centered( &
            foldnorm_loglikelihood(fn_data,1.3_dp,0.1_dp+h,1.2_dp), &
            foldnorm_loglikelihood(fn_data,1.3_dp,0.1_dp-h,1.2_dp),h),tol,tol,failures)
        call check_close('foldnorm score scale',q3(3),centered( &
            foldnorm_loglikelihood(fn_data,1.3_dp,0.1_dp,1.2_dp+h), &
            foldnorm_loglikelihood(fn_data,1.3_dp,0.1_dp,1.2_dp-h),h),tol,tol,failures)

        q3 = foldcauchy_score(fc_data,0.8_dp,0.1_dp,1.2_dp)
        call check_close('foldcauchy score c',q3(1),centered( &
            foldcauchy_loglikelihood(fc_data,0.8_dp+h,0.1_dp,1.2_dp), &
            foldcauchy_loglikelihood(fc_data,0.8_dp-h,0.1_dp,1.2_dp),h),tol,tol,failures)
        call check_close('foldcauchy score loc',q3(2),centered( &
            foldcauchy_loglikelihood(fc_data,0.8_dp,0.1_dp+h,1.2_dp), &
            foldcauchy_loglikelihood(fc_data,0.8_dp,0.1_dp-h,1.2_dp),h),tol,tol,failures)
        call check_close('foldcauchy score scale',q3(3),centered( &
            foldcauchy_loglikelihood(fc_data,0.8_dp,0.1_dp,1.2_dp+h), &
            foldcauchy_loglikelihood(fc_data,0.8_dp,0.1_dp,1.2_dp-h),h),tol,tol,failures)

        q3 = recipinvgauss_score(ri_data,0.7_dp,0.1_dp,1.2_dp)
        call check_close('recipinvgauss score mu',q3(1),centered( &
            recipinvgauss_loglikelihood(ri_data,0.7_dp+h,0.1_dp,1.2_dp), &
            recipinvgauss_loglikelihood(ri_data,0.7_dp-h,0.1_dp,1.2_dp),h),tol,tol,failures)
        call check_close('recipinvgauss score loc',q3(2),centered( &
            recipinvgauss_loglikelihood(ri_data,0.7_dp,0.1_dp+h,1.2_dp), &
            recipinvgauss_loglikelihood(ri_data,0.7_dp,0.1_dp-h,1.2_dp),h),tol,tol,failures)
        call check_close('recipinvgauss score scale',q3(3),centered( &
            recipinvgauss_loglikelihood(ri_data,0.7_dp,0.1_dp,1.2_dp+h), &
            recipinvgauss_loglikelihood(ri_data,0.7_dp,0.1_dp,1.2_dp-h),h),tol,tol,failures)

        q4 = truncpareto_score(tp_data,1.7_dp,5.0_dp,0.1_dp,1.2_dp)
        call check_close('truncpareto score b',q4(1),centered( &
            truncpareto_loglikelihood(tp_data,1.7_dp+h,5.0_dp,0.1_dp,1.2_dp), &
            truncpareto_loglikelihood(tp_data,1.7_dp-h,5.0_dp,0.1_dp,1.2_dp),h),tol,tol,failures)
        call check_close('truncpareto score c',q4(2),centered( &
            truncpareto_loglikelihood(tp_data,1.7_dp,5.0_dp+h,0.1_dp,1.2_dp), &
            truncpareto_loglikelihood(tp_data,1.7_dp,5.0_dp-h,0.1_dp,1.2_dp),h),tol,tol,failures)
        call check_close('truncpareto score loc',q4(3),centered( &
            truncpareto_loglikelihood(tp_data,1.7_dp,5.0_dp,0.1_dp+h,1.2_dp), &
            truncpareto_loglikelihood(tp_data,1.7_dp,5.0_dp,0.1_dp-h,1.2_dp),h),tol,tol,failures)
        call check_close('truncpareto score scale',q4(4),centered( &
            truncpareto_loglikelihood(tp_data,1.7_dp,5.0_dp,0.1_dp,1.2_dp+h), &
            truncpareto_loglikelihood(tp_data,1.7_dp,5.0_dp,0.1_dp,1.2_dp-h),h),tol,tol,failures)
    end subroutine test_scores

    subroutine test_rng_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: fn_data(4) = [0.4_dp,1.0_dp,2.0_dp,3.2_dp]
        real(dp), parameter :: fc_data(4) = [0.4_dp,0.9_dp,2.0_dp,4.0_dp]
        real(dp), parameter :: ri_data(4) = [0.3_dp,0.7_dp,2.0_dp,5.0_dp]
        real(dp), parameter :: tp_data(4) = [1.6_dp,2.0_dp,3.4_dp,5.0_dp]
        real(dp) :: actual, expected
        type(rng_state) :: reference, state
        type(fit_result) :: result

        call rng_seed(state,246813)
        call rng_seed(reference,246813)
        actual = foldnorm_rvs(state,1.3_dp,0.1_dp,1.2_dp)
        expected = foldnorm_ppf(rng_uniform(reference),1.3_dp,0.1_dp,1.2_dp)
        call check_close('foldnorm rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = foldcauchy_rvs(state,0.8_dp,0.1_dp,1.2_dp)
        expected = foldcauchy_ppf(rng_uniform(reference),0.8_dp,0.1_dp,1.2_dp)
        call check_close('foldcauchy rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = recipinvgauss_rvs(state,0.7_dp,0.1_dp,1.2_dp)
        expected = recipinvgauss_ppf(rng_uniform(reference),0.7_dp,0.1_dp,1.2_dp)
        call check_close('recipinvgauss rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = truncpareto_rvs(state,1.7_dp,5.0_dp,0.1_dp,1.2_dp)
        expected = truncpareto_ppf(rng_uniform(reference),1.7_dp,5.0_dp,0.1_dp,1.2_dp)
        call check_close('truncpareto rvs',actual,expected,0.0_dp,0.0_dp,failures)

        call foldnorm_fit(fn_data,[1.3_dp,0.1_dp,1.2_dp],[1.3_dp,0.1_dp,1.2_dp],result)
        call check_true('foldnorm fixed fit',result%success,failures)
        call foldcauchy_fit(fc_data,[0.8_dp,0.1_dp,1.2_dp],[0.8_dp,0.1_dp,1.2_dp],result)
        call check_true('foldcauchy fixed fit',result%success,failures)
        call recipinvgauss_fit(ri_data,[0.7_dp,0.1_dp,1.2_dp],[0.7_dp,0.1_dp,1.2_dp],result)
        call check_true('recipinvgauss fixed fit',result%success,failures)
        call truncpareto_fit(tp_data,[1.7_dp,5.0_dp,0.1_dp,1.2_dp], &
            [1.7_dp,5.0_dp,0.1_dp,1.2_dp],result)
        call check_true('truncpareto fixed fit',result%success,failures)

        call check_close('C ABI foldnorm cdf',scifort_foldnorm_cdf_f64(1.0_dp,1.3_dp,0.0_dp,1.0_dp), &
            foldnorm_cdf(1.0_dp,1.3_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI foldcauchy pdf', &
            scifort_foldcauchy_pdf_f64(1.0_dp,0.8_dp,0.0_dp,1.0_dp), &
            foldcauchy_pdf(1.0_dp,0.8_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI recipinvgauss ppf', &
            scifort_recipinvgauss_ppf_f64(0.9_dp,0.7_dp,0.0_dp,1.0_dp), &
            recipinvgauss_ppf(0.9_dp,0.7_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI truncpareto cdf', &
            scifort_truncpareto_cdf_f64(2.0_dp,1.7_dp,5.0_dp,0.0_dp,1.0_dp), &
            truncpareto_cdf(2.0_dp,1.7_dp,5.0_dp),0.0_dp,0.0_dp,failures)
    end subroutine test_rng_fit_c_api

    pure function centered(plus, minus, h) result(value)
        real(dp), intent(in) :: plus !! objective at theta+h
        real(dp), intent(in) :: minus !! objective at theta-h
        real(dp), intent(in) :: h !! positive finite-difference step
        real(dp) :: value
        value = (plus - minus) / (2.0_dp * h)
    end function centered

end program test_folded_reciprocal_four
