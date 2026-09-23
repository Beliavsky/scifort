! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_scipy_transform_four
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
        print '(a)', 'test_scipy_transform_four: PASS'
    else
        print '(a,1x,i0)', 'test_scipy_transform_four: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: dg_x(5) = [-4.0_dp, -1.0_dp, -0.2_dp, 0.5_dp, 3.0_dp]
        real(dp), parameter :: dg_pdf(5) = [0.047588984334845479_dp, 0.15765650437455672_dp, &
            0.043299858075380365_dp, 0.10556504558240531_dp, 0.088998064181960851_dp]
        real(dp), parameter :: dg_cdf(5) = [0.064066659566475823_dp, 0.40477288194092603_dp, &
            0.49599557343049178_dp, 0.52687173824053724_dp, 0.86919092165630707_dp]
        real(dp), parameter :: dg_ppf(5) = [-16.704764935077005_dp, -3.3867280007307086_dp, &
            0.0_dp, 3.3867280007307117_dp, 16.704764935045972_dp]
        real(dp), parameter :: la_x(5) = [-3.0_dp, -0.5_dp, 0.0_dp, 0.7_dp, 3.0_dp]
        real(dp), parameter :: la_pdf(5) = [0.074833712855048806_dp, 0.32566092260229224_dp, &
            0.43701799485861187_dp, 0.13295016681578331_dp, 0.0026643879592227722_dp]
        real(dp), parameter :: la_cdf(5) = [0.12721731185358298_dp, 0.55362356842389682_dp, &
            0.74293059125964001_dp, 0.92179401952012752_dp, 0.99843271296516312_dp]
        real(dp), parameter :: la_ppf(5) = [-22.981208434179042_dp, -3.4092351437296533_dp, &
            -0.67319069259168296_dp, 0.55539760903746538_dp, 7.3277067060618606_dp]
        real(dp), parameter :: tn_x(5) = [-1.0_dp, -0.5_dp, 0.0_dp, 0.8_dp, 1.8_dp]
        real(dp), parameter :: tn_pdf(5) = [0.28064982833997104_dp, 0.40834309074872988_dp, &
            0.46271334160244992_dp, 0.33599884744576874_dp, 0.091570368351757039_dp]
        real(dp), parameter :: tn_cdf(5) = [0.050552754307559615_dp, 0.22439377405232433_dp, &
            0.44646157583653989_dp, 0.7806661912310825_dp, 0.98471295536835834_dp]
        real(dp), parameter :: tn_ppf(5) = [-1.1999955600421406_dp, -0.83703057864282104_dp, &
            0.11596477140566201_dp, 1.2320342147680767_dp, 1.9999840312836441_dp]
        real(dp), parameter :: lu_x(5) = [0.25_dp, 0.4_dp, 1.0_dp, 2.0_dp, 4.5_dp]
        real(dp), parameter :: lu_pdf(5) = [1.2426698691192237_dp, 0.77666866819951474_dp, &
            0.31066746727980593_dp, 0.15533373363990299_dp, 0.06903721495106796_dp]
        real(dp), parameter :: lu_cdf(5) = [0.069323441926606943_dp, 0.21533827903669653_dp, &
            0.5_dp, 0.71533827903669644_dp, 0.96726791544928881_dp]
        real(dp), parameter :: lu_ppf(5) = [0.20000064377620108_dp, 0.27594593229224301_dp, &
            1.0_dp, 3.6238983183884779_dp, 4.9999839056467774_dp]
        integer :: i

        do i = 1, 5
            call check_close('dgamma pdf', dgamma_pdf(dg_x(i), 2.3_dp), dg_pdf(i), &
                3.0e-12_dp, 3.0e-11_dp, failures)
            call check_close('dgamma cdf', dgamma_cdf(dg_x(i), 2.3_dp), dg_cdf(i), &
                3.0e-12_dp, 3.0e-11_dp, failures)
            call check_close('dgamma ppf', dgamma_ppf(probs(i), 2.3_dp), dg_ppf(i), &
                2.0e-9_dp, 2.0e-9_dp, failures)
            call check_close('laplace_asymmetric pdf', &
                laplace_asymmetric_pdf(la_x(i), 1.7_dp), la_pdf(i), &
                3.0e-12_dp, 3.0e-11_dp, failures)
            call check_close('laplace_asymmetric cdf', &
                laplace_asymmetric_cdf(la_x(i), 1.7_dp), la_cdf(i), &
                3.0e-12_dp, 3.0e-11_dp, failures)
            call check_close('laplace_asymmetric ppf', &
                laplace_asymmetric_ppf(probs(i), 1.7_dp), la_ppf(i), &
                3.0e-11_dp, 3.0e-10_dp, failures)
            call check_close('truncnorm pdf', truncnorm_pdf(tn_x(i), -1.2_dp, 2.0_dp), &
                tn_pdf(i), 3.0e-12_dp, 3.0e-11_dp, failures)
            call check_close('truncnorm cdf', truncnorm_cdf(tn_x(i), -1.2_dp, 2.0_dp), &
                tn_cdf(i), 3.0e-12_dp, 3.0e-11_dp, failures)
            call check_close('truncnorm ppf', truncnorm_ppf(probs(i), -1.2_dp, 2.0_dp), &
                tn_ppf(i), 3.0e-11_dp, 3.0e-10_dp, failures)
            call check_close('loguniform pdf', loguniform_pdf(lu_x(i), 0.2_dp, 5.0_dp), &
                lu_pdf(i), 3.0e-12_dp, 3.0e-11_dp, failures)
            call check_close('loguniform cdf', loguniform_cdf(lu_x(i), 0.2_dp, 5.0_dp), &
                lu_cdf(i), 3.0e-12_dp, 3.0e-11_dp, failures)
            call check_close('loguniform ppf', loguniform_ppf(probs(i), 0.2_dp, 5.0_dp), &
                lu_ppf(i), 3.0e-12_dp, 3.0e-11_dp, failures)
        end do
        call check_close('dgamma logsf', dgamma_logsf(3.0_dp, 2.3_dp), &
            -2.0340164360701976_dp, 3.0e-12_dp, 3.0e-11_dp, failures)
        call check_close('laplace_asymmetric logsf', laplace_asymmetric_logsf(3.0_dp, 1.7_dp), &
            -6.4584091576303546_dp, 3.0e-12_dp, 3.0e-11_dp, failures)
        call check_close('truncnorm logsf', truncnorm_logsf(1.8_dp, -1.2_dp, 2.0_dp), &
            -4.1807495653860487_dp, 3.0e-12_dp, 3.0e-11_dp, failures)
        call check_close('loguniform logcdf', loguniform_logcdf(0.25_dp, 0.2_dp, 5.0_dp), &
            -2.6689721626465714_dp, 3.0e-12_dp, 3.0e-11_dp, failures)
    end subroutine test_reference

    subroutine test_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp), parameter :: tol = 8.0e-5_dp
        real(dp), parameter :: dg_data(4) = [-2.0_dp, -0.8_dp, 0.5_dp, 1.7_dp]
        real(dp), parameter :: la_data(4) = [-2.0_dp, -0.4_dp, 0.5_dp, 2.0_dp]
        real(dp), parameter :: tn_data(4) = [-0.94_dp, -0.16_dp, 0.62_dp, 1.66_dp]
        real(dp), parameter :: lu_data(4) = [0.46_dp, 0.82_dp, 1.9_dp, 3.7_dp]
        real(dp) :: q3(3), q4(4)

        q3 = dgamma_score(dg_data, 2.3_dp, 0.1_dp, 1.2_dp)
        call check_close('dgamma score a', q3(1), centered( &
            dgamma_loglikelihood(dg_data,2.3_dp+h,0.1_dp,1.2_dp), &
            dgamma_loglikelihood(dg_data,2.3_dp-h,0.1_dp,1.2_dp),h), tol,tol,failures)
        call check_close('dgamma score loc', q3(2), centered( &
            dgamma_loglikelihood(dg_data,2.3_dp,0.1_dp+h,1.2_dp), &
            dgamma_loglikelihood(dg_data,2.3_dp,0.1_dp-h,1.2_dp),h), tol,tol,failures)
        call check_close('dgamma score scale', q3(3), centered( &
            dgamma_loglikelihood(dg_data,2.3_dp,0.1_dp,1.2_dp+h), &
            dgamma_loglikelihood(dg_data,2.3_dp,0.1_dp,1.2_dp-h),h), tol,tol,failures)

        q3 = laplace_asymmetric_score(la_data, 1.7_dp, 0.1_dp, 1.2_dp)
        call check_close('laplace_asymmetric score kappa', q3(1), centered( &
            laplace_asymmetric_loglikelihood(la_data,1.7_dp+h,0.1_dp,1.2_dp), &
            laplace_asymmetric_loglikelihood(la_data,1.7_dp-h,0.1_dp,1.2_dp),h),tol,tol,failures)
        call check_close('laplace_asymmetric score loc', q3(2), centered( &
            laplace_asymmetric_loglikelihood(la_data,1.7_dp,0.1_dp+h,1.2_dp), &
            laplace_asymmetric_loglikelihood(la_data,1.7_dp,0.1_dp-h,1.2_dp),h),tol,tol,failures)
        call check_close('laplace_asymmetric score scale', q3(3), centered( &
            laplace_asymmetric_loglikelihood(la_data,1.7_dp,0.1_dp,1.2_dp+h), &
            laplace_asymmetric_loglikelihood(la_data,1.7_dp,0.1_dp,1.2_dp-h),h),tol,tol,failures)

        q4 = truncnorm_score(tn_data, -1.2_dp, 2.0_dp, 0.1_dp, 1.3_dp)
        call check_close('truncnorm score a', q4(1), centered( &
            truncnorm_loglikelihood(tn_data,-1.2_dp+h,2.0_dp,0.1_dp,1.3_dp), &
            truncnorm_loglikelihood(tn_data,-1.2_dp-h,2.0_dp,0.1_dp,1.3_dp),h),tol,tol,failures)
        call check_close('truncnorm score b', q4(2), centered( &
            truncnorm_loglikelihood(tn_data,-1.2_dp,2.0_dp+h,0.1_dp,1.3_dp), &
            truncnorm_loglikelihood(tn_data,-1.2_dp,2.0_dp-h,0.1_dp,1.3_dp),h),tol,tol,failures)
        call check_close('truncnorm score loc', q4(3), centered( &
            truncnorm_loglikelihood(tn_data,-1.2_dp,2.0_dp,0.1_dp+h,1.3_dp), &
            truncnorm_loglikelihood(tn_data,-1.2_dp,2.0_dp,0.1_dp-h,1.3_dp),h),tol,tol,failures)
        call check_close('truncnorm score scale', q4(4), centered( &
            truncnorm_loglikelihood(tn_data,-1.2_dp,2.0_dp,0.1_dp,1.3_dp+h), &
            truncnorm_loglikelihood(tn_data,-1.2_dp,2.0_dp,0.1_dp,1.3_dp-h),h),tol,tol,failures)

        q4 = loguniform_score(lu_data, 0.2_dp, 5.0_dp, 0.1_dp, 1.2_dp)
        call check_close('loguniform score a', q4(1), centered( &
            loguniform_loglikelihood(lu_data,0.2_dp+h,5.0_dp,0.1_dp,1.2_dp), &
            loguniform_loglikelihood(lu_data,0.2_dp-h,5.0_dp,0.1_dp,1.2_dp),h),tol,tol,failures)
        call check_close('loguniform score b', q4(2), centered( &
            loguniform_loglikelihood(lu_data,0.2_dp,5.0_dp+h,0.1_dp,1.2_dp), &
            loguniform_loglikelihood(lu_data,0.2_dp,5.0_dp-h,0.1_dp,1.2_dp),h),tol,tol,failures)
        call check_close('loguniform score loc', q4(3), centered( &
            loguniform_loglikelihood(lu_data,0.2_dp,5.0_dp,0.1_dp+h,1.2_dp), &
            loguniform_loglikelihood(lu_data,0.2_dp,5.0_dp,0.1_dp-h,1.2_dp),h),tol,tol,failures)
        call check_close('loguniform score scale', q4(4), centered( &
            loguniform_loglikelihood(lu_data,0.2_dp,5.0_dp,0.1_dp,1.2_dp+h), &
            loguniform_loglikelihood(lu_data,0.2_dp,5.0_dp,0.1_dp,1.2_dp-h),h),tol,tol,failures)
    end subroutine test_scores

    subroutine test_rng_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: dg_data(4) = [-2.0_dp, -0.8_dp, 0.5_dp, 1.7_dp]
        real(dp), parameter :: la_data(4) = [-2.0_dp, -0.4_dp, 0.5_dp, 2.0_dp]
        real(dp), parameter :: tn_data(4) = [-0.94_dp, -0.16_dp, 0.62_dp, 1.66_dp]
        real(dp), parameter :: lu_data(4) = [0.46_dp, 0.82_dp, 1.9_dp, 3.7_dp]
        real(dp) :: actual, expected
        type(rng_state) :: reference, state
        type(fit_result) :: result

        call rng_seed(state, 987654)
        call rng_seed(reference, 987654)
        actual = dgamma_rvs(state, 2.3_dp, 0.1_dp, 1.2_dp)
        expected = dgamma_ppf(rng_uniform(reference), 2.3_dp, 0.1_dp, 1.2_dp)
        call check_close('dgamma rvs', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = laplace_asymmetric_rvs(state, 1.7_dp, 0.1_dp, 1.2_dp)
        expected = laplace_asymmetric_ppf(rng_uniform(reference), 1.7_dp, 0.1_dp, 1.2_dp)
        call check_close('laplace_asymmetric rvs', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = truncnorm_rvs(state, -1.2_dp, 2.0_dp, 0.1_dp, 1.3_dp)
        expected = truncnorm_ppf(rng_uniform(reference), -1.2_dp, 2.0_dp, 0.1_dp, 1.3_dp)
        call check_close('truncnorm rvs', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = loguniform_rvs(state, 0.2_dp, 5.0_dp, 0.1_dp, 1.2_dp)
        expected = loguniform_ppf(rng_uniform(reference), 0.2_dp, 5.0_dp, 0.1_dp, 1.2_dp)
        call check_close('loguniform rvs', actual, expected, 0.0_dp, 0.0_dp, failures)

        call dgamma_fit(dg_data, [2.3_dp,0.1_dp,1.2_dp], [2.3_dp,0.1_dp,1.2_dp], result)
        call check_true('dgamma fixed fit', result%success, failures)
        call laplace_asymmetric_fit(la_data, [1.7_dp,0.1_dp,1.2_dp], &
            [1.7_dp,0.1_dp,1.2_dp], result)
        call check_true('laplace_asymmetric fixed fit', result%success, failures)
        call truncnorm_fit(tn_data, [-1.2_dp,2.0_dp,0.1_dp,1.3_dp], &
            [-1.2_dp,2.0_dp,0.1_dp,1.3_dp], result)
        call check_true('truncnorm fixed fit', result%success, failures)
        call loguniform_fit(lu_data, [0.2_dp,5.0_dp,0.1_dp,1.2_dp], &
            [0.2_dp,5.0_dp,0.1_dp,1.2_dp], result)
        call check_true('loguniform fixed fit', result%success, failures)

        call check_close('C ABI dgamma cdf', scifort_dgamma_cdf_f64(-1.0_dp,2.3_dp,0.0_dp,1.0_dp), &
            dgamma_cdf(-1.0_dp,2.3_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI laplace asymmetric pdf', &
            scifort_laplace_asymmetric_pdf_f64(0.7_dp,1.7_dp,0.0_dp,1.0_dp), &
            laplace_asymmetric_pdf(0.7_dp,1.7_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI truncnorm ppf', &
            scifort_truncnorm_ppf_f64(0.9_dp,-1.2_dp,2.0_dp,0.0_dp,1.0_dp), &
            truncnorm_ppf(0.9_dp,-1.2_dp,2.0_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI loguniform cdf', &
            scifort_loguniform_cdf_f64(2.0_dp,0.2_dp,5.0_dp,0.0_dp,1.0_dp), &
            loguniform_cdf(2.0_dp,0.2_dp,5.0_dp), 0.0_dp, 0.0_dp, failures)
    end subroutine test_rng_fit_c_api

    pure function centered(plus, minus, h) result(value)
        real(dp), intent(in) :: plus !! objective at theta+h
        real(dp), intent(in) :: minus !! objective at theta-h
        real(dp), intent(in) :: h !! positive finite-difference step
        real(dp) :: value
        value = (plus - minus) / (2.0_dp * h)
    end function centered

end program test_scipy_transform_four
