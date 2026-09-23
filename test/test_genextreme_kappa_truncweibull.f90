! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_genextreme_kappa_truncweibull
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_c_api
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state, rng_uniform
    use scifort_stats
    use test_genextreme_kappa_truncweibull_reference
    use test_checks, only : check_close, check_true
    implicit none

    integer :: failures

    failures = 0
    call test_reference(failures)
    call test_identities_endpoints(failures)
    call test_scores(failures)
    call test_rng_fit_c_api(failures)

    if (failures == 0) then
        print '(a)', 'test_genextreme_kappa_truncweibull: PASS'
    else
        print '(a,1x,i0)', 'test_genextreme_kappa_truncweibull: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        integer :: i

        do i = 1, size(reference_probs)
            call check_close('genextreme pdf',genextreme_pdf(ge_x(i),0.35_dp),ge_pdf(i), &
                1.0e-11_dp,1.0e-10_dp,failures)
            call check_close('genextreme logpdf',genextreme_logpdf(ge_x(i),0.35_dp),ge_logpdf(i), &
                1.0e-11_dp,1.0e-10_dp,failures)
            call check_close('genextreme cdf',genextreme_cdf(ge_x(i),0.35_dp),ge_cdf(i), &
                1.0e-11_dp,1.0e-10_dp,failures)
            call check_close('genextreme sf',genextreme_sf(ge_x(i),0.35_dp),ge_sf(i), &
                1.0e-11_dp,1.0e-10_dp,failures)
            call check_close('genextreme logcdf',genextreme_logcdf(ge_x(i),0.35_dp),ge_logcdf(i), &
                1.0e-11_dp,1.0e-10_dp,failures)
            call check_close('genextreme logsf',genextreme_logsf(ge_x(i),0.35_dp),ge_logsf(i), &
                1.0e-11_dp,1.0e-10_dp,failures)
            call check_close('genextreme ppf',genextreme_ppf(reference_probs(i),0.35_dp),ge_ppf(i), &
                3.0e-10_dp,3.0e-9_dp,failures)
            call check_close('genextreme isf',genextreme_isf(reference_probs(i),0.35_dp),ge_isf(i), &
                3.0e-10_dp,3.0e-9_dp,failures)

            call check_close('kappa3 pdf',kappa3_pdf(k3_x(i),2.4_dp),k3_pdf(i), &
                1.0e-11_dp,1.0e-10_dp,failures)
            call check_close('kappa3 logpdf',kappa3_logpdf(k3_x(i),2.4_dp),k3_logpdf(i), &
                1.0e-11_dp,1.0e-10_dp,failures)
            call check_close('kappa3 cdf',kappa3_cdf(k3_x(i),2.4_dp),k3_cdf(i), &
                1.0e-11_dp,1.0e-10_dp,failures)
            call check_close('kappa3 sf',kappa3_sf(k3_x(i),2.4_dp),k3_sf(i), &
                1.0e-11_dp,1.0e-10_dp,failures)
            call check_close('kappa3 logcdf',kappa3_logcdf(k3_x(i),2.4_dp),k3_logcdf(i), &
                1.0e-11_dp,1.0e-10_dp,failures)
            call check_close('kappa3 logsf',kappa3_logsf(k3_x(i),2.4_dp),k3_logsf(i), &
                1.0e-11_dp,1.0e-10_dp,failures)
            call check_close('kappa3 ppf',kappa3_ppf(reference_probs(i),2.4_dp),k3_ppf(i), &
                3.0e-9_dp,3.0e-9_dp,failures)
            call check_close('kappa3 isf',kappa3_isf(reference_probs(i),2.4_dp),k3_isf(i), &
                3.0e-9_dp,3.0e-9_dp,failures)

            call check_close('kappa4 pdf',kappa4_pdf(k4_x(i),0.5_dp,0.3_dp),k4_pdf(i), &
                1.0e-11_dp,1.0e-10_dp,failures)
            call check_close('kappa4 logpdf',kappa4_logpdf(k4_x(i),0.5_dp,0.3_dp),k4_logpdf(i), &
                1.0e-11_dp,1.0e-10_dp,failures)
            call check_close('kappa4 cdf',kappa4_cdf(k4_x(i),0.5_dp,0.3_dp),k4_cdf(i), &
                1.0e-11_dp,1.0e-10_dp,failures)
            call check_close('kappa4 sf',kappa4_sf(k4_x(i),0.5_dp,0.3_dp),k4_sf(i), &
                1.0e-11_dp,1.0e-10_dp,failures)
            call check_close('kappa4 logcdf',kappa4_logcdf(k4_x(i),0.5_dp,0.3_dp),k4_logcdf(i), &
                1.0e-11_dp,1.0e-10_dp,failures)
            call check_close('kappa4 logsf',kappa4_logsf(k4_x(i),0.5_dp,0.3_dp),k4_logsf(i), &
                1.0e-11_dp,1.0e-10_dp,failures)
            call check_close('kappa4 ppf',kappa4_ppf(reference_probs(i),0.5_dp,0.3_dp),k4_ppf(i), &
                3.0e-9_dp,3.0e-9_dp,failures)
            call check_close('kappa4 isf',kappa4_isf(reference_probs(i),0.5_dp,0.3_dp),k4_isf(i), &
                3.0e-9_dp,3.0e-9_dp,failures)

            call check_close('truncweibull pdf', &
                truncweibull_min_pdf(tw_x(i),1.7_dp,0.2_dp,2.4_dp),tw_pdf(i), &
                1.0e-11_dp,1.0e-10_dp,failures)
            call check_close('truncweibull logpdf', &
                truncweibull_min_logpdf(tw_x(i),1.7_dp,0.2_dp,2.4_dp),tw_logpdf(i), &
                1.0e-11_dp,1.0e-10_dp,failures)
            call check_close('truncweibull cdf', &
                truncweibull_min_cdf(tw_x(i),1.7_dp,0.2_dp,2.4_dp),tw_cdf(i), &
                1.0e-11_dp,1.0e-10_dp,failures)
            call check_close('truncweibull sf', &
                truncweibull_min_sf(tw_x(i),1.7_dp,0.2_dp,2.4_dp),tw_sf(i), &
                1.0e-11_dp,1.0e-10_dp,failures)
            call check_close('truncweibull logcdf', &
                truncweibull_min_logcdf(tw_x(i),1.7_dp,0.2_dp,2.4_dp),tw_logcdf(i), &
                1.0e-11_dp,1.0e-10_dp,failures)
            call check_close('truncweibull logsf', &
                truncweibull_min_logsf(tw_x(i),1.7_dp,0.2_dp,2.4_dp),tw_logsf(i), &
                1.0e-11_dp,1.0e-10_dp,failures)
            call check_close('truncweibull ppf', &
                truncweibull_min_ppf(reference_probs(i),1.7_dp,0.2_dp,2.4_dp),tw_ppf(i), &
                3.0e-10_dp,3.0e-9_dp,failures)
            call check_close('truncweibull isf', &
                truncweibull_min_isf(reference_probs(i),1.7_dp,0.2_dp,2.4_dp),tw_isf(i), &
                3.0e-10_dp,3.0e-9_dp,failures)
        end do
    end subroutine test_reference

    subroutine test_identities_endpoints(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: xvals(5) = [-0.4_dp,0.0_dp,0.4_dp,0.8_dp,1.3_dp]
        real(dp), parameter :: pvals(5) = [1.0e-12_dp,1.0e-6_dp,0.2_dp,0.8_dp,1.0_dp-1.0e-12_dp]
        real(dp) :: a3(3), lower, upper
        integer :: i

        do i = 1, size(xvals)
            call check_close('genextreme-gumbel pdf',genextreme_pdf(xvals(i),0.0_dp), &
                gumbel_r_pdf(xvals(i)),2.0e-14_dp,2.0e-13_dp,failures)
            call check_close('kappa4-genextreme pdf',kappa4_pdf(xvals(i),0.0_dp,0.3_dp), &
                genextreme_pdf(xvals(i),0.3_dp),2.0e-14_dp,2.0e-13_dp,failures)
            call check_close('kappa4-genpareto cdf',kappa4_cdf(xvals(i),1.0_dp,0.3_dp), &
                genpareto_cdf(xvals(i),-0.3_dp),2.0e-14_dp,2.0e-13_dp,failures)
            call check_close('kappa4-logistic cdf',kappa4_cdf(xvals(i),-1.0_dp,0.0_dp), &
                logistic_cdf(xvals(i)),2.0e-14_dp,2.0e-13_dp,failures)
        end do
        do i = 1, size(pvals)
            call check_close('genextreme cdf-ppf', &
                genextreme_cdf(genextreme_ppf(pvals(i),0.35_dp),0.35_dp),pvals(i), &
                3.0e-12_dp,3.0e-10_dp,failures)
            call check_close('genextreme sf-isf', &
                genextreme_sf(genextreme_isf(pvals(i),0.35_dp),0.35_dp),pvals(i), &
                3.0e-12_dp,3.0e-10_dp,failures)
            call check_close('kappa3 cdf-ppf',kappa3_cdf(kappa3_ppf(pvals(i),2.4_dp),2.4_dp), &
                pvals(i),3.0e-12_dp,3.0e-10_dp,failures)
            call check_close('kappa3 sf-isf',kappa3_sf(kappa3_isf(pvals(i),2.4_dp),2.4_dp), &
                pvals(i),3.0e-12_dp,3.0e-10_dp,failures)
            call check_close('kappa4 cdf-ppf', &
                kappa4_cdf(kappa4_ppf(pvals(i),0.5_dp,0.3_dp),0.5_dp,0.3_dp),pvals(i), &
                3.0e-12_dp,3.0e-10_dp,failures)
            call check_close('kappa4 sf-isf', &
                kappa4_sf(kappa4_isf(pvals(i),0.5_dp,0.3_dp),0.5_dp,0.3_dp),pvals(i), &
                3.0e-12_dp,3.0e-10_dp,failures)
            call check_close('truncweibull cdf-ppf', &
                truncweibull_min_cdf(truncweibull_min_ppf(pvals(i),1.7_dp,0.2_dp,2.4_dp), &
                1.7_dp,0.2_dp,2.4_dp),pvals(i),3.0e-12_dp,3.0e-10_dp,failures)
            call check_close('truncweibull sf-isf', &
                truncweibull_min_sf(truncweibull_min_isf(pvals(i),1.7_dp,0.2_dp,2.4_dp), &
                1.7_dp,0.2_dp,2.4_dp),pvals(i),3.0e-12_dp,3.0e-10_dp,failures)
        end do

        upper = 0.2_dp + 1.3_dp / 0.35_dp
        call check_close('genextreme upper endpoint',genextreme_cdf(upper,0.35_dp,0.2_dp,1.3_dp), &
            1.0_dp,0.0_dp,0.0_dp,failures)
        lower = 0.1_dp + 1.2_dp * 0.2_dp
        upper = 0.1_dp + 1.2_dp * 2.4_dp
        call check_close('truncweibull lower endpoint', &
            truncweibull_min_cdf(lower,1.7_dp,0.2_dp,2.4_dp,0.1_dp,1.2_dp), &
            0.0_dp,0.0_dp,0.0_dp,failures)
        call check_close('truncweibull upper endpoint', &
            truncweibull_min_cdf(upper,1.7_dp,0.2_dp,2.4_dp,0.1_dp,1.2_dp), &
            1.0_dp,0.0_dp,0.0_dp,failures)
        call check_true('genextreme invalid scale', &
            ieee_is_nan(genextreme_pdf(0.0_dp,0.3_dp,0.0_dp,-1.0_dp)),failures)
        call check_true('kappa3 invalid shape',ieee_is_nan(kappa3_pdf(1.0_dp,0.0_dp)),failures)
        call check_true('kappa4 invalid scale',ieee_is_nan(kappa4_pdf(0.0_dp,0.5_dp,0.3_dp, &
            0.0_dp,-1.0_dp)),failures)
        call check_true('truncweibull invalid bounds', &
            ieee_is_nan(truncweibull_min_pdf(1.0_dp,1.7_dp,2.4_dp,0.2_dp)),failures)
        call check_true('truncweibull singular endpoint', &
            .not. ieee_is_finite(truncweibull_min_pdf(0.0_dp,0.7_dp,0.0_dp,2.0_dp)),failures)

        a3 = kappa3_pdf([0.3_dp,1.0_dp,3.0_dp],2.4_dp)
        call check_close('kappa3 elemental 1',a3(1),kappa3_pdf(0.3_dp,2.4_dp),0.0_dp,0.0_dp,failures)
        call check_close('kappa3 elemental 3',a3(3),kappa3_pdf(3.0_dp,2.4_dp),0.0_dp,0.0_dp,failures)
    end subroutine test_identities_endpoints

    subroutine test_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: fd = 1.0e-6_dp, tol = 7.0e-4_dp
        real(dp), parameter :: ge_data(4) = [-1.0_dp,0.2_dp,1.0_dp,2.0_dp]
        real(dp), parameter :: k3_data(4) = [0.5_dp,0.9_dp,1.8_dp,4.0_dp]
        real(dp), parameter :: k4_data(4) = [-0.2_dp,0.5_dp,1.5_dp,3.0_dp]
        real(dp), parameter :: tw_data(4) = [0.6_dp,1.0_dp,1.8_dp,2.6_dp]
        real(dp) :: s3(3), s4(4), s5(5)

        s3 = genextreme_score(ge_data,0.35_dp,0.1_dp,1.2_dp)
        call check_close('genextreme score c',s3(1),centered( &
            genextreme_loglikelihood(ge_data,0.35_dp+fd,0.1_dp,1.2_dp), &
            genextreme_loglikelihood(ge_data,0.35_dp-fd,0.1_dp,1.2_dp),fd),tol,tol,failures)
        call check_close('genextreme score loc',s3(2),centered( &
            genextreme_loglikelihood(ge_data,0.35_dp,0.1_dp+fd,1.2_dp), &
            genextreme_loglikelihood(ge_data,0.35_dp,0.1_dp-fd,1.2_dp),fd),tol,tol,failures)
        call check_close('genextreme score scale',s3(3),centered( &
            genextreme_loglikelihood(ge_data,0.35_dp,0.1_dp,1.2_dp+fd), &
            genextreme_loglikelihood(ge_data,0.35_dp,0.1_dp,1.2_dp-fd),fd),tol,tol,failures)

        s3 = kappa3_score(k3_data,2.4_dp,0.1_dp,1.2_dp)
        call check_close('kappa3 score a',s3(1),centered( &
            kappa3_loglikelihood(k3_data,2.4_dp+fd,0.1_dp,1.2_dp), &
            kappa3_loglikelihood(k3_data,2.4_dp-fd,0.1_dp,1.2_dp),fd),tol,tol,failures)
        call check_close('kappa3 score loc',s3(2),centered( &
            kappa3_loglikelihood(k3_data,2.4_dp,0.1_dp+fd,1.2_dp), &
            kappa3_loglikelihood(k3_data,2.4_dp,0.1_dp-fd,1.2_dp),fd),tol,tol,failures)
        call check_close('kappa3 score scale',s3(3),centered( &
            kappa3_loglikelihood(k3_data,2.4_dp,0.1_dp,1.2_dp+fd), &
            kappa3_loglikelihood(k3_data,2.4_dp,0.1_dp,1.2_dp-fd),fd),tol,tol,failures)

        s4 = kappa4_score(k4_data,0.5_dp,0.3_dp,0.1_dp,1.2_dp)
        call check_close('kappa4 score h',s4(1),centered( &
            kappa4_loglikelihood(k4_data,0.5_dp+fd,0.3_dp,0.1_dp,1.2_dp), &
            kappa4_loglikelihood(k4_data,0.5_dp-fd,0.3_dp,0.1_dp,1.2_dp),fd),tol,tol,failures)
        call check_close('kappa4 score k',s4(2),centered( &
            kappa4_loglikelihood(k4_data,0.5_dp,0.3_dp+fd,0.1_dp,1.2_dp), &
            kappa4_loglikelihood(k4_data,0.5_dp,0.3_dp-fd,0.1_dp,1.2_dp),fd),tol,tol,failures)
        call check_close('kappa4 score loc',s4(3),centered( &
            kappa4_loglikelihood(k4_data,0.5_dp,0.3_dp,0.1_dp+fd,1.2_dp), &
            kappa4_loglikelihood(k4_data,0.5_dp,0.3_dp,0.1_dp-fd,1.2_dp),fd),tol,tol,failures)
        call check_close('kappa4 score scale',s4(4),centered( &
            kappa4_loglikelihood(k4_data,0.5_dp,0.3_dp,0.1_dp,1.2_dp+fd), &
            kappa4_loglikelihood(k4_data,0.5_dp,0.3_dp,0.1_dp,1.2_dp-fd),fd),tol,tol,failures)

        s5 = truncweibull_min_score(tw_data,1.7_dp,0.2_dp,2.4_dp,0.1_dp,1.2_dp)
        call check_close('truncweibull score c',s5(1),centered( &
            truncweibull_min_loglikelihood(tw_data,1.7_dp+fd,0.2_dp,2.4_dp,0.1_dp,1.2_dp), &
            truncweibull_min_loglikelihood(tw_data,1.7_dp-fd,0.2_dp,2.4_dp,0.1_dp,1.2_dp),fd), &
            tol,tol,failures)
        call check_close('truncweibull score a',s5(2),centered( &
            truncweibull_min_loglikelihood(tw_data,1.7_dp,0.2_dp+fd,2.4_dp,0.1_dp,1.2_dp), &
            truncweibull_min_loglikelihood(tw_data,1.7_dp,0.2_dp-fd,2.4_dp,0.1_dp,1.2_dp),fd), &
            tol,tol,failures)
        call check_close('truncweibull score b',s5(3),centered( &
            truncweibull_min_loglikelihood(tw_data,1.7_dp,0.2_dp,2.4_dp+fd,0.1_dp,1.2_dp), &
            truncweibull_min_loglikelihood(tw_data,1.7_dp,0.2_dp,2.4_dp-fd,0.1_dp,1.2_dp),fd), &
            tol,tol,failures)
        call check_close('truncweibull score loc',s5(4),centered( &
            truncweibull_min_loglikelihood(tw_data,1.7_dp,0.2_dp,2.4_dp,0.1_dp+fd,1.2_dp), &
            truncweibull_min_loglikelihood(tw_data,1.7_dp,0.2_dp,2.4_dp,0.1_dp-fd,1.2_dp),fd), &
            tol,tol,failures)
        call check_close('truncweibull score scale',s5(5),centered( &
            truncweibull_min_loglikelihood(tw_data,1.7_dp,0.2_dp,2.4_dp,0.1_dp,1.2_dp+fd), &
            truncweibull_min_loglikelihood(tw_data,1.7_dp,0.2_dp,2.4_dp,0.1_dp,1.2_dp-fd),fd), &
            tol,tol,failures)
    end subroutine test_scores

    subroutine test_rng_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: ge_data(4) = [-1.0_dp,0.2_dp,1.0_dp,2.0_dp]
        real(dp), parameter :: k3_data(4) = [0.5_dp,0.9_dp,1.8_dp,4.0_dp]
        real(dp), parameter :: k4_data(4) = [-0.2_dp,0.5_dp,1.5_dp,3.0_dp]
        real(dp), parameter :: tw_data(4) = [0.6_dp,1.0_dp,1.8_dp,2.6_dp]
        real(dp) :: actual, expected
        type(rng_state) :: reference, state
        type(fit_result) :: result

        call rng_seed(state,24680)
        call rng_seed(reference,24680)
        actual = genextreme_rvs(state,0.35_dp,0.1_dp,1.2_dp)
        expected = genextreme_ppf(rng_uniform(reference),0.35_dp,0.1_dp,1.2_dp)
        call check_close('genextreme rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = kappa3_rvs(state,2.4_dp,0.1_dp,1.2_dp)
        expected = kappa3_ppf(rng_uniform(reference),2.4_dp,0.1_dp,1.2_dp)
        call check_close('kappa3 rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = kappa4_rvs(state,0.5_dp,0.3_dp,0.1_dp,1.2_dp)
        expected = kappa4_ppf(rng_uniform(reference),0.5_dp,0.3_dp,0.1_dp,1.2_dp)
        call check_close('kappa4 rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = truncweibull_min_rvs(state,1.7_dp,0.2_dp,2.4_dp,0.1_dp,1.2_dp)
        expected = truncweibull_min_ppf(rng_uniform(reference),1.7_dp,0.2_dp,2.4_dp,0.1_dp,1.2_dp)
        call check_close('truncweibull rvs',actual,expected,0.0_dp,0.0_dp,failures)

        call rng_seed(state,13579)
        call rng_seed(reference,13579)
        actual = kappa3_rvs(state,-1.0_dp)
        call check_true('invalid kappa3 rvs NaN',ieee_is_nan(actual),failures)
        call check_close('invalid rvs preserves state',rng_uniform(state),rng_uniform(reference), &
            0.0_dp,0.0_dp,failures)

        call genextreme_fit(ge_data,[0.35_dp,0.1_dp,1.2_dp],[0.35_dp,0.1_dp,1.2_dp],result)
        call check_true('genextreme fixed fit',result%success,failures)
        call kappa3_fit(k3_data,[2.4_dp,0.1_dp,1.2_dp],[2.4_dp,0.1_dp,1.2_dp],result)
        call check_true('kappa3 fixed fit',result%success,failures)
        call kappa4_fit(k4_data,[0.5_dp,0.3_dp,0.1_dp,1.2_dp], &
            [0.5_dp,0.3_dp,0.1_dp,1.2_dp],result)
        call check_true('kappa4 fixed fit',result%success,failures)
        call truncweibull_min_fit(tw_data,[1.7_dp,0.2_dp,2.4_dp,0.1_dp,1.2_dp], &
            [1.7_dp,0.2_dp,2.4_dp,0.1_dp,1.2_dp],result)
        call check_true('truncweibull fixed fit',result%success,failures)

        call check_close('C ABI genextreme cdf', &
            scifort_genextreme_cdf_f64(1.0_dp,0.35_dp,0.0_dp,1.0_dp), &
            genextreme_cdf(1.0_dp,0.35_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI kappa3 pdf',scifort_kappa3_pdf_f64(1.0_dp,2.4_dp,0.0_dp,1.0_dp), &
            kappa3_pdf(1.0_dp,2.4_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI kappa4 ppf', &
            scifort_kappa4_ppf_f64(0.9_dp,0.5_dp,0.3_dp,0.0_dp,1.0_dp), &
            kappa4_ppf(0.9_dp,0.5_dp,0.3_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI truncweibull cdf', &
            scifort_truncweibull_min_cdf_f64(1.0_dp,1.7_dp,0.2_dp,2.4_dp,0.0_dp,1.0_dp), &
            truncweibull_min_cdf(1.0_dp,1.7_dp,0.2_dp,2.4_dp),0.0_dp,0.0_dp,failures)
    end subroutine test_rng_fit_c_api

    pure function centered(plus, minus, h) result(value)
        real(dp), intent(in) :: plus !! objective at theta+h
        real(dp), intent(in) :: minus !! objective at theta-h
        real(dp), intent(in) :: h !! positive finite-difference step
        real(dp) :: value
        value = (plus - minus) / (2.0_dp * h)
    end function centered

end program test_genextreme_kappa_truncweibull
