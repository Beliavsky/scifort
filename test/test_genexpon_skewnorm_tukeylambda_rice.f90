! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_genexpon_skewnorm_tukeylambda_rice
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_positive_inf, ieee_value
    use scifort_c_api
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state, rng_uniform
    use scifort_stats
    use test_checks, only : check_close, check_true
    use test_genexpon_skewnorm_tukeylambda_rice_reference
    implicit none

    integer :: failures

    failures = 0
    call test_reference(failures)
    call test_identities_endpoints(failures)
    call test_scores(failures)
    call test_rng_fit_c_api(failures)

    if (failures == 0) then
        print '(a)', 'test_genexpon_skewnorm_tukeylambda_rice: PASS'
    else
        print '(a,1x,i0)', 'test_genexpon_skewnorm_tukeylambda_rice: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        integer :: i

        do i = 1, size(reference_probs)
            call check_close('genexpon pdf',genexpon_pdf(gx_x(i),1.3_dp,0.8_dp,1.7_dp), &
                gx_pdf(i),2.0e-13_dp,2.0e-12_dp,failures)
            call check_close('genexpon logpdf',genexpon_logpdf(gx_x(i),1.3_dp,0.8_dp,1.7_dp), &
                gx_logpdf(i),2.0e-13_dp,2.0e-12_dp,failures)
            call check_close('genexpon cdf',genexpon_cdf(gx_x(i),1.3_dp,0.8_dp,1.7_dp), &
                gx_cdf(i),2.0e-13_dp,2.0e-12_dp,failures)
            call check_close('genexpon sf',genexpon_sf(gx_x(i),1.3_dp,0.8_dp,1.7_dp), &
                gx_sf(i),2.0e-13_dp,2.0e-12_dp,failures)
            call check_close('genexpon logcdf',genexpon_logcdf(gx_x(i),1.3_dp,0.8_dp,1.7_dp), &
                gx_logcdf(i),2.0e-13_dp,2.0e-12_dp,failures)
            call check_close('genexpon logsf',genexpon_logsf(gx_x(i),1.3_dp,0.8_dp,1.7_dp), &
                gx_logsf(i),2.0e-13_dp,2.0e-12_dp,failures)
            call check_close('genexpon ppf',genexpon_ppf(reference_probs(i),1.3_dp,0.8_dp,1.7_dp), &
                gx_ppf(i),2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('genexpon isf',genexpon_isf(reference_probs(i),1.3_dp,0.8_dp,1.7_dp), &
                gx_isf(i),2.0e-11_dp,2.0e-10_dp,failures)

            call check_close('skewnorm pdf',skewnorm_pdf(sn_x(i),3.5_dp),sn_pdf(i), &
                2.0e-12_dp,2.0e-10_dp,failures)
            call check_close('skewnorm logpdf',skewnorm_logpdf(sn_x(i),3.5_dp),sn_logpdf(i), &
                2.0e-12_dp,2.0e-10_dp,failures)
            call check_close('skewnorm cdf',skewnorm_cdf(sn_x(i),3.5_dp),sn_cdf(i), &
                2.0e-14_dp,1.0e-8_dp,failures)
            call check_close('skewnorm sf',skewnorm_sf(sn_x(i),3.5_dp),sn_sf(i), &
                2.0e-14_dp,1.0e-8_dp,failures)
            call check_close('skewnorm logcdf',skewnorm_logcdf(sn_x(i),3.5_dp),sn_logcdf(i), &
                1.0e-8_dp,1.0e-8_dp,failures)
            call check_close('skewnorm logsf',skewnorm_logsf(sn_x(i),3.5_dp),sn_logsf(i), &
                1.0e-8_dp,1.0e-8_dp,failures)
            call check_close('skewnorm ppf',skewnorm_ppf(reference_probs(i),3.5_dp),sn_ppf(i), &
                2.0e-9_dp,2.0e-8_dp,failures)
            call check_close('skewnorm isf',skewnorm_isf(reference_probs(i),3.5_dp),sn_isf(i), &
                2.0e-9_dp,2.0e-8_dp,failures)

            call check_close('tukeylambda pdf',tukeylambda_pdf(tl_x(i),0.14_dp),tl_pdf(i), &
                2.0e-12_dp,2.0e-11_dp,failures)
            call check_close('tukeylambda logpdf',tukeylambda_logpdf(tl_x(i),0.14_dp), &
                tl_logpdf(i),2.0e-12_dp,2.0e-11_dp,failures)
            call check_close('tukeylambda cdf',tukeylambda_cdf(tl_x(i),0.14_dp),tl_cdf(i), &
                2.0e-12_dp,2.0e-11_dp,failures)
            call check_close('tukeylambda sf',tukeylambda_sf(tl_x(i),0.14_dp),tl_sf(i), &
                2.0e-12_dp,2.0e-11_dp,failures)
            call check_close('tukeylambda logcdf',tukeylambda_logcdf(tl_x(i),0.14_dp), &
                tl_logcdf(i),2.0e-12_dp,2.0e-11_dp,failures)
            call check_close('tukeylambda logsf',tukeylambda_logsf(tl_x(i),0.14_dp), &
                tl_logsf(i),2.0e-12_dp,2.0e-11_dp,failures)
            call check_close('tukeylambda ppf',tukeylambda_ppf(reference_probs(i),0.14_dp), &
                tl_ppf(i),2.0e-12_dp,2.0e-11_dp,failures)
            call check_close('tukeylambda isf',tukeylambda_isf(reference_probs(i),0.14_dp), &
                tl_isf(i),2.0e-12_dp,2.0e-11_dp,failures)

            call check_close('rice pdf',rice_pdf(ri_x(i),2.0_dp),ri_pdf(i), &
                2.0e-12_dp,2.0e-10_dp,failures)
            call check_close('rice logpdf',rice_logpdf(ri_x(i),2.0_dp),ri_logpdf(i), &
                2.0e-12_dp,2.0e-10_dp,failures)
            call check_close('rice cdf',rice_cdf(ri_x(i),2.0_dp),ri_cdf(i), &
                2.0e-12_dp,2.0e-9_dp,failures)
            call check_close('rice sf',rice_sf(ri_x(i),2.0_dp),ri_sf(i), &
                2.0e-12_dp,2.0e-9_dp,failures)
            call check_close('rice logcdf',rice_logcdf(ri_x(i),2.0_dp),ri_logcdf(i), &
                2.0e-11_dp,2.0e-9_dp,failures)
            call check_close('rice logsf',rice_logsf(ri_x(i),2.0_dp),ri_logsf(i), &
                2.0e-11_dp,2.0e-9_dp,failures)
            call check_close('rice ppf',rice_ppf(reference_probs(i),2.0_dp),ri_ppf(i), &
                2.0e-9_dp,2.0e-8_dp,failures)
            call check_close('rice isf',rice_isf(reference_probs(i),2.0_dp),ri_isf(i), &
                2.0e-9_dp,2.0e-8_dp,failures)
        end do
    end subroutine test_reference

    subroutine test_identities_endpoints(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: xvals(5) = [-2.0_dp,-0.5_dp,0.0_dp,0.5_dp,2.0_dp]
        real(dp), parameter :: pvals(5) = [1.0e-10_dp,1.0e-5_dp,0.2_dp,0.8_dp,1.0_dp-1.0e-10_dp]
        integer :: i

        do i = 1, size(xvals)
            call check_close('skewnorm a=0 pdf',skewnorm_pdf(xvals(i),0.0_dp), &
                normal_pdf(xvals(i)),2.0e-15_dp,2.0e-15_dp,failures)
            call check_close('skewnorm a=0 cdf',skewnorm_cdf(xvals(i),0.0_dp), &
                normal_cdf(xvals(i)),2.0e-15_dp,2.0e-15_dp,failures)
            call check_close('tukey lambda=0 pdf',tukeylambda_pdf(xvals(i),0.0_dp), &
                logistic_pdf(xvals(i)),2.0e-12_dp,2.0e-11_dp,failures)
            call check_close('tukey lambda=0 cdf',tukeylambda_cdf(xvals(i),0.0_dp), &
                logistic_cdf(xvals(i)),2.0e-12_dp,2.0e-11_dp,failures)
            if (xvals(i) > 0.0_dp) then
                call check_close('rice b=0 pdf',rice_pdf(xvals(i),0.0_dp), &
                    rayleigh_pdf(xvals(i)),2.0e-13_dp,2.0e-12_dp,failures)
                call check_close('rice b=0 cdf',rice_cdf(xvals(i),0.0_dp), &
                    rayleigh_cdf(xvals(i)),2.0e-13_dp,2.0e-12_dp,failures)
            end if
        end do
        call check_close('tukey lambda=1 center pdf',tukeylambda_pdf(0.0_dp,1.0_dp), &
            0.5_dp,2.0e-14_dp,2.0e-14_dp,failures)
        call check_close('tukey lambda=1 lower',tukeylambda_cdf(-1.0_dp,1.0_dp), &
            0.0_dp,0.0_dp,0.0_dp,failures)
        call check_close('tukey lambda=1 upper',tukeylambda_cdf(1.0_dp,1.0_dp), &
            1.0_dp,0.0_dp,0.0_dp,failures)

        do i = 1, size(pvals)
            call check_close('genexpon cdf-ppf', &
                genexpon_cdf(genexpon_ppf(pvals(i),1.3_dp,0.8_dp,1.7_dp),1.3_dp,0.8_dp,1.7_dp), &
                pvals(i),2.0e-11_dp,2.0e-8_dp,failures)
            call check_close('skewnorm cdf-ppf', &
                skewnorm_cdf(skewnorm_ppf(pvals(i),3.5_dp),3.5_dp), &
                pvals(i),2.0e-11_dp,2.0e-8_dp,failures)
            call check_close('tukeylambda sf-isf', &
                tukeylambda_sf(tukeylambda_isf(pvals(i),0.14_dp),0.14_dp), &
                pvals(i),2.0e-11_dp,2.0e-8_dp,failures)
            call check_close('rice sf-isf',rice_sf(rice_isf(pvals(i),2.0_dp),2.0_dp), &
                pvals(i),2.0e-10_dp,3.0e-8_dp,failures)
        end do

        call check_true('genexpon invalid shape',ieee_is_nan(genexpon_pdf(1.0_dp,0.0_dp,1.0_dp,1.0_dp)), &
            failures)
        call check_true('skewnorm invalid shape',ieee_is_nan(skewnorm_pdf(1.0_dp,ieee_value(1.0_dp,ieee_positive_inf))),failures)
        call check_true('tukeylambda invalid shape', &
            ieee_is_nan(tukeylambda_pdf(0.0_dp,ieee_value(1.0_dp,ieee_positive_inf))),failures)
        call check_true('rice invalid shape',ieee_is_nan(rice_pdf(1.0_dp,-1.0_dp)),failures)
        call check_true('rice upper quantile infinite', &
            .not. ieee_is_finite(rice_ppf(1.0_dp,2.0_dp)),failures)
    end subroutine test_identities_endpoints

    subroutine test_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: fd = 1.0e-6_dp, tol = 2.0e-3_dp
        real(dp), parameter :: gx_data(4) = [0.2_dp,0.5_dp,1.0_dp,1.8_dp]
        real(dp), parameter :: sn_data(4) = [-0.4_dp,0.2_dp,0.8_dp,1.5_dp]
        real(dp), parameter :: tl_data(4) = [-1.0_dp,-0.2_dp,0.5_dp,1.1_dp]
        real(dp), parameter :: ri_data(4) = [0.4_dp,1.2_dp,2.5_dp,3.4_dp]
        real(dp) :: s3(3), s5(5)

        s5 = genexpon_score(gx_data,1.3_dp,0.8_dp,1.7_dp,0.1_dp,1.2_dp)
        call check_close('genexpon score a',s5(1),centered( &
            genexpon_loglikelihood(gx_data,1.3_dp+fd,0.8_dp,1.7_dp,0.1_dp,1.2_dp), &
            genexpon_loglikelihood(gx_data,1.3_dp-fd,0.8_dp,1.7_dp,0.1_dp,1.2_dp),fd),tol,tol,failures)
        call check_close('genexpon score b',s5(2),centered( &
            genexpon_loglikelihood(gx_data,1.3_dp,0.8_dp+fd,1.7_dp,0.1_dp,1.2_dp), &
            genexpon_loglikelihood(gx_data,1.3_dp,0.8_dp-fd,1.7_dp,0.1_dp,1.2_dp),fd),tol,tol,failures)
        call check_close('genexpon score c',s5(3),centered( &
            genexpon_loglikelihood(gx_data,1.3_dp,0.8_dp,1.7_dp+fd,0.1_dp,1.2_dp), &
            genexpon_loglikelihood(gx_data,1.3_dp,0.8_dp,1.7_dp-fd,0.1_dp,1.2_dp),fd),tol,tol,failures)
        call check_close('genexpon score loc',s5(4),centered( &
            genexpon_loglikelihood(gx_data,1.3_dp,0.8_dp,1.7_dp,0.1_dp+fd,1.2_dp), &
            genexpon_loglikelihood(gx_data,1.3_dp,0.8_dp,1.7_dp,0.1_dp-fd,1.2_dp),fd),tol,tol,failures)
        call check_close('genexpon score scale',s5(5),centered( &
            genexpon_loglikelihood(gx_data,1.3_dp,0.8_dp,1.7_dp,0.1_dp,1.2_dp+fd), &
            genexpon_loglikelihood(gx_data,1.3_dp,0.8_dp,1.7_dp,0.1_dp,1.2_dp-fd),fd),tol,tol,failures)

        s3 = skewnorm_score(sn_data,2.2_dp,0.1_dp,1.3_dp)
        call check_score3('skewnorm',s3, &
            centered(skewnorm_loglikelihood(sn_data,2.2_dp+fd,0.1_dp,1.3_dp), &
                skewnorm_loglikelihood(sn_data,2.2_dp-fd,0.1_dp,1.3_dp),fd), &
            centered(skewnorm_loglikelihood(sn_data,2.2_dp,0.1_dp+fd,1.3_dp), &
                skewnorm_loglikelihood(sn_data,2.2_dp,0.1_dp-fd,1.3_dp),fd), &
            centered(skewnorm_loglikelihood(sn_data,2.2_dp,0.1_dp,1.3_dp+fd), &
                skewnorm_loglikelihood(sn_data,2.2_dp,0.1_dp,1.3_dp-fd),fd),tol,failures)

        s3 = tukeylambda_score(tl_data,0.14_dp,0.1_dp,1.2_dp)
        call check_score3('tukeylambda',s3, &
            centered(tukeylambda_loglikelihood(tl_data,0.14_dp+fd,0.1_dp,1.2_dp), &
                tukeylambda_loglikelihood(tl_data,0.14_dp-fd,0.1_dp,1.2_dp),fd), &
            centered(tukeylambda_loglikelihood(tl_data,0.14_dp,0.1_dp+fd,1.2_dp), &
                tukeylambda_loglikelihood(tl_data,0.14_dp,0.1_dp-fd,1.2_dp),fd), &
            centered(tukeylambda_loglikelihood(tl_data,0.14_dp,0.1_dp,1.2_dp+fd), &
                tukeylambda_loglikelihood(tl_data,0.14_dp,0.1_dp,1.2_dp-fd),fd),tol,failures)

        s3 = rice_score(ri_data,1.7_dp,0.1_dp,1.1_dp)
        call check_score3('rice',s3, &
            centered(rice_loglikelihood(ri_data,1.7_dp+fd,0.1_dp,1.1_dp), &
                rice_loglikelihood(ri_data,1.7_dp-fd,0.1_dp,1.1_dp),fd), &
            centered(rice_loglikelihood(ri_data,1.7_dp,0.1_dp+fd,1.1_dp), &
                rice_loglikelihood(ri_data,1.7_dp,0.1_dp-fd,1.1_dp),fd), &
            centered(rice_loglikelihood(ri_data,1.7_dp,0.1_dp,1.1_dp+fd), &
                rice_loglikelihood(ri_data,1.7_dp,0.1_dp,1.1_dp-fd),fd),tol,failures)
    end subroutine test_scores

    subroutine test_rng_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: gx_data(4) = [0.2_dp,0.5_dp,1.0_dp,1.8_dp]
        real(dp), parameter :: sn_data(4) = [-0.4_dp,0.2_dp,0.8_dp,1.5_dp]
        real(dp), parameter :: tl_data(4) = [-1.0_dp,-0.2_dp,0.5_dp,1.1_dp]
        real(dp), parameter :: ri_data(4) = [0.4_dp,1.2_dp,2.5_dp,3.4_dp]
        real(dp) :: actual, expected
        type(rng_state) :: reference, state
        type(fit_result) :: result

        call rng_seed(state,24680)
        call rng_seed(reference,24680)
        actual = genexpon_rvs(state,1.3_dp,0.8_dp,1.7_dp,0.1_dp,1.2_dp)
        expected = genexpon_ppf(rng_uniform(reference),1.3_dp,0.8_dp,1.7_dp,0.1_dp,1.2_dp)
        call check_close('genexpon rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = skewnorm_rvs(state,2.2_dp,0.1_dp,1.2_dp)
        expected = skewnorm_ppf(rng_uniform(reference),2.2_dp,0.1_dp,1.2_dp)
        call check_close('skewnorm rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = tukeylambda_rvs(state,0.14_dp,0.1_dp,1.2_dp)
        expected = tukeylambda_ppf(rng_uniform(reference),0.14_dp,0.1_dp,1.2_dp)
        call check_close('tukeylambda rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = rice_rvs(state,2.0_dp,0.1_dp,1.2_dp)
        expected = rice_ppf(rng_uniform(reference),2.0_dp,0.1_dp,1.2_dp)
        call check_close('rice rvs',actual,expected,0.0_dp,0.0_dp,failures)

        call rng_seed(state,13579)
        call rng_seed(reference,13579)
        actual = genexpon_rvs(state,0.0_dp,1.0_dp,1.0_dp)
        call check_true('invalid genexpon rvs NaN',ieee_is_nan(actual),failures)
        call check_close('invalid rvs preserves state',rng_uniform(state),rng_uniform(reference), &
            0.0_dp,0.0_dp,failures)

        call genexpon_fit(gx_data,[1.3_dp,0.8_dp,1.7_dp,0.1_dp,1.2_dp], &
            [1.3_dp,0.8_dp,1.7_dp,0.1_dp,1.2_dp],result)
        call check_true('genexpon fixed fit',result%success,failures)
        call skewnorm_fit(sn_data,[2.2_dp,0.1_dp,1.2_dp],[2.2_dp,0.1_dp,1.2_dp],result)
        call check_true('skewnorm fixed fit',result%success,failures)
        call tukeylambda_fit(tl_data,[0.14_dp,0.1_dp,1.2_dp],[0.14_dp,0.1_dp,1.2_dp],result)
        call check_true('tukeylambda fixed fit',result%success,failures)
        call rice_fit(ri_data,[2.0_dp,0.1_dp,1.2_dp],[2.0_dp,0.1_dp,1.2_dp],result)
        call check_true('rice fixed fit',result%success,failures)

        call check_close('C ABI genexpon pdf', &
            scifort_genexpon_pdf_f64(0.8_dp,1.3_dp,0.8_dp,1.7_dp,0.0_dp,1.0_dp), &
            genexpon_pdf(0.8_dp,1.3_dp,0.8_dp,1.7_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI genexpon cdf', &
            scifort_genexpon_cdf_f64(0.8_dp,1.3_dp,0.8_dp,1.7_dp,0.0_dp,1.0_dp), &
            genexpon_cdf(0.8_dp,1.3_dp,0.8_dp,1.7_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI genexpon ppf', &
            scifort_genexpon_ppf_f64(0.9_dp,1.3_dp,0.8_dp,1.7_dp,0.0_dp,1.0_dp), &
            genexpon_ppf(0.9_dp,1.3_dp,0.8_dp,1.7_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI skewnorm pdf',scifort_skewnorm_pdf_f64(0.5_dp,3.5_dp,0.0_dp,1.0_dp), &
            skewnorm_pdf(0.5_dp,3.5_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI skewnorm cdf',scifort_skewnorm_cdf_f64(0.5_dp,3.5_dp,0.0_dp,1.0_dp), &
            skewnorm_cdf(0.5_dp,3.5_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI skewnorm ppf',scifort_skewnorm_ppf_f64(0.9_dp,3.5_dp,0.0_dp,1.0_dp), &
            skewnorm_ppf(0.9_dp,3.5_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI tukeylambda pdf', &
            scifort_tukeylambda_pdf_f64(0.5_dp,0.14_dp,0.0_dp,1.0_dp), &
            tukeylambda_pdf(0.5_dp,0.14_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI tukeylambda cdf', &
            scifort_tukeylambda_cdf_f64(0.5_dp,0.14_dp,0.0_dp,1.0_dp), &
            tukeylambda_cdf(0.5_dp,0.14_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI tukeylambda ppf', &
            scifort_tukeylambda_ppf_f64(0.9_dp,0.14_dp,0.0_dp,1.0_dp), &
            tukeylambda_ppf(0.9_dp,0.14_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI rice pdf',scifort_rice_pdf_f64(2.0_dp,2.0_dp,0.0_dp,1.0_dp), &
            rice_pdf(2.0_dp,2.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI rice cdf',scifort_rice_cdf_f64(2.0_dp,2.0_dp,0.0_dp,1.0_dp), &
            rice_cdf(2.0_dp,2.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI rice ppf',scifort_rice_ppf_f64(0.9_dp,2.0_dp,0.0_dp,1.0_dp), &
            rice_ppf(0.9_dp,2.0_dp),0.0_dp,0.0_dp,failures)
    end subroutine test_rng_fit_c_api

    subroutine check_score3(label, actual, e1, e2, e3, tol, failures)
        character(len=*), intent(in) :: label !! distribution label for diagnostics
        real(dp), intent(in) :: actual(3) !! analytic score vector
        real(dp), intent(in) :: e1 !! finite-difference first component
        real(dp), intent(in) :: e2 !! finite-difference second component
        real(dp), intent(in) :: e3 !! finite-difference third component
        real(dp), intent(in) :: tol !! absolute and relative tolerance
        integer, intent(inout) :: failures !! running count of failed checks
        call check_close(trim(label)//' score 1',actual(1),e1,tol,tol,failures)
        call check_close(trim(label)//' score 2',actual(2),e2,tol,tol,failures)
        call check_close(trim(label)//' score 3',actual(3),e3,tol,tol,failures)
    end subroutine check_score3

    pure function centered(plus, minus, h) result(value)
        real(dp), intent(in) :: plus !! objective at theta+h
        real(dp), intent(in) :: minus !! objective at theta-h
        real(dp), intent(in) :: h !! positive finite-difference step
        real(dp) :: value
        value = (plus - minus) / (2.0_dp * h)
    end function centered

end program test_genexpon_skewnorm_tukeylambda_rice
