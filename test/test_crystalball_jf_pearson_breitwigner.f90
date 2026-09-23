! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_crystalball_jf_pearson_breitwigner
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_positive_inf, ieee_value
    use scifort_c_api
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state, rng_uniform
    use scifort_stats
    use test_crystalball_jf_pearson_breitwigner_reference
    use test_checks, only : check_close, check_true
    implicit none

    integer :: failures

    failures = 0
    call test_reference(failures)
    call test_identities_endpoints(failures)
    call test_scores(failures)
    call test_rng_fit_c_api(failures)

    if (failures == 0) then
        print '(a)', 'test_crystalball_jf_pearson_breitwigner: PASS'
    else
        print '(a,1x,i0)', 'test_crystalball_jf_pearson_breitwigner: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        integer :: i

        do i = 1, size(reference_probs)
            call check_close('crystalball pdf',crystalball_pdf(cb_x(i),1.7_dp,3.5_dp),cb_pdf(i), &
                2.0e-12_dp,2.0e-11_dp,failures)
            call check_close('crystalball logpdf',crystalball_logpdf(cb_x(i),1.7_dp,3.5_dp),cb_logpdf(i), &
                2.0e-12_dp,2.0e-11_dp,failures)
            call check_close('crystalball cdf',crystalball_cdf(cb_x(i),1.7_dp,3.5_dp),cb_cdf(i), &
                2.0e-12_dp,2.0e-11_dp,failures)
            call check_close('crystalball sf',crystalball_sf(cb_x(i),1.7_dp,3.5_dp),cb_sf(i), &
                2.0e-12_dp,2.0e-11_dp,failures)
            call check_close('crystalball logcdf',crystalball_logcdf(cb_x(i),1.7_dp,3.5_dp),cb_logcdf(i), &
                2.0e-12_dp,2.0e-11_dp,failures)
            call check_close('crystalball logsf',crystalball_logsf(cb_x(i),1.7_dp,3.5_dp),cb_logsf(i), &
                2.0e-12_dp,2.0e-11_dp,failures)
            call check_close('crystalball ppf',crystalball_ppf(reference_probs(i),1.7_dp,3.5_dp),cb_ppf(i), &
                2.0e-9_dp,2.0e-10_dp,failures)
            call check_close('crystalball isf',crystalball_isf(reference_probs(i),1.7_dp,3.5_dp),cb_isf(i), &
                2.0e-9_dp,2.0e-10_dp,failures)

            call check_close('jf skew t pdf',jf_skew_t_pdf(jf_x(i),2.3_dp,4.1_dp),jf_pdf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('jf skew t logpdf',jf_skew_t_logpdf(jf_x(i),2.3_dp,4.1_dp),jf_logpdf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('jf skew t cdf',jf_skew_t_cdf(jf_x(i),2.3_dp,4.1_dp),jf_cdf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('jf skew t sf',jf_skew_t_sf(jf_x(i),2.3_dp,4.1_dp),jf_sf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('jf skew t logcdf',jf_skew_t_logcdf(jf_x(i),2.3_dp,4.1_dp),jf_logcdf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('jf skew t logsf',jf_skew_t_logsf(jf_x(i),2.3_dp,4.1_dp),jf_logsf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('jf skew t ppf',jf_skew_t_ppf(reference_probs(i),2.3_dp,4.1_dp),jf_ppf(i), &
                5.0e-9_dp,5.0e-9_dp,failures)
            call check_close('jf skew t isf',jf_skew_t_isf(reference_probs(i),2.3_dp,4.1_dp),jf_isf(i), &
                5.0e-9_dp,5.0e-9_dp,failures)

            call check_close('pearson3 pdf',pearson3_pdf(p3_x(i),-1.2_dp),p3_pdf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('pearson3 logpdf',pearson3_logpdf(p3_x(i),-1.2_dp),p3_logpdf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('pearson3 cdf',pearson3_cdf(p3_x(i),-1.2_dp),p3_cdf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('pearson3 sf',pearson3_sf(p3_x(i),-1.2_dp),p3_sf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('pearson3 logcdf',pearson3_logcdf(p3_x(i),-1.2_dp),p3_logcdf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('pearson3 logsf',pearson3_logsf(p3_x(i),-1.2_dp),p3_logsf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('pearson3 ppf',pearson3_ppf(reference_probs(i),-1.2_dp),p3_ppf(i), &
                5.0e-9_dp,5.0e-9_dp,failures)
            call check_close('pearson3 isf',pearson3_isf(reference_probs(i),-1.2_dp),p3_isf(i), &
                5.0e-9_dp,5.0e-9_dp,failures)

            call check_close('rel breitwigner pdf',rel_breitwigner_pdf(rb_x(i),5.0_dp),rb_pdf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('rel breitwigner logpdf',rel_breitwigner_logpdf(rb_x(i),5.0_dp),rb_logpdf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('rel breitwigner cdf',rel_breitwigner_cdf(rb_x(i),5.0_dp),rb_cdf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('rel breitwigner sf',rel_breitwigner_sf(rb_x(i),5.0_dp),rb_sf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('rel breitwigner logcdf',rel_breitwigner_logcdf(rb_x(i),5.0_dp),rb_logcdf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('rel breitwigner logsf',rel_breitwigner_logsf(rb_x(i),5.0_dp),rb_logsf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('rel breitwigner ppf',rel_breitwigner_ppf(reference_probs(i),5.0_dp),rb_ppf(i), &
                2.0e-8_dp,2.0e-9_dp,failures)
            call check_close('rel breitwigner isf',rel_breitwigner_isf(reference_probs(i),5.0_dp),rb_isf(i), &
                2.0e-8_dp,2.0e-9_dp,failures)
        end do
    end subroutine test_reference

    subroutine test_identities_endpoints(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: xvals(5) = [-3.0_dp,-1.0_dp,0.0_dp,1.0_dp,3.0_dp]
        real(dp), parameter :: pvals(5) = [1.0e-12_dp,1.0e-6_dp,0.2_dp,0.8_dp,1.0_dp-1.0e-12_dp]
        real(dp) :: boundary, pinf
        integer :: i

        do i = 1, size(xvals)
            call check_close('jf symmetric t pdf',jf_skew_t_pdf(xvals(i),2.3_dp,2.3_dp), &
                t_pdf(xvals(i),4.6_dp),3.0e-13_dp,3.0e-12_dp,failures)
            call check_close('jf symmetric t cdf',jf_skew_t_cdf(xvals(i),2.3_dp,2.3_dp), &
                t_cdf(xvals(i),4.6_dp),3.0e-13_dp,3.0e-12_dp,failures)
            call check_close('pearson3 zero skew pdf',pearson3_pdf(xvals(i),0.0_dp), &
                normal_pdf(xvals(i)),3.0e-15_dp,3.0e-15_dp,failures)
            call check_close('pearson3 zero skew cdf',pearson3_cdf(xvals(i),0.0_dp), &
                normal_cdf(xvals(i)),0.0_dp,0.0_dp,failures)
            call check_close('pearson3 reflection pdf',pearson3_pdf(xvals(i),-1.2_dp), &
                pearson3_pdf(-xvals(i),1.2_dp),3.0e-13_dp,3.0e-12_dp,failures)
            call check_close('pearson3 reflection cdf',pearson3_cdf(xvals(i),-1.2_dp), &
                pearson3_sf(-xvals(i),1.2_dp),3.0e-13_dp,3.0e-12_dp,failures)
        end do

        do i = 1, size(pvals)
            call check_close('crystalball cdf-ppf', &
                crystalball_cdf(crystalball_ppf(pvals(i),1.7_dp,3.5_dp),1.7_dp,3.5_dp), &
                pvals(i),5.0e-12_dp,2.0e-8_dp,failures)
            call check_close('crystalball sf-isf', &
                crystalball_sf(crystalball_isf(pvals(i),1.7_dp,3.5_dp),1.7_dp,3.5_dp), &
                pvals(i),5.0e-12_dp,2.0e-8_dp,failures)
            call check_close('jf skew t cdf-ppf', &
                jf_skew_t_cdf(jf_skew_t_ppf(pvals(i),2.3_dp,4.1_dp),2.3_dp,4.1_dp), &
                pvals(i),5.0e-12_dp,2.0e-8_dp,failures)
            call check_close('jf skew t sf-isf', &
                jf_skew_t_sf(jf_skew_t_isf(pvals(i),2.3_dp,4.1_dp),2.3_dp,4.1_dp), &
                pvals(i),5.0e-12_dp,2.0e-8_dp,failures)
            call check_close('pearson3 cdf-ppf', &
                pearson3_cdf(pearson3_ppf(pvals(i),1.2_dp),1.2_dp), &
                pvals(i),5.0e-12_dp,2.0e-8_dp,failures)
            call check_close('pearson3 sf-isf', &
                pearson3_sf(pearson3_isf(pvals(i),-1.2_dp),-1.2_dp), &
                pvals(i),5.0e-12_dp,2.0e-8_dp,failures)
            call check_close('rel breitwigner cdf-ppf', &
                rel_breitwigner_cdf(rel_breitwigner_ppf(pvals(i),5.0_dp),5.0_dp), &
                pvals(i),2.0e-11_dp,2.0e-7_dp,failures)
            call check_close('rel breitwigner sf-isf', &
                rel_breitwigner_sf(rel_breitwigner_isf(pvals(i),5.0_dp),5.0_dp), &
                pvals(i),2.0e-11_dp,2.0e-7_dp,failures)
        end do

        boundary = -1.7_dp
        call check_close('crystalball breakpoint pdf continuity', &
            crystalball_pdf(boundary-1.0e-9_dp,1.7_dp,3.5_dp), &
            crystalball_pdf(boundary+1.0e-9_dp,1.7_dp,3.5_dp),1.0e-9_dp,1.0e-8_dp,failures)
        call check_close('crystalball cdf+sf',crystalball_cdf(0.3_dp,1.7_dp,3.5_dp)+ &
            crystalball_sf(0.3_dp,1.7_dp,3.5_dp),1.0_dp,2.0e-15_dp,2.0e-15_dp,failures)
        call check_close('rel breitwigner lower cdf',rel_breitwigner_cdf(0.0_dp,5.0_dp), &
            0.0_dp,0.0_dp,0.0_dp,failures)
        call check_close('rel breitwigner lower sf',rel_breitwigner_sf(0.0_dp,5.0_dp), &
            1.0_dp,0.0_dp,0.0_dp,failures)
        call check_close('rel breitwigner cdf+sf',rel_breitwigner_cdf(7.0_dp,5.0_dp)+ &
            rel_breitwigner_sf(7.0_dp,5.0_dp),1.0_dp,2.0e-14_dp,2.0e-14_dp,failures)
        call check_true('jf extreme right logsf finite', &
            ieee_is_finite(jf_skew_t_logsf(1.0e150_dp,2.3_dp,4.1_dp)) .and. &
            jf_skew_t_logsf(1.0e150_dp,2.3_dp,4.1_dp) < -1000.0_dp,failures)
        call check_true('jf extreme left logcdf finite', &
            ieee_is_finite(jf_skew_t_logcdf(-1.0e150_dp,2.3_dp,4.1_dp)) .and. &
            jf_skew_t_logcdf(-1.0e150_dp,2.3_dp,4.1_dp) < -1000.0_dp,failures)

        call check_true('crystalball invalid beta',ieee_is_nan(crystalball_pdf(0.0_dp,0.0_dp,3.5_dp)),failures)
        call check_true('crystalball invalid m',ieee_is_nan(crystalball_pdf(0.0_dp,1.7_dp,1.0_dp)),failures)
        call check_true('jf skew t invalid a',ieee_is_nan(jf_skew_t_pdf(0.0_dp,0.0_dp,4.1_dp)),failures)
        call check_true('jf skew t invalid b',ieee_is_nan(jf_skew_t_pdf(0.0_dp,2.3_dp,-1.0_dp)),failures)
        pinf = ieee_value(0.0_dp,ieee_positive_inf)
        call check_true('pearson3 invalid skew',ieee_is_nan(pearson3_pdf(0.0_dp,pinf)),failures)
        call check_true('rel breitwigner invalid rho',ieee_is_nan(rel_breitwigner_pdf(1.0_dp,0.0_dp)),failures)
        call check_true('rel breitwigner upper quantile infinite', &
            .not. ieee_is_finite(rel_breitwigner_ppf(1.0_dp,5.0_dp)),failures)
        call check_true('rel breitwigner endpoint score undefined', &
            all(ieee_is_nan(rel_breitwigner_score([0.0_dp,1.0_dp],5.0_dp))),failures)
    end subroutine test_identities_endpoints

    subroutine test_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: fd = 1.0e-6_dp, tol = 1.5e-3_dp
        real(dp), parameter :: cb_data(4) = [-2.5_dp,-1.0_dp,0.4_dp,2.0_dp]
        real(dp), parameter :: jf_data(4) = [-2.0_dp,-0.5_dp,0.5_dp,2.0_dp]
        real(dp), parameter :: p3_data(4) = [-1.0_dp,0.2_dp,1.0_dp,2.0_dp]
        real(dp), parameter :: rb_data(4) = [0.2_dp,1.5_dp,4.5_dp,8.0_dp]
        real(dp) :: s3(3), s4(4)

        s4 = crystalball_score(cb_data,1.7_dp,3.5_dp,0.1_dp,1.2_dp)
        call check_close('crystalball score beta',s4(1),centered( &
            crystalball_loglikelihood(cb_data,1.7_dp+fd,3.5_dp,0.1_dp,1.2_dp), &
            crystalball_loglikelihood(cb_data,1.7_dp-fd,3.5_dp,0.1_dp,1.2_dp),fd),tol,tol,failures)
        call check_close('crystalball score m',s4(2),centered( &
            crystalball_loglikelihood(cb_data,1.7_dp,3.5_dp+fd,0.1_dp,1.2_dp), &
            crystalball_loglikelihood(cb_data,1.7_dp,3.5_dp-fd,0.1_dp,1.2_dp),fd),tol,tol,failures)
        call check_close('crystalball score loc',s4(3),centered( &
            crystalball_loglikelihood(cb_data,1.7_dp,3.5_dp,0.1_dp+fd,1.2_dp), &
            crystalball_loglikelihood(cb_data,1.7_dp,3.5_dp,0.1_dp-fd,1.2_dp),fd),tol,tol,failures)
        call check_close('crystalball score scale',s4(4),centered( &
            crystalball_loglikelihood(cb_data,1.7_dp,3.5_dp,0.1_dp,1.2_dp+fd), &
            crystalball_loglikelihood(cb_data,1.7_dp,3.5_dp,0.1_dp,1.2_dp-fd),fd),tol,tol,failures)

        s4 = jf_skew_t_score(jf_data,2.3_dp,4.1_dp,0.1_dp,1.2_dp)
        call check_close('jf skew t score a',s4(1),centered( &
            jf_skew_t_loglikelihood(jf_data,2.3_dp+fd,4.1_dp,0.1_dp,1.2_dp), &
            jf_skew_t_loglikelihood(jf_data,2.3_dp-fd,4.1_dp,0.1_dp,1.2_dp),fd),tol,tol,failures)
        call check_close('jf skew t score b',s4(2),centered( &
            jf_skew_t_loglikelihood(jf_data,2.3_dp,4.1_dp+fd,0.1_dp,1.2_dp), &
            jf_skew_t_loglikelihood(jf_data,2.3_dp,4.1_dp-fd,0.1_dp,1.2_dp),fd),tol,tol,failures)
        call check_close('jf skew t score loc',s4(3),centered( &
            jf_skew_t_loglikelihood(jf_data,2.3_dp,4.1_dp,0.1_dp+fd,1.2_dp), &
            jf_skew_t_loglikelihood(jf_data,2.3_dp,4.1_dp,0.1_dp-fd,1.2_dp),fd),tol,tol,failures)
        call check_close('jf skew t score scale',s4(4),centered( &
            jf_skew_t_loglikelihood(jf_data,2.3_dp,4.1_dp,0.1_dp,1.2_dp+fd), &
            jf_skew_t_loglikelihood(jf_data,2.3_dp,4.1_dp,0.1_dp,1.2_dp-fd),fd),tol,tol,failures)

        s3 = pearson3_score(p3_data,1.2_dp,0.1_dp,1.2_dp)
        call check_close('pearson3 score skew',s3(1),centered( &
            pearson3_loglikelihood(p3_data,1.2_dp+fd,0.1_dp,1.2_dp), &
            pearson3_loglikelihood(p3_data,1.2_dp-fd,0.1_dp,1.2_dp),fd),tol,tol,failures)
        call check_close('pearson3 score loc',s3(2),centered( &
            pearson3_loglikelihood(p3_data,1.2_dp,0.1_dp+fd,1.2_dp), &
            pearson3_loglikelihood(p3_data,1.2_dp,0.1_dp-fd,1.2_dp),fd),tol,tol,failures)
        call check_close('pearson3 score scale',s3(3),centered( &
            pearson3_loglikelihood(p3_data,1.2_dp,0.1_dp,1.2_dp+fd), &
            pearson3_loglikelihood(p3_data,1.2_dp,0.1_dp,1.2_dp-fd),fd),tol,tol,failures)

        s3 = rel_breitwigner_score(rb_data,5.0_dp,0.1_dp,1.2_dp)
        call check_close('rel breitwigner score rho',s3(1),centered( &
            rel_breitwigner_loglikelihood(rb_data,5.0_dp+fd,0.1_dp,1.2_dp), &
            rel_breitwigner_loglikelihood(rb_data,5.0_dp-fd,0.1_dp,1.2_dp),fd),tol,tol,failures)
        call check_close('rel breitwigner score loc',s3(2),centered( &
            rel_breitwigner_loglikelihood(rb_data,5.0_dp,0.1_dp+fd,1.2_dp), &
            rel_breitwigner_loglikelihood(rb_data,5.0_dp,0.1_dp-fd,1.2_dp),fd),tol,tol,failures)
        call check_close('rel breitwigner score scale',s3(3),centered( &
            rel_breitwigner_loglikelihood(rb_data,5.0_dp,0.1_dp,1.2_dp+fd), &
            rel_breitwigner_loglikelihood(rb_data,5.0_dp,0.1_dp,1.2_dp-fd),fd),tol,tol,failures)
    end subroutine test_scores

    subroutine test_rng_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: cb_data(4) = [-2.5_dp,-1.0_dp,0.4_dp,2.0_dp]
        real(dp), parameter :: jf_data(4) = [-2.0_dp,-0.5_dp,0.5_dp,2.0_dp]
        real(dp), parameter :: p3_data(4) = [-1.0_dp,0.2_dp,1.0_dp,2.0_dp]
        real(dp), parameter :: rb_data(4) = [0.2_dp,1.5_dp,4.5_dp,8.0_dp]
        real(dp) :: actual, expected
        type(rng_state) :: reference, state
        type(fit_result) :: result

        call rng_seed(state,24680)
        call rng_seed(reference,24680)
        actual = crystalball_rvs(state,1.7_dp,3.5_dp,0.1_dp,1.2_dp)
        expected = crystalball_ppf(rng_uniform(reference),1.7_dp,3.5_dp,0.1_dp,1.2_dp)
        call check_close('crystalball rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = jf_skew_t_rvs(state,2.3_dp,4.1_dp,0.1_dp,1.2_dp)
        expected = jf_skew_t_ppf(rng_uniform(reference),2.3_dp,4.1_dp,0.1_dp,1.2_dp)
        call check_close('jf skew t rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = pearson3_rvs(state,1.2_dp,0.1_dp,1.2_dp)
        expected = pearson3_ppf(rng_uniform(reference),1.2_dp,0.1_dp,1.2_dp)
        call check_close('pearson3 rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = rel_breitwigner_rvs(state,5.0_dp,0.1_dp,1.2_dp)
        expected = rel_breitwigner_ppf(rng_uniform(reference),5.0_dp,0.1_dp,1.2_dp)
        call check_close('rel breitwigner rvs',actual,expected,0.0_dp,0.0_dp,failures)

        call rng_seed(state,13579)
        call rng_seed(reference,13579)
        actual = crystalball_rvs(state,0.0_dp,3.5_dp)
        call check_true('invalid crystalball rvs NaN',ieee_is_nan(actual),failures)
        call check_close('invalid rvs preserves state',rng_uniform(state),rng_uniform(reference), &
            0.0_dp,0.0_dp,failures)

        call crystalball_fit(cb_data,[1.7_dp,3.5_dp,0.1_dp,1.2_dp], &
            [1.7_dp,3.5_dp,0.1_dp,1.2_dp],result)
        call check_true('crystalball fixed fit',result%success,failures)
        call jf_skew_t_fit(jf_data,[2.3_dp,4.1_dp,0.1_dp,1.2_dp], &
            [2.3_dp,4.1_dp,0.1_dp,1.2_dp],result)
        call check_true('jf skew t fixed fit',result%success,failures)
        call pearson3_fit(p3_data,[1.2_dp,0.1_dp,1.2_dp],[1.2_dp,0.1_dp,1.2_dp],result)
        call check_true('pearson3 fixed fit',result%success,failures)
        call rel_breitwigner_fit(rb_data,[5.0_dp,0.1_dp,1.2_dp],[5.0_dp,0.1_dp,1.2_dp],result)
        call check_true('rel breitwigner fixed fit',result%success,failures)

        call check_close('C ABI crystalball pdf', &
            scifort_crystalball_pdf_f64(-1.0_dp,1.7_dp,3.5_dp,0.0_dp,1.0_dp), &
            crystalball_pdf(-1.0_dp,1.7_dp,3.5_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI crystalball cdf', &
            scifort_crystalball_cdf_f64(0.5_dp,1.7_dp,3.5_dp,0.0_dp,1.0_dp), &
            crystalball_cdf(0.5_dp,1.7_dp,3.5_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI crystalball ppf', &
            scifort_crystalball_ppf_f64(0.9_dp,1.7_dp,3.5_dp,0.0_dp,1.0_dp), &
            crystalball_ppf(0.9_dp,1.7_dp,3.5_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI jf skew t pdf', &
            scifort_jf_skew_t_pdf_f64(-1.0_dp,2.3_dp,4.1_dp,0.0_dp,1.0_dp), &
            jf_skew_t_pdf(-1.0_dp,2.3_dp,4.1_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI jf skew t cdf', &
            scifort_jf_skew_t_cdf_f64(0.5_dp,2.3_dp,4.1_dp,0.0_dp,1.0_dp), &
            jf_skew_t_cdf(0.5_dp,2.3_dp,4.1_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI jf skew t ppf', &
            scifort_jf_skew_t_ppf_f64(0.9_dp,2.3_dp,4.1_dp,0.0_dp,1.0_dp), &
            jf_skew_t_ppf(0.9_dp,2.3_dp,4.1_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI pearson3 pdf', &
            scifort_pearson3_pdf_f64(0.2_dp,-1.2_dp,0.0_dp,1.0_dp), &
            pearson3_pdf(0.2_dp,-1.2_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI pearson3 cdf', &
            scifort_pearson3_cdf_f64(0.2_dp,-1.2_dp,0.0_dp,1.0_dp), &
            pearson3_cdf(0.2_dp,-1.2_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI pearson3 ppf', &
            scifort_pearson3_ppf_f64(0.9_dp,-1.2_dp,0.0_dp,1.0_dp), &
            pearson3_ppf(0.9_dp,-1.2_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI rel breitwigner pdf', &
            scifort_rel_breitwigner_pdf_f64(5.0_dp,5.0_dp,0.0_dp,1.0_dp), &
            rel_breitwigner_pdf(5.0_dp,5.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI rel breitwigner cdf', &
            scifort_rel_breitwigner_cdf_f64(5.0_dp,5.0_dp,0.0_dp,1.0_dp), &
            rel_breitwigner_cdf(5.0_dp,5.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI rel breitwigner ppf', &
            scifort_rel_breitwigner_ppf_f64(0.9_dp,5.0_dp,0.0_dp,1.0_dp), &
            rel_breitwigner_ppf(0.9_dp,5.0_dp),0.0_dp,0.0_dp,failures)
    end subroutine test_rng_fit_c_api

    pure function centered(plus, minus, h) result(value)
        real(dp), intent(in) :: plus !! objective at theta+h
        real(dp), intent(in) :: minus !! objective at theta-h
        real(dp), intent(in) :: h !! positive finite-difference step
        real(dp) :: value
        value = (plus - minus) / (2.0_dp * h)
    end function centered

end program test_crystalball_jf_pearson_breitwigner
