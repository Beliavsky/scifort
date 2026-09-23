! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_gengamma_halfgennorm_argus_erlang
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_c_api
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state, rng_uniform
    use scifort_stats
    use test_gengamma_halfgennorm_argus_erlang_reference
    use test_checks, only : check_close, check_true
    implicit none

    integer :: failures

    failures = 0
    call test_reference(failures)
    call test_identities_endpoints(failures)
    call test_scores(failures)
    call test_rng_fit_c_api(failures)

    if (failures == 0) then
        print '(a)', 'test_gengamma_halfgennorm_argus_erlang: PASS'
    else
        print '(a,1x,i0)', 'test_gengamma_halfgennorm_argus_erlang: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        integer :: i

        do i = 1, size(reference_probs)
            call check_close('gengamma pdf',gengamma_pdf(gg_x(i),2.3_dp,-1.3_dp),gg_pdf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('gengamma logpdf',gengamma_logpdf(gg_x(i),2.3_dp,-1.3_dp), &
                gg_logpdf(i),2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('gengamma cdf',gengamma_cdf(gg_x(i),2.3_dp,-1.3_dp),gg_cdf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('gengamma sf',gengamma_sf(gg_x(i),2.3_dp,-1.3_dp),gg_sf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('gengamma logcdf',gengamma_logcdf(gg_x(i),2.3_dp,-1.3_dp), &
                gg_logcdf(i),2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('gengamma logsf',gengamma_logsf(gg_x(i),2.3_dp,-1.3_dp), &
                gg_logsf(i),2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('gengamma ppf',gengamma_ppf(reference_probs(i),2.3_dp,-1.3_dp), &
                gg_ppf(i),5.0e-9_dp,5.0e-9_dp,failures)
            call check_close('gengamma isf',gengamma_isf(reference_probs(i),2.3_dp,-1.3_dp), &
                gg_isf(i),5.0e-9_dp,5.0e-9_dp,failures)

            call check_close('halfgennorm pdf',halfgennorm_pdf(hg_x(i),1.7_dp),hg_pdf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('halfgennorm logpdf',halfgennorm_logpdf(hg_x(i),1.7_dp), &
                hg_logpdf(i),2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('halfgennorm cdf',halfgennorm_cdf(hg_x(i),1.7_dp),hg_cdf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('halfgennorm sf',halfgennorm_sf(hg_x(i),1.7_dp),hg_sf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('halfgennorm logcdf',halfgennorm_logcdf(hg_x(i),1.7_dp), &
                hg_logcdf(i),2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('halfgennorm logsf',halfgennorm_logsf(hg_x(i),1.7_dp), &
                hg_logsf(i),2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('halfgennorm ppf',halfgennorm_ppf(reference_probs(i),1.7_dp), &
                hg_ppf(i),5.0e-9_dp,5.0e-9_dp,failures)
            call check_close('halfgennorm isf',halfgennorm_isf(reference_probs(i),1.7_dp), &
                hg_isf(i),5.0e-9_dp,5.0e-9_dp,failures)

            call check_close('argus pdf',argus_pdf(ar_x(i),2.3_dp),ar_pdf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('argus logpdf',argus_logpdf(ar_x(i),2.3_dp),ar_logpdf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('argus cdf',argus_cdf(ar_x(i),2.3_dp),ar_cdf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('argus sf',argus_sf(ar_x(i),2.3_dp),ar_sf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('argus logcdf',argus_logcdf(ar_x(i),2.3_dp),ar_logcdf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('argus logsf',argus_logsf(ar_x(i),2.3_dp),ar_logsf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('argus ppf',argus_ppf(reference_probs(i),2.3_dp),ar_ppf(i), &
                2.0e-9_dp,2.0e-9_dp,failures)
            call check_close('argus isf',argus_isf(reference_probs(i),2.3_dp),ar_isf(i), &
                2.0e-9_dp,2.0e-9_dp,failures)

            call check_close('erlang pdf',erlang_pdf(er_x(i),4.0_dp),er_pdf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('erlang logpdf',erlang_logpdf(er_x(i),4.0_dp),er_logpdf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('erlang cdf',erlang_cdf(er_x(i),4.0_dp),er_cdf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('erlang sf',erlang_sf(er_x(i),4.0_dp),er_sf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('erlang logcdf',erlang_logcdf(er_x(i),4.0_dp),er_logcdf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('erlang logsf',erlang_logsf(er_x(i),4.0_dp),er_logsf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('erlang ppf',erlang_ppf(reference_probs(i),4.0_dp),er_ppf(i), &
                5.0e-9_dp,5.0e-9_dp,failures)
            call check_close('erlang isf',erlang_isf(reference_probs(i),4.0_dp),er_isf(i), &
                5.0e-9_dp,5.0e-9_dp,failures)
        end do
    end subroutine test_reference

    subroutine test_identities_endpoints(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: xvals(5) = [0.1_dp,0.3_dp,0.8_dp,1.7_dp,4.0_dp]
        real(dp), parameter :: pvals(5) = [1.0e-12_dp,1.0e-6_dp,0.2_dp,0.8_dp,1.0_dp-1.0e-12_dp]
        real(dp), parameter :: invsqrt2 = 0.707106781186547524400844362104849039_dp
        real(dp) :: a3(3)
        integer :: i

        do i = 1, size(xvals)
            call check_close('gengamma-gamma pdf',gengamma_pdf(xvals(i),2.3_dp,1.0_dp), &
                gamma_pdf(xvals(i),2.3_dp),3.0e-14_dp,3.0e-13_dp,failures)
            call check_close('gengamma-weibull cdf',gengamma_cdf(xvals(i),1.0_dp,1.7_dp), &
                weibull_cdf(xvals(i),1.7_dp),3.0e-14_dp,3.0e-13_dp,failures)
            call check_close('halfgennorm-exponential pdf',halfgennorm_pdf(xvals(i),1.0_dp), &
                exponential_pdf(xvals(i)),3.0e-14_dp,3.0e-13_dp,failures)
            call check_close('halfgennorm-halfnormal cdf',halfgennorm_cdf(xvals(i),2.0_dp), &
                halfnorm_cdf(xvals(i),0.0_dp,invsqrt2),3.0e-13_dp,3.0e-12_dp,failures)
            call check_close('erlang-gamma pdf integer',erlang_pdf(xvals(i),4.0_dp), &
                gamma_pdf(xvals(i),4.0_dp),0.0_dp,0.0_dp,failures)
            call check_close('erlang-gamma cdf noninteger',erlang_cdf(xvals(i),2.5_dp), &
                gamma_cdf(xvals(i),2.5_dp),0.0_dp,0.0_dp,failures)
        end do

        do i = 1, size(pvals)
            call check_close('gengamma cdf-ppf positive c', &
                gengamma_cdf(gengamma_ppf(pvals(i),2.3_dp,1.3_dp),2.3_dp,1.3_dp), &
                pvals(i),5.0e-12_dp,5.0e-9_dp,failures)
            call check_close('gengamma sf-isf negative c', &
                gengamma_sf(gengamma_isf(pvals(i),2.3_dp,-1.3_dp),2.3_dp,-1.3_dp), &
                pvals(i),5.0e-12_dp,5.0e-9_dp,failures)
            call check_close('halfgennorm cdf-ppf', &
                halfgennorm_cdf(halfgennorm_ppf(pvals(i),1.7_dp),1.7_dp), &
                pvals(i),5.0e-12_dp,5.0e-9_dp,failures)
            call check_close('halfgennorm sf-isf', &
                halfgennorm_sf(halfgennorm_isf(pvals(i),1.7_dp),1.7_dp), &
                pvals(i),5.0e-12_dp,5.0e-9_dp,failures)
            call check_close('argus cdf-ppf',argus_cdf(argus_ppf(pvals(i),2.3_dp),2.3_dp), &
                pvals(i),3.0e-11_dp,3.0e-8_dp,failures)
            call check_close('argus sf-isf',argus_sf(argus_isf(pvals(i),2.3_dp),2.3_dp), &
                pvals(i),3.0e-11_dp,3.0e-8_dp,failures)
            call check_close('erlang cdf-ppf',erlang_cdf(erlang_ppf(pvals(i),4.0_dp),4.0_dp), &
                pvals(i),5.0e-12_dp,5.0e-9_dp,failures)
        end do

        call check_true('gengamma singular lower density', &
            .not. ieee_is_finite(gengamma_pdf(0.0_dp,0.5_dp,1.0_dp)),failures)
        call check_close('gengamma finite lower density',gengamma_pdf(0.0_dp,1.0_dp,1.0_dp), &
            1.0_dp,0.0_dp,0.0_dp,failures)
        call check_close('gengamma zero lower density',gengamma_pdf(0.0_dp,2.0_dp,1.0_dp), &
            0.0_dp,0.0_dp,0.0_dp,failures)
        call check_close('gengamma negative-c lower density',gengamma_pdf(0.0_dp,2.3_dp,-1.3_dp), &
            0.0_dp,0.0_dp,0.0_dp,failures)
        call check_close('gengamma lower cdf',gengamma_cdf(0.1_dp,2.3_dp,-1.3_dp,0.1_dp,1.2_dp), &
            0.0_dp,0.0_dp,0.0_dp,failures)
        call check_close('halfgennorm lower cdf',halfgennorm_cdf(0.1_dp,1.7_dp,0.1_dp,1.2_dp), &
            0.0_dp,0.0_dp,0.0_dp,failures)
        call check_close('argus lower cdf',argus_cdf(0.1_dp,2.3_dp,0.1_dp,1.2_dp), &
            0.0_dp,0.0_dp,0.0_dp,failures)
        call check_close('argus upper cdf',argus_cdf(1.3_dp,2.3_dp,0.1_dp,1.2_dp), &
            1.0_dp,0.0_dp,0.0_dp,failures)
        call check_close('erlang lower cdf',erlang_cdf(0.1_dp,4.0_dp,0.1_dp,1.2_dp), &
            0.0_dp,0.0_dp,0.0_dp,failures)

        call check_true('gengamma invalid a',ieee_is_nan(gengamma_pdf(1.0_dp,-1.0_dp,1.3_dp)),failures)
        call check_true('gengamma invalid c',ieee_is_nan(gengamma_pdf(1.0_dp,2.3_dp,0.0_dp)),failures)
        call check_true('halfgennorm invalid beta',ieee_is_nan(halfgennorm_pdf(1.0_dp,0.0_dp)),failures)
        call check_true('argus invalid chi',ieee_is_nan(argus_pdf(0.5_dp,-1.0_dp)),failures)
        call check_true('erlang invalid a',ieee_is_nan(erlang_pdf(1.0_dp,0.0_dp)),failures)
        call check_close('argus lower endpoint density',argus_pdf(0.0_dp,2.3_dp), &
            0.0_dp,0.0_dp,0.0_dp,failures)
        call check_close('argus upper endpoint density',argus_pdf(1.0_dp,2.3_dp), &
            0.0_dp,0.0_dp,0.0_dp,failures)
        call check_true('gengamma upper quantile infinite',.not. ieee_is_finite(gengamma_ppf(1.0_dp,2.3_dp,1.3_dp)), &
            failures)

        a3 = argus_pdf([0.2_dp,0.5_dp,0.9_dp],2.3_dp)
        call check_close('argus elemental 1',a3(1),argus_pdf(0.2_dp,2.3_dp),0.0_dp,0.0_dp,failures)
        call check_close('argus elemental 3',a3(3),argus_pdf(0.9_dp,2.3_dp),0.0_dp,0.0_dp,failures)
    end subroutine test_identities_endpoints

    subroutine test_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: fd = 1.0e-6_dp, tol = 8.0e-4_dp
        real(dp), parameter :: gg_data(4) = [0.4_dp,0.8_dp,1.8_dp,3.0_dp]
        real(dp), parameter :: hg_data(4) = [0.4_dp,0.8_dp,1.8_dp,3.0_dp]
        real(dp), parameter :: ar_data(4) = [0.25_dp,0.55_dp,0.85_dp,1.1_dp]
        real(dp), parameter :: er_data(4) = [0.4_dp,0.8_dp,1.8_dp,3.0_dp]
        real(dp) :: s3(3), s4(4)

        s4 = gengamma_score(gg_data,2.3_dp,-1.3_dp,0.1_dp,1.2_dp)
        call check_close('gengamma score a',s4(1),centered( &
            gengamma_loglikelihood(gg_data,2.3_dp+fd,-1.3_dp,0.1_dp,1.2_dp), &
            gengamma_loglikelihood(gg_data,2.3_dp-fd,-1.3_dp,0.1_dp,1.2_dp),fd),tol,tol,failures)
        call check_close('gengamma score c',s4(2),centered( &
            gengamma_loglikelihood(gg_data,2.3_dp,-1.3_dp+fd,0.1_dp,1.2_dp), &
            gengamma_loglikelihood(gg_data,2.3_dp,-1.3_dp-fd,0.1_dp,1.2_dp),fd),tol,tol,failures)
        call check_close('gengamma score loc',s4(3),centered( &
            gengamma_loglikelihood(gg_data,2.3_dp,-1.3_dp,0.1_dp+fd,1.2_dp), &
            gengamma_loglikelihood(gg_data,2.3_dp,-1.3_dp,0.1_dp-fd,1.2_dp),fd),tol,tol,failures)
        call check_close('gengamma score scale',s4(4),centered( &
            gengamma_loglikelihood(gg_data,2.3_dp,-1.3_dp,0.1_dp,1.2_dp+fd), &
            gengamma_loglikelihood(gg_data,2.3_dp,-1.3_dp,0.1_dp,1.2_dp-fd),fd),tol,tol,failures)

        s3 = halfgennorm_score(hg_data,1.7_dp,0.1_dp,1.2_dp)
        call check_close('halfgennorm score beta',s3(1),centered( &
            halfgennorm_loglikelihood(hg_data,1.7_dp+fd,0.1_dp,1.2_dp), &
            halfgennorm_loglikelihood(hg_data,1.7_dp-fd,0.1_dp,1.2_dp),fd),tol,tol,failures)
        call check_close('halfgennorm score loc',s3(2),centered( &
            halfgennorm_loglikelihood(hg_data,1.7_dp,0.1_dp+fd,1.2_dp), &
            halfgennorm_loglikelihood(hg_data,1.7_dp,0.1_dp-fd,1.2_dp),fd),tol,tol,failures)
        call check_close('halfgennorm score scale',s3(3),centered( &
            halfgennorm_loglikelihood(hg_data,1.7_dp,0.1_dp,1.2_dp+fd), &
            halfgennorm_loglikelihood(hg_data,1.7_dp,0.1_dp,1.2_dp-fd),fd),tol,tol,failures)

        s3 = argus_score(ar_data,2.3_dp,0.1_dp,1.2_dp)
        call check_close('argus score chi',s3(1),centered( &
            argus_loglikelihood(ar_data,2.3_dp+fd,0.1_dp,1.2_dp), &
            argus_loglikelihood(ar_data,2.3_dp-fd,0.1_dp,1.2_dp),fd),tol,tol,failures)
        call check_close('argus score loc',s3(2),centered( &
            argus_loglikelihood(ar_data,2.3_dp,0.1_dp+fd,1.2_dp), &
            argus_loglikelihood(ar_data,2.3_dp,0.1_dp-fd,1.2_dp),fd),tol,tol,failures)
        call check_close('argus score scale',s3(3),centered( &
            argus_loglikelihood(ar_data,2.3_dp,0.1_dp,1.2_dp+fd), &
            argus_loglikelihood(ar_data,2.3_dp,0.1_dp,1.2_dp-fd),fd),tol,tol,failures)

        s3 = erlang_score(er_data,4.0_dp,0.1_dp,1.2_dp)
        call check_close('erlang score a',s3(1),centered( &
            erlang_loglikelihood(er_data,4.0_dp+fd,0.1_dp,1.2_dp), &
            erlang_loglikelihood(er_data,4.0_dp-fd,0.1_dp,1.2_dp),fd),tol,tol,failures)
        call check_close('erlang score loc',s3(2),centered( &
            erlang_loglikelihood(er_data,4.0_dp,0.1_dp+fd,1.2_dp), &
            erlang_loglikelihood(er_data,4.0_dp,0.1_dp-fd,1.2_dp),fd),tol,tol,failures)
        call check_close('erlang score scale',s3(3),centered( &
            erlang_loglikelihood(er_data,4.0_dp,0.1_dp,1.2_dp+fd), &
            erlang_loglikelihood(er_data,4.0_dp,0.1_dp,1.2_dp-fd),fd),tol,tol,failures)
    end subroutine test_scores

    subroutine test_rng_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: gg_data(4) = [0.4_dp,0.8_dp,1.8_dp,3.0_dp]
        real(dp), parameter :: hg_data(4) = [0.4_dp,0.8_dp,1.8_dp,3.0_dp]
        real(dp), parameter :: ar_data(4) = [0.25_dp,0.55_dp,0.85_dp,1.1_dp]
        real(dp), parameter :: er_data(4) = [0.4_dp,0.8_dp,1.8_dp,3.0_dp]
        real(dp) :: actual, expected
        type(rng_state) :: reference, state
        type(fit_result) :: result

        call rng_seed(state,24680)
        call rng_seed(reference,24680)
        actual = gengamma_rvs(state,2.3_dp,-1.3_dp,0.1_dp,1.2_dp)
        expected = gengamma_ppf(rng_uniform(reference),2.3_dp,-1.3_dp,0.1_dp,1.2_dp)
        call check_close('gengamma rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = halfgennorm_rvs(state,1.7_dp,0.1_dp,1.2_dp)
        expected = halfgennorm_ppf(rng_uniform(reference),1.7_dp,0.1_dp,1.2_dp)
        call check_close('halfgennorm rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = argus_rvs(state,2.3_dp,0.1_dp,1.2_dp)
        expected = argus_ppf(rng_uniform(reference),2.3_dp,0.1_dp,1.2_dp)
        call check_close('argus rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = erlang_rvs(state,4.0_dp,0.1_dp,1.2_dp)
        expected = erlang_ppf(rng_uniform(reference),4.0_dp,0.1_dp,1.2_dp)
        call check_close('erlang rvs',actual,expected,0.0_dp,0.0_dp,failures)

        call rng_seed(state,13579)
        call rng_seed(reference,13579)
        actual = argus_rvs(state,-1.0_dp)
        call check_true('invalid argus rvs NaN',ieee_is_nan(actual),failures)
        call check_close('invalid rvs preserves state',rng_uniform(state),rng_uniform(reference), &
            0.0_dp,0.0_dp,failures)

        call gengamma_fit(gg_data,[2.3_dp,-1.3_dp,0.1_dp,1.2_dp], &
            [2.3_dp,-1.3_dp,0.1_dp,1.2_dp],result)
        call check_true('gengamma fixed fit',result%success,failures)
        call halfgennorm_fit(hg_data,[1.7_dp,0.1_dp,1.2_dp],[1.7_dp,0.1_dp,1.2_dp],result)
        call check_true('halfgennorm fixed fit',result%success,failures)
        call argus_fit(ar_data,[2.3_dp,0.1_dp,1.2_dp],[2.3_dp,0.1_dp,1.2_dp],result)
        call check_true('argus fixed fit',result%success,failures)
        call erlang_fit(er_data,[4.0_dp,0.1_dp,1.2_dp],[4.0_dp,0.1_dp,1.2_dp],result)
        call check_true('erlang fixed fit',result%success,failures)

        call check_close('C ABI gengamma cdf', &
            scifort_gengamma_cdf_f64(0.8_dp,2.3_dp,-1.3_dp,0.0_dp,1.0_dp), &
            gengamma_cdf(0.8_dp,2.3_dp,-1.3_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI halfgennorm pdf', &
            scifort_halfgennorm_pdf_f64(0.8_dp,1.7_dp,0.0_dp,1.0_dp), &
            halfgennorm_pdf(0.8_dp,1.7_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI argus ppf',scifort_argus_ppf_f64(0.9_dp,2.3_dp,0.0_dp,1.0_dp), &
            argus_ppf(0.9_dp,2.3_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI erlang cdf',scifort_erlang_cdf_f64(3.0_dp,4.0_dp,0.0_dp,1.0_dp), &
            erlang_cdf(3.0_dp,4.0_dp),0.0_dp,0.0_dp,failures)
    end subroutine test_rng_fit_c_api

    pure function centered(plus, minus, h) result(value)
        real(dp), intent(in) :: plus !! objective at theta+h
        real(dp), intent(in) :: minus !! objective at theta-h
        real(dp), intent(in) :: h !! positive finite-difference step
        real(dp) :: value
        value = (plus - minus) / (2.0_dp * h)
    end function centered

end program test_gengamma_halfgennorm_argus_erlang
