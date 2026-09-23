! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_burr_mielke_gibrat_wrapcauchy
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_c_api
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state, rng_uniform
    use scifort_stats
    use test_burr_mielke_gibrat_wrapcauchy_reference
    use test_checks, only : check_close, check_true
    implicit none

    integer :: failures

    failures = 0
    call test_reference(failures)
    call test_identities_endpoints(failures)
    call test_scores(failures)
    call test_rng_fit_c_api(failures)

    if (failures == 0) then
        print '(a)', 'test_burr_mielke_gibrat_wrapcauchy: PASS'
    else
        print '(a,1x,i0)', 'test_burr_mielke_gibrat_wrapcauchy: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        integer :: i

        do i = 1, size(reference_probs)
            call check_close('burr pdf',burr_pdf(bu_x(i),2.3_dp,1.7_dp),bu_pdf(i), &
                8.0e-12_dp,8.0e-11_dp,failures)
            call check_close('burr logpdf',burr_logpdf(bu_x(i),2.3_dp,1.7_dp),bu_logpdf(i), &
                8.0e-12_dp,8.0e-11_dp,failures)
            call check_close('burr cdf',burr_cdf(bu_x(i),2.3_dp,1.7_dp),bu_cdf(i), &
                8.0e-12_dp,8.0e-11_dp,failures)
            call check_close('burr sf',burr_sf(bu_x(i),2.3_dp,1.7_dp),bu_sf(i), &
                8.0e-12_dp,8.0e-11_dp,failures)
            call check_close('burr logcdf',burr_logcdf(bu_x(i),2.3_dp,1.7_dp),bu_logcdf(i), &
                8.0e-12_dp,8.0e-11_dp,failures)
            call check_close('burr logsf',burr_logsf(bu_x(i),2.3_dp,1.7_dp),bu_logsf(i), &
                8.0e-12_dp,8.0e-11_dp,failures)
            call check_close('burr ppf',burr_ppf(reference_probs(i),2.3_dp,1.7_dp),bu_ppf(i), &
                2.0e-9_dp,2.0e-9_dp,failures)
            call check_close('burr isf',burr_isf(reference_probs(i),2.3_dp,1.7_dp),bu_isf(i), &
                2.0e-9_dp,2.0e-9_dp,failures)

            call check_close('mielke pdf',mielke_pdf(mi_x(i),2.4_dp,1.6_dp),mi_pdf(i), &
                8.0e-12_dp,8.0e-11_dp,failures)
            call check_close('mielke logpdf',mielke_logpdf(mi_x(i),2.4_dp,1.6_dp),mi_logpdf(i), &
                8.0e-12_dp,8.0e-11_dp,failures)
            call check_close('mielke cdf',mielke_cdf(mi_x(i),2.4_dp,1.6_dp),mi_cdf(i), &
                8.0e-12_dp,8.0e-11_dp,failures)
            call check_close('mielke sf',mielke_sf(mi_x(i),2.4_dp,1.6_dp),mi_sf(i), &
                8.0e-12_dp,8.0e-11_dp,failures)
            call check_close('mielke logcdf',mielke_logcdf(mi_x(i),2.4_dp,1.6_dp),mi_logcdf(i), &
                8.0e-12_dp,8.0e-11_dp,failures)
            call check_close('mielke logsf',mielke_logsf(mi_x(i),2.4_dp,1.6_dp),mi_logsf(i), &
                8.0e-12_dp,8.0e-11_dp,failures)
            call check_close('mielke ppf',mielke_ppf(reference_probs(i),2.4_dp,1.6_dp),mi_ppf(i), &
                3.0e-8_dp,3.0e-9_dp,failures)
            call check_close('mielke isf',mielke_isf(reference_probs(i),2.4_dp,1.6_dp),mi_isf(i), &
                3.0e-8_dp,3.0e-9_dp,failures)

            call check_close('gibrat pdf',gibrat_pdf(gi_x(i)),gi_pdf(i), &
                8.0e-12_dp,8.0e-11_dp,failures)
            call check_close('gibrat logpdf',gibrat_logpdf(gi_x(i)),gi_logpdf(i), &
                8.0e-12_dp,8.0e-11_dp,failures)
            call check_close('gibrat cdf',gibrat_cdf(gi_x(i)),gi_cdf(i), &
                8.0e-12_dp,8.0e-11_dp,failures)
            call check_close('gibrat sf',gibrat_sf(gi_x(i)),gi_sf(i), &
                8.0e-12_dp,8.0e-11_dp,failures)
            call check_close('gibrat logcdf',gibrat_logcdf(gi_x(i)),gi_logcdf(i), &
                8.0e-12_dp,8.0e-11_dp,failures)
            call check_close('gibrat logsf',gibrat_logsf(gi_x(i)),gi_logsf(i), &
                8.0e-12_dp,8.0e-11_dp,failures)
            call check_close('gibrat ppf',gibrat_ppf(reference_probs(i)),gi_ppf(i), &
                2.0e-10_dp,2.0e-9_dp,failures)
            call check_close('gibrat isf',gibrat_isf(reference_probs(i)),gi_isf(i), &
                2.0e-10_dp,2.0e-9_dp,failures)

            call check_close('wrapcauchy pdf',wrapcauchy_pdf(wc_x(i),0.4_dp),wc_pdf(i), &
                8.0e-12_dp,8.0e-11_dp,failures)
            call check_close('wrapcauchy logpdf',wrapcauchy_logpdf(wc_x(i),0.4_dp),wc_logpdf(i), &
                8.0e-12_dp,8.0e-11_dp,failures)
            call check_close('wrapcauchy cdf',wrapcauchy_cdf(wc_x(i),0.4_dp),wc_cdf(i), &
                8.0e-12_dp,8.0e-11_dp,failures)
            call check_close('wrapcauchy sf',wrapcauchy_sf(wc_x(i),0.4_dp),wc_sf(i), &
                8.0e-12_dp,8.0e-11_dp,failures)
            call check_close('wrapcauchy logcdf',wrapcauchy_logcdf(wc_x(i),0.4_dp),wc_logcdf(i), &
                8.0e-12_dp,8.0e-11_dp,failures)
            call check_close('wrapcauchy logsf',wrapcauchy_logsf(wc_x(i),0.4_dp),wc_logsf(i), &
                8.0e-12_dp,8.0e-11_dp,failures)
            call check_close('wrapcauchy ppf',wrapcauchy_ppf(reference_probs(i),0.4_dp),wc_ppf(i), &
                2.0e-10_dp,2.0e-9_dp,failures)
            call check_close('wrapcauchy isf',wrapcauchy_isf(reference_probs(i),0.4_dp),wc_isf(i), &
                2.0e-10_dp,2.0e-9_dp,failures)
        end do
    end subroutine test_reference

    subroutine test_identities_endpoints(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: xvals(5) = [0.08_dp,0.3_dp,0.9_dp,2.0_dp,12.0_dp]
        real(dp), parameter :: pvals(5) = [1.0e-12_dp,1.0e-6_dp,0.2_dp,0.8_dp,1.0_dp-1.0e-12_dp]
        real(dp) :: a(3), upper
        integer :: i

        do i = 1, size(xvals)
            call check_close('mielke-burr pdf',mielke_pdf(xvals(i),2.4_dp,1.6_dp), &
                burr_pdf(xvals(i),1.6_dp,1.5_dp),5.0e-12_dp,5.0e-11_dp,failures)
            call check_close('mielke-burr cdf',mielke_cdf(xvals(i),2.4_dp,1.6_dp), &
                burr_cdf(xvals(i),1.6_dp,1.5_dp),5.0e-12_dp,5.0e-11_dp,failures)
            call check_close('gibrat-lognormal pdf',gibrat_pdf(xvals(i)), &
                lognormal_pdf(xvals(i),1.0_dp),0.0_dp,0.0_dp,failures)
            call check_close('gibrat-lognormal cdf',gibrat_cdf(xvals(i)), &
                lognormal_cdf(xvals(i),1.0_dp),0.0_dp,0.0_dp,failures)
        end do
        do i = 1, size(pvals)
            call check_close('burr cdf-ppf',burr_cdf(burr_ppf(pvals(i),2.3_dp,1.7_dp), &
                2.3_dp,1.7_dp),pvals(i),2.0e-12_dp,2.0e-10_dp,failures)
            call check_close('mielke cdf-ppf',mielke_cdf(mielke_ppf(pvals(i),2.4_dp,1.6_dp), &
                2.4_dp,1.6_dp),pvals(i),2.0e-12_dp,2.0e-10_dp,failures)
            call check_close('gibrat cdf-ppf',gibrat_cdf(gibrat_ppf(pvals(i))),pvals(i), &
                2.0e-12_dp,2.0e-10_dp,failures)
            call check_close('wrapcauchy cdf-ppf',wrapcauchy_cdf(wrapcauchy_ppf(pvals(i),0.4_dp), &
                0.4_dp),pvals(i),2.0e-12_dp,2.0e-10_dp,failures)
            call check_close('burr sf-isf',burr_sf(burr_isf(pvals(i),2.3_dp,1.7_dp), &
                2.3_dp,1.7_dp),pvals(i),2.0e-12_dp,2.0e-10_dp,failures)
            call check_close('mielke sf-isf',mielke_sf(mielke_isf(pvals(i),2.4_dp,1.6_dp), &
                2.4_dp,1.6_dp),pvals(i),2.0e-12_dp,2.0e-10_dp,failures)
            call check_close('gibrat sf-isf',gibrat_sf(gibrat_isf(pvals(i))),pvals(i), &
                2.0e-12_dp,2.0e-10_dp,failures)
            call check_close('wrapcauchy sf-isf',wrapcauchy_sf(wrapcauchy_isf(pvals(i),0.4_dp), &
                0.4_dp),pvals(i),2.0e-12_dp,2.0e-10_dp,failures)
        end do

        call check_close('burr lower endpoint',burr_cdf(0.0_dp,2.3_dp,1.7_dp), &
            0.0_dp,0.0_dp,0.0_dp,failures)
        call check_close('burr ppf zero',burr_ppf(0.0_dp,2.3_dp,1.7_dp,0.3_dp,1.2_dp), &
            0.3_dp,0.0_dp,0.0_dp,failures)
        call check_true('burr ppf one inf',.not. ieee_is_finite(burr_ppf(1.0_dp,2.3_dp,1.7_dp)), &
            failures)
        call check_true('mielke invalid shape',ieee_is_nan(mielke_pdf(1.0_dp,0.0_dp,1.6_dp)),failures)
        call check_true('gibrat invalid scale',ieee_is_nan(gibrat_cdf(1.0_dp,0.0_dp,-1.0_dp)),failures)
        call check_true('wrapcauchy invalid c',ieee_is_nan(wrapcauchy_pdf(1.0_dp,1.0_dp)),failures)

        upper = 0.2_dp + 1.3_dp * (2.0_dp * scifort_pi)
        call check_close('wrapcauchy lower endpoint',wrapcauchy_cdf(0.2_dp,0.4_dp,0.2_dp,1.3_dp), &
            0.0_dp,0.0_dp,0.0_dp,failures)
        call check_close('wrapcauchy upper endpoint',wrapcauchy_cdf(upper,0.4_dp,0.2_dp,1.3_dp), &
            1.0_dp,0.0_dp,0.0_dp,failures)
        call check_close('wrapcauchy ppf one',wrapcauchy_ppf(1.0_dp,0.4_dp,0.2_dp,1.3_dp), &
            upper,4.0e-15_dp,4.0e-15_dp,failures)

        a = burr_pdf([0.3_dp,1.0_dp,3.0_dp],2.3_dp,1.7_dp)
        call check_close('burr elemental 1',a(1),burr_pdf(0.3_dp,2.3_dp,1.7_dp),0.0_dp,0.0_dp,failures)
        call check_close('burr elemental 3',a(3),burr_pdf(3.0_dp,2.3_dp,1.7_dp),0.0_dp,0.0_dp,failures)
    end subroutine test_identities_endpoints

    subroutine test_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp), parameter :: tol = 3.0e-4_dp
        real(dp), parameter :: bu_data(4) = [0.5_dp,0.9_dp,1.8_dp,4.0_dp]
        real(dp), parameter :: mi_data(4) = [0.5_dp,0.9_dp,1.8_dp,4.0_dp]
        real(dp), parameter :: gi_data(4) = [0.5_dp,0.9_dp,1.8_dp,4.0_dp]
        real(dp), parameter :: wc_data(4) = [0.6_dp,1.8_dp,4.0_dp,6.5_dp]
        real(dp) :: q2(2), q3(3), q4(4)

        q4 = burr_score(bu_data,2.3_dp,1.7_dp,0.1_dp,1.2_dp)
        call check_close('burr score c',q4(1),centered( &
            burr_loglikelihood(bu_data,2.3_dp+h,1.7_dp,0.1_dp,1.2_dp), &
            burr_loglikelihood(bu_data,2.3_dp-h,1.7_dp,0.1_dp,1.2_dp),h),tol,tol,failures)
        call check_close('burr score d',q4(2),centered( &
            burr_loglikelihood(bu_data,2.3_dp,1.7_dp+h,0.1_dp,1.2_dp), &
            burr_loglikelihood(bu_data,2.3_dp,1.7_dp-h,0.1_dp,1.2_dp),h),tol,tol,failures)
        call check_close('burr score loc',q4(3),centered( &
            burr_loglikelihood(bu_data,2.3_dp,1.7_dp,0.1_dp+h,1.2_dp), &
            burr_loglikelihood(bu_data,2.3_dp,1.7_dp,0.1_dp-h,1.2_dp),h),tol,tol,failures)
        call check_close('burr score scale',q4(4),centered( &
            burr_loglikelihood(bu_data,2.3_dp,1.7_dp,0.1_dp,1.2_dp+h), &
            burr_loglikelihood(bu_data,2.3_dp,1.7_dp,0.1_dp,1.2_dp-h),h),tol,tol,failures)

        q4 = mielke_score(mi_data,2.4_dp,1.6_dp,0.1_dp,1.2_dp)
        call check_close('mielke score k',q4(1),centered( &
            mielke_loglikelihood(mi_data,2.4_dp+h,1.6_dp,0.1_dp,1.2_dp), &
            mielke_loglikelihood(mi_data,2.4_dp-h,1.6_dp,0.1_dp,1.2_dp),h),tol,tol,failures)
        call check_close('mielke score s',q4(2),centered( &
            mielke_loglikelihood(mi_data,2.4_dp,1.6_dp+h,0.1_dp,1.2_dp), &
            mielke_loglikelihood(mi_data,2.4_dp,1.6_dp-h,0.1_dp,1.2_dp),h),tol,tol,failures)
        call check_close('mielke score loc',q4(3),centered( &
            mielke_loglikelihood(mi_data,2.4_dp,1.6_dp,0.1_dp+h,1.2_dp), &
            mielke_loglikelihood(mi_data,2.4_dp,1.6_dp,0.1_dp-h,1.2_dp),h),tol,tol,failures)
        call check_close('mielke score scale',q4(4),centered( &
            mielke_loglikelihood(mi_data,2.4_dp,1.6_dp,0.1_dp,1.2_dp+h), &
            mielke_loglikelihood(mi_data,2.4_dp,1.6_dp,0.1_dp,1.2_dp-h),h),tol,tol,failures)

        q2 = gibrat_score(gi_data,0.1_dp,1.2_dp)
        call check_close('gibrat score loc',q2(1),centered( &
            gibrat_loglikelihood(gi_data,0.1_dp+h,1.2_dp), &
            gibrat_loglikelihood(gi_data,0.1_dp-h,1.2_dp),h),tol,tol,failures)
        call check_close('gibrat score scale',q2(2),centered( &
            gibrat_loglikelihood(gi_data,0.1_dp,1.2_dp+h), &
            gibrat_loglikelihood(gi_data,0.1_dp,1.2_dp-h),h),tol,tol,failures)

        q3 = wrapcauchy_score(wc_data,0.4_dp,0.1_dp,1.2_dp)
        call check_close('wrapcauchy score c',q3(1),centered( &
            wrapcauchy_loglikelihood(wc_data,0.4_dp+h,0.1_dp,1.2_dp), &
            wrapcauchy_loglikelihood(wc_data,0.4_dp-h,0.1_dp,1.2_dp),h),tol,tol,failures)
        call check_close('wrapcauchy score loc',q3(2),centered( &
            wrapcauchy_loglikelihood(wc_data,0.4_dp,0.1_dp+h,1.2_dp), &
            wrapcauchy_loglikelihood(wc_data,0.4_dp,0.1_dp-h,1.2_dp),h),tol,tol,failures)
        call check_close('wrapcauchy score scale',q3(3),centered( &
            wrapcauchy_loglikelihood(wc_data,0.4_dp,0.1_dp,1.2_dp+h), &
            wrapcauchy_loglikelihood(wc_data,0.4_dp,0.1_dp,1.2_dp-h),h),tol,tol,failures)
    end subroutine test_scores

    subroutine test_rng_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: bu_data(4) = [0.5_dp,0.9_dp,1.8_dp,4.0_dp]
        real(dp), parameter :: mi_data(4) = [0.5_dp,0.9_dp,1.8_dp,4.0_dp]
        real(dp), parameter :: gi_data(4) = [0.5_dp,0.9_dp,1.8_dp,4.0_dp]
        real(dp), parameter :: wc_data(4) = [0.6_dp,1.8_dp,4.0_dp,6.5_dp]
        real(dp) :: actual, expected
        type(rng_state) :: reference, state
        type(fit_result) :: result

        call rng_seed(state,86420)
        call rng_seed(reference,86420)
        actual = burr_rvs(state,2.3_dp,1.7_dp,0.1_dp,1.2_dp)
        expected = burr_ppf(rng_uniform(reference),2.3_dp,1.7_dp,0.1_dp,1.2_dp)
        call check_close('burr rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = mielke_rvs(state,2.4_dp,1.6_dp,0.1_dp,1.2_dp)
        expected = mielke_ppf(rng_uniform(reference),2.4_dp,1.6_dp,0.1_dp,1.2_dp)
        call check_close('mielke rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = gibrat_rvs(state,0.1_dp,1.2_dp)
        expected = gibrat_ppf(rng_uniform(reference),0.1_dp,1.2_dp)
        call check_close('gibrat rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = wrapcauchy_rvs(state,0.4_dp,0.1_dp,1.2_dp)
        expected = wrapcauchy_ppf(rng_uniform(reference),0.4_dp,0.1_dp,1.2_dp)
        call check_close('wrapcauchy rvs',actual,expected,0.0_dp,0.0_dp,failures)

        call rng_seed(state,11111)
        call rng_seed(reference,11111)
        actual = burr_rvs(state,-1.0_dp,1.7_dp)
        call check_true('invalid burr rvs NaN',ieee_is_nan(actual),failures)
        call check_close('invalid rvs preserves state',rng_uniform(state),rng_uniform(reference), &
            0.0_dp,0.0_dp,failures)

        call burr_fit(bu_data,[2.3_dp,1.7_dp,0.1_dp,1.2_dp], &
            [2.3_dp,1.7_dp,0.1_dp,1.2_dp],result)
        call check_true('burr fixed fit',result%success,failures)
        call mielke_fit(mi_data,[2.4_dp,1.6_dp,0.1_dp,1.2_dp], &
            [2.4_dp,1.6_dp,0.1_dp,1.2_dp],result)
        call check_true('mielke fixed fit',result%success,failures)
        call gibrat_fit(gi_data,[0.1_dp,1.2_dp],[0.1_dp,1.2_dp],result)
        call check_true('gibrat fixed fit',result%success,failures)
        call wrapcauchy_fit(wc_data,[0.4_dp,0.1_dp,1.2_dp], &
            [0.4_dp,0.1_dp,1.2_dp],result)
        call check_true('wrapcauchy fixed fit',result%success,failures)

        call check_close('C ABI burr cdf', &
            scifort_burr_cdf_f64(1.0_dp,2.3_dp,1.7_dp,0.0_dp,1.0_dp), &
            burr_cdf(1.0_dp,2.3_dp,1.7_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI mielke pdf', &
            scifort_mielke_pdf_f64(1.0_dp,2.4_dp,1.6_dp,0.0_dp,1.0_dp), &
            mielke_pdf(1.0_dp,2.4_dp,1.6_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI gibrat ppf',scifort_gibrat_ppf_f64(0.9_dp,0.0_dp,1.0_dp), &
            gibrat_ppf(0.9_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI wrapcauchy cdf', &
            scifort_wrapcauchy_cdf_f64(0.8_dp,0.4_dp,0.0_dp,1.0_dp), &
            wrapcauchy_cdf(0.8_dp,0.4_dp),0.0_dp,0.0_dp,failures)
    end subroutine test_rng_fit_c_api

    pure function centered(plus, minus, h) result(value)
        real(dp), intent(in) :: plus !! objective at theta+h
        real(dp), intent(in) :: minus !! objective at theta-h
        real(dp), intent(in) :: h !! positive finite-difference step
        real(dp) :: value
        value = (plus - minus) / (2.0_dp * h)
    end function centered

end program test_burr_mielke_gibrat_wrapcauchy
