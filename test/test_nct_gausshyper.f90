! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_nct_gausshyper
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_c_api
    use scifort_fit, only : fit_result
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state, rng_uniform
    use scifort_stats
    use test_checks, only : check_close, check_true
    use test_nct_gausshyper_reference
    implicit none

    integer :: failures

    failures = 0
    call test_reference(failures)
    call test_identities_and_tails(failures)
    call test_scores(failures)
    call test_rng_fit_c_api(failures)

    if (failures == 0) then
        print '(a)', 'test_nct_gausshyper: PASS'
    else
        print '(a,1x,i0)', 'test_nct_gausshyper: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        integer :: i

        do i = 1, size(nct_x)
            call check_close('nct pdf',nct_pdf(nct_x(i),5.0_dp,2.0_dp), &
                nct_ref_pdf(i),2.0e-13_dp,3.0e-11_dp,failures)
            call check_close('nct logpdf',nct_logpdf(nct_x(i),5.0_dp,2.0_dp), &
                nct_ref_logpdf(i),2.0e-12_dp,3.0e-11_dp,failures)
            call check_close('nct cdf',nct_cdf(nct_x(i),5.0_dp,2.0_dp), &
                nct_ref_cdf(i),2.0e-13_dp,3.0e-10_dp,failures)
            call check_close('nct sf',nct_sf(nct_x(i),5.0_dp,2.0_dp), &
                nct_ref_sf(i),2.0e-13_dp,3.0e-10_dp,failures)
            call check_close('nct logcdf',nct_logcdf(nct_x(i),5.0_dp,2.0_dp), &
                nct_ref_logcdf(i),3.0e-11_dp,3.0e-10_dp,failures)
            call check_close('nct logsf',nct_logsf(nct_x(i),5.0_dp,2.0_dp), &
                nct_ref_logsf(i),3.0e-11_dp,3.0e-10_dp,failures)
        end do
        do i = 1, size(reference_probs)
            call check_close('nct ppf',nct_ppf(reference_probs(i),5.0_dp,2.0_dp), &
                nct_ref_ppf(i),3.0e-9_dp,3.0e-9_dp,failures)
            call check_close('nct isf',nct_isf(reference_probs(i),5.0_dp,2.0_dp), &
                nct_ref_isf(i),3.0e-9_dp,3.0e-9_dp,failures)
        end do

        do i = 1, size(gh_x)
            call check_close('gausshyper pdf',gausshyper_pdf(gh_x(i),2.0_dp,3.0_dp,1.5_dp,2.0_dp), &
                gh_pdf(i),3.0e-13_dp,3.0e-11_dp,failures)
            call check_close('gausshyper logpdf',gausshyper_logpdf(gh_x(i),2.0_dp,3.0_dp,1.5_dp,2.0_dp), &
                gh_logpdf(i),3.0e-12_dp,3.0e-11_dp,failures)
            call check_close('gausshyper cdf',gausshyper_cdf(gh_x(i),2.0_dp,3.0_dp,1.5_dp,2.0_dp), &
                gh_cdf(i),3.0e-12_dp,3.0e-10_dp,failures)
            call check_close('gausshyper sf',gausshyper_sf(gh_x(i),2.0_dp,3.0_dp,1.5_dp,2.0_dp), &
                gh_sf(i),3.0e-12_dp,3.0e-10_dp,failures)
            call check_close('gausshyper logcdf',gausshyper_logcdf(gh_x(i),2.0_dp,3.0_dp,1.5_dp,2.0_dp), &
                gh_logcdf(i),3.0e-11_dp,3.0e-10_dp,failures)
            call check_close('gausshyper logsf',gausshyper_logsf(gh_x(i),2.0_dp,3.0_dp,1.5_dp,2.0_dp), &
                gh_logsf(i),3.0e-11_dp,3.0e-10_dp,failures)
        end do
        do i = 1, size(reference_probs)
            call check_close('gausshyper ppf',gausshyper_ppf(reference_probs(i),2.0_dp,3.0_dp,1.5_dp,2.0_dp), &
                gh_ppf(i),3.0e-10_dp,3.0e-9_dp,failures)
            call check_close('gausshyper isf',gausshyper_isf(reference_probs(i),2.0_dp,3.0_dp,1.5_dp,2.0_dp), &
                gh_isf(i),3.0e-10_dp,3.0e-9_dp,failures)
        end do
    end subroutine test_reference

    subroutine test_identities_and_tails(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        integer :: i
        real(dp), parameter :: xs(7) = [-4.0_dp,-2.0_dp,-0.5_dp,0.0_dp,0.7_dp,2.0_dp,5.0_dp]
        real(dp) :: x

        do i = 1, size(xs)
            x = xs(i)
            call check_close('nct central pdf',nct_pdf(x,7.0_dp,0.0_dp), &
                t_pdf(x,7.0_dp),2.0e-13_dp,3.0e-11_dp,failures)
            call check_close('nct central cdf',nct_cdf(x,7.0_dp,0.0_dp), &
                t_cdf(x,7.0_dp),2.0e-13_dp,3.0e-10_dp,failures)
        end do
        call check_close('nct tails',nct_cdf(1.7_dp,5.0_dp,2.0_dp)+nct_sf(1.7_dp,5.0_dp,2.0_dp), &
            1.0_dp,3.0e-14_dp,3.0e-14_dp,failures)
        call check_close('nct shifted scaled',nct_cdf(2.4_dp,5.0_dp,1.2_dp,0.4_dp,2.0_dp), &
            nct_cdf(1.0_dp,5.0_dp,1.2_dp),3.0e-14_dp,3.0e-14_dp,failures)

        call check_close('gausshyper beta c=0 pdf',gausshyper_pdf(0.37_dp,2.0_dp,3.0_dp,0.0_dp,8.0_dp), &
            beta_pdf(0.37_dp,2.0_dp,3.0_dp),3.0e-13_dp,3.0e-11_dp,failures)
        call check_close('gausshyper beta c=0 cdf',gausshyper_cdf(0.37_dp,2.0_dp,3.0_dp,0.0_dp,8.0_dp), &
            beta_cdf(0.37_dp,2.0_dp,3.0_dp),3.0e-12_dp,3.0e-10_dp,failures)
        call check_close('gausshyper beta z=0 pdf',gausshyper_pdf(0.37_dp,2.0_dp,3.0_dp,4.0_dp,0.0_dp), &
            beta_pdf(0.37_dp,2.0_dp,3.0_dp),3.0e-13_dp,3.0e-11_dp,failures)
        call check_close('gausshyper beta z=0 cdf',gausshyper_cdf(0.37_dp,2.0_dp,3.0_dp,4.0_dp,0.0_dp), &
            beta_cdf(0.37_dp,2.0_dp,3.0_dp),3.0e-12_dp,3.0e-10_dp,failures)
        call check_close('gausshyper tails', &
            gausshyper_cdf(0.63_dp,2.0_dp,3.0_dp,1.5_dp,2.0_dp)+ &
            gausshyper_sf(0.63_dp,2.0_dp,3.0_dp,1.5_dp,2.0_dp), &
            1.0_dp,3.0e-13_dp,3.0e-13_dp,failures)

        ! Independent defining-integral checks.  SciPy 1.17/Boost loses the
        ! noncentral-t tail below and SciPy's generic gausshyper integration
        ! loses the two small upper tails below.
        call check_close('nct direct deep logsf',nct_logsf(10.0_dp,30.0_dp,-5.0_dp), &
            -57.214203664073054_dp,2.0e-9_dp,3.0e-11_dp,failures)
        call check_close('nct small-df cdf',nct_cdf(10.0_dp,0.2_dp,6.0_dp), &
            0.2504117582414005_dp,3.0e-13_dp,3.0e-12_dp,failures)
        call check_close('nct small-df pdf',nct_pdf(10.0_dp,0.2_dp,6.0_dp), &
            0.01449116494956138_dp,3.0e-14_dp,3.0e-12_dp,failures)
        call check_close('gausshyper direct deep logsf', &
            gausshyper_logsf(0.8046156546172839_dp,6.958263543888631_dp, &
                39.14506398459203_dp,10.510642249857504_dp,6.881072658767391_dp), &
            -63.29983233932789_dp,3.0e-9_dp,3.0e-11_dp,failures)
        call check_close('gausshyper concentrated logsf', &
            gausshyper_logsf(0.5_dp,50.0_dp,80.0_dp,100.0_dp,5.0_dp), &
            -55.36602535445925_dp,3.0e-8_dp,3.0e-10_dp,failures)
        call check_close('gausshyper hard logpdf', &
            gausshyper_logpdf(0.884632573494884_dp,1.1758527186186538_dp, &
                45.84409749288438_dp,3.803538990491404_dp,2.024498237278596_dp), &
            -96.02140974354578_dp,3.0e-9_dp,3.0e-11_dp,failures)

        call check_true('nct invalid df',ieee_is_nan(nct_pdf(0.0_dp,0.0_dp,1.0_dp)),failures)
        call check_true('gausshyper invalid a', &
            ieee_is_nan(gausshyper_pdf(0.5_dp,0.0_dp,2.0_dp,1.0_dp,1.0_dp)),failures)
        call check_true('gausshyper invalid z', &
            ieee_is_nan(gausshyper_pdf(0.5_dp,2.0_dp,2.0_dp,1.0_dp,-1.0_dp)),failures)
    end subroutine test_identities_and_tails

    subroutine test_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: h = 2.0e-6_dp
        real(dp), parameter :: nd(4) = [-1.0_dp,0.2_dp,1.3_dp,3.0_dp]
        real(dp), parameter :: gd(4) = [0.15_dp,0.30_dp,0.55_dp,0.90_dp]
        real(dp) :: ns(4), gs(6), fd

        ns = nct_score(nd,5.0_dp,1.5_dp,0.1_dp,1.2_dp)
        fd = centered(nct_loglikelihood(nd,5.0_dp+h,1.5_dp,0.1_dp,1.2_dp), &
            nct_loglikelihood(nd,5.0_dp-h,1.5_dp,0.1_dp,1.2_dp),h)
        call check_close('nct score df',ns(1),fd,2.0e-5_dp,2.0e-5_dp,failures)
        fd = centered(nct_loglikelihood(nd,5.0_dp,1.5_dp+h,0.1_dp,1.2_dp), &
            nct_loglikelihood(nd,5.0_dp,1.5_dp-h,0.1_dp,1.2_dp),h)
        call check_close('nct score nc',ns(2),fd,2.0e-5_dp,2.0e-5_dp,failures)
        fd = centered(nct_loglikelihood(nd,5.0_dp,1.5_dp,0.1_dp+h,1.2_dp), &
            nct_loglikelihood(nd,5.0_dp,1.5_dp,0.1_dp-h,1.2_dp),h)
        call check_close('nct score loc',ns(3),fd,2.0e-5_dp,2.0e-5_dp,failures)
        fd = centered(nct_loglikelihood(nd,5.0_dp,1.5_dp,0.1_dp,1.2_dp+h), &
            nct_loglikelihood(nd,5.0_dp,1.5_dp,0.1_dp,1.2_dp-h),h)
        call check_close('nct score scale',ns(4),fd,2.0e-5_dp,2.0e-5_dp,failures)

        gs = gausshyper_score(gd,2.0_dp,3.0_dp,1.5_dp,2.0_dp,0.0_dp,1.2_dp)
        fd = centered(gausshyper_loglikelihood(gd,2.0_dp+h,3.0_dp,1.5_dp,2.0_dp,0.0_dp,1.2_dp), &
            gausshyper_loglikelihood(gd,2.0_dp-h,3.0_dp,1.5_dp,2.0_dp,0.0_dp,1.2_dp),h)
        call check_close('gausshyper score a',gs(1),fd,4.0e-5_dp,4.0e-5_dp,failures)
        fd = centered(gausshyper_loglikelihood(gd,2.0_dp,3.0_dp+h,1.5_dp,2.0_dp,0.0_dp,1.2_dp), &
            gausshyper_loglikelihood(gd,2.0_dp,3.0_dp-h,1.5_dp,2.0_dp,0.0_dp,1.2_dp),h)
        call check_close('gausshyper score b',gs(2),fd,4.0e-5_dp,4.0e-5_dp,failures)
        fd = centered(gausshyper_loglikelihood(gd,2.0_dp,3.0_dp,1.5_dp+h,2.0_dp,0.0_dp,1.2_dp), &
            gausshyper_loglikelihood(gd,2.0_dp,3.0_dp,1.5_dp-h,2.0_dp,0.0_dp,1.2_dp),h)
        call check_close('gausshyper score c',gs(3),fd,4.0e-5_dp,4.0e-5_dp,failures)
        fd = centered(gausshyper_loglikelihood(gd,2.0_dp,3.0_dp,1.5_dp,2.0_dp+h,0.0_dp,1.2_dp), &
            gausshyper_loglikelihood(gd,2.0_dp,3.0_dp,1.5_dp,2.0_dp-h,0.0_dp,1.2_dp),h)
        call check_close('gausshyper score z',gs(4),fd,4.0e-5_dp,4.0e-5_dp,failures)
        fd = centered(gausshyper_loglikelihood(gd,2.0_dp,3.0_dp,1.5_dp,2.0_dp,h,1.2_dp), &
            gausshyper_loglikelihood(gd,2.0_dp,3.0_dp,1.5_dp,2.0_dp,-h,1.2_dp),h)
        call check_close('gausshyper score loc',gs(5),fd,4.0e-5_dp,4.0e-5_dp,failures)
        fd = centered(gausshyper_loglikelihood(gd,2.0_dp,3.0_dp,1.5_dp,2.0_dp,0.0_dp,1.2_dp+h), &
            gausshyper_loglikelihood(gd,2.0_dp,3.0_dp,1.5_dp,2.0_dp,0.0_dp,1.2_dp-h),h)
        call check_close('gausshyper score scale',gs(6),fd,4.0e-5_dp,4.0e-5_dp,failures)
    end subroutine test_scores

    subroutine test_rng_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: nd(4) = [-1.0_dp,0.2_dp,1.3_dp,3.0_dp]
        real(dp), parameter :: gd(4) = [0.15_dp,0.30_dp,0.55_dp,0.90_dp]
        type(rng_state) :: state, reference
        type(fit_result) :: result
        real(dp) :: actual, expected

        call rng_seed(state,24680); call rng_seed(reference,24680)
        actual = nct_rvs(state,5.0_dp,1.5_dp)
        expected = nct_ppf(rng_uniform(reference),5.0_dp,1.5_dp)
        call check_close('nct rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = gausshyper_rvs(state,2.0_dp,3.0_dp,1.5_dp,2.0_dp)
        expected = gausshyper_ppf(rng_uniform(reference),2.0_dp,3.0_dp,1.5_dp,2.0_dp)
        call check_close('gausshyper rvs',actual,expected,0.0_dp,0.0_dp,failures)

        call rng_seed(state,13579); call rng_seed(reference,13579)
        actual = gausshyper_rvs(state,0.0_dp,3.0_dp,1.5_dp,2.0_dp)
        call check_true('invalid gausshyper rvs NaN',ieee_is_nan(actual),failures)
        call check_close('invalid rvs preserves state',rng_uniform(state),rng_uniform(reference), &
            0.0_dp,0.0_dp,failures)

        call nct_fit(nd,[5.0_dp,1.5_dp,0.1_dp,1.2_dp], &
            [5.0_dp,1.5_dp,0.1_dp,1.2_dp],result)
        call check_true('nct fixed fit',result%success,failures)
        call gausshyper_fit(gd,[2.0_dp,3.0_dp,1.5_dp,2.0_dp,0.0_dp,1.2_dp], &
            [2.0_dp,3.0_dp,1.5_dp,2.0_dp,0.0_dp,1.2_dp],result)
        call check_true('gausshyper fixed fit',result%success,failures)

        call check_close('C ABI nct pdf',scifort_nct_pdf_f64(0.3_dp,5.0_dp,1.5_dp,0.1_dp,1.2_dp), &
            nct_pdf(0.3_dp,5.0_dp,1.5_dp,0.1_dp,1.2_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI nct cdf',scifort_nct_cdf_f64(0.3_dp,5.0_dp,1.5_dp,0.1_dp,1.2_dp), &
            nct_cdf(0.3_dp,5.0_dp,1.5_dp,0.1_dp,1.2_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI nct ppf',scifort_nct_ppf_f64(0.7_dp,5.0_dp,1.5_dp,0.1_dp,1.2_dp), &
            nct_ppf(0.7_dp,5.0_dp,1.5_dp,0.1_dp,1.2_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI gausshyper pdf', &
            scifort_gausshyper_pdf_f64(0.3_dp,2.0_dp,3.0_dp,1.5_dp,2.0_dp,0.0_dp,1.0_dp), &
            gausshyper_pdf(0.3_dp,2.0_dp,3.0_dp,1.5_dp,2.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI gausshyper cdf', &
            scifort_gausshyper_cdf_f64(0.3_dp,2.0_dp,3.0_dp,1.5_dp,2.0_dp,0.0_dp,1.0_dp), &
            gausshyper_cdf(0.3_dp,2.0_dp,3.0_dp,1.5_dp,2.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI gausshyper ppf', &
            scifort_gausshyper_ppf_f64(0.7_dp,2.0_dp,3.0_dp,1.5_dp,2.0_dp,0.0_dp,1.0_dp), &
            gausshyper_ppf(0.7_dp,2.0_dp,3.0_dp,1.5_dp,2.0_dp),0.0_dp,0.0_dp,failures)
    end subroutine test_rng_fit_c_api

    pure function centered(plus_value, minus_value, h) result(value)
        real(dp), intent(in) :: plus_value !! function value at theta+h
        real(dp), intent(in) :: minus_value !! function value at theta-h
        real(dp), intent(in) :: h !! positive centered-difference step
        real(dp) :: value
        value = (plus_value-minus_value)/(2.0_dp*h)
    end function centered

end program test_nct_gausshyper
