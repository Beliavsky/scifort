! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_irwinhall_ksone_ncx2_ncf
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_c_api
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state, rng_uniform
    use scifort_stats
    use test_checks, only : check_close, check_true
    use test_irwinhall_ksone_ncx2_ncf_reference
    implicit none

    integer :: failures

    failures = 0
    call test_reference(failures)
    call test_identities_endpoints(failures)
    call test_scores(failures)
    call test_rng_fit_c_api(failures)

    if (failures == 0) then
        print '(a)', 'test_irwinhall_ksone_ncx2_ncf: PASS'
    else
        print '(a,1x,i0)', 'test_irwinhall_ksone_ncx2_ncf: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        integer :: i

        do i = 1, size(reference_probs)
            call check_close('irwinhall pdf',irwinhall_pdf(ih_x(i),5.0_dp),ih_pdf(i), &
                5.0e-13_dp,5.0e-12_dp,failures)
            call check_close('irwinhall logpdf',irwinhall_logpdf(ih_x(i),5.0_dp),ih_logpdf(i), &
                5.0e-13_dp,5.0e-12_dp,failures)
            call check_close('irwinhall cdf',irwinhall_cdf(ih_x(i),5.0_dp),ih_cdf(i), &
                5.0e-13_dp,5.0e-12_dp,failures)
            call check_close('irwinhall sf',irwinhall_sf(ih_x(i),5.0_dp),ih_sf(i), &
                5.0e-13_dp,5.0e-12_dp,failures)
            call check_close('irwinhall logcdf',irwinhall_logcdf(ih_x(i),5.0_dp),ih_logcdf(i), &
                2.0e-12_dp,2.0e-11_dp,failures)
            call check_close('irwinhall logsf',irwinhall_logsf(ih_x(i),5.0_dp),ih_logsf(i), &
                2.0e-12_dp,2.0e-11_dp,failures)
            call check_close('irwinhall ppf',irwinhall_ppf(reference_probs(i),5.0_dp),ih_ppf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('irwinhall isf',irwinhall_isf(reference_probs(i),5.0_dp),ih_isf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)

            call check_close('ksone pdf',ksone_pdf(ko_x(i),20.0_dp),ko_pdf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('ksone logpdf',ksone_logpdf(ko_x(i),20.0_dp),ko_logpdf(i), &
                2.0e-11_dp,2.0e-10_dp,failures)
            call check_close('ksone cdf',ksone_cdf(ko_x(i),20.0_dp),ko_cdf(i), &
                2.0e-12_dp,2.0e-11_dp,failures)
            call check_close('ksone sf',ksone_sf(ko_x(i),20.0_dp),ko_sf(i), &
                2.0e-12_dp,2.0e-11_dp,failures)
            call check_close('ksone logcdf',ksone_logcdf(ko_x(i),20.0_dp),ko_logcdf(i), &
                3.0e-11_dp,3.0e-10_dp,failures)
            call check_close('ksone logsf',ksone_logsf(ko_x(i),20.0_dp),ko_logsf(i), &
                3.0e-11_dp,3.0e-10_dp,failures)
            call check_close('ksone ppf',ksone_ppf(reference_probs(i),20.0_dp),ko_ppf(i), &
                2.0e-10_dp,2.0e-9_dp,failures)
            call check_close('ksone isf',ksone_isf(reference_probs(i),20.0_dp),ko_isf(i), &
                2.0e-10_dp,2.0e-9_dp,failures)

            call check_close('ncx2 pdf',ncx2_pdf(nx_x(i),5.0_dp,3.0_dp),nx_pdf(i), &
                3.0e-12_dp,3.0e-11_dp,failures)
            call check_close('ncx2 logpdf',ncx2_logpdf(nx_x(i),5.0_dp,3.0_dp),nx_logpdf(i), &
                3.0e-12_dp,3.0e-11_dp,failures)
            call check_close('ncx2 cdf',ncx2_cdf(nx_x(i),5.0_dp,3.0_dp),nx_cdf(i), &
                3.0e-12_dp,3.0e-11_dp,failures)
            call check_close('ncx2 sf',ncx2_sf(nx_x(i),5.0_dp,3.0_dp),nx_sf(i), &
                3.0e-12_dp,3.0e-11_dp,failures)
            call check_close('ncx2 logcdf',ncx2_logcdf(nx_x(i),5.0_dp,3.0_dp),nx_logcdf(i), &
                3.0e-11_dp,3.0e-10_dp,failures)
            call check_close('ncx2 logsf',ncx2_logsf(nx_x(i),5.0_dp,3.0_dp),nx_logsf(i), &
                3.0e-11_dp,3.0e-10_dp,failures)
            call check_close('ncx2 ppf',ncx2_ppf(reference_probs(i),5.0_dp,3.0_dp),nx_ppf(i), &
                3.0e-10_dp,3.0e-9_dp,failures)
            call check_close('ncx2 isf',ncx2_isf(reference_probs(i),5.0_dp,3.0_dp),nx_isf(i), &
                3.0e-10_dp,3.0e-9_dp,failures)

            call check_close('ncf pdf',ncf_pdf(nf_x(i),5.0_dp,12.0_dp,3.0_dp),nf_pdf(i), &
                5.0e-12_dp,5.0e-11_dp,failures)
            call check_close('ncf logpdf',ncf_logpdf(nf_x(i),5.0_dp,12.0_dp,3.0_dp),nf_logpdf(i), &
                5.0e-12_dp,5.0e-11_dp,failures)
            call check_close('ncf cdf',ncf_cdf(nf_x(i),5.0_dp,12.0_dp,3.0_dp),nf_cdf(i), &
                5.0e-12_dp,5.0e-11_dp,failures)
            call check_close('ncf sf',ncf_sf(nf_x(i),5.0_dp,12.0_dp,3.0_dp),nf_sf(i), &
                5.0e-12_dp,5.0e-11_dp,failures)
            call check_close('ncf logcdf',ncf_logcdf(nf_x(i),5.0_dp,12.0_dp,3.0_dp),nf_logcdf(i), &
                5.0e-11_dp,5.0e-10_dp,failures)
            call check_close('ncf logsf',ncf_logsf(nf_x(i),5.0_dp,12.0_dp,3.0_dp),nf_logsf(i), &
                5.0e-11_dp,5.0e-10_dp,failures)
            call check_close('ncf ppf',ncf_ppf(reference_probs(i),5.0_dp,12.0_dp,3.0_dp),nf_ppf(i), &
                5.0e-9_dp,5.0e-9_dp,failures)
            call check_close('ncf isf',ncf_isf(reference_probs(i),5.0_dp,12.0_dp,3.0_dp),nf_isf(i), &
                5.0e-9_dp,5.0e-9_dp,failures)
        end do
        call check_close('ksone tiny ppf',ksone_ppf(1.0e-10_dp,250.0_dp),ko_tiny_ppf(1), &
            5.0e-24_dp,5.0e-13_dp,failures)
        call check_close('ncx2 tiny ppf',ncx2_ppf(1.0e-9_dp,0.5_dp,0.0_dp),nx_tiny_ppf(1), &
            5.0e-49_dp,5.0e-13_dp,failures)
        call check_close('ncx2 far-tail pdf',ncx2_pdf(900.0_dp,80.0_dp,100.0_dp),nx_far_pdf(1), &
            1.0e-84_dp,5.0e-12_dp,failures)
        call check_close('ncf far-tail sf',ncf_sf(378227639.49492204_dp,5.0_dp,2.5_dp,100.0_dp), &
            nf_far_sf(1),2.0e-21_dp,5.0e-12_dp,failures)
    end subroutine test_reference

    subroutine test_identities_endpoints(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: pvals(5) = [1.0e-9_dp,1.0e-5_dp,0.2_dp,0.8_dp,1.0_dp-1.0e-9_dp]
        real(dp), parameter :: xvals(4) = [0.15_dp,0.4_dp,0.7_dp,0.93_dp]
        integer :: i

        do i = 1, size(xvals)
            call check_close('irwinhall n1 uniform pdf',irwinhall_pdf(xvals(i),1.0_dp), &
                uniform_pdf(xvals(i)),2.0e-15_dp,2.0e-15_dp,failures)
            call check_close('irwinhall n1 uniform cdf',irwinhall_cdf(xvals(i),1.0_dp), &
                uniform_cdf(xvals(i)),2.0e-15_dp,2.0e-15_dp,failures)
            call check_close('ksone n1 uniform pdf',ksone_pdf(xvals(i),1.0_dp), &
                uniform_pdf(xvals(i)),2.0e-13_dp,2.0e-12_dp,failures)
            call check_close('ksone n1 uniform cdf',ksone_cdf(xvals(i),1.0_dp), &
                uniform_cdf(xvals(i)),2.0e-13_dp,2.0e-12_dp,failures)
            call check_close('ncx2 central pdf',ncx2_pdf(3.0_dp*xvals(i),5.0_dp,0.0_dp), &
                chi2_pdf(3.0_dp*xvals(i),5.0_dp),3.0e-13_dp,3.0e-12_dp,failures)
            call check_close('ncx2 central cdf',ncx2_cdf(3.0_dp*xvals(i),5.0_dp,0.0_dp), &
                chi2_cdf(3.0_dp*xvals(i),5.0_dp),3.0e-13_dp,3.0e-12_dp,failures)
            call check_close('ncf central pdf',ncf_pdf(2.0_dp*xvals(i),5.0_dp,12.0_dp,0.0_dp), &
                f_pdf(2.0_dp*xvals(i),5.0_dp,12.0_dp),3.0e-13_dp,3.0e-12_dp,failures)
            call check_close('ncf central cdf',ncf_cdf(2.0_dp*xvals(i),5.0_dp,12.0_dp,0.0_dp), &
                f_cdf(2.0_dp*xvals(i),5.0_dp,12.0_dp),3.0e-13_dp,3.0e-12_dp,failures)
        end do
        do i = 1, size(pvals)
            call check_close('irwinhall cdf-ppf',irwinhall_cdf(irwinhall_ppf(pvals(i),5.0_dp),5.0_dp), &
                pvals(i),3.0e-11_dp,3.0e-9_dp,failures)
            call check_close('irwinhall sf-isf',irwinhall_sf(irwinhall_isf(pvals(i),5.0_dp),5.0_dp), &
                pvals(i),3.0e-11_dp,3.0e-9_dp,failures)
            call check_close('ksone cdf-ppf',ksone_cdf(ksone_ppf(pvals(i),20.0_dp),20.0_dp), &
                pvals(i),3.0e-10_dp,3.0e-8_dp,failures)
            call check_close('ksone sf-isf',ksone_sf(ksone_isf(pvals(i),20.0_dp),20.0_dp), &
                pvals(i),3.0e-10_dp,3.0e-8_dp,failures)
            call check_close('ncx2 cdf-ppf',ncx2_cdf(ncx2_ppf(pvals(i),5.0_dp,3.0_dp),5.0_dp,3.0_dp), &
                pvals(i),3.0e-10_dp,3.0e-8_dp,failures)
            call check_close('ncx2 sf-isf',ncx2_sf(ncx2_isf(pvals(i),5.0_dp,3.0_dp),5.0_dp,3.0_dp), &
                pvals(i),3.0e-10_dp,3.0e-8_dp,failures)
            call check_close('ncf cdf-ppf',ncf_cdf(ncf_ppf(pvals(i),5.0_dp,12.0_dp,3.0_dp), &
                5.0_dp,12.0_dp,3.0_dp),pvals(i),5.0e-10_dp,5.0e-8_dp,failures)
            call check_close('ncf sf-isf',ncf_sf(ncf_isf(pvals(i),5.0_dp,12.0_dp,3.0_dp), &
                5.0_dp,12.0_dp,3.0_dp),pvals(i),5.0e-10_dp,5.0e-8_dp,failures)
        end do
        call check_close('irwinhall lower cdf',irwinhall_cdf(0.0_dp,5.0_dp),0.0_dp,0.0_dp,0.0_dp,failures)
        call check_close('irwinhall upper cdf',irwinhall_cdf(5.0_dp,5.0_dp),1.0_dp,0.0_dp,0.0_dp,failures)
        call check_close('ksone lower cdf',ksone_cdf(0.0_dp,20.0_dp),0.0_dp,0.0_dp,0.0_dp,failures)
        call check_close('ksone upper cdf',ksone_cdf(1.0_dp,20.0_dp),1.0_dp,0.0_dp,0.0_dp,failures)
        call check_true('irwinhall invalid noninteger n',ieee_is_nan(irwinhall_pdf(1.0_dp,2.5_dp)),failures)
        call check_true('ksone invalid noninteger n',ieee_is_nan(ksone_pdf(0.2_dp,3.5_dp)),failures)
        call check_true('ncx2 invalid nc',ieee_is_nan(ncx2_pdf(1.0_dp,5.0_dp,-1.0_dp)),failures)
        call check_true('ncf invalid df',ieee_is_nan(ncf_pdf(1.0_dp,0.0_dp,12.0_dp,3.0_dp)),failures)
        call check_true('ncx2 ppf one infinite',.not. ieee_is_finite(ncx2_ppf(1.0_dp,5.0_dp,3.0_dp)),failures)
        call check_true('ncf ppf one infinite',.not. ieee_is_finite(ncf_ppf(1.0_dp,5.0_dp,12.0_dp,3.0_dp)),failures)
    end subroutine test_identities_endpoints

    subroutine test_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: h = 1.0e-6_dp, tol = 7.0e-4_dp
        real(dp), parameter :: ihd(4) = [0.6_dp,1.5_dp,2.6_dp,4.2_dp]
        real(dp), parameter :: kod(4) = [0.06_dp,0.13_dp,0.21_dp,0.31_dp]
        real(dp), parameter :: nxd(4) = [1.1_dp,3.0_dp,7.5_dp,12.0_dp]
        real(dp), parameter :: nfd(4) = [0.3_dp,0.8_dp,1.7_dp,3.2_dp]
        real(dp) :: s2(2), s4(4), s5(5)

        s2 = irwinhall_score(ihd,5.0_dp,0.1_dp,1.2_dp)
        call check_close('irwinhall score loc',s2(1),centered( &
            irwinhall_loglikelihood(ihd,5.0_dp,0.1_dp+h,1.2_dp), &
            irwinhall_loglikelihood(ihd,5.0_dp,0.1_dp-h,1.2_dp),h),tol,tol,failures)
        call check_close('irwinhall score scale',s2(2),centered( &
            irwinhall_loglikelihood(ihd,5.0_dp,0.1_dp,1.2_dp+h), &
            irwinhall_loglikelihood(ihd,5.0_dp,0.1_dp,1.2_dp-h),h),tol,tol,failures)

        s2 = ksone_score(kod,20.0_dp,0.01_dp,1.1_dp)
        call check_close('ksone score loc',s2(1),centered( &
            ksone_loglikelihood(kod,20.0_dp,0.01_dp+h,1.1_dp), &
            ksone_loglikelihood(kod,20.0_dp,0.01_dp-h,1.1_dp),h),tol,tol,failures)
        call check_close('ksone score scale',s2(2),centered( &
            ksone_loglikelihood(kod,20.0_dp,0.01_dp,1.1_dp+h), &
            ksone_loglikelihood(kod,20.0_dp,0.01_dp,1.1_dp-h),h),tol,tol,failures)

        s4 = ncx2_score(nxd,5.5_dp,2.5_dp,0.2_dp,1.1_dp)
        call check_close('ncx2 score df',s4(1),centered( &
            ncx2_loglikelihood(nxd,5.5_dp+h,2.5_dp,0.2_dp,1.1_dp), &
            ncx2_loglikelihood(nxd,5.5_dp-h,2.5_dp,0.2_dp,1.1_dp),h),tol,tol,failures)
        call check_close('ncx2 score nc',s4(2),centered( &
            ncx2_loglikelihood(nxd,5.5_dp,2.5_dp+h,0.2_dp,1.1_dp), &
            ncx2_loglikelihood(nxd,5.5_dp,2.5_dp-h,0.2_dp,1.1_dp),h),tol,tol,failures)
        call check_close('ncx2 score loc',s4(3),centered( &
            ncx2_loglikelihood(nxd,5.5_dp,2.5_dp,0.2_dp+h,1.1_dp), &
            ncx2_loglikelihood(nxd,5.5_dp,2.5_dp,0.2_dp-h,1.1_dp),h),tol,tol,failures)
        call check_close('ncx2 score scale',s4(4),centered( &
            ncx2_loglikelihood(nxd,5.5_dp,2.5_dp,0.2_dp,1.1_dp+h), &
            ncx2_loglikelihood(nxd,5.5_dp,2.5_dp,0.2_dp,1.1_dp-h),h),tol,tol,failures)

        s5 = ncf_score(nfd,5.5_dp,11.5_dp,2.7_dp,0.1_dp,1.2_dp)
        call check_close('ncf score dfn',s5(1),centered( &
            ncf_loglikelihood(nfd,5.5_dp+h,11.5_dp,2.7_dp,0.1_dp,1.2_dp), &
            ncf_loglikelihood(nfd,5.5_dp-h,11.5_dp,2.7_dp,0.1_dp,1.2_dp),h),tol,tol,failures)
        call check_close('ncf score dfd',s5(2),centered( &
            ncf_loglikelihood(nfd,5.5_dp,11.5_dp+h,2.7_dp,0.1_dp,1.2_dp), &
            ncf_loglikelihood(nfd,5.5_dp,11.5_dp-h,2.7_dp,0.1_dp,1.2_dp),h),tol,tol,failures)
        call check_close('ncf score nc',s5(3),centered( &
            ncf_loglikelihood(nfd,5.5_dp,11.5_dp,2.7_dp+h,0.1_dp,1.2_dp), &
            ncf_loglikelihood(nfd,5.5_dp,11.5_dp,2.7_dp-h,0.1_dp,1.2_dp),h),tol,tol,failures)
        call check_close('ncf score loc',s5(4),centered( &
            ncf_loglikelihood(nfd,5.5_dp,11.5_dp,2.7_dp,0.1_dp+h,1.2_dp), &
            ncf_loglikelihood(nfd,5.5_dp,11.5_dp,2.7_dp,0.1_dp-h,1.2_dp),h),tol,tol,failures)
        call check_close('ncf score scale',s5(5),centered( &
            ncf_loglikelihood(nfd,5.5_dp,11.5_dp,2.7_dp,0.1_dp,1.2_dp+h), &
            ncf_loglikelihood(nfd,5.5_dp,11.5_dp,2.7_dp,0.1_dp,1.2_dp-h),h),tol,tol,failures)
    end subroutine test_scores

    subroutine test_rng_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: ihd(4) = [0.6_dp,1.5_dp,2.6_dp,4.2_dp]
        real(dp), parameter :: kod(4) = [0.06_dp,0.13_dp,0.21_dp,0.31_dp]
        real(dp), parameter :: nxd(4) = [1.1_dp,3.0_dp,7.5_dp,12.0_dp]
        real(dp), parameter :: nfd(4) = [0.3_dp,0.8_dp,1.7_dp,3.2_dp]
        real(dp) :: actual, expected
        type(rng_state) :: reference, state
        type(fit_result) :: result

        call rng_seed(state,24680)
        call rng_seed(reference,24680)
        actual = irwinhall_rvs(state,5.0_dp,0.1_dp,1.2_dp)
        expected = irwinhall_ppf(rng_uniform(reference),5.0_dp,0.1_dp,1.2_dp)
        call check_close('irwinhall rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = ksone_rvs(state,20.0_dp,0.01_dp,1.1_dp)
        expected = ksone_ppf(rng_uniform(reference),20.0_dp,0.01_dp,1.1_dp)
        call check_close('ksone rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = ncx2_rvs(state,5.5_dp,2.5_dp,0.2_dp,1.1_dp)
        expected = ncx2_ppf(rng_uniform(reference),5.5_dp,2.5_dp,0.2_dp,1.1_dp)
        call check_close('ncx2 rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = ncf_rvs(state,5.5_dp,11.5_dp,2.7_dp,0.1_dp,1.2_dp)
        expected = ncf_ppf(rng_uniform(reference),5.5_dp,11.5_dp,2.7_dp,0.1_dp,1.2_dp)
        call check_close('ncf rvs',actual,expected,0.0_dp,0.0_dp,failures)

        call rng_seed(state,13579)
        call rng_seed(reference,13579)
        actual = ksone_rvs(state,2.5_dp)
        call check_true('invalid ksone rvs NaN',ieee_is_nan(actual),failures)
        call check_close('invalid rvs preserves state',rng_uniform(state),rng_uniform(reference), &
            0.0_dp,0.0_dp,failures)

        call irwinhall_fit(ihd,[5.0_dp,0.1_dp,1.2_dp],[5.0_dp,0.1_dp,1.2_dp],result)
        call check_true('irwinhall fixed fit',result%success,failures)
        call ksone_fit(kod,[20.0_dp,0.01_dp,1.1_dp],[20.0_dp,0.01_dp,1.1_dp],result)
        call check_true('ksone fixed fit',result%success,failures)
        call ncx2_fit(nxd,[5.5_dp,2.5_dp,0.2_dp,1.1_dp],[5.5_dp,2.5_dp,0.2_dp,1.1_dp],result)
        call check_true('ncx2 fixed fit',result%success,failures)
        call ncf_fit(nfd,[5.5_dp,11.5_dp,2.7_dp,0.1_dp,1.2_dp], &
            [5.5_dp,11.5_dp,2.7_dp,0.1_dp,1.2_dp],result)
        call check_true('ncf fixed fit',result%success,failures)

        call check_close('C ABI irwinhall pdf',scifort_irwinhall_pdf_f64(2.3_dp,5.0_dp,0.0_dp,1.0_dp), &
            irwinhall_pdf(2.3_dp,5.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI irwinhall cdf',scifort_irwinhall_cdf_f64(2.3_dp,5.0_dp,0.0_dp,1.0_dp), &
            irwinhall_cdf(2.3_dp,5.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI irwinhall ppf',scifort_irwinhall_ppf_f64(0.4_dp,5.0_dp,0.0_dp,1.0_dp), &
            irwinhall_ppf(0.4_dp,5.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI ksone pdf',scifort_ksone_pdf_f64(0.17_dp,20.0_dp,0.0_dp,1.0_dp), &
            ksone_pdf(0.17_dp,20.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI ksone cdf',scifort_ksone_cdf_f64(0.17_dp,20.0_dp,0.0_dp,1.0_dp), &
            ksone_cdf(0.17_dp,20.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI ksone ppf',scifort_ksone_ppf_f64(0.4_dp,20.0_dp,0.0_dp,1.0_dp), &
            ksone_ppf(0.4_dp,20.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI ncx2 pdf',scifort_ncx2_pdf_f64(8.0_dp,5.0_dp,3.0_dp,0.0_dp,1.0_dp), &
            ncx2_pdf(8.0_dp,5.0_dp,3.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI ncx2 cdf',scifort_ncx2_cdf_f64(8.0_dp,5.0_dp,3.0_dp,0.0_dp,1.0_dp), &
            ncx2_cdf(8.0_dp,5.0_dp,3.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI ncx2 ppf',scifort_ncx2_ppf_f64(0.4_dp,5.0_dp,3.0_dp,0.0_dp,1.0_dp), &
            ncx2_ppf(0.4_dp,5.0_dp,3.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI ncf pdf',scifort_ncf_pdf_f64(1.7_dp,5.0_dp,12.0_dp,3.0_dp,0.0_dp,1.0_dp), &
            ncf_pdf(1.7_dp,5.0_dp,12.0_dp,3.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI ncf cdf',scifort_ncf_cdf_f64(1.7_dp,5.0_dp,12.0_dp,3.0_dp,0.0_dp,1.0_dp), &
            ncf_cdf(1.7_dp,5.0_dp,12.0_dp,3.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI ncf ppf',scifort_ncf_ppf_f64(0.4_dp,5.0_dp,12.0_dp,3.0_dp,0.0_dp,1.0_dp), &
            ncf_ppf(0.4_dp,5.0_dp,12.0_dp,3.0_dp),0.0_dp,0.0_dp,failures)
    end subroutine test_rng_fit_c_api

    pure function centered(fp, fm, h) result(value)
        real(dp), intent(in) :: fp !! function value at positive perturbation
        real(dp), intent(in) :: fm !! function value at negative perturbation
        real(dp), intent(in) :: h !! perturbation magnitude
        real(dp) :: value
        value = (fp - fm) / (2.0_dp*h)
    end function centered

end program test_irwinhall_ksone_ncx2_ncf
