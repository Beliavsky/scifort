! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_dpareto_vonmises_kstwobign
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_c_api
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state, rng_uniform
    use scifort_stats
    use test_checks, only : check_close, check_true
    use test_dpareto_vonmises_kstwobign_reference
    implicit none

    integer :: failures

    failures = 0
    call test_reference(failures)
    call test_identities_endpoints(failures)
    call test_scores(failures)
    call test_rng_fit_c_api(failures)

    if (failures == 0) then
        print '(a)', 'test_dpareto_vonmises_kstwobign: PASS'
    else
        print '(a,1x,i0)', 'test_dpareto_vonmises_kstwobign: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        integer :: i

        do i = 1, size(reference_probs)
            call check_close('dpareto pdf',dpareto_lognorm_pdf(dp_x(i),0.2_dp,0.7_dp,1.3_dp,2.1_dp), &
                dp_pdf(i),3.0e-13_dp,3.0e-12_dp,failures)
            call check_close('dpareto logpdf',dpareto_lognorm_logpdf(dp_x(i),0.2_dp,0.7_dp,1.3_dp,2.1_dp), &
                dp_logpdf(i),3.0e-13_dp,3.0e-12_dp,failures)
            call check_close('dpareto cdf',dpareto_lognorm_cdf(dp_x(i),0.2_dp,0.7_dp,1.3_dp,2.1_dp), &
                dp_cdf(i),3.0e-13_dp,3.0e-12_dp,failures)
            call check_close('dpareto sf',dpareto_lognorm_sf(dp_x(i),0.2_dp,0.7_dp,1.3_dp,2.1_dp), &
                dp_sf(i),3.0e-13_dp,3.0e-12_dp,failures)
            call check_close('dpareto logcdf',dpareto_lognorm_logcdf(dp_x(i),0.2_dp,0.7_dp,1.3_dp,2.1_dp), &
                dp_logcdf(i),3.0e-12_dp,3.0e-12_dp,failures)
            call check_close('dpareto logsf',dpareto_lognorm_logsf(dp_x(i),0.2_dp,0.7_dp,1.3_dp,2.1_dp), &
                dp_logsf(i),3.0e-12_dp,3.0e-12_dp,failures)
            call check_close('dpareto ppf',dpareto_lognorm_ppf(reference_probs(i),0.2_dp,0.7_dp,1.3_dp,2.1_dp), &
                dp_ppf(i),2.0e-10_dp,2.0e-10_dp,failures)
            call check_close('dpareto isf',dpareto_lognorm_isf(reference_probs(i),0.2_dp,0.7_dp,1.3_dp,2.1_dp), &
                dp_isf(i),2.0e-10_dp,2.0e-10_dp,failures)

            call check_close('vonmises pdf',vonmises_pdf(vm_x(i),3.0_dp),vm_pdf(i), &
                3.0e-14_dp,3.0e-13_dp,failures)
            call check_close('vonmises logpdf',vonmises_logpdf(vm_x(i),3.0_dp),vm_logpdf(i), &
                3.0e-14_dp,3.0e-13_dp,failures)
            call check_close('vonmises cdf',vonmises_cdf(vm_x(i),3.0_dp),vm_cdf(i), &
                3.0e-13_dp,3.0e-12_dp,failures)
            call check_close('vonmises sf',vonmises_sf(vm_x(i),3.0_dp),vm_sf(i), &
                3.0e-13_dp,3.0e-12_dp,failures)
            call check_close('vonmises logcdf',vonmises_logcdf(vm_x(i),3.0_dp),vm_logcdf(i), &
                3.0e-12_dp,3.0e-12_dp,failures)
            call check_close('vonmises logsf',vonmises_logsf(vm_x(i),3.0_dp),vm_logsf(i), &
                3.0e-12_dp,3.0e-12_dp,failures)
            call check_close('vonmises ppf',vonmises_ppf(reference_probs(i),3.0_dp),vm_ppf(i), &
                5.0e-12_dp,5.0e-11_dp,failures)
            call check_close('vonmises isf',vonmises_isf(reference_probs(i),3.0_dp),vm_isf(i), &
                5.0e-12_dp,5.0e-11_dp,failures)

            call check_close('vonmises_line pdf',vonmises_line_pdf(vl_x(i),3.0_dp),vl_pdf(i), &
                3.0e-14_dp,3.0e-13_dp,failures)
            call check_close('vonmises_line logpdf',vonmises_line_logpdf(vl_x(i),3.0_dp),vl_logpdf(i), &
                3.0e-14_dp,3.0e-13_dp,failures)
            call check_close('vonmises_line cdf',vonmises_line_cdf(vl_x(i),3.0_dp),vl_cdf(i), &
                3.0e-13_dp,3.0e-12_dp,failures)
            call check_close('vonmises_line sf',vonmises_line_sf(vl_x(i),3.0_dp),vl_sf(i), &
                3.0e-13_dp,3.0e-12_dp,failures)
            call check_close('vonmises_line logcdf',vonmises_line_logcdf(vl_x(i),3.0_dp),vl_logcdf(i), &
                3.0e-12_dp,3.0e-12_dp,failures)
            call check_close('vonmises_line logsf',vonmises_line_logsf(vl_x(i),3.0_dp),vl_logsf(i), &
                3.0e-12_dp,3.0e-12_dp,failures)
            call check_close('vonmises_line ppf',vonmises_line_ppf(reference_probs(i),3.0_dp),vl_ppf(i), &
                5.0e-12_dp,5.0e-11_dp,failures)
            call check_close('vonmises_line isf',vonmises_line_isf(reference_probs(i),3.0_dp),vl_isf(i), &
                5.0e-12_dp,5.0e-11_dp,failures)

            call check_close('kstwobign pdf',kstwobign_pdf(ks_x(i)),ks_pdf(i), &
                3.0e-12_dp,3.0e-11_dp,failures)
            call check_close('kstwobign logpdf',kstwobign_logpdf(ks_x(i)),ks_logpdf(i), &
                3.0e-12_dp,3.0e-11_dp,failures)
            call check_close('kstwobign cdf',kstwobign_cdf(ks_x(i)),ks_cdf(i), &
                3.0e-13_dp,3.0e-12_dp,failures)
            call check_close('kstwobign sf',kstwobign_sf(ks_x(i)),ks_sf(i), &
                3.0e-13_dp,3.0e-12_dp,failures)
            call check_close('kstwobign logcdf',kstwobign_logcdf(ks_x(i)),ks_logcdf(i), &
                3.0e-12_dp,3.0e-12_dp,failures)
            call check_close('kstwobign logsf',kstwobign_logsf(ks_x(i)),ks_logsf(i), &
                3.0e-12_dp,3.0e-12_dp,failures)
            call check_close('kstwobign ppf',kstwobign_ppf(reference_probs(i)),ks_ppf(i), &
                3.0e-12_dp,3.0e-11_dp,failures)
            call check_close('kstwobign isf',kstwobign_isf(reference_probs(i)),ks_isf(i), &
                3.0e-12_dp,3.0e-11_dp,failures)
        end do
    end subroutine test_reference

    subroutine test_identities_endpoints(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: pvals(5) = [1.0e-10_dp,1.0e-5_dp,0.2_dp,0.8_dp,1.0_dp-1.0e-10_dp]
        real(dp), parameter :: xvals(4) = [0.2_dp,0.6_dp,1.4_dp,3.0_dp]
        real(dp) :: pi
        integer :: i

        pi = acos(-1.0_dp)
        do i = 1, size(xvals)
            call check_close('symmetric dpareto reciprocal', &
                dpareto_lognorm_cdf(xvals(i),0.0_dp,0.8_dp,1.5_dp,1.5_dp), &
                dpareto_lognorm_sf(1.0_dp/xvals(i),0.0_dp,0.8_dp,1.5_dp,1.5_dp), &
                3.0e-13_dp,3.0e-12_dp,failures)
        end do
        do i = 1, size(vm_x)
            call check_close('vonmises line/circular cdf',vonmises_line_cdf(vm_x(i),3.0_dp), &
                vonmises_cdf(vm_x(i),3.0_dp),3.0e-13_dp,3.0e-12_dp,failures)
            call check_close('vonmises line/circular pdf',vonmises_line_pdf(vm_x(i),3.0_dp), &
                vonmises_pdf(vm_x(i),3.0_dp),3.0e-14_dp,3.0e-13_dp,failures)
        end do
        call check_close('vonmises periodic cdf',vonmises_cdf(0.4_dp+2.0_dp*pi,2.0_dp), &
            vonmises_cdf(0.4_dp,2.0_dp)+1.0_dp,3.0e-13_dp,3.0e-12_dp,failures)
        call check_close('vonmises line uniform pdf',vonmises_line_pdf(0.4_dp,0.0_dp), &
            1.0_dp/(2.0_dp*pi),2.0e-15_dp,2.0e-15_dp,failures)
        call check_close('vonmises line uniform cdf',vonmises_line_cdf(0.4_dp,0.0_dp), &
            (0.4_dp+pi)/(2.0_dp*pi),2.0e-15_dp,2.0e-15_dp,failures)

        do i = 1, size(pvals)
            call check_close('dpareto cdf-ppf', &
                dpareto_lognorm_cdf(dpareto_lognorm_ppf(pvals(i),0.2_dp,0.7_dp,1.3_dp,2.1_dp), &
                    0.2_dp,0.7_dp,1.3_dp,2.1_dp),pvals(i),2.0e-11_dp,2.0e-9_dp,failures)
            call check_close('vonmises cdf-ppf',vonmises_cdf(vonmises_ppf(pvals(i),3.0_dp),3.0_dp), &
                pvals(i),2.0e-11_dp,2.0e-9_dp,failures)
            call check_close('vonmises line cdf-ppf', &
                vonmises_line_cdf(vonmises_line_ppf(pvals(i),3.0_dp),3.0_dp), &
                pvals(i),2.0e-11_dp,2.0e-9_dp,failures)
            call check_close('kstwobign cdf-ppf',kstwobign_cdf(kstwobign_ppf(pvals(i))), &
                pvals(i),2.0e-11_dp,2.0e-9_dp,failures)
            call check_close('kstwobign sf-isf',kstwobign_sf(kstwobign_isf(pvals(i))), &
                pvals(i),2.0e-11_dp,2.0e-9_dp,failures)
        end do
        call check_close('kstwobign complement below switch',kstwobign_cdf(0.81_dp)+kstwobign_sf(0.81_dp), &
            1.0_dp,2.0e-15_dp,2.0e-15_dp,failures)
        call check_close('kstwobign complement above switch',kstwobign_cdf(0.83_dp)+kstwobign_sf(0.83_dp), &
            1.0_dp,2.0e-15_dp,2.0e-15_dp,failures)
        call check_close('dpareto lower cdf',dpareto_lognorm_cdf(0.0_dp,0.2_dp,0.7_dp,1.3_dp,2.1_dp), &
            0.0_dp,0.0_dp,0.0_dp,failures)
        call check_true('dpareto ppf one infinite', &
            .not. ieee_is_finite(dpareto_lognorm_ppf(1.0_dp,0.2_dp,0.7_dp,1.3_dp,2.1_dp)),failures)
        call check_true('dpareto invalid shape', &
            ieee_is_nan(dpareto_lognorm_pdf(1.0_dp,0.2_dp,0.0_dp,1.3_dp,2.1_dp)),failures)
        call check_true('vonmises invalid shape',ieee_is_nan(vonmises_pdf(0.0_dp,-1.0_dp)),failures)
        call check_true('kstwobign ppf one infinite',.not. ieee_is_finite(kstwobign_ppf(1.0_dp)),failures)
    end subroutine test_identities_endpoints

    subroutine test_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: h = 1.0e-6_dp, tol = 3.0e-4_dp
        real(dp), parameter :: dd(4) = [0.35_dp,0.8_dp,1.6_dp,3.0_dp]
        real(dp), parameter :: vd(4) = [-0.8_dp,-0.1_dp,0.5_dp,1.1_dp]
        real(dp), parameter :: kd(4) = [0.35_dp,0.7_dp,1.1_dp,1.7_dp]
        real(dp) :: s2(2), s3(3), s6(6)

        s6 = dpareto_lognorm_score(dd,0.2_dp,0.7_dp,1.3_dp,2.1_dp,0.1_dp,1.2_dp)
        call check_close('dpareto score u',s6(1),centered( &
            dpareto_lognorm_loglikelihood(dd,0.2_dp+h,0.7_dp,1.3_dp,2.1_dp,0.1_dp,1.2_dp), &
            dpareto_lognorm_loglikelihood(dd,0.2_dp-h,0.7_dp,1.3_dp,2.1_dp,0.1_dp,1.2_dp),h),tol,tol,failures)
        call check_close('dpareto score s',s6(2),centered( &
            dpareto_lognorm_loglikelihood(dd,0.2_dp,0.7_dp+h,1.3_dp,2.1_dp,0.1_dp,1.2_dp), &
            dpareto_lognorm_loglikelihood(dd,0.2_dp,0.7_dp-h,1.3_dp,2.1_dp,0.1_dp,1.2_dp),h),tol,tol,failures)
        call check_close('dpareto score a',s6(3),centered( &
            dpareto_lognorm_loglikelihood(dd,0.2_dp,0.7_dp,1.3_dp+h,2.1_dp,0.1_dp,1.2_dp), &
            dpareto_lognorm_loglikelihood(dd,0.2_dp,0.7_dp,1.3_dp-h,2.1_dp,0.1_dp,1.2_dp),h),tol,tol,failures)
        call check_close('dpareto score b',s6(4),centered( &
            dpareto_lognorm_loglikelihood(dd,0.2_dp,0.7_dp,1.3_dp,2.1_dp+h,0.1_dp,1.2_dp), &
            dpareto_lognorm_loglikelihood(dd,0.2_dp,0.7_dp,1.3_dp,2.1_dp-h,0.1_dp,1.2_dp),h),tol,tol,failures)
        call check_close('dpareto score loc',s6(5),centered( &
            dpareto_lognorm_loglikelihood(dd,0.2_dp,0.7_dp,1.3_dp,2.1_dp,0.1_dp+h,1.2_dp), &
            dpareto_lognorm_loglikelihood(dd,0.2_dp,0.7_dp,1.3_dp,2.1_dp,0.1_dp-h,1.2_dp),h),tol,tol,failures)
        call check_close('dpareto score scale',s6(6),centered( &
            dpareto_lognorm_loglikelihood(dd,0.2_dp,0.7_dp,1.3_dp,2.1_dp,0.1_dp,1.2_dp+h), &
            dpareto_lognorm_loglikelihood(dd,0.2_dp,0.7_dp,1.3_dp,2.1_dp,0.1_dp,1.2_dp-h),h),tol,tol,failures)

        s3 = vonmises_score(vd,2.5_dp,0.1_dp,1.2_dp)
        call check_score3('vonmises',s3, &
            centered(vonmises_loglikelihood(vd,2.5_dp+h,0.1_dp,1.2_dp), &
                vonmises_loglikelihood(vd,2.5_dp-h,0.1_dp,1.2_dp),h), &
            centered(vonmises_loglikelihood(vd,2.5_dp,0.1_dp+h,1.2_dp), &
                vonmises_loglikelihood(vd,2.5_dp,0.1_dp-h,1.2_dp),h), &
            centered(vonmises_loglikelihood(vd,2.5_dp,0.1_dp,1.2_dp+h), &
                vonmises_loglikelihood(vd,2.5_dp,0.1_dp,1.2_dp-h),h),tol,failures)

        s3 = vonmises_line_score(vd,2.5_dp,0.1_dp,1.2_dp)
        call check_score3('vonmises_line',s3, &
            centered(vonmises_line_loglikelihood(vd,2.5_dp+h,0.1_dp,1.2_dp), &
                vonmises_line_loglikelihood(vd,2.5_dp-h,0.1_dp,1.2_dp),h), &
            centered(vonmises_line_loglikelihood(vd,2.5_dp,0.1_dp+h,1.2_dp), &
                vonmises_line_loglikelihood(vd,2.5_dp,0.1_dp-h,1.2_dp),h), &
            centered(vonmises_line_loglikelihood(vd,2.5_dp,0.1_dp,1.2_dp+h), &
                vonmises_line_loglikelihood(vd,2.5_dp,0.1_dp,1.2_dp-h),h),tol,failures)

        s2 = kstwobign_score(kd,0.1_dp,1.2_dp)
        call check_close('kstwobign score loc',s2(1),centered( &
            kstwobign_loglikelihood(kd,0.1_dp+h,1.2_dp), &
            kstwobign_loglikelihood(kd,0.1_dp-h,1.2_dp),h),tol,tol,failures)
        call check_close('kstwobign score scale',s2(2),centered( &
            kstwobign_loglikelihood(kd,0.1_dp,1.2_dp+h), &
            kstwobign_loglikelihood(kd,0.1_dp,1.2_dp-h),h),tol,tol,failures)
    end subroutine test_scores

    subroutine test_rng_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: dd(4) = [0.35_dp,0.8_dp,1.6_dp,3.0_dp]
        real(dp), parameter :: vd(4) = [-0.8_dp,-0.1_dp,0.5_dp,1.1_dp]
        real(dp), parameter :: kd(4) = [0.35_dp,0.7_dp,1.1_dp,1.7_dp]
        real(dp) :: actual, expected
        type(rng_state) :: reference, state
        type(fit_result) :: result

        call rng_seed(state,24680)
        call rng_seed(reference,24680)
        actual = dpareto_lognorm_rvs(state,0.2_dp,0.7_dp,1.3_dp,2.1_dp,0.1_dp,1.2_dp)
        expected = dpareto_lognorm_ppf(rng_uniform(reference),0.2_dp,0.7_dp,1.3_dp,2.1_dp,0.1_dp,1.2_dp)
        call check_close('dpareto rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = vonmises_rvs(state,2.5_dp,0.1_dp,1.2_dp)
        expected = vonmises_ppf(rng_uniform(reference),2.5_dp,0.1_dp,1.2_dp)
        call check_close('vonmises rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = vonmises_line_rvs(state,2.5_dp,0.1_dp,1.2_dp)
        expected = vonmises_line_ppf(rng_uniform(reference),2.5_dp,0.1_dp,1.2_dp)
        call check_close('vonmises_line rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = kstwobign_rvs(state,0.1_dp,1.2_dp)
        expected = kstwobign_ppf(rng_uniform(reference),0.1_dp,1.2_dp)
        call check_close('kstwobign rvs',actual,expected,0.0_dp,0.0_dp,failures)

        call rng_seed(state,13579)
        call rng_seed(reference,13579)
        actual = dpareto_lognorm_rvs(state,0.2_dp,0.0_dp,1.3_dp,2.1_dp)
        call check_true('invalid dpareto rvs NaN',ieee_is_nan(actual),failures)
        call check_close('invalid rvs preserves state',rng_uniform(state),rng_uniform(reference), &
            0.0_dp,0.0_dp,failures)

        call dpareto_lognorm_fit(dd,[0.2_dp,0.7_dp,1.3_dp,2.1_dp,0.1_dp,1.2_dp], &
            [0.2_dp,0.7_dp,1.3_dp,2.1_dp,0.1_dp,1.2_dp],result)
        call check_true('dpareto fixed fit',result%success,failures)
        call vonmises_fit(vd,[2.5_dp,0.1_dp,1.2_dp],[2.5_dp,0.1_dp,1.2_dp],result)
        call check_true('vonmises fixed fit',result%success,failures)
        call vonmises_line_fit(vd,[2.5_dp,0.1_dp,1.2_dp],[2.5_dp,0.1_dp,1.2_dp],result)
        call check_true('vonmises_line fixed fit',result%success,failures)
        call kstwobign_fit(kd,[0.1_dp,1.2_dp],[0.1_dp,1.2_dp],result)
        call check_true('kstwobign fixed fit',result%success,failures)

        call check_close('C ABI dpareto pdf', &
            scifort_dpareto_lognorm_pdf_f64(1.2_dp,0.2_dp,0.7_dp,1.3_dp,2.1_dp,0.1_dp,1.2_dp), &
            dpareto_lognorm_pdf(1.2_dp,0.2_dp,0.7_dp,1.3_dp,2.1_dp,0.1_dp,1.2_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI dpareto cdf', &
            scifort_dpareto_lognorm_cdf_f64(1.2_dp,0.2_dp,0.7_dp,1.3_dp,2.1_dp,0.1_dp,1.2_dp), &
            dpareto_lognorm_cdf(1.2_dp,0.2_dp,0.7_dp,1.3_dp,2.1_dp,0.1_dp,1.2_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI dpareto ppf', &
            scifort_dpareto_lognorm_ppf_f64(0.8_dp,0.2_dp,0.7_dp,1.3_dp,2.1_dp,0.1_dp,1.2_dp), &
            dpareto_lognorm_ppf(0.8_dp,0.2_dp,0.7_dp,1.3_dp,2.1_dp,0.1_dp,1.2_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI vonmises pdf',scifort_vonmises_pdf_f64(0.4_dp,2.5_dp,0.1_dp,1.2_dp), &
            vonmises_pdf(0.4_dp,2.5_dp,0.1_dp,1.2_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI vonmises cdf',scifort_vonmises_cdf_f64(0.4_dp,2.5_dp,0.1_dp,1.2_dp), &
            vonmises_cdf(0.4_dp,2.5_dp,0.1_dp,1.2_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI vonmises ppf',scifort_vonmises_ppf_f64(0.8_dp,2.5_dp,0.1_dp,1.2_dp), &
            vonmises_ppf(0.8_dp,2.5_dp,0.1_dp,1.2_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI vonmises_line pdf', &
            scifort_vonmises_line_pdf_f64(0.4_dp,2.5_dp,0.1_dp,1.2_dp), &
            vonmises_line_pdf(0.4_dp,2.5_dp,0.1_dp,1.2_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI vonmises_line cdf', &
            scifort_vonmises_line_cdf_f64(0.4_dp,2.5_dp,0.1_dp,1.2_dp), &
            vonmises_line_cdf(0.4_dp,2.5_dp,0.1_dp,1.2_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI vonmises_line ppf', &
            scifort_vonmises_line_ppf_f64(0.8_dp,2.5_dp,0.1_dp,1.2_dp), &
            vonmises_line_ppf(0.8_dp,2.5_dp,0.1_dp,1.2_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI kstwobign pdf',scifort_kstwobign_pdf_f64(1.0_dp,0.1_dp,1.2_dp), &
            kstwobign_pdf(1.0_dp,0.1_dp,1.2_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI kstwobign cdf',scifort_kstwobign_cdf_f64(1.0_dp,0.1_dp,1.2_dp), &
            kstwobign_cdf(1.0_dp,0.1_dp,1.2_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI kstwobign ppf',scifort_kstwobign_ppf_f64(0.8_dp,0.1_dp,1.2_dp), &
            kstwobign_ppf(0.8_dp,0.1_dp,1.2_dp),0.0_dp,0.0_dp,failures)
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
        value = (plus-minus)/(2.0_dp*h)
    end function centered

end program test_dpareto_vonmises_kstwobign
