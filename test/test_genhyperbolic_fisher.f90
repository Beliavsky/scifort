! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_genhyperbolic_fisher
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_c_api
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state, rng_uniform
    use scifort_special, only : log_besselk_scaled
    use scifort_stats
    use test_genhyperbolic_fisher_reference
    use test_checks, only : check_close, check_true
    implicit none

    integer :: failures

    failures = 0
    call test_reference(failures)
    call test_identities_endpoints(failures)
    call test_scores(failures)
    call test_rng_fit_c_api(failures)

    if (failures == 0) then
        print '(a)', 'test_genhyperbolic_fisher: PASS'
    else
        print '(a,1x,i0)', 'test_genhyperbolic_fisher: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        integer :: i

        do i = 1, size(gh_x)
            call check_close('GH pdf', genhyperbolic_pdf(gh_x(i),1.2_dp,2.5_dp,0.7_dp), &
                gh_pdf(i),3.0e-13_dp,3.0e-11_dp,failures)
            call check_close('GH logpdf',genhyperbolic_logpdf(gh_x(i),1.2_dp,2.5_dp,0.7_dp), &
                gh_logpdf(i),3.0e-12_dp,3.0e-11_dp,failures)
            call check_close('GH cdf',genhyperbolic_cdf(gh_x(i),1.2_dp,2.5_dp,0.7_dp), &
                gh_cdf(i),3.0e-12_dp,3.0e-10_dp,failures)
            call check_close('GH sf',genhyperbolic_sf(gh_x(i),1.2_dp,2.5_dp,0.7_dp), &
                gh_sf(i),3.0e-12_dp,3.0e-10_dp,failures)
            call check_close('GH logcdf',genhyperbolic_logcdf(gh_x(i),1.2_dp,2.5_dp,0.7_dp), &
                gh_logcdf(i),3.0e-10_dp,3.0e-9_dp,failures)
            call check_close('GH logsf',genhyperbolic_logsf(gh_x(i),1.2_dp,2.5_dp,0.7_dp), &
                gh_logsf(i),3.0e-10_dp,3.0e-9_dp,failures)
        end do
        do i = 1, size(reference_probs)
            call check_close('GH ppf',genhyperbolic_ppf(reference_probs(i),1.2_dp,2.5_dp,0.7_dp), &
                gh_ppf(i),2.0e-8_dp,3.0e-9_dp,failures)
            call check_close('GH isf',genhyperbolic_isf(reference_probs(i),1.2_dp,2.5_dp,0.7_dp), &
                gh_isf(i),2.0e-8_dp,3.0e-9_dp,failures)
        end do

        do i = 1, size(fh_x)
            call check_close('Fisher pmf',nchypergeom_fisher_pmf(fh_x(i),30.0_dp,12.0_dp,7.0_dp,2.3_dp), &
                fh_pmf(i),3.0e-15_dp,3.0e-13_dp,failures)
            if (.not. ieee_is_nan(fh_logpmf(i))) then
                call check_close('Fisher logpmf',nchypergeom_fisher_logpmf(fh_x(i),30.0_dp,12.0_dp,7.0_dp,2.3_dp), &
                    fh_logpmf(i),3.0e-13_dp,3.0e-13_dp,failures)
            end if
            call check_close('Fisher cdf',nchypergeom_fisher_cdf(fh_x(i),30.0_dp,12.0_dp,7.0_dp,2.3_dp), &
                fh_cdf(i),3.0e-15_dp,3.0e-13_dp,failures)
            call check_close('Fisher sf',nchypergeom_fisher_sf(fh_x(i),30.0_dp,12.0_dp,7.0_dp,2.3_dp), &
                fh_sf(i),3.0e-15_dp,3.0e-13_dp,failures)
            call check_close('Fisher logcdf',nchypergeom_fisher_logcdf(fh_x(i),30.0_dp,12.0_dp,7.0_dp,2.3_dp), &
                fh_logcdf(i),3.0e-13_dp,3.0e-13_dp,failures)
            call check_close('Fisher logsf',nchypergeom_fisher_logsf(fh_x(i),30.0_dp,12.0_dp,7.0_dp,2.3_dp), &
                fh_logsf(i),3.0e-13_dp,3.0e-13_dp,failures)
        end do
        do i = 1, size(reference_probs)
            call check_close('Fisher ppf',nchypergeom_fisher_ppf(reference_probs(i),30.0_dp,12.0_dp,7.0_dp,2.3_dp), &
                fh_ppf(i),0.0_dp,0.0_dp,failures)
            call check_close('Fisher isf',nchypergeom_fisher_isf(reference_probs(i),30.0_dp,12.0_dp,7.0_dp,2.3_dp), &
                fh_isf(i),0.0_dp,0.0_dp,failures)
        end do

        do i = 1, size(bk_scaled_nu)
            call check_close('scaled Bessel K',log_besselk_scaled(bk_scaled_nu(i),bk_scaled_x(i)), &
                bk_scaled_log(i),3.0e-12_dp,3.0e-12_dp,failures)
        end do
    end subroutine test_reference

    subroutine test_identities_endpoints(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        integer :: k
        real(dp) :: total

        call check_close('GH tails',genhyperbolic_cdf(1.3_dp,1.2_dp,2.5_dp,0.7_dp)+ &
            genhyperbolic_sf(1.3_dp,1.2_dp,2.5_dp,0.7_dp),1.0_dp,3.0e-13_dp,3.0e-13_dp,failures)
        call check_close('Fisher tails',nchypergeom_fisher_cdf(3,30,12,7,2.3_dp)+ &
            nchypergeom_fisher_sf(3,30,12,7,2.3_dp),1.0_dp,3.0e-15_dp,3.0e-14_dp,failures)

        do k = -4, 4
            call check_close('GH NIG pdf identity',genhyperbolic_pdf(real(k,dp)/2.0_dp,-0.5_dp,2.5_dp,0.7_dp), &
                norminvgauss_pdf(real(k,dp)/2.0_dp,2.5_dp,0.7_dp),3.0e-12_dp,3.0e-10_dp,failures)
            call check_close('GH NIG cdf identity',genhyperbolic_cdf(real(k,dp)/2.0_dp,-0.5_dp,2.5_dp,0.7_dp), &
                norminvgauss_cdf(real(k,dp)/2.0_dp,2.5_dp,0.7_dp),3.0e-11_dp,3.0e-9_dp,failures)
        end do

        total = 0.0_dp
        do k = 0, 7
            total = total + nchypergeom_fisher_pmf(k,30,12,7,2.3_dp)
            call check_close('Fisher hypergeom pmf identity',nchypergeom_fisher_pmf(k,30,12,7,1.0_dp), &
                hypergeom_pmf(k,30,12,7),3.0e-15_dp,3.0e-13_dp,failures)
            call check_close('Fisher hypergeom cdf identity',nchypergeom_fisher_cdf(k,30,12,7,1.0_dp), &
                hypergeom_cdf(k,30,12,7),3.0e-15_dp,3.0e-13_dp,failures)
        end do
        call check_close('Fisher normalization',total,1.0_dp,3.0e-15_dp,3.0e-14_dp,failures)
        call check_close('Fisher shifted',nchypergeom_fisher_cdf(3.5_dp,30.0_dp,12.0_dp,7.0_dp,2.3_dp,0.5_dp), &
            nchypergeom_fisher_cdf(3.0_dp,30.0_dp,12.0_dp,7.0_dp,2.3_dp),0.0_dp,0.0_dp,failures)

        call check_true('GH invalid open boundary',ieee_is_nan(genhyperbolic_pdf(0.0_dp,1.0_dp,2.0_dp,2.0_dp)),failures)
        call check_true('GH negative-p closed boundary', &
            ieee_is_finite(genhyperbolic_pdf(0.0_dp,-1.0_dp,2.0_dp,2.0_dp)),failures)
        call check_close('GH skewed lower log tail', &
            genhyperbolic_logcdf(-20.0_dp,-2.0_dp,2.5_dp,2.4_dp), &
            -108.05810925056922_dp,5.0e-9_dp,5.0e-11_dp,failures)
        call check_close('GH skewed upper log tail', &
            genhyperbolic_logsf(20.0_dp,2.0_dp,5.0_dp,-4.8_dp), &
            -198.27997111944788_dp,1.0e-7_dp,1.0e-9_dp,failures)

        call check_true('Fisher invalid odds',ieee_is_nan(nchypergeom_fisher_pmf(2,30,12,7,0.0_dp)),failures)
    end subroutine test_identities_endpoints

    subroutine test_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: h = 2.0e-6_dp
        real(dp), parameter :: gd(4) = [-1.2_dp,-0.2_dp,0.7_dp,2.1_dp]
        integer, parameter :: fd_data(5) = [1,2,3,4,5]
        real(dp) :: s5(5), s1(1), fd

        s5 = genhyperbolic_score(gd,1.2_dp,2.5_dp,0.7_dp,0.1_dp,1.2_dp)
        fd = centered(genhyperbolic_loglikelihood(gd,1.2_dp+h,2.5_dp,0.7_dp,0.1_dp,1.2_dp), &
            genhyperbolic_loglikelihood(gd,1.2_dp-h,2.5_dp,0.7_dp,0.1_dp,1.2_dp),h)
        call check_close('GH score p',s5(1),fd,8.0e-6_dp,8.0e-6_dp,failures)
        fd = centered(genhyperbolic_loglikelihood(gd,1.2_dp,2.5_dp+h,0.7_dp,0.1_dp,1.2_dp), &
            genhyperbolic_loglikelihood(gd,1.2_dp,2.5_dp-h,0.7_dp,0.1_dp,1.2_dp),h)
        call check_close('GH score a',s5(2),fd,8.0e-6_dp,8.0e-6_dp,failures)
        fd = centered(genhyperbolic_loglikelihood(gd,1.2_dp,2.5_dp,0.7_dp+h,0.1_dp,1.2_dp), &
            genhyperbolic_loglikelihood(gd,1.2_dp,2.5_dp,0.7_dp-h,0.1_dp,1.2_dp),h)
        call check_close('GH score b',s5(3),fd,8.0e-6_dp,8.0e-6_dp,failures)
        fd = centered(genhyperbolic_loglikelihood(gd,1.2_dp,2.5_dp,0.7_dp,0.1_dp+h,1.2_dp), &
            genhyperbolic_loglikelihood(gd,1.2_dp,2.5_dp,0.7_dp,0.1_dp-h,1.2_dp),h)
        call check_close('GH score loc',s5(4),fd,8.0e-6_dp,8.0e-6_dp,failures)
        fd = centered(genhyperbolic_loglikelihood(gd,1.2_dp,2.5_dp,0.7_dp,0.1_dp,1.2_dp+h), &
            genhyperbolic_loglikelihood(gd,1.2_dp,2.5_dp,0.7_dp,0.1_dp,1.2_dp-h),h)
        call check_close('GH score scale',s5(5),fd,8.0e-6_dp,8.0e-6_dp,failures)

        s1 = nchypergeom_fisher_score(fd_data,30,12,7,2.3_dp)
        fd = centered(nchypergeom_fisher_loglikelihood(fd_data,30,12,7,2.3_dp+h), &
            nchypergeom_fisher_loglikelihood(fd_data,30,12,7,2.3_dp-h),h)
        call check_close('Fisher score odds',s1(1),fd,2.0e-7_dp,2.0e-7_dp,failures)
    end subroutine test_scores

    subroutine test_rng_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: gd(4) = [-1.2_dp,-0.2_dp,0.7_dp,2.1_dp]
        integer, parameter :: fd_data(5) = [1,2,3,4,5]
        real(dp) :: actual, expected
        type(rng_state) :: state, reference
        type(fit_result) :: result

        call rng_seed(state,24680); call rng_seed(reference,24680)
        actual = genhyperbolic_rvs(state,1.2_dp,2.5_dp,0.7_dp)
        expected = genhyperbolic_ppf(rng_uniform(reference),1.2_dp,2.5_dp,0.7_dp)
        call check_close('GH rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = nchypergeom_fisher_rvs(state,30.0_dp,12.0_dp,7.0_dp,2.3_dp)
        expected = nchypergeom_fisher_ppf(rng_uniform(reference),30.0_dp,12.0_dp,7.0_dp,2.3_dp)
        call check_close('Fisher rvs',actual,expected,0.0_dp,0.0_dp,failures)

        call rng_seed(state,13579); call rng_seed(reference,13579)
        actual = genhyperbolic_rvs(state,1.0_dp,2.0_dp,2.0_dp)
        call check_true('invalid GH rvs NaN',ieee_is_nan(actual),failures)
        call check_close('invalid rvs preserves state',rng_uniform(state),rng_uniform(reference), &
            0.0_dp,0.0_dp,failures)

        call genhyperbolic_fit(gd,[1.2_dp,2.5_dp,0.7_dp,0.0_dp,1.0_dp], &
            [1.2_dp,2.5_dp,0.7_dp,0.0_dp,1.0_dp],result)
        call check_true('GH fixed fit',result%success,failures)
        call nchypergeom_fisher_fit(fd_data,[30.0_dp,12.0_dp,7.0_dp,2.3_dp,0.0_dp], &
            [30.0_dp,12.0_dp,7.0_dp,2.3_dp,0.0_dp],result)
        call check_true('Fisher fixed fit',result%success,failures)

        call check_close('C ABI GH pdf',scifort_genhyperbolic_pdf_f64(0.3_dp,1.2_dp,2.5_dp,0.7_dp,0.0_dp,1.0_dp), &
            genhyperbolic_pdf(0.3_dp,1.2_dp,2.5_dp,0.7_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI GH cdf',scifort_genhyperbolic_cdf_f64(0.3_dp,1.2_dp,2.5_dp,0.7_dp,0.0_dp,1.0_dp), &
            genhyperbolic_cdf(0.3_dp,1.2_dp,2.5_dp,0.7_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI GH ppf',scifort_genhyperbolic_ppf_f64(0.7_dp,1.2_dp,2.5_dp,0.7_dp,0.0_dp,1.0_dp), &
            genhyperbolic_ppf(0.7_dp,1.2_dp,2.5_dp,0.7_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI Fisher pmf', &
            scifort_nchypergeom_fisher_pmf_f64(2.0_dp,30.0_dp,12.0_dp,7.0_dp,2.3_dp,0.0_dp), &
            nchypergeom_fisher_pmf(2.0_dp,30.0_dp,12.0_dp,7.0_dp,2.3_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI Fisher cdf', &
            scifort_nchypergeom_fisher_cdf_f64(2.0_dp,30.0_dp,12.0_dp,7.0_dp,2.3_dp,0.0_dp), &
            nchypergeom_fisher_cdf(2.0_dp,30.0_dp,12.0_dp,7.0_dp,2.3_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI Fisher ppf', &
            scifort_nchypergeom_fisher_ppf_f64(0.7_dp,30.0_dp,12.0_dp,7.0_dp,2.3_dp,0.0_dp), &
            nchypergeom_fisher_ppf(0.7_dp,30.0_dp,12.0_dp,7.0_dp,2.3_dp),0.0_dp,0.0_dp,failures)
    end subroutine test_rng_fit_c_api

    pure function centered(plus_value,minus_value,h) result(value)
        real(dp), intent(in) :: plus_value !! function value at theta+h
        real(dp), intent(in) :: minus_value !! function value at theta-h
        real(dp), intent(in) :: h !! positive centered-difference step
        real(dp) :: value
        value = (plus_value-minus_value)/(2.0_dp*h)
    end function centered

end program test_genhyperbolic_fisher
