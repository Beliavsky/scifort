! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_gig_nig_skellam
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_c_api
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state, rng_uniform
    use scifort_special, only : besselk_log_derivative_x, log_besselk, log_besselk_order_derivative
    use scifort_stats
    use test_gig_nig_skellam_reference
    use test_checks, only : check_close, check_true
    implicit none

    integer :: failures

    failures = 0
    call test_reference(failures)
    call test_identities_endpoints(failures)
    call test_scores(failures)
    call test_rng_fit_c_api(failures)

    if (failures == 0) then
        print '(a)', 'test_gig_nig_skellam: PASS'
    else
        print '(a,1x,i0)', 'test_gig_nig_skellam: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        integer :: i
        do i=1,size(gig_x)
            call check_close('GIG pdf',geninvgauss_pdf(gig_x(i),1.2_dp,2.3_dp), &
                gig_pdf(i),2.0e-13_dp,2.0e-11_dp,failures)
            call check_close('GIG logpdf',geninvgauss_logpdf(gig_x(i),1.2_dp,2.3_dp), &
                gig_logpdf(i),2.0e-12_dp,2.0e-11_dp,failures)
            call check_close('GIG cdf',geninvgauss_cdf(gig_x(i),1.2_dp,2.3_dp), &
                gig_cdf(i),2.0e-12_dp,2.0e-10_dp,failures)
            call check_close('GIG sf',geninvgauss_sf(gig_x(i),1.2_dp,2.3_dp), &
                gig_sf(i),2.0e-12_dp,2.0e-10_dp,failures)
            call check_close('GIG logcdf',geninvgauss_logcdf(gig_x(i),1.2_dp,2.3_dp), &
                gig_logcdf(i),2.0e-10_dp,2.0e-10_dp,failures)
            if (i < size(gig_x)) then
                call check_close('GIG logsf',geninvgauss_logsf(gig_x(i),1.2_dp,2.3_dp), &
                    gig_logsf(i),2.0e-10_dp,2.0e-10_dp,failures)
            else
                ! Direct high-accuracy quadrature of the defining density.
                ! SciPy 1.17.0 differs here by about 6.3e-6 in log survival.
                call check_close('GIG logsf high precision', &
                    geninvgauss_logsf(gig_x(i),1.2_dp,2.3_dp), &
                    -15.328256967854234_dp,2.0e-13_dp,2.0e-13_dp,failures)
            end if

            call check_close('NIG pdf',norminvgauss_pdf(nig_x(i),2.5_dp,0.7_dp), &
                nig_pdf(i),3.0e-13_dp,2.0e-11_dp,failures)
            call check_close('NIG logpdf',norminvgauss_logpdf(nig_x(i),2.5_dp,0.7_dp), &
                nig_logpdf(i),3.0e-12_dp,2.0e-11_dp,failures)
            call check_close('NIG cdf',norminvgauss_cdf(nig_x(i),2.5_dp,0.7_dp), &
                nig_cdf(i),3.0e-12_dp,3.0e-9_dp,failures)
            call check_close('NIG sf',norminvgauss_sf(nig_x(i),2.5_dp,0.7_dp), &
                nig_sf(i),3.0e-12_dp,3.0e-9_dp,failures)
            call check_close('NIG logcdf',norminvgauss_logcdf(nig_x(i),2.5_dp,0.7_dp), &
                nig_logcdf(i),5.0e-7_dp,3.0e-8_dp,failures)
            call check_close('NIG logsf',norminvgauss_logsf(nig_x(i),2.5_dp,0.7_dp), &
                nig_logsf(i),5.0e-7_dp,3.0e-8_dp,failures)

            call check_close('Skellam pmf',skellam_pmf(sk_x(i),3.2_dp,1.7_dp), &
                sk_pmf(i),3.0e-14_dp,3.0e-13_dp,failures)
            call check_close('Skellam logpmf',skellam_logpmf(sk_x(i),3.2_dp,1.7_dp), &
                sk_logpmf(i),3.0e-13_dp,3.0e-13_dp,failures)
            call check_close('Skellam cdf',skellam_cdf(sk_x(i),3.2_dp,1.7_dp), &
                sk_cdf(i),3.0e-14_dp,3.0e-12_dp,failures)
            call check_close('Skellam sf',skellam_sf(sk_x(i),3.2_dp,1.7_dp), &
                sk_sf(i),3.0e-14_dp,3.0e-12_dp,failures)
            call check_close('Skellam logcdf',skellam_logcdf(sk_x(i),3.2_dp,1.7_dp), &
                sk_logcdf(i),3.0e-12_dp,3.0e-12_dp,failures)
            call check_close('Skellam logsf',skellam_logsf(sk_x(i),3.2_dp,1.7_dp), &
                sk_logsf(i),3.0e-12_dp,3.0e-12_dp,failures)
        end do
        do i=1,size(reference_probs)
            call check_close('GIG ppf',geninvgauss_ppf(reference_probs(i),1.2_dp,2.3_dp), &
                gig_ppf(i),2.0e-8_dp,2.0e-9_dp,failures)
            call check_close('GIG isf',geninvgauss_isf(reference_probs(i),1.2_dp,2.3_dp), &
                gig_isf(i),2.0e-8_dp,2.0e-9_dp,failures)
            call check_close('NIG ppf',norminvgauss_ppf(reference_probs(i),2.5_dp,0.7_dp), &
                nig_ppf(i),2.0e-8_dp,2.0e-9_dp,failures)
            call check_close('NIG isf',norminvgauss_isf(reference_probs(i),2.5_dp,0.7_dp), &
                nig_isf(i),2.0e-8_dp,2.0e-9_dp,failures)
            call check_close('Skellam ppf',skellam_ppf(reference_probs(i),3.2_dp,1.7_dp), &
                sk_ppf(i),0.0_dp,0.0_dp,failures)
            call check_close('Skellam isf',skellam_isf(reference_probs(i),3.2_dp,1.7_dp), &
                sk_isf(i),0.0_dp,0.0_dp,failures)
        end do
        do i=1,size(bk_nu)
            call check_close('log Bessel K',log_besselk(bk_nu(i),bk_x(i)), &
                bk_log(i),2.0e-13_dp,2.0e-13_dp,failures)
        end do
    end subroutine test_reference

    subroutine test_identities_endpoints(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        integer :: k
        real(dp) :: total

        call check_close('GIG tails',geninvgauss_cdf(2.0_dp,1.2_dp,2.3_dp)+ &
            geninvgauss_sf(2.0_dp,1.2_dp,2.3_dp),1.0_dp,2.0e-13_dp,2.0e-13_dp,failures)
        call check_close('NIG tails',norminvgauss_cdf(1.0_dp,2.5_dp,0.7_dp)+ &
            norminvgauss_sf(1.0_dp,2.5_dp,0.7_dp),1.0_dp,2.0e-13_dp,2.0e-13_dp,failures)
        call check_close('Skellam tails',skellam_cdf(2,3.2_dp,1.7_dp)+ &
            skellam_sf(2,3.2_dp,1.7_dp),1.0_dp,2.0e-14_dp,2.0e-14_dp,failures)

        total=0.0_dp
        do k=-30,35
            total=total+skellam_pmf(k,3.2_dp,1.7_dp)
        end do
        call check_close('Skellam normalization',total,1.0_dp,2.0e-14_dp,2.0e-14_dp,failures)
        call check_close('Skellam symmetry',skellam_pmf(4,3.2_dp,1.7_dp), &
            skellam_pmf(-4,1.7_dp,3.2_dp),2.0e-15_dp,2.0e-14_dp,failures)
        call check_close('Skellam shifted',skellam_cdf(2.5_dp,3.2_dp,1.7_dp,0.5_dp), &
            skellam_cdf(2.0_dp,3.2_dp,1.7_dp),0.0_dp,0.0_dp,failures)

        ! High-precision stress references for the shared Bessel-K kernel and
        ! direct log-tail quadrature.  Constants were evaluated independently
        ! with 50-80 digit arithmetic from the defining integral formulas.
        call check_close('Bessel K wide-order log',log_besselk(100.0_dp,1000.0_dp), &
            -998.234850360727953_dp,2.0e-12_dp,2.0e-14_dp,failures)
        call check_close('Bessel K wide-order nu derivative', &
            log_besselk_order_derivative(100.0_dp,1000.0_dp), &
            0.099784626620319048_dp,2.0e-11_dp,2.0e-11_dp,failures)
        call check_close('Bessel K wide-order x derivative', &
            besselk_log_derivative_x(100.0_dp,1000.0_dp), &
            -1.00548249467332335_dp,2.0e-12_dp,2.0e-12_dp,failures)
        call check_close('GIG deep logcdf',geninvgauss_logcdf(1.0e-3_dp,1.2_dp,2.3_dp), &
            -1163.75841193668282_dp,1.0e-8_dp,1.0e-11_dp,failures)
        call check_close('GIG deep logsf',geninvgauss_logsf(1000.0_dp,1.2_dp,2.3_dp), &
            -1147.17771592627017_dp,1.0e-8_dp,1.0e-11_dp,failures)
        call check_close('GIG difficult isf',geninvgauss_isf(1.0e-6_dp,-0.5_dp,0.01_dp), &
            1048.9614348032378_dp,1.0e-6_dp,1.0e-10_dp,failures)
        call check_close('NIG deep logcdf',norminvgauss_logcdf(-300.0_dp,2.5_dp,0.7_dp), &
            -967.784848851011178_dp,1.0e-8_dp,1.0e-11_dp,failures)
        call check_close('NIG deep logsf',norminvgauss_logsf(300.0_dp,2.5_dp,0.7_dp), &
            -547.210690919970967_dp,1.0e-8_dp,1.0e-11_dp,failures)
        call check_close('NIG difficult isf',norminvgauss_isf(1.0e-6_dp,1.0_dp,0.99_dp), &
            753.9034692822579_dp,1.0e-6_dp,1.0e-10_dp,failures)

        call check_true('GIG invalid b',ieee_is_nan(geninvgauss_pdf(1.0_dp,1.2_dp,0.0_dp)),failures)
        call check_true('NIG invalid shapes',ieee_is_nan(norminvgauss_pdf(0.0_dp,1.0_dp,1.0_dp)),failures)
        call check_true('Skellam invalid mean',ieee_is_nan(skellam_pmf(0,0.0_dp,1.0_dp)),failures)
    end subroutine test_identities_endpoints

    subroutine test_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: h=2.0e-6_dp
        real(dp), parameter :: gd(4)=[0.45_dp,0.8_dp,1.4_dp,2.8_dp]
        real(dp), parameter :: nd(4)=[-1.2_dp,-0.2_dp,0.7_dp,2.1_dp]
        real(dp), parameter :: sd(5)=[-2.0_dp,0.0_dp,1.0_dp,3.0_dp,5.0_dp]
        real(dp) :: s4(4),s2(2),fd

        s4=geninvgauss_score(gd,1.2_dp,2.3_dp,0.1_dp,1.1_dp)
        fd=centered(geninvgauss_loglikelihood(gd,1.2_dp+h,2.3_dp,0.1_dp,1.1_dp), &
            geninvgauss_loglikelihood(gd,1.2_dp-h,2.3_dp,0.1_dp,1.1_dp),h)
        call check_close('GIG score p',s4(1),fd,2.0e-6_dp,2.0e-6_dp,failures)
        fd=centered(geninvgauss_loglikelihood(gd,1.2_dp,2.3_dp+h,0.1_dp,1.1_dp), &
            geninvgauss_loglikelihood(gd,1.2_dp,2.3_dp-h,0.1_dp,1.1_dp),h)
        call check_close('GIG score b',s4(2),fd,2.0e-6_dp,2.0e-6_dp,failures)
        fd=centered(geninvgauss_loglikelihood(gd,1.2_dp,2.3_dp,0.1_dp+h,1.1_dp), &
            geninvgauss_loglikelihood(gd,1.2_dp,2.3_dp,0.1_dp-h,1.1_dp),h)
        call check_close('GIG score loc',s4(3),fd,2.0e-6_dp,2.0e-6_dp,failures)
        fd=centered(geninvgauss_loglikelihood(gd,1.2_dp,2.3_dp,0.1_dp,1.1_dp+h), &
            geninvgauss_loglikelihood(gd,1.2_dp,2.3_dp,0.1_dp,1.1_dp-h),h)
        call check_close('GIG score scale',s4(4),fd,2.0e-6_dp,2.0e-6_dp,failures)

        s4=norminvgauss_score(nd,2.5_dp,0.7_dp,0.1_dp,1.2_dp)
        fd=centered(norminvgauss_loglikelihood(nd,2.5_dp+h,0.7_dp,0.1_dp,1.2_dp), &
            norminvgauss_loglikelihood(nd,2.5_dp-h,0.7_dp,0.1_dp,1.2_dp),h)
        call check_close('NIG score a',s4(1),fd,3.0e-6_dp,3.0e-6_dp,failures)
        fd=centered(norminvgauss_loglikelihood(nd,2.5_dp,0.7_dp+h,0.1_dp,1.2_dp), &
            norminvgauss_loglikelihood(nd,2.5_dp,0.7_dp-h,0.1_dp,1.2_dp),h)
        call check_close('NIG score b',s4(2),fd,3.0e-6_dp,3.0e-6_dp,failures)
        fd=centered(norminvgauss_loglikelihood(nd,2.5_dp,0.7_dp,0.1_dp+h,1.2_dp), &
            norminvgauss_loglikelihood(nd,2.5_dp,0.7_dp,0.1_dp-h,1.2_dp),h)
        call check_close('NIG score loc',s4(3),fd,3.0e-6_dp,3.0e-6_dp,failures)
        fd=centered(norminvgauss_loglikelihood(nd,2.5_dp,0.7_dp,0.1_dp,1.2_dp+h), &
            norminvgauss_loglikelihood(nd,2.5_dp,0.7_dp,0.1_dp,1.2_dp-h),h)
        call check_close('NIG score scale',s4(4),fd,3.0e-6_dp,3.0e-6_dp,failures)

        s2=skellam_score(sd,3.2_dp,1.7_dp)
        fd=centered(skellam_loglikelihood(sd,3.2_dp+h,1.7_dp), &
            skellam_loglikelihood(sd,3.2_dp-h,1.7_dp),h)
        call check_close('Skellam score mu1',s2(1),fd,2.0e-6_dp,2.0e-6_dp,failures)
        fd=centered(skellam_loglikelihood(sd,3.2_dp,1.7_dp+h), &
            skellam_loglikelihood(sd,3.2_dp,1.7_dp-h),h)
        call check_close('Skellam score mu2',s2(2),fd,2.0e-6_dp,2.0e-6_dp,failures)
    end subroutine test_scores

    subroutine test_rng_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: gd(4)=[0.45_dp,0.8_dp,1.4_dp,2.8_dp]
        real(dp), parameter :: nd(4)=[-1.2_dp,-0.2_dp,0.7_dp,2.1_dp]
        real(dp), parameter :: sd(4)=[-2.0_dp,0.0_dp,2.0_dp,5.0_dp]
        real(dp) :: actual,expected
        type(rng_state) :: state,reference
        type(fit_result) :: result

        call rng_seed(state,24680); call rng_seed(reference,24680)
        actual=geninvgauss_rvs(state,1.2_dp,2.3_dp)
        expected=geninvgauss_ppf(rng_uniform(reference),1.2_dp,2.3_dp)
        call check_close('GIG rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual=norminvgauss_rvs(state,2.5_dp,0.7_dp)
        expected=norminvgauss_ppf(rng_uniform(reference),2.5_dp,0.7_dp)
        call check_close('NIG rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual=skellam_rvs(state,3.2_dp,1.7_dp)
        expected=skellam_ppf(rng_uniform(reference),3.2_dp,1.7_dp)
        call check_close('Skellam rvs',actual,expected,0.0_dp,0.0_dp,failures)

        call rng_seed(state,13579); call rng_seed(reference,13579)
        actual=norminvgauss_rvs(state,1.0_dp,1.0_dp)
        call check_true('invalid NIG rvs NaN',ieee_is_nan(actual),failures)
        call check_close('invalid rvs preserves state',rng_uniform(state),rng_uniform(reference), &
            0.0_dp,0.0_dp,failures)

        call geninvgauss_fit(gd,[1.2_dp,2.3_dp,0.0_dp,1.0_dp], &
            [1.2_dp,2.3_dp,0.0_dp,1.0_dp],result)
        call check_true('GIG fixed fit',result%success,failures)
        call norminvgauss_fit(nd,[2.5_dp,0.7_dp,0.0_dp,1.0_dp], &
            [2.5_dp,0.7_dp,0.0_dp,1.0_dp],result)
        call check_true('NIG fixed fit',result%success,failures)
        call skellam_fit(sd,[3.2_dp,1.7_dp,0.0_dp],[3.2_dp,1.7_dp,0.0_dp],result)
        call check_true('Skellam fixed fit',result%success,failures)

        call check_close('C ABI GIG pdf',scifort_geninvgauss_pdf_f64(1.3_dp,1.2_dp,2.3_dp,0.0_dp,1.0_dp), &
            geninvgauss_pdf(1.3_dp,1.2_dp,2.3_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI GIG cdf',scifort_geninvgauss_cdf_f64(1.3_dp,1.2_dp,2.3_dp,0.0_dp,1.0_dp), &
            geninvgauss_cdf(1.3_dp,1.2_dp,2.3_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI GIG ppf',scifort_geninvgauss_ppf_f64(0.7_dp,1.2_dp,2.3_dp,0.0_dp,1.0_dp), &
            geninvgauss_ppf(0.7_dp,1.2_dp,2.3_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI NIG pdf',scifort_norminvgauss_pdf_f64(0.3_dp,2.5_dp,0.7_dp,0.0_dp,1.0_dp), &
            norminvgauss_pdf(0.3_dp,2.5_dp,0.7_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI NIG cdf',scifort_norminvgauss_cdf_f64(0.3_dp,2.5_dp,0.7_dp,0.0_dp,1.0_dp), &
            norminvgauss_cdf(0.3_dp,2.5_dp,0.7_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI NIG ppf',scifort_norminvgauss_ppf_f64(0.7_dp,2.5_dp,0.7_dp,0.0_dp,1.0_dp), &
            norminvgauss_ppf(0.7_dp,2.5_dp,0.7_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI Skellam pmf',scifort_skellam_pmf_f64(2.0_dp,3.2_dp,1.7_dp,0.0_dp), &
            skellam_pmf(2.0_dp,3.2_dp,1.7_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI Skellam cdf',scifort_skellam_cdf_f64(2.0_dp,3.2_dp,1.7_dp,0.0_dp), &
            skellam_cdf(2.0_dp,3.2_dp,1.7_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI Skellam ppf',scifort_skellam_ppf_f64(0.7_dp,3.2_dp,1.7_dp,0.0_dp), &
            skellam_ppf(0.7_dp,3.2_dp,1.7_dp),0.0_dp,0.0_dp,failures)
    end subroutine test_rng_fit_c_api

    pure function centered(plus_value,minus_value,h) result(value)
        real(dp), intent(in) :: plus_value !! function value at theta+h
        real(dp), intent(in) :: minus_value !! function value at theta-h
        real(dp), intent(in) :: h !! positive centered-difference step
        real(dp) :: value
        value=(plus_value-minus_value)/(2.0_dp*h)
    end function centered

end program test_gig_nig_skellam
