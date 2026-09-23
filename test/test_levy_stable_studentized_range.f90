! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_levy_stable_studentized_range
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_c_api
    use scifort_fit, only : fit_result
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state, rng_uniform
    use scifort_stats
    use test_checks, only : check_close, check_true
    use test_levy_stable_studentized_range_reference
    implicit none

    integer :: failures

    failures = 0
    call test_reference(failures)
    call test_structure(failures)
    call test_scores(failures)
    call test_rng_fit_c_api(failures)

    if (failures == 0) then
        print '(a)', 'test_levy_stable_studentized_range: PASS'
    else
        print '(a,1x,i0)', 'test_levy_stable_studentized_range: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        integer :: i

        do i = 1, size(ls_x_ref)
            call check_close('levy_stable pdf', levy_stable_pdf(ls_x_ref(i),ls_alpha_ref(i),ls_beta_ref(i)), &
                ls_pdf_ref(i), 5.0e-13_dp, 5.0e-11_dp, failures)
            call check_close('levy_stable logpdf', levy_stable_logpdf(ls_x_ref(i),ls_alpha_ref(i),ls_beta_ref(i)), &
                ls_logpdf_ref(i), 2.0e-12_dp, 5.0e-11_dp, failures)
            call check_close('levy_stable cdf', levy_stable_cdf(ls_x_ref(i),ls_alpha_ref(i),ls_beta_ref(i)), &
                ls_cdf_ref(i), 5.0e-13_dp, 5.0e-11_dp, failures)
            call check_close('levy_stable sf', levy_stable_sf(ls_x_ref(i),ls_alpha_ref(i),ls_beta_ref(i)), &
                ls_sf_ref(i), 5.0e-13_dp, 5.0e-11_dp, failures)
            call check_close('levy_stable logcdf', levy_stable_logcdf(ls_x_ref(i),ls_alpha_ref(i),ls_beta_ref(i)), &
                ls_logcdf_ref(i), 3.0e-12_dp, 8.0e-11_dp, failures)
            call check_close('levy_stable logsf', levy_stable_logsf(ls_x_ref(i),ls_alpha_ref(i),ls_beta_ref(i)), &
                ls_logsf_ref(i), 3.0e-12_dp, 8.0e-11_dp, failures)
        end do
        do i = 1, size(ls_prob_ref)
            call check_close('levy_stable ppf', levy_stable_ppf(ls_prob_ref(i),ls_qalpha_ref(i),ls_qbeta_ref(i)), &
                ls_ppf_ref(i), 3.0e-10_dp, 3.0e-10_dp, failures)
            call check_close('levy_stable isf', levy_stable_isf(ls_prob_ref(i),ls_qalpha_ref(i),ls_qbeta_ref(i)), &
                ls_isf_ref(i), 3.0e-10_dp, 3.0e-10_dp, failures)
        end do

        do i = 1, size(sr_x_ref)
            call check_close('studentized_range pdf', studentized_range_pdf(sr_x_ref(i),sr_k_ref(i),sr_df_ref(i)), &
                sr_pdf_ref(i), 3.0e-11_dp, 5.0e-7_dp, failures)
            call check_close('studentized_range logpdf', studentized_range_logpdf(sr_x_ref(i),sr_k_ref(i),sr_df_ref(i)), &
                sr_logpdf_ref(i), 3.0e-10_dp, 5.0e-7_dp, failures)
            call check_close('studentized_range cdf', studentized_range_cdf(sr_x_ref(i),sr_k_ref(i),sr_df_ref(i)), &
                sr_cdf_ref(i), 5.0e-11_dp, 5.0e-10_dp, failures)
            call check_close('studentized_range sf', studentized_range_sf(sr_x_ref(i),sr_k_ref(i),sr_df_ref(i)), &
                sr_sf_ref(i), 5.0e-11_dp, 5.0e-10_dp, failures)
            call check_close('studentized_range logcdf', studentized_range_logcdf(sr_x_ref(i),sr_k_ref(i),sr_df_ref(i)), &
                sr_logcdf_ref(i), 5.0e-10_dp, 5.0e-9_dp, failures)
            call check_close('studentized_range logsf', studentized_range_logsf(sr_x_ref(i),sr_k_ref(i),sr_df_ref(i)), &
                sr_logsf_ref(i), 5.0e-9_dp, 5.0e-8_dp, failures)
        end do
        do i = 1, size(sr_prob_ref)
            call check_close('studentized_range ppf', studentized_range_ppf(sr_prob_ref(i),sr_qk_ref(i),sr_qdf_ref(i)), &
                sr_ppf_ref(i), 5.0e-10_dp, 5.0e-10_dp, failures)
            call check_close('studentized_range isf', studentized_range_isf(sr_prob_ref(i),sr_qk_ref(i),sr_qdf_ref(i)), &
                sr_isf_ref(i), 5.0e-10_dp, 5.0e-10_dp, failures)
        end do
    end subroutine test_reference

    subroutine test_structure(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp) :: shift, x, p

        call check_close('levy_stable alpha2 normal pdf', levy_stable_pdf(0.7_dp,2.0_dp,0.9_dp), &
            normal_pdf(0.7_dp/sqrt(2.0_dp))/sqrt(2.0_dp), 3.0e-15_dp, 3.0e-15_dp, failures)
        call check_close('levy_stable cauchy pdf', levy_stable_pdf(0.7_dp,1.0_dp,0.0_dp), &
            cauchy_pdf(0.7_dp), 3.0e-15_dp, 3.0e-15_dp, failures)
        call check_close('levy_stable levy pdf', levy_stable_pdf(1.7_dp,0.5_dp,1.0_dp), &
            levy_pdf(1.7_dp), 3.0e-15_dp, 3.0e-15_dp, failures)
        call check_close('levy_stable symmetry', levy_stable_sf(1.3_dp,1.4_dp,0.6_dp), &
            levy_stable_cdf(-1.3_dp,1.4_dp,-0.6_dp), 3.0e-15_dp, 3.0e-15_dp, failures)

        shift = 2.0_dp * 0.4_dp * 1.7_dp * log(1.7_dp) / acos(-1.0_dp)
        x = 0.3_dp + shift + 1.7_dp * 0.8_dp
        call check_close('levy_stable alpha1 S1 scale shift pdf', &
            levy_stable_pdf(x,1.0_dp,0.4_dp,0.3_dp,1.7_dp), &
            levy_stable_pdf(0.8_dp,1.0_dp,0.4_dp)/1.7_dp, 5.0e-13_dp, 5.0e-12_dp, failures)
        call check_close('levy_stable alpha1 S1 scale shift cdf', &
            levy_stable_cdf(x,1.0_dp,0.4_dp,0.3_dp,1.7_dp), &
            levy_stable_cdf(0.8_dp,1.0_dp,0.4_dp), 5.0e-13_dp, 5.0e-12_dp, failures)
        call check_true('levy_stable invalid alpha', ieee_is_nan(levy_stable_pdf(0.0_dp,0.0_dp,0.0_dp)), failures)
        call check_true('levy_stable invalid beta', ieee_is_nan(levy_stable_cdf(0.0_dp,1.2_dp,1.1_dp)), failures)

        call check_close('studentized_range below support cdf',studentized_range_cdf(-0.1_dp,5.0_dp,10.0_dp), &
            0.0_dp,0.0_dp,0.0_dp,failures)
        call check_close('studentized_range below support pdf',studentized_range_pdf(-0.1_dp,5.0_dp,10.0_dp), &
            0.0_dp,0.0_dp,0.0_dp,failures)
        call check_close('studentized_range k2 dfInf pdf0',studentized_range_pdf(0.0_dp,2.0_dp,1.0e5_dp), &
            1.0_dp/sqrt(acos(-1.0_dp)),3.0e-15_dp,3.0e-15_dp,failures)
        call check_close('studentized_range shift scale cdf',studentized_range_cdf(5.2_dp,5.0_dp,10.0_dp,0.4_dp,1.6_dp), &
            studentized_range_cdf(3.0_dp,5.0_dp,10.0_dp),5.0e-12_dp,5.0e-11_dp,failures)
        call check_close('studentized_range shift scale pdf',studentized_range_pdf(5.2_dp,5.0_dp,10.0_dp,0.4_dp,1.6_dp), &
            studentized_range_pdf(3.0_dp,5.0_dp,10.0_dp)/1.6_dp,5.0e-12_dp,5.0e-11_dp,failures)
        call check_true('studentized_range invalid k',ieee_is_nan(studentized_range_pdf(1.0_dp,1.0_dp,10.0_dp)),failures)
        call check_true('studentized_range invalid df',ieee_is_nan(studentized_range_cdf(1.0_dp,3.0_dp,0.0_dp)),failures)

        p = 1.0e-4_dp
        x = levy_stable_ppf(p,1.4_dp,0.35_dp)
        call check_close('levy_stable cdf-ppf',levy_stable_cdf(x,1.4_dp,0.35_dp),p,5.0e-11_dp,5.0e-9_dp,failures)
        x = studentized_range_isf(p,5.0_dp,10.0_dp)
        call check_close('studentized_range sf-isf',studentized_range_sf(x,5.0_dp,10.0_dp),p,5.0e-11_dp,5.0e-8_dp,failures)
    end subroutine test_structure

    subroutine test_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: ld(4) = [-1.2_dp,-0.1_dp,0.8_dp,2.0_dp]
        real(dp), parameter :: sd(4) = [1.1_dp,2.0_dp,3.2_dp,4.5_dp]
        real(dp), parameter :: h = 2.0e-5_dp
        real(dp) :: score_l(4), score_s(4), fd

        score_l = levy_stable_score(ld,1.4_dp,0.25_dp,0.1_dp,1.2_dp)
        fd = centered(levy_stable_loglikelihood(ld,1.4_dp+h,0.25_dp,0.1_dp,1.2_dp), &
            levy_stable_loglikelihood(ld,1.4_dp-h,0.25_dp,0.1_dp,1.2_dp),h)
        call check_close('levy_stable score alpha',score_l(1),fd,2.0e-5_dp,2.0e-5_dp,failures)
        fd = centered(levy_stable_loglikelihood(ld,1.4_dp,0.25_dp+h,0.1_dp,1.2_dp), &
            levy_stable_loglikelihood(ld,1.4_dp,0.25_dp-h,0.1_dp,1.2_dp),h)
        call check_close('levy_stable score beta',score_l(2),fd,2.0e-5_dp,2.0e-5_dp,failures)

        score_s = studentized_range_score(sd,5.0_dp,10.0_dp,0.2_dp,1.1_dp)
        fd = centered(studentized_range_loglikelihood(sd,5.0_dp+h,10.0_dp,0.2_dp,1.1_dp), &
            studentized_range_loglikelihood(sd,5.0_dp-h,10.0_dp,0.2_dp,1.1_dp),h)
        call check_close('studentized_range score k',score_s(1),fd,2.0e-5_dp,2.0e-5_dp,failures)
        fd = centered(studentized_range_loglikelihood(sd,5.0_dp,10.0_dp+h,0.2_dp,1.1_dp), &
            studentized_range_loglikelihood(sd,5.0_dp,10.0_dp-h,0.2_dp,1.1_dp),h)
        call check_close('studentized_range score df',score_s(2),fd,2.0e-5_dp,2.0e-5_dp,failures)
    end subroutine test_scores

    subroutine test_rng_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: ld(4) = [-1.2_dp,-0.1_dp,0.8_dp,2.0_dp]
        real(dp), parameter :: sd(4) = [1.1_dp,2.0_dp,3.2_dp,4.5_dp]
        type(rng_state) :: state, reference
        type(fit_result) :: result
        real(dp) :: actual, expected

        call rng_seed(state,24680)
        call rng_seed(reference,24680)
        actual = levy_stable_rvs(state,1.4_dp,0.25_dp,0.1_dp,1.2_dp)
        expected = levy_stable_ppf(rng_uniform(reference),1.4_dp,0.25_dp,0.1_dp,1.2_dp)
        call check_close('levy_stable rvs',actual,expected,0.0_dp,0.0_dp,failures)

        call rng_seed(state,24680)
        call rng_seed(reference,24680)
        actual = studentized_range_rvs(state,5.0_dp,10.0_dp,0.2_dp,1.1_dp)
        expected = studentized_range_ppf(rng_uniform(reference),5.0_dp,10.0_dp,0.2_dp,1.1_dp)
        call check_close('studentized_range rvs',actual,expected,0.0_dp,0.0_dp,failures)

        call levy_stable_fit(ld,[1.4_dp,0.25_dp,0.1_dp,1.2_dp],[1.4_dp,0.25_dp,0.1_dp,1.2_dp],result)
        call check_true('levy_stable fixed fit',result%success,failures)
        call studentized_range_fit(sd,[5.0_dp,10.0_dp,0.2_dp,1.1_dp],[5.0_dp,10.0_dp,0.2_dp,1.1_dp],result)
        call check_true('studentized_range fixed fit',result%success,failures)

        call check_close('C ABI levy_stable pdf',scifort_levy_stable_pdf_f64(0.7_dp,1.4_dp,0.25_dp,0.0_dp,1.0_dp), &
            levy_stable_pdf(0.7_dp,1.4_dp,0.25_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI levy_stable cdf',scifort_levy_stable_cdf_f64(0.7_dp,1.4_dp,0.25_dp,0.0_dp,1.0_dp), &
            levy_stable_cdf(0.7_dp,1.4_dp,0.25_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI levy_stable ppf',scifort_levy_stable_ppf_f64(0.7_dp,1.4_dp,0.25_dp,0.0_dp,1.0_dp), &
            levy_stable_ppf(0.7_dp,1.4_dp,0.25_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI studentized_range pdf',scifort_studentized_range_pdf_f64(3.0_dp,5.0_dp,10.0_dp,0.0_dp,1.0_dp), &
            studentized_range_pdf(3.0_dp,5.0_dp,10.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI studentized_range cdf',scifort_studentized_range_cdf_f64(3.0_dp,5.0_dp,10.0_dp,0.0_dp,1.0_dp), &
            studentized_range_cdf(3.0_dp,5.0_dp,10.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI studentized_range ppf',scifort_studentized_range_ppf_f64(0.7_dp,5.0_dp,10.0_dp,0.0_dp,1.0_dp), &
            studentized_range_ppf(0.7_dp,5.0_dp,10.0_dp),0.0_dp,0.0_dp,failures)
    end subroutine test_rng_fit_c_api

    pure function centered(plus_value, minus_value, h) result(value)
        real(dp), intent(in) :: plus_value !! function value at theta+h
        real(dp), intent(in) :: minus_value !! function value at theta-h
        real(dp), intent(in) :: h !! positive centered-difference step
        real(dp) :: value
        value = (plus_value - minus_value) / (2.0_dp*h)
    end function centered

end program test_levy_stable_studentized_range
