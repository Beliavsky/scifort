! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Compact full-stack regression coverage for hyperbolic-secant and half-logistic additions.

program test_hypsecant_halflogistic
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_c_api
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state, rng_uniform
    use scifort_stats
    use test_checks, only : check_close, check_true
    implicit none

    integer :: failures

    failures = 0
    call test_reference_values(failures)
    call test_identities_and_tails(failures)
    call test_likelihood_scores(failures)
    call test_random_fit_c_api(failures)

    if (failures == 0) then
        print '(a)', 'test_hypsecant_halflogistic: PASS'
    else
        print '(a,1x,i0)', 'test_hypsecant_halflogistic: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference_values(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: xh(5) = [-3.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 3.0_dp]
        real(dp), parameter :: hypsecant_pdf_ref(5) = [ &
            0.03161706127175161_dp, 0.2062820820908705_dp, &
            0.3183098861837907_dp, 0.2062820820908705_dp, &
            0.03161706127175161_dp]
        real(dp), parameter :: hypsecant_cdf_ref(5) = [ &
            0.03166928263726926_dp, 0.22441701432858502_dp, 0.5_dp, &
            0.775582985671415_dp, 0.9683307173627308_dp]
        real(dp), parameter :: xl(5) = [0.0_dp, 0.2_dp, 1.0_dp, 2.0_dp, 5.0_dp]
        real(dp), parameter :: halflogistic_pdf_ref(5) = [ &
            0.5_dp, 0.4950331454237199_dp, 0.39322386648296365_dp, &
            0.20998717080701307_dp, 0.01329611334158031_dp]
        real(dp), parameter :: halflogistic_cdf_ref(5) = [ &
            0.0_dp, 0.09966799462495582_dp, 0.46211715726000974_dp, &
            0.7615941559557649_dp, 0.9866142981514303_dp]
        integer :: i
        character(len=56) :: label

        do i = 1, size(xh)
            write(label, '(a,i0)') 'hypsecant pdf ref ', i
            call check_close(trim(label), hypsecant_pdf(xh(i)), hypsecant_pdf_ref(i), &
                8.0e-15_dp, 8.0e-14_dp, failures)
            write(label, '(a,i0)') 'hypsecant cdf ref ', i
            call check_close(trim(label), hypsecant_cdf(xh(i)), hypsecant_cdf_ref(i), &
                8.0e-15_dp, 8.0e-14_dp, failures)
            write(label, '(a,i0)') 'halflogistic pdf ref ', i
            call check_close(trim(label), halflogistic_pdf(xl(i)), &
                halflogistic_pdf_ref(i), 8.0e-15_dp, 8.0e-14_dp, failures)
            write(label, '(a,i0)') 'halflogistic cdf ref ', i
            call check_close(trim(label), halflogistic_cdf(xl(i)), &
                halflogistic_cdf_ref(i), 8.0e-15_dp, 8.0e-14_dp, failures)
        end do

        call check_close('hypsecant ppf 0.1', hypsecant_ppf(0.1_dp), &
            -1.842730034701113_dp, 2.0e-14_dp, 2.0e-13_dp, failures)
        call check_close('hypsecant ppf 0.9', hypsecant_ppf(0.9_dp), &
            1.8427300347011126_dp, 2.0e-14_dp, 2.0e-13_dp, failures)
        call check_close('halflogistic ppf 0.1', halflogistic_ppf(0.1_dp), &
            0.20067069546215116_dp, 2.0e-14_dp, 2.0e-13_dp, failures)
        call check_close('halflogistic ppf 0.9', halflogistic_ppf(0.9_dp), &
            2.9444389791664407_dp, 2.0e-14_dp, 2.0e-13_dp, failures)
    end subroutine test_reference_values

    subroutine test_identities_and_tails(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: probs(5) = [ &
            1.0e-10_dp, 0.1_dp, 0.5_dp, 0.9_dp, 1.0_dp - 1.0e-10_dp]
        real(dp), parameter :: points(5) = [-4.0_dp, -0.5_dp, 0.0_dp, 1.0_dp, 8.0_dp]
        integer :: i
        real(dp) :: x

        do i = 1, size(points)
            call check_close('hypsecant cdf+sf', &
                hypsecant_cdf(points(i), 0.2_dp, 1.3_dp) + &
                hypsecant_sf(points(i), 0.2_dp, 1.3_dp), &
                1.0_dp, 4.0e-15_dp, 4.0e-15_dp, failures)
        end do

        do i = 1, size(probs)
            x = hypsecant_ppf(probs(i), -0.2_dp, 1.7_dp)
            call check_close('hypsecant ppf/cdf', hypsecant_cdf(x, -0.2_dp, 1.7_dp), &
                probs(i), 3.0e-13_dp, 3.0e-8_dp, failures)
            x = halflogistic_isf(probs(i), -0.2_dp, 1.7_dp)
            call check_close('halflogistic isf/sf', &
                halflogistic_sf(x, -0.2_dp, 1.7_dp), &
                probs(i), 3.0e-13_dp, 3.0e-8_dp, failures)
        end do

        call check_close('hypsecant deep logcdf', hypsecant_logcdf(-1000.0_dp), &
            log(2.0_dp / acos(-1.0_dp)) - 1000.0_dp, &
            2.0e-13_dp, 2.0e-13_dp, failures)
        call check_close('halflogistic deep logsf', halflogistic_logsf(1000.0_dp), &
            log(2.0_dp) - 1000.0_dp, 2.0e-13_dp, 2.0e-13_dp, failures)
        call check_true('hypsecant invalid scale nan', &
            ieee_is_nan(hypsecant_pdf(0.0_dp, scale=0.0_dp)), failures)
        call check_true('halflogistic invalid scale nan', &
            ieee_is_nan(halflogistic_cdf(1.0_dp, scale=0.0_dp)), failures)
    end subroutine test_identities_and_tails

    subroutine test_likelihood_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: xh(5) = [-1.7_dp, -0.4_dp, 0.2_dp, 1.1_dp, 2.3_dp]
        real(dp), parameter :: xl(5) = [0.3_dp, 0.6_dp, 1.0_dp, 1.8_dp, 3.0_dp]
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp), parameter :: tol = 1.0e-5_dp
        real(dp) :: sh(2)
        real(dp) :: sl(2)

        call check_close('hypsecant loglike sum', &
            hypsecant_loglikelihood(xh, 0.1_dp, 1.4_dp), &
            sum(hypsecant_logpdf(xh, 0.1_dp, 1.4_dp)), &
            3.0e-15_dp, 3.0e-15_dp, failures)
        call check_close('halflogistic loglike sum', &
            halflogistic_loglikelihood(xl, 0.1_dp, 1.4_dp), &
            sum(halflogistic_logpdf(xl, 0.1_dp, 1.4_dp)), &
            3.0e-15_dp, 3.0e-15_dp, failures)

        sh = hypsecant_score(xh, 0.1_dp, 1.4_dp)
        call check_close('hypsecant score loc', sh(1), centered( &
            hypsecant_loglikelihood(xh, 0.1_dp + h, 1.4_dp), &
            hypsecant_loglikelihood(xh, 0.1_dp - h, 1.4_dp), h), &
            tol, tol, failures)
        call check_close('hypsecant score scale', sh(2), centered( &
            hypsecant_loglikelihood(xh, 0.1_dp, 1.4_dp + h), &
            hypsecant_loglikelihood(xh, 0.1_dp, 1.4_dp - h), h), &
            tol, tol, failures)

        sl = halflogistic_score(xl, 0.1_dp, 1.4_dp)
        call check_close('halflogistic score loc', sl(1), centered( &
            halflogistic_loglikelihood(xl, 0.1_dp + h, 1.4_dp), &
            halflogistic_loglikelihood(xl, 0.1_dp - h, 1.4_dp), h), &
            tol, tol, failures)
        call check_close('halflogistic score scale', sl(2), centered( &
            halflogistic_loglikelihood(xl, 0.1_dp, 1.4_dp + h), &
            halflogistic_loglikelihood(xl, 0.1_dp, 1.4_dp - h), h), &
            tol, tol, failures)
    end subroutine test_likelihood_scores

    subroutine test_random_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: xh(5) = [-1.7_dp, -0.4_dp, 0.2_dp, 1.1_dp, 2.3_dp]
        real(dp), parameter :: xl(5) = [0.3_dp, 0.6_dp, 1.0_dp, 1.8_dp, 3.0_dp]
        real(dp) :: actual
        real(dp) :: expected
        type(fit_result) :: result
        type(rng_state) :: reference
        type(rng_state) :: state

        call rng_seed(state, 97531)
        call rng_seed(reference, 97531)
        actual = hypsecant_rvs(state, -0.2_dp, 1.8_dp)
        expected = hypsecant_ppf(rng_uniform(reference), -0.2_dp, 1.8_dp)
        call check_close('rvs hypsecant', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = halflogistic_rvs(state, 0.1_dp, 1.3_dp)
        expected = halflogistic_ppf(rng_uniform(reference), 0.1_dp, 1.3_dp)
        call check_close('rvs halflogistic', actual, expected, 0.0_dp, 0.0_dp, failures)

        call hypsecant_fit(xh, [0.1_dp, 1.4_dp], [0.1_dp, 1.4_dp], result)
        call check_true('fit hypsecant fixed success', result%success, failures)
        call halflogistic_fit(xl, [0.1_dp, 1.4_dp], [0.1_dp, 1.4_dp], result)
        call check_true('fit halflogistic fixed success', result%success, failures)

        call check_close('C ABI hypsecant cdf', &
            scifort_hypsecant_cdf_f64(0.2_dp, -0.2_dp, 1.8_dp), &
            hypsecant_cdf(0.2_dp, -0.2_dp, 1.8_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI halflogistic ppf', &
            scifort_halflogistic_ppf_f64(0.4_dp, -0.2_dp, 1.8_dp), &
            halflogistic_ppf(0.4_dp, -0.2_dp, 1.8_dp), 0.0_dp, 0.0_dp, failures)
    end subroutine test_random_fit_c_api

    pure function centered(plus, minus, h) result(value)
        real(dp), intent(in) :: plus !! objective at theta + h
        real(dp), intent(in) :: minus !! objective at theta - h
        real(dp), intent(in) :: h !! positive centered-difference half step
        real(dp) :: value

        value = (plus - minus) / (2.0_dp * h)
    end function centered

end program test_hypsecant_halflogistic
