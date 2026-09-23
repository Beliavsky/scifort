! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Full-stack regression coverage for inverse gamma, inverse Gaussian, Levy,
! and log-Laplace distributions.

program test_inverse_positive
    use scifort_c_api
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state, rng_uniform
    use scifort_stats
    use test_checks, only : check_close, check_true
    implicit none

    integer :: failures

    failures = 0
    call test_reference_values(failures)
    call test_tail_identities(failures)
    call test_likelihood_scores(failures)
    call test_random_fit_c_api(failures)

    if (failures == 0) then
        print '(a)', 'test_inverse_positive: PASS'
    else
        print '(a,1x,i0)', 'test_inverse_positive: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference_values(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: x(5) = [0.1_dp, 0.5_dp, 1.0_dp, 2.0_dp, 5.0_dp]
        real(dp), parameter :: invgamma_pdf_ref(5) = [ &
            0.10799881274785479_dp, 1.1518072856146782_dp, &
            0.2767383316137298_dp, 0.0403284540865239_dp, &
            0.002203483935860386_dp]
        real(dp), parameter :: invgamma_cdf_ref(5) = [ &
            0.0012497305630313773_dp, 0.5494159513527802_dp, &
            0.8491450360846096_dp, 0.9625657732472964_dp, &
            0.9953295932358704_dp]
        real(dp), parameter :: invgauss_pdf_ref(5) = [ &
            0.1504496747647393_dp, 0.6855816625980469_dp, &
            0.3665158103980622_dp, 0.13995353963889726_dp, &
            0.02447967712281798_dp]
        real(dp), parameter :: invgauss_cdf_ref(5) = [ &
            0.002777346394023719_dp, 0.26807883088722934_dp, &
            0.522236848726317_dp, 0.750396685145634_dp, &
            0.933771635802996_dp]
        real(dp), parameter :: levy_pdf_ref(5) = [ &
            0.08500366602520341_dp, 0.4151074974205947_dp, &
            0.24197072451914337_dp, 0.1098478223669306_dp, &
            0.03228684517430723_dp]
        real(dp), parameter :: levy_cdf_ref(5) = [ &
            0.001565402258002548_dp, 0.15729920705028516_dp, &
            0.31731050786291415_dp, 0.4795001221869535_dp, &
            0.6547208460185769_dp]
        real(dp), parameter :: loglaplace_pdf_ref(5) = [ &
            0.14264038732150022_dp, 0.5169142597486657_dp, 0.9_dp, &
            0.12922856493716645_dp, 0.009934053562520749_dp]
        real(dp), parameter :: loglaplace_cdf_ref(5) = [ &
            0.007924465962305567_dp, 0.14358729437462936_dp, 0.5_dp, &
            0.8564127056253706_dp, 0.9724054067707757_dp]
        integer :: i
        character(len=56) :: label

        do i = 1, size(x)
            write(label, '(a,i0)') 'invgamma pdf ref ', i
            call check_close(trim(label), invgamma_pdf(x(i), 2.5_dp), &
                invgamma_pdf_ref(i), 2.0e-14_dp, 2.0e-13_dp, failures)
            write(label, '(a,i0)') 'invgamma cdf ref ', i
            call check_close(trim(label), invgamma_cdf(x(i), 2.5_dp), &
                invgamma_cdf_ref(i), 2.0e-13_dp, 2.0e-12_dp, failures)
            write(label, '(a,i0)') 'invgauss pdf ref ', i
            call check_close(trim(label), invgauss_pdf(x(i), 1.7_dp), &
                invgauss_pdf_ref(i), 2.0e-14_dp, 2.0e-13_dp, failures)
            write(label, '(a,i0)') 'invgauss cdf ref ', i
            call check_close(trim(label), invgauss_cdf(x(i), 1.7_dp), &
                invgauss_cdf_ref(i), 3.0e-14_dp, 3.0e-13_dp, failures)
            write(label, '(a,i0)') 'levy pdf ref ', i
            call check_close(trim(label), levy_pdf(x(i)), levy_pdf_ref(i), &
                2.0e-14_dp, 2.0e-13_dp, failures)
            write(label, '(a,i0)') 'levy cdf ref ', i
            call check_close(trim(label), levy_cdf(x(i)), levy_cdf_ref(i), &
                2.0e-14_dp, 2.0e-13_dp, failures)
            write(label, '(a,i0)') 'loglaplace pdf ref ', i
            call check_close(trim(label), loglaplace_pdf(x(i), 1.8_dp), &
                loglaplace_pdf_ref(i), 2.0e-14_dp, 2.0e-13_dp, failures)
            write(label, '(a,i0)') 'loglaplace cdf ref ', i
            call check_close(trim(label), loglaplace_cdf(x(i), 1.8_dp), &
                loglaplace_cdf_ref(i), 2.0e-14_dp, 2.0e-13_dp, failures)
        end do

        call check_close('invgamma ppf 0.1', invgamma_ppf(0.1_dp, 2.5_dp), &
            0.21653559100205366_dp, 3.0e-13_dp, 3.0e-12_dp, failures)
        call check_close('invgamma ppf 0.9', invgamma_ppf(0.9_dp, 2.5_dp), &
            1.2419984352016977_dp, 3.0e-13_dp, 3.0e-12_dp, failures)
        call check_close('invgauss ppf 0.1', invgauss_ppf(0.1_dp, 1.7_dp), &
            0.2771802643381196_dp, 4.0e-13_dp, 4.0e-12_dp, failures)
        call check_close('invgauss ppf 0.9', invgauss_ppf(0.9_dp, 1.7_dp), &
            3.938673759695976_dp, 4.0e-12_dp, 4.0e-12_dp, failures)
        call check_close('levy ppf 0.1', levy_ppf(0.1_dp), &
            0.3696115094681948_dp, 3.0e-13_dp, 3.0e-12_dp, failures)
        call check_close('levy ppf 0.9', levy_ppf(0.9_dp), &
            63.328117677016756_dp, 3.0e-11_dp, 3.0e-12_dp, failures)
        call check_close('loglaplace ppf 0.1', loglaplace_ppf(0.1_dp, 1.8_dp), &
            0.4089623530229582_dp, 2.0e-14_dp, 2.0e-13_dp, failures)
        call check_close('loglaplace ppf 0.9', loglaplace_ppf(0.9_dp, 1.8_dp), &
            2.445212848097689_dp, 2.0e-13_dp, 2.0e-13_dp, failures)
    end subroutine test_reference_values

    subroutine test_tail_identities(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: probs(5) = [ &
            1.0e-8_dp, 0.1_dp, 0.5_dp, 0.9_dp, 1.0_dp - 1.0e-8_dp]
        integer :: i
        real(dp) :: x

        do i = 1, size(probs)
            x = invgamma_ppf(probs(i), 2.2_dp, -0.3_dp, 1.4_dp)
            call check_close('invgamma ppf/cdf', invgamma_cdf(x, 2.2_dp, -0.3_dp, 1.4_dp), &
                probs(i), 3.0e-11_dp, 3.0e-8_dp, failures)
            x = invgauss_ppf(probs(i), 1.3_dp, -0.3_dp, 1.4_dp)
            call check_close('invgauss ppf/cdf', invgauss_cdf(x, 1.3_dp, -0.3_dp, 1.4_dp), &
                probs(i), 3.0e-11_dp, 3.0e-8_dp, failures)
            x = levy_isf(probs(i), -0.3_dp, 1.4_dp)
            call check_close('levy isf/sf', levy_sf(x, -0.3_dp, 1.4_dp), &
                probs(i), 3.0e-11_dp, 3.0e-8_dp, failures)
            x = loglaplace_isf(probs(i), 1.7_dp, -0.3_dp, 1.4_dp)
            call check_close('loglaplace isf/sf', loglaplace_sf(x, 1.7_dp, -0.3_dp, 1.4_dp), &
                probs(i), 3.0e-11_dp, 3.0e-8_dp, failures)
        end do
    end subroutine test_tail_identities

    subroutine test_likelihood_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: xi(5) = [0.4_dp, 0.7_dp, 1.1_dp, 1.8_dp, 3.0_dp]
        real(dp), parameter :: xl(5) = [0.4_dp, 0.8_dp, 1.6_dp, 2.7_dp, 4.0_dp]
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp), parameter :: tol = 3.0e-5_dp
        real(dp) :: s2(2)
        real(dp) :: s3(3)

        s3 = invgamma_score(xi, 2.3_dp, 0.1_dp, 1.2_dp)
        call check_close('invgamma score a', s3(1), centered( &
            invgamma_loglikelihood(xi, 2.3_dp + h, 0.1_dp, 1.2_dp), &
            invgamma_loglikelihood(xi, 2.3_dp - h, 0.1_dp, 1.2_dp), h), tol, tol, failures)
        call check_close('invgamma score loc', s3(2), centered( &
            invgamma_loglikelihood(xi, 2.3_dp, 0.1_dp + h, 1.2_dp), &
            invgamma_loglikelihood(xi, 2.3_dp, 0.1_dp - h, 1.2_dp), h), tol, tol, failures)
        call check_close('invgamma score scale', s3(3), centered( &
            invgamma_loglikelihood(xi, 2.3_dp, 0.1_dp, 1.2_dp + h), &
            invgamma_loglikelihood(xi, 2.3_dp, 0.1_dp, 1.2_dp - h), h), tol, tol, failures)

        s3 = invgauss_score(xi, 1.4_dp, 0.1_dp, 1.2_dp)
        call check_close('invgauss score shape', s3(1), centered( &
            invgauss_loglikelihood(xi, 1.4_dp + h, 0.1_dp, 1.2_dp), &
            invgauss_loglikelihood(xi, 1.4_dp - h, 0.1_dp, 1.2_dp), h), tol, tol, failures)
        call check_close('invgauss score loc', s3(2), centered( &
            invgauss_loglikelihood(xi, 1.4_dp, 0.1_dp + h, 1.2_dp), &
            invgauss_loglikelihood(xi, 1.4_dp, 0.1_dp - h, 1.2_dp), h), tol, tol, failures)
        call check_close('invgauss score scale', s3(3), centered( &
            invgauss_loglikelihood(xi, 1.4_dp, 0.1_dp, 1.2_dp + h), &
            invgauss_loglikelihood(xi, 1.4_dp, 0.1_dp, 1.2_dp - h), h), tol, tol, failures)

        s2 = levy_score(xi, 0.1_dp, 1.2_dp)
        call check_close('levy score loc', s2(1), centered( &
            levy_loglikelihood(xi, 0.1_dp + h, 1.2_dp), &
            levy_loglikelihood(xi, 0.1_dp - h, 1.2_dp), h), tol, tol, failures)
        call check_close('levy score scale', s2(2), centered( &
            levy_loglikelihood(xi, 0.1_dp, 1.2_dp + h), &
            levy_loglikelihood(xi, 0.1_dp, 1.2_dp - h), h), tol, tol, failures)

        s3 = loglaplace_score(xl, 1.8_dp, 0.1_dp, 1.2_dp)
        call check_close('loglaplace score c', s3(1), centered( &
            loglaplace_loglikelihood(xl, 1.8_dp + h, 0.1_dp, 1.2_dp), &
            loglaplace_loglikelihood(xl, 1.8_dp - h, 0.1_dp, 1.2_dp), h), tol, tol, failures)
        call check_close('loglaplace score loc', s3(2), centered( &
            loglaplace_loglikelihood(xl, 1.8_dp, 0.1_dp + h, 1.2_dp), &
            loglaplace_loglikelihood(xl, 1.8_dp, 0.1_dp - h, 1.2_dp), h), tol, tol, failures)
        call check_close('loglaplace score scale', s3(3), centered( &
            loglaplace_loglikelihood(xl, 1.8_dp, 0.1_dp, 1.2_dp + h), &
            loglaplace_loglikelihood(xl, 1.8_dp, 0.1_dp, 1.2_dp - h), h), tol, tol, failures)
    end subroutine test_likelihood_scores

    subroutine test_random_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: xi(5) = [0.4_dp, 0.7_dp, 1.1_dp, 1.8_dp, 3.0_dp]
        real(dp), parameter :: xl(5) = [0.4_dp, 0.8_dp, 1.6_dp, 2.7_dp, 4.0_dp]
        real(dp) :: actual
        real(dp) :: expected
        type(fit_result) :: result
        type(rng_state) :: reference
        type(rng_state) :: state

        call rng_seed(state, 86420)
        call rng_seed(reference, 86420)
        actual = invgamma_rvs(state, 2.2_dp, -0.2_dp, 1.3_dp)
        expected = invgamma_ppf(rng_uniform(reference), 2.2_dp, -0.2_dp, 1.3_dp)
        call check_close('rvs invgamma', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = invgauss_rvs(state, 1.2_dp, -0.2_dp, 1.3_dp)
        expected = invgauss_ppf(rng_uniform(reference), 1.2_dp, -0.2_dp, 1.3_dp)
        call check_close('rvs invgauss', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = levy_rvs(state, -0.2_dp, 1.3_dp)
        expected = levy_ppf(rng_uniform(reference), -0.2_dp, 1.3_dp)
        call check_close('rvs levy', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = loglaplace_rvs(state, 1.7_dp, -0.2_dp, 1.3_dp)
        expected = loglaplace_ppf(rng_uniform(reference), 1.7_dp, -0.2_dp, 1.3_dp)
        call check_close('rvs loglaplace', actual, expected, 0.0_dp, 0.0_dp, failures)

        call invgamma_fit(xi, [2.3_dp, 0.1_dp, 1.2_dp], &
            [2.3_dp, 0.1_dp, 1.2_dp], result)
        call check_true('fit invgamma fixed success', result%success, failures)
        call invgauss_fit(xi, [1.4_dp, 0.1_dp, 1.2_dp], &
            [1.4_dp, 0.1_dp, 1.2_dp], result)
        call check_true('fit invgauss fixed success', result%success, failures)
        call levy_fit(xi, [0.1_dp, 1.2_dp], [0.1_dp, 1.2_dp], result)
        call check_true('fit levy fixed success', result%success, failures)
        call loglaplace_fit(xl, [1.8_dp, 0.1_dp, 1.2_dp], &
            [1.8_dp, 0.1_dp, 1.2_dp], result)
        call check_true('fit loglaplace fixed success', result%success, failures)

        call check_close('C ABI invgamma cdf', &
            scifort_invgamma_cdf_f64(1.2_dp, 2.2_dp, -0.2_dp, 1.3_dp), &
            invgamma_cdf(1.2_dp, 2.2_dp, -0.2_dp, 1.3_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI invgauss ppf', &
            scifort_invgauss_ppf_f64(0.4_dp, 1.2_dp, -0.2_dp, 1.3_dp), &
            invgauss_ppf(0.4_dp, 1.2_dp, -0.2_dp, 1.3_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI levy pdf', scifort_levy_pdf_f64(1.2_dp, -0.2_dp, 1.3_dp), &
            levy_pdf(1.2_dp, -0.2_dp, 1.3_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI loglaplace cdf', &
            scifort_loglaplace_cdf_f64(1.2_dp, 1.7_dp, -0.2_dp, 1.3_dp), &
            loglaplace_cdf(1.2_dp, 1.7_dp, -0.2_dp, 1.3_dp), 0.0_dp, 0.0_dp, failures)
    end subroutine test_random_fit_c_api

    pure function centered(plus, minus, h) result(value)
        real(dp), intent(in) :: plus !! objective at theta + h
        real(dp), intent(in) :: minus !! objective at theta - h
        real(dp), intent(in) :: h !! positive centered-difference half step
        real(dp) :: value

        value = (plus - minus) / (2.0_dp * h)
    end function centered

end program test_inverse_positive
