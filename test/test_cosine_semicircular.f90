! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Compact full-stack regression coverage for cosine and semicircular additions.

program test_cosine_semicircular
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
    call test_identities_and_round_trips(failures)
    call test_likelihood_scores(failures)
    call test_random_fit_c_api(failures)

    if (failures == 0) then
        print '(a)', 'test_cosine_semicircular: PASS'
    else
        print '(a,1x,i0)', 'test_cosine_semicircular: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference_values(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: xc(5) = [-2.5_dp, -1.0_dp, 0.0_dp, 1.0_dp, 2.5_dp]
        real(dp), parameter :: cosine_pdf_ref(5) = [ &
            0.03164897655108783_dp, 0.24514672583475897_dp, &
            0.3183098861837907_dp, 0.24514672583475897_dp, &
            0.03164897655108783_dp]
        real(dp), parameter :: cosine_cdf_ref(5) = [ &
            0.00686284223331188_dp, 0.20692079020752274_dp, &
            0.5_dp, 0.7930792097924773_dp, 0.9931371577666881_dp]
        real(dp), parameter :: xs(5) = [-0.8_dp, -0.3_dp, 0.0_dp, 0.4_dp, 0.9_dp]
        real(dp), parameter :: semi_pdf_ref(5) = [ &
            0.38197186342054873_dp, 0.6072965572585683_dp, &
            0.6366197723675814_dp, 0.5834716591559995_dp, &
            0.2774961253210154_dp]
        real(dp), parameter :: semi_cdf_ref(5) = [ &
            0.052044019330913904_dp, 0.31191883239053647_dp, &
            0.5_dp, 0.7476842122656545_dp, 0.9813069632657507_dp]
        integer :: i
        character(len=48) :: label

        do i = 1, size(xc)
            write(label, '(a,i0)') 'cosine pdf ref ', i
            call check_close(trim(label), cosine_pdf(xc(i)), cosine_pdf_ref(i), &
                8.0e-15_dp, 8.0e-14_dp, failures)
            write(label, '(a,i0)') 'cosine cdf ref ', i
            call check_close(trim(label), cosine_cdf(xc(i)), cosine_cdf_ref(i), &
                8.0e-15_dp, 8.0e-14_dp, failures)
            write(label, '(a,i0)') 'semicircular pdf ref ', i
            call check_close(trim(label), semicircular_pdf(xs(i)), semi_pdf_ref(i), &
                8.0e-15_dp, 8.0e-14_dp, failures)
            write(label, '(a,i0)') 'semicircular cdf ref ', i
            call check_close(trim(label), semicircular_cdf(xs(i)), semi_cdf_ref(i), &
                8.0e-15_dp, 8.0e-14_dp, failures)
        end do
    end subroutine test_reference_values

    subroutine test_identities_and_round_trips(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: probs(5) = [ &
            1.0e-8_dp, 0.1_dp, 0.5_dp, 0.9_dp, 1.0_dp - 1.0e-8_dp]
        real(dp), parameter :: points(5) = [-1.1_dp, -0.4_dp, 0.0_dp, 0.7_dp, 1.4_dp]
        integer :: i
        real(dp) :: x

        do i = 1, size(points)
            call check_close('cosine cdf+sf', &
                cosine_cdf(points(i), 0.2_dp, 1.3_dp) + cosine_sf(points(i), 0.2_dp, 1.3_dp), &
                1.0_dp, 4.0e-15_dp, 4.0e-15_dp, failures)
            call check_close('semicircular cdf+sf', &
                semicircular_cdf(points(i), 0.2_dp, 1.6_dp) + &
                semicircular_sf(points(i), 0.2_dp, 1.6_dp), &
                1.0_dp, 4.0e-15_dp, 4.0e-15_dp, failures)
        end do

        do i = 1, size(probs)
            x = cosine_ppf(probs(i), -0.2_dp, 1.7_dp)
            call check_close('cosine ppf/cdf', cosine_cdf(x, -0.2_dp, 1.7_dp), &
                probs(i), 2.0e-13_dp, 2.0e-9_dp, failures)
            x = semicircular_isf(probs(i), -0.2_dp, 1.7_dp)
            call check_close('semicircular isf/sf', semicircular_sf(x, -0.2_dp, 1.7_dp), &
                probs(i), 2.0e-13_dp, 2.0e-9_dp, failures)
        end do

        call check_true('cosine invalid scale nan', &
            ieee_is_nan(cosine_pdf(0.0_dp, scale=0.0_dp)), failures)
        call check_true('semicircular invalid scale nan', &
            ieee_is_nan(semicircular_cdf(0.0_dp, scale=0.0_dp)), failures)
    end subroutine test_identities_and_round_trips

    subroutine test_likelihood_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: xc(5) = [-1.2_dp, -0.4_dp, 0.2_dp, 0.9_dp, 1.6_dp]
        real(dp), parameter :: xs(5) = [-0.8_dp, -0.3_dp, 0.2_dp, 0.7_dp, 1.0_dp]
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp), parameter :: tol = 1.0e-5_dp
        real(dp) :: sc(2)
        real(dp) :: ss(2)

        call check_close('cosine loglike sum', cosine_loglikelihood(xc, 0.1_dp, 1.4_dp), &
            sum(cosine_logpdf(xc, 0.1_dp, 1.4_dp)), 3.0e-15_dp, 3.0e-15_dp, failures)
        call check_close('semicircular loglike sum', &
            semicircular_loglikelihood(xs, 0.1_dp, 1.4_dp), &
            sum(semicircular_logpdf(xs, 0.1_dp, 1.4_dp)), &
            3.0e-15_dp, 3.0e-15_dp, failures)

        sc = cosine_score(xc, 0.1_dp, 1.4_dp)
        call check_close('cosine score loc', sc(1), centered( &
            cosine_loglikelihood(xc, 0.1_dp + h, 1.4_dp), &
            cosine_loglikelihood(xc, 0.1_dp - h, 1.4_dp), h), tol, tol, failures)
        call check_close('cosine score scale', sc(2), centered( &
            cosine_loglikelihood(xc, 0.1_dp, 1.4_dp + h), &
            cosine_loglikelihood(xc, 0.1_dp, 1.4_dp - h), h), tol, tol, failures)

        ss = semicircular_score(xs, 0.1_dp, 1.4_dp)
        call check_close('semicircular score loc', ss(1), centered( &
            semicircular_loglikelihood(xs, 0.1_dp + h, 1.4_dp), &
            semicircular_loglikelihood(xs, 0.1_dp - h, 1.4_dp), h), tol, tol, failures)
        call check_close('semicircular score scale', ss(2), centered( &
            semicircular_loglikelihood(xs, 0.1_dp, 1.4_dp + h), &
            semicircular_loglikelihood(xs, 0.1_dp, 1.4_dp - h), h), tol, tol, failures)
    end subroutine test_likelihood_scores

    subroutine test_random_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: xc(5) = [-1.2_dp, -0.4_dp, 0.2_dp, 0.9_dp, 1.6_dp]
        real(dp), parameter :: xs(5) = [-0.8_dp, -0.3_dp, 0.2_dp, 0.7_dp, 1.0_dp]
        real(dp) :: actual
        real(dp) :: expected
        type(fit_result) :: result
        type(rng_state) :: reference
        type(rng_state) :: state

        call rng_seed(state, 97531)
        call rng_seed(reference, 97531)
        actual = cosine_rvs(state, -0.2_dp, 1.8_dp)
        expected = cosine_ppf(rng_uniform(reference), -0.2_dp, 1.8_dp)
        call check_close('rvs cosine', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = semicircular_rvs(state, 0.1_dp, 1.3_dp)
        expected = semicircular_ppf(rng_uniform(reference), 0.1_dp, 1.3_dp)
        call check_close('rvs semicircular', actual, expected, 0.0_dp, 0.0_dp, failures)

        call cosine_fit(xc, [0.1_dp, 1.4_dp], [0.1_dp, 1.4_dp], result)
        call check_true('fit cosine fixed success', result%success, failures)
        call semicircular_fit(xs, [0.1_dp, 1.4_dp], [0.1_dp, 1.4_dp], result)
        call check_true('fit semicircular fixed success', result%success, failures)

        call check_close('C ABI cosine cdf', &
            scifort_cosine_cdf_f64(0.4_dp, -0.2_dp, 1.8_dp), &
            cosine_cdf(0.4_dp, -0.2_dp, 1.8_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI semicircular ppf', &
            scifort_semicircular_ppf_f64(0.4_dp, -0.2_dp, 1.8_dp), &
            semicircular_ppf(0.4_dp, -0.2_dp, 1.8_dp), 0.0_dp, 0.0_dp, failures)
    end subroutine test_random_fit_c_api

    pure function centered(plus, minus, h) result(value)
        real(dp), intent(in) :: plus !! objective at theta + h
        real(dp), intent(in) :: minus !! objective at theta - h
        real(dp), intent(in) :: h !! positive centered-difference half step
        real(dp) :: value

        value = (plus - minus) / (2.0_dp * h)
    end function centered

end program test_cosine_semicircular
