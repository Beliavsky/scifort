! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Compact full-stack regression coverage for the chi and Maxwell additions.

program test_chi_maxwell
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
        print '(a)', 'test_chi_maxwell: PASS'
    else
        print '(a,1x,i0)', 'test_chi_maxwell: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference_values(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: x(5) = [0.1_dp, 0.5_dp, 1.0_dp, 2.0_dp, 4.0_dp]
        real(dp), parameter :: chi_pdf_ref(5) = [ &
            0.017478997285437642_dp, 0.23913887645004395_dp, &
            0.53399937211163817_dp, 0.38712392836095172_dp, &
            0.0031176985239521351_dp]
        real(dp), parameter :: chi_cdf_ref(5) = [ &
            0.00064874971154539802_dp, 0.04673104184090323_dp, &
            0.24693359869598669_dp, 0.77922479651655252_dp, &
            0.99918880790448228_dp]
        real(dp), parameter :: max_pdf_ref(5) = [ &
            0.0079390509495402377_dp, 0.17603266338214973_dp, &
            0.48394144903828673_dp, 0.4319277321055045_dp, &
            0.0042825672244763318_dp]
        real(dp), parameter :: max_cdf_ref(5) = [ &
            0.00026516505865561009_dp, 0.030859595783726757_dp, &
            0.19874804309879915_dp, 0.73853587005088883_dp, &
            0.99886601571021472_dp]
        integer :: i
        character(len=48) :: label

        do i = 1, size(x)
            write(label, '(a,i0)') 'chi pdf ref ', i
            call check_close(trim(label), chi_pdf(x(i), 2.7_dp), chi_pdf_ref(i), &
                8.0e-15_dp, 8.0e-14_dp, failures)
            write(label, '(a,i0)') 'chi cdf ref ', i
            call check_close(trim(label), chi_cdf(x(i), 2.7_dp), chi_cdf_ref(i), &
                8.0e-15_dp, 8.0e-14_dp, failures)
            write(label, '(a,i0)') 'maxwell pdf ref ', i
            call check_close(trim(label), maxwell_pdf(x(i)), max_pdf_ref(i), &
                8.0e-15_dp, 8.0e-14_dp, failures)
            write(label, '(a,i0)') 'maxwell cdf ref ', i
            call check_close(trim(label), maxwell_cdf(x(i)), max_cdf_ref(i), &
                8.0e-15_dp, 8.0e-14_dp, failures)
        end do

        call check_close('chi df1 endpoint pdf', chi_pdf(0.0_dp, 1.0_dp), &
            0.797884560802865406_dp, 3.0e-15_dp, 3.0e-14_dp, failures)
        call check_close('maxwell endpoint pdf', maxwell_pdf(0.0_dp), &
            0.0_dp, 0.0_dp, 0.0_dp, failures)
    end subroutine test_reference_values

    subroutine test_identities_and_round_trips(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: probs(5) = [ &
            1.0e-8_dp, 0.1_dp, 0.5_dp, 0.9_dp, 1.0_dp - 1.0e-8_dp]
        real(dp), parameter :: points(4) = [0.2_dp, 0.7_dp, 1.3_dp, 3.0_dp]
        integer :: i
        real(dp) :: x

        do i = 1, size(points)
            call check_close('maxwell equals chi3 pdf', maxwell_pdf(points(i), -0.2_dp, 1.7_dp), &
                chi_pdf(points(i), 3.0_dp, -0.2_dp, 1.7_dp), 0.0_dp, 0.0_dp, failures)
            call check_close('maxwell equals chi3 sf', maxwell_sf(points(i), -0.2_dp, 1.7_dp), &
                chi_sf(points(i), 3.0_dp, -0.2_dp, 1.7_dp), 0.0_dp, 0.0_dp, failures)
        end do

        do i = 1, size(probs)
            x = chi_ppf(probs(i), 2.7_dp, -0.2_dp, 1.7_dp)
            call check_close('chi ppf/cdf', chi_cdf(x, 2.7_dp, -0.2_dp, 1.7_dp), &
                probs(i), 8.0e-14_dp, 5.0e-10_dp, failures)
            x = maxwell_isf(probs(i), -0.2_dp, 1.7_dp)
            call check_close('maxwell isf/sf', maxwell_sf(x, -0.2_dp, 1.7_dp), &
                probs(i), 8.0e-14_dp, 5.0e-10_dp, failures)
        end do

        call check_true('chi invalid df nan', ieee_is_nan(chi_pdf(1.0_dp, 0.0_dp)), failures)
        call check_true('maxwell invalid scale nan', &
            ieee_is_nan(maxwell_cdf(1.0_dp, scale=0.0_dp)), failures)
    end subroutine test_identities_and_round_trips

    subroutine test_likelihood_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: xc(5) = [0.2_dp, 0.7_dp, 1.1_dp, 1.8_dp, 2.6_dp]
        real(dp), parameter :: xm(5) = [0.3_dp, 0.8_dp, 1.4_dp, 2.0_dp, 3.0_dp]
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp), parameter :: tol = 1.0e-5_dp
        real(dp) :: sc(3)
        real(dp) :: sm(2)

        call check_close('chi loglike sum', chi_loglikelihood(xc, 2.7_dp, -0.2_dp, 1.4_dp), &
            sum(chi_logpdf(xc, 2.7_dp, -0.2_dp, 1.4_dp)), &
            3.0e-15_dp, 3.0e-15_dp, failures)
        call check_close('maxwell loglike sum', maxwell_loglikelihood(xm, -0.2_dp, 1.4_dp), &
            sum(maxwell_logpdf(xm, -0.2_dp, 1.4_dp)), 3.0e-15_dp, 3.0e-15_dp, failures)

        sc = chi_score(xc, 2.7_dp, -0.2_dp, 1.4_dp)
        call check_close('chi score df', sc(1), centered( &
            chi_loglikelihood(xc, 2.7_dp + h, -0.2_dp, 1.4_dp), &
            chi_loglikelihood(xc, 2.7_dp - h, -0.2_dp, 1.4_dp), h), tol, tol, failures)
        call check_close('chi score loc', sc(2), centered( &
            chi_loglikelihood(xc, 2.7_dp, -0.2_dp + h, 1.4_dp), &
            chi_loglikelihood(xc, 2.7_dp, -0.2_dp - h, 1.4_dp), h), tol, tol, failures)
        call check_close('chi score scale', sc(3), centered( &
            chi_loglikelihood(xc, 2.7_dp, -0.2_dp, 1.4_dp + h), &
            chi_loglikelihood(xc, 2.7_dp, -0.2_dp, 1.4_dp - h), h), tol, tol, failures)

        sm = maxwell_score(xm, -0.2_dp, 1.4_dp)
        call check_close('maxwell score loc', sm(1), centered( &
            maxwell_loglikelihood(xm, -0.2_dp + h, 1.4_dp), &
            maxwell_loglikelihood(xm, -0.2_dp - h, 1.4_dp), h), tol, tol, failures)
        call check_close('maxwell score scale', sm(2), centered( &
            maxwell_loglikelihood(xm, -0.2_dp, 1.4_dp + h), &
            maxwell_loglikelihood(xm, -0.2_dp, 1.4_dp - h), h), tol, tol, failures)
    end subroutine test_likelihood_scores

    subroutine test_random_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: xc(5) = [0.2_dp, 0.7_dp, 1.1_dp, 1.8_dp, 2.6_dp]
        real(dp), parameter :: xm(5) = [0.3_dp, 0.8_dp, 1.4_dp, 2.0_dp, 3.0_dp]
        real(dp) :: actual
        real(dp) :: expected
        type(fit_result) :: result
        type(rng_state) :: reference
        type(rng_state) :: state

        call rng_seed(state, 86420)
        call rng_seed(reference, 86420)
        actual = chi_rvs(state, 2.7_dp, -0.2_dp, 1.8_dp)
        expected = chi_ppf(rng_uniform(reference), 2.7_dp, -0.2_dp, 1.8_dp)
        call check_close('rvs chi', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = maxwell_rvs(state, 0.1_dp, 1.3_dp)
        expected = maxwell_ppf(rng_uniform(reference), 0.1_dp, 1.3_dp)
        call check_close('rvs maxwell', actual, expected, 0.0_dp, 0.0_dp, failures)

        call chi_fit(xc, [2.7_dp, -0.2_dp, 1.4_dp], &
            [2.7_dp, -0.2_dp, 1.4_dp], result)
        call check_true('fit chi fixed success', result%success, failures)
        call maxwell_fit(xm, [-0.2_dp, 1.4_dp], [-0.2_dp, 1.4_dp], result)
        call check_true('fit maxwell fixed success', result%success, failures)

        call check_close('C ABI chi cdf', &
            scifort_chi_cdf_f64(1.2_dp, 2.7_dp, -0.2_dp, 1.8_dp), &
            chi_cdf(1.2_dp, 2.7_dp, -0.2_dp, 1.8_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI maxwell ppf', &
            scifort_maxwell_ppf_f64(0.4_dp, -0.2_dp, 1.8_dp), &
            maxwell_ppf(0.4_dp, -0.2_dp, 1.8_dp), 0.0_dp, 0.0_dp, failures)
    end subroutine test_random_fit_c_api

    pure function centered(plus, minus, h) result(value)
        real(dp), intent(in) :: plus !! objective at theta + h
        real(dp), intent(in) :: minus !! objective at theta - h
        real(dp), intent(in) :: h !! positive centered-difference half step
        real(dp) :: value

        value = (plus - minus) / (2.0_dp * h)
    end function centered

end program test_chi_maxwell
