! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Compact full-stack regression coverage for Anglit and Moyal additions.

program test_anglit_moyal
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
        print '(a)', 'test_anglit_moyal: PASS'
    else
        print '(a,1x,i0)', 'test_anglit_moyal: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference_values(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: xa(5) = [-0.7_dp, -0.3_dp, 0.0_dp, 0.3_dp, 0.7_dp]
        real(dp), parameter :: anglit_pdf_ref(5) = [ &
            0.16996714290024104_dp, 0.8253356149096783_dp, 1.0_dp, &
            0.8253356149096783_dp, 0.16996714290024104_dp]
        real(dp), parameter :: anglit_cdf_ref(5) = [ &
            0.007275135005769911_dp, 0.2176787633024823_dp, &
            0.4999999999999999_dp, 0.7823212366975176_dp, &
            0.99272486499423_dp]
        real(dp), parameter :: xm(5) = [-2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 4.0_dp]
        real(dp), parameter :: moyal_pdf_ref(5) = [ &
            0.02695823175881603_dp, 0.1689623369069973_dp, &
            0.24197072451914337_dp, 0.20131624406488796_dp, &
            0.0534987840888832_dp]
        real(dp), parameter :: moyal_cdf_ref(5) = [ &
            0.006562191672591345_dp, 0.0992047504111147_dp, &
            0.31731050786291415_dp, 0.5441624293623031_dp, &
            0.892346789695769_dp]
        integer :: i
        character(len=48) :: label

        do i = 1, size(xa)
            write(label, '(a,i0)') 'anglit pdf ref ', i
            call check_close(trim(label), anglit_pdf(xa(i)), anglit_pdf_ref(i), &
                8.0e-15_dp, 8.0e-14_dp, failures)
            write(label, '(a,i0)') 'anglit cdf ref ', i
            call check_close(trim(label), anglit_cdf(xa(i)), anglit_cdf_ref(i), &
                8.0e-15_dp, 8.0e-14_dp, failures)
            write(label, '(a,i0)') 'moyal pdf ref ', i
            call check_close(trim(label), moyal_pdf(xm(i)), moyal_pdf_ref(i), &
                8.0e-15_dp, 8.0e-14_dp, failures)
            write(label, '(a,i0)') 'moyal cdf ref ', i
            call check_close(trim(label), moyal_cdf(xm(i)), moyal_cdf_ref(i), &
                8.0e-15_dp, 8.0e-14_dp, failures)
        end do

        call check_close('moyal ppf 0.1', moyal_ppf(0.1_dp), &
            -0.9953027993634334_dp, 2.0e-14_dp, 2.0e-13_dp, failures)
        call check_close('moyal ppf 0.9', moyal_ppf(0.9_dp), &
            4.148329427633069_dp, 2.0e-14_dp, 2.0e-13_dp, failures)
    end subroutine test_reference_values

    subroutine test_identities_and_round_trips(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: probs(5) = [ &
            1.0e-10_dp, 0.1_dp, 0.5_dp, 0.9_dp, 1.0_dp - 1.0e-10_dp]
        real(dp), parameter :: points(5) = [-1.0_dp, -0.4_dp, 0.0_dp, 0.7_dp, 2.0_dp]
        integer :: i
        real(dp) :: x

        do i = 1, size(points)
            call check_close('anglit cdf+sf', &
                anglit_cdf(points(i), 0.2_dp, 1.3_dp) + anglit_sf(points(i), 0.2_dp, 1.3_dp), &
                1.0_dp, 4.0e-15_dp, 4.0e-15_dp, failures)
            call check_close('moyal cdf+sf', &
                moyal_cdf(points(i), 0.2_dp, 1.3_dp) + moyal_sf(points(i), 0.2_dp, 1.3_dp), &
                1.0_dp, 4.0e-15_dp, 4.0e-15_dp, failures)
        end do

        do i = 1, size(probs)
            x = anglit_ppf(probs(i), -0.2_dp, 1.7_dp)
            call check_close('anglit ppf/cdf', anglit_cdf(x, -0.2_dp, 1.7_dp), &
                probs(i), 2.0e-13_dp, 2.0e-8_dp, failures)
            x = moyal_isf(probs(i), -0.2_dp, 1.7_dp)
            call check_close('moyal isf/sf', moyal_sf(x, -0.2_dp, 1.7_dp), &
                probs(i), 2.0e-13_dp, 2.0e-8_dp, failures)
        end do

        call check_close('moyal far-tail logsf', moyal_logsf(100.0_dp), &
            log(moyal_sf(100.0_dp)), 1.0e-13_dp, 1.0e-13_dp, failures)
        call check_true('anglit invalid scale nan', &
            ieee_is_nan(anglit_pdf(0.0_dp, scale=0.0_dp)), failures)
        call check_true('moyal invalid scale nan', &
            ieee_is_nan(moyal_cdf(0.0_dp, scale=0.0_dp)), failures)
    end subroutine test_identities_and_round_trips

    subroutine test_likelihood_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: xa(5) = [-0.7_dp, -0.3_dp, 0.0_dp, 0.3_dp, 0.7_dp]
        real(dp), parameter :: xm(5) = [-1.4_dp, -0.5_dp, 0.2_dp, 1.1_dp, 2.4_dp]
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp), parameter :: tol = 1.0e-5_dp
        real(dp) :: sa(2)
        real(dp) :: sm(2)

        call check_close('anglit loglike sum', anglit_loglikelihood(xa, 0.0_dp, 1.2_dp), &
            sum(anglit_logpdf(xa, 0.0_dp, 1.2_dp)), 3.0e-15_dp, 3.0e-15_dp, failures)
        call check_close('moyal loglike sum', moyal_loglikelihood(xm, 0.1_dp, 1.4_dp), &
            sum(moyal_logpdf(xm, 0.1_dp, 1.4_dp)), 3.0e-15_dp, 3.0e-15_dp, failures)

        sa = anglit_score(xa, 0.0_dp, 1.2_dp)
        call check_close('anglit score loc', sa(1), centered( &
            anglit_loglikelihood(xa, h, 1.2_dp), &
            anglit_loglikelihood(xa, -h, 1.2_dp), h), tol, tol, failures)
        call check_close('anglit score scale', sa(2), centered( &
            anglit_loglikelihood(xa, 0.0_dp, 1.2_dp + h), &
            anglit_loglikelihood(xa, 0.0_dp, 1.2_dp - h), h), tol, tol, failures)

        sm = moyal_score(xm, 0.1_dp, 1.4_dp)
        call check_close('moyal score loc', sm(1), centered( &
            moyal_loglikelihood(xm, 0.1_dp + h, 1.4_dp), &
            moyal_loglikelihood(xm, 0.1_dp - h, 1.4_dp), h), tol, tol, failures)
        call check_close('moyal score scale', sm(2), centered( &
            moyal_loglikelihood(xm, 0.1_dp, 1.4_dp + h), &
            moyal_loglikelihood(xm, 0.1_dp, 1.4_dp - h), h), tol, tol, failures)
    end subroutine test_likelihood_scores

    subroutine test_random_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: xa(5) = [-0.7_dp, -0.3_dp, 0.0_dp, 0.3_dp, 0.7_dp]
        real(dp), parameter :: xm(5) = [-1.4_dp, -0.5_dp, 0.2_dp, 1.1_dp, 2.4_dp]
        real(dp) :: actual
        real(dp) :: expected
        type(fit_result) :: result
        type(rng_state) :: reference
        type(rng_state) :: state

        call rng_seed(state, 86420)
        call rng_seed(reference, 86420)
        actual = anglit_rvs(state, -0.2_dp, 1.8_dp)
        expected = anglit_ppf(rng_uniform(reference), -0.2_dp, 1.8_dp)
        call check_close('rvs anglit', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = moyal_rvs(state, 0.1_dp, 1.3_dp)
        expected = moyal_ppf(rng_uniform(reference), 0.1_dp, 1.3_dp)
        call check_close('rvs moyal', actual, expected, 0.0_dp, 0.0_dp, failures)

        call anglit_fit(xa, [0.0_dp, 1.2_dp], [0.0_dp, 1.2_dp], result)
        call check_true('fit anglit fixed success', result%success, failures)
        call moyal_fit(xm, [0.1_dp, 1.4_dp], [0.1_dp, 1.4_dp], result)
        call check_true('fit moyal fixed success', result%success, failures)

        call check_close('C ABI anglit cdf', &
            scifort_anglit_cdf_f64(0.2_dp, -0.2_dp, 1.8_dp), &
            anglit_cdf(0.2_dp, -0.2_dp, 1.8_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI moyal ppf', &
            scifort_moyal_ppf_f64(0.4_dp, -0.2_dp, 1.8_dp), &
            moyal_ppf(0.4_dp, -0.2_dp, 1.8_dp), 0.0_dp, 0.0_dp, failures)
    end subroutine test_random_fit_c_api

    pure function centered(plus, minus, h) result(value)
        real(dp), intent(in) :: plus !! objective at theta + h
        real(dp), intent(in) :: minus !! objective at theta - h
        real(dp), intent(in) :: h !! positive centered-difference half step
        real(dp) :: value

        value = (plus - minus) / (2.0_dp * h)
    end function centered

end program test_anglit_moyal
