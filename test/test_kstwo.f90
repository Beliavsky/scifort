! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_kstwo
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_c_api
    use scifort_fit, only : fit_result
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state, rng_uniform
    use scifort_stats
    use test_checks, only : check_close, check_true
    use test_kstwo_reference
    implicit none

    integer :: failures

    failures = 0
    call test_reference(failures)
    call test_structure(failures)
    call test_score(failures)
    call test_rng_fit_c_api(failures)

    if (failures == 0) then
        print '(a)', 'test_kstwo: PASS'
    else
        print '(a,1x,i0)', 'test_kstwo: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        integer :: i
        real(dp) :: n

        do i = 1, size(kstwo_x_ref)
            n = real(kstwo_n_ref(i), dp)
            call check_close('kstwo pdf', kstwo_pdf(kstwo_x_ref(i), n), &
                kstwo_pdf_ref(i), 5.0e-8_dp, 5.0e-8_dp, failures)
            call check_close('kstwo logpdf', kstwo_logpdf(kstwo_x_ref(i), n), &
                kstwo_logpdf_ref(i), 5.0e-8_dp, 5.0e-8_dp, failures)
            call check_close('kstwo cdf', kstwo_cdf(kstwo_x_ref(i), n), &
                kstwo_cdf_ref(i), 2.0e-12_dp, 2.0e-12_dp, failures)
            call check_close('kstwo sf', kstwo_sf(kstwo_x_ref(i), n), &
                kstwo_sf_ref(i), 2.0e-12_dp, 2.0e-12_dp, failures)
            call check_close('kstwo logcdf', kstwo_logcdf(kstwo_x_ref(i), n), &
                kstwo_logcdf_ref(i), 2.0e-11_dp, 2.0e-11_dp, failures)
            call check_close('kstwo logsf', kstwo_logsf(kstwo_x_ref(i), n), &
                kstwo_logsf_ref(i), 2.0e-11_dp, 2.0e-11_dp, failures)
        end do
        do i = 1, size(kstwo_prob_ref)
            n = real(kstwo_qn_ref(i), dp)
            call check_close('kstwo ppf', kstwo_ppf(kstwo_prob_ref(i), n), &
                kstwo_ppf_ref(i), 3.0e-9_dp, 3.0e-9_dp, failures)
            call check_close('kstwo isf', kstwo_isf(kstwo_prob_ref(i), n), &
                kstwo_isf_ref(i), 3.0e-9_dp, 3.0e-9_dp, failures)
        end do
    end subroutine test_reference

    subroutine test_structure(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp) :: x, p

        call check_close('kstwo n=1 cdf', kstwo_cdf(0.75_dp, 1.0_dp), &
            0.5_dp, 1.0e-14_dp, 1.0e-14_dp, failures)
        call check_close('kstwo n=1 pdf', kstwo_pdf(0.75_dp, 1.0_dp), &
            2.0_dp, 1.0e-14_dp, 1.0e-14_dp, failures)
        call check_close('kstwo ppf lower endpoint', kstwo_ppf(0.0_dp, 20.0_dp), &
            0.025_dp, 0.0_dp, 0.0_dp, failures)
        call check_close('kstwo ppf upper endpoint', kstwo_ppf(1.0_dp, 20.0_dp), &
            1.0_dp, 0.0_dp, 0.0_dp, failures)
        call check_close('kstwo shift/scale cdf', kstwo_cdf(0.7_dp, 20.0_dp, 0.1_dp, 2.0_dp), &
            kstwo_cdf(0.3_dp, 20.0_dp), 2.0e-14_dp, 2.0e-14_dp, failures)
        call check_close('kstwo shift/scale pdf', kstwo_pdf(0.7_dp, 20.0_dp, 0.1_dp, 2.0_dp), &
            0.5_dp * kstwo_pdf(0.3_dp, 20.0_dp), 2.0e-12_dp, 2.0e-12_dp, failures)
        call check_close('kstwo one-sided tail identity', kstwo_sf(0.6_dp, 20.0_dp), &
            2.0_dp * ksone_sf(0.6_dp, 20.0_dp), 1.0e-15_dp, 1.0e-15_dp, failures)

        p = 1.0e-8_dp
        x = kstwo_isf(p, 200.0_dp)
        call check_close('kstwo deep isf log-roundtrip', kstwo_logsf(x, 200.0_dp), &
            log(p), 2.0e-12_dp, 2.0e-12_dp, failures)
        call check_true('kstwo invalid noninteger n', ieee_is_nan(kstwo_cdf(0.2_dp, 2.5_dp)), failures)
        call check_true('kstwo invalid n zero', ieee_is_nan(kstwo_pdf(0.2_dp, 0.0_dp)), failures)
    end subroutine test_structure

    subroutine test_score(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: data(4) = [0.131_dp, 0.173_dp, 0.227_dp, 0.311_dp]
        real(dp), parameter :: n = 20.0_dp
        real(dp), parameter :: h = 2.0e-6_dp
        real(dp) :: fd
        real(dp) :: score(2)

        score = kstwo_score(data, n, 0.01_dp, 1.1_dp)
        fd = centered(kstwo_loglikelihood(data, n, 0.01_dp + h, 1.1_dp), &
            kstwo_loglikelihood(data, n, 0.01_dp - h, 1.1_dp), h)
        call check_close('kstwo score loc', score(1), fd, 5.0e-3_dp, 5.0e-4_dp, failures)
        fd = centered(kstwo_loglikelihood(data, n, 0.01_dp, 1.1_dp + h), &
            kstwo_loglikelihood(data, n, 0.01_dp, 1.1_dp - h), h)
        call check_close('kstwo score scale', score(2), fd, 5.0e-3_dp, 5.0e-4_dp, failures)
    end subroutine test_score

    subroutine test_rng_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: data(4) = [0.131_dp, 0.173_dp, 0.227_dp, 0.311_dp]
        type(rng_state) :: reference
        type(rng_state) :: state
        type(fit_result) :: result
        real(dp) :: actual, expected

        call rng_seed(state, 24680)
        call rng_seed(reference, 24680)
        actual = kstwo_rvs(state, 20.0_dp, 0.01_dp, 1.1_dp)
        expected = kstwo_ppf(rng_uniform(reference), 20.0_dp, 0.01_dp, 1.1_dp)
        call check_close('kstwo rvs', actual, expected, 0.0_dp, 0.0_dp, failures)

        call rng_seed(state, 13579)
        call rng_seed(reference, 13579)
        actual = kstwo_rvs(state, 2.5_dp)
        call check_true('invalid kstwo rvs NaN', ieee_is_nan(actual), failures)
        call check_close('invalid kstwo rvs preserves state', rng_uniform(state), &
            rng_uniform(reference), 0.0_dp, 0.0_dp, failures)

        call kstwo_fit(data, [20.0_dp, 0.01_dp, 1.1_dp], &
            [20.0_dp, 0.01_dp, 1.1_dp], result)
        call check_true('kstwo fixed fit', result%success, failures)

        call check_close('C ABI kstwo pdf', scifort_kstwo_pdf_f64(0.2_dp, 20.0_dp, 0.0_dp, 1.0_dp), &
            kstwo_pdf(0.2_dp, 20.0_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI kstwo cdf', scifort_kstwo_cdf_f64(0.2_dp, 20.0_dp, 0.0_dp, 1.0_dp), &
            kstwo_cdf(0.2_dp, 20.0_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI kstwo ppf', scifort_kstwo_ppf_f64(0.7_dp, 20.0_dp, 0.0_dp, 1.0_dp), &
            kstwo_ppf(0.7_dp, 20.0_dp), 0.0_dp, 0.0_dp, failures)
    end subroutine test_rng_fit_c_api

    pure function centered(plus_value, minus_value, h) result(value)
        real(dp), intent(in) :: plus_value !! function value at theta+h
        real(dp), intent(in) :: minus_value !! function value at theta-h
        real(dp), intent(in) :: h !! positive centered-difference step
        real(dp) :: value
        value = (plus_value - minus_value) / (2.0_dp * h)
    end function centered

end program test_kstwo
