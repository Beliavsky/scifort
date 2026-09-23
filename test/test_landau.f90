! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_landau
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_c_api
    use scifort_fit, only : fit_result
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state, rng_uniform
    use scifort_stats
    use test_checks, only : check_close, check_true
    use test_landau_reference
    implicit none

    integer :: failures

    failures = 0
    call test_reference(failures)
    call test_structure_and_tails(failures)
    call test_score(failures)
    call test_rng_fit_c_api(failures)

    if (failures == 0) then
        print '(a)', 'test_landau: PASS'
    else
        print '(a,1x,i0)', 'test_landau: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        integer :: i

        do i = 1, size(landau_x_ref)
            call check_close('landau pdf', landau_pdf(landau_x_ref(i)), &
                landau_pdf_ref(i), 3.0e-13_dp, 2.0e-8_dp, failures)
            call check_close('landau logpdf', landau_logpdf(landau_x_ref(i)), &
                landau_logpdf_ref(i), 2.0e-8_dp, 2.0e-9_dp, failures)
            call check_close('landau cdf', landau_cdf(landau_x_ref(i)), &
                landau_cdf_ref(i), 3.0e-13_dp, 2.0e-8_dp, failures)
            call check_close('landau sf', landau_sf(landau_x_ref(i)), &
                landau_sf_ref(i), 3.0e-13_dp, 2.0e-8_dp, failures)
            call check_close('landau logcdf', landau_logcdf(landau_x_ref(i)), &
                landau_logcdf_ref(i), 2.0e-8_dp, 2.0e-9_dp, failures)
            call check_close('landau logsf', landau_logsf(landau_x_ref(i)), &
                landau_logsf_ref(i), 2.0e-10_dp, 2.0e-9_dp, failures)
        end do
        do i = 1, size(landau_prob_ref)
            call check_close('landau ppf', landau_ppf(landau_prob_ref(i)), &
                landau_ppf_ref(i), 2.0e-8_dp, 3.0e-10_dp, failures)
            call check_close('landau isf', landau_isf(landau_prob_ref(i)), &
                landau_isf_ref(i), 2.0e-8_dp, 3.0e-10_dp, failures)
        end do
    end subroutine test_reference

    subroutine test_structure_and_tails(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp) :: x
        real(dp) :: p

        x = 1.7_dp
        call check_close('landau tails sum', landau_cdf(x) + landau_sf(x), &
            1.0_dp, 2.0e-14_dp, 2.0e-14_dp, failures)
        call check_close('landau shifted scaled cdf', landau_cdf(2.3_dp, 0.3_dp, 2.0_dp), &
            landau_cdf(1.0_dp), 2.0e-14_dp, 2.0e-14_dp, failures)
        call check_close('landau shifted scaled pdf', landau_pdf(2.3_dp, 0.3_dp, 2.0_dp), &
            0.5_dp * landau_pdf(1.0_dp), 2.0e-14_dp, 2.0e-14_dp, failures)

        p = 1.0e-100_dp
        x = landau_ppf(p)
        call check_true('landau deep lower ppf finite', ieee_is_finite(x), failures)
        call check_close('landau deep lower roundtrip', landau_logcdf(x), log(p), &
            2.0e-9_dp, 2.0e-10_dp, failures)
        p = 1.0e-12_dp
        x = landau_isf(p)
        call check_true('landau deep upper isf finite', ieee_is_finite(x), failures)
        call check_close('landau deep upper roundtrip', landau_logsf(x), log(p), &
            2.0e-9_dp, 2.0e-10_dp, failures)

        call check_true('landau invalid scale', ieee_is_nan(landau_pdf(0.0_dp, 0.0_dp, 0.0_dp)), failures)
        call check_true('landau ppf zero infinite', .not. ieee_is_finite(landau_ppf(0.0_dp)), failures)
        call check_true('landau ppf one infinite', .not. ieee_is_finite(landau_ppf(1.0_dp)), failures)
    end subroutine test_structure_and_tails

    subroutine test_score(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: data(4) = [-1.0_dp, 0.2_dp, 1.5_dp, 5.0_dp]
        real(dp), parameter :: h = 2.0e-6_dp
        real(dp) :: fd
        real(dp) :: score(2)

        score = landau_score(data, 0.1_dp, 1.2_dp)
        fd = centered(landau_loglikelihood(data, 0.1_dp + h, 1.2_dp), &
            landau_loglikelihood(data, 0.1_dp - h, 1.2_dp), h)
        call check_close('landau score loc', score(1), fd, 2.0e-6_dp, 2.0e-6_dp, failures)
        fd = centered(landau_loglikelihood(data, 0.1_dp, 1.2_dp + h), &
            landau_loglikelihood(data, 0.1_dp, 1.2_dp - h), h)
        call check_close('landau score scale', score(2), fd, 2.0e-6_dp, 2.0e-6_dp, failures)
    end subroutine test_score

    subroutine test_rng_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: data(4) = [-1.0_dp, 0.2_dp, 1.5_dp, 5.0_dp]
        type(rng_state) :: reference
        type(rng_state) :: state
        type(fit_result) :: result
        real(dp) :: actual
        real(dp) :: expected

        call rng_seed(state, 24680)
        call rng_seed(reference, 24680)
        actual = landau_rvs(state, 0.1_dp, 1.2_dp)
        expected = landau_ppf(rng_uniform(reference), 0.1_dp, 1.2_dp)
        call check_close('landau rvs', actual, expected, 0.0_dp, 0.0_dp, failures)

        call rng_seed(state, 13579)
        call rng_seed(reference, 13579)
        actual = landau_rvs(state, 0.0_dp, 0.0_dp)
        call check_true('invalid landau rvs NaN', ieee_is_nan(actual), failures)
        call check_close('invalid landau rvs preserves state', rng_uniform(state), &
            rng_uniform(reference), 0.0_dp, 0.0_dp, failures)

        call landau_fit(data, [0.1_dp, 1.2_dp], [0.1_dp, 1.2_dp], result)
        call check_true('landau fixed fit', result%success, failures)

        call check_close('C ABI landau pdf', scifort_landau_pdf_f64(0.3_dp, 0.1_dp, 1.2_dp), &
            landau_pdf(0.3_dp, 0.1_dp, 1.2_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI landau cdf', scifort_landau_cdf_f64(0.3_dp, 0.1_dp, 1.2_dp), &
            landau_cdf(0.3_dp, 0.1_dp, 1.2_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI landau ppf', scifort_landau_ppf_f64(0.7_dp, 0.1_dp, 1.2_dp), &
            landau_ppf(0.7_dp, 0.1_dp, 1.2_dp), 0.0_dp, 0.0_dp, failures)
    end subroutine test_rng_fit_c_api

    pure function centered(plus_value, minus_value, h) result(value)
        real(dp), intent(in) :: plus_value !! function value at theta+h
        real(dp), intent(in) :: minus_value !! function value at theta-h
        real(dp), intent(in) :: h !! positive centered-difference step
        real(dp) :: value

        value = (plus_value - minus_value) / (2.0_dp * h)
    end function centered

end program test_landau
