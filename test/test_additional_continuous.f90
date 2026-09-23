! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Compact cross-layer regression coverage for half-Cauchy and Lomax.

program test_additional_continuous
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
    call test_round_trips(failures)
    call test_likelihood_scores(failures)
    call test_random_fit_c_api(failures)

    if (failures == 0) then
        print '(a)', 'test_additional_continuous: PASS'
    else
        print '(a,1x,i0)', 'test_additional_continuous: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference_values(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: x(5) = [0.1_dp, 0.5_dp, 1.0_dp, 3.0_dp, 10.0_dp]
        real(dp), parameter :: hpdf(5) = [ &
            0.630316606304536009_dp, 0.509295817894065084_dp, &
            0.318309886183790691_dp, 0.0636619772367581355_dp, &
            0.0063031660630453604_dp]
        real(dp), parameter :: hcdf(5) = [ &
            0.0634510348611071473_dp, 0.295167235300866526_dp, &
            0.5_dp, 0.795167235300866637_dp, 0.936548965138892964_dp]
        real(dp), parameter :: lpdf(5) = [ &
            1.67931425674992485_dp, 0.603429699134749420_dp, &
            0.233522563952417744_dp, 0.0237099077717004454_dp, &
            0.000841650866731334267_dp]
        real(dp), parameter :: lcdf(5) = [ &
            0.196849703293513895_dp, 0.606458891868641770_dp, &
            0.796936900910941048_dp, 0.958765377788347029_dp, &
            0.995974713246067567_dp]
        integer :: i
        character(len=48) :: label

        do i = 1, size(x)
            write(label, '(a,i0)') 'halfcauchy pdf ref ', i
            call check_close(trim(label), halfcauchy_pdf(x(i)), hpdf(i), &
                5.0e-15_dp, 5.0e-14_dp, failures)
            write(label, '(a,i0)') 'halfcauchy cdf ref ', i
            call check_close(trim(label), halfcauchy_cdf(x(i)), hcdf(i), &
                5.0e-15_dp, 5.0e-14_dp, failures)
            write(label, '(a,i0)') 'lomax pdf ref ', i
            call check_close(trim(label), lomax_pdf(x(i), 2.3_dp), lpdf(i), &
                5.0e-15_dp, 5.0e-14_dp, failures)
            write(label, '(a,i0)') 'lomax cdf ref ', i
            call check_close(trim(label), lomax_cdf(x(i), 2.3_dp), lcdf(i), &
                5.0e-15_dp, 5.0e-14_dp, failures)
        end do

        call check_close('halfcauchy endpoint pdf', halfcauchy_pdf(0.0_dp), &
            0.636619772367581382_dp, 4.0e-15_dp, 4.0e-14_dp, failures)
        call check_close('lomax endpoint pdf', lomax_pdf(0.0_dp, 2.3_dp), &
            2.3_dp, 0.0_dp, 0.0_dp, failures)
    end subroutine test_reference_values

    subroutine test_round_trips(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: probs(5) = [1.0e-8_dp, 0.1_dp, 0.5_dp, 0.9_dp, 1.0_dp - 1.0e-8_dp]
        integer :: i
        real(dp) :: x

        do i = 1, size(probs)
            x = halfcauchy_ppf(probs(i), -0.2_dp, 1.7_dp)
            call check_close('halfcauchy ppf/cdf', halfcauchy_cdf(x, -0.2_dp, 1.7_dp), &
                probs(i), 8.0e-15_dp, 8.0e-10_dp, failures)
            x = lomax_isf(probs(i), 2.3_dp, -0.2_dp, 1.7_dp)
            call check_close('lomax isf/sf', lomax_sf(x, 2.3_dp, -0.2_dp, 1.7_dp), &
                probs(i), 8.0e-15_dp, 8.0e-10_dp, failures)
        end do

        call check_close('halfcauchy tail complement', &
            halfcauchy_cdf(1.2_dp) + halfcauchy_sf(1.2_dp), &
            1.0_dp, 4.0e-15_dp, 4.0e-15_dp, failures)
        call check_close('lomax tail complement', lomax_cdf(1.2_dp, 2.3_dp) + &
            lomax_sf(1.2_dp, 2.3_dp), 1.0_dp, 4.0e-15_dp, 4.0e-15_dp, failures)
        call check_true('halfcauchy invalid scale nan', &
            ieee_is_nan(halfcauchy_cdf(1.0_dp, scale=0.0_dp)), failures)
        call check_true('lomax invalid shape nan', ieee_is_nan(lomax_pdf(1.0_dp, 0.0_dp)), failures)
    end subroutine test_round_trips

    subroutine test_likelihood_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: xh(5) = [0.1_dp, 0.4_dp, 0.8_dp, 1.2_dp, 2.0_dp]
        real(dp), parameter :: xl(5) = [0.1_dp, 0.4_dp, 0.8_dp, 1.4_dp, 2.4_dp]
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp), parameter :: tol = 1.0e-5_dp
        real(dp) :: sh(2)
        real(dp) :: sl(3)

        call check_close('halfcauchy loglike sum', &
            halfcauchy_loglikelihood(xh, -0.2_dp, 1.4_dp), &
            sum(halfcauchy_logpdf(xh, -0.2_dp, 1.4_dp)), &
            3.0e-15_dp, 3.0e-15_dp, failures)
        call check_close('lomax loglike sum', lomax_loglikelihood(xl, 2.3_dp, -0.2_dp, 1.4_dp), &
            sum(lomax_logpdf(xl, 2.3_dp, -0.2_dp, 1.4_dp)), &
            3.0e-15_dp, 3.0e-15_dp, failures)

        sh = halfcauchy_score(xh, -0.2_dp, 1.4_dp)
        call check_close('halfcauchy score loc', sh(1), centered( &
            halfcauchy_loglikelihood(xh, -0.2_dp + h, 1.4_dp), &
            halfcauchy_loglikelihood(xh, -0.2_dp - h, 1.4_dp), h), tol, tol, failures)
        call check_close('halfcauchy score scale', sh(2), centered( &
            halfcauchy_loglikelihood(xh, -0.2_dp, 1.4_dp + h), &
            halfcauchy_loglikelihood(xh, -0.2_dp, 1.4_dp - h), h), tol, tol, failures)

        sl = lomax_score(xl, 2.3_dp, -0.2_dp, 1.4_dp)
        call check_close('lomax score shape', sl(1), centered( &
            lomax_loglikelihood(xl, 2.3_dp + h, -0.2_dp, 1.4_dp), &
            lomax_loglikelihood(xl, 2.3_dp - h, -0.2_dp, 1.4_dp), h), tol, tol, failures)
        call check_close('lomax score loc', sl(2), centered( &
            lomax_loglikelihood(xl, 2.3_dp, -0.2_dp + h, 1.4_dp), &
            lomax_loglikelihood(xl, 2.3_dp, -0.2_dp - h, 1.4_dp), h), tol, tol, failures)
        call check_close('lomax score scale', sl(3), centered( &
            lomax_loglikelihood(xl, 2.3_dp, -0.2_dp, 1.4_dp + h), &
            lomax_loglikelihood(xl, 2.3_dp, -0.2_dp, 1.4_dp - h), h), tol, tol, failures)
    end subroutine test_likelihood_scores

    subroutine test_random_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: xh(5) = [0.1_dp, 0.4_dp, 0.8_dp, 1.2_dp, 2.0_dp]
        real(dp), parameter :: xl(5) = [0.1_dp, 0.4_dp, 0.8_dp, 1.4_dp, 2.4_dp]
        real(dp) :: actual
        real(dp) :: expected
        type(fit_result) :: result
        type(rng_state) :: reference
        type(rng_state) :: state

        call rng_seed(state, 86420)
        call rng_seed(reference, 86420)
        actual = halfcauchy_rvs(state, -0.2_dp, 1.8_dp)
        expected = halfcauchy_ppf(rng_uniform(reference), -0.2_dp, 1.8_dp)
        call check_close('rvs halfcauchy', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = lomax_rvs(state, 2.3_dp, -0.2_dp, 1.8_dp)
        expected = lomax_ppf(rng_uniform(reference), 2.3_dp, -0.2_dp, 1.8_dp)
        call check_close('rvs lomax', actual, expected, 0.0_dp, 0.0_dp, failures)

        call halfcauchy_fit(xh, [-0.2_dp, 1.4_dp], [-0.2_dp, 1.4_dp], result)
        call check_true('fit halfcauchy fixed success', result%success, failures)
        call lomax_fit(xl, [2.3_dp, -0.2_dp, 1.4_dp], [2.3_dp, -0.2_dp, 1.4_dp], result)
        call check_true('fit lomax fixed success', result%success, failures)

        call check_close('C ABI halfcauchy cdf', &
            scifort_halfcauchy_cdf_f64(0.4_dp, -0.2_dp, 1.8_dp), &
            halfcauchy_cdf(0.4_dp, -0.2_dp, 1.8_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI lomax ppf', scifort_lomax_ppf_f64(0.4_dp, 2.3_dp, -0.2_dp, 1.8_dp), &
            lomax_ppf(0.4_dp, 2.3_dp, -0.2_dp, 1.8_dp), 0.0_dp, 0.0_dp, failures)
    end subroutine test_random_fit_c_api

    pure function centered(plus, minus, h) result(value)
        real(dp), intent(in) :: plus !! objective at theta + h
        real(dp), intent(in) :: minus !! objective at theta - h
        real(dp), intent(in) :: h !! positive centered-difference half step
        real(dp) :: value

        value = (plus - minus) / (2.0_dp * h)
    end function centered

end program test_additional_continuous
