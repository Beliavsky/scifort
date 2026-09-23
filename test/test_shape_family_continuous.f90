! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Full-stack regression coverage for alpha, fatigue-life, generalized-logistic,
! and generalized-normal distributions.

program test_shape_family_continuous
    use scifort_c_api
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state, rng_uniform
    use scifort_stats
    use test_checks, only : check_close, check_true
    use test_shape_family_reference
    implicit none

    integer :: failures

    failures = 0
    call test_reference_values(failures)
    call test_tail_identities(failures)
    call test_likelihood_scores(failures)
    call test_random_fit_c_api(failures)

    if (failures == 0) then
        print '(a)', 'test_shape_family_continuous: PASS'
    else
        print '(a,1x,i0)', 'test_shape_family_continuous: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference_values(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        integer :: i
        character(len=64) :: label

        do i = 1, 5
            write(label, '(a,i0)') 'alpha pdf ref ', i
            call check_close(trim(label), alpha_pdf(al_x(i), reference_shape), &
                al_pdf(i), 4.0e-13_dp, 4.0e-12_dp, failures)
            write(label, '(a,i0)') 'alpha cdf ref ', i
            call check_close(trim(label), alpha_cdf(al_x(i), reference_shape), &
                al_cdf(i), 4.0e-13_dp, 4.0e-12_dp, failures)
            write(label, '(a,i0)') 'alpha sf ref ', i
            call check_close(trim(label), alpha_sf(al_x(i), reference_shape), &
                al_sf(i), 4.0e-13_dp, 4.0e-12_dp, failures)
            write(label, '(a,i0)') 'alpha logcdf ref ', i
            call check_close(trim(label), alpha_logcdf(al_x(i), reference_shape), &
                al_logcdf(i), 4.0e-13_dp, 4.0e-12_dp, failures)
            write(label, '(a,i0)') 'alpha logsf ref ', i
            call check_close(trim(label), alpha_logsf(al_x(i), reference_shape), &
                al_logsf(i), 4.0e-13_dp, 4.0e-12_dp, failures)

            write(label, '(a,i0)') 'fatiguelife pdf ref ', i
            call check_close(trim(label), fatiguelife_pdf(fl_x(i), reference_shape), &
                fl_pdf(i), 4.0e-13_dp, 4.0e-12_dp, failures)
            write(label, '(a,i0)') 'fatiguelife cdf ref ', i
            call check_close(trim(label), fatiguelife_cdf(fl_x(i), reference_shape), &
                fl_cdf(i), 4.0e-13_dp, 4.0e-12_dp, failures)
            write(label, '(a,i0)') 'fatiguelife sf ref ', i
            call check_close(trim(label), fatiguelife_sf(fl_x(i), reference_shape), &
                fl_sf(i), 4.0e-13_dp, 4.0e-12_dp, failures)
            write(label, '(a,i0)') 'fatiguelife logcdf ref ', i
            call check_close(trim(label), fatiguelife_logcdf(fl_x(i), reference_shape), &
                fl_logcdf(i), 4.0e-13_dp, 4.0e-12_dp, failures)
            write(label, '(a,i0)') 'fatiguelife logsf ref ', i
            call check_close(trim(label), fatiguelife_logsf(fl_x(i), reference_shape), &
                fl_logsf(i), 4.0e-13_dp, 4.0e-12_dp, failures)

            write(label, '(a,i0)') 'genlogistic pdf ref ', i
            call check_close(trim(label), genlogistic_pdf(gl_x(i), reference_shape), &
                gl_pdf(i), 4.0e-13_dp, 4.0e-12_dp, failures)
            write(label, '(a,i0)') 'genlogistic cdf ref ', i
            call check_close(trim(label), genlogistic_cdf(gl_x(i), reference_shape), &
                gl_cdf(i), 4.0e-13_dp, 4.0e-12_dp, failures)
            write(label, '(a,i0)') 'genlogistic sf ref ', i
            call check_close(trim(label), genlogistic_sf(gl_x(i), reference_shape), &
                gl_sf(i), 4.0e-13_dp, 4.0e-12_dp, failures)
            write(label, '(a,i0)') 'genlogistic logcdf ref ', i
            call check_close(trim(label), genlogistic_logcdf(gl_x(i), reference_shape), &
                gl_logcdf(i), 4.0e-13_dp, 4.0e-12_dp, failures)
            write(label, '(a,i0)') 'genlogistic logsf ref ', i
            call check_close(trim(label), genlogistic_logsf(gl_x(i), reference_shape), &
                gl_logsf(i), 4.0e-13_dp, 4.0e-12_dp, failures)

            write(label, '(a,i0)') 'gennorm pdf ref ', i
            call check_close(trim(label), gennorm_pdf(gn_x(i), reference_shape), &
                gn_pdf(i), 4.0e-12_dp, 4.0e-11_dp, failures)
            write(label, '(a,i0)') 'gennorm cdf ref ', i
            call check_close(trim(label), gennorm_cdf(gn_x(i), reference_shape), &
                gn_cdf(i), 4.0e-12_dp, 4.0e-11_dp, failures)
            write(label, '(a,i0)') 'gennorm sf ref ', i
            call check_close(trim(label), gennorm_sf(gn_x(i), reference_shape), &
                gn_sf(i), 4.0e-12_dp, 4.0e-11_dp, failures)
            write(label, '(a,i0)') 'gennorm logcdf ref ', i
            call check_close(trim(label), gennorm_logcdf(gn_x(i), reference_shape), &
                gn_logcdf(i), 4.0e-12_dp, 4.0e-11_dp, failures)
            write(label, '(a,i0)') 'gennorm logsf ref ', i
            call check_close(trim(label), gennorm_logsf(gn_x(i), reference_shape), &
                gn_logsf(i), 4.0e-12_dp, 4.0e-11_dp, failures)
        end do

        do i = 1, 5
            write(label, '(a,i0)') 'alpha ppf ref ', i
            call check_close(trim(label), alpha_ppf(reference_probs(i), reference_shape), &
                al_ppf(i), 2.0e-11_dp, 2.0e-10_dp, failures)
            write(label, '(a,i0)') 'alpha isf ref ', i
            call check_close(trim(label), alpha_isf(reference_probs(i), reference_shape), &
                al_isf(i), 2.0e-11_dp, 2.0e-10_dp, failures)
            write(label, '(a,i0)') 'fatiguelife ppf ref ', i
            call check_close(trim(label), fatiguelife_ppf(reference_probs(i), reference_shape), &
                fl_ppf(i), 2.0e-11_dp, 2.0e-10_dp, failures)
            write(label, '(a,i0)') 'fatiguelife isf ref ', i
            call check_close(trim(label), fatiguelife_isf(reference_probs(i), reference_shape), &
                fl_isf(i), 2.0e-11_dp, 2.0e-10_dp, failures)
            write(label, '(a,i0)') 'genlogistic ppf ref ', i
            call check_close(trim(label), genlogistic_ppf(reference_probs(i), reference_shape), &
                gl_ppf(i), 2.0e-11_dp, 2.0e-10_dp, failures)
            write(label, '(a,i0)') 'genlogistic isf ref ', i
            call check_close(trim(label), genlogistic_isf(reference_probs(i), reference_shape), &
                gl_isf(i), 2.0e-11_dp, 2.0e-10_dp, failures)
            write(label, '(a,i0)') 'gennorm ppf ref ', i
            call check_close(trim(label), gennorm_ppf(reference_probs(i), reference_shape), &
                gn_ppf(i), 2.0e-10_dp, 2.0e-9_dp, failures)
            write(label, '(a,i0)') 'gennorm isf ref ', i
            call check_close(trim(label), gennorm_isf(reference_probs(i), reference_shape), &
                gn_isf(i), 2.0e-10_dp, 2.0e-9_dp, failures)
        end do
    end subroutine test_reference_values

    subroutine test_tail_identities(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: probs(5) = [1.0e-8_dp, 0.1_dp, 0.5_dp, 0.9_dp, &
            1.0_dp - 1.0e-8_dp]
        integer :: i
        real(dp) :: x

        do i = 1, size(probs)
            x = alpha_ppf(probs(i), 1.6_dp, -0.3_dp, 1.4_dp)
            call check_close('alpha ppf/cdf', alpha_cdf(x, 1.6_dp, -0.3_dp, 1.4_dp), &
                probs(i), 2.0e-10_dp, 2.0e-8_dp, failures)
            x = fatiguelife_isf(probs(i), 1.6_dp, -0.3_dp, 1.4_dp)
            call check_close('fatiguelife isf/sf', fatiguelife_sf(x, 1.6_dp, -0.3_dp, 1.4_dp), &
                probs(i), 2.0e-10_dp, 2.0e-8_dp, failures)
            x = genlogistic_ppf(probs(i), 1.6_dp, -0.3_dp, 1.4_dp)
            call check_close('genlogistic ppf/cdf', genlogistic_cdf(x, 1.6_dp, -0.3_dp, 1.4_dp), &
                probs(i), 2.0e-11_dp, 2.0e-9_dp, failures)
            x = gennorm_isf(probs(i), 1.6_dp, -0.3_dp, 1.4_dp)
            call check_close('gennorm isf/sf', gennorm_sf(x, 1.6_dp, -0.3_dp, 1.4_dp), &
                probs(i), 2.0e-9_dp, 2.0e-7_dp, failures)
        end do
    end subroutine test_tail_identities

    subroutine test_likelihood_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: xa(5) = [0.3_dp, 0.6_dp, 1.1_dp, 1.8_dp, 3.0_dp]
        real(dp), parameter :: xf(5) = [0.35_dp, 0.7_dp, 1.2_dp, 2.0_dp, 3.5_dp]
        real(dp), parameter :: xl(5) = [-2.0_dp, -0.5_dp, 0.4_dp, 1.5_dp, 3.0_dp]
        real(dp), parameter :: xn(5) = [-2.0_dp, -0.6_dp, 0.4_dp, 1.4_dp, 2.7_dp]
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp), parameter :: tol = 7.0e-5_dp
        real(dp) :: s3(3)

        s3 = alpha_score(xa, 1.7_dp, 0.1_dp, 1.2_dp)
        call compare_score('alpha', s3, &
            centered(alpha_loglikelihood(xa, 1.7_dp + h, 0.1_dp, 1.2_dp), &
                alpha_loglikelihood(xa, 1.7_dp - h, 0.1_dp, 1.2_dp), h), &
            centered(alpha_loglikelihood(xa, 1.7_dp, 0.1_dp + h, 1.2_dp), &
                alpha_loglikelihood(xa, 1.7_dp, 0.1_dp - h, 1.2_dp), h), &
            centered(alpha_loglikelihood(xa, 1.7_dp, 0.1_dp, 1.2_dp + h), &
                alpha_loglikelihood(xa, 1.7_dp, 0.1_dp, 1.2_dp - h), h), tol, failures)

        s3 = fatiguelife_score(xf, 1.7_dp, 0.1_dp, 1.2_dp)
        call compare_score('fatiguelife', s3, &
            centered(fatiguelife_loglikelihood(xf, 1.7_dp + h, 0.1_dp, 1.2_dp), &
                fatiguelife_loglikelihood(xf, 1.7_dp - h, 0.1_dp, 1.2_dp), h), &
            centered(fatiguelife_loglikelihood(xf, 1.7_dp, 0.1_dp + h, 1.2_dp), &
                fatiguelife_loglikelihood(xf, 1.7_dp, 0.1_dp - h, 1.2_dp), h), &
            centered(fatiguelife_loglikelihood(xf, 1.7_dp, 0.1_dp, 1.2_dp + h), &
                fatiguelife_loglikelihood(xf, 1.7_dp, 0.1_dp, 1.2_dp - h), h), tol, failures)

        s3 = genlogistic_score(xl, 1.7_dp, 0.1_dp, 1.2_dp)
        call compare_score('genlogistic', s3, &
            centered(genlogistic_loglikelihood(xl, 1.7_dp + h, 0.1_dp, 1.2_dp), &
                genlogistic_loglikelihood(xl, 1.7_dp - h, 0.1_dp, 1.2_dp), h), &
            centered(genlogistic_loglikelihood(xl, 1.7_dp, 0.1_dp + h, 1.2_dp), &
                genlogistic_loglikelihood(xl, 1.7_dp, 0.1_dp - h, 1.2_dp), h), &
            centered(genlogistic_loglikelihood(xl, 1.7_dp, 0.1_dp, 1.2_dp + h), &
                genlogistic_loglikelihood(xl, 1.7_dp, 0.1_dp, 1.2_dp - h), h), tol, failures)

        s3 = gennorm_score(xn, 1.7_dp, 0.1_dp, 1.2_dp)
        call compare_score('gennorm', s3, &
            centered(gennorm_loglikelihood(xn, 1.7_dp + h, 0.1_dp, 1.2_dp), &
                gennorm_loglikelihood(xn, 1.7_dp - h, 0.1_dp, 1.2_dp), h), &
            centered(gennorm_loglikelihood(xn, 1.7_dp, 0.1_dp + h, 1.2_dp), &
                gennorm_loglikelihood(xn, 1.7_dp, 0.1_dp - h, 1.2_dp), h), &
            centered(gennorm_loglikelihood(xn, 1.7_dp, 0.1_dp, 1.2_dp + h), &
                gennorm_loglikelihood(xn, 1.7_dp, 0.1_dp, 1.2_dp - h), h), tol, failures)
    end subroutine test_likelihood_scores

    subroutine compare_score(name, score, dshape, dloc, dscale, tol, failures)
        character(len=*), intent(in) :: name !! family name used in check labels
        real(dp), intent(in) :: score(3) !! analytic score vector
        real(dp), intent(in) :: dshape !! finite-difference shape derivative
        real(dp), intent(in) :: dloc !! finite-difference location derivative
        real(dp), intent(in) :: dscale !! finite-difference scale derivative
        real(dp), intent(in) :: tol !! absolute and relative comparison tolerance
        integer, intent(inout) :: failures !! running count of failed checks
        call check_close(name // ' score shape', score(1), dshape, tol, tol, failures)
        call check_close(name // ' score loc', score(2), dloc, tol, tol, failures)
        call check_close(name // ' score scale', score(3), dscale, tol, tol, failures)
    end subroutine compare_score

    subroutine test_random_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: xa(5) = [0.3_dp, 0.6_dp, 1.1_dp, 1.8_dp, 3.0_dp]
        real(dp), parameter :: xf(5) = [0.35_dp, 0.7_dp, 1.2_dp, 2.0_dp, 3.5_dp]
        real(dp), parameter :: xl(5) = [-2.0_dp, -0.5_dp, 0.4_dp, 1.5_dp, 3.0_dp]
        real(dp), parameter :: xn(5) = [-2.0_dp, -0.6_dp, 0.4_dp, 1.4_dp, 2.7_dp]
        real(dp) :: actual
        real(dp) :: expected
        type(fit_result) :: result
        type(rng_state) :: reference
        type(rng_state) :: state

        call rng_seed(state, 86420)
        call rng_seed(reference, 86420)
        actual = alpha_rvs(state, 1.6_dp, -0.2_dp, 1.3_dp)
        expected = alpha_ppf(rng_uniform(reference), 1.6_dp, -0.2_dp, 1.3_dp)
        call check_close('rvs alpha', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = fatiguelife_rvs(state, 1.6_dp, -0.2_dp, 1.3_dp)
        expected = fatiguelife_ppf(rng_uniform(reference), 1.6_dp, -0.2_dp, 1.3_dp)
        call check_close('rvs fatiguelife', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = genlogistic_rvs(state, 1.6_dp, -0.2_dp, 1.3_dp)
        expected = genlogistic_ppf(rng_uniform(reference), 1.6_dp, -0.2_dp, 1.3_dp)
        call check_close('rvs genlogistic', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = gennorm_rvs(state, 1.6_dp, -0.2_dp, 1.3_dp)
        expected = gennorm_ppf(rng_uniform(reference), 1.6_dp, -0.2_dp, 1.3_dp)
        call check_close('rvs gennorm', actual, expected, 0.0_dp, 0.0_dp, failures)

        call alpha_fit(xa, [1.7_dp, 0.1_dp, 1.2_dp], [1.7_dp, 0.1_dp, 1.2_dp], result)
        call check_true('fit alpha fixed success', result%success, failures)
        call fatiguelife_fit(xf, [1.7_dp, 0.1_dp, 1.2_dp], [1.7_dp, 0.1_dp, 1.2_dp], result)
        call check_true('fit fatiguelife fixed success', result%success, failures)
        call genlogistic_fit(xl, [1.7_dp, 0.1_dp, 1.2_dp], [1.7_dp, 0.1_dp, 1.2_dp], result)
        call check_true('fit genlogistic fixed success', result%success, failures)
        call gennorm_fit(xn, [1.7_dp, 0.1_dp, 1.2_dp], [1.7_dp, 0.1_dp, 1.2_dp], result)
        call check_true('fit gennorm fixed success', result%success, failures)

        call check_close('C ABI alpha cdf', &
            scifort_alpha_cdf_f64(0.8_dp, 1.7_dp, -0.2_dp, 1.3_dp), &
            alpha_cdf(0.8_dp, 1.7_dp, -0.2_dp, 1.3_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI fatiguelife ppf', &
            scifort_fatiguelife_ppf_f64(0.4_dp, 1.7_dp, -0.2_dp, 1.3_dp), &
            fatiguelife_ppf(0.4_dp, 1.7_dp, -0.2_dp, 1.3_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI genlogistic pdf', &
            scifort_genlogistic_pdf_f64(0.3_dp, 1.7_dp, -0.2_dp, 1.3_dp), &
            genlogistic_pdf(0.3_dp, 1.7_dp, -0.2_dp, 1.3_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI gennorm cdf', &
            scifort_gennorm_cdf_f64(0.3_dp, 1.7_dp, -0.2_dp, 1.3_dp), &
            gennorm_cdf(0.3_dp, 1.7_dp, -0.2_dp, 1.3_dp), 0.0_dp, 0.0_dp, failures)
    end subroutine test_random_fit_c_api

    pure function centered(plus, minus, h) result(value)
        real(dp), intent(in) :: plus !! objective at theta+h
        real(dp), intent(in) :: minus !! objective at theta-h
        real(dp), intent(in) :: h !! positive finite-difference step
        real(dp) :: value
        value = (plus - minus) / (2.0_dp * h)
    end function centered

end program test_shape_family_continuous
