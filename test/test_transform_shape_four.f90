! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_transform_shape_four
    use scifort_c_api
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state, rng_uniform
    use scifort_stats
    use test_checks, only : check_close, check_true
    implicit none

    integer :: failures
    real(dp), parameter :: probs(5) = [1.0e-6_dp, 0.1_dp, 0.5_dp, 0.9_dp, 1.0_dp - 1.0e-6_dp]

    failures = 0
    call test_reference(failures)
    call test_scores(failures)
    call test_rng_fit_c_api(failures)

    if (failures == 0) then
        print '(a)', 'test_transform_shape_four: PASS'
    else
        print '(a,1x,i0)', 'test_transform_shape_four: FAIL', failures
        error stop 1
    end if
contains
    subroutine test_reference(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: gh_x(5) = [0.05_dp, 0.2_dp, 0.5_dp, 0.9_dp, 1.2_dp]
        real(dp), parameter :: gh_pdf(5) = [0.51779931585316696_dp, 0.57469960519501351_dp, &
            0.70076166800219752_dp, 0.84721054127966644_dp, 0.79209842448979595_dp]
        real(dp), parameter :: gh_cdf(5) = [0.025442492082888356_dp, 0.10731579131039697_dp, &
            0.29834515514273607_dp, 0.61078981610128924_dp, 0.86401881007491099_dp]
        real(dp), parameter :: gh_ppf(5) = [1.9999986000662773e-6_dp, 0.1872140300593772_dp, &
            0.76648134746861463_dp, 1.2466958284595799_dp, 1.4285159428655991_dp]
        real(dp), parameter :: ep_x(5) = [0.05_dp, 0.2_dp, 0.7_dp, 1.2_dp, 1.8_dp]
        real(dp), parameter :: ep_pdf(5) = [0.52910732565872798_dp, 0.79580200304483995_dp, &
            0.91270486385073957_dp, 0.38010458027267058_dp, 0.0069135481774016162_dp]
        real(dp), parameter :: ep_cdf(5) = [0.020353113980914453_dp, 0.12308439898287141_dp, &
            0.5834175596441149_dp, 0.92206136919445869_dp, 0.99947916765605305_dp]
        real(dp), parameter :: ep_ppf(5) = [2.424462017082642e-5_dp, 0.17034987098604404_dp, &
            0.61058699072227951_dp, 1.1466514897361162_dp, 2.1442858111002856_dp]
        real(dp), parameter :: ew_x(5) = [0.05_dp, 0.2_dp, 0.7_dp, 1.2_dp, 2.0_dp]
        real(dp), parameter :: ew_pdf(5) = [0.037840368766994217_dp, 0.23891550865439065_dp, &
            0.75989611322929473_dp, 0.61458628477614896_dp, 0.1477046751237355_dp]
        real(dp), parameter :: ew_cdf(5) = [0.00079741173289160693_dp, 0.020741982712062623_dp, &
            0.29720934392573206_dp, 0.66124753401613157_dp, 0.94608744260009037_dp]
        real(dp), parameter :: ew_ppf(5) = [0.0030129936245567421_dp, 0.40433698112404437_dp, &
            0.96446515902749763_dp, 1.763891760952184_dp, 4.7528197669283196_dp]
        real(dp), parameter :: pl_x(5) = [0.1_dp, 0.3_dp, 0.8_dp, 1.5_dp, 3.0_dp]
        real(dp), parameter :: pl_pdf(5) = [0.11873591426841582_dp, 0.77643145063031715_dp, &
            0.70232582434651536_dp, 0.2426595130815658_dp, 0.0282858111503261_dp]
        real(dp), parameter :: pl_cdf(5) = [0.0029978380539853917_dp, 0.097589295596525558_dp, &
            0.5237498896007694_dp, 0.83061583391986837_dp, 0.97529060650666177_dp]
        real(dp), parameter :: pl_ppf(5) = [0.020904814580402234_dp, 0.30309228296237739_dp, &
            0.76689986660244824_dp, 1.8778731920516332_dp, 19.593799974461465_dp]
        integer :: i
        do i = 1, 5
            call check_close('genhalflogistic pdf', genhalflogistic_pdf(gh_x(i), 0.7_dp), &
                gh_pdf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            call check_close('genhalflogistic cdf', genhalflogistic_cdf(gh_x(i), 0.7_dp), &
                gh_cdf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            call check_close('genhalflogistic ppf', genhalflogistic_ppf(probs(i), 0.7_dp), &
                gh_ppf(i), 3.0e-10_dp, 3.0e-9_dp, failures)
            call check_close('exponpow pdf', exponpow_pdf(ep_x(i), 1.3_dp), ep_pdf(i), &
                5.0e-12_dp, 5.0e-11_dp, failures)
            call check_close('exponpow cdf', exponpow_cdf(ep_x(i), 1.3_dp), ep_cdf(i), &
                5.0e-12_dp, 5.0e-11_dp, failures)
            call check_close('exponpow ppf', exponpow_ppf(probs(i), 1.3_dp), ep_ppf(i), &
                3.0e-10_dp, 3.0e-9_dp, failures)
            call check_close('exponweib pdf', exponweib_pdf(ew_x(i), 1.4_dp, 1.7_dp), &
                ew_pdf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            call check_close('exponweib cdf', exponweib_cdf(ew_x(i), 1.4_dp, 1.7_dp), &
                ew_cdf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            call check_close('exponweib ppf', exponweib_ppf(probs(i), 1.4_dp, 1.7_dp), &
                ew_ppf(i), 3.0e-10_dp, 3.0e-9_dp, failures)
            call check_close('powerlognorm pdf', powerlognorm_pdf(pl_x(i), 1.5_dp, 0.8_dp), &
                pl_pdf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            call check_close('powerlognorm cdf', powerlognorm_cdf(pl_x(i), 1.5_dp, 0.8_dp), &
                pl_cdf(i), 5.0e-12_dp, 5.0e-11_dp, failures)
            call check_close('powerlognorm ppf', powerlognorm_ppf(probs(i), 1.5_dp, 0.8_dp), &
                pl_ppf(i), 3.0e-9_dp, 3.0e-8_dp, failures)
        end do
    end subroutine test_reference

    subroutine test_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: x(4) = [0.2_dp, 0.5_dp, 0.8_dp, 1.2_dp]
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp), parameter :: tol = 3.0e-4_dp
        real(dp) :: q3(3), q4(4)
        q3 = genhalflogistic_score(x, 0.7_dp, 0.05_dp, 1.1_dp)
        call check_close('genhalflogistic score c', q3(1), centered( &
            genhalflogistic_loglikelihood(x, 0.7_dp+h, 0.05_dp, 1.1_dp), &
            genhalflogistic_loglikelihood(x, 0.7_dp-h, 0.05_dp, 1.1_dp), h), tol, tol, failures)
        call check_close('genhalflogistic score loc', q3(2), centered( &
            genhalflogistic_loglikelihood(x, 0.7_dp, 0.05_dp+h, 1.1_dp), &
            genhalflogistic_loglikelihood(x, 0.7_dp, 0.05_dp-h, 1.1_dp), h), tol, tol, failures)
        q3 = exponpow_score(x, 1.3_dp, 0.05_dp, 1.1_dp)
        call check_close('exponpow score b', q3(1), centered( &
            exponpow_loglikelihood(x, 1.3_dp+h, 0.05_dp, 1.1_dp), &
            exponpow_loglikelihood(x, 1.3_dp-h, 0.05_dp, 1.1_dp), h), tol, tol, failures)
        q4 = exponweib_score(x, 1.4_dp, 1.7_dp, 0.05_dp, 1.1_dp)
        call check_close('exponweib score a', q4(1), centered( &
            exponweib_loglikelihood(x, 1.4_dp+h, 1.7_dp, 0.05_dp, 1.1_dp), &
            exponweib_loglikelihood(x, 1.4_dp-h, 1.7_dp, 0.05_dp, 1.1_dp), h), tol, tol, failures)
        call check_close('exponweib score c', q4(2), centered( &
            exponweib_loglikelihood(x, 1.4_dp, 1.7_dp+h, 0.05_dp, 1.1_dp), &
            exponweib_loglikelihood(x, 1.4_dp, 1.7_dp-h, 0.05_dp, 1.1_dp), h), tol, tol, failures)
        q4 = powerlognorm_score(x, 1.5_dp, 0.8_dp, 0.05_dp, 1.1_dp)
        call check_close('powerlognorm score c', q4(1), centered( &
            powerlognorm_loglikelihood(x, 1.5_dp+h, 0.8_dp, 0.05_dp, 1.1_dp), &
            powerlognorm_loglikelihood(x, 1.5_dp-h, 0.8_dp, 0.05_dp, 1.1_dp), h), tol, tol, failures)
        call check_close('powerlognorm score s', q4(2), centered( &
            powerlognorm_loglikelihood(x, 1.5_dp, 0.8_dp+h, 0.05_dp, 1.1_dp), &
            powerlognorm_loglikelihood(x, 1.5_dp, 0.8_dp-h, 0.05_dp, 1.1_dp), h), tol, tol, failures)
    end subroutine test_scores

    subroutine test_rng_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: x(4) = [0.2_dp, 0.5_dp, 0.8_dp, 1.2_dp]
        real(dp) :: actual, expected
        type(rng_state) :: reference, state
        type(fit_result) :: result
        call rng_seed(state, 123456)
        call rng_seed(reference, 123456)
        actual = genhalflogistic_rvs(state, 0.7_dp, 0.05_dp, 1.1_dp)
        expected = genhalflogistic_ppf(rng_uniform(reference), 0.7_dp, 0.05_dp, 1.1_dp)
        call check_close('genhalflogistic rvs', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = exponpow_rvs(state, 1.3_dp, 0.05_dp, 1.1_dp)
        expected = exponpow_ppf(rng_uniform(reference), 1.3_dp, 0.05_dp, 1.1_dp)
        call check_close('exponpow rvs', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = exponweib_rvs(state, 1.4_dp, 1.7_dp, 0.05_dp, 1.1_dp)
        expected = exponweib_ppf(rng_uniform(reference), 1.4_dp, 1.7_dp, 0.05_dp, 1.1_dp)
        call check_close('exponweib rvs', actual, expected, 0.0_dp, 0.0_dp, failures)
        actual = powerlognorm_rvs(state, 1.5_dp, 0.8_dp, 0.05_dp, 1.1_dp)
        expected = powerlognorm_ppf(rng_uniform(reference), 1.5_dp, 0.8_dp, 0.05_dp, 1.1_dp)
        call check_close('powerlognorm rvs', actual, expected, 0.0_dp, 0.0_dp, failures)
        call genhalflogistic_fit(x, [0.7_dp,0.05_dp,1.1_dp], [0.7_dp,0.05_dp,1.1_dp], result)
        call check_true('genhalflogistic fixed fit', result%success, failures)
        call exponweib_fit(x, [1.4_dp,1.7_dp,0.05_dp,1.1_dp], &
            [1.4_dp,1.7_dp,0.05_dp,1.1_dp], result)
        call check_true('exponweib fixed fit', result%success, failures)
        call check_close('C ABI genhalflogistic cdf', &
            scifort_genhalflogistic_cdf_f64(0.5_dp,0.7_dp,0.0_dp,1.0_dp), &
            genhalflogistic_cdf(0.5_dp,0.7_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI exponpow pdf', scifort_exponpow_pdf_f64(0.7_dp,1.3_dp,0.0_dp,1.0_dp), &
            exponpow_pdf(0.7_dp,1.3_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI exponweib ppf', &
            scifort_exponweib_ppf_f64(0.5_dp,1.4_dp,1.7_dp,0.0_dp,1.0_dp), &
            exponweib_ppf(0.5_dp,1.4_dp,1.7_dp), 0.0_dp, 0.0_dp, failures)
        call check_close('C ABI powerlognorm cdf', &
            scifort_powerlognorm_cdf_f64(0.8_dp,1.5_dp,0.8_dp,0.0_dp,1.0_dp), &
            powerlognorm_cdf(0.8_dp,1.5_dp,0.8_dp), 0.0_dp, 0.0_dp, failures)
    end subroutine test_rng_fit_c_api

    pure function centered(plus, minus, h) result(value)
        real(dp), intent(in) :: plus !! objective at theta+h
        real(dp), intent(in) :: minus !! objective at theta-h
        real(dp), intent(in) :: h !! positive finite-difference step
        real(dp) :: value
        value = (plus - minus) / (2.0_dp * h)
    end function centered
end program test_transform_shape_four
