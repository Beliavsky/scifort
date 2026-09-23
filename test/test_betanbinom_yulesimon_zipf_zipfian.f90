! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_betanbinom_yulesimon_zipf_zipfian
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_c_api
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state, rng_uniform
    use scifort_special, only : hurwitz_zeta, hurwitz_zeta_derivative
    use scifort_stats
    use test_betanbinom_yulesimon_zipf_zipfian_reference
    use test_checks, only : check_close, check_true
    implicit none

    integer :: failures

    failures = 0
    call test_reference(failures)
    call test_identities_endpoints(failures)
    call test_scores(failures)
    call test_rng_fit_c_api(failures)

    if (failures == 0) then
        print '(a)', 'test_betanbinom_yulesimon_zipf_zipfian: PASS'
    else
        print '(a,1x,i0)', 'test_betanbinom_yulesimon_zipf_zipfian: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        integer :: i

        do i = 1, size(reference_probs)
            call check_close('betanbinom pmf',betanbinom_pmf(bn_x(i),4.0_dp,2.7_dp,1.4_dp), &
                bn_pmf(i),3.0e-14_dp,4.0e-13_dp,failures)
            call check_close('betanbinom logpmf',betanbinom_logpmf(bn_x(i),4.0_dp,2.7_dp,1.4_dp), &
                bn_logpmf(i),3.0e-13_dp,4.0e-13_dp,failures)
            call check_close('betanbinom cdf',betanbinom_cdf(bn_x(i),4.0_dp,2.7_dp,1.4_dp), &
                bn_cdf(i),4.0e-14_dp,5.0e-13_dp,failures)
            call check_close('betanbinom sf',betanbinom_sf(bn_x(i),4.0_dp,2.7_dp,1.4_dp), &
                bn_sf(i),4.0e-14_dp,5.0e-10_dp,failures)
            call check_close('betanbinom logcdf',betanbinom_logcdf(bn_x(i),4.0_dp,2.7_dp,1.4_dp), &
                bn_logcdf(i),4.0e-13_dp,5.0e-13_dp,failures)
            call check_close('betanbinom logsf',betanbinom_logsf(bn_x(i),4.0_dp,2.7_dp,1.4_dp), &
                bn_logsf(i),3.0e-10_dp,5.0e-10_dp,failures)
            call check_close('betanbinom ppf',betanbinom_ppf(reference_probs(i),4.0_dp,2.7_dp,1.4_dp), &
                bn_ppf(i),0.0_dp,0.0_dp,failures)
            call check_close('betanbinom isf',betanbinom_isf(reference_probs(i),4.0_dp,2.7_dp,1.4_dp), &
                bn_isf(i),0.0_dp,0.0_dp,failures)

            call check_close('yulesimon pmf',yulesimon_pmf(ys_x(i),2.5_dp),ys_pmf(i), &
                3.0e-14_dp,3.0e-13_dp,failures)
            call check_close('yulesimon logpmf',yulesimon_logpmf(ys_x(i),2.5_dp),ys_logpmf(i), &
                3.0e-13_dp,3.0e-13_dp,failures)
            call check_close('yulesimon cdf',yulesimon_cdf(ys_x(i),2.5_dp),ys_cdf(i), &
                3.0e-14_dp,3.0e-13_dp,failures)
            call check_close('yulesimon sf',yulesimon_sf(ys_x(i),2.5_dp),ys_sf(i), &
                3.0e-14_dp,3.0e-13_dp,failures)
            call check_close('yulesimon logcdf',yulesimon_logcdf(ys_x(i),2.5_dp),ys_logcdf(i), &
                3.0e-13_dp,3.0e-13_dp,failures)
            call check_close('yulesimon logsf',yulesimon_logsf(ys_x(i),2.5_dp),ys_logsf(i), &
                3.0e-13_dp,3.0e-13_dp,failures)
            call check_close('yulesimon ppf',yulesimon_ppf(reference_probs(i),2.5_dp),ys_ppf(i), &
                0.0_dp,0.0_dp,failures)
            call check_close('yulesimon isf',yulesimon_isf(reference_probs(i),2.5_dp),ys_isf(i), &
                0.0_dp,0.0_dp,failures)

            call check_close('zipf pmf',zipf_pmf(zp_x(i),2.3_dp),zp_pmf(i), &
                3.0e-14_dp,3.0e-13_dp,failures)
            call check_close('zipf logpmf',zipf_logpmf(zp_x(i),2.3_dp),zp_logpmf(i), &
                3.0e-13_dp,3.0e-13_dp,failures)
            call check_close('zipf cdf',zipf_cdf(zp_x(i),2.3_dp),zp_cdf(i), &
                3.0e-14_dp,3.0e-13_dp,failures)
            call check_close('zipf sf',zipf_sf(zp_x(i),2.3_dp),zp_sf(i), &
                3.0e-14_dp,3.0e-13_dp,failures)
            call check_close('zipf logcdf',zipf_logcdf(zp_x(i),2.3_dp),zp_logcdf(i), &
                3.0e-13_dp,3.0e-13_dp,failures)
            call check_close('zipf logsf',zipf_logsf(zp_x(i),2.3_dp),zp_logsf(i), &
                3.0e-13_dp,3.0e-13_dp,failures)
            call check_close('zipf ppf',zipf_ppf(reference_probs(i),2.3_dp),zp_ppf(i), &
                0.0_dp,0.0_dp,failures)
            call check_close('zipf isf',zipf_isf(reference_probs(i),2.3_dp),zp_isf(i), &
                0.0_dp,0.0_dp,failures)

            call check_close('zipfian pmf',zipfian_pmf(zi_x(i),1.2_dp,50.0_dp),zi_pmf(i), &
                3.0e-14_dp,3.0e-13_dp,failures)
            call check_close('zipfian logpmf',zipfian_logpmf(zi_x(i),1.2_dp,50.0_dp),zi_logpmf(i), &
                3.0e-13_dp,3.0e-13_dp,failures)
            call check_close('zipfian cdf',zipfian_cdf(zi_x(i),1.2_dp,50.0_dp),zi_cdf(i), &
                3.0e-14_dp,3.0e-13_dp,failures)
            call check_close('zipfian sf',zipfian_sf(zi_x(i),1.2_dp,50.0_dp),zi_sf(i), &
                3.0e-14_dp,3.0e-13_dp,failures)
            call check_close('zipfian logcdf',zipfian_logcdf(zi_x(i),1.2_dp,50.0_dp),zi_logcdf(i), &
                3.0e-13_dp,3.0e-13_dp,failures)
            call check_close('zipfian logsf',zipfian_logsf(zi_x(i),1.2_dp,50.0_dp),zi_logsf(i), &
                3.0e-13_dp,3.0e-13_dp,failures)
            call check_close('zipfian ppf',zipfian_ppf(reference_probs(i),1.2_dp,50.0_dp),zi_ppf(i), &
                0.0_dp,0.0_dp,failures)
            call check_close('zipfian isf',zipfian_isf(reference_probs(i),1.2_dp,50.0_dp),zi_isf(i), &
                0.0_dp,0.0_dp,failures)
        end do
    end subroutine test_reference

    subroutine test_identities_endpoints(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp) :: norm
        integer :: k

        call check_close('zeta two',hurwitz_zeta(2.0_dp,1.0_dp), &
            1.6449340668482264365_dp,2.0e-15_dp,2.0e-15_dp,failures)
        call check_close('zeta derivative two',hurwitz_zeta_derivative(2.0_dp,1.0_dp), &
            -0.9375482543158437537_dp,3.0e-15_dp,3.0e-15_dp,failures)

        do k = 1, 20
            call check_close('Yule-Simon alpha one',yulesimon_pmf(k,1.0_dp), &
                1.0_dp/(real(k,dp)*real(k+1,dp)),3.0e-15_dp,3.0e-14_dp,failures)
        end do
        norm = sum([(zipfian_pmf(k,0.0_dp,20),k=1,20)])
        call check_close('zipfian a zero normalization',norm,1.0_dp,3.0e-15_dp,3.0e-15_dp,failures)
        call check_close('zipfian a zero uniform',zipfian_pmf(7,0.0_dp,20),0.05_dp, &
            3.0e-15_dp,3.0e-15_dp,failures)
        call check_close('betanbinom tails',betanbinom_cdf(10,4,2.7_dp,1.4_dp) + &
            betanbinom_sf(10,4,2.7_dp,1.4_dp),1.0_dp,2.0e-15_dp,2.0e-15_dp,failures)
        call check_close('zipf tails',zipf_cdf(100,2.3_dp)+zipf_sf(100,2.3_dp), &
            1.0_dp,3.0e-15_dp,3.0e-15_dp,failures)

        call check_close('betanbinom ppf zero',betanbinom_ppf(0.0_dp,4.0_dp,2.7_dp,1.4_dp), &
            -1.0_dp,0.0_dp,0.0_dp,failures)
        call check_close('yulesimon ppf zero',yulesimon_ppf(0.0_dp,2.5_dp),0.0_dp,0.0_dp,0.0_dp,failures)
        call check_close('zipf ppf zero',zipf_ppf(0.0_dp,2.3_dp),0.0_dp,0.0_dp,0.0_dp,failures)
        call check_close('zipfian ppf one',zipfian_ppf(1.0_dp,1.2_dp,50.0_dp),50.0_dp,0.0_dp,0.0_dp,failures)

        call check_true('betanbinom invalid a',ieee_is_nan(betanbinom_pmf(2,4,0.0_dp,1.4_dp)),failures)
        call check_true('yulesimon invalid alpha',ieee_is_nan(yulesimon_pmf(2,0.0_dp)),failures)
        call check_true('zipf invalid a',ieee_is_nan(zipf_pmf(2,1.0_dp)),failures)
        call check_true('zipfian invalid n',ieee_is_nan(zipfian_pmf(2,1.2_dp,0)),failures)

        call check_close('betanbinom shifted',betanbinom_pmf(3.5_dp,4.0_dp,2.7_dp,1.4_dp,0.5_dp), &
            betanbinom_pmf(3.0_dp,4.0_dp,2.7_dp,1.4_dp),0.0_dp,0.0_dp,failures)
        call check_close('yulesimon shifted',yulesimon_cdf(4.5_dp,2.5_dp,0.5_dp), &
            yulesimon_cdf(4.0_dp,2.5_dp),0.0_dp,0.0_dp,failures)
        call check_close('zipf shifted',zipf_sf(7.5_dp,2.3_dp,0.5_dp), &
            zipf_sf(7.0_dp,2.3_dp),0.0_dp,0.0_dp,failures)
        call check_close('zipfian shifted',zipfian_ppf(0.7_dp,1.2_dp,50.0_dp,0.5_dp), &
            zipfian_ppf(0.7_dp,1.2_dp,50.0_dp)+0.5_dp,0.0_dp,0.0_dp,failures)
    end subroutine test_identities_endpoints

    subroutine test_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp), parameter :: bd(5) = [0.0_dp,1.0_dp,2.0_dp,5.0_dp,9.0_dp]
        real(dp), parameter :: yd(5) = [1.0_dp,1.0_dp,2.0_dp,3.0_dp,8.0_dp]
        real(dp), parameter :: zd(5) = [1.0_dp,1.0_dp,2.0_dp,4.0_dp,10.0_dp]
        real(dp), parameter :: zid(5) = [1.0_dp,2.0_dp,4.0_dp,10.0_dp,30.0_dp]
        real(dp) :: s1(1), s2(2)

        s2 = betanbinom_score(bd,4.0_dp,2.7_dp,1.4_dp)
        call check_close('betanbinom score a',s2(1),centered( &
            betanbinom_loglikelihood(bd,4.0_dp,2.7_dp+h,1.4_dp), &
            betanbinom_loglikelihood(bd,4.0_dp,2.7_dp-h,1.4_dp),h),5.0e-8_dp,5.0e-8_dp,failures)
        call check_close('betanbinom score b',s2(2),centered( &
            betanbinom_loglikelihood(bd,4.0_dp,2.7_dp,1.4_dp+h), &
            betanbinom_loglikelihood(bd,4.0_dp,2.7_dp,1.4_dp-h),h),5.0e-8_dp,5.0e-8_dp,failures)
        s1 = yulesimon_score(yd,2.5_dp)
        call check_close('yulesimon score',s1(1),centered(yulesimon_loglikelihood(yd,2.5_dp+h), &
            yulesimon_loglikelihood(yd,2.5_dp-h),h),5.0e-8_dp,5.0e-8_dp,failures)
        s1 = zipf_score(zd,2.3_dp)
        call check_close('zipf score',s1(1),centered(zipf_loglikelihood(zd,2.3_dp+h), &
            zipf_loglikelihood(zd,2.3_dp-h),h),5.0e-8_dp,5.0e-8_dp,failures)
        s1 = zipfian_score(zid,1.2_dp,50.0_dp)
        call check_close('zipfian score',s1(1),centered(zipfian_loglikelihood(zid,1.2_dp+h,50.0_dp), &
            zipfian_loglikelihood(zid,1.2_dp-h,50.0_dp),h),5.0e-8_dp,5.0e-8_dp,failures)
    end subroutine test_scores

    subroutine test_rng_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: bd(4) = [0.0_dp,1.0_dp,3.0_dp,7.0_dp]
        real(dp), parameter :: yd(4) = [1.0_dp,1.0_dp,2.0_dp,5.0_dp]
        real(dp), parameter :: zd(4) = [1.0_dp,1.0_dp,3.0_dp,9.0_dp]
        real(dp), parameter :: zid(4) = [1.0_dp,2.0_dp,5.0_dp,20.0_dp]
        real(dp) :: actual, expected
        type(rng_state) :: reference, state
        type(fit_result) :: result

        call rng_seed(state,24680)
        call rng_seed(reference,24680)
        actual = betanbinom_rvs(state,4.0_dp,2.7_dp,1.4_dp)
        expected = betanbinom_ppf(rng_uniform(reference),4.0_dp,2.7_dp,1.4_dp)
        call check_close('betanbinom rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = yulesimon_rvs(state,2.5_dp)
        expected = yulesimon_ppf(rng_uniform(reference),2.5_dp)
        call check_close('yulesimon rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = zipf_rvs(state,2.3_dp)
        expected = zipf_ppf(rng_uniform(reference),2.3_dp)
        call check_close('zipf rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = zipfian_rvs(state,1.2_dp,50.0_dp)
        expected = zipfian_ppf(rng_uniform(reference),1.2_dp,50.0_dp)
        call check_close('zipfian rvs',actual,expected,0.0_dp,0.0_dp,failures)

        call rng_seed(state,13579)
        call rng_seed(reference,13579)
        actual = zipf_rvs(state,1.0_dp)
        call check_true('invalid zipf rvs NaN',ieee_is_nan(actual),failures)
        call check_close('invalid rvs preserves state',rng_uniform(state),rng_uniform(reference), &
            0.0_dp,0.0_dp,failures)

        call betanbinom_fit(bd,[4.0_dp,2.7_dp,1.4_dp,0.0_dp],[4.0_dp,2.7_dp,1.4_dp,0.0_dp],result)
        call check_true('betanbinom fixed fit',result%success,failures)
        call yulesimon_fit(yd,[2.5_dp,0.0_dp],[2.5_dp,0.0_dp],result)
        call check_true('yulesimon fixed fit',result%success,failures)
        call zipf_fit(zd,[2.3_dp,0.0_dp],[2.3_dp,0.0_dp],result)
        call check_true('zipf fixed fit',result%success,failures)
        call zipfian_fit(zid,[1.2_dp,50.0_dp,0.0_dp],[1.2_dp,50.0_dp,0.0_dp],result)
        call check_true('zipfian fixed fit',result%success,failures)

        call check_close('C ABI betanbinom pmf',scifort_betanbinom_pmf_f64(3.0_dp,4.0_dp,2.7_dp,1.4_dp,0.0_dp), &
            betanbinom_pmf(3.0_dp,4.0_dp,2.7_dp,1.4_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI betanbinom cdf',scifort_betanbinom_cdf_f64(3.0_dp,4.0_dp,2.7_dp,1.4_dp,0.0_dp), &
            betanbinom_cdf(3.0_dp,4.0_dp,2.7_dp,1.4_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI betanbinom ppf',scifort_betanbinom_ppf_f64(0.7_dp,4.0_dp,2.7_dp,1.4_dp,0.0_dp), &
            betanbinom_ppf(0.7_dp,4.0_dp,2.7_dp,1.4_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI yulesimon pmf',scifort_yulesimon_pmf_f64(3.0_dp,2.5_dp,0.0_dp), &
            yulesimon_pmf(3.0_dp,2.5_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI yulesimon cdf',scifort_yulesimon_cdf_f64(3.0_dp,2.5_dp,0.0_dp), &
            yulesimon_cdf(3.0_dp,2.5_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI yulesimon ppf',scifort_yulesimon_ppf_f64(0.7_dp,2.5_dp,0.0_dp), &
            yulesimon_ppf(0.7_dp,2.5_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI zipf pmf',scifort_zipf_pmf_f64(3.0_dp,2.3_dp,0.0_dp), &
            zipf_pmf(3.0_dp,2.3_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI zipf cdf',scifort_zipf_cdf_f64(3.0_dp,2.3_dp,0.0_dp), &
            zipf_cdf(3.0_dp,2.3_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI zipf ppf',scifort_zipf_ppf_f64(0.7_dp,2.3_dp,0.0_dp), &
            zipf_ppf(0.7_dp,2.3_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI zipfian pmf',scifort_zipfian_pmf_f64(3.0_dp,1.2_dp,50.0_dp,0.0_dp), &
            zipfian_pmf(3.0_dp,1.2_dp,50.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI zipfian cdf',scifort_zipfian_cdf_f64(3.0_dp,1.2_dp,50.0_dp,0.0_dp), &
            zipfian_cdf(3.0_dp,1.2_dp,50.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI zipfian ppf',scifort_zipfian_ppf_f64(0.7_dp,1.2_dp,50.0_dp,0.0_dp), &
            zipfian_ppf(0.7_dp,1.2_dp,50.0_dp),0.0_dp,0.0_dp,failures)
    end subroutine test_rng_fit_c_api

    pure function centered(plus_value, minus_value, h) result(value)
        real(dp), intent(in) :: plus_value !! function value at theta+h
        real(dp), intent(in) :: minus_value !! function value at theta-h
        real(dp), intent(in) :: h !! positive centered-difference step
        real(dp) :: value
        value = (plus_value-minus_value)/(2.0_dp*h)
    end function centered
end program test_betanbinom_yulesimon_zipf_zipfian
