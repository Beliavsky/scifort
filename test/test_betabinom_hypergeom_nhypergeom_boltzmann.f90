! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_betabinom_hypergeom_nhypergeom_boltzmann
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_c_api
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state, rng_uniform
    use scifort_stats
    use test_betabinom_hypergeom_nhypergeom_boltzmann_reference
    use test_checks, only : check_close, check_true
    implicit none

    integer :: failures

    failures = 0
    call test_reference(failures)
    call test_identities_endpoints(failures)
    call test_scores(failures)
    call test_rng_fit_c_api(failures)

    if (failures == 0) then
        print '(a)', 'test_betabinom_hypergeom_nhypergeom_boltzmann: PASS'
    else
        print '(a,1x,i0)', 'test_betabinom_hypergeom_nhypergeom_boltzmann: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        integer :: i

        do i = 1, size(reference_probs)
            call check_close('betabinom pmf',betabinom_pmf(bb_x(i),12.0_dp,2.5_dp,4.0_dp),bb_pmf(i), &
                2.0e-14_dp,2.0e-13_dp,failures)
            call check_close('betabinom logpmf',betabinom_logpmf(bb_x(i),12.0_dp,2.5_dp,4.0_dp),bb_logpmf(i), &
                2.0e-13_dp,2.0e-13_dp,failures)
            call check_close('betabinom cdf',betabinom_cdf(bb_x(i),12.0_dp,2.5_dp,4.0_dp),bb_cdf(i), &
                2.0e-14_dp,2.0e-13_dp,failures)
            call check_close('betabinom sf',betabinom_sf(bb_x(i),12.0_dp,2.5_dp,4.0_dp),bb_sf(i), &
                2.0e-14_dp,2.0e-13_dp,failures)
            call check_close('betabinom logcdf',betabinom_logcdf(bb_x(i),12.0_dp,2.5_dp,4.0_dp),bb_logcdf(i), &
                2.0e-13_dp,2.0e-13_dp,failures)
            call check_close('betabinom logsf',betabinom_logsf(bb_x(i),12.0_dp,2.5_dp,4.0_dp),bb_logsf(i), &
                2.0e-13_dp,2.0e-13_dp,failures)
            call check_close('betabinom ppf',betabinom_ppf(reference_probs(i),12.0_dp,2.5_dp,4.0_dp),bb_ppf(i), &
                0.0_dp,0.0_dp,failures)
            call check_close('betabinom isf',betabinom_isf(reference_probs(i),12.0_dp,2.5_dp,4.0_dp),bb_isf(i), &
                0.0_dp,0.0_dp,failures)

            call check_close('hypergeom pmf',hypergeom_pmf(hg_x(i),40.0_dp,15.0_dp,12.0_dp),hg_pmf(i), &
                2.0e-14_dp,3.0e-13_dp,failures)
            call check_close('hypergeom logpmf',hypergeom_logpmf(hg_x(i),40.0_dp,15.0_dp,12.0_dp),hg_logpmf(i), &
                2.0e-13_dp,3.0e-13_dp,failures)
            call check_close('hypergeom cdf',hypergeom_cdf(hg_x(i),40.0_dp,15.0_dp,12.0_dp),hg_cdf(i), &
                2.0e-14_dp,3.0e-13_dp,failures)
            call check_close('hypergeom sf',hypergeom_sf(hg_x(i),40.0_dp,15.0_dp,12.0_dp),hg_sf(i), &
                2.0e-14_dp,3.0e-13_dp,failures)
            call check_close('hypergeom logcdf',hypergeom_logcdf(hg_x(i),40.0_dp,15.0_dp,12.0_dp),hg_logcdf(i), &
                2.0e-13_dp,3.0e-13_dp,failures)
            call check_close('hypergeom logsf',hypergeom_logsf(hg_x(i),40.0_dp,15.0_dp,12.0_dp),hg_logsf(i), &
                2.0e-13_dp,3.0e-13_dp,failures)
            call check_close('hypergeom ppf',hypergeom_ppf(reference_probs(i),40.0_dp,15.0_dp,12.0_dp),hg_ppf(i), &
                0.0_dp,0.0_dp,failures)
            call check_close('hypergeom isf',hypergeom_isf(reference_probs(i),40.0_dp,15.0_dp,12.0_dp),hg_isf(i), &
                0.0_dp,0.0_dp,failures)

            call check_close('nhypergeom pmf',nhypergeom_pmf(nh_x(i),40.0_dp,15.0_dp,8.0_dp),nh_pmf(i), &
                2.0e-14_dp,3.0e-13_dp,failures)
            call check_close('nhypergeom logpmf',nhypergeom_logpmf(nh_x(i),40.0_dp,15.0_dp,8.0_dp),nh_logpmf(i), &
                2.0e-13_dp,3.0e-13_dp,failures)
            call check_close('nhypergeom cdf',nhypergeom_cdf(nh_x(i),40.0_dp,15.0_dp,8.0_dp),nh_cdf(i), &
                2.0e-14_dp,3.0e-13_dp,failures)
            call check_close('nhypergeom sf',nhypergeom_sf(nh_x(i),40.0_dp,15.0_dp,8.0_dp),nh_sf(i), &
                2.0e-14_dp,3.0e-13_dp,failures)
            call check_close('nhypergeom logcdf',nhypergeom_logcdf(nh_x(i),40.0_dp,15.0_dp,8.0_dp),nh_logcdf(i), &
                2.0e-13_dp,3.0e-13_dp,failures)
            call check_close('nhypergeom logsf',nhypergeom_logsf(nh_x(i),40.0_dp,15.0_dp,8.0_dp),nh_logsf(i), &
                2.0e-13_dp,3.0e-13_dp,failures)
            call check_close('nhypergeom ppf',nhypergeom_ppf(reference_probs(i),40.0_dp,15.0_dp,8.0_dp),nh_ppf(i), &
                0.0_dp,0.0_dp,failures)
            call check_close('nhypergeom isf',nhypergeom_isf(reference_probs(i),40.0_dp,15.0_dp,8.0_dp),nh_isf(i), &
                0.0_dp,0.0_dp,failures)

            call check_close('boltzmann pmf',boltzmann_pmf(bo_x(i),0.35_dp,12.0_dp),bo_pmf(i), &
                3.0e-15_dp,3.0e-14_dp,failures)
            call check_close('boltzmann logpmf',boltzmann_logpmf(bo_x(i),0.35_dp,12.0_dp),bo_logpmf(i), &
                3.0e-14_dp,3.0e-14_dp,failures)
            call check_close('boltzmann cdf',boltzmann_cdf(bo_x(i),0.35_dp,12.0_dp),bo_cdf(i), &
                3.0e-15_dp,3.0e-14_dp,failures)
            call check_close('boltzmann sf',boltzmann_sf(bo_x(i),0.35_dp,12.0_dp),bo_sf(i), &
                3.0e-15_dp,3.0e-14_dp,failures)
            call check_close('boltzmann logcdf',boltzmann_logcdf(bo_x(i),0.35_dp,12.0_dp),bo_logcdf(i), &
                3.0e-14_dp,3.0e-14_dp,failures)
            call check_close('boltzmann logsf',boltzmann_logsf(bo_x(i),0.35_dp,12.0_dp),bo_logsf(i), &
                3.0e-14_dp,3.0e-14_dp,failures)
            call check_close('boltzmann ppf',boltzmann_ppf(reference_probs(i),0.35_dp,12.0_dp),bo_ppf(i), &
                0.0_dp,0.0_dp,failures)
            call check_close('boltzmann isf',boltzmann_isf(reference_probs(i),0.35_dp,12.0_dp),bo_isf(i), &
                0.0_dp,0.0_dp,failures)
        end do

        call check_close('betabinom stress sf',betabinom_sf(150.0_dp,200.0_dp,0.2_dp,3.5_dp),bb_stress_sf(1), &
            2.0e-14_dp,2.0e-10_dp,failures)
        call check_close('betabinom stress isf',betabinom_isf(1.0e-10_dp,200.0_dp,0.2_dp,3.5_dp), &
            bb_stress_isf(1),0.0_dp,0.0_dp,failures)
        call check_close('hypergeom stress sf',hypergeom_sf(140.0_dp,10000.0_dp,2500.0_dp,400.0_dp),hg_stress_sf(1), &
            2.0e-16_dp,5.0e-10_dp,failures)
        call check_close('hypergeom stress ppf',hypergeom_ppf(1.0_dp-1.0e-10_dp,10000.0_dp,2500.0_dp,400.0_dp), &
            hg_stress_ppf(1),0.0_dp,0.0_dp,failures)
        call check_close('nhypergeom stress sf',nhypergeom_sf(60.0_dp,500.0_dp,120.0_dp,100.0_dp),nh_stress_sf(1), &
            5.0e-19_dp,2.0e-12_dp,failures)
        call check_close('nhypergeom stress isf',nhypergeom_isf(1.0e-10_dp,500.0_dp,120.0_dp,100.0_dp), &
            nh_stress_isf(1),0.0_dp,0.0_dp,failures)
        call check_close('boltzmann stress ppf',boltzmann_ppf(0.999999_dp,1.0e-8_dp,100000.0_dp),bo_stress_ppf(1), &
            0.0_dp,0.0_dp,failures)
        call check_close('boltzmann stress sf',boltzmann_sf(99990.0_dp,1.0e-8_dp,100000.0_dp),bo_stress_sf(1), &
            5.0e-18_dp,2.0e-12_dp,failures)
    end subroutine test_reference

    subroutine test_identities_endpoints(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: a = 2.3_dp
        real(dp), parameter :: b = 4.1_dp
        real(dp), parameter :: lambda = 0.35_dp
        real(dp) :: norm, rhs
        integer :: k

        call check_close('betabinom Bernoulli zero',betabinom_pmf(0,1,a,b),bernoulli_pmf(0,a/(a+b)), &
            3.0e-15_dp,3.0e-14_dp,failures)
        call check_close('betabinom Bernoulli one',betabinom_pmf(1,1,a,b),bernoulli_pmf(1,a/(a+b)), &
            3.0e-15_dp,3.0e-14_dp,failures)
        do k = 0, 12
            call check_close('betabinom symmetric',betabinom_pmf(k,12,2.5_dp,2.5_dp), &
                betabinom_pmf(12-k,12,2.5_dp,2.5_dp),2.0e-14_dp,2.0e-13_dp,failures)
        end do
        norm = sum([(betabinom_pmf(k,12,2.5_dp,4.0_dp),k=0,12)])
        call check_close('betabinom normalization',norm,1.0_dp,3.0e-14_dp,3.0e-14_dp,failures)

        do k = 0, 12
            call check_close('hypergeom complement',hypergeom_pmf(k,40,15,12), &
                hypergeom_pmf(12-k,40,25,12),2.0e-14_dp,3.0e-13_dp,failures)
        end do
        norm = sum([(hypergeom_pmf(k,40,15,12),k=0,12)])
        call check_close('hypergeom normalization',norm,1.0_dp,3.0e-14_dp,3.0e-14_dp,failures)

        do k = 0, 12
            rhs = hypergeom_pmf(k,40,15,k+8-1) * real(40-15-(8-1),dp) / real(40-(k+8-1),dp)
            call check_close('nhypergeom-hypergeom identity',nhypergeom_pmf(k,40,15,8),rhs, &
                4.0e-14_dp,4.0e-13_dp,failures)
        end do
        norm = sum([(nhypergeom_pmf(k,40,15,8),k=0,15)])
        call check_close('nhypergeom normalization',norm,1.0_dp,4.0e-14_dp,4.0e-14_dp,failures)

        do k = 0, 11
            call check_close('boltzmann-planck truncation',boltzmann_pmf(k,lambda,12), &
                planck_pmf(k,lambda)/planck_cdf(11,lambda),3.0e-15_dp,3.0e-14_dp,failures)
        end do
        norm = sum([(boltzmann_pmf(k,lambda,12),k=0,11)])
        call check_close('boltzmann normalization',norm,1.0_dp,4.0e-15_dp,4.0e-15_dp,failures)

        call check_close('betabinom ppf zero',betabinom_ppf(0.0_dp,12.0_dp,2.5_dp,4.0_dp),-1.0_dp,0.0_dp,0.0_dp,failures)
        call check_close('betabinom ppf one',betabinom_ppf(1.0_dp,12.0_dp,2.5_dp,4.0_dp),12.0_dp,0.0_dp,0.0_dp,failures)
        call check_close('hypergeom ppf zero',hypergeom_ppf(0.0_dp,40.0_dp,15.0_dp,12.0_dp),-1.0_dp,0.0_dp,0.0_dp,failures)
        call check_close('hypergeom ppf one',hypergeom_ppf(1.0_dp,40.0_dp,15.0_dp,12.0_dp),12.0_dp,0.0_dp,0.0_dp,failures)
        call check_close('nhypergeom ppf zero',nhypergeom_ppf(0.0_dp,40.0_dp,15.0_dp,8.0_dp),-1.0_dp,0.0_dp,0.0_dp,failures)
        call check_close('nhypergeom ppf one',nhypergeom_ppf(1.0_dp,40.0_dp,15.0_dp,8.0_dp),15.0_dp,0.0_dp,0.0_dp,failures)
        call check_close('boltzmann ppf zero',boltzmann_ppf(0.0_dp,0.35_dp,12.0_dp),-1.0_dp,0.0_dp,0.0_dp,failures)
        call check_close('boltzmann ppf one',boltzmann_ppf(1.0_dp,0.35_dp,12.0_dp),11.0_dp,0.0_dp,0.0_dp,failures)

        call check_true('betabinom invalid a',ieee_is_nan(betabinom_pmf(2,12,0.0_dp,4.0_dp)),failures)
        call check_true('hypergeom invalid n',ieee_is_nan(hypergeom_pmf(2,40,41,12)),failures)
        call check_true('nhypergeom invalid r',ieee_is_nan(nhypergeom_pmf(2,40,15,26)),failures)
        call check_true('boltzmann invalid lambda',ieee_is_nan(boltzmann_pmf(2,0.0_dp,12)),failures)

        call check_close('betabinom shifted pmf',betabinom_pmf(4.5_dp,12.0_dp,2.5_dp,4.0_dp,0.5_dp), &
            betabinom_pmf(4.0_dp,12.0_dp,2.5_dp,4.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('hypergeom shifted cdf',hypergeom_cdf(4.5_dp,40.0_dp,15.0_dp,12.0_dp,0.5_dp), &
            hypergeom_cdf(4.0_dp,40.0_dp,15.0_dp,12.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('nhypergeom shifted sf',nhypergeom_sf(4.5_dp,40.0_dp,15.0_dp,8.0_dp,0.5_dp), &
            nhypergeom_sf(4.0_dp,40.0_dp,15.0_dp,8.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('boltzmann shifted ppf',boltzmann_ppf(0.4_dp,0.35_dp,12.0_dp,0.5_dp), &
            boltzmann_ppf(0.4_dp,0.35_dp,12.0_dp)+0.5_dp,0.0_dp,0.0_dp,failures)
    end subroutine test_identities_endpoints

    subroutine test_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp), parameter :: bd(5) = [0.0_dp,1.0_dp,3.0_dp,5.0_dp,8.0_dp]
        real(dp), parameter :: zd(5) = [0.0_dp,1.0_dp,1.0_dp,3.0_dp,7.0_dp]
        real(dp) :: s1(1), s2(2)

        s2 = betabinom_score(bd,10.0_dp,2.3_dp,4.1_dp)
        call check_close('betabinom score a',s2(1),centered( &
            betabinom_loglikelihood(bd,10.0_dp,2.3_dp+h,4.1_dp), &
            betabinom_loglikelihood(bd,10.0_dp,2.3_dp-h,4.1_dp),h),5.0e-8_dp,5.0e-8_dp,failures)
        call check_close('betabinom score b',s2(2),centered( &
            betabinom_loglikelihood(bd,10.0_dp,2.3_dp,4.1_dp+h), &
            betabinom_loglikelihood(bd,10.0_dp,2.3_dp,4.1_dp-h),h),5.0e-8_dp,5.0e-8_dp,failures)
        call check_true('hypergeom score zero size',size(hypergeom_score([1,3,5],40,15,12)) == 0,failures)
        call check_true('nhypergeom score zero size',size(nhypergeom_score([0,2,4],40,15,8)) == 0,failures)
        s1 = boltzmann_score(zd,0.35_dp,12.0_dp)
        call check_close('boltzmann score lambda',s1(1),centered( &
            boltzmann_loglikelihood(zd,0.35_dp+h,12.0_dp), &
            boltzmann_loglikelihood(zd,0.35_dp-h,12.0_dp),h),5.0e-8_dp,5.0e-8_dp,failures)
    end subroutine test_scores

    subroutine test_rng_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: bbd(4) = [0.0_dp,1.0_dp,3.0_dp,5.0_dp]
        real(dp), parameter :: hgd(4) = [1.0_dp,3.0_dp,5.0_dp,7.0_dp]
        real(dp), parameter :: nhd(4) = [0.0_dp,2.0_dp,4.0_dp,7.0_dp]
        real(dp), parameter :: bod(4) = [0.0_dp,1.0_dp,2.0_dp,4.0_dp]
        real(dp) :: actual, expected
        type(rng_state) :: reference, state
        type(fit_result) :: result

        call rng_seed(state,24680)
        call rng_seed(reference,24680)
        actual = betabinom_rvs(state,12.0_dp,2.5_dp,4.0_dp)
        expected = betabinom_ppf(rng_uniform(reference),12.0_dp,2.5_dp,4.0_dp)
        call check_close('betabinom rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = hypergeom_rvs(state,40.0_dp,15.0_dp,12.0_dp)
        expected = hypergeom_ppf(rng_uniform(reference),40.0_dp,15.0_dp,12.0_dp)
        call check_close('hypergeom rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = nhypergeom_rvs(state,40.0_dp,15.0_dp,8.0_dp)
        expected = nhypergeom_ppf(rng_uniform(reference),40.0_dp,15.0_dp,8.0_dp)
        call check_close('nhypergeom rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = boltzmann_rvs(state,0.35_dp,12.0_dp)
        expected = boltzmann_ppf(rng_uniform(reference),0.35_dp,12.0_dp)
        call check_close('boltzmann rvs',actual,expected,0.0_dp,0.0_dp,failures)

        call rng_seed(state,13579)
        call rng_seed(reference,13579)
        actual = betabinom_rvs(state,12.0_dp,0.0_dp,4.0_dp)
        call check_true('invalid betabinom rvs NaN',ieee_is_nan(actual),failures)
        call check_close('invalid rvs preserves state',rng_uniform(state),rng_uniform(reference), &
            0.0_dp,0.0_dp,failures)

        call betabinom_fit(bbd,[12.0_dp,2.5_dp,4.0_dp,0.0_dp],[12.0_dp,2.5_dp,4.0_dp,0.0_dp],result)
        call check_true('betabinom fixed fit',result%success,failures)
        call hypergeom_fit(hgd,[40.0_dp,15.0_dp,12.0_dp,0.0_dp],[40.0_dp,15.0_dp,12.0_dp,0.0_dp],result)
        call check_true('hypergeom fixed fit',result%success,failures)
        call nhypergeom_fit(nhd,[40.0_dp,15.0_dp,8.0_dp,0.0_dp],[40.0_dp,15.0_dp,8.0_dp,0.0_dp],result)
        call check_true('nhypergeom fixed fit',result%success,failures)
        call boltzmann_fit(bod,[0.35_dp,12.0_dp,0.0_dp],[0.35_dp,12.0_dp,0.0_dp],result)
        call check_true('boltzmann fixed fit',result%success,failures)

        call check_close('C ABI betabinom pmf',scifort_betabinom_pmf_f64(4.0_dp,12.0_dp,2.5_dp,4.0_dp,0.0_dp), &
            betabinom_pmf(4.0_dp,12.0_dp,2.5_dp,4.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI betabinom cdf',scifort_betabinom_cdf_f64(4.0_dp,12.0_dp,2.5_dp,4.0_dp,0.0_dp), &
            betabinom_cdf(4.0_dp,12.0_dp,2.5_dp,4.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI betabinom ppf',scifort_betabinom_ppf_f64(0.4_dp,12.0_dp,2.5_dp,4.0_dp,0.0_dp), &
            betabinom_ppf(0.4_dp,12.0_dp,2.5_dp,4.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI hypergeom pmf',scifort_hypergeom_pmf_f64(4.0_dp,40.0_dp,15.0_dp,12.0_dp,0.0_dp), &
            hypergeom_pmf(4.0_dp,40.0_dp,15.0_dp,12.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI hypergeom cdf',scifort_hypergeom_cdf_f64(4.0_dp,40.0_dp,15.0_dp,12.0_dp,0.0_dp), &
            hypergeom_cdf(4.0_dp,40.0_dp,15.0_dp,12.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI hypergeom ppf',scifort_hypergeom_ppf_f64(0.4_dp,40.0_dp,15.0_dp,12.0_dp,0.0_dp), &
            hypergeom_ppf(0.4_dp,40.0_dp,15.0_dp,12.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI nhypergeom pmf',scifort_nhypergeom_pmf_f64(4.0_dp,40.0_dp,15.0_dp,8.0_dp,0.0_dp), &
            nhypergeom_pmf(4.0_dp,40.0_dp,15.0_dp,8.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI nhypergeom cdf',scifort_nhypergeom_cdf_f64(4.0_dp,40.0_dp,15.0_dp,8.0_dp,0.0_dp), &
            nhypergeom_cdf(4.0_dp,40.0_dp,15.0_dp,8.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI nhypergeom ppf',scifort_nhypergeom_ppf_f64(0.4_dp,40.0_dp,15.0_dp,8.0_dp,0.0_dp), &
            nhypergeom_ppf(0.4_dp,40.0_dp,15.0_dp,8.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI boltzmann pmf',scifort_boltzmann_pmf_f64(4.0_dp,0.35_dp,12.0_dp,0.0_dp), &
            boltzmann_pmf(4.0_dp,0.35_dp,12.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI boltzmann cdf',scifort_boltzmann_cdf_f64(4.0_dp,0.35_dp,12.0_dp,0.0_dp), &
            boltzmann_cdf(4.0_dp,0.35_dp,12.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI boltzmann ppf',scifort_boltzmann_ppf_f64(0.4_dp,0.35_dp,12.0_dp,0.0_dp), &
            boltzmann_ppf(0.4_dp,0.35_dp,12.0_dp),0.0_dp,0.0_dp,failures)
    end subroutine test_rng_fit_c_api

    pure function centered(fp, fm, h) result(value)
        real(dp), intent(in) :: fp !! function value at positive perturbation
        real(dp), intent(in) :: fm !! function value at negative perturbation
        real(dp), intent(in) :: h !! perturbation magnitude
        real(dp) :: value
        value = (fp - fm) / (2.0_dp*h)
    end function centered

end program test_betabinom_hypergeom_nhypergeom_boltzmann
