! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_randint_planck_dlaplace_logser
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_c_api
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state, rng_uniform
    use scifort_stats
    use test_checks, only : check_close, check_true
    use test_randint_planck_dlaplace_logser_reference
    implicit none

    integer :: failures

    failures = 0
    call test_reference(failures)
    call test_identities_endpoints(failures)
    call test_scores(failures)
    call test_rng_fit_c_api(failures)

    if (failures == 0) then
        print '(a)', 'test_randint_planck_dlaplace_logser: PASS'
    else
        print '(a,1x,i0)', 'test_randint_planck_dlaplace_logser: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        integer :: i

        do i = 1, size(reference_probs)
            call check_close('randint pmf',randint_pmf(ri_x(i),2.0_dp,9.0_dp),ri_pmf(i), &
                2.0e-15_dp,2.0e-15_dp,failures)
            call check_close('randint logpmf',randint_logpmf(ri_x(i),2.0_dp,9.0_dp),ri_logpmf(i), &
                2.0e-15_dp,2.0e-15_dp,failures)
            call check_close('randint cdf',randint_cdf(ri_x(i),2.0_dp,9.0_dp),ri_cdf(i), &
                2.0e-15_dp,2.0e-15_dp,failures)
            call check_close('randint sf',randint_sf(ri_x(i),2.0_dp,9.0_dp),ri_sf(i), &
                2.0e-15_dp,2.0e-15_dp,failures)
            call check_close('randint logcdf',randint_logcdf(ri_x(i),2.0_dp,9.0_dp),ri_logcdf(i), &
                2.0e-15_dp,2.0e-15_dp,failures)
            call check_close('randint logsf',randint_logsf(ri_x(i),2.0_dp,9.0_dp),ri_logsf(i), &
                2.0e-15_dp,2.0e-15_dp,failures)
            call check_close('randint ppf',randint_ppf(reference_probs(i),2.0_dp,9.0_dp),ri_ppf(i), &
                0.0_dp,0.0_dp,failures)
            call check_close('randint isf',randint_isf(reference_probs(i),2.0_dp,9.0_dp),ri_isf(i), &
                0.0_dp,0.0_dp,failures)

            call check_close('planck pmf',planck_pmf(pl_x(i),0.7_dp),pl_pmf(i), &
                3.0e-15_dp,3.0e-14_dp,failures)
            call check_close('planck logpmf',planck_logpmf(pl_x(i),0.7_dp),pl_logpmf(i), &
                3.0e-15_dp,3.0e-14_dp,failures)
            call check_close('planck cdf',planck_cdf(pl_x(i),0.7_dp),pl_cdf(i), &
                3.0e-15_dp,3.0e-14_dp,failures)
            call check_close('planck sf',planck_sf(pl_x(i),0.7_dp),pl_sf(i), &
                3.0e-15_dp,3.0e-14_dp,failures)
            call check_close('planck logcdf',planck_logcdf(pl_x(i),0.7_dp),pl_logcdf(i), &
                3.0e-15_dp,3.0e-14_dp,failures)
            call check_close('planck logsf',planck_logsf(pl_x(i),0.7_dp),pl_logsf(i), &
                3.0e-15_dp,3.0e-14_dp,failures)
            call check_close('planck ppf',planck_ppf(reference_probs(i),0.7_dp),pl_ppf(i), &
                0.0_dp,0.0_dp,failures)
            call check_close('planck isf',planck_isf(reference_probs(i),0.7_dp),pl_isf(i), &
                0.0_dp,0.0_dp,failures)

            call check_close('dlaplace pmf',dlaplace_pmf(dl_x(i),0.8_dp),dl_pmf(i), &
                3.0e-15_dp,3.0e-14_dp,failures)
            call check_close('dlaplace logpmf',dlaplace_logpmf(dl_x(i),0.8_dp),dl_logpmf(i), &
                3.0e-15_dp,3.0e-14_dp,failures)
            call check_close('dlaplace cdf',dlaplace_cdf(dl_x(i),0.8_dp),dl_cdf(i), &
                3.0e-15_dp,3.0e-14_dp,failures)
            call check_close('dlaplace sf',dlaplace_sf(dl_x(i),0.8_dp),dl_sf(i), &
                3.0e-15_dp,3.0e-14_dp,failures)
            call check_close('dlaplace logcdf',dlaplace_logcdf(dl_x(i),0.8_dp),dl_logcdf(i), &
                3.0e-15_dp,3.0e-14_dp,failures)
            call check_close('dlaplace logsf',dlaplace_logsf(dl_x(i),0.8_dp),dl_logsf(i), &
                3.0e-15_dp,3.0e-14_dp,failures)
            call check_close('dlaplace ppf',dlaplace_ppf(reference_probs(i),0.8_dp),dl_ppf(i), &
                0.0_dp,0.0_dp,failures)
            call check_close('dlaplace isf',dlaplace_isf(reference_probs(i),0.8_dp),dl_isf(i), &
                0.0_dp,0.0_dp,failures)

            call check_close('logser pmf',logser_pmf(ls_x(i),0.8_dp),ls_pmf(i), &
                8.0e-15_dp,8.0e-13_dp,failures)
            call check_close('logser logpmf',logser_logpmf(ls_x(i),0.8_dp),ls_logpmf(i), &
                8.0e-14_dp,8.0e-13_dp,failures)
            call check_close('logser cdf',logser_cdf(ls_x(i),0.8_dp),ls_cdf(i), &
                8.0e-15_dp,8.0e-13_dp,failures)
            call check_close('logser sf',logser_sf(ls_x(i),0.8_dp),ls_sf(i), &
                8.0e-15_dp,8.0e-13_dp,failures)
            call check_close('logser logcdf',logser_logcdf(ls_x(i),0.8_dp),ls_logcdf(i), &
                8.0e-14_dp,8.0e-13_dp,failures)
            call check_close('logser logsf',logser_logsf(ls_x(i),0.8_dp),ls_logsf(i), &
                8.0e-14_dp,8.0e-13_dp,failures)
            call check_close('logser ppf',logser_ppf(reference_probs(i),0.8_dp),ls_ppf(i), &
                0.0_dp,0.0_dp,failures)
            call check_close('logser isf',logser_isf(reference_probs(i),0.8_dp),ls_isf(i), &
                0.0_dp,0.0_dp,failures)
        end do

        call check_close('planck far isf',planck_isf(1.0e-12_dp,0.05_dp),pl_far_isf(1), &
            0.0_dp,0.0_dp,failures)
        call check_close('dlaplace far isf',dlaplace_isf(1.0e-12_dp,0.05_dp),dl_far_isf(1), &
            0.0_dp,0.0_dp,failures)
        call check_close('logser far isf',logser_isf(1.0e-10_dp,0.99_dp),ls_far_isf(1), &
            0.0_dp,0.0_dp,failures)
        call check_close('logser far sf',logser_sf(1000.0_dp,0.99_dp),ls_far_sf(1), &
            2.0e-18_dp,2.0e-11_dp,failures)
    end subroutine test_reference

    subroutine test_identities_endpoints(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: lambda = 0.7_dp
        real(dp) :: gp, total
        integer :: k

        total = 0.0_dp
        do k = 2, 8
            total = total + randint_pmf(k,2,9)
        end do
        call check_close('randint normalization',total,1.0_dp,2.0e-15_dp,2.0e-15_dp,failures)
        call check_close('randint ppf zero',randint_ppf(0.0_dp,2.0_dp,9.0_dp),1.0_dp, &
            0.0_dp,0.0_dp,failures)
        call check_close('randint ppf one',randint_ppf(1.0_dp,2.0_dp,9.0_dp),8.0_dp, &
            0.0_dp,0.0_dp,failures)

        gp = 1.0_dp - exp(-lambda)
        do k = 0, 8
            call check_close('planck-geometric pmf',planck_pmf(k,lambda),geometric_pmf(k+1,gp), &
                3.0e-15_dp,3.0e-14_dp,failures)
            call check_close('planck-geometric cdf',planck_cdf(k,lambda),geometric_cdf(k+1,gp), &
                3.0e-15_dp,3.0e-14_dp,failures)
        end do
        do k = 0, 8
            call check_close('dlaplace symmetry pmf',dlaplace_pmf(-k,0.8_dp),dlaplace_pmf(k,0.8_dp), &
                2.0e-15_dp,2.0e-15_dp,failures)
            call check_close('dlaplace reflected tails',dlaplace_cdf(-k-1,0.8_dp), &
                dlaplace_sf(k,0.8_dp),3.0e-15_dp,3.0e-14_dp,failures)
        end do
        call check_close('logser cdf+sf',logser_cdf(100.0_dp,0.99_dp)+logser_sf(100.0_dp,0.99_dp), &
            1.0_dp,3.0e-15_dp,3.0e-15_dp,failures)
        call check_close('logser mass-tail identity',sum([(logser_pmf(k,0.8_dp),k=1,100)]) + &
            logser_sf(100.0_dp,0.8_dp),1.0_dp,2.0e-14_dp,2.0e-14_dp,failures)

        call check_true('randint invalid bounds',ieee_is_nan(randint_pmf(2.0_dp,3.0_dp,3.0_dp)),failures)
        call check_true('randint invalid noninteger',ieee_is_nan(randint_pmf(2.0_dp,2.5_dp,9.0_dp)),failures)
        call check_true('planck invalid lambda',ieee_is_nan(planck_pmf(0.0_dp,0.0_dp)),failures)
        call check_true('dlaplace invalid a',ieee_is_nan(dlaplace_pmf(0.0_dp,0.0_dp)),failures)
        call check_true('logser invalid p',ieee_is_nan(logser_pmf(1.0_dp,1.0_dp)),failures)
    end subroutine test_identities_endpoints

    subroutine test_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp), parameter :: pd(5) = [0.0_dp,1.0_dp,1.0_dp,3.0_dp,5.0_dp]
        real(dp), parameter :: dd(5) = [-3.0_dp,-1.0_dp,0.0_dp,2.0_dp,4.0_dp]
        real(dp), parameter :: ld(5) = [1.0_dp,1.0_dp,2.0_dp,3.0_dp,7.0_dp]
        real(dp) :: s1(1)

        call check_true('randint score has zero size',size(randint_score([2,4,8],2,9)) == 0,failures)

        s1 = planck_score(pd,0.7_dp)
        call check_close('planck score lambda',s1(1),centered( &
            planck_loglikelihood(pd,0.7_dp+h),planck_loglikelihood(pd,0.7_dp-h),h), &
            5.0e-8_dp,5.0e-8_dp,failures)
        s1 = dlaplace_score(dd,0.8_dp)
        call check_close('dlaplace score a',s1(1),centered( &
            dlaplace_loglikelihood(dd,0.8_dp+h),dlaplace_loglikelihood(dd,0.8_dp-h),h), &
            5.0e-8_dp,5.0e-8_dp,failures)
        s1 = logser_score(ld,0.8_dp)
        call check_close('logser score p',s1(1),centered( &
            logser_loglikelihood(ld,0.8_dp+h),logser_loglikelihood(ld,0.8_dp-h),h), &
            5.0e-8_dp,5.0e-8_dp,failures)
    end subroutine test_scores

    subroutine test_rng_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: rid(4) = [2.0_dp,3.0_dp,5.0_dp,8.0_dp]
        real(dp), parameter :: pld(4) = [0.0_dp,1.0_dp,2.0_dp,4.0_dp]
        real(dp), parameter :: dld(4) = [-2.0_dp,-1.0_dp,0.0_dp,3.0_dp]
        real(dp), parameter :: lsd(4) = [1.0_dp,1.0_dp,2.0_dp,5.0_dp]
        real(dp) :: actual, expected
        type(rng_state) :: reference, state
        type(fit_result) :: result

        call rng_seed(state,24680)
        call rng_seed(reference,24680)
        actual = randint_rvs(state,2.0_dp,9.0_dp)
        expected = randint_ppf(rng_uniform(reference),2.0_dp,9.0_dp)
        call check_close('randint rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = planck_rvs(state,0.7_dp)
        expected = planck_ppf(rng_uniform(reference),0.7_dp)
        call check_close('planck rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = dlaplace_rvs(state,0.8_dp)
        expected = dlaplace_ppf(rng_uniform(reference),0.8_dp)
        call check_close('dlaplace rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = logser_rvs(state,0.8_dp)
        expected = logser_ppf(rng_uniform(reference),0.8_dp)
        call check_close('logser rvs',actual,expected,0.0_dp,0.0_dp,failures)

        call rng_seed(state,13579)
        call rng_seed(reference,13579)
        actual = logser_rvs(state,1.0_dp)
        call check_true('invalid logser rvs NaN',ieee_is_nan(actual),failures)
        call check_close('invalid rvs preserves state',rng_uniform(state),rng_uniform(reference), &
            0.0_dp,0.0_dp,failures)

        call randint_fit(rid,[2.0_dp,9.0_dp,0.0_dp],[2.0_dp,9.0_dp,0.0_dp],result)
        call check_true('randint fixed fit',result%success,failures)
        call planck_fit(pld,[0.7_dp,0.0_dp],[0.7_dp,0.0_dp],result)
        call check_true('planck fixed fit',result%success,failures)
        call dlaplace_fit(dld,[0.8_dp,0.0_dp],[0.8_dp,0.0_dp],result)
        call check_true('dlaplace fixed fit',result%success,failures)
        call logser_fit(lsd,[0.8_dp,0.0_dp],[0.8_dp,0.0_dp],result)
        call check_true('logser fixed fit',result%success,failures)

        call check_close('C ABI randint pmf',scifort_randint_pmf_f64(4.0_dp,2.0_dp,9.0_dp,0.0_dp), &
            randint_pmf(4.0_dp,2.0_dp,9.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI randint cdf',scifort_randint_cdf_f64(4.0_dp,2.0_dp,9.0_dp,0.0_dp), &
            randint_cdf(4.0_dp,2.0_dp,9.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI randint ppf',scifort_randint_ppf_f64(0.4_dp,2.0_dp,9.0_dp,0.0_dp), &
            randint_ppf(0.4_dp,2.0_dp,9.0_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI planck pmf',scifort_planck_pmf_f64(3.0_dp,0.7_dp,0.0_dp), &
            planck_pmf(3.0_dp,0.7_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI planck cdf',scifort_planck_cdf_f64(3.0_dp,0.7_dp,0.0_dp), &
            planck_cdf(3.0_dp,0.7_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI planck ppf',scifort_planck_ppf_f64(0.4_dp,0.7_dp,0.0_dp), &
            planck_ppf(0.4_dp,0.7_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI dlaplace pmf',scifort_dlaplace_pmf_f64(-2.0_dp,0.8_dp,0.0_dp), &
            dlaplace_pmf(-2.0_dp,0.8_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI dlaplace cdf',scifort_dlaplace_cdf_f64(-2.0_dp,0.8_dp,0.0_dp), &
            dlaplace_cdf(-2.0_dp,0.8_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI dlaplace ppf',scifort_dlaplace_ppf_f64(0.4_dp,0.8_dp,0.0_dp), &
            dlaplace_ppf(0.4_dp,0.8_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI logser pmf',scifort_logser_pmf_f64(4.0_dp,0.8_dp,0.0_dp), &
            logser_pmf(4.0_dp,0.8_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI logser cdf',scifort_logser_cdf_f64(4.0_dp,0.8_dp,0.0_dp), &
            logser_cdf(4.0_dp,0.8_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI logser ppf',scifort_logser_ppf_f64(0.4_dp,0.8_dp,0.0_dp), &
            logser_ppf(0.4_dp,0.8_dp),0.0_dp,0.0_dp,failures)
    end subroutine test_rng_fit_c_api

    pure function centered(fp, fm, h) result(value)
        real(dp), intent(in) :: fp !! function value at positive perturbation
        real(dp), intent(in) :: fm !! function value at negative perturbation
        real(dp), intent(in) :: h !! perturbation magnitude
        real(dp) :: value
        value = (fp - fm) / (2.0_dp*h)
    end function centered

end program test_randint_planck_dlaplace_logser
