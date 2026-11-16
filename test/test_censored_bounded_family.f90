! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_bounded_family
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        uniform_censored_loglikelihood, uniform_censored_nnlf, uniform_fit, &
        uniform_logpdf, uniform_logcdf, uniform_logsf, uniform_ppf, uniform_loglikelihood, &
        powerlaw_censored_loglikelihood, powerlaw_censored_nnlf, powerlaw_fit, &
        powerlaw_logpdf, powerlaw_pdf, powerlaw_logcdf, powerlaw_logsf, powerlaw_ppf, &
        powerlaw_loglikelihood, arcsine_censored_loglikelihood, arcsine_censored_nnlf, arcsine_fit, &
        arcsine_logpdf, arcsine_logcdf, arcsine_logsf, arcsine_ppf, arcsine_loglikelihood, &
        beta_censored_loglikelihood, fit_status_invalid_input, fit_status_no_finite_objective, &
        fit_status_max_iter
    use scifort_math, only : quiet_nan, positive_infinity
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
    implicit none
    type(censored_data) :: data, single, empty, other
    type(fit_result) :: fit, old
    real(dp) :: inf, nan, actual, expected, x, a, p, previous, small
    real(dp) :: reference(3)
    integer :: family, i, j, checks, failures
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    small = nearest(0.0_dp,1.0_dp)
    reference = [-10.596634733096073_dp,-10.086476648213639_dp,-10.496569053424526_dp]
    data = make_censored_data(uncensored=[0.1_dp,0.3_dp,0.6_dp],left=[0.05_dp], &
        right=[0.8_dp],interval=reshape([0.2_dp,0.7_dp,0.25_dp,0.75_dp],[2,2]))
    do family = 1, 3
        actual = likelihood(data,family,0.6_dp)
        call close(actual,reference(family),6.0e-14_dp,'mixed mpmath reference')
        select case (family)
        case (1)
            expected = uniform_censored_nnlf(data)
        case (2)
            expected = powerlaw_censored_nnlf(data,0.6_dp)
        case (3)
            expected = arcsine_censored_nnlf(data)
        end select
        call close(expected,-actual,0.0_dp,'negative likelihood')
        do j = 1, 3
            a = 0.3_dp
            if (j == 2) a = 1.0_dp
            if (j == 3) a = 3.0_dp
            previous = -inf
            do i = 1, 12
                x = real(i,dp)/13.0_dp
                single = make_censored_data(uncensored=[x])
                select case (family)
                case (1)
                    expected = uniform_logpdf(x)
                case (2)
                    expected = powerlaw_logpdf(x,a)
                case (3)
                    expected = arcsine_logpdf(x)
                end select
                call close(likelihood(single,family,a),expected,3.0e-14_dp,'scalar density')
                single = make_censored_data(left=[x])
                select case (family)
                case (1)
                    expected = uniform_logcdf(x)
                case (2)
                    expected = powerlaw_logcdf(x,a)
                case (3)
                    expected = arcsine_logcdf(x)
                end select
                actual = likelihood(single,family,a)
                call close(actual,expected,3.0e-14_dp,'scalar CDF')
                call check(actual >= previous .and. actual <= 0.0_dp,'CDF monotonicity and range')
                previous = actual
                single = make_censored_data(right=[x])
                select case (family)
                case (1)
                    expected = uniform_logsf(x)
                case (2)
                    expected = powerlaw_logsf(x,a)
                case (3)
                    expected = arcsine_logsf(x)
                end select
                call close(likelihood(single,family,a),expected,3.0e-14_dp,'scalar SF')
            end do
            do i = 1, 9
                p = real(i,dp)/10.0_dp
                select case (family)
                case (1)
                    x = uniform_ppf(p)
                case (2)
                    x = powerlaw_ppf(p,a)
                case (3)
                    x = arcsine_ppf(p)
                end select
                single = make_censored_data(left=[x])
                call close(exp(likelihood(single,family,a)),p,3.0e-14_dp,'CDF/PPF round trip')
            end do
        end do
        single = make_censored_data(uncensored=[0.25_dp,0.5_dp,0.75_dp])
        select case (family)
        case (1)
            expected = uniform_loglikelihood(single%uncensored)
            call uniform_fit(single,[0.0_dp,1.0_dp],[0.0_dp,1.0_dp],fit)
            call uniform_fit(single%uncensored,[0.0_dp,1.0_dp],[0.0_dp,1.0_dp],old)
        case (2)
            expected = powerlaw_loglikelihood(single%uncensored,0.6_dp)
            call powerlaw_fit(single,[0.6_dp,0.0_dp,1.0_dp],[0.6_dp,0.0_dp,1.0_dp],fit)
            call powerlaw_fit(single%uncensored,[0.6_dp,0.0_dp,1.0_dp],[0.6_dp,0.0_dp,1.0_dp],old)
        case (3)
            expected = arcsine_loglikelihood(single%uncensored)
            call arcsine_fit(single,[0.0_dp,1.0_dp],[0.0_dp,1.0_dp],fit)
            call arcsine_fit(single%uncensored,[0.0_dp,1.0_dp],[0.0_dp,1.0_dp],old)
        end select
        call check(fit%success .and. old%success,'fixed fits converge')
        call close(fit%nllf,-expected,3.0e-14_dp,'array likelihood agreement')
        call close(fit%nllf,old%nllf,3.0e-14_dp,'array fit unchanged')
        call close(likelihood(empty,family,1.0_dp),0.0_dp,0.0_dp,'valid empty sample')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,loc=nan)),'NaN location')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,loc=inf)),'infinite location')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,scale=0.0_dp)),'zero width')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,scale=-1.0_dp)),'negative width')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,scale=nan)),'NaN width')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,scale=inf)),'infinite width')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,loc=huge(1.0_dp), &
            scale=huge(1.0_dp))),'overflowing upper edge rejected')
        single = make_censored_data(uncensored=[nan])
        call check(ieee_is_nan(likelihood(single,family,1.0_dp)),'invalid data')
        single = make_censored_data(uncensored=[-0.1_dp])
        call check(likelihood(single,family,1.0_dp) == -inf,'density below support')
        single = make_censored_data(uncensored=[1.1_dp])
        call check(likelihood(single,family,1.0_dp) == -inf,'density above support')
        single = make_censored_data(left=[0.0_dp])
        call check(likelihood(single,family,1.0_dp) == -inf,'lower CDF endpoint')
        single = make_censored_data(right=[1.0_dp])
        call check(likelihood(single,family,1.0_dp) == -inf,'upper SF endpoint')
        single = make_censored_data(left=[1.1_dp],right=[-0.1_dp])
        call close(likelihood(single,family,1.0_dp),0.0_dp,0.0_dp,'unit outside-support tails')
        single = make_censored_data(interval=reshape([-inf,0.5_dp],[1,2]))
        other = make_censored_data(left=[0.5_dp])
        call close(likelihood(single,family,0.6_dp),likelihood(other,family,0.6_dp), &
            0.0_dp,'infinite lower endpoint')
        single = make_censored_data(interval=reshape([0.5_dp,inf],[1,2]))
        other = make_censored_data(right=[0.5_dp])
        call close(likelihood(single,family,0.6_dp),likelihood(other,family,0.6_dp), &
            0.0_dp,'infinite upper endpoint')
        single = make_censored_data(interval=reshape([-1.0_dp,2.0_dp],[1,2]))
        call close(likelihood(single,family,0.6_dp),0.0_dp,0.0_dp,'support-clipped full interval')
        single = make_censored_data(interval=reshape([-1.0_dp,0.0_dp],[1,2]))
        call check(likelihood(single,family,0.6_dp) == -inf,'impossible interval')
        single = make_censored_data(interval=reshape([0.5_dp,0.5_dp],[1,2]))
        other = make_censored_data(uncensored=[0.5_dp])
        call close(likelihood(single,family,0.6_dp),likelihood(other,family,0.6_dp), &
            0.0_dp,'equal endpoints exact')
        single = make_censored_data(interval=reshape([0.5_dp,nearest(0.5_dp,1.0_dp)],[1,2]))
        call close(likelihood(single,family,0.6_dp), &
            likelihood(other,family,0.6_dp)+log(epsilon(1.0_dp)/2.0_dp), &
            4.0e-14_dp,'adjacent interval')
    end do
    call check(ieee_is_nan(powerlaw_censored_loglikelihood(empty,0.0_dp)),'zero shape')
    call check(ieee_is_nan(powerlaw_censored_loglikelihood(empty,nan)),'NaN shape')
    call check(ieee_is_nan(powerlaw_censored_loglikelihood(empty,inf)),'infinite shape')
    single = make_censored_data(uncensored=[0.0_dp])
    call check(powerlaw_censored_loglikelihood(single,0.5_dp) == inf,'singular powerlaw endpoint')
    call check(powerlaw_pdf(0.0_dp,0.5_dp) == inf,'scalar powerlaw endpoint density corrected')
    call check(powerlaw_logpdf(0.0_dp,0.5_dp) == inf,'scalar endpoint log density corrected')
    call close(powerlaw_censored_loglikelihood(single,1.0_dp),0.0_dp,0.0_dp,'unit lower density')
    call check(powerlaw_censored_loglikelihood(single,2.0_dp) == -inf,'zero lower density')
    call check(arcsine_censored_loglikelihood(single) == inf,'arcsine lower singularity')
    single = make_censored_data(uncensored=[1.0_dp])
    call close(powerlaw_censored_loglikelihood(single,2.0_dp),log(2.0_dp),0.0_dp,'upper density')
    call close(uniform_censored_loglikelihood(single),0.0_dp,0.0_dp,'uniform included endpoint')
    call check(arcsine_censored_loglikelihood(single) == inf,'arcsine upper singularity')
    single = make_censored_data(uncensored=[0.25_dp,0.5_dp,1.0_dp])
    call powerlaw_fit(single,[0.1_dp,0.0_dp,1.0_dp],[10.0_dp,0.0_dp,1.0_dp],fit)
    call check(fit%success,'analytic shape MLE converges')
    call close(fit%params(1),1.0_dp/log(2.0_dp),4.0e-6_dp,'analytic powerlaw shape MLE')
    call powerlaw_fit(data,[0.1_dp,0.0_dp,1.0_dp],[10.0_dp,0.0_dp,1.0_dp],fit)
    call check(fit%success,'mixed shape fit converges')
    call close(fit%params(1),0.6711914062499998_dp,2.0e-4_dp,'SciPy fitted shape')
    call close(fit%nllf,powerlaw_censored_nnlf(data,fit%params(1),fit%params(2),fit%params(3)), &
        0.0_dp,'fit objective agreement')
    single = make_censored_data(uncensored=[2.0_dp,3.0_dp,5.0_dp])
    call uniform_fit(single,[1.0_dp,3.0_dp],[2.0_dp,6.0_dp],fit)
    call check(fit%success,'uniform bounded all-free local fit converges')
    call close(fit%nllf,uniform_censored_nnlf(single,fit%params(1),fit%params(2)), &
        0.0_dp,'uniform free fit objective agreement')
    call uniform_fit(single,[2.0_dp,3.0_dp],[2.0_dp,6.0_dp],fit)
    call check(fit%success,'uniform fixed-location MLE converges')
    call close(fit%params(2),3.0_dp,4.0e-6_dp,'uniform analytic width MLE')
    call powerlaw_fit(data,[0.3_dp,-0.2_dp,1.0_dp],[3.0_dp,0.0_dp,2.0_dp],fit)
    call check(fit%success,'powerlaw all parameters free')
    call arcsine_fit(data,[-0.2_dp,1.0_dp],[0.0_dp,2.0_dp],fit)
    call check(fit%success,'arcsine location and width free')
    call close(fit%nllf,arcsine_censored_nnlf(data,fit%params(1),fit%params(2)), &
        0.0_dp,'arcsine fit objective agreement')
    call uniform_fit(empty,[0.0_dp,1.0_dp],[0.0_dp,1.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'uniform empty fit rejected')
    call powerlaw_fit(empty,[1.0_dp,0.0_dp,1.0_dp],[1.0_dp,0.0_dp,1.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'powerlaw empty fit rejected')
    call arcsine_fit(empty,[0.0_dp,1.0_dp],[0.0_dp,1.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'arcsine empty fit rejected')
    single = make_censored_data(uncensored=[2.0_dp])
    call uniform_fit(single,[0.0_dp,1.0_dp],[0.0_dp,1.0_dp],fit)
    call check(fit%status == fit_status_no_finite_objective,'uniform impossible fit')
    call powerlaw_fit(single,[1.0_dp,0.0_dp,1.0_dp],[1.0_dp,0.0_dp,1.0_dp],fit)
    call check(fit%status == fit_status_no_finite_objective,'powerlaw impossible fit')
    call arcsine_fit(single,[0.0_dp,1.0_dp],[0.0_dp,1.0_dp],fit)
    call check(fit%status == fit_status_no_finite_objective,'arcsine impossible fit')
    call powerlaw_fit(data,[0.1_dp,0.0_dp,1.0_dp],[10.0_dp,0.0_dp,1.0_dp],fit, &
        guess=[1.0_dp,0.0_dp,1.0_dp],max_iter=1,tolerance=1.0e-10_dp)
    call check(fit%status == fit_status_max_iter,'powerlaw optional controls')
    call arcsine_fit(data,[-0.2_dp,1.0_dp],[0.0_dp,2.0_dp],fit, &
        guess=[-0.1_dp,1.5_dp],max_iter=1,tolerance=1.0e-10_dp)
    call check(fit%status == fit_status_max_iter,'arcsine optional controls')
    call uniform_fit(data,[-0.2_dp,1.0_dp],[0.0_dp,2.0_dp],fit, &
        guess=[-0.1_dp,1.5_dp],max_iter=1,tolerance=1.0e-10_dp)
    call check(fit%status == fit_status_max_iter,'uniform optional controls')
    call close(uniform_censored_loglikelihood(data),powerlaw_censored_loglikelihood(data,1.0_dp), &
        0.0_dp,'uniform powerlaw identity')
    call close(arcsine_censored_loglikelihood(data),beta_censored_loglikelihood(data,0.5_dp,0.5_dp), &
        3.0e-14_dp,'arcsine beta identity')
    call close(powerlaw_censored_loglikelihood(data,0.6_dp), &
        beta_censored_loglikelihood(data,0.6_dp,1.0_dp),3.0e-14_dp,'powerlaw beta identity')
    single = make_censored_data(interval=reshape([2.0e-100_dp,1.0_dp],[1,2]))
    call close(uniform_censored_loglikelihood(single,1.0e-100_dp)/(-2.0e-100_dp), &
        1.0_dp,5.0e-14_dp,'uniform near-unit interval ratio')
    expected = -4.0_dp*1.0e-50_dp/acos(-1.0_dp)
    call close(arcsine_censored_loglikelihood(single,1.0e-100_dp)/expected, &
        1.0_dp,2.0e-13_dp,'arcsine near-unit interval ratio')
    single = make_censored_data(right=[0.5_dp])
    call close(powerlaw_censored_loglikelihood(single,small), &
        log(small)+log(log(2.0_dp)),2.0e-13_dp,'smallest shape SF')
    single = make_censored_data(uncensored=[1.0e-300_dp])
    call close(powerlaw_censored_loglikelihood(single,0.5_dp,scale=1.0e300_dp), &
        -log(2.0_dp),3.0e-13_dp,'underflowed ratio density')
    call close(arcsine_censored_loglikelihood(single,scale=1.0e300_dp), &
        -log(acos(-1.0_dp)),3.0e-13_dp,'underflowed arcsine ratio density')
    single = make_censored_data(left=[0.01_dp])
    call check(powerlaw_censored_loglikelihood(single,huge(1.0_dp)) == -inf, &
        'unrepresentable log CDF')
    single = make_censored_data(left=[1.0_dp])
    call close(powerlaw_censored_loglikelihood(single,1.0e308_dp,1.0e-300_dp), &
        -1.0e8_dp,1.0e-7_dp,'large shape at compensated upper edge')
    single = make_censored_data(left=[1.0e300_dp])
    actual = powerlaw_censored_loglikelihood(single,1.0e300_dp,1.0e-300_dp,1.0e300_dp)
    call close(actual/(-1.0e-300_dp),1.0_dp,3.0e-13_dp,'compensated underflowed upper distance')
    single = make_censored_data(right=[1.0e300_dp])
    call close(powerlaw_censored_loglikelihood(single,1.0e300_dp,1.0e-300_dp,1.0e300_dp), &
        log(1.0e-300_dp),3.0e-13_dp,'finite SF log at underflowed upper distance')
    write(*,'(a,i0,a)') 'Bounded family censored: ',checks,' checks'
    if (failures /= 0) error stop 'Bounded family censoring failures'
contains
    pure function likelihood(data, family, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        integer, intent(in) :: family !! 1 uniform, 2 powerlaw, 3 arcsine
        real(dp), intent(in) :: a !! power-function shape, ignored by other families
        real(dp), intent(in), optional :: loc !! lower support endpoint
        real(dp), intent(in), optional :: scale !! positive support width
        real(dp) :: value
        select case (family)
        case (1)
            value = uniform_censored_loglikelihood(data,loc,scale)
        case (2)
            value = powerlaw_censored_loglikelihood(data,a,loc,scale)
        case default
            value = arcsine_censored_loglikelihood(data,loc,scale)
        end select
    end function likelihood
    subroutine check(condition,label)
        implicit none
        logical, intent(in) :: condition !! true when expected behavior holds
        character(*), intent(in) :: label !! failure diagnostic
        checks = checks+1
        if (condition) return
        failures = failures+1
        write(*,'(a)') 'FAIL: '//label
    end subroutine check
    subroutine close(value,reference,tolerance,label)
        implicit none
        real(dp), intent(in) :: value !! computed finite result
        real(dp), intent(in) :: reference !! expected finite value
        real(dp), intent(in) :: tolerance !! absolute tolerance for this regime
        character(*), intent(in) :: label !! failure diagnostic
        call check(ieee_is_finite(value),label//' finite')
        if (.not. ieee_is_finite(value)) return
        call check(abs(value-reference) <= tolerance,label)
        if (abs(value-reference) > tolerance) write(*,*) value,reference
    end subroutine close
end program test_censored_bounded_family
