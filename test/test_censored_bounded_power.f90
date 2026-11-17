! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_bounded_power
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        loguniform_censored_loglikelihood, loguniform_censored_nnlf, loguniform_fit, &
        loguniform_logpdf, loguniform_logcdf, loguniform_logsf, loguniform_ppf, &
        truncpareto_censored_loglikelihood, truncpareto_censored_nnlf, truncpareto_fit, &
        truncpareto_logpdf, truncpareto_logcdf, truncpareto_logsf, truncpareto_ppf, &
        genhalflogistic_censored_loglikelihood, genhalflogistic_censored_nnlf, &
        genhalflogistic_fit, genhalflogistic_logpdf, genhalflogistic_logcdf, &
        genhalflogistic_logsf, genhalflogistic_ppf, fit_status_invalid_input
    use scifort_math, only : positive_infinity, quiet_nan
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    implicit none
    type(censored_data) :: data, one, other, empty
    type(fit_result) :: fit, old
    real(dp) :: x, p, shape, actual, expected, previous, inf, nan, lower(4), upper(4), t
    integer :: family, i, j, checks, failures
    real(dp), parameter :: scipy_shape(3) = &
        [0.93061758996822452_dp,2.1125005039405504_dp,0.38009762494385113_dp]
    real(dp), parameter :: scipy_objective(3) = &
        [10.427146663355687_dp,9.2205606831064895_dp,11.999191342035351_dp]
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    empty = make_censored_data()
    data = make_censored_data(uncensored=[1.2_dp,1.4_dp,1.6_dp],left=[1.1_dp], &
        right=[1.7_dp],interval=reshape([1.25_dp,1.5_dp,1.35_dp,1.55_dp],[2,2]))
    do family = 1, 3
        do j = 1, 3
            shape = real(j,dp)/4.0_dp
            if (family == 2) shape = real(2*j-3,dp)
            previous = -inf
            do i = 1, 19
                p = real(i,dp)/20.0_dp
                select case (family)
                case (1)
                    x = shape+(3.0_dp-shape)*p
                    expected = loguniform_logpdf(x,shape,3.0_dp)
                case (2)
                    x = 1.0_dp+2.0_dp*p
                    expected = truncpareto_logpdf(x,shape,3.0_dp)
                case (3)
                    x = p/shape
                    expected = genhalflogistic_logpdf(x,shape)
                end select
                one = make_censored_data(uncensored=[x])
                call close(likelihood(one,family,shape),expected,2.0e-12_dp,'scalar density')
                one = make_censored_data(left=[x])
                select case (family)
                case (1)
                    expected = loguniform_logcdf(x,shape,3.0_dp)
                case (2)
                    expected = truncpareto_logcdf(x,shape,3.0_dp)
                case (3)
                    expected = genhalflogistic_logcdf(x,shape)
                end select
                actual = likelihood(one,family,shape)
                call close(actual,expected,2.0e-12_dp,'scalar CDF')
                call check(actual >= previous .and. actual <= 0.0_dp,'CDF monotonicity')
                previous = actual
                one = make_censored_data(right=[x])
                select case (family)
                case (1)
                    expected = loguniform_logsf(x,shape,3.0_dp)
                case (2)
                    expected = truncpareto_logsf(x,shape,3.0_dp)
                case (3)
                    expected = genhalflogistic_logsf(x,shape)
                end select
                call close(likelihood(one,family,shape),expected,2.0e-12_dp,'scalar SF')
                one = make_censored_data(interval=reshape([x,nearest(x,1.0_dp)],[1,2]))
                other = make_censored_data(uncensored=[x])
                expected = likelihood(other,family,shape)+log(nearest(x,1.0_dp)-x)
                call close(likelihood(one,family,shape),expected,3.0e-13_dp,'adjacent interval')
            end do
            do i = 1, 9
                p = real(i,dp)/10.0_dp
                select case (family)
                case (1)
                    x = loguniform_ppf(p,shape,3.0_dp)
                case (2)
                    x = truncpareto_ppf(p,shape,3.0_dp)
                case (3)
                    x = genhalflogistic_ppf(p,shape)
                end select
                one = make_censored_data(left=[x])
                call close(likelihood(one,family,shape),log(p),3.0e-12_dp,'CDF/PPF round trip')
            end do
        end do
        call close(likelihood(empty,family,0.5_dp),0.0_dp,0.0_dp,'empty sample')
        call check(ieee_is_nan(likelihood(empty,family,nan)),'NaN shape')
        call check(ieee_is_nan(likelihood(empty,family,inf)),'infinite shape')
        call check(ieee_is_nan(likelihood(empty,family,0.0_dp)),'zero shape')
        if (family /= 2) call check(ieee_is_nan(likelihood(empty,family,-1.0_dp)), &
            'negative shape')
        call check(ieee_is_nan(likelihood(empty,family,0.5_dp,scale=0.0_dp)),'zero scale')
        call check(ieee_is_nan(likelihood(empty,family,0.5_dp,loc=nan)),'NaN location')
        one = make_censored_data(uncensored=[nan])
        call check(ieee_is_nan(likelihood(one,family,0.5_dp)),'NaN observation')
        one = make_censored_data(uncensored=[inf])
        call check(ieee_is_nan(likelihood(one,family,0.5_dp)),'infinite observation')
        one = make_censored_data(uncensored=[-1.0_dp])
        call check(likelihood(one,family,0.5_dp) == -inf,'outside support density')
        one = make_censored_data(left=[-1.0_dp])
        call check(likelihood(one,family,0.5_dp) == -inf,'outside support CDF')
        one = make_censored_data(right=[-1.0_dp])
        call close(likelihood(one,family,0.5_dp),0.0_dp,0.0_dp,'outside support SF')
        one = make_censored_data(interval=reshape([-inf,1.5_dp],[1,2]))
        other = make_censored_data(left=[1.5_dp])
        call close(likelihood(one,family,0.5_dp),likelihood(other,family,0.5_dp), &
            0.0_dp,'infinite interval endpoint')
        one = make_censored_data(interval=reshape([1.5_dp,1.5_dp],[1,2]))
        other = make_censored_data(uncensored=[1.5_dp])
        call close(likelihood(one,family,0.5_dp),likelihood(other,family,0.5_dp), &
            0.0_dp,'equal interval endpoints')
        lower = [0.5_dp,3.0_dp,0.0_dp,1.0_dp]
        call fitting(data,family,lower,lower,fit)
        call check(fit%success,'fixed parameter fit')
        call close(fit%nllf,-likelihood(data,family,0.5_dp),0.0_dp,'fixed fit objective')
        one = make_censored_data(uncensored=[1.2_dp,1.4_dp,1.6_dp])
        call fitting(one,family,lower,lower,fit)
        select case (family)
        case (1)
            call loguniform_fit(one%uncensored,lower,lower,old)
            actual = loguniform_censored_nnlf(data,0.5_dp,3.0_dp)
        case (2)
            call truncpareto_fit(one%uncensored,lower,lower,old)
            actual = truncpareto_censored_nnlf(data,0.5_dp,3.0_dp)
        case (3)
            call genhalflogistic_fit(one%uncensored,lower([1,3,4]),lower([1,3,4]),old)
            actual = genhalflogistic_censored_nnlf(data,0.5_dp)
        end select
        call check(old%success,'original array overload')
        call close(fit%nllf,old%nllf,2.0e-12_dp,'array/censored fit agreement')
        call close(actual,-likelihood(data,family,0.5_dp),0.0_dp,'negative likelihood')
        call fitting(empty,family,lower,lower,fit)
        call check(fit%status == fit_status_invalid_input,'empty fit rejected')
        upper = lower
        if (family == 1) then
            lower(1) = 0.1_dp
            upper(1) = 1.09_dp
        else if (family == 2) then
            lower(1) = 0.01_dp
            upper(1) = 5.0_dp
        else
            lower(1) = 0.05_dp
            upper(1) = 0.58_dp
        end if
        call fitting(data,family,lower,upper,fit)
        call check(fit%success,'free shape censored fit')
        call close(fit%nllf,-likelihood(data,family,fit%params(1)),0.0_dp,'free fit objective')
        call close(fit%params(1),scipy_shape(family),3.0e-5_dp,'SciPy shape optimum')
        call close(fit%nllf,scipy_objective(family),3.0e-10_dp,'SciPy optimal objective')
    end do
    call check(ieee_is_nan(loguniform_censored_loglikelihood(empty,1.0_dp,1.0_dp)), &
        'equal log-uniform shapes')
    call check(ieee_is_nan(truncpareto_censored_loglikelihood(empty,1.0_dp,1.0_dp)), &
        'unit Pareto cutoff')
    x = 4.575698893664817e200_dp
    p = huge(1.0_dp)/x
    call check(ieee_is_nan(loguniform_censored_loglikelihood(empty,1.0_dp,x,scale=p)), &
        'rounded quotient admits overflowing log-uniform cutoff')
    call check(ieee_is_nan(truncpareto_censored_loglikelihood(empty,1.0_dp,x,scale=p)), &
        'rounded quotient admits overflowing Pareto cutoff')
    call close(loguniform_censored_loglikelihood(empty,1.0_dp,x, &
        scale=nearest(p,-1.0_dp)),0.0_dp,0.0_dp,'neighboring finite log-uniform cutoff')
    call close(truncpareto_censored_loglikelihood(empty,1.0_dp,x, &
        scale=nearest(p,-1.0_dp)),0.0_dp,0.0_dp,'neighboring finite Pareto cutoff')
    call check(ieee_is_nan(loguniform_censored_loglikelihood(empty,1.0_dp,2.0_dp, &
        scale=huge(1.0_dp))),'overflowed cutoffs')
    call check(ieee_is_nan(truncpareto_censored_loglikelihood(empty,1.0_dp,2.0_dp, &
        loc=huge(1.0_dp))),'collapsed translated support')
    call check(ieee_is_nan(genhalflogistic_censored_loglikelihood(empty,0.5_dp, &
        scale=huge(1.0_dp))),'overflowed generalized half-logistic width')
    one = make_censored_data(uncensored=[1.0_dp])
    call check(genhalflogistic_censored_loglikelihood(one,2.0_dp,scale=2.0_dp) == inf, &
        'singular upper density')
    call check(genhalflogistic_censored_loglikelihood(one,0.5_dp,scale=0.5_dp) == -inf, &
        'vanishing upper density')
    call close(genhalflogistic_censored_loglikelihood(one,1.0_dp),log(2.0_dp), &
        1.0e-14_dp,'finite upper density')
    one = make_censored_data(uncensored=[0.25_dp])
    call close(genhalflogistic_censored_loglikelihood(one,1.0_dp), &
        log(2.0_dp)-2.0_dp*log(1.75_dp),1.0e-14_dp,'unit shape density identity')
    one = make_censored_data(left=[0.25_dp])
    call close(genhalflogistic_censored_loglikelihood(one,1.0_dp), &
        log(0.25_dp/1.75_dp),1.0e-14_dp,'unit shape CDF identity')
    t = nearest(0.0_dp,1.0_dp)
    call close(loguniform_censored_loglikelihood(empty,1.49_dp,2.49_dp,scale=t), &
        0.0_dp,0.0_dp,'single-rounded subnormal cutoff products')
    one = make_censored_data(left=[t])
    call close(genhalflogistic_censored_loglikelihood(one,t,scale=t), &
        log((1.0_dp-exp(-1.0_dp))/(1.0_dp+exp(-1.0_dp))),3.0e-13_dp, &
        'smallest shape and scale')
    one = make_censored_data(left=[1.0e-300_dp])
    call close(genhalflogistic_censored_loglikelihood(one,1.0e300_dp,scale=1.0e300_dp), &
        log(1.0e-300_dp)-log(1.0e300_dp)-log(2.0_dp),3.0e-13_dp, &
        'underflowed transformed hazard')
    one = make_censored_data(uncensored=[1.5_dp])
    actual = loguniform_censored_loglikelihood(one,1.0_dp,2.0_dp)
    call close(truncpareto_censored_loglikelihood(one,1.0e-300_dp,2.0_dp),actual, &
        3.0e-13_dp,'positive power zero limit')
    call close(truncpareto_censored_loglikelihood(one,-1.0e-300_dp,2.0_dp),actual, &
        3.0e-13_dp,'negative power zero limit')
    call check(ieee_is_nan(loguniform_censored_loglikelihood(empty,1.0_dp,nan)), &
        'NaN log-uniform upper shape')
    call check(ieee_is_nan(truncpareto_censored_loglikelihood(empty,1.0_dp,inf)), &
        'infinite Pareto upper shape')
    call check(ieee_is_nan(loguniform_censored_loglikelihood(empty,t,2.0_dp*t,scale=t)), &
        'underflowed lower cutoff')
    call check(ieee_is_nan(genhalflogistic_censored_loglikelihood(empty,2.0_dp,scale=t)), &
        'underflowed generalized half-logistic width')
    call check(ieee_is_nan(loguniform_censored_loglikelihood(empty,1.0_dp,2.0_dp, &
        loc=1.0e300_dp)),'collapsed log-uniform support')
    call check(ieee_is_nan(genhalflogistic_censored_loglikelihood(empty,1.0_dp, &
        loc=1.0e300_dp)),'collapsed generalized half-logistic support')
    one = make_censored_data(interval=reshape([0.0_dp,4.0_dp],[1,2]))
    do family = 1, 3
        call close(likelihood(one,family,0.5_dp),0.0_dp,0.0_dp,'whole support interval')
    end do
    one = make_censored_data(uncensored=[1.0_dp])
    call close(loguniform_censored_loglikelihood(one,1.0_dp,nearest(1.0_dp,1.0_dp)), &
        -log(log(nearest(1.0_dp,1.0_dp))),2.0e-14_dp,'adjacent support cutoffs')
    one = make_censored_data(interval=reshape([1.5_dp,nearest(1.5_dp,1.0_dp)],[1,2]))
    actual = loguniform_censored_loglikelihood(one,1.0_dp,2.0_dp)
    call close(truncpareto_censored_loglikelihood(one,1.0e-300_dp,2.0_dp),actual, &
        3.0e-13_dp,'tiny positive power interval')
    call close(truncpareto_censored_loglikelihood(one,-1.0e-300_dp,2.0_dp),actual, &
        3.0e-13_dp,'tiny negative power interval')
    call loguniform_fit(data,[0.1_dp,2.0_dp,-0.1_dp,0.9_dp], &
        [1.0_dp,4.0_dp,0.1_dp,1.1_dp],fit,max_iter=10000)
    call check(fit%success,'four free log-uniform parameters')
    call close(fit%nllf,loguniform_censored_nnlf(data,fit%params(1),fit%params(2), &
        fit%params(3),fit%params(4)),0.0_dp,'four parameter objective')
    call truncpareto_fit(data,[0.1_dp,2.0_dp,-0.1_dp,0.9_dp], &
        [3.0_dp,4.0_dp,0.1_dp,1.1_dp],fit,max_iter=10000)
    call check(fit%success,'four free truncated Pareto parameters')
    call close(fit%nllf,truncpareto_censored_nnlf(data,fit%params(1),fit%params(2), &
        fit%params(3),fit%params(4)),0.0_dp,'four parameter Pareto objective')
    call truncpareto_fit(data,[-5.0_dp,3.0_dp,0.0_dp,1.0_dp], &
        [-0.1_dp,3.0_dp,0.0_dp,1.0_dp],fit,max_iter=10000)
    call check(fit%success,'negative Pareto shape fit')
    call close(fit%nllf,truncpareto_censored_nnlf(data,fit%params(1),3.0_dp), &
        0.0_dp,'negative Pareto fit objective')
    print '(a,i0)', 'Censored bounded-power checks: ',checks
    if (failures /= 0) error stop 'censored bounded-power failures'
contains
    function likelihood(sample, family, shape, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: sample !! exact and censored observations
        integer, intent(in) :: family !! 1 loguniform, 2 truncpareto, 3 genhalflogistic
        real(dp), intent(in) :: shape !! first shape; Pareto permits either nonzero sign
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        select case (family)
        case (1)
            value = loguniform_censored_loglikelihood(sample,shape,3.0_dp,loc,scale)
        case (2)
            value = truncpareto_censored_loglikelihood(sample,shape,3.0_dp,loc,scale)
        case (3)
            value = genhalflogistic_censored_loglikelihood(sample,shape,loc,scale)
        case default
            value = quiet_nan(0.0_dp)
        end select
    end function likelihood
    subroutine fitting(sample, family, lower, upper, result)
        implicit none
        type(censored_data), intent(in) :: sample !! exact and censored observations
        integer, intent(in) :: family !! family selector as in likelihood
        real(dp), intent(in) :: lower(4) !! shape, second shape, location, scale bounds
        real(dp), intent(in) :: upper(4) !! upper bounds; equality fixes a parameter
        type(fit_result), intent(out) :: result !! parameters and fit diagnostics
        select case (family)
        case (1)
            call loguniform_fit(sample,lower,upper,result,max_iter=10000)
        case (2)
            call truncpareto_fit(sample,lower,upper,result,max_iter=10000)
        case (3)
            call genhalflogistic_fit(sample,lower([1,3,4]),upper([1,3,4]), &
                result,max_iter=10000)
        end select
    end subroutine fitting
    subroutine check(condition, label)
        implicit none
        logical, intent(in) :: condition !! assertion that must hold
        character(*), intent(in) :: label !! diagnostic description
        checks = checks+1
        if (condition) return
        failures = failures+1
        print '(a)', 'FAIL: '//label
    end subroutine check
    subroutine close(actual, expected, tolerance, label)
        implicit none
        real(dp), intent(in) :: actual !! computed finite value
        real(dp), intent(in) :: expected !! independently expected finite value
        real(dp), intent(in) :: tolerance !! scaled absolute tolerance
        character(*), intent(in) :: label !! diagnostic description
        call check(ieee_is_finite(actual),label//' finite')
        if (.not. ieee_is_finite(actual)) return
        call check(abs(actual-expected) <= tolerance*max(1.0_dp,abs(expected)),label)
        if (abs(actual-expected) > tolerance*max(1.0_dp,abs(expected))) &
            print '(2es25.16)', actual, expected
    end subroutine close
end program test_censored_bounded_power
