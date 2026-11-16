! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_reflected
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        dgamma_censored_loglikelihood, dgamma_censored_nnlf, dgamma_fit, &
        dgamma_logpdf, dgamma_logcdf, dgamma_logsf, dgamma_ppf, &
        dweibull_censored_loglikelihood, dweibull_censored_nnlf, dweibull_fit, &
        dweibull_logpdf, dweibull_logcdf, dweibull_logsf, dweibull_ppf, &
        loggamma_censored_loglikelihood, loggamma_censored_nnlf, loggamma_fit, &
        loggamma_logpdf, loggamma_logcdf, loggamma_logsf, loggamma_ppf, &
        laplace_censored_loglikelihood, &
        fit_status_invalid_input, fit_status_no_finite_objective, fit_status_max_iter
    use scifort_log_gamma, only : stirling_remainder
    use scifort_math, only : quiet_nan, positive_infinity
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
    implicit none
    type(censored_data) :: data, single, empty, other
    type(fit_result) :: fit, old
    real(dp) :: inf, nan, actual, expected, x, c, p, previous, reference(3), lower(3), upper(3), scipy_shape(3)
    integer :: family, i, j, checks, failures
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    reference = [-16.344657035802172_dp,-22.754745701778361_dp,-41.125493018751435_dp]
    scipy_shape = [1.8682617187500021_dp,0.8721679687499999_dp,4.708984375000009_dp]
    data = make_censored_data(uncensored=[0.5_dp,1.2_dp,2.5_dp],left=[0.3_dp], &
        right=[3.5_dp],interval=reshape([1.1_dp,4.0_dp,1.4_dp,4.4_dp],[2,2]))
    do family = 1, 3
        actual = likelihood(data,family,1.7_dp,scale=1.3_dp)
        call close(actual,reference(family),3.0e-13_dp,'mixed mpmath reference')
        select case (family)
        case (1)
            expected = dgamma_censored_nnlf(data,1.7_dp,scale=1.3_dp)
        case (2)
            expected = dweibull_censored_nnlf(data,1.7_dp,scale=1.3_dp)
        case (3)
            expected = loggamma_censored_nnlf(data,1.7_dp,scale=1.3_dp)
        end select
        call close(expected,-actual,0.0_dp,'negative likelihood')
        do j = 1, 3
            c = 0.4_dp
            if (j == 2) c = 1.0_dp
            if (j == 3) c = 3.0_dp
            previous = -inf
            do i = 1, 12
                x = real(i-7,dp)/2.0_dp
                single = make_censored_data(uncensored=[x])
                select case (family)
                case (1)
                    expected = dgamma_logpdf(x,c)
                case (2)
                    expected = dweibull_logpdf(x,c)
                case (3)
                    expected = loggamma_logpdf(x,c)
                end select
                if (ieee_is_finite(expected)) then
                    call close(likelihood(single,family,c),expected,6.0e-14_dp,'scalar density')
                else
                    call check(likelihood(single,family,c) == expected,'scalar center density')
                end if
                single = make_censored_data(left=[x])
                select case (family)
                case (1)
                    expected = dgamma_logcdf(x,c)
                case (2)
                    expected = dweibull_logcdf(x,c)
                case (3)
                    expected = loggamma_logcdf(x,c)
                end select
                actual = likelihood(single,family,c)
                call close(actual,expected,6.0e-14_dp,'scalar CDF')
                call check(actual >= previous .and. actual <= 0.0_dp,'CDF monotonicity and range')
                previous = actual
                single = make_censored_data(right=[x])
                select case (family)
                case (1)
                    expected = dgamma_logsf(x,c)
                case (2)
                    expected = dweibull_logsf(x,c)
                case (3)
                    expected = loggamma_logsf(x,c)
                end select
                call close(likelihood(single,family,c),expected,6.0e-14_dp,'scalar SF')
            end do
            do i = 1, 9
                p = real(i,dp)/10.0_dp
                select case (family)
                case (1)
                    x = dgamma_ppf(p,c)
                case (2)
                    x = dweibull_ppf(p,c)
                case (3)
                    x = loggamma_ppf(p,c)
                end select
                single = make_censored_data(left=[x])
                call close(exp(likelihood(single,family,c)),p,5.0e-14_dp,'CDF/PPF round trip')
            end do
        end do
        single = make_censored_data(uncensored=[0.5_dp,1.0_dp,2.0_dp])
        lower = [1.7_dp,0.0_dp,1.3_dp]
        upper = lower
        call fitting(single,family,lower,upper,fit)
        select case (family)
        case (1)
            call dgamma_fit(single%uncensored,lower,upper,old)
        case (2)
            call dweibull_fit(single%uncensored,lower,upper,old)
        case (3)
            call loggamma_fit(single%uncensored,lower,upper,old)
        end select
        call check(fit%success .and. old%success,'fixed fits converge')
        call close(fit%nllf,old%nllf,2.0e-13_dp,'array fit compatibility')
        call close(fit%nllf,-likelihood(single,family,1.7_dp,scale=1.3_dp), &
            0.0_dp,'fit objective agreement')
        call close(likelihood(empty,family,1.0_dp),0.0_dp,0.0_dp,'empty valid sample')
        call check(ieee_is_nan(likelihood(empty,family,nan)),'NaN shape')
        call check(ieee_is_nan(likelihood(empty,family,inf)),'infinite shape')
        call check(ieee_is_nan(likelihood(empty,family,0.0_dp)),'zero shape')
        call check(ieee_is_nan(likelihood(empty,family,-1.0_dp)),'negative shape')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,loc=nan)),'NaN location')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,loc=inf)),'infinite location')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,scale=nan)),'NaN scale')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,scale=inf)),'infinite scale')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,scale=0.0_dp)),'zero scale')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,scale=-1.0_dp)),'negative scale')
        single = make_censored_data(uncensored=[nan])
        call check(ieee_is_nan(likelihood(single,family,1.0_dp)),'NaN exact value')
        single = make_censored_data(uncensored=[inf])
        call check(ieee_is_nan(likelihood(single,family,1.0_dp)),'infinite exact value')
        single = make_censored_data(interval=reshape([1.0_dp,0.0_dp],[1,2]))
        call check(ieee_is_nan(likelihood(single,family,1.0_dp)),'reversed interval')
        single = make_censored_data(interval=reshape([-inf,1.0_dp],[1,2]))
        other = make_censored_data(left=[1.0_dp])
        call close(likelihood(single,family,1.0_dp),likelihood(other,family,1.0_dp), &
            0.0_dp,'infinite lower interval endpoint')
        single = make_censored_data(interval=reshape([1.0_dp,inf],[1,2]))
        other = make_censored_data(right=[1.0_dp])
        call close(likelihood(single,family,1.0_dp),likelihood(other,family,1.0_dp), &
            0.0_dp,'infinite upper interval endpoint')
        single = make_censored_data(interval=reshape([1.0_dp,1.0_dp],[1,2]))
        other = make_censored_data(uncensored=[1.0_dp])
        call close(likelihood(single,family,1.0_dp),likelihood(other,family,1.0_dp), &
            0.0_dp,'equal endpoints represent exact data')
        lower = [0.3_dp,-0.1_dp,0.8_dp]
        upper = [3.0_dp,0.1_dp,2.0_dp]
        call fitting(data,family,lower,upper,fit)
        call check(fit%success,'all-free bounded fit')
        if (allocated(fit%params)) then
            call check(all(fit%params >= lower) .and. all(fit%params <= upper),'fit within bounds')
            call close(fit%nllf,-likelihood(data,family,fit%params(1), &
                fit%params(2),fit%params(3)),0.0_dp,'free fit objective')
        end if
        lower = [1.7_dp,0.0_dp,1.3_dp]
        upper = [1.7_dp,0.0_dp,3.0_dp]
        call fitting(data,family,lower,upper,fit)
        call check(fit%success,'partly fixed fit')
        lower = [0.01_dp,0.0_dp,1.3_dp]
        upper = [10.0_dp,0.0_dp,1.3_dp]
        call fitting(data,family,lower,upper,fit)
        call check(fit%success,'fixed support shape fit')
        call check(fit%nllf <= -reference(family),'shape fit improves objective')
        call close(fit%params(1),scipy_shape(family),2.0e-4_dp,'SciPy constrained shape fit')
        call fitting(empty,family,lower,upper,fit)
        call check(fit%status == fit_status_invalid_input,'empty fit rejected')
        lower = [0.0_dp,0.0_dp,1.3_dp]
        upper = lower
        call fitting(data,family,lower,upper,fit)
        call check(.not. fit%success,'invalid fixed shape rejected')
    end do
    do family = 1, 2
        single = make_censored_data(left=[0.0_dp])
        call close(likelihood(single,family,0.4_dp),-log(2.0_dp),0.0_dp,'median CDF')
        single = make_censored_data(uncensored=[0.0_dp])
        call check(likelihood(single,family,0.4_dp) == inf,'singular center density')
        call check(likelihood(single,family,2.0_dp) == -inf,'zero center density')
        call close(likelihood(single,family,1.0_dp),-log(2.0_dp),0.0_dp,'unit center density')
        lower = [0.4_dp,0.0_dp,1.0_dp]
        call fitting(single,family,lower,lower,fit)
        call check(fit%status == fit_status_no_finite_objective,'singular fit objective')
        single = make_censored_data(interval=reshape([-1.0e-300_dp,2.0e-300_dp],[1,2]))
        other = make_censored_data(interval=reshape([-2.0e-300_dp,1.0e-300_dp],[1,2]))
        call close(likelihood(single,family,0.4_dp,scale=1.0e300_dp), &
            likelihood(other,family,0.4_dp,scale=1.0e300_dp),0.0_dp,'reflected narrow interval')
        single = make_censored_data(interval=reshape([-230.0_dp,231.0_dp],[1,2]))
        actual = likelihood(single,family,1.0_dp)
        expected = -0.5_dp*(exp(-230.0_dp)+exp(-231.0_dp))
        call close(actual/expected,1.0_dp,5.0e-14_dp,'both omitted tails retained')
    end do
    do family = 1, 2
        call close(likelihood(data,family,1.0_dp,scale=1.3_dp), &
            laplace_censored_loglikelihood(data,scale=1.3_dp),2.0e-13_dp,'Laplace identity')
    end do
    single = make_censored_data(interval=reshape([0.0_dp,nearest(0.0_dp,1.0_dp)],[1,2]))
    expected = log(nearest(0.0_dp,1.0_dp))-log(2.0_dp)
    call close(dweibull_censored_loglikelihood(single,1.0_dp),expected,3.0e-13_dp, &
        'smallest reflected Weibull interval')
    single = make_censored_data(uncensored=[-1.0e300_dp])
    actual = loggamma_censored_loglikelihood(single,1.0e-300_dp,scale=1.0e-300_dp)
    call close(actual/(-1.0e300_dp),1.0_dp,3.0e-13_dp,'tiny shape rescues overflowing z')
    single = make_censored_data(left=[-1.0e300_dp])
    actual = loggamma_censored_loglikelihood(single,1.0e-300_dp,scale=1.0e-300_dp)
    call close(actual/(-1.0e300_dp),1.0_dp,3.0e-13_dp,'overflowing z left tail')
    single = make_censored_data(right=[-1000.0_dp])
    actual = loggamma_censored_loglikelihood(single,1.0e-300_dp)
    call close(actual,log(1.0e-300_dp)+log(1000.0_dp-0.57721566490153286_dp), &
        3.0e-13_dp,'tiny log gamma upper tail')
    single = make_censored_data(interval=reshape([0.0_dp,nearest(0.0_dp,1.0_dp)],[1,2]))
    expected = log(nearest(0.0_dp,1.0_dp))-1.0_dp
    call close(loggamma_censored_loglikelihood(single,1.0_dp),expected,3.0e-13_dp, &
        'smallest log gamma interval')
    single = make_censored_data(uncensored=[log(1000.0_dp)])
    expected = 0.5_dp*log(1000.0_dp)-0.91893853320467274_dp &
        -1.0_dp/12000.0_dp+1.0_dp/(360.0_dp*1000.0_dp**3)
    call close(loggamma_censored_loglikelihood(single,1000.0_dp),expected,3.0e-13_dp, &
        'log gamma mode stable cancellation')
    call close(stirling_remainder(1.0e308_dp)/(1.0_dp/12.0_dp/1.0e308_dp), &
        1.0_dp,1.0e-14_dp,'large-shape Stirling remainder')
    single = make_censored_data(uncensored=[710.0_dp])
    actual = loggamma_censored_loglikelihood(single,1.0e308_dp)
    expected = -4.302034083277817e307_dp
    call close(actual/expected,1.0_dp,5.0e-13_dp,'density beyond representable exp(z)')
    single = make_censored_data(right=[710.0_dp])
    actual = loggamma_censored_loglikelihood(single,1.0e308_dp)
    call close(actual/expected,1.0_dp,5.0e-13_dp,'survival beyond representable exp(z)')
    single = make_censored_data(uncensored=[-1.0e308_dp])
    expected = -2.0_dp-log(2.0_dp)-log(1.0e308_dp)
    do family = 1, 2
        call close(likelihood(single,family,1.0_dp,loc=1.0e308_dp,scale=1.0e308_dp), &
            expected,3.0e-13_dp,'overflowing unscaled reflected distance')
    end do
    write(*,'(a,i0,a)') 'Reflected and log gamma censored: ',checks,' checks'
    if (failures /= 0) error stop 'Reflected and log gamma censoring failures'
contains
    pure function likelihood(data, family, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        integer, intent(in) :: family !! 1 double gamma, 2 double Weibull, 3 log gamma
        real(dp), intent(in) :: c !! positive shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive scale, default one
        real(dp) :: value
        select case (family)
        case (1)
            value = dgamma_censored_loglikelihood(data,c,loc,scale)
        case (2)
            value = dweibull_censored_loglikelihood(data,c,loc,scale)
        case (3)
            value = loggamma_censored_loglikelihood(data,c,loc,scale)
        end select
    end function likelihood
    subroutine fitting(data, family, lower, upper, result)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        integer, intent(in) :: family !! distribution identifier used by this test
        real(dp), intent(in) :: lower(3) !! shape, location, scale lower bounds
        real(dp), intent(in) :: upper(3) !! shape, location, scale upper bounds
        type(fit_result), intent(out) :: result !! bounded fit diagnostics
        select case (family)
        case (1)
            call dgamma_fit(data,lower,upper,result,max_iter=10000,tolerance=1.0e-7_dp)
        case (2)
            call dweibull_fit(data,lower,upper,result,max_iter=10000,tolerance=1.0e-7_dp)
        case (3)
            call loggamma_fit(data,lower,upper,result,max_iter=10000,tolerance=1.0e-7_dp)
        end select
    end subroutine fitting
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
end program test_censored_reflected
