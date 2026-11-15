! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_half_families
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        halfcauchy_censored_loglikelihood, halfcauchy_censored_nnlf, halfcauchy_fit, &
        halflogistic_censored_loglikelihood, halflogistic_censored_nnlf, halflogistic_fit, &
        halfcauchy_logpdf, halfcauchy_logcdf, halfcauchy_logsf, halfcauchy_loglikelihood, &
        halflogistic_logpdf, halflogistic_logcdf, halflogistic_logsf, halflogistic_loglikelihood, &
        cauchy_censored_loglikelihood, logistic_censored_loglikelihood, &
        fit_status_invalid_input, fit_status_no_finite_objective, fit_status_max_iter
    use scifort_math, only : quiet_nan, positive_infinity
    use scifort_constants, only : scifort_pi
    implicit none
    type(censored_data) :: data, other, empty
    type(fit_result) :: fit
    real(dp), parameter :: mixed(2) = [-14.499835370067348_dp,-14.320770635057839_dp]
    real(dp), parameter :: fitted_scale(2) = [1.4831054687500016_dp,1.575110174798021_dp]
    real(dp) :: inf, nan, logtwo, logpi, width, base, value, expected, x, a, b, logorigin
    integer :: family, checks, failures, i
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    logtwo = log(2.0_dp)
    logpi = log(scifort_pi)
    width = nearest(0.0_dp,1.0_dp)
    do family = 1, 2
        logorigin = -logtwo
        if (family == 1) logorigin = logtwo-logpi
        data = make_censored_data(uncensored=[0.5_dp,1.2_dp,2.5_dp],left=[0.3_dp], &
            right=[3.5_dp],interval=reshape([1.1_dp,4.0_dp,1.4_dp,4.4_dp],[2,2]))
        call close(ll(data,0.0_dp,1.3_dp),mixed(family),4.0e-14_dp,'mixed mpmath reference')
        if (family == 1) then
            value = halfcauchy_censored_nnlf(data,0.0_dp,1.3_dp)
        else
            value = halflogistic_censored_nnlf(data,0.0_dp,1.3_dp)
        end if
        call close(value,-mixed(family),4.0e-14_dp,'negative log-likelihood')
        call fit_data(data,[0.0_dp,0.05_dp],[0.0_dp,10.0_dp],fit)
        call check(fit%success,'fixed-location mixed fit converges')
        call close(fit%params(2),fitted_scale(family),4.0e-4_dp,'SciPy fitted scale')
        call close(fit%nllf,-ll(data,fit%params(1),fit%params(2)),0.0_dp,'fit objective')
        call fit_data(data,[0.0_dp,1.3_dp],[0.0_dp,1.3_dp],fit)
        call check(fit%success,'both parameters fixed')
        call close(fit%nllf,-mixed(family),4.0e-14_dp,'fixed objective')
        base = fit%nllf
        call fit_data(data,[-3.0_dp,0.05_dp],[0.29_dp,10.0_dp],fit)
        call check(fit%success,'both parameters free')
        call check(fit%nllf <= base,'free fit improves fixed objective')
        call fit_data(data,[-3.0_dp,1.3_dp],[0.29_dp,1.3_dp],fit)
        call check(fit%success .and. fit%params(2) == 1.3_dp,'fixed scale fit')
        call check(fit%nllf <= base,'free location improves objective')
        call fit_data(data,[-3.0_dp,0.05_dp],[0.29_dp,10.0_dp],fit, &
            guess=[0.0_dp,1.3_dp],max_iter=1)
        call check(fit%status == fit_status_max_iter,'iteration limit forwarded')
        other = data
        deallocate(other%left)
        if (family == 1) then
            expected = cauchy_censored_loglikelihood(other,0.0_dp,1.3_dp)
        else
            expected = logistic_censored_loglikelihood(other,0.0_dp,1.3_dp)
        end if
        call close(ll(other,0.0_dp,1.3_dp),expected+6.0_dp*logtwo, &
            4.0e-14_dp,'positive-side parent identity')
        data = make_censored_data(uncensored=[0.0_dp,0.5_dp,1.0_dp,3.0_dp])
        if (family == 1) then
            value = halfcauchy_loglikelihood(data%uncensored,0.0_dp,1.3_dp)
            call halfcauchy_fit(data%uncensored,[0.0_dp,0.05_dp],[0.0_dp,10.0_dp],fit)
        else
            value = halflogistic_loglikelihood(data%uncensored,0.0_dp,1.3_dp)
            call halflogistic_fit(data%uncensored,[0.0_dp,0.05_dp],[0.0_dp,10.0_dp],fit)
        end if
        call close(ll(data,0.0_dp,1.3_dp),value,2.0e-14_dp,'array likelihood agreement')
        base = fit%nllf
        call fit_data(data,[0.0_dp,0.05_dp],[0.0_dp,10.0_dp],fit)
        call check(fit%success,'exact-data fit converges')
        call close(fit%nllf,base,2.0e-12_dp,'array fit agreement')
        do i = 1, 10
            x = 0.5_dp*real(i,dp)
            data = make_censored_data(uncensored=[x])
            call close(ll(data,0.1_dp,1.3_dp),scalar(x,1),8.0e-15_dp,'scalar density')
            data = make_censored_data(left=[x])
            call close(ll(data,0.1_dp,1.3_dp),scalar(x,2),8.0e-15_dp,'scalar CDF')
            data = make_censored_data(right=[x])
            call close(ll(data,0.1_dp,1.3_dp),scalar(x,3),8.0e-15_dp,'scalar SF')
        end do
        data = make_censored_data(uncensored=[0.0_dp])
        call close(ll(data),logorigin,4.0e-16_dp,'density at support endpoint')
        data = make_censored_data(uncensored=[-1.0_dp])
        call check(ll(data) == -inf,'density outside support')
        data = make_censored_data(left=[0.0_dp])
        call check(ll(data) == -inf,'zero CDF at endpoint')
        data = make_censored_data(left=[-1.0_dp])
        call check(ll(data) == -inf,'zero CDF below support')
        data = make_censored_data(right=[0.0_dp,-1.0_dp])
        call close(ll(data),0.0_dp,0.0_dp,'unit SF at/below support')
        call interval(-1.0_dp,0.0_dp)
        call check(ll(data) == -inf,'interval outside support')
        call interval(0.0_dp,0.0_dp)
        call close(ll(data),logorigin,4.0e-16_dp,'equal endpoint is exact observation')
        call interval(-1.0_dp,1.0_dp)
        base = ll(data)
        call interval(0.0_dp,1.0_dp)
        call close(ll(data),base,0.0_dp,'clip interval to support')
        data = make_censored_data(left=[1.0_dp])
        call close(ll(data),base,0.0_dp,'clipped interval equals CDF')
        call interval(0.0_dp,2.0_dp)
        call check(ll(data) >= base .and. ll(data) <= 0.0_dp,'nested interval mass')
        call interval(0.0_dp,width)
        call close(ll(data),log(width)+logorigin,2.0e-13_dp,'smallest width at endpoint')
        call interval(-width,width)
        call close(ll(data),log(width)+logorigin,2.0e-13_dp,'smallest clipped width')
        call interval(1.0e-300_dp,2.0e-300_dp)
        call close(ll(data,scale=1.0e300_dp),log(1.0e-300_dp)-log(1.0e300_dp)+logorigin, &
            4.0e-13_dp,'standardized width underflows')
        data = make_censored_data(right=[1.0e-300_dp])
        expected = -exp(logorigin)*1.0e-300_dp
        call close(ll(data)/expected,1.0_dp,3.0e-13_dp,'near-one SF at boundary')
        data = make_censored_data(left=[1.0e-300_dp])
        call close(ll(data),log(1.0e-300_dp)+logorigin,2.0e-13_dp,'tiny boundary CDF')
        a = 1000.0_dp
        b = nearest(a,1.0_dp)
        call interval(a,b)
        if (family == 1) then
            expected = logtwo-logpi-log(1.0_dp+a*a)+log(b-a)
        else
            expected = logtwo-a+log(b-a)
        end if
        call close(ll(data),expected,4.0e-13_dp,'adjacent far-tail interval')
        call interval(-1.0e308_dp,1.0e308_dp)
        if (family == 1) then
            expected = -logtwo
        else
            expected = log(tanh(0.5_dp))
        end if
        call close(ll(data,scale=1.0e308_dp),expected,4.0e-14_dp,'overflowing raw width')
        data = make_censored_data(uncensored=[1.0e308_dp])
        if (family == 1) then
            expected = logtwo-logpi-log(5.0_dp)-log(1.0e308_dp)
        else
            expected = logtwo-2.0_dp-2.0_dp*log(1.0_dp+exp(-2.0_dp))-log(1.0e308_dp)
        end if
        call close(ll(data,-1.0e308_dp,1.0e308_dp),expected,2.0e-13_dp,'raw displacement overflow')
        if (family == 1) then
            expected = logtwo-logpi+log(1.0e-300_dp)-2.0_dp*log(1.0e308_dp)
            call close(ll(data,scale=1.0e-300_dp),expected,8.0e-13_dp,'finite heavy-tail log density')
            data = make_censored_data(left=[1.0e300_dp])
            expected = -2.0_dp/(scifort_pi*1.0e300_dp)
        else
            call check(ll(data,scale=1.0e-300_dp) == -inf,'unrepresentable exponential log density')
            call interval(1.0e308_dp,nearest(1.0e308_dp,1.0_dp))
            call check(ll(data,scale=1.0e-300_dp) == -inf,'unrepresentable interval log mass')
            data = make_censored_data(left=[700.0_dp])
            expected = -2.0_dp*exp(-700.0_dp)
        end if
        call close(ll(data)/expected,1.0_dp,3.0e-13_dp,'near-one CDF keeps tiny log probability')
        if (family == 1) then
            call interval(1.0e-300_dp,1.0e300_dp)
            expected = -4.0_dp/(scifort_pi*1.0e300_dp)
        else
            call interval(1.0e-300_dp,700.0_dp)
            expected = -0.5e-300_dp-2.0_dp*exp(-700.0_dp)
        end if
        call close(ll(data)/expected,1.0_dp,3.0e-13_dp,'almost whole support retains excluded mass')
        call interval(-inf,1.0_dp)
        base = ll(data)
        data = make_censored_data(left=[1.0_dp])
        call close(ll(data),base,0.0_dp,'infinite lower endpoint')
        call interval(0.0_dp,inf)
        call close(ll(data),0.0_dp,0.0_dp,'whole support interval')
        call close(ll(empty),0.0_dp,0.0_dp,'empty valid sample')
        call check(ieee_is_nan(ll(empty,scale=0.0_dp)),'invalid scale even for empty sample')
        call check(ieee_is_nan(ll(data,scale=-1.0_dp)),'negative scale')
        call check(ieee_is_nan(ll(data,scale=inf)),'infinite scale')
        call check(ieee_is_nan(ll(data,loc=inf)),'infinite location')
        call check(ieee_is_nan(ll(data,scale=nan)),'NaN scale')
        call check(ieee_is_nan(ll(data,loc=nan)),'NaN location')
        data = make_censored_data(uncensored=[nan])
        call check(ieee_is_nan(ll(data)),'NaN observation')
        data = make_censored_data(uncensored=[inf])
        call check(ieee_is_nan(ll(data)),'infinite exact observation')
        call interval(1.0_dp,0.0_dp)
        call check(ieee_is_nan(ll(data)),'reversed interval')
        call interval(-inf,inf)
        call check(ieee_is_nan(ll(data)),'two infinite endpoints')
        call fit_data(data,[-3.0_dp,0.05_dp],[0.29_dp,10.0_dp],fit)
        call check(fit%status == fit_status_invalid_input,'invalid fit data')
        call fit_data(empty,[-3.0_dp,0.05_dp],[0.29_dp,10.0_dp],fit)
        call check(fit%status == fit_status_invalid_input,'empty fit data')
        data = make_censored_data(uncensored=[0.0_dp])
        call fit_data(data,[0.0_dp,0.0_dp],[0.0_dp,0.0_dp],fit)
        call check(fit%status == fit_status_no_finite_objective,'invalid fixed scale objective')
        data = make_censored_data(left=[0.0_dp])
        call fit_data(data,[0.0_dp,0.1_dp],[0.0_dp,2.0_dp],fit)
        call check(fit%status == fit_status_no_finite_objective,'impossible left-censoring objective')
        data = make_censored_data(right=[3.0_dp])
        call fit_data(data,[-5.0_dp,1.0_dp],[8.0_dp,1.0_dp],fit)
        call check(fit%success,'all-right-censored bounded fit')
        call close(fit%nllf,0.0_dp,0.0_dp,'one-sided fit reaches unit-probability region')
    end do
    write(*,'(a,i0,a,i0)') 'Half-family censored checks: ',checks,', failures: ',failures
    if (failures /= 0) error stop 1
contains
    pure function ll(sample,loc,scale) result(value)
        implicit none
        type(censored_data), intent(in) :: sample !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        if (family == 1) then
            value = halfcauchy_censored_loglikelihood(sample,loc,scale)
        else
            value = halflogistic_censored_loglikelihood(sample,loc,scale)
        end if
    end function ll
    function scalar(x,which) result(value)
        implicit none
        real(dp), intent(in) :: x !! point above support endpoint 0.1
        integer, intent(in) :: which !! select log density, CDF, or SF with 1, 2, or 3
        real(dp) :: value
        if (family == 1) then
            select case (which)
            case (1)
                value = halfcauchy_logpdf(x,0.1_dp,1.3_dp)
            case (2)
                value = halfcauchy_logcdf(x,0.1_dp,1.3_dp)
            case (3)
                value = halfcauchy_logsf(x,0.1_dp,1.3_dp)
            end select
        else
            select case (which)
            case (1)
                value = halflogistic_logpdf(x,0.1_dp,1.3_dp)
            case (2)
                value = halflogistic_logcdf(x,0.1_dp,1.3_dp)
            case (3)
                value = halflogistic_logsf(x,0.1_dp,1.3_dp)
            end select
        end if
    end function scalar
    subroutine fit_data(sample,lower,upper,result,guess,max_iter)
        implicit none
        type(censored_data), intent(in) :: sample !! exact and censored observations
        real(dp), intent(in) :: lower(2) !! lower location/scale bounds
        real(dp), intent(in) :: upper(2) !! upper location/scale bounds
        type(fit_result), intent(out) :: result !! fitted parameters and diagnostics
        real(dp), intent(in), optional :: guess(2) !! initial location and scale
        integer, intent(in), optional :: max_iter !! maximum coordinate-search sweeps
        if (family == 1) then
            call halfcauchy_fit(sample,lower,upper,result,guess,max_iter)
        else
            call halflogistic_fit(sample,lower,upper,result,guess,max_iter)
        end if
    end subroutine fit_data
    subroutine interval(lower,upper)
        implicit none
        real(dp), intent(in) :: lower !! lower interval endpoint
        real(dp), intent(in) :: upper !! upper interval endpoint
        data = make_censored_data(interval=reshape([lower,upper],[1,2]))
    end subroutine interval
    subroutine check(condition,label)
        implicit none
        logical, intent(in) :: condition !! true when the expected behavior holds
        character(*), intent(in) :: label !! test description for failure diagnostics
        checks = checks+1
        if (condition) return
        failures = failures+1
        write(*,'(a,i0,a)') 'FAIL family ',family,': '//label
    end subroutine check
    subroutine close(value,reference,tolerance,label)
        implicit none
        real(dp), intent(in) :: value !! computed finite value
        real(dp), intent(in) :: reference !! expected finite result
        real(dp), intent(in) :: tolerance !! absolute error allowed for this regime
        character(*), intent(in) :: label !! test description for failure diagnostics
        call check(ieee_is_finite(value),label//' is finite')
        if (.not. ieee_is_finite(value)) return
        call check(abs(value-reference) <= tolerance,label)
        if (abs(value-reference) > tolerance) write(*,*) value,reference
    end subroutine close
end program test_censored_half_families
