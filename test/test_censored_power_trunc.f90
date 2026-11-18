! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_power_trunc
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        invgamma_censored_loglikelihood, weibull_censored_loglikelihood, &
        fit_status_invalid_input, &
        gengamma_censored_loglikelihood, gengamma_censored_nnlf, gengamma_fit, &
        gengamma_logpdf, gengamma_logcdf, gengamma_logsf, gengamma_ppf, gengamma_isf, &
        truncnorm_censored_loglikelihood, truncnorm_censored_nnlf, truncnorm_fit, &
        truncnorm_logpdf, truncnorm_logcdf, truncnorm_logsf, truncnorm_ppf, truncnorm_isf, &
        truncweibull_min_censored_loglikelihood, truncweibull_min_censored_nnlf, truncweibull_min_fit, &
        truncweibull_min_logpdf, truncweibull_min_logcdf, truncweibull_min_logsf, truncweibull_min_ppf, truncweibull_min_isf
    use scifort_math, only : positive_infinity, quiet_nan
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    implicit none
    type(censored_data) :: data, one, other, empty
    type(fit_result) :: fit, old
    real(dp) :: par(5), lower(5), upper(5), x, expected, actual, previous, p, inf, nan, sf
    real(dp), parameter :: scipy_parameter(3) = &
        [1.7716401277384342_dp,0.09999999778083914_dp,1.6519866276629518_dp]
    real(dp), parameter :: scipy_nllf(3) = &
        [10.216323674112758_dp,13.024419603423851_dp,10.531775107244417_dp]
    integer :: family, i, j, k, n, checks, failures
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    empty = make_censored_data()
    data = make_censored_data(uncensored=[0.5_dp,0.8_dp,1.1_dp],left=[0.4_dp], &
        right=[1.7_dp],interval=reshape([0.6_dp,0.9_dp,0.7_dp,1.0_dp],[2,2]))
    do family = 1, 3
        n = 4
        if (family == 3) n = 5
        do j = 1, 4
            call parameters(family,par)
            if (family == 1) par(2) = (-1.0_dp)**j*real(j,dp)/2.0_dp
            if (family == 3) par(1) = real(j,dp)/2.0_dp
            previous = -inf
            do i = 1, 25
                x = 0.15_dp+real(i,dp)/10.0_dp
                one = make_censored_data(uncensored=[x])
                call close(ll(one,family,par),scalar(x,family,par,1),3e-12_dp,'density')
                one = make_censored_data(left=[x])
                actual = ll(one,family,par)
                call close(actual,scalar(x,family,par,2),3e-12_dp,'CDF')
                call check(actual >= previous .and. actual <= 0.0_dp,'CDF monotonicity')
                previous = actual
                one = make_censored_data(right=[x])
                call close(ll(one,family,par),scalar(x,family,par,3),3e-12_dp,'SF')
                one = make_censored_data(interval=reshape([x,nearest(x,1.0_dp)],[1,2]))
                expected = scalar(x,family,par,1)+log(nearest(x,1.0_dp)-x)
                call close(ll(one,family,par),expected,3e-12_dp,'adjacent interval')
                one = make_censored_data(interval=reshape([x,x+0.1_dp],[1,2]))
                if (scalar(x,family,par,2) < log(0.5_dp)) then
                    expected = log(exp(scalar(x+0.1_dp,family,par,2))-exp(scalar(x,family,par,2)))
                else
                    expected = log(exp(scalar(x,family,par,3))-exp(scalar(x+0.1_dp,family,par,3)))
                end if
                call close(ll(one,family,par),expected,3e-11_dp,'ordinary interval')
            end do
            do i = 1, 9
                p = real(i,dp)/10.0_dp
                x = quantile(p,family,par,.false.)
                one = make_censored_data(left=[x])
                call close(ll(one,family,par),log(p),5e-11_dp,'CDF/PPF')
                x = quantile(p,family,par,.true.)
                one = make_censored_data(right=[x])
                call close(ll(one,family,par),log(p),5e-11_dp,'SF/ISF')
            end do
        end do
        call parameters(family,par)
        do k = -1, 1
            par(n) = 10.0_dp**(200*k)
            par(n-1) = -0.3_dp*par(n)
            x = 0.7_dp*par(n)
            one = make_censored_data(uncensored=[x])
            expected = scalar(1.0_dp,family,par,1)-log(par(n))
            call close(ll(one,family,par),expected,3e-12_dp,'affine density')
            one = make_censored_data(interval=reshape([x,nearest(x,1.0_dp)],[1,2]))
            expected = expected+log(nearest(x,1.0_dp)-x)
            call close(ll(one,family,par),expected,3e-12_dp,'affine adjacent interval')
        end do
        call parameters(family,par)
        call close(ll(empty,family,par),0.0_dp,0.0_dp,'valid empty')
        par(n) = 0.0_dp
        call check(ieee_is_nan(ll(empty,family,par)),'invalid empty scale')
        par(n) = inf
        call check(ieee_is_nan(ll(empty,family,par)),'infinite scale')
        call parameters(family,par)
        par(1) = nan
        call check(ieee_is_nan(ll(empty,family,par)),'NaN shape')
        call parameters(family,par)
        one = make_censored_data(uncensored=[inf])
        call check(ieee_is_nan(ll(one,family,par)),'infinite exact invalid')
        one = make_censored_data(left=[nan])
        call check(ieee_is_nan(ll(one,family,par)),'NaN data invalid')
        one = make_censored_data(interval=reshape([1.0_dp,0.9_dp],[1,2]))
        call check(ieee_is_nan(ll(one,family,par)),'reversed interval')
        one = make_censored_data(interval=reshape([-inf,inf],[1,2]))
        call check(ieee_is_nan(ll(one,family,par)),'whole infinite interval invalid')
        one = make_censored_data(interval=reshape([0.7_dp,0.7_dp],[1,2]))
        call close(ll(one,family,par),scalar(0.7_dp,family,par,1),3e-12_dp,'equal endpoints exact')
        one = make_censored_data(interval=reshape([-inf,0.7_dp],[1,2]))
        call close(ll(one,family,par),scalar(0.7_dp,family,par,2),3e-12_dp,'left interval')
        one = make_censored_data(interval=reshape([0.7_dp,inf],[1,2]))
        call close(ll(one,family,par),scalar(0.7_dp,family,par,3),3e-12_dp,'right interval')
        lower = par
        upper = par
        call fit_family(data,family,lower,upper,fit)
        call check(fit%success,'fixed fit converged')
        call close(fit%nllf,-ll(data,family,par),3e-13_dp,'fixed NNLF')
        call fit_family(empty,family,lower,upper,fit)
        call check(fit%status == fit_status_invalid_input,'empty fit rejected')
        one = make_censored_data(uncensored=[0.5_dp,0.8_dp,1.1_dp])
        call fit_family(one,family,lower,upper,fit)
        call fit_array([0.5_dp,0.8_dp,1.1_dp],family,lower,upper,old)
        call check(fit%success .and. old%success,'old generic retained')
        call close(fit%nllf,old%nllf,3e-12_dp,'array/exact fit agreement')
        call parameters(family,par)
        if (family == 1) par(2) = -1.3_dp
        lower = par
        upper = par
        j = 1
        if (family == 2) j = 3
        lower(j) = 0.3_dp
        upper(j) = 5.0_dp
        if (family == 2) then
            lower(j) = -0.1_dp
            upper(j) = 0.1_dp
        end if
        call fit_family(data,family,lower,upper,fit)
        call check(fit%success,'one-parameter convergence')
        call close(fit%params(j),scipy_parameter(family),3e-5_dp,'SciPy optimum parameter')
        call close(fit%nllf,scipy_nllf(family),2e-8_dp,'SciPy optimum likelihood')
        lower = par
        upper = par
        if (family == 1) then
            lower(:4) = [0.3_dp,-2.0_dp,-0.1_dp,0.7_dp]
            upper(:4) = [5.0_dp,-0.5_dp,0.1_dp,1.3_dp]
        else if (family == 2) then
            lower(:4) = [-2.0_dp,2.8_dp,-0.1_dp,0.7_dp]
            upper(:4) = [-0.5_dp,4.0_dp,0.1_dp,1.3_dp]
        else
            lower = [0.5_dp,0.0_dp,2.8_dp,-0.1_dp,0.7_dp]
            upper = [3.0_dp,0.2_dp,4.0_dp,0.1_dp,1.3_dp]
        end if
        call fit_family(data,family,lower,upper,fit)
        call check(ieee_is_finite(fit%nllf),'all-free finite candidate')
        call check(fit%nllf <= -ll(data,family,par),'all-free candidate improves')
    end do
    one = make_censored_data(uncensored=[0.7_dp],left=[0.5_dp],right=[2.0_dp], &
        interval=reshape([0.8_dp,1.3_dp],[1,2]))
    call close(gengamma_censored_loglikelihood(one,1.7_dp,-1.0_dp), &
        invgamma_censored_loglikelihood(one,1.7_dp),2e-12_dp,'inverse gamma identity')
    call close(gengamma_censored_loglikelihood(one,1.0_dp,1.3_dp), &
        weibull_censored_loglikelihood(one,1.3_dp),2e-12_dp,'Weibull identity')
    call check(ieee_is_nan(gengamma_censored_loglikelihood(empty,0.0_dp,1.0_dp)),'gamma zero shape')
    call check(ieee_is_nan(gengamma_censored_loglikelihood(empty,1.0_dp,0.0_dp)),'zero power')
    call check(ieee_is_nan(gengamma_censored_loglikelihood(empty,1.0_dp,1e7_dp)),'power range')
    call check(ieee_is_nan(truncnorm_censored_loglikelihood(empty,1.0_dp,1.0_dp)),'collapsed normal')
    call check(ieee_is_nan(truncnorm_censored_loglikelihood(empty,-101.0_dp,1.0_dp)),'normal range')
    call check(ieee_is_nan(truncweibull_min_censored_loglikelihood(empty,1.0_dp,-1.0_dp,2.0_dp)), &
        'negative Weibull truncation')
    one = make_censored_data(left=[0.0_dp],right=[3.0_dp])
    call close(truncnorm_censored_loglikelihood(one,0.0_dp,3.0_dp),-inf,0.0_dp,'normal endpoints')
    one = make_censored_data(interval=reshape([-1.0_dp,4.0_dp],[1,2]))
    call close(truncnorm_censored_loglikelihood(one,0.0_dp,3.0_dp),0.0_dp,0.0_dp,'normal support')
    call close(truncweibull_min_censored_loglikelihood(one,1.3_dp,0.1_dp,3.0_dp), &
        0.0_dp,0.0_dp,'Weibull support')
    one = make_censored_data(uncensored=[0.0_dp])
    call close(truncweibull_min_censored_loglikelihood(one,0.5_dp,0.0_dp,3.0_dp), &
        inf,0.0_dp,'Weibull infinite endpoint density')
    call close(truncweibull_min_censored_loglikelihood(one,2.0_dp,0.0_dp,3.0_dp), &
        -inf,0.0_dp,'Weibull zero endpoint density')
    one = make_censored_data(uncensored=[1.0_dp])
    expected = -0.5_dp-0.5_dp*log(2.0_dp*acos(-1.0_dp))
    call close(truncnorm_censored_loglikelihood(one,-100.0_dp,100.0_dp), &
        expected,2e-13_dp,'near-untruncated normal')
    do family = 2, 3
        one = make_censored_data(interval=reshape([1.0_dp,nearest(1.0_dp,1.0_dp)],[1,2]))
        x = nearest(nearest(1.0_dp,1.0_dp),1.0_dp)
        if (family == 2) then
            actual = truncnorm_censored_loglikelihood(one,1.0_dp,x)
        else
            actual = truncweibull_min_censored_loglikelihood(one,2.0_dp,1.0_dp,x)
        end if
        call close(actual,-log(2.0_dp),2e-12_dp,'two-ULP window half mass')
    end do
    call check(ieee_is_nan(truncnorm_censored_loglikelihood(empty,1.0_dp,x,loc=0.1_dp)), &
        'narrow translated window unsupported')
    call check(ieee_is_nan(truncnorm_censored_loglikelihood(empty,1.3_dp,3.0_dp, &
        scale=nearest(0.0_dp,1.0_dp))),'subnormal affine normal bounds rejected')
    call check(ieee_is_nan(truncweibull_min_censored_loglikelihood(empty,1.3_dp,1.3_dp,3.0_dp, &
        scale=nearest(0.0_dp,1.0_dp))),'subnormal affine Weibull bounds rejected')
    call check(ieee_is_nan(truncnorm_censored_loglikelihood(empty,-1.3_dp,3.0_dp,loc=1e10_dp)), &
        'poorly resolved normal support rejected')
    call check(ieee_is_nan(truncweibull_min_censored_loglikelihood(empty,1.3_dp,0.1_dp,3.0_dp, &
        loc=1e10_dp)),'poorly resolved Weibull support rejected')
    x = nearest(1e-100_dp,-1.0_dp)
    one = make_censored_data(uncensored=[x])
    expected = log(1e6_dp)-log(x)+1e6_dp*log(1.0_dp+(x-1e-100_dp)/1e-100_dp)
    call close(truncweibull_min_censored_loglikelihood(one,1e6_dp,0.0_dp,1e-100_dp), &
        expected,2e-12_dp,'tiny normalizer cancels in log density')
    one = make_censored_data(left=[0.7_dp])
    expected = -tiny(1.0_dp)*(-1e6_dp*log(0.7_dp)-0.5772156649015328606_dp)
    actual = gengamma_censored_loglikelihood(one,tiny(1.0_dp),1e6_dp)
    call close(actual/expected,1.0_dp,1e-12_dp,'minimum gamma shape small log-CDF relative')
    one = make_censored_data(uncensored=[1.0_dp])
    call close(gengamma_censored_loglikelihood(one,1.0_dp,1e-6_dp), &
        log(1e-6_dp)-1.0_dp,2e-13_dp,'minimum power density')
    call close(gengamma_censored_loglikelihood(one,1.0_dp,-1e6_dp), &
        log(1e6_dp)-1.0_dp,2e-13_dp,'maximum negative power density')
    call close(gengamma_censored_nnlf(data,1.7_dp,-1.3_dp), &
        -gengamma_censored_loglikelihood(data,1.7_dp,-1.3_dp),0.0_dp,'gamma NNLF')
    call close(truncnorm_censored_nnlf(data,-1.0_dp,3.0_dp), &
        -truncnorm_censored_loglikelihood(data,-1.0_dp,3.0_dp),0.0_dp,'normal NNLF')
    call close(truncweibull_min_censored_nnlf(data,1.3_dp,0.1_dp,3.0_dp), &
        -truncweibull_min_censored_loglikelihood(data,1.3_dp,0.1_dp,3.0_dp),0.0_dp,'Weibull NNLF')
    write(*,'(a,i0,a,i0)') 'power/truncation censoring checks: ',checks,', failures: ',failures
    if (failures /= 0) error stop 1
contains
    subroutine parameters(family, par)
        implicit none
        integer, intent(in) :: family !! generalized gamma=1, normal=2, Weibull=3
        real(dp), intent(out) :: par(5) !! shapes followed by loc and scale
        select case (family)
        case (1)
            par = [1.7_dp,1.3_dp,0.0_dp,1.0_dp,0.0_dp]
        case (2)
            par = [-1.0_dp,3.0_dp,0.0_dp,1.0_dp,0.0_dp]
        case (3)
            par = [1.3_dp,0.1_dp,3.0_dp,0.0_dp,1.0_dp]
        end select
    end subroutine parameters
    pure function ll(data, family, par) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! observations
        integer, intent(in) :: family !! generalized gamma=1, normal=2, Weibull=3
        real(dp), intent(in) :: par(5) !! family shapes followed by loc and scale
        real(dp) :: value
        select case (family)
        case (1)
            value = gengamma_censored_loglikelihood(data,par(1),par(2),par(3),par(4))
        case (2)
            value = truncnorm_censored_loglikelihood(data,par(1),par(2),par(3),par(4))
        case (3)
            value = truncweibull_min_censored_loglikelihood(data,par(1),par(2),par(3),par(4),par(5))
        end select
    end function ll

    pure function scalar(x, family, par, part) result(value)
        implicit none
        real(dp), intent(in) :: x !! standardized observation or limit
        integer, intent(in) :: part !! density=1, CDF=2, SF=3
        integer, intent(in) :: family !! generalized gamma=1, normal=2, Weibull=3
        real(dp), intent(in) :: par(5) !! family shapes followed by loc and scale
        real(dp) :: value
        select case (family)
        case (1)
            select case (part)
            case (1)
                value = gengamma_logpdf(x,par(1),par(2))
            case (2)
                value = gengamma_logcdf(x,par(1),par(2))
            case (3)
                value = gengamma_logsf(x,par(1),par(2))
            end select
        case (2)
            select case (part)
            case (1)
                value = truncnorm_logpdf(x,par(1),par(2))
            case (2)
                value = truncnorm_logcdf(x,par(1),par(2))
            case (3)
                value = truncnorm_logsf(x,par(1),par(2))
            end select
        case (3)
            select case (part)
            case (1)
                value = truncweibull_min_logpdf(x,par(1),par(2),par(3))
            case (2)
                value = truncweibull_min_logcdf(x,par(1),par(2),par(3))
            case (3)
                value = truncweibull_min_logsf(x,par(1),par(2),par(3))
            end select
        end select
    end function scalar

    pure function quantile(p, family, par, upper) result(value)
        implicit none
        real(dp), intent(in) :: p !! probability in (0,1)
        logical, intent(in) :: upper !! true for survival quantile
        integer, intent(in) :: family !! generalized gamma=1, normal=2, Weibull=3
        real(dp), intent(in) :: par(5) !! family shapes followed by loc and scale
        real(dp) :: value
        select case (family)
        case (1)
            if (upper) then
                value = gengamma_isf(p,par(1),par(2))
            else
                value = gengamma_ppf(p,par(1),par(2))
            end if
        case (2)
            if (upper) then
                value = truncnorm_isf(p,par(1),par(2))
            else
                value = truncnorm_ppf(p,par(1),par(2))
            end if
        case (3)
            if (upper) then
                value = truncweibull_min_isf(p,par(1),par(2),par(3))
            else
                value = truncweibull_min_ppf(p,par(1),par(2),par(3))
            end if
        end select
    end function quantile

    subroutine fit_family(data, family, lower, upper, result)
        implicit none
        type(censored_data), intent(in) :: data !! observations
        real(dp), intent(in) :: lower(5), upper(5) !! shape, loc, scale bounds; unused tail ignored
        type(fit_result), intent(out) :: result !! fit and diagnostics
        integer, intent(in) :: family !! generalized gamma=1, normal=2, Weibull=3
        select case (family)
        case (1)
            call gengamma_fit(data,lower(:4),upper(:4),result)
        case (2)
            call truncnorm_fit(data,lower(:4),upper(:4),result)
        case (3)
            call truncweibull_min_fit(data,lower(:5),upper(:5),result)
        end select
    end subroutine fit_family

    subroutine fit_array(data, family, lower, upper, result)
        implicit none
        real(dp), intent(in) :: data(:) !! observations
        real(dp), intent(in) :: lower(5), upper(5) !! shape, loc, scale bounds; unused tail ignored
        type(fit_result), intent(out) :: result !! fit and diagnostics
        integer, intent(in) :: family !! generalized gamma=1, normal=2, Weibull=3
        select case (family)
        case (1)
            call gengamma_fit(data,lower(:4),upper(:4),result)
        case (2)
            call truncnorm_fit(data,lower(:4),upper(:4),result)
        case (3)
            call truncweibull_min_fit(data,lower(:5),upper(:5),result)
        end select
    end subroutine fit_array

    subroutine check(ok, label)
        implicit none
        logical, intent(in) :: ok !! true when the assertion passes
        character(*), intent(in) :: label !! assertion description
        checks = checks+1
        if (ok) return
        failures = failures+1
        write(*,*) 'FAIL: ',label
    end subroutine check

    subroutine close(actual, expected, tol, label)
        implicit none
        real(dp), intent(in) :: actual !! observed result
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tol !! scaled tolerance with unit absolute floor
        character(*), intent(in) :: label !! assertion description
        logical :: ok
        ok = actual == expected
        if (ieee_is_finite(actual) .and. ieee_is_finite(expected)) &
            ok = abs(actual-expected) <= tol*max(1.0_dp,abs(expected))
        call check(ok,label)
        if (.not. ok) write(*,*) actual,expected
    end subroutine close
end program test_censored_power_trunc
