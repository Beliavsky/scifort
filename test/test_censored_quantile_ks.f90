! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_quantile_ks
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        fit_status_invalid_input, &
        tukeylambda_censored_loglikelihood, tukeylambda_censored_nnlf, tukeylambda_fit, &
        tukeylambda_logpdf, tukeylambda_logcdf, tukeylambda_logsf, tukeylambda_ppf, &
        ksone_censored_loglikelihood, ksone_censored_nnlf, ksone_fit, &
        ksone_logpdf, ksone_logcdf, ksone_logsf, ksone_ppf, &
        vonmises_censored_loglikelihood, vonmises_censored_nnlf, vonmises_fit, &
        vonmises_line_logpdf, vonmises_line_logcdf, vonmises_line_logsf, vonmises_line_ppf, &
        logistic_censored_loglikelihood, uniform_censored_loglikelihood, &
        vonmises_line_censored_loglikelihood
    use scifort_math, only : quiet_nan, positive_infinity
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
    implicit none
    type(censored_data) :: data, one, empty
    type(fit_result) :: fit, old
    real(dp) :: par(3), lower(3), upper(3), x, expected, actual, previous, inf, nan, q
    real(dp) :: best, candidate, scales(3)
    integer :: family, i, j, k, checks, failures
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    empty = make_censored_data()
    data = make_censored_data(uncensored=[0.2_dp,0.4_dp,0.6_dp],left=[0.1_dp], &
        right=[0.8_dp],interval=reshape([0.3_dp,0.5_dp,0.35_dp,0.55_dp],[2,2]))
    do family = 1, 3
        call parameters(family,par)
        previous = -inf
        do i = 1, 19
            x = real(i,dp)/20.0_dp
            if (family /= 2) x = 4.0_dp*x-2.0_dp
            one = make_censored_data(uncensored=[x])
            call close(ll(one,family,par),scalar(x,family,par,1),3e-11_dp,'density')
            one = make_censored_data(left=[x])
            actual = ll(one,family,par)
            call close(actual,scalar(x,family,par,2),3e-11_dp,'CDF')
            call check(actual >= previous .and. actual <= 0.0_dp,'CDF monotone')
            previous = actual
            one = make_censored_data(right=[x])
            call close(ll(one,family,par),scalar(x,family,par,3),3e-11_dp,'SF')
            ! Avoid the defined midpoint density at the first KSone knot.
            if (family /= 2 .or. abs(x-1.0_dp/par(1)) > 1e-12_dp) then
                one = make_censored_data(interval=reshape([x,nearest(x,1.0_dp)],[1,2]))
                expected = scalar(x,family,par,1)+log(nearest(x,1.0_dp)-x)
                call close(ll(one,family,par),expected,3e-11_dp,'adjacent interval')
            end if
            one = make_censored_data(interval=reshape([x,x+0.01_dp],[1,2]))
            if (scalar(x,family,par,2) < log(0.5_dp)) then
                expected = log(exp(scalar(x+0.01_dp,family,par,2))-exp(scalar(x,family,par,2)))
            else
                expected = log(exp(scalar(x,family,par,3))-exp(scalar(x+0.01_dp,family,par,3)))
            end if
            call close(ll(one,family,par),expected,2e-9_dp,'ordinary interval')
        end do
        do i = 1, 3
            q = real(i,dp)/4.0_dp
            x = quantile(q,family,par)
            one = make_censored_data(left=[x])
            call close(ll(one,family,par),log(q),3e-10_dp,'CDF/PPF')
        end do
        do k = -1, 1
            par(3) = 10.0_dp**(200*k)
            par(2) = -0.3_dp*par(3)
            x = 0.3_dp*par(3)
            expected = scalar(0.6_dp,family,par,1)-log(par(3))
            one = make_censored_data(uncensored=[x])
            call close(ll(one,family,par),expected,3e-11_dp,'affine density')
            one = make_censored_data(interval=reshape([x,nearest(x,1.0_dp)],[1,2]))
            call close(ll(one,family,par),expected+log(nearest(x,1.0_dp)-x),3e-11_dp,'affine adjacent')
        end do
        call parameters(family,par)
        call close(ll(empty,family,par),0.0_dp,0.0_dp,'valid empty')
        do k = 1, 3
            lower = par
            lower(k) = nan
            call check(ieee_is_nan(ll(empty,family,lower)),'NaN parameter')
        end do
        lower = par
        lower(3) = 0.0_dp
        call check(ieee_is_nan(ll(empty,family,lower)),'zero scale')
        lower(3) = inf
        call check(ieee_is_nan(ll(empty,family,lower)),'infinite scale')
        one = make_censored_data(left=[nan])
        call check(ieee_is_nan(ll(one,family,par)),'NaN observation')
        one = make_censored_data(uncensored=[inf])
        call check(ieee_is_nan(ll(one,family,par)),'infinite exact invalid')
        one = make_censored_data(interval=reshape([0.7_dp,0.6_dp],[1,2]))
        call check(ieee_is_nan(ll(one,family,par)),'reversed interval')
        one = make_censored_data(interval=reshape([-inf,inf],[1,2]))
        call check(ieee_is_nan(ll(one,family,par)),'all infinite invalid')
        one = make_censored_data(interval=reshape([0.7_dp,0.7_dp],[1,2]))
        call close(ll(one,family,par),scalar(0.7_dp,family,par,1),3e-11_dp,'equal endpoints exact')
        one = make_censored_data(interval=reshape([-inf,0.7_dp],[1,2]))
        call close(ll(one,family,par),scalar(0.7_dp,family,par,2),3e-11_dp,'left interval')
        one = make_censored_data(interval=reshape([0.7_dp,inf],[1,2]))
        call close(ll(one,family,par),scalar(0.7_dp,family,par,3),3e-11_dp,'right interval')
        lower = par
        upper = par
        call fitter(data,family,lower,upper,fit)
        call check(fit%success,'fixed fit success')
        call close(fit%nllf,-ll(data,family,par),3e-12_dp,'fixed NNLF')
        call fitter(empty,family,lower,upper,fit)
        call check(fit%status == fit_status_invalid_input,'empty fit invalid')
        one = make_censored_data(uncensored=[0.2_dp,0.4_dp,0.6_dp])
        call fitter(one,family,lower,upper,fit)
        call array_fit([0.2_dp,0.4_dp,0.6_dp],family,lower,upper,old)
        call check(fit%success .and. old%success,'array generic retained')
        call close(fit%nllf,old%nllf,3e-11_dp,'array/exact fit agreement')
        if (family /= 2) then
            lower(3) = 0.3_dp
            upper(3) = 3.0_dp
            call fitter(data,family,lower,upper,fit)
            call check(fit%success,'scale fit converges')
            call check(fit%nllf <= -ll(data,family,par),'scale fit improves')
            if (family == 1) then
                call close(fit%params(3),0.4088112186967253_dp,3e-5_dp,'independent Tukey scale')
                call close(fit%nllf,12.307544101328782_dp,2e-8_dp,'independent Tukey fit NNLF')
            else
                call close(fit%params(3),0.491565248524064_dp,3e-5_dp,'independent branch scale')
                call close(fit%nllf,12.428763478977562_dp,2e-8_dp,'independent branch fit NNLF')
            end if
        else
            lower(1) = 3.0_dp
            upper(1) = 7.0_dp
            best = inf
            do j = 3, 7
                par(1) = real(j,dp)
                candidate = -ll(data,family,par)
                best = min(best,candidate)
            end do
            call fitter(data,family,lower,upper,fit)
            call check(fit%success,'integer-shape fit converges')
            call check(fit%params(1) == real(nint(fit%params(1)),dp),'integer shape retained')
            call close(fit%nllf,best,3e-11_dp,'integer fit matches enumeration')
        end if
    end do
    call check(ieee_is_nan(tukeylambda_censored_loglikelihood(empty,1e-8_dp)),'small shape limit')
    call check(ieee_is_nan(tukeylambda_censored_loglikelihood(empty,6.0_dp)),'large shape limit')
    call check(ieee_is_nan(ksone_censored_loglikelihood(empty,2.5_dp)),'noninteger shape')
    call check(ieee_is_nan(ksone_censored_loglikelihood(empty,101.0_dp)),'large sample limit')
    call check(ieee_is_nan(vonmises_censored_loglikelihood(empty,-1.0_dp)),'negative kappa')
    call check(ieee_is_nan(vonmises_censored_loglikelihood(empty,101.0_dp)),'large kappa')
    call close(tukeylambda_censored_loglikelihood(data,0.0_dp), &
        logistic_censored_loglikelihood(data),3e-11_dp,'logistic identity')
    call close(tukeylambda_censored_loglikelihood(data,1.0_dp), &
        uniform_censored_loglikelihood(data,loc=-1.0_dp,scale=2.0_dp),3e-11_dp,'uniform lambda one')
    call close(tukeylambda_censored_loglikelihood(data,2.0_dp,scale=2.0_dp), &
        uniform_censored_loglikelihood(data,loc=-1.0_dp,scale=2.0_dp),3e-11_dp,'uniform lambda two')
    call close(ksone_censored_loglikelihood(data,1.0_dp), &
        uniform_censored_loglikelihood(data),3e-11_dp,'uniform n one')
    call close(vonmises_censored_loglikelihood(data,100.0_dp), &
        vonmises_line_censored_loglikelihood(data,100.0_dp),0.0_dp,'fixed-branch line identity')
    call close(tukeylambda_censored_loglikelihood(data,0.14_dp), &
        -tukeylambda_censored_nnlf(data,0.14_dp),0.0_dp,'Tukey NNLF')
    call close(ksone_censored_loglikelihood(data,10.0_dp), &
        -ksone_censored_nnlf(data,10.0_dp),0.0_dp,'KS NNLF')
    call close(vonmises_censored_loglikelihood(data,1.0_dp), &
        -vonmises_censored_nnlf(data,1.0_dp),0.0_dp,'von Mises NNLF')
    one = make_censored_data(uncensored=[1.0_dp])
    call close(tukeylambda_censored_loglikelihood(one,1.0_dp),-log(2.0_dp),0.0_dp,'uniform endpoint PDF')
    call close(ksone_censored_loglikelihood(one,1.0_dp),0.0_dp,0.0_dp,'n one endpoint PDF')
    call close(ksone_censored_loglikelihood(one,2.0_dp),-inf,0.0_dp,'KS endpoint convention')
    one = make_censored_data(left=[0.0_dp])
    call close(ksone_censored_loglikelihood(one,10.0_dp),-inf,0.0_dp,'KS support CDF')
    one = make_censored_data(left=[400.0_dp])
    expected = -exp(-400.0_dp)
    call close(tukeylambda_censored_loglikelihood(one,0.0_dp)/expected,1.0_dp,3e-11_dp,'tiny logistic complement')
    do i = 1, 7
        q = real(i,dp)/8.0_dp
        x = q/10.0_dp
        one = make_censored_data(left=[x])
        call close(ksone_censored_loglikelihood(one,10.0_dp), &
            log(x)+9.0_dp*log(1.0_dp+x),3e-12_dp,'small KS polynomial')
    end do
    x = nearest(0.0_dp,1.0_dp)
    one = make_censored_data(interval=reshape([x,2.0_dp*x],[1,2]))
    call close(tukeylambda_censored_loglikelihood(one,-0.3_dp), &
        -2.3_dp*log(2.0_dp)+log(x),3e-11_dp,'subnormal Tukey interval')
    call close(ksone_censored_loglikelihood(one,100.0_dp),log(x),3e-11_dp,'subnormal KS interval')
    one = make_censored_data(right=[x])
    call close(ksone_censored_loglikelihood(one,100.0_dp)/(-x),1.0_dp,0.0_dp,'subnormal KS complement')
    one = make_censored_data(interval=reshape([x,2.0_dp*x],[1,2]))
    call close(vonmises_censored_loglikelihood(one,1.0_dp), &
        vonmises_line_logpdf(0.0_dp,1.0_dp)+log(x),3e-11_dp,'subnormal branch interval')
    x = 0.1_dp
    expected = log((1.0_dp+x)**8*(1.0_dp+10.0_dp*x)-1.0_dp)+log(nearest(x,1.0_dp)-x)
    one = make_censored_data(interval=reshape([x,nearest(x,1.0_dp)],[1,2]))
    call close(ksone_censored_loglikelihood(one,10.0_dp),expected,3e-11_dp,'right first-knot interval')
    x = nearest(0.1_dp,-1.0_dp)
    expected = 8.0_dp*log(1.1_dp)+log(2.0_dp)+log(0.1_dp-x)
    one = make_censored_data(interval=reshape([x,0.1_dp],[1,2]))
    call close(ksone_censored_loglikelihood(one,10.0_dp),expected,3e-11_dp,'left first-knot interval')
    scales = [1e-200_dp,1.0_dp,1e200_dp]
    do k = 1, 3
        x = scales(k)
        one = make_censored_data(interval=reshape([nearest(x,-1.0_dp),x],[1,2]))
        expected = 10.0_dp*(log(x-nearest(x,-1.0_dp))-log(x))
        call close(ksone_censored_loglikelihood(one,10.0_dp,scale=x),expected,3e-11_dp,'upper endpoint interval')
    end do
    write(*,'(a,i0,a,i0)') 'Quantile/KS/branch censored checks: ',checks,', failures: ',failures
    if (failures /= 0) stop 1
contains
    subroutine parameters(family,par)
        implicit none
        integer, intent(in) :: family !! distribution selector
        real(dp), intent(out) :: par(3) !! shape, location, scale
        par = [0.14_dp,0.0_dp,1.0_dp]
        if (family == 2) par(1) = 10.0_dp
        if (family == 3) par(1) = 1.0_dp
    end subroutine parameters
    pure function ll(data,family,par) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! observations
        integer, intent(in) :: family !! distribution selector
        real(dp), intent(in) :: par(3) !! shape, loc, scale
        real(dp) :: value
        select case (family)
        case (1)
            value = tukeylambda_censored_loglikelihood(data,par(1),par(2),par(3))
        case (2)
            value = ksone_censored_loglikelihood(data,par(1),par(2),par(3))
        case (3)
            value = vonmises_censored_loglikelihood(data,par(1),par(2),par(3))
        end select
    end function ll

    pure function scalar(x,family,par,part) result(value)
        implicit none
        real(dp), intent(in) :: x !! standardized coordinate
        integer, intent(in) :: family, part !! family and PDF/CDF/SF selectors
        real(dp), intent(in) :: par(3) !! parameters
        real(dp) :: value
        select case (family)
        case (1)
            select case (part)
            case (1)
                value = tukeylambda_logpdf(x,par(1))
            case (2)
                value = tukeylambda_logcdf(x,par(1))
            case (3)
                value = tukeylambda_logsf(x,par(1))
            end select
        case (2)
            select case (part)
            case (1)
                value = ksone_logpdf(x,par(1))
            case (2)
                value = ksone_logcdf(x,par(1))
            case (3)
                value = ksone_logsf(x,par(1))
            end select
        case (3)
            select case (part)
            case (1)
                value = vonmises_line_logpdf(x,par(1))
            case (2)
                value = vonmises_line_logcdf(x,par(1))
            case (3)
                value = vonmises_line_logsf(x,par(1))
            end select
        end select
    end function scalar

    pure function quantile(q,family,par) result(value)
        implicit none
        real(dp), intent(in) :: q !! probability in (0,1)
        integer, intent(in) :: family !! distribution selector
        real(dp), intent(in) :: par(3) !! parameters
        real(dp) :: value
        select case (family)
        case (1)
            value = tukeylambda_ppf(q,par(1))
        case (2)
            value = ksone_ppf(q,par(1))
        case (3)
            value = vonmises_line_ppf(q,par(1))
        end select
    end function quantile

    subroutine fitter(data,family,lower,upper,result)
        implicit none
        type(censored_data), intent(in) :: data !! observations
        integer, intent(in) :: family !! distribution selector
        real(dp), intent(in) :: lower(3), upper(3) !! shape, loc, scale bounds
        type(fit_result), intent(out) :: result !! fit diagnostics
        select case (family)
        case (1)
            call tukeylambda_fit(data,lower,upper,result)
        case (2)
            call ksone_fit(data,lower,upper,result)
        case (3)
            call vonmises_fit(data,lower,upper,result)
        end select
    end subroutine fitter

    subroutine array_fit(data,family,lower,upper,result)
        implicit none
        real(dp), intent(in) :: data(:) !! observations
        integer, intent(in) :: family !! distribution selector
        real(dp), intent(in) :: lower(3), upper(3) !! shape, loc, scale bounds
        type(fit_result), intent(out) :: result !! fit diagnostics
        select case (family)
        case (1)
            call tukeylambda_fit(data,lower,upper,result)
        case (2)
            call ksone_fit(data,lower,upper,result)
        case (3)
            call vonmises_fit(data,lower,upper,result)
        end select
    end subroutine array_fit

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
end program test_censored_quantile_ks
