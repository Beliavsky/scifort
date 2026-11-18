! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_bessel
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        fit_status_invalid_input, &
        geninvgauss_censored_loglikelihood, geninvgauss_censored_nnlf, geninvgauss_fit, &
        geninvgauss_logpdf, geninvgauss_logcdf, geninvgauss_logsf, geninvgauss_ppf, &
        norminvgauss_censored_loglikelihood, norminvgauss_censored_nnlf, norminvgauss_fit, &
        norminvgauss_logpdf, norminvgauss_logcdf, norminvgauss_logsf, norminvgauss_ppf, &
        genhyperbolic_censored_loglikelihood, genhyperbolic_censored_nnlf, genhyperbolic_fit, &
        genhyperbolic_logpdf, genhyperbolic_logcdf, genhyperbolic_logsf, genhyperbolic_ppf
    use scifort_math, only : quiet_nan, positive_infinity
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
    implicit none
    type(censored_data) :: data, one, empty
    type(fit_result) :: fit, old
    real(dp) :: par(5), lower(5), upper(5), x, actual, expected, previous, inf, nan, q
    integer :: family, i, n, k, checks, failures
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    empty = make_censored_data()
    data = make_censored_data(uncensored=[0.2_dp,0.4_dp,0.6_dp],left=[0.1_dp], &
        right=[0.8_dp],interval=reshape([0.3_dp,0.5_dp,0.35_dp,0.55_dp],[2,2]))
    do family = 1, 3
        call parameters(family,par,n)
        previous = -inf
        do i = 1, 7
            x = real(i,dp)/2.0_dp
            if (family /= 1) x = x-2.0_dp
            one = make_censored_data(uncensored=[x])
            call close(ll(one,family,par),scalar(x,family,par,1),3e-11_dp,'density')
            one = make_censored_data(left=[x])
            actual = ll(one,family,par)
            call close(actual,scalar(x,family,par,2),3e-11_dp,'CDF')
            call check(actual >= previous .and. actual <= 0.0_dp,'CDF monotone')
            previous = actual
            one = make_censored_data(right=[x])
            call close(ll(one,family,par),scalar(x,family,par,3),3e-11_dp,'SF')
            one = make_censored_data(interval=reshape([x,nearest(x,1.0_dp)],[1,2]))
            expected = scalar(x,family,par,1)+log(nearest(x,1.0_dp)-x)
            call close(ll(one,family,par),expected,3e-11_dp,'adjacent interval')
            one = make_censored_data(interval=reshape([x,x+0.1_dp],[1,2]))
            expected = log(exp(scalar(x+0.1_dp,family,par,2))-exp(scalar(x,family,par,2)))
            call close(ll(one,family,par),expected,1e-9_dp,'ordinary interval')
        end do
        do i = 1, 1
            q = 0.5_dp
            x = quantile(q,family,par)
            one = make_censored_data(left=[x])
            call close(ll(one,family,par),log(q),3e-10_dp,'CDF/PPF')
        end do
        do k = -1, 1
            par(n) = 10.0_dp**(200*k)
            par(n-1) = -0.3_dp*par(n)
            x = 0.3_dp*par(n)
            one = make_censored_data(uncensored=[x])
            expected = scalar(0.6_dp,family,par,1)-log(par(n))
            call close(ll(one,family,par),expected,3e-11_dp,'affine PDF')
            one = make_censored_data(interval=reshape([x,nearest(x,1.0_dp)],[1,2]))
            call close(ll(one,family,par),expected+log(nearest(x,1.0_dp)-x),3e-11_dp,'affine interval')
        end do
        call parameters(family,par,n)
        do k = 1, 3
            par(n-1) = 10.0_dp**(4*k)
            x = par(n-1)+0.7_dp
            one = make_censored_data(interval=reshape([x,nearest(x,1.0_dp)],[1,2]))
            call check(ieee_is_finite(ll(one,family,par)),'large-location interval finite')
        end do
        call parameters(family,par,n)
        x = nearest(0.0_dp,1.0_dp)
        one = make_censored_data(interval=reshape([x,2.0_dp*x],[1,2]))
        if (family == 1) then
            call close(ll(one,family,par),-inf,0.0_dp,'unrepresentable subnormal log mass')
        else
            expected = scalar(0.0_dp,family,par,1)+log(x)
            call close(ll(one,family,par),expected,3e-11_dp,'subnormal interval')
        end if
        x = 1.0e9_dp
        one = make_censored_data(interval=reshape([x,nearest(x,1.0_dp)],[1,2]))
        call check(ieee_is_finite(ll(one,family,par)),'remote adjacent interval finite')
        call close(ll(empty,family,par),0.0_dp,0.0_dp,'valid empty')
        do k = 1, n
            lower = par
            lower(k) = nan
            call check(ieee_is_nan(ll(empty,family,lower)),'NaN parameter')
        end do
        lower = par
        lower(n) = 0.0_dp
        call check(ieee_is_nan(ll(empty,family,lower)),'invalid scale')
        lower(n) = inf
        call check(ieee_is_nan(ll(empty,family,lower)),'infinite scale')
        one = make_censored_data(left=[nan])
        call check(ieee_is_nan(ll(one,family,par)),'NaN data')
        one = make_censored_data(uncensored=[inf])
        call check(ieee_is_nan(ll(one,family,par)),'infinite exact invalid')
        one = make_censored_data(interval=reshape([1.0_dp,0.9_dp],[1,2]))
        call check(ieee_is_nan(ll(one,family,par)),'reversed interval')
        one = make_censored_data(interval=reshape([-inf,inf],[1,2]))
        call check(ieee_is_nan(ll(one,family,par)),'all infinite invalid')
        one = make_censored_data(interval=reshape([0.7_dp,0.7_dp],[1,2]))
        call close(ll(one,family,par),scalar(0.7_dp,family,par,1),3e-11_dp,'equal endpoint exact')
        one = make_censored_data(interval=reshape([-inf,0.7_dp],[1,2]))
        call close(ll(one,family,par),scalar(0.7_dp,family,par,2),3e-11_dp,'left interval')
        one = make_censored_data(interval=reshape([0.7_dp,inf],[1,2]))
        call close(ll(one,family,par),scalar(0.7_dp,family,par,3),3e-11_dp,'right interval')
        lower = par
        upper = par
        call fitter(data,family,lower,upper,fit)
        call check(fit%success,'fixed fit success')
        call close(fit%nllf,-ll(data,family,par),3e-12_dp,'fixed fit NNLF')
        call fitter(empty,family,lower,upper,fit)
        call check(fit%status == fit_status_invalid_input,'empty fit invalid')
        one = make_censored_data(uncensored=[0.2_dp,0.4_dp,0.6_dp])
        call fitter(one,family,lower,upper,fit)
        call array_fit([0.2_dp,0.4_dp,0.6_dp],family,lower,upper,old)
        call check(fit%success .and. old%success,'array generic retained')
        call close(fit%nllf,old%nllf,3e-11_dp,'array/exact agreement')
        if (family == 1) then
        lower(n) = 0.3_dp
        upper(n) = 3.0_dp
        call fitter(data,family,lower,upper,fit)
        call check(fit%success,'scale fit converges')
        call check(fit%nllf <= -ll(data,family,par),'scale fit improves')
        end if
    end do
    one = make_censored_data(left=[0.00644270263257127_dp])
    call close(geninvgauss_censored_loglikelihood(one,1.3957936504673851_dp, &
        6.9007004756020205_dp),-542.06273954714593454206824146_dp,3e-12_dp, &
        'independent concentrated lower tail')
    one = make_censored_data(left=[0.0_dp])
    call close(geninvgauss_censored_loglikelihood(one,0.7_dp,1.2_dp),-inf,0.0_dp,'GIG endpoint CDF')
    one = make_censored_data(right=[0.0_dp])
    call close(geninvgauss_censored_loglikelihood(one,0.7_dp,1.2_dp),0.0_dp,0.0_dp,'GIG endpoint SF')
    one = make_censored_data(uncensored=[-1.0_dp])
    call close(geninvgauss_censored_loglikelihood(one,0.7_dp,1.2_dp),-inf,0.0_dp,'GIG support PDF')
    call check(ieee_is_nan(geninvgauss_censored_loglikelihood(empty,6.0_dp,1.0_dp)),'order limit')
    call check(ieee_is_nan(norminvgauss_censored_loglikelihood(empty,0.4_dp,0.0_dp)),'tail limit')
    call check(ieee_is_nan(genhyperbolic_censored_loglikelihood(empty,0.0_dp,1.0_dp,0.9_dp)),'skew limit')
    call close(norminvgauss_censored_loglikelihood(data,1.5_dp,0.3_dp), &
        genhyperbolic_censored_loglikelihood(data,-0.5_dp,1.5_dp,0.3_dp),0.0_dp,'NIG/GH identity')
    call close(geninvgauss_censored_loglikelihood(data,0.7_dp,1.2_dp), &
        -geninvgauss_censored_nnlf(data,0.7_dp,1.2_dp),0.0_dp,'GIG NNLF')
    call close(norminvgauss_censored_loglikelihood(data,1.5_dp,0.3_dp), &
        -norminvgauss_censored_nnlf(data,1.5_dp,0.3_dp),0.0_dp,'NIG NNLF')
    call close(genhyperbolic_censored_loglikelihood(data,0.7_dp,1.5_dp,0.3_dp), &
        -genhyperbolic_censored_nnlf(data,0.7_dp,1.5_dp,0.3_dp),0.0_dp,'GH NNLF')
    do i = 1, 15
        x = 10.0_dp**(real(i-8,dp))
        one = make_censored_data(left=[x])
        data = make_censored_data(right=[1.0_dp/x])
        call close(geninvgauss_censored_loglikelihood(one,0.7_dp,1.2_dp), &
            geninvgauss_censored_loglikelihood(data,-0.7_dp,1.2_dp),3e-11_dp,'GIG reciprocal')
        one = make_censored_data(left=[x])
        data = make_censored_data(right=[-x])
        call close(genhyperbolic_censored_loglikelihood(one,0.7_dp,1.5_dp,0.3_dp), &
            genhyperbolic_censored_loglikelihood(data,0.7_dp,1.5_dp,-0.3_dp),3e-11_dp,'GH reflection')
    end do
    write(*,'(a,i0,a,i0)') 'Censored Bessel checks: ',checks,', failures: ',failures
    if (failures /= 0) stop 1
contains
    subroutine parameters(family,par,n)
        implicit none
        integer, intent(in) :: family !! distribution selector
        real(dp), intent(out) :: par(5) !! shapes, location, scale
        integer, intent(out) :: n !! parameter count
        par = 0.0_dp
        n = 4
        select case (family)
        case (1)
            par(:4) = [0.7_dp,1.2_dp,0.0_dp,1.0_dp]
        case (2)
            par(:4) = [1.5_dp,0.3_dp,0.0_dp,1.0_dp]
        case (3)
            n = 5
            par = [0.7_dp,1.5_dp,0.3_dp,0.0_dp,1.0_dp]
        end select
    end subroutine parameters
    pure function ll(data,family,par) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! observations
        integer, intent(in) :: family !! distribution selector
        real(dp), intent(in) :: par(5) !! shapes, location, scale
        real(dp) :: value
        select case (family)
        case (1)
            value = geninvgauss_censored_loglikelihood(data,par(1),par(2),par(3),par(4))
        case (2)
            value = norminvgauss_censored_loglikelihood(data,par(1),par(2),par(3),par(4))
        case (3)
            value = genhyperbolic_censored_loglikelihood(data,par(1),par(2),par(3),par(4),par(5))
        end select
    end function ll

    pure function scalar(x,family,par,part) result(value)
        implicit none
        real(dp), intent(in) :: x !! standardized coordinate
        integer, intent(in) :: family, part !! family selector and PDF/CDF/SF selector
        real(dp), intent(in) :: par(5) !! parameters
        real(dp) :: value
        select case (family)
        case (1)
            select case (part)
            case (1)
                value = geninvgauss_logpdf(x,par(1),par(2))
            case (2)
                value = geninvgauss_logcdf(x,par(1),par(2))
            case (3)
                value = geninvgauss_logsf(x,par(1),par(2))
            end select
        case (2)
            select case (part)
            case (1)
                value = norminvgauss_logpdf(x,par(1),par(2))
            case (2)
                value = norminvgauss_logcdf(x,par(1),par(2))
            case (3)
                value = norminvgauss_logsf(x,par(1),par(2))
            end select
        case (3)
            select case (part)
            case (1)
                value = genhyperbolic_logpdf(x,par(1),par(2),par(3))
            case (2)
                value = genhyperbolic_logcdf(x,par(1),par(2),par(3))
            case (3)
                value = genhyperbolic_logsf(x,par(1),par(2),par(3))
            end select
        end select
    end function scalar

    pure function quantile(q,family,par) result(value)
        implicit none
        real(dp), intent(in) :: q !! probability in (0,1)
        integer, intent(in) :: family !! distribution selector
        real(dp), intent(in) :: par(5) !! parameters
        real(dp) :: value
        select case (family)
        case (1)
            value = geninvgauss_ppf(q,par(1),par(2))
        case (2)
            value = norminvgauss_ppf(q,par(1),par(2))
        case (3)
            value = genhyperbolic_ppf(q,par(1),par(2),par(3))
        end select
    end function quantile

    subroutine fitter(data,family,lower,upper,result)
        implicit none
        type(censored_data), intent(in) :: data !! observations
        integer, intent(in) :: family !! distribution selector
        real(dp), intent(in) :: lower(5), upper(5) !! bounds in shapes, loc, scale order
        type(fit_result), intent(out) :: result !! fit diagnostics
        select case (family)
        case (1)
            call geninvgauss_fit(data,lower(:4),upper(:4),result)
        case (2)
            call norminvgauss_fit(data,lower(:4),upper(:4),result)
        case (3)
            call genhyperbolic_fit(data,lower(:5),upper(:5),result)
        end select
    end subroutine fitter

    subroutine array_fit(data,family,lower,upper,result)
        implicit none
        real(dp), intent(in) :: data(:) !! observations
        integer, intent(in) :: family !! distribution selector
        real(dp), intent(in) :: lower(5), upper(5) !! bounds in shapes, loc, scale order
        type(fit_result), intent(out) :: result !! fit diagnostics
        select case (family)
        case (1)
            call geninvgauss_fit(data,lower(:4),upper(:4),result)
        case (2)
            call norminvgauss_fit(data,lower(:4),upper(:4),result)
        case (3)
            call genhyperbolic_fit(data,lower(:5),upper(:5),result)
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
end program test_censored_bessel
