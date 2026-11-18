! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_three_tail
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        fit_status_invalid_input, &
        kstwobign_censored_loglikelihood, kstwobign_censored_nnlf, kstwobign_fit, &
        kstwobign_logpdf, kstwobign_logcdf, kstwobign_logsf, kstwobign_ppf, kstwobign_isf, &
        wrapcauchy_censored_loglikelihood, wrapcauchy_censored_nnlf, wrapcauchy_fit, &
        wrapcauchy_logpdf, wrapcauchy_logcdf, wrapcauchy_logsf, wrapcauchy_ppf, wrapcauchy_isf, &
        rel_breitwigner_censored_loglikelihood, rel_breitwigner_censored_nnlf, rel_breitwigner_fit, &
        rel_breitwigner_logpdf, rel_breitwigner_logcdf, rel_breitwigner_logsf, rel_breitwigner_ppf, rel_breitwigner_isf
    use scifort_math, only : positive_infinity, quiet_nan
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    implicit none
    type(censored_data) :: data, one, empty
    type(fit_result) :: fit, old
    real(dp) :: par(6), lower(6), upper(6), x, actual, expected, previous, p, inf, nan
    integer :: checks, failures, family, i, j, k, n
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    empty = make_censored_data()
    data = make_censored_data(uncensored=[0.2_dp,0.4_dp,0.6_dp],left=[0.1_dp], &
        right=[0.8_dp],interval=reshape([0.3_dp,0.5_dp,0.35_dp,0.55_dp],[2,2]))
    do family = 1, 3
        n = 3
        if (family == 1) n = 2
        do j = 1, 4
            call parameters(family,par)
            if (family == 1) par(2) = 1.0_dp
            if (family == 2) par(1) = real(j,dp)/5.0_dp
            if (family == 3) par(1) = real(j,dp)/2.0_dp
            previous = -inf
            do i = 1, 25
                x = 0.05_dp+real(i,dp)/10.0_dp
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
            x = 0.3_dp*par(n)
            one = make_censored_data(uncensored=[x])
            expected = scalar(0.6_dp,family,par,1)-log(par(n))
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
        one = make_censored_data(uncensored=[0.2_dp,0.4_dp,0.6_dp])
        call fit_family(one,family,lower,upper,fit)
        call fit_array([0.2_dp,0.4_dp,0.6_dp],family,lower,upper,old)
        call check(fit%success .and. old%success,'old generic retained')
        call close(fit%nllf,old%nllf,3e-12_dp,'array/exact fit agreement')
        lower(n) = 0.3_dp
        upper(n) = 3.0_dp
        call fit_family(data,family,lower,upper,fit)
        call check(fit%success,'one-parameter fit converges')
        call check(fit%nllf <= -ll(data,family,par),'one-parameter fit improves')
        select case (family)
        case (1)
            x = 0.36696009974403854_dp
            expected = 28.498269735035024_dp
        case (2)
            x = 0.3000000072924005_dp
            expected = 12.550326038133834_dp
        case (3)
            x = 0.46747001313866393_dp
            expected = 9.354988140535767_dp
        end select
        call close(fit%params(n),x,3e-5_dp,'independent scale fit')
        call close(fit%nllf,expected,2e-8_dp,'independent fitted NNLF')
        lower = par-0.1_dp
        upper = par+0.1_dp
        select case (family)
        case (1)
            call kstwobign_fit(data,lower(:2),upper(:2),fit,guess=par(:2),max_iter=20)
        case (2)
            call wrapcauchy_fit(data,lower(:3),upper(:3),fit,guess=par(:3),max_iter=20)
        case (3)
            call rel_breitwigner_fit(data,lower(:3),upper(:3),fit,guess=par(:3),max_iter=20)
        end select
        call check(ieee_is_finite(fit%nllf),'all-free finite candidate')
        call check(fit%nllf <= -ll(data,family,par),'all-free improves initial point')
    end do
    one = make_censored_data(left=[0.0_dp])
    call close(kstwobign_censored_loglikelihood(one),-inf,0.0_dp,'Kolmogorov support CDF')
    call close(wrapcauchy_censored_loglikelihood(one,0.4_dp),-inf,0.0_dp,'angular support CDF')
    call close(rel_breitwigner_censored_loglikelihood(one,1.0_dp),-inf,0.0_dp,'resonance support CDF')
    one = make_censored_data(uncensored=[0.0_dp])
    call close(kstwobign_censored_loglikelihood(one),-inf,0.0_dp,'zero Kolmogorov density')
    call close(wrapcauchy_censored_loglikelihood(one,0.4_dp), &
        wrapcauchy_logpdf(0.0_dp,0.4_dp),3e-13_dp,'finite angular endpoint density')
    call close(rel_breitwigner_censored_loglikelihood(one,1.0_dp), &
        rel_breitwigner_logpdf(0.0_dp,1.0_dp),3e-13_dp,'finite resonance endpoint density')
    empty = make_censored_data()
    call check(ieee_is_nan(wrapcauchy_censored_loglikelihood(empty,0.0_dp)),'angular lower shape')
    call check(ieee_is_nan(wrapcauchy_censored_loglikelihood(empty,1.0_dp)),'angular upper shape')
    call check(ieee_is_nan(rel_breitwigner_censored_loglikelihood(empty,0.01_dp)),'resonance shape')
    one = make_censored_data(interval=reshape([6.0_dp,0.2_dp],[1,2]))
    call check(ieee_is_nan(wrapcauchy_censored_loglikelihood(one,0.4_dp)),'crossing arc not an interval')
    call close(kstwobign_censored_loglikelihood(data),-kstwobign_censored_nnlf(data), &
        0.0_dp,'Kolmogorov NNLF')
    call close(wrapcauchy_censored_loglikelihood(data,0.4_dp), &
        -wrapcauchy_censored_nnlf(data,0.4_dp),0.0_dp,'angular NNLF')
    call close(rel_breitwigner_censored_loglikelihood(data,1.0_dp), &
        -rel_breitwigner_censored_nnlf(data,1.0_dp),0.0_dp,'resonance NNLF')
    one = make_censored_data(left=[18.0_dp])
    expected = -2.0_dp*exp(-648.0_dp)
    call close(kstwobign_censored_loglikelihood(one)/expected,1.0_dp,3e-12_dp,'tiny Kolmogorov log-CDF')
    one = make_censored_data(right=[1e100_dp])
    expected = log(2.0_dp)-2.0e200_dp
    call close(kstwobign_censored_loglikelihood(one),expected,3e-13_dp,'large Kolmogorov SF')
    one = make_censored_data(left=[1e-150_dp])
    expected = 0.5_dp*log(2.0_dp*acos(-1.0_dp))-log(1e-150_dp) &
        -(acos(-1.0_dp)**2/8.0_dp)*1e300_dp
    call close(kstwobign_censored_loglikelihood(one),expected,3e-13_dp,'small Kolmogorov CDF')
    x = nearest(0.0_dp,1.0_dp)
    one = make_censored_data(interval=reshape([x,2.0_dp*x],[1,2]))
    expected = wrapcauchy_logpdf(0.0_dp,0.4_dp)+log(x)
    call close(wrapcauchy_censored_loglikelihood(one,0.4_dp),expected,3e-13_dp,'subnormal angular interval')
    expected = rel_breitwigner_logpdf(0.0_dp,1.0_dp)+log(x)
    call close(rel_breitwigner_censored_loglikelihood(one,1.0_dp),expected,3e-13_dp,'subnormal resonance interval')
    one = make_censored_data(right=[1e100_dp])
    call check(ieee_is_finite(rel_breitwigner_censored_loglikelihood(one,1.0_dp)), &
        'far resonance log-SF retained')
    one = make_censored_data(left=[acos(-1.0_dp)])
    call close(wrapcauchy_censored_loglikelihood(one,0.99_dp),-log(2.0_dp), &
        3e-13_dp,'angular branch midpoint')
    one = make_censored_data(interval=reshape([-1.0_dp,7.0_dp],[1,2]))
    call close(wrapcauchy_censored_loglikelihood(one,0.4_dp),0.0_dp,0.0_dp,'whole angular support')
    call check(ieee_is_nan(wrapcauchy_censored_loglikelihood(empty,0.999_dp)),'angular shape cap')
    call check(ieee_is_nan(wrapcauchy_censored_loglikelihood(empty,0.4_dp,loc=11.0_dp)), &
        'angular location cancellation guard')
    one = make_censored_data(left=[huge(1.0_dp)])
    call close(rel_breitwigner_censored_loglikelihood(one,1.0_dp, &
        loc=-huge(1.0_dp),scale=huge(1.0_dp)),rel_breitwigner_logcdf(2.0_dp,1.0_dp), &
        3e-12_dp,'overflowed physical resonance distance')
    x = 1e12_dp+1.0_dp
    one = make_censored_data(interval=reshape([x,nearest(x,1.0_dp)],[1,2]))
    call close(rel_breitwigner_censored_loglikelihood(one,1.0_dp,loc=1e12_dp), &
        -9.210035685388604_dp,3e-13_dp,'large-location narrow resonance interval')
    x = 1e10_dp+0.8_dp
    one = make_censored_data(interval=reshape([x,nearest(x,1.0_dp)],[1,2]))
    call close(kstwobign_censored_loglikelihood(one,loc=1e10_dp), &
        -12.683043840480472_dp,3e-13_dp,'large-location narrow Kolmogorov interval')
    write(*,'(a,i0,a,i0)') 'three-tail censoring: ',checks,' checks; failures: ',failures
    if (failures /= 0) error stop 1
contains
    subroutine parameters(family,par)
        implicit none
        integer, intent(in) :: family !! Kolmogorov=1, wrapped Cauchy=2, Breit-Wigner=3
        real(dp), intent(out) :: par(6) !! shapes, location, scale, with unused tail zero
        par = 0.0_dp
        select case (family)
        case (1)
            par(:2) = [0.0_dp,1.0_dp]
        case (2)
            par(:3) = [0.4_dp,0.0_dp,1.0_dp]
        case (3)
            par(:3) = [1.0_dp,0.0_dp,1.0_dp]
        end select
    end subroutine parameters

    pure function ll(data,family,par) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        integer, intent(in) :: family !! Kolmogorov=1, wrapped Cauchy=2, Breit-Wigner=3
        real(dp), intent(in) :: par(6) !! shapes, loc, scale; unused tail ignored
        real(dp) :: value
        select case (family)
        case (1)
            value = kstwobign_censored_loglikelihood(data,par(1),par(2))
        case (2)
            value = wrapcauchy_censored_loglikelihood(data,par(1),par(2),par(3))
        case (3)
            value = rel_breitwigner_censored_loglikelihood(data,par(1),par(2),par(3))
        end select
    end function ll

    pure function scalar(x,family,par,part) result(value)
        implicit none
        real(dp), intent(in) :: x !! standardized observation
        integer, intent(in) :: part !! density=1, CDF=2, SF=3
        integer, intent(in) :: family !! Kolmogorov=1, wrapped Cauchy=2, Breit-Wigner=3
        real(dp), intent(in) :: par(6) !! shapes, loc, scale; unused tail ignored
        real(dp) :: value
        select case (family)
        case (1)
            select case (part)
            case (1)
                value = kstwobign_logpdf(x)
            case (2)
                value = kstwobign_logcdf(x)
            case (3)
                value = kstwobign_logsf(x)
            end select
        case (2)
            select case (part)
            case (1)
                value = wrapcauchy_logpdf(x,par(1))
            case (2)
                value = wrapcauchy_logcdf(x,par(1))
            case (3)
                value = wrapcauchy_logsf(x,par(1))
            end select
        case (3)
            select case (part)
            case (1)
                value = rel_breitwigner_logpdf(x,par(1))
            case (2)
                value = rel_breitwigner_logcdf(x,par(1))
            case (3)
                value = rel_breitwigner_logsf(x,par(1))
            end select
        end select
    end function scalar

    pure function quantile(p,family,par,upper) result(value)
        implicit none
        real(dp), intent(in) :: p !! probability in (0,1)
        logical, intent(in) :: upper !! select survival quantile
        integer, intent(in) :: family !! Kolmogorov=1, wrapped Cauchy=2, Breit-Wigner=3
        real(dp), intent(in) :: par(6) !! shapes, loc, scale; unused tail ignored
        real(dp) :: value
        select case (family)
        case (1)
            if (upper) then
                value = kstwobign_isf(p)
            else
                value = kstwobign_ppf(p)
            end if
        case (2)
            if (upper) then
                value = wrapcauchy_isf(p,par(1))
            else
                value = wrapcauchy_ppf(p,par(1))
            end if
        case (3)
            if (upper) then
                value = rel_breitwigner_isf(p,par(1))
            else
                value = rel_breitwigner_ppf(p,par(1))
            end if
        end select
    end function quantile

    subroutine fit_family(data,family,lower,upper,result)
        implicit none
        type(censored_data), intent(in) :: data !! observations
        real(dp), intent(in) :: lower(6), upper(6) !! bounds in shape, loc, scale order
        type(fit_result), intent(out) :: result !! fit and diagnostics
        integer, intent(in) :: family !! Kolmogorov=1, wrapped Cauchy=2, Breit-Wigner=3
        select case (family)
        case (1)
            call kstwobign_fit(data,lower(:2),upper(:2),result)
        case (2)
            call wrapcauchy_fit(data,lower(:3),upper(:3),result)
        case (3)
            call rel_breitwigner_fit(data,lower(:3),upper(:3),result)
        end select
    end subroutine fit_family

    subroutine fit_array(data,family,lower,upper,result)
        implicit none
        real(dp), intent(in) :: data(:) !! exact observations
        real(dp), intent(in) :: lower(6), upper(6) !! bounds in shape, loc, scale order
        type(fit_result), intent(out) :: result !! fit and diagnostics
        integer, intent(in) :: family !! Kolmogorov=1, wrapped Cauchy=2, Breit-Wigner=3
        select case (family)
        case (1)
            call kstwobign_fit(data,lower(:2),upper(:2),result)
        case (2)
            call wrapcauchy_fit(data,lower(:3),upper(:3),result)
        case (3)
            call rel_breitwigner_fit(data,lower(:3),upper(:3),result)
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
end program test_censored_three_tail
