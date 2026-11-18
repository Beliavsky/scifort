! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_shape
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        fit_status_invalid_input, gumbel_r_censored_loglikelihood, &
        alpha_censored_loglikelihood, alpha_censored_nnlf, alpha_fit, &
        alpha_logpdf, alpha_logcdf, alpha_logsf, alpha_ppf, alpha_isf, &
        genextreme_censored_loglikelihood, genextreme_censored_nnlf, genextreme_fit, &
        genextreme_logpdf, genextreme_logcdf, genextreme_logsf, genextreme_ppf, genextreme_isf, &
        crystalball_censored_loglikelihood, crystalball_censored_nnlf, crystalball_fit, &
        crystalball_logpdf, crystalball_logcdf, crystalball_logsf, crystalball_ppf, crystalball_isf
    use scifort_math, only : positive_infinity, quiet_nan
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    implicit none
    type(censored_data) :: data, one, empty
    type(fit_result) :: fit, old
    real(dp) :: par(4), lower(4), upper(4), x, actual, expected, previous, p, inf, nan
    real(dp), parameter :: scipy_scale(3) = &
        [1.2307260337339929_dp,0.90425515523373112_dp,1.1798418087757312_dp]
    real(dp), parameter :: scipy_nllf(3) = &
        [10.399301854809206_dp,13.013831218213928_dp,14.678207227692237_dp]
    integer :: checks, failures, family, i, j, k, n
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    empty = make_censored_data()
    data = make_censored_data(uncensored=[0.5_dp,0.8_dp,1.1_dp],left=[0.4_dp], &
        right=[1.7_dp],interval=reshape([0.6_dp,0.9_dp,0.7_dp,1.0_dp],[2,2]))
    do family = 1, 3
        n = 3
        if (family == 3) n = 4
        do j = 1, 4
            call parameters(family,par)
            if (family == 1) par(1) = real(j,dp)/2.0_dp
            if (family == 2) par(1) = real(j-2,dp)/10.0_dp
            if (family == 3) par(1) = real(j,dp)/2.0_dp
            previous = -inf
            do i = 1, 25
                x = 0.05_dp+real(i,dp)/10.0_dp
                if (family == 3) x = real(i,dp)/3.0_dp-6.0_dp
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
        lower(n) = 0.3_dp
        upper(n) = 3.0_dp
        call fit_family(data,family,lower,upper,fit)
        call check(fit%success,'one-parameter fit converges')
        call check(fit%nllf <= -ll(data,family,par),'one-parameter fit improves')
        call close(fit%params(n),scipy_scale(family),3e-5_dp,'SciPy scale fit')
        call close(fit%nllf,scipy_nllf(family),2e-8_dp,'SciPy fitted likelihood')
        lower = [0.3_dp,-0.1_dp,0.7_dp,0.0_dp]
        upper = [5.0_dp,0.1_dp,1.3_dp,0.0_dp]
        if (family == 2) then
            lower(1) = -0.3_dp
            upper(1) = 0.3_dp
        else if (family == 3) then
            lower = [0.3_dp,1.5_dp,-0.1_dp,0.7_dp]
            upper = [3.0_dp,4.0_dp,0.1_dp,1.3_dp]
        end if
        call fit_family(data,family,lower,upper,fit)
        call check(ieee_is_finite(fit%nllf),'all-free finite candidate')
        call check(fit%nllf <= -ll(data,family,par),'all-free candidate improves')
    end do
    call close(genextreme_censored_loglikelihood(data,0.0_dp), &
        gumbel_r_censored_loglikelihood(data),2e-12_dp,'Gumbel identity')
    one = make_censored_data(uncensored=[1.0_dp])
    call close(genextreme_censored_loglikelihood(one,1.0_dp),0.0_dp,0.0_dp,'GEV finite endpoint')
    one = make_censored_data(uncensored=[0.5_dp])
    call close(genextreme_censored_loglikelihood(one,2.0_dp),inf,0.0_dp,'GEV singular endpoint')
    one = make_censored_data(left=[-1.0_dp])
    call close(genextreme_censored_loglikelihood(one,-1.0_dp),-inf,0.0_dp,'GEV lower support')
    call close(alpha_censored_loglikelihood(one,1.0_dp),-inf,0.0_dp,'alpha lower support')
    one = make_censored_data(right=[1000.0_dp])
    call close(genextreme_censored_loglikelihood(one,0.0_dp),-1000.0_dp,1e-14_dp,'GEV tiny SF')
    one = make_censored_data(right=[1e300_dp])
    ! Phi(a) normalizer is fixed independently through the scalar alpha density at z=1.
    expected = alpha_logpdf(1.0_dp,1.0_dp)-0.5_dp-log(1e300_dp)
    call close(alpha_censored_loglikelihood(one,1.0_dp),expected,2e-13_dp,'alpha tiny SF')
    one = make_censored_data(left=[1e300_dp])
    call close(alpha_censored_loglikelihood(one,1.0_dp)/(-exp(expected)), &
        1.0_dp,2e-13_dp,'alpha tiny log-CDF relative')
    one = make_censored_data(interval=reshape([1e-100_dp,nearest(1e-100_dp,1.0_dp)],[1,2]))
    expected = -0.5_dp/nearest(1e-100_dp,1.0_dp)**2
    call close(alpha_censored_loglikelihood(one,1.0_dp),expected,3e-13_dp,'alpha steep interval')
    one = make_censored_data(left=[-1e-6_dp])
    x = nearest(1.0_dp,1.0_dp)
    expected = -sqrt(2.0_dp*acos(-1.0_dp))*0.5_dp*erfc(-1e-6_dp/sqrt(2.0_dp)) &
        *(x-1.0_dp)*1e-6_dp/x*exp(0.5e-12_dp)
    actual = crystalball_censored_loglikelihood(one,1e-6_dp,x)
    call close(actual/expected,1.0_dp,2e-13_dp,'Crystal Ball tiny log-CDF relative')
    do family = 1, 3
        call parameters(family,par)
        n = 3
        if (family == 3) n = 4
        par(n) = 1e-300_dp
        if (family == 2) par(1) = -1e-6_dp
        x = 1e300_dp
        if (family == 3) x = -x
        one = make_censored_data(uncensored=[x],left=[x],right=[x], &
            interval=reshape([x,nearest(x,1.0_dp)],[1,2]))
        call check(ieee_is_finite(ll(one,family,par)),'overflowed standardized coordinate')
    end do
    one = make_censored_data(right=[40.0_dp])
    call check(ieee_is_finite(crystalball_censored_loglikelihood(one,2.0_dp,3.0_dp)), &
        'crystalball extreme Gaussian SF')
    do i = 1, 50
        par = [0.1_dp+real(i,dp)/7.0_dp,2.7_dp,0.0_dp,1.0_dp]
        par(4) = 10.0_dp**(real(i-25,dp)*8.0_dp)
        par(3) = 0.2_dp*par(4)
        x = par(3)-par(1)*par(4)
        one = make_censored_data(interval= &
            reshape([nearest(x,-1.0_dp),nearest(x,1.0_dp)],[1,2]))
        expected = crystalball_logpdf(-par(1),par(1),par(2)) &
            -log(par(4))+log(nearest(x,1.0_dp)-nearest(x,-1.0_dp))
        call close(ll(one,3,par),expected,3e-12_dp,'adjacent affine transition')
    end do
    x = 1.5e154_dp
    one = make_censored_data(right=[x])
    expected = -x*(0.5_dp*x)
    call close(crystalball_censored_loglikelihood(one,2.0_dp,3.0_dp), &
        expected,1e-13_dp,'Crystal Ball largest finite log tails')
    one = make_censored_data(uncensored=[x])
    call close(crystalball_censored_loglikelihood(one,2.0_dp,3.0_dp), &
        expected,1e-13_dp,'Crystal Ball largest finite log density')
    one = make_censored_data(interval=reshape([x,nearest(x,1.0_dp)],[1,2]))
    call close(crystalball_censored_loglikelihood(one,2.0_dp,3.0_dp), &
        expected,1e-13_dp,'Crystal Ball largest finite interval')
    x = 1.0_dp/1.5e154_dp
    one = make_censored_data(left=[x])
    call close(alpha_censored_loglikelihood(one,1.0_dp),expected,3e-13_dp, &
        'alpha largest finite log tail')
    call check(ieee_is_nan(alpha_censored_loglikelihood(empty,0.0_dp)),'alpha shape range')
    call check(ieee_is_nan(genextreme_censored_loglikelihood(empty,1e-8_dp)),'GEV shape range')
    call check(ieee_is_nan(crystalball_censored_loglikelihood(empty,1.0_dp,1.0_dp)), &
        'crystalball exponent range')
    call close(alpha_censored_nnlf(data,1.3_dp), &
        -alpha_censored_loglikelihood(data,1.3_dp),0.0_dp,'alpha NNLF')
    call close(genextreme_censored_nnlf(data,1.3_dp), &
        -genextreme_censored_loglikelihood(data,1.3_dp),0.0_dp,'genextreme NNLF')
    call close(crystalball_censored_nnlf(data,1.3_dp,2.7_dp), &
        -crystalball_censored_loglikelihood(data,1.3_dp,2.7_dp),0.0_dp,'crystalball NNLF')
    write(*,'(a,i0,a,i0)') 'shape censoring checks: ',checks,', failures: ',failures
    if (failures /= 0) error stop 1
contains
    subroutine parameters(family, par)
        implicit none
        integer, intent(in) :: family !! alpha=1, GEV=2, Crystal Ball=3
        real(dp), intent(out) :: par(4) !! shapes, loc, scale; unused tail ignored
        par = [1.7_dp,0.0_dp,1.0_dp,0.0_dp]
        if (family == 2) par(1) = 0.2_dp
        if (family == 3) par = [1.5_dp,2.7_dp,0.0_dp,1.0_dp]
    end subroutine parameters
    pure function ll(data, family, par) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! observations
        integer, intent(in) :: family !! alpha=1, GEV=2, Crystal Ball=3
        real(dp), intent(in) :: par(4) !! shapes, loc, scale; unused tail ignored
        real(dp) :: value
        select case (family)
        case (1)
            value = alpha_censored_loglikelihood(data,par(1),par(2),par(3))
        case (2)
            value = genextreme_censored_loglikelihood(data,par(1),par(2),par(3))
        case (3)
            value = crystalball_censored_loglikelihood(data,par(1),par(2),par(3),par(4))
        end select
    end function ll

    pure function scalar(x, family, par, part) result(value)
        implicit none
        real(dp), intent(in) :: x !! standardized observation
        integer, intent(in) :: part !! density=1, CDF=2, SF=3
        integer, intent(in) :: family !! alpha=1, GEV=2, Crystal Ball=3
        real(dp), intent(in) :: par(4) !! shapes, loc, scale; unused tail ignored
        real(dp) :: value
        select case (family)
        case (1)
            select case (part)
            case (1)
                value = alpha_logpdf(x,par(1))
            case (2)
                value = alpha_logcdf(x,par(1))
            case (3)
                value = alpha_logsf(x,par(1))
            end select
        case (2)
            select case (part)
            case (1)
                value = genextreme_logpdf(x,par(1))
            case (2)
                value = genextreme_logcdf(x,par(1))
            case (3)
                value = genextreme_logsf(x,par(1))
            end select
        case (3)
            select case (part)
            case (1)
                value = crystalball_logpdf(x,par(1),par(2))
            case (2)
                value = crystalball_logcdf(x,par(1),par(2))
            case (3)
                value = crystalball_logsf(x,par(1),par(2))
            end select
        end select
    end function scalar

    pure function quantile(p, family, par, upper) result(value)
        implicit none
        real(dp), intent(in) :: p !! probability in (0,1)
        logical, intent(in) :: upper !! true for survival quantile
        integer, intent(in) :: family !! alpha=1, GEV=2, Crystal Ball=3
        real(dp), intent(in) :: par(4) !! shapes, loc, scale; unused tail ignored
        real(dp) :: value
        select case (family)
        case (1)
            if (upper) then
                value = alpha_isf(p,par(1))
            else
                value = alpha_ppf(p,par(1))
            end if
        case (2)
            if (upper) then
                value = genextreme_isf(p,par(1))
            else
                value = genextreme_ppf(p,par(1))
            end if
        case (3)
            if (upper) then
                value = crystalball_isf(p,par(1),par(2))
            else
                value = crystalball_ppf(p,par(1),par(2))
            end if
        end select
    end function quantile

    subroutine fit_family(data, family, lower, upper, result)
        implicit none
        type(censored_data), intent(in) :: data !! observations
        real(dp), intent(in) :: lower(4), upper(4) !! parameter bounds in shape, loc, scale order
        type(fit_result), intent(out) :: result !! fit and diagnostics
        integer, intent(in) :: family !! alpha=1, GEV=2, Crystal Ball=3
        select case (family)
        case (1)
            call alpha_fit(data,lower(:3),upper(:3),result)
        case (2)
            call genextreme_fit(data,lower(:3),upper(:3),result)
        case (3)
            call crystalball_fit(data,lower(:4),upper(:4),result)
        end select
    end subroutine fit_family

    subroutine fit_array(data, family, lower, upper, result)
        implicit none
        real(dp), intent(in) :: data(:) !! exact observations
        real(dp), intent(in) :: lower(4), upper(4) !! parameter bounds in shape, loc, scale order
        type(fit_result), intent(out) :: result !! fit and diagnostics
        integer, intent(in) :: family !! alpha=1, GEV=2, Crystal Ball=3
        select case (family)
        case (1)
            call alpha_fit(data,lower(:3),upper(:3),result)
        case (2)
            call genextreme_fit(data,lower(:3),upper(:3),result)
        case (3)
            call crystalball_fit(data,lower(:4),upper(:4),result)
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
end program test_censored_shape
