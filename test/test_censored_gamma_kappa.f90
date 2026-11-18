! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_gamma_kappa
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        fit_status_invalid_input, normal_censored_loglikelihood, &
        genextreme_censored_loglikelihood, exponential_censored_loglikelihood, &
        uniform_censored_loglikelihood, logistic_censored_loglikelihood, &
        argus_censored_loglikelihood, argus_censored_nnlf, argus_fit, &
        argus_logpdf, argus_logcdf, argus_logsf, argus_ppf, argus_isf, &
        pearson3_censored_loglikelihood, pearson3_censored_nnlf, pearson3_fit, &
        pearson3_logpdf, pearson3_logcdf, pearson3_logsf, pearson3_ppf, pearson3_isf, &
        kappa4_censored_loglikelihood, kappa4_censored_nnlf, kappa4_fit, &
        kappa4_logpdf, kappa4_logcdf, kappa4_logsf, kappa4_ppf, kappa4_isf
    use scifort_math, only : positive_infinity, quiet_nan
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    implicit none
    type(censored_data) :: data, one, empty
    type(fit_result) :: fit, old
    real(dp) :: par(4), lower(4), upper(4), x, actual, expected, previous, p, inf, nan
    real(dp), parameter :: scipy_scale(3) = &
        [0.8676833205360297_dp,0.5688502460453513_dp,0.47467771562661354_dp]
    real(dp), parameter :: scipy_nllf(3) = &
        [13.678940903065282_dp,13.368238839719634_dp,10.672577136502522_dp]
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
        if (family == 3) n = 4
        do j = 1, 4
            call parameters(family,par)
            if (family == 1) par(1) = real(j,dp)/2.0_dp
            if (family == 2) par(1) = (-1.0_dp)**j*real(j,dp)/2.0_dp
            if (family == 3) par(1) = real(j-2,dp)/2.0_dp
            previous = -inf
            do i = 1, 25
                x = 0.01_dp+real(i,dp)/40.0_dp
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
        call close(fit%params(n),scipy_scale(family),3e-5_dp,'SciPy scale fit')
        call close(fit%nllf,scipy_nllf(family),2e-8_dp,'SciPy fitted likelihood')
        lower = [0.3_dp,-0.1_dp,0.9_dp,0.0_dp]
        upper = [5.0_dp,0.05_dp,1.3_dp,0.0_dp]
        if (family == 2) then
            lower(:3) = [0.2_dp,-0.1_dp,0.5_dp]
            upper(:3) = [2.0_dp,0.1_dp,1.5_dp]
        else if (family == 3) then
            lower = [-0.5_dp,-0.2_dp,-0.1_dp,0.5_dp]
            upper = [0.7_dp,0.3_dp,0.1_dp,1.5_dp]
        end if
        call fit_family(data,family,lower,upper,fit)
        call check(ieee_is_finite(fit%nllf),'all-free finite candidate')
        call check(fit%nllf <= -ll(data,family,par),'all-free candidate improves')
    end do
    call close(pearson3_censored_loglikelihood(data,0.0_dp), &
        normal_censored_loglikelihood(data),2e-12_dp,'Pearson normal identity')
    call close(pearson3_censored_loglikelihood(data,2.0_dp), &
        exponential_censored_loglikelihood(data,loc=-1.0_dp),2e-12_dp,'Pearson exponential identity')
    call close(kappa4_censored_loglikelihood(data,0.0_dp,0.2_dp), &
        genextreme_censored_loglikelihood(data,0.2_dp),2e-12_dp,'kappa GEV identity')
    call close(kappa4_censored_loglikelihood(data,1.0_dp,1.0_dp), &
        uniform_censored_loglikelihood(data),2e-12_dp,'kappa uniform identity')
    call close(kappa4_censored_loglikelihood(data,-1.0_dp,0.0_dp), &
        logistic_censored_loglikelihood(data),2e-12_dp,'kappa logistic identity')
    call close(kappa4_censored_loglikelihood(data,1.0_dp,0.0_dp), &
        exponential_censored_loglikelihood(data),2e-12_dp,'kappa exponential identity')
    x = 1e-100_dp
    one = make_censored_data(left=[x])
    expected = argus_logpdf(x,2.0_dp)+log(x)-log(2.0_dp)
    call close(argus_censored_loglikelihood(one,2.0_dp),expected,2e-13_dp,'ARGUS tiny CDF')
    one = make_censored_data(right=[x])
    call close(argus_censored_loglikelihood(one,2.0_dp)/(-exp(expected)), &
        1.0_dp,2e-12_dp,'ARGUS tiny log-SF relative')
    do i = 1, 9
        x = real(i,dp)/10.0_dp
        one = make_censored_data(left=[x])
        actual = pearson3_censored_loglikelihood(one,1.3_dp)
        one = make_censored_data(right=[-x])
        call close(pearson3_censored_loglikelihood(one,-1.3_dp),actual,2e-13_dp, &
            'Pearson reflection identity')
    end do
    one = make_censored_data(left=[1e300_dp])
    call close(pearson3_censored_loglikelihood(one,0.0_dp,scale=1e-300_dp), &
        0.0_dp,0.0_dp,'Pearson normal overflowed coordinate')
    one = make_censored_data(right=[1000.0_dp])
    call close(kappa4_censored_loglikelihood(one,-1.0_dp,0.0_dp), &
        -1000.0_dp,2e-13_dp,'kappa underflowed survival')
    one = make_censored_data(uncensored=[0.0_dp,1.0_dp])
    call close(argus_censored_loglikelihood(one,2.0_dp),-inf,0.0_dp,'ARGUS endpoint densities')
    one = make_censored_data(interval=reshape([-1.0_dp,2.0_dp],[1,2]))
    call close(argus_censored_loglikelihood(one,2.0_dp),0.0_dp,0.0_dp,'ARGUS support interval')
    call close(kappa4_censored_loglikelihood(one,1.0_dp,1.0_dp),0.0_dp,0.0_dp,'kappa support interval')
    one = make_censored_data(uncensored=[-1.0_dp])
    call close(pearson3_censored_loglikelihood(one,2.0_dp),0.0_dp,0.0_dp,'Pearson finite endpoint')
    one = make_censored_data(uncensored=[-0.5_dp])
    call close(pearson3_censored_loglikelihood(one,4.0_dp),inf,0.0_dp,'Pearson infinite endpoint')
    one = make_censored_data(uncensored=[0.0_dp])
    call close(kappa4_censored_loglikelihood(one,2.0_dp,1.0_dp),-inf,0.0_dp,'kappa outside support')
    one = make_censored_data(uncensored=[log(2.0_dp)])
    call close(kappa4_censored_loglikelihood(one,2.0_dp,0.0_dp),inf,0.0_dp,'kappa infinite endpoint')
    call check(ieee_is_nan(argus_censored_loglikelihood(empty,0.0_dp)),'ARGUS shape range')
    call check(ieee_is_nan(pearson3_censored_loglikelihood(empty,1e-4_dp)),'Pearson shape range')
    call check(ieee_is_nan(kappa4_censored_loglikelihood(empty,1e-8_dp,1.0_dp)),'kappa shape range')
    call check(ieee_is_nan(kappa4_censored_loglikelihood(empty,10.0_dp,10.0_dp)), &
        'kappa compressed support rejected')
    call check(ieee_is_nan(kappa4_censored_loglikelihood(empty,1e-6_dp,-1.0_dp)), &
        'kappa compressed outer endpoint rejected')
    one = make_censored_data(uncensored=[-1e308_dp])
    call close(kappa4_censored_loglikelihood(one,-1.0_dp,0.0_dp), &
        -1e308_dp,2e-13_dp,'kappa finite largest log density')
    x = 1.5e154_dp
    one = make_censored_data(interval=reshape([x,nearest(x,1.0_dp)],[1,2]))
    call close(pearson3_censored_loglikelihood(one,0.0_dp), &
        -x*(0.5_dp*x),2e-13_dp,'Pearson normal largest finite interval')
    one = make_censored_data(uncensored=[-1.0_dp])
    call close(kappa4_censored_loglikelihood(one,-1.0_dp,-1.0_dp), &
        0.0_dp,0.0_dp,'kappa finite lower endpoint')
    call close(kappa4_censored_loglikelihood(one,-2.0_dp,-1.0_dp), &
        inf,0.0_dp,'kappa singular lower endpoint')
    one = make_censored_data(uncensored=[0.5_dp])
    call close(kappa4_censored_loglikelihood(one,0.0_dp,2.0_dp), &
        inf,0.0_dp,'kappa singular upper endpoint')
    call close(argus_censored_nnlf(data,1.3_dp), &
        -argus_censored_loglikelihood(data,1.3_dp),0.0_dp,'argus NNLF')
    call close(pearson3_censored_nnlf(data,1.3_dp), &
        -pearson3_censored_loglikelihood(data,1.3_dp),0.0_dp,'pearson3 NNLF')
    call close(kappa4_censored_nnlf(data,-0.3_dp,0.2_dp), &
        -kappa4_censored_loglikelihood(data,-0.3_dp,0.2_dp),0.0_dp,'kappa4 NNLF')
    write(*,'(a,i0,a,i0)') 'gamma/kappa censoring checks: ',checks,', failures: ',failures
    if (failures /= 0) error stop 1
contains
    subroutine parameters(family, par)
        implicit none
        integer, intent(in) :: family !! ARGUS=1, Pearson=2, kappa=3
        real(dp), intent(out) :: par(4) !! shapes, loc, scale; unused tail ignored
        par = [2.0_dp,0.0_dp,1.0_dp,0.0_dp]
        if (family == 2) par(1) = 1.0_dp
        if (family == 3) par = [0.3_dp,0.2_dp,0.0_dp,1.0_dp]
    end subroutine parameters
    pure function ll(data, family, par) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! observations
        integer, intent(in) :: family !! ARGUS=1, Pearson=2, kappa=3
        real(dp), intent(in) :: par(4) !! shapes, loc, scale; unused tail ignored
        real(dp) :: value
        select case (family)
        case (1)
            value = argus_censored_loglikelihood(data,par(1),par(2),par(3))
        case (2)
            value = pearson3_censored_loglikelihood(data,par(1),par(2),par(3))
        case (3)
            value = kappa4_censored_loglikelihood(data,par(1),par(2),par(3),par(4))
        end select
    end function ll

    pure function scalar(x, family, par, part) result(value)
        implicit none
        real(dp), intent(in) :: x !! standardized observation
        integer, intent(in) :: part !! density=1, CDF=2, SF=3
        integer, intent(in) :: family !! ARGUS=1, Pearson=2, kappa=3
        real(dp), intent(in) :: par(4) !! shapes, loc, scale; unused tail ignored
        real(dp) :: value
        select case (family)
        case (1)
            select case (part)
            case (1)
                value = argus_logpdf(x,par(1))
            case (2)
                value = argus_logcdf(x,par(1))
            case (3)
                value = argus_logsf(x,par(1))
            end select
        case (2)
            select case (part)
            case (1)
                value = pearson3_logpdf(x,par(1))
            case (2)
                value = pearson3_logcdf(x,par(1))
            case (3)
                value = pearson3_logsf(x,par(1))
            end select
        case (3)
            select case (part)
            case (1)
                value = kappa4_logpdf(x,par(1),par(2))
            case (2)
                value = kappa4_logcdf(x,par(1),par(2))
            case (3)
                value = kappa4_logsf(x,par(1),par(2))
            end select
        end select
    end function scalar

    pure function quantile(p, family, par, upper) result(value)
        implicit none
        real(dp), intent(in) :: p !! probability in (0,1)
        logical, intent(in) :: upper !! true for survival quantile
        integer, intent(in) :: family !! ARGUS=1, Pearson=2, kappa=3
        real(dp), intent(in) :: par(4) !! shapes, loc, scale; unused tail ignored
        real(dp) :: value
        select case (family)
        case (1)
            if (upper) then
                value = argus_isf(p,par(1))
            else
                value = argus_ppf(p,par(1))
            end if
        case (2)
            if (upper) then
                value = pearson3_isf(p,par(1))
            else
                value = pearson3_ppf(p,par(1))
            end if
        case (3)
            if (upper) then
                value = kappa4_isf(p,par(1),par(2))
            else
                value = kappa4_ppf(p,par(1),par(2))
            end if
        end select
    end function quantile

    subroutine fit_family(data, family, lower, upper, result)
        implicit none
        type(censored_data), intent(in) :: data !! observations
        real(dp), intent(in) :: lower(4), upper(4) !! parameter bounds in shape, loc, scale order
        type(fit_result), intent(out) :: result !! fit and diagnostics
        integer, intent(in) :: family !! ARGUS=1, Pearson=2, kappa=3
        select case (family)
        case (1)
            call argus_fit(data,lower(:3),upper(:3),result)
        case (2)
            call pearson3_fit(data,lower(:3),upper(:3),result)
        case (3)
            call kappa4_fit(data,lower(:4),upper(:4),result)
        end select
    end subroutine fit_family

    subroutine fit_array(data, family, lower, upper, result)
        implicit none
        real(dp), intent(in) :: data(:) !! exact observations
        real(dp), intent(in) :: lower(4), upper(4) !! parameter bounds in shape, loc, scale order
        type(fit_result), intent(out) :: result !! fit and diagnostics
        integer, intent(in) :: family !! ARGUS=1, Pearson=2, kappa=3
        select case (family)
        case (1)
            call argus_fit(data,lower(:3),upper(:3),result)
        case (2)
            call pearson3_fit(data,lower(:3),upper(:3),result)
        case (3)
            call kappa4_fit(data,lower(:4),upper(:4),result)
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
end program test_censored_gamma_kappa
