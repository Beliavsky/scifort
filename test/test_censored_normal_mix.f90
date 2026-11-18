! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_normal_mix
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        fit_status_invalid_input, normal_censored_loglikelihood, &
        exponnorm_censored_loglikelihood, exponnorm_censored_nnlf, exponnorm_fit, &
        exponnorm_logpdf, exponnorm_logcdf, exponnorm_logsf, exponnorm_ppf, exponnorm_isf, &
        skewnorm_censored_loglikelihood, skewnorm_censored_nnlf, skewnorm_fit, &
        skewnorm_logpdf, skewnorm_logcdf, skewnorm_logsf, skewnorm_ppf, skewnorm_isf, &
        dpareto_lognorm_censored_loglikelihood, dpareto_lognorm_censored_nnlf, dpareto_lognorm_fit, &
        dpareto_lognorm_logpdf, dpareto_lognorm_logcdf, dpareto_lognorm_logsf, dpareto_lognorm_ppf, dpareto_lognorm_isf
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
        if (family == 3) n = 6
        do j = 1, 4
            call parameters(family,par)
            if (family == 1) par(1) = real(j,dp)/2.0_dp
            if (family == 2) par(1) = real(j-2,dp)/2.0_dp
            if (family == 3) par(1) = real(j-2,dp)/5.0_dp
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
            x = 0.3000000072924005_dp
            expected = 10.081147658879386_dp
        case (2)
            x = 0.5118115365121739_dp
            expected = 10.105222635536382_dp
        case (3)
            x = 0.3000000072924005_dp
            expected = 10.008137885354882_dp
        end select
        call close(fit%params(n),x,3e-5_dp,'independent scale fit')
        call close(fit%nllf,expected,2e-8_dp,'independent fitted NNLF')
    end do
    call close(skewnorm_censored_loglikelihood(data,0.0_dp), &
        normal_censored_loglikelihood(data),2e-12_dp,'normal identity')
    one = make_censored_data(left=[0.0_dp])
    do j = -10, 10
        expected = log(0.5_dp-atan(real(j,dp))/acos(-1.0_dp))
        call close(skewnorm_censored_loglikelihood(one,real(j,dp)),expected,3e-13_dp, &
            'skew-normal origin identity')
    end do
    one = make_censored_data(right=[1000.0_dp])
    call close(exponnorm_censored_loglikelihood(one,1.0_dp),-999.5_dp,2e-13_dp,'exponential tail')
    one = make_censored_data(interval=reshape([1000.0_dp,nearest(1000.0_dp,1.0_dp)],[1,2]))
    expected = -999.5_dp+log(nearest(1000.0_dp,1.0_dp)-1000.0_dp)
    call close(exponnorm_censored_loglikelihood(one,1.0_dp),expected,2e-13_dp,'tiny exponential interval')
    one = make_censored_data(left=[-1000.0_dp])
    call check(ieee_is_finite(exponnorm_censored_loglikelihood(one,100.0_dp)),'tiny EMG CDF')
    call check(ieee_is_finite(skewnorm_censored_loglikelihood(one,10.0_dp)),'tiny skew CDF')
    one = make_censored_data(uncensored=[0.0_dp])
    call close(dpareto_lognorm_censored_loglikelihood(one,0.0_dp,1.0_dp,2.0_dp,3.0_dp), &
        -inf,0.0_dp,'positive support zero density')
    one = make_censored_data(left=[0.0_dp])
    call close(dpareto_lognorm_censored_loglikelihood(one,0.0_dp,1.0_dp,2.0_dp,3.0_dp), &
        -inf,0.0_dp,'positive support CDF')
    one = make_censored_data(right=[0.0_dp])
    call close(dpareto_lognorm_censored_loglikelihood(one,0.0_dp,1.0_dp,2.0_dp,3.0_dp), &
        0.0_dp,0.0_dp,'positive support SF')
    call check(ieee_is_nan(exponnorm_censored_loglikelihood(empty,0.001_dp)),'EMG shape range')
    call check(ieee_is_nan(skewnorm_censored_loglikelihood(empty,11.0_dp)),'skew shape range')
    call check(ieee_is_nan(dpareto_lognorm_censored_loglikelihood(empty,0.0_dp,0.001_dp, &
        2.0_dp,3.0_dp)),'lognormal shape range')
    one = make_censored_data(uncensored=[1.8e154_dp])
    call close(skewnorm_censored_loglikelihood(one,1.0_dp),-1.62e308_dp,3e-13_dp, &
        'large finite Gaussian log density')
    one = make_censored_data(uncensored=[1.9e154_dp])
    call close(skewnorm_censored_loglikelihood(one,1.0_dp),-inf,0.0_dp,'Gaussian log range')
    one = make_censored_data(left=[-huge(1.0_dp)])
    call close(skewnorm_censored_loglikelihood(one,-10.0_dp),-inf,0.0_dp,'large skew log range')
    call close(exponnorm_censored_loglikelihood(one,0.01_dp),-inf,0.0_dp,'large EMG log range')
    one = make_censored_data(right=[huge(1.0_dp)])
    call close(exponnorm_censored_loglikelihood(one,100.0_dp),-huge(1.0_dp)/100.0_dp, &
        3e-13_dp,'finite largest-magnitude exponential tail')
    one = make_censored_data(left=[-1e300_dp])
    call close(exponnorm_censored_loglikelihood(one,1.0_dp,scale=1e-300_dp),-inf, &
        0.0_dp,'overflowed standard coordinate')
    call close(skewnorm_censored_loglikelihood(one,1.0_dp,scale=1e-300_dp),-inf, &
        0.0_dp,'overflowed skew coordinate')
    call close(dpareto_lognorm_censored_loglikelihood(data,0.1_dp,0.7_dp,2.0_dp,3.0_dp), &
        dpareto_lognorm_censored_loglikelihood(data,0.1_dp+log(2.0_dp), &
            0.7_dp,2.0_dp,3.0_dp,scale=0.5_dp),3e-12_dp,'u and scale identity')
    one = make_censored_data(left=[0.4_dp])
    empty = make_censored_data(right=[2.5_dp])
    call close(dpareto_lognorm_censored_loglikelihood(one,0.1_dp,0.7_dp,2.0_dp,3.0_dp), &
        dpareto_lognorm_censored_loglikelihood(empty,-0.1_dp,0.7_dp,3.0_dp,2.0_dp), &
        3e-12_dp,'reciprocal tail identity')
    x = nearest(0.0_dp,1.0_dp)
    one = make_censored_data(left=[x])
    expected = log(0.4_dp)+0.5_dp*(3.0_dp*0.7_dp)**2+3.0_dp*(log(x)-0.1_dp)
    call close(dpareto_lognorm_censored_loglikelihood(one,0.1_dp,0.7_dp,2.0_dp,3.0_dp), &
        expected,3e-13_dp,'subnormal Pareto CDF')
    one = make_censored_data(interval=reshape([x,2.0_dp*x],[1,2]))
    expected = expected+3.0_dp*log(2.0_dp)+log(1.0_dp-0.125_dp)
    call close(dpareto_lognorm_censored_loglikelihood(one,0.1_dp,0.7_dp,2.0_dp,3.0_dp), &
        expected,3e-13_dp,'subnormal Pareto interval')
    one = make_censored_data(left=[huge(1.0_dp)])
    expected = dpareto_lognorm_logcdf(2.0_dp,0.1_dp,0.7_dp,2.0_dp,3.0_dp)
    call close(dpareto_lognorm_censored_loglikelihood(one,0.1_dp,0.7_dp,2.0_dp,3.0_dp, &
        loc=-huge(1.0_dp),scale=huge(1.0_dp)),expected,3e-12_dp,'overflowed positive distance')
    one = make_censored_data(left=[700.0_dp])
    expected = -exp(0.5_dp-700.0_dp)
    call close(exponnorm_censored_loglikelihood(one,1.0_dp)/expected,1.0_dp, &
        3e-11_dp,'tiny EMG log-CDF relative accuracy')
    one = make_censored_data(right=[-35.0_dp])
    expected = -3.1198801213181618e-270_dp
    call close(exponnorm_censored_loglikelihood(one,1.0_dp)/expected,1.0_dp, &
        3e-11_dp,'tiny EMG log-SF relative accuracy')
    one = make_censored_data(left=[35.0_dp])
    expected = -erfc(35.0_dp/sqrt(2.0_dp))
    call close(skewnorm_censored_loglikelihood(one,2.0_dp)/expected,1.0_dp, &
        3e-11_dp,'tiny skew log-CDF relative accuracy')
    one = make_censored_data(right=[-35.0_dp])
    call close(skewnorm_censored_loglikelihood(one,-2.0_dp)/expected,1.0_dp, &
        3e-11_dp,'tiny skew log-SF relative accuracy')
    one = make_censored_data(left=[1e100_dp])
    expected = -0.4_dp*exp(0.3_dp+0.5_dp*(3.0_dp*0.7_dp)**2)*1e-300_dp
    call close(dpareto_lognorm_censored_loglikelihood(one,0.1_dp,0.7_dp,3.0_dp,2.0_dp) &
        /expected,1.0_dp,3e-11_dp,'tiny Pareto log-CDF relative accuracy')
    one = make_censored_data(right=[1e-100_dp])
    expected = -0.4_dp*exp(-0.3_dp+0.5_dp*(3.0_dp*0.7_dp)**2)*1e-300_dp
    call close(dpareto_lognorm_censored_loglikelihood(one,0.1_dp,0.7_dp,2.0_dp,3.0_dp) &
        /expected,1.0_dp,3e-11_dp,'tiny Pareto log-SF relative accuracy')
    empty = make_censored_data()
    do family = 1, 3
        call parameters(family,par)
        n = 3
        if (family == 3) n = 6
        lower = par-0.1_dp
        upper = par+0.1_dp
        select case (family)
        case (1)
            call exponnorm_fit(data,lower(:3),upper(:3),fit,guess=par(:3),max_iter=20)
        case (2)
            call skewnorm_fit(data,lower(:3),upper(:3),fit,guess=par(:3),max_iter=20)
        case (3)
            call dpareto_lognorm_fit(data,lower,upper,fit,guess=par,max_iter=20)
        end select
        call check(ieee_is_finite(fit%nllf),'all-free finite candidate')
        call check(fit%nllf <= -ll(data,family,par),'all-free improves initial candidate')
        par(n-1) = inf
        call check(ieee_is_nan(ll(empty,family,par)),'infinite location')
    end do
    call close(exponnorm_censored_nnlf(data,1.3_dp), &
        -exponnorm_censored_loglikelihood(data,1.3_dp),0.0_dp,'exponnorm NNLF')
    call close(skewnorm_censored_nnlf(data,1.3_dp), &
        -skewnorm_censored_loglikelihood(data,1.3_dp),0.0_dp,'skewnorm NNLF')
    call close(dpareto_lognorm_censored_nnlf(data,0.1_dp,0.7_dp,2.0_dp,3.0_dp), &
        -dpareto_lognorm_censored_loglikelihood(data,0.1_dp,0.7_dp,2.0_dp,3.0_dp),0.0_dp,'dpareto_lognorm NNLF')
    write(*,'(a,i0,a,i0)') 'normal-mixture censoring checks: ',checks,', failures: ',failures
    if (failures /= 0) error stop 1
contains
    subroutine parameters(family, par)
        implicit none
        integer, intent(in) :: family !! EMG=1, skew normal=2, double Pareto-lognormal=3
        real(dp), intent(out) :: par(6) !! shapes, loc, scale; unused tail ignored
        par = [1.7_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp]
        if (family == 2) par(1) = 1.0_dp
        if (family == 3) par = [0.1_dp,0.7_dp,2.0_dp,3.0_dp,0.0_dp,1.0_dp]
    end subroutine parameters
    pure function ll(data, family, par) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! observations
        integer, intent(in) :: family !! EMG=1, skew normal=2, double Pareto-lognormal=3
        real(dp), intent(in) :: par(6) !! shapes, loc, scale; unused tail ignored
        real(dp) :: value
        select case (family)
        case (1)
            value = exponnorm_censored_loglikelihood(data,par(1),par(2),par(3))
        case (2)
            value = skewnorm_censored_loglikelihood(data,par(1),par(2),par(3))
        case (3)
            value = dpareto_lognorm_censored_loglikelihood(data,par(1),par(2),par(3),par(4),par(5),par(6))
        end select
    end function ll

    pure function scalar(x, family, par, part) result(value)
        implicit none
        real(dp), intent(in) :: x !! standardized observation
        integer, intent(in) :: part !! density=1, CDF=2, SF=3
        integer, intent(in) :: family !! EMG=1, skew normal=2, double Pareto-lognormal=3
        real(dp), intent(in) :: par(6) !! shapes, loc, scale; unused tail ignored
        real(dp) :: value
        select case (family)
        case (1)
            select case (part)
            case (1)
                value = exponnorm_logpdf(x,par(1))
            case (2)
                value = exponnorm_logcdf(x,par(1))
            case (3)
                value = exponnorm_logsf(x,par(1))
            end select
        case (2)
            select case (part)
            case (1)
                value = skewnorm_logpdf(x,par(1))
            case (2)
                value = skewnorm_logcdf(x,par(1))
            case (3)
                value = skewnorm_logsf(x,par(1))
            end select
        case (3)
            select case (part)
            case (1)
                value = dpareto_lognorm_logpdf(x,par(1),par(2),par(3),par(4))
            case (2)
                value = dpareto_lognorm_logcdf(x,par(1),par(2),par(3),par(4))
            case (3)
                value = dpareto_lognorm_logsf(x,par(1),par(2),par(3),par(4))
            end select
        end select
    end function scalar

    pure function quantile(p, family, par, upper) result(value)
        implicit none
        real(dp), intent(in) :: p !! probability in (0,1)
        logical, intent(in) :: upper !! true for survival quantile
        integer, intent(in) :: family !! EMG=1, skew normal=2, double Pareto-lognormal=3
        real(dp), intent(in) :: par(6) !! shapes, loc, scale; unused tail ignored
        real(dp) :: value
        select case (family)
        case (1)
            if (upper) then
                value = exponnorm_isf(p,par(1))
            else
                value = exponnorm_ppf(p,par(1))
            end if
        case (2)
            if (upper) then
                value = skewnorm_isf(p,par(1))
            else
                value = skewnorm_ppf(p,par(1))
            end if
        case (3)
            if (upper) then
                value = dpareto_lognorm_isf(p,par(1),par(2),par(3),par(4))
            else
                value = dpareto_lognorm_ppf(p,par(1),par(2),par(3),par(4))
            end if
        end select
    end function quantile

    subroutine fit_family(data, family, lower, upper, result)
        implicit none
        type(censored_data), intent(in) :: data !! observations
        real(dp), intent(in) :: lower(6), upper(6) !! parameter bounds in shape, loc, scale order
        type(fit_result), intent(out) :: result !! fit and diagnostics
        integer, intent(in) :: family !! EMG=1, skew normal=2, double Pareto-lognormal=3
        select case (family)
        case (1)
            call exponnorm_fit(data,lower(:3),upper(:3),result)
        case (2)
            call skewnorm_fit(data,lower(:3),upper(:3),result)
        case (3)
            call dpareto_lognorm_fit(data,lower(:6),upper(:6),result)
        end select
    end subroutine fit_family

    subroutine fit_array(data, family, lower, upper, result)
        implicit none
        real(dp), intent(in) :: data(:) !! exact observations
        real(dp), intent(in) :: lower(6), upper(6) !! parameter bounds in shape, loc, scale order
        type(fit_result), intent(out) :: result !! fit and diagnostics
        integer, intent(in) :: family !! EMG=1, skew normal=2, double Pareto-lognormal=3
        select case (family)
        case (1)
            call exponnorm_fit(data,lower(:3),upper(:3),result)
        case (2)
            call skewnorm_fit(data,lower(:3),upper(:3),result)
        case (3)
            call dpareto_lognorm_fit(data,lower(:6),upper(:6),result)
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
end program test_censored_normal_mix
