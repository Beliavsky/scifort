! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_bounded_final
    use scifort_gausshyper, only : gausshyper_logkernel_integral
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        fit_status_invalid_input, beta_censored_loglikelihood, uniform_censored_loglikelihood, &
        gausshyper_censored_loglikelihood, gausshyper_censored_nnlf, gausshyper_fit, &
        gausshyper_logpdf, gausshyper_logcdf, gausshyper_logsf, gausshyper_ppf, gausshyper_isf, &
        irwinhall_censored_loglikelihood, irwinhall_censored_nnlf, irwinhall_fit, &
        irwinhall_logpdf, irwinhall_logcdf, irwinhall_logsf, irwinhall_ppf, irwinhall_isf, &
        vonmises_line_censored_loglikelihood, vonmises_line_censored_nnlf, vonmises_line_fit, &
        vonmises_line_logpdf, vonmises_line_logcdf, vonmises_line_logsf, vonmises_line_ppf, vonmises_line_isf
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
        if (family == 1) n = 6
        do j = 1, 4
            call parameters(family,par)
            if (family == 1) par(1) = real(j,dp)/2.0_dp
            if (family == 2) par(1) = real(j,dp)
            if (family == 3) par(1) = real(j,dp)/2.0_dp
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
        select case (family)
        case (1)
            x = 1.2233448777774485_dp
            expected = 9.704537797653995_dp
        case (2)
            x = 0.33700479821178475_dp
            expected = 12.807592741295426_dp
        case (3)
            x = 0.7485579168664227_dp
            expected = 12.366078083348853_dp
        end select
        call close(fit%params(n),x,3e-5_dp,'independent fitted scale')
        call close(fit%nllf,expected,2e-8_dp,'independent fitted NNLF')
    end do
    call close(gausshyper_censored_loglikelihood(data,2.0_dp,3.0_dp,0.0_dp,1.0_dp), &
        beta_censored_loglikelihood(data,2.0_dp,3.0_dp),2e-12_dp,'beta identity')
    call close(irwinhall_censored_loglikelihood(data,1.0_dp), &
        uniform_censored_loglikelihood(data),2e-12_dp,'uniform sum identity')
    call close(vonmises_line_censored_loglikelihood(data,0.0_dp), &
        uniform_censored_loglikelihood(data,loc=-acos(-1.0_dp),scale=2*acos(-1.0_dp)), &
        2e-12_dp,'uniform line identity')
    one = make_censored_data(left=[1e-200_dp])
    expected = 32.0_dp*log(1e-200_dp)-log_gamma(33.0_dp)
    call close(irwinhall_censored_loglikelihood(one,32.0_dp),expected,3e-13_dp,'tiny spline CDF')
    one = make_censored_data(interval=reshape([1e-200_dp,2e-200_dp],[1,2]))
    expected = 32.0_dp*log(2e-200_dp)-log_gamma(33.0_dp)+log(1.0_dp-2.0_dp**(-32))
    call close(irwinhall_censored_loglikelihood(one,32.0_dp),expected,3e-13_dp,'tiny spline interval')
    x = nearest(0.0_dp,1.0_dp)
    one = make_censored_data(left=[x])
    expected = 32.0_dp*log(x)-log_gamma(33.0_dp)
    call close(irwinhall_censored_loglikelihood(one,32.0_dp),expected,3e-13_dp, &
        'smallest positive spline CDF')
    one = make_censored_data(interval=reshape([x,2.0_dp*x],[1,2]))
    expected = 32.0_dp*log(2.0_dp*x)-log_gamma(33.0_dp)+log(1.0_dp-2.0_dp**(-32))
    call close(irwinhall_censored_loglikelihood(one,32.0_dp),expected,3e-13_dp, &
        'subnormal spline interval')
    one = make_censored_data(uncensored=[0.0_dp,1.0_dp])
    call close(irwinhall_censored_loglikelihood(one,1.0_dp),0.0_dp,0.0_dp,'uniform endpoints')
    call close(irwinhall_censored_loglikelihood(one,2.0_dp),-inf,0.0_dp,'sum zero endpoint')
    one = make_censored_data(uncensored=[0.0_dp])
    call close(gausshyper_censored_loglikelihood(one,0.5_dp,2.0_dp,0.0_dp,1.0_dp), &
        inf,0.0_dp,'singular beta endpoint')
    call check(ieee_is_nan(irwinhall_censored_loglikelihood(empty,1.5_dp)),'integral shape')
    call check(ieee_is_nan(irwinhall_censored_loglikelihood(empty,33.0_dp)),'sum range')
    call check(ieee_is_nan(vonmises_line_censored_loglikelihood(empty,101.0_dp)),'concentration range')
    call check(ieee_is_nan(gausshyper_censored_loglikelihood(empty,0.01_dp,2.0_dp,0.0_dp,1.0_dp)), &
        'beta shape range')
    call check(ieee_is_nan(gausshyper_censored_loglikelihood(empty,1.0_dp,2.0_dp,0.0_dp,-0.9_dp)), &
        'tilt range')
    lower = [1.2_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp]
    upper = [3.8_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp]
    call fit_family(data,2,lower,upper,fit)
    call check(fit%success,'integer censored fit')
    call check(fit%params(1) == real(nint(fit%params(1)),dp),'integer fit result')
    x = nearest(1.0_dp,-1.0_dp)
    one = make_censored_data(interval=reshape([nearest(x,-1.0_dp),x],[1,2]))
    expected = log((1.0_dp-nearest(x,-1.0_dp))**0.1_dp-(1.0_dp-x)**0.1_dp)
    call close(gausshyper_censored_loglikelihood(one,1.0_dp,0.1_dp,0.0_dp,1.0_dp), &
        expected,3e-12_dp,'reflected endpoint interval')
    x = nearest(32.0_dp,-1.0_dp)
    one = make_censored_data(interval=reshape([nearest(x,-1.0_dp),x],[1,2]))
    expected = 32.0_dp*log(32.0_dp-nearest(x,-1.0_dp))-log_gamma(33.0_dp) &
        +log(1.0_dp-2.0_dp**(-32))
    call close(irwinhall_censored_loglikelihood(one,32.0_dp),expected,3e-13_dp, &
        'reflected spline endpoint interval')
    call close(gausshyper_logkernel_integral(0.0_dp,1.0_dp,2.0_dp,3.0_dp,0.0_dp,1.0_dp), &
        -log(12.0_dp),3e-13_dp,'kernel integral normalization')
    call close(gausshyper_logkernel_integral(0.4_dp,0.4_dp,2.0_dp,3.0_dp,0.0_dp,1.0_dp), &
        -inf,0.0_dp,'kernel zero width')
    call check(ieee_is_nan(gausshyper_logkernel_integral(-0.1_dp,0.5_dp, &
        2.0_dp,3.0_dp,0.0_dp,1.0_dp)),'kernel invalid endpoint')
    do family = 1, 3
        call parameters(family,par)
        n = 3
        if (family == 1) n = 6
        lower = par
        upper = par
        do j = 1, n
            lower(j) = par(j)-0.1_dp
            upper(j) = par(j)+0.1_dp
        end do
        if (family == 2) then
            lower(1) = 2.0_dp
            upper(1) = 4.0_dp
        end if
        select case (family)
        case (1)
            call gausshyper_fit(data,lower,upper,fit,guess=par,max_iter=20)
        case (2)
            call irwinhall_fit(data,lower(:3),upper(:3),fit,guess=par(:3),max_iter=20)
        case (3)
            call vonmises_line_fit(data,lower(:3),upper(:3),fit,guess=par(:3),max_iter=20)
        end select
        call check(ieee_is_finite(fit%nllf),'all-free finite fit')
        call check(fit%nllf <= -ll(data,family,par),'all-free fit improves initial candidate')
        call check(ieee_is_nan(ll(empty,family,[par(:n-2),inf,1.0_dp, &
            (0.0_dp,i=1,6-n)])),'infinite location rejected')
    end do
    call close(gausshyper_censored_nnlf(data,1.3_dp,2.7_dp,0.5_dp,1.0_dp), &
        -gausshyper_censored_loglikelihood(data,1.3_dp,2.7_dp,0.5_dp,1.0_dp),0.0_dp,'gausshyper NNLF')
    call close(irwinhall_censored_nnlf(data,3.0_dp), &
        -irwinhall_censored_loglikelihood(data,3.0_dp),0.0_dp,'irwinhall NNLF')
    call close(vonmises_line_censored_nnlf(data,1.3_dp), &
        -vonmises_line_censored_loglikelihood(data,1.3_dp),0.0_dp,'vonmises_line NNLF')
    write(*,'(a,i0,a,i0)') 'bounded censoring checks: ',checks,', failures: ',failures
    if (failures /= 0) error stop 1
contains
    subroutine parameters(family, par)
        implicit none
        integer, intent(in) :: family !! Gauss hypergeometric=1, Irwin-Hall=2, line von Mises=3
        real(dp), intent(out) :: par(6) !! shapes, loc, scale; unused tail ignored
        par = [1.7_dp,2.7_dp,0.5_dp,1.0_dp,0.0_dp,1.0_dp]
        if (family == 2) par = [3.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp]
        if (family == 3) par = [2.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp]
    end subroutine parameters
    pure function ll(data, family, par) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! observations
        integer, intent(in) :: family !! Gauss hypergeometric=1, Irwin-Hall=2, line von Mises=3
        real(dp), intent(in) :: par(6) !! shapes, loc, scale; unused tail ignored
        real(dp) :: value
        select case (family)
        case (1)
            value = gausshyper_censored_loglikelihood(data,par(1),par(2),par(3),par(4),par(5),par(6))
        case (2)
            value = irwinhall_censored_loglikelihood(data,par(1),par(2),par(3))
        case (3)
            value = vonmises_line_censored_loglikelihood(data,par(1),par(2),par(3))
        end select
    end function ll

    pure function scalar(x, family, par, part) result(value)
        implicit none
        real(dp), intent(in) :: x !! standardized observation
        integer, intent(in) :: part !! density=1, CDF=2, SF=3
        integer, intent(in) :: family !! Gauss hypergeometric=1, Irwin-Hall=2, line von Mises=3
        real(dp), intent(in) :: par(6) !! shapes, loc, scale; unused tail ignored
        real(dp) :: value
        select case (family)
        case (1)
            select case (part)
            case (1)
                value = gausshyper_logpdf(x,par(1),par(2),par(3),par(4))
            case (2)
                value = gausshyper_logcdf(x,par(1),par(2),par(3),par(4))
            case (3)
                value = gausshyper_logsf(x,par(1),par(2),par(3),par(4))
            end select
        case (2)
            select case (part)
            case (1)
                value = irwinhall_logpdf(x,par(1))
            case (2)
                value = irwinhall_logcdf(x,par(1))
            case (3)
                value = irwinhall_logsf(x,par(1))
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

    pure function quantile(p, family, par, upper) result(value)
        implicit none
        real(dp), intent(in) :: p !! probability in (0,1)
        logical, intent(in) :: upper !! true for survival quantile
        integer, intent(in) :: family !! Gauss hypergeometric=1, Irwin-Hall=2, line von Mises=3
        real(dp), intent(in) :: par(6) !! shapes, loc, scale; unused tail ignored
        real(dp) :: value
        select case (family)
        case (1)
            if (upper) then
                value = gausshyper_isf(p,par(1),par(2),par(3),par(4))
            else
                value = gausshyper_ppf(p,par(1),par(2),par(3),par(4))
            end if
        case (2)
            if (upper) then
                value = irwinhall_isf(p,par(1))
            else
                value = irwinhall_ppf(p,par(1))
            end if
        case (3)
            if (upper) then
                value = vonmises_line_isf(p,par(1))
            else
                value = vonmises_line_ppf(p,par(1))
            end if
        end select
    end function quantile

    subroutine fit_family(data, family, lower, upper, result)
        implicit none
        type(censored_data), intent(in) :: data !! observations
        real(dp), intent(in) :: lower(6), upper(6) !! parameter bounds in shape, loc, scale order
        type(fit_result), intent(out) :: result !! fit and diagnostics
        integer, intent(in) :: family !! Gauss hypergeometric=1, Irwin-Hall=2, line von Mises=3
        select case (family)
        case (1)
            call gausshyper_fit(data,lower(:6),upper(:6),result)
        case (2)
            call irwinhall_fit(data,lower(:3),upper(:3),result)
        case (3)
            call vonmises_line_fit(data,lower(:3),upper(:3),result)
        end select
    end subroutine fit_family

    subroutine fit_array(data, family, lower, upper, result)
        implicit none
        real(dp), intent(in) :: data(:) !! exact observations
        real(dp), intent(in) :: lower(6), upper(6) !! parameter bounds in shape, loc, scale order
        type(fit_result), intent(out) :: result !! fit and diagnostics
        integer, intent(in) :: family !! Gauss hypergeometric=1, Irwin-Hall=2, line von Mises=3
        select case (family)
        case (1)
            call gausshyper_fit(data,lower(:6),upper(:6),result)
        case (2)
            call irwinhall_fit(data,lower(:3),upper(:3),result)
        case (3)
            call vonmises_line_fit(data,lower(:3),upper(:3),result)
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
end program test_censored_bounded_final
