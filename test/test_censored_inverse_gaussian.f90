! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_inverse_gaussian
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        fit_status_invalid_input, &
        invgauss_censored_loglikelihood, invgauss_censored_nnlf, invgauss_fit, &
        invgauss_logpdf, invgauss_logcdf, invgauss_logsf, invgauss_ppf, invgauss_isf, &
        wald_censored_loglikelihood, wald_censored_nnlf, wald_fit, &
        wald_logpdf, wald_logcdf, wald_logsf, wald_ppf, wald_isf, &
        recipinvgauss_censored_loglikelihood, recipinvgauss_censored_nnlf, recipinvgauss_fit, &
        recipinvgauss_logpdf, recipinvgauss_logcdf, recipinvgauss_logsf, recipinvgauss_ppf, recipinvgauss_isf
    use scifort_math, only : positive_infinity, quiet_nan
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    implicit none
    type(censored_data) :: one, data, other, empty
    type(fit_result) :: fit, old
    real(dp) :: x, shape, actual, expected, previous, p, scale, loc, inf, nan
    real(dp) :: lower(3), upper(3), tiny_shape
    real(dp), parameter :: scipy_parameter(3) = &
        [0.3349859947937545_dp,0.47823955115344424_dp,4.999999874519365_dp]
    real(dp), parameter :: scipy_nllf(3) = &
        [10.735705684804694_dp,9.221866539671925_dp,11.23530622645976_dp]
    integer :: family, i, j, k, checks, failures
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    empty = make_censored_data()
    data = make_censored_data(uncensored=[0.2_dp,0.4_dp,0.6_dp],left=[0.1_dp], &
        right=[0.7_dp],interval=reshape([0.25_dp,0.5_dp,0.35_dp,0.55_dp],[2,2]))
    do family = 1, 3
        do j = 1, 4
            shape = real(j,dp)/3.0_dp
            previous = -inf
            do i = 1, 25
                x = real(i,dp)/5.0_dp
                one = make_censored_data(uncensored=[x])
                call close(ll(one,family,shape),scalar(x,family,shape,1),2e-12_dp,'density')
                one = make_censored_data(left=[x])
                actual = ll(one,family,shape)
                call close(actual,scalar(x,family,shape,2),2e-12_dp,'CDF')
                call check(actual >= previous .and. actual <= 0.0_dp,'CDF monotonicity')
                previous = actual
                one = make_censored_data(right=[x])
                call close(ll(one,family,shape),scalar(x,family,shape,3),2e-12_dp,'SF')
                one = make_censored_data(interval=reshape([x,nearest(x,1.0_dp)],[1,2]))
                expected = scalar(x,family,shape,1)+log(nearest(x,1.0_dp)-x)
                call close(ll(one,family,shape),expected,3e-13_dp,'adjacent interval')
                one = make_censored_data(interval=reshape([x,x+0.1_dp],[1,2]))
                if (scalar(x,family,shape,2) < log(0.5_dp)) then
                    expected = log(exp(scalar(x+0.1_dp,family,shape,2)) &
                        -exp(scalar(x,family,shape,2)))
                else
                    expected = log(exp(scalar(x,family,shape,3)) &
                        -exp(scalar(x+0.1_dp,family,shape,3)))
                end if
                call close(ll(one,family,shape),expected,2e-11_dp,'ordinary interval')
            end do
            do i = 1, 9
                p = real(i,dp)/10.0_dp
                x = quantile(p,family,shape,.false.)
                one = make_censored_data(left=[x])
                call close(ll(one,family,shape),log(p),3e-12_dp,'CDF/PPF')
                x = quantile(p,family,shape,.true.)
                one = make_censored_data(right=[x])
                call close(ll(one,family,shape),log(p),3e-12_dp,'SF/ISF')
            end do
        end do
        do k = -1, 1
            scale = 10.0_dp**(300*k)
            loc = -0.3_dp*scale
            x = 0.7_dp*scale
            one = make_censored_data(uncensored=[x])
            other = make_censored_data(uncensored=[1.0_dp])
            call close(ll(one,family,1.3_dp,loc,scale), &
                ll(other,family,1.3_dp)-log(scale),3e-13_dp,'affine density')
            one = make_censored_data(interval=reshape([x,nearest(x,1.0_dp)],[1,2]))
            expected = ll(other,family,1.3_dp)-log(scale)+log(nearest(x,1.0_dp)-x)
            call close(ll(one,family,1.3_dp,loc,scale),expected,3e-13_dp,'affine interval')
        end do
        call close(ll(empty,family,1.0_dp),0.0_dp,0.0_dp,'valid empty')
        if (family /= 2) then
            call check(ieee_is_nan(ll(empty,family,0.0_dp)),'zero shape')
            call check(ieee_is_nan(ll(empty,family,nan)),'NaN shape')
            call check(ieee_is_nan(ll(empty,family,-1.0_dp)),'negative shape')
            call check(ieee_is_nan(ll(empty,family,inf)),'infinite shape')
        end if
        call check(ieee_is_nan(ll(empty,family,1.0_dp,scale=0.0_dp)),'zero scale')
        call check(ieee_is_nan(ll(empty,family,1.0_dp,loc=inf)),'infinite location')
        call check(ieee_is_nan(ll(empty,family,1.0_dp,scale=nan)),'NaN scale')
        one = make_censored_data(uncensored=[0.0_dp])
        call close(ll(one,family,1.3_dp),-inf,0.0_dp,'endpoint density')
        one = make_censored_data(left=[0.0_dp])
        call close(ll(one,family,1.3_dp),-inf,0.0_dp,'endpoint CDF')
        one = make_censored_data(right=[0.0_dp])
        call close(ll(one,family,1.3_dp),0.0_dp,0.0_dp,'endpoint SF')
        one = make_censored_data(uncensored=[inf])
        call check(ieee_is_nan(ll(one,family,1.3_dp)),'infinite exact invalid')
        one = make_censored_data(left=[nan])
        call check(ieee_is_nan(ll(one,family,1.3_dp)),'NaN data invalid')
        one = make_censored_data(interval=reshape([0.7_dp,0.7_dp],[1,2]))
        call close(ll(one,family,1.3_dp),scalar(0.7_dp,family,1.3_dp,1),2e-13_dp,'equal endpoint')
        one = make_censored_data(interval=reshape([-inf,0.7_dp],[1,2]))
        call close(ll(one,family,1.3_dp),scalar(0.7_dp,family,1.3_dp,2),2e-13_dp,'left interval')
        one = make_censored_data(interval=reshape([0.7_dp,inf],[1,2]))
        call close(ll(one,family,1.3_dp),scalar(0.7_dp,family,1.3_dp,3),2e-13_dp,'right interval')
        one = make_censored_data(interval=reshape([0.8_dp,0.7_dp],[1,2]))
        call check(ieee_is_nan(ll(one,family,1.3_dp)),'reversed interval')
        one = make_censored_data(interval=reshape([-inf,inf],[1,2]))
        call check(ieee_is_nan(ll(one,family,1.3_dp)),'unbounded interval invalid')
        lower = [1.3_dp,0.0_dp,1.0_dp]
        upper = lower
        call fit_family(data,family,lower,upper,fit)
        call check(fit%success,'fixed fit success')
        call close(fit%nllf,-ll(data,family,1.3_dp),2e-13_dp,'fixed fit NNLF')
        call fit_family(empty,family,lower,upper,fit)
        call check(fit%status == fit_status_invalid_input,'empty fit rejected')
        one = make_censored_data(uncensored=[0.2_dp,0.4_dp,0.6_dp])
        call fit_family(one,family,lower,upper,fit)
        call fit_array([0.2_dp,0.4_dp,0.6_dp],family,lower,upper,old)
        call check(old%success .and. fit%success,'original fit generic')
        call close(fit%nllf,old%nllf,2e-13_dp,'exact/array fit agreement')
        lower = [0.3_dp,0.0_dp,1.3_dp]
        upper = [5.0_dp,0.0_dp,1.3_dp]
        if (family == 2) then
            lower(3) = 0.3_dp
            upper(3) = 5.0_dp
        end if
        call fit_family(data,family,lower,upper,fit)
        call check(fit%success,'one-parameter fit converged')
        j = 1
        if (family == 2) j = 2
        call close(fit%params(j),scipy_parameter(family),3e-5_dp,'SciPy optimum parameter')
        call close(fit%nllf,scipy_nllf(family),2e-8_dp,'SciPy optimum likelihood')
        upper = [5.0_dp,0.05_dp,2.0_dp]
        lower = [0.3_dp,-0.2_dp,0.5_dp]
        call fit_family(data,family,lower,upper,fit)
        call check(ieee_is_finite(fit%nllf),'all-free finite candidate')
        call check(fit%nllf <= -ll(data,family,1.0_dp,0.0_dp,1.0_dp), &
            'all-free candidate improvement')
    end do
    do i = 1, 30
        x = 0.1_dp*real(i,dp)
        one = make_censored_data(right=[x])
        call close(wald_censored_loglikelihood(one),invgauss_censored_loglikelihood(one,1.0_dp), &
            0.0_dp,'Wald identity')
        other = make_censored_data(left=[1.0_dp/x])
        call close(recipinvgauss_censored_loglikelihood(other,1.7_dp), &
            invgauss_censored_loglikelihood(one,1.7_dp),2e-13_dp,'reciprocal tail identity')
    end do
    one = make_censored_data(right=[1e300_dp])
    expected = -0.5e300_dp
    call close(wald_censored_loglikelihood(one),expected,2e-13_dp,'huge right tail')
    one = make_censored_data(left=[1e-300_dp])
    call close(wald_censored_loglikelihood(one),expected,2e-13_dp,'tiny left tail')
    ! Large shapes approach Levy; retain its very small complementary log probability.
    one = make_censored_data(left=[1e300_dp])
    expected = -sqrt(2.0_dp/acos(-1.0_dp))*1e-150_dp
    actual = invgauss_censored_loglikelihood(one,1e300_dp)
    call close(actual/expected,1.0_dp,3e-13_dp,'tiny complementary log relative')
    tiny_shape = nearest(0.0_dp,1.0_dp)
    one = make_censored_data(uncensored=[tiny_shape])
    actual = invgauss_censored_loglikelihood(one,tiny_shape)
    expected = -0.5_dp*log(2.0_dp*acos(-1.0_dp))-1.5_dp*log(tiny_shape)
    call close(actual,expected,2e-13_dp,'minimum subnormal shape at mode coordinate')
    shape = 2.0_dp**(-1000)
    x = nearest(shape,1.0_dp)
    one = make_censored_data(uncensored=[x])
    p = epsilon(1.0_dp)
    expected = -0.5_dp*(p/sqrt(shape))**2/(1.0_dp+p) &
        -0.5_dp*log(2.0_dp*acos(-1.0_dp))-1.5_dp*log(x)
    call close(invgauss_censored_loglikelihood(one,shape),expected,6e-13_dp, &
        'tiny shape adjacent to center')
    one = make_censored_data(interval=reshape([0.7_dp,nearest(0.7_dp,1.0_dp)],[1,2]))
    call close(invgauss_censored_loglikelihood(one,1e-100_dp),-3.5e199_dp, &
        3e-13_dp,'steep falling adjacent interval')
    expected = -0.5_dp/(nearest(0.7_dp,1.0_dp)*1e-200_dp)
    call close(recipinvgauss_censored_loglikelihood(one,1e-100_dp),expected, &
        3e-13_dp,'steep rising adjacent interval')
    call close(invgauss_censored_nnlf(data,1.3_dp),-ll(data,1,1.3_dp),0.0_dp,'invgauss NNLF')
    call close(wald_censored_nnlf(data),-ll(data,2,1.0_dp),0.0_dp,'Wald NNLF')
    call close(recipinvgauss_censored_nnlf(data,1.3_dp),-ll(data,3,1.3_dp),0.0_dp,'reciprocal NNLF')
    write(*,'(a,i0,a,i0)') 'inverse Gaussian censoring checks: ',checks,', failures: ',failures
    if (failures /= 0) error stop 1
contains
    pure function ll(data, family, shape, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! observations
        integer, intent(in) :: family !! inverse Gaussian=1, Wald=2, reciprocal=3
        real(dp), intent(in) :: shape !! positive shape, ignored for Wald
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        select case (family)
        case (1)
            value = invgauss_censored_loglikelihood(data,shape,loc,scale)
        case (2)
            value = wald_censored_loglikelihood(data,loc,scale)
        case (3)
            value = recipinvgauss_censored_loglikelihood(data,shape,loc,scale)
        end select
    end function ll

    pure function scalar(x, family, shape, part) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or limit
        integer, intent(in) :: family !! inverse Gaussian=1, Wald=2, reciprocal=3
        real(dp), intent(in) :: shape !! positive shape, ignored for Wald
        integer, intent(in) :: part !! density=1, CDF=2, SF=3
        real(dp) :: value
        select case (family)
        case (1)
            select case (part)
            case (1)
                value = invgauss_logpdf(x,shape)
            case (2)
                value = invgauss_logcdf(x,shape)
            case (3)
                value = invgauss_logsf(x,shape)
            end select
        case (2)
            select case (part)
            case (1)
                value = wald_logpdf(x)
            case (2)
                value = wald_logcdf(x)
            case (3)
                value = wald_logsf(x)
            end select
        case (3)
            select case (part)
            case (1)
                value = recipinvgauss_logpdf(x,shape)
            case (2)
                value = recipinvgauss_logcdf(x,shape)
            case (3)
                value = recipinvgauss_logsf(x,shape)
            end select
        end select
    end function scalar

    pure function quantile(p, family, shape, upper) result(value)
        implicit none
        real(dp), intent(in) :: p !! probability in (0,1)
        integer, intent(in) :: family !! inverse Gaussian=1, Wald=2, reciprocal=3
        real(dp), intent(in) :: shape !! positive shape, ignored for Wald
        logical, intent(in) :: upper !! true for survival quantile
        real(dp) :: value
        select case (family)
        case (1)
            if (upper) then
                value = invgauss_isf(p,shape)
            else
                value = invgauss_ppf(p,shape)
            end if
        case (2)
            if (upper) then
                value = wald_isf(p)
            else
                value = wald_ppf(p)
            end if
        case (3)
            if (upper) then
                value = recipinvgauss_isf(p,shape)
            else
                value = recipinvgauss_ppf(p,shape)
            end if
        end select
    end function quantile

    subroutine fit_family(data, family, lower, upper, result)
        implicit none
        type(censored_data), intent(in) :: data !! observations
        integer, intent(in) :: family !! inverse Gaussian=1, Wald=2, reciprocal=3
        real(dp), intent(in) :: lower(3) !! shape, loc, scale lower bounds; Wald ignores shape
        real(dp), intent(in) :: upper(3) !! upper bounds; equality fixes parameters
        type(fit_result), intent(out) :: result !! fit and diagnostics
        select case (family)
        case (1)
            call invgauss_fit(data,lower,upper,result)
        case (2)
            call wald_fit(data,lower(2:3),upper(2:3),result)
        case (3)
            call recipinvgauss_fit(data,lower,upper,result)
        end select
    end subroutine fit_family

    subroutine fit_array(data, family, lower, upper, result)
        implicit none
        real(dp), intent(in) :: data(:) !! observations
        integer, intent(in) :: family !! inverse Gaussian=1, Wald=2, reciprocal=3
        real(dp), intent(in) :: lower(3) !! shape, loc, scale lower bounds; Wald ignores shape
        real(dp), intent(in) :: upper(3) !! upper bounds; equality fixes parameters
        type(fit_result), intent(out) :: result !! fit and diagnostics
        select case (family)
        case (1)
            call invgauss_fit(data,lower,upper,result)
        case (2)
            call wald_fit(data,lower(2:3),upper(2:3),result)
        case (3)
            call recipinvgauss_fit(data,lower,upper,result)
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
end program test_censored_inverse_gaussian

