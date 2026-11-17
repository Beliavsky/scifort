! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_beta_t
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        normal_censored_loglikelihood, cauchy_censored_loglikelihood, &
        betaprime_censored_loglikelihood, fit_status_invalid_input, &
        t_censored_loglikelihood, t_censored_nnlf, t_fit, &
        t_logpdf, t_logcdf, t_logsf, t_ppf, t_isf, &
        f_censored_loglikelihood, f_censored_nnlf, f_fit, &
        f_logpdf, f_logcdf, f_logsf, f_ppf, f_isf, &
        jf_skew_t_censored_loglikelihood, jf_skew_t_censored_nnlf, jf_skew_t_fit, &
        jf_skew_t_logpdf, jf_skew_t_logcdf, jf_skew_t_logsf, jf_skew_t_ppf, jf_skew_t_isf
    use scifort_math, only : positive_infinity, quiet_nan
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    implicit none
    type(censored_data) :: one, data, other, empty
    type(fit_result) :: fit, old
    real(dp) :: x, a, b, actual, expected, previous, p, scale, loc, inf, nan, c
    real(dp) :: lower(4), upper(4)
    real(dp), parameter :: scipy_shape(3) = &
        [5.0_dp,1.1741350701452165_dp,2.2316074118819835_dp]
    real(dp), parameter :: scipy_objective(3) = &
        [13.595649458888033_dp,12.162483724958465_dp,13.427914058395741_dp]
    real(dp), parameter :: minimum_mixed(3) = &
        [-709.4260379494452_dp,-710.1191851300052_dp,-709.4260379494452_dp]
    real(dp), parameter :: minimum_interval(3) = &
        [-708.3477977876847_dp,-709.0409449682446_dp,-708.3477977876847_dp]
    real(dp), parameter :: maximum_mixed(3) = &
        [-2.3052331443245633_dp,3.9093749540976286_dp,-2.3052331443245633_dp]
    real(dp), parameter :: maximum_interval(3) = &
        [-745.3590107045859_dp,-30.747984073899637_dp,-745.3590107045859_dp]
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
            a = real(j,dp)/3.0_dp
            b = 1.7_dp
            previous = -inf
            do i = 1, 25
                x = -3.0_dp+real(i,dp)/4.0_dp
                if (family == 2) x = real(i,dp)/5.0_dp
                one = make_censored_data(uncensored=[x])
                call close(ll(one,family,a,b),scalar(x,family,a,b,1),2e-12_dp,'density')
                one = make_censored_data(left=[x])
                actual = ll(one,family,a,b)
                call close(actual,scalar(x,family,a,b,2),2e-12_dp,'CDF')
                call check(actual >= previous .and. actual <= 0.0_dp,'CDF monotonicity')
                previous = actual
                one = make_censored_data(right=[x])
                call close(ll(one,family,a,b),scalar(x,family,a,b,3),2e-12_dp,'SF')
                one = make_censored_data(interval=reshape([x,nearest(x,1.0_dp)],[1,2]))
                expected = scalar(x,family,a,b,1)+log(nearest(x,1.0_dp)-x)
                call close(ll(one,family,a,b),expected,3e-13_dp,'adjacent interval')
            end do
            do i = 1, 9
                p = real(i,dp)/10.0_dp
                x = quantile(p,family,a,b,.false.)
                one = make_censored_data(left=[x])
                call close(ll(one,family,a,b),log(p),3e-12_dp,'CDF/PPF')
                x = quantile(p,family,a,b,.true.)
                one = make_censored_data(right=[x])
                call close(ll(one,family,a,b),log(p),3e-12_dp,'SF/ISF')
            end do
        end do
        do k = -1, 1
            scale = 10.0_dp**(200*k)
            loc = -0.3_dp*scale
            x = 0.7_dp*scale
            one = make_censored_data(uncensored=[x])
            other = make_censored_data(uncensored=[1.0_dp])
            call close(ll(one,family,1.3_dp,1.7_dp,loc,scale), &
                ll(other,family,1.3_dp,1.7_dp)-log(scale),3e-13_dp,'affine density')
            one = make_censored_data(interval=reshape([x,nearest(x,1.0_dp)],[1,2]))
            expected = ll(other,family,1.3_dp,1.7_dp)-log(scale)+log(nearest(x,1.0_dp)-x)
            call close(ll(one,family,1.3_dp,1.7_dp,loc,scale),expected,3e-13_dp,'affine interval')
        end do
        call close(ll(empty,family,1.0_dp,1.7_dp),0.0_dp,0.0_dp,'valid empty')
        call check(ieee_is_nan(ll(empty,family,0.0_dp,1.7_dp)),'zero first shape')
        call check(ieee_is_nan(ll(empty,family,nan,1.7_dp)),'NaN first shape')
        call check(ieee_is_nan(ll(empty,family,-1.0_dp,1.7_dp)),'negative first shape')
        call check(ieee_is_nan(ll(empty,family,1.0_dp,1.7_dp,scale=0.0_dp)),'zero scale')
        call check(ieee_is_nan(ll(empty,family,1.0_dp,1.7_dp,loc=inf)),'infinite location')
        if (family > 1) then
            call check(ieee_is_nan(ll(empty,family,1.0_dp,0.0_dp)),'zero second shape')
            call check(ieee_is_nan(ll(empty,family,1.0_dp,inf)),'infinite second shape')
        end if
        call check(ieee_is_nan(ll(empty,family,nearest(0.0_dp,1.0_dp),1.0_dp)), &
            'subnormal beta shape unsupported')
        call check(ieee_is_nan(ll(empty,family,2.0e6_dp,1.0_dp)),'shape sum above cap')
        one = make_censored_data(uncensored=[nan])
        call check(ieee_is_nan(ll(one,family,1.0_dp,1.0_dp)),'NaN exact')
        one = make_censored_data(uncensored=[inf])
        call check(ieee_is_nan(ll(one,family,1.0_dp,1.0_dp)),'infinite exact')
        one = make_censored_data(interval=reshape([0.3_dp,0.3_dp],[1,2]))
        other = make_censored_data(uncensored=[0.3_dp])
        call close(ll(one,family,1.3_dp,1.7_dp),ll(other,family,1.3_dp,1.7_dp),0.0_dp, &
            'equal interval exact')
        one = make_censored_data(interval=reshape([-inf,0.3_dp],[1,2]))
        other = make_censored_data(left=[0.3_dp])
        call close(ll(one,family,1.3_dp,1.7_dp),ll(other,family,1.3_dp,1.7_dp),0.0_dp, &
            'left-infinite interval')
        one = make_censored_data(interval=reshape([0.3_dp,inf],[1,2]))
        other = make_censored_data(right=[0.3_dp])
        call close(ll(one,family,1.3_dp,1.7_dp),ll(other,family,1.3_dp,1.7_dp),0.0_dp, &
            'right-infinite interval')
        a = tiny(1.0_dp)
        if (family < 3) a = 2.0_dp*a
        call close(ll(empty,family,a,a),0.0_dp,0.0_dp,'minimum normal beta shapes')
        one = make_censored_data(uncensored=[0.7_dp],left=[0.7_dp],right=[0.7_dp])
        call close(ll(one,family,a,a),minimum_mixed(family),3e-13_dp,'minimum mixed reference')
        one = make_censored_data(interval=reshape([0.7_dp,2.0_dp],[1,2]))
        call check(ieee_is_finite(ll(one,family,a,a)),'tiny shape wide interval')
        call close(ll(one,family,a,a),minimum_interval(family),3e-13_dp,'minimum interval reference')
        c = 1.0e6_dp
        x = 0.0_dp
        if (family == 2) x = 1.0_dp
        if (family == 3) c = 5.0e5_dp
        one = make_censored_data(uncensored=[x],left=[x],right=[x])
        call close(ll(one,family,c,c),maximum_mixed(family),2e-11_dp,'maximum shape reference')
        one = make_censored_data(interval=reshape([x,nearest(x,1.0_dp)],[1,2]))
        call close(ll(one,family,c,c),maximum_interval(family),3e-13_dp,'maximum adjacent reference')
        lower = [1.3_dp,1.7_dp,0.0_dp,1.3_dp]
        call fit_family(data,family,lower,lower,fit)
        call check(fit%success,'fixed fit')
        call close(fit%nllf,-ll(data,family,1.3_dp,1.7_dp,scale=1.3_dp),0.0_dp,'fixed objective')
        select case (family)
        case (1)
            actual = t_censored_nnlf(data,1.3_dp,scale=1.3_dp)
        case (2)
            actual = f_censored_nnlf(data,1.3_dp,1.7_dp,scale=1.3_dp)
        case (3)
            actual = jf_skew_t_censored_nnlf(data,1.3_dp,1.7_dp,scale=1.3_dp)
        end select
        call close(actual,fit%nllf,0.0_dp,'NNLF')
        call fit_array([0.2_dp,0.4_dp,0.6_dp],family,lower,lower,old)
        call check(old%success,'array fit preserved')
        one = make_censored_data(uncensored=[0.2_dp,0.4_dp,0.6_dp])
        call fit_family(one,family,lower,lower,fit)
        call close(fit%nllf,old%nllf,2e-12_dp,'array/censored agreement')
        lower(1) = 0.3_dp
        upper = lower
        upper(1) = 5.0_dp
        call fit_family(data,family,lower,upper,fit)
        call check(fit%success,'shape fit')
        call close(fit%params(1),scipy_shape(family),3e-5_dp,'SciPy optimum shape')
        call close(fit%nllf,scipy_objective(family),2e-8_dp,'SciPy optimum objective')
        call check(fit%nllf <= -ll(data,family,1.3_dp,1.7_dp,scale=1.3_dp),'shape improves')
        lower = [0.3_dp,0.3_dp,-0.2_dp,0.5_dp]
        upper = [5.0_dp,5.0_dp,0.05_dp,2.0_dp]
        expected = -ll(data,family,2.65_dp,2.65_dp,-0.075_dp,1.25_dp)
        call fit_family(data,family,lower,upper,fit)
        call check(ieee_is_finite(fit%nllf),'all-free finite candidate')
        call check(fit%nllf <= expected,'all-free improves midpoint')
        call fit_family(empty,family,lower,upper,fit)
        call check(fit%status == fit_status_invalid_input,'empty fit invalid')
    end do
    call close(t_censored_loglikelihood(data,1.0_dp),cauchy_censored_loglikelihood(data), &
        3e-13_dp,'t/Cauchy')
    call close(t_censored_loglikelihood(data,inf),normal_censored_loglikelihood(data), &
        0.0_dp,'normal limit')
    call close(t_censored_loglikelihood(data,3.0_dp), &
        jf_skew_t_censored_loglikelihood(data,1.5_dp,1.5_dp),0.0_dp,'symmetric skew-t')
    call close(f_censored_loglikelihood(data,4.0_dp,6.0_dp), &
        betaprime_censored_loglikelihood(data,2.0_dp,3.0_dp,scale=1.5_dp),3e-13_dp,'F/beta-prime')
    one = make_censored_data(uncensored=[0.0_dp])
    call check(f_censored_loglikelihood(one,1.0_dp,3.0_dp) == inf,'F singular endpoint')
    call close(f_censored_loglikelihood(one,2.0_dp,3.0_dp),0.0_dp,3e-13_dp,'F finite endpoint')
    call check(f_censored_loglikelihood(one,3.0_dp,3.0_dp) == -inf,'F zero endpoint')
    one = make_censored_data(left=[0.7_dp])
    call close(f_censored_loglikelihood(one,1e-300_dp,1e-20_dp)/1e-280_dp,-1.0_dp, &
        3e-12_dp,'F tiny complementary log CDF')
    call close(jf_skew_t_censored_loglikelihood(one,1e-300_dp,1e-20_dp)/1e-280_dp,-1.0_dp, &
        3e-12_dp,'skew-t tiny complementary log CDF')
    one = make_censored_data(right=[0.7_dp])
    call close(f_censored_loglikelihood(one,1e-20_dp,1e-300_dp)/1e-280_dp,-1.0_dp, &
        3e-12_dp,'F tiny complementary log SF')
    call close(jf_skew_t_censored_loglikelihood(one,1e-20_dp,1e-300_dp)/1e-280_dp,-1.0_dp, &
        3e-12_dp,'skew-t tiny complementary log SF')
    one = make_censored_data(uncensored=[-1e300_dp])
    other = make_censored_data(uncensored=[1e300_dp])
    call close(jf_skew_t_censored_loglikelihood(one,0.3_dp,2.4_dp), &
        jf_skew_t_censored_loglikelihood(other,2.4_dp,0.3_dp),3e-13_dp,'reflected far density')
    one = make_censored_data(left=[-1e300_dp])
    other = make_censored_data(right=[1e300_dp])
    call close(jf_skew_t_censored_loglikelihood(one,0.3_dp,2.4_dp), &
        jf_skew_t_censored_loglikelihood(other,2.4_dp,0.3_dp),3e-13_dp,'reflected far tail')
    one = make_censored_data(interval=reshape([-2e300_dp,-1e300_dp],[1,2]))
    other = make_censored_data(interval=reshape([1e300_dp,2e300_dp],[1,2]))
    call close(jf_skew_t_censored_loglikelihood(one,0.3_dp,2.4_dp), &
        jf_skew_t_censored_loglikelihood(other,2.4_dp,0.3_dp),3e-13_dp,'reflected far interval')
    call close(jf_skew_t_censored_loglikelihood(one,1e-300_dp,2e-300_dp), &
        jf_skew_t_censored_loglikelihood(other,2e-300_dp,1e-300_dp),3e-13_dp,'reflected tiny interval')
    write(*,'(a,i0,a,i0)') 'Beta/t censored checks: ',checks,', failures: ',failures
    if (failures /= 0) error stop 1
contains
    pure function ll(data, family, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! observations
        integer, intent(in) :: family !! t=1, F=2, Jones-Faddy=3
        real(dp), intent(in) :: a, b !! family shape parameters; b ignored for t
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive scale, default one
        real(dp) :: value
        select case (family)
        case (1)
            value = t_censored_loglikelihood(data,a,loc,scale)
        case (2)
            value = f_censored_loglikelihood(data,a,b,loc,scale)
        case (3)
            value = jf_skew_t_censored_loglikelihood(data,a,b,loc,scale)
        end select
    end function ll

    pure function scalar(x, family, a, b, part) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        integer, intent(in) :: family !! t=1, F=2, Jones-Faddy=3
        real(dp), intent(in) :: a, b !! family shapes; b ignored for t
        integer, intent(in) :: part !! density=1, CDF=2, SF=3
        real(dp) :: value
        select case (family)
        case (1)
            select case (part)
            case (1)
                value = t_logpdf(x,a)
            case (2)
                value = t_logcdf(x,a)
            case (3)
                value = t_logsf(x,a)
            end select
        case (2)
            select case (part)
            case (1)
                value = f_logpdf(x,a,b)
            case (2)
                value = f_logcdf(x,a,b)
            case (3)
                value = f_logsf(x,a,b)
            end select
        case (3)
            select case (part)
            case (1)
                value = jf_skew_t_logpdf(x,a,b)
            case (2)
                value = jf_skew_t_logcdf(x,a,b)
            case (3)
                value = jf_skew_t_logsf(x,a,b)
            end select
        end select
    end function scalar

    pure function quantile(p, family, a, b, upper) result(value)
        implicit none
        real(dp), intent(in) :: p !! probability strictly between zero and one
        integer, intent(in) :: family !! t=1, F=2, Jones-Faddy=3
        real(dp), intent(in) :: a, b !! family shapes; b ignored for t
        logical, intent(in) :: upper !! true for the survival quantile
        real(dp) :: value
        select case (family)
        case (1)
            if (upper) then
                value = t_isf(p,a)
            else
                value = t_ppf(p,a)
            end if
        case (2)
            if (upper) then
                value = f_isf(p,a,b)
            else
                value = f_ppf(p,a,b)
            end if
        case (3)
            if (upper) then
                value = jf_skew_t_isf(p,a,b)
            else
                value = jf_skew_t_ppf(p,a,b)
            end if
        end select
    end function quantile

    subroutine fit_family(data, family, lower, upper, result)
        implicit none
        type(censored_data), intent(in) :: data !! observations
        integer, intent(in) :: family !! t=1, F=2, Jones-Faddy=3
        real(dp), intent(in) :: lower(4) !! shape, second shape, loc, scale lower bounds
        real(dp), intent(in) :: upper(4) !! upper bounds; equality fixes parameters
        type(fit_result), intent(out) :: result !! fit and diagnostics
        select case (family)
        case (1)
            call t_fit(data,lower([1,3,4]),upper([1,3,4]),result)
        case (2)
            call f_fit(data,lower,upper,result)
        case (3)
            call jf_skew_t_fit(data,lower,upper,result)
        end select
    end subroutine fit_family

    subroutine fit_array(data, family, lower, upper, result)
        implicit none
        real(dp), intent(in) :: data(:) !! observations
        integer, intent(in) :: family !! t=1, F=2, Jones-Faddy=3
        real(dp), intent(in) :: lower(4) !! shape, second shape, loc, scale lower bounds
        real(dp), intent(in) :: upper(4) !! upper bounds; equality fixes parameters
        type(fit_result), intent(out) :: result !! fit and diagnostics
        select case (family)
        case (1)
            call t_fit(data,lower([1,3,4]),upper([1,3,4]),result)
        case (2)
            call f_fit(data,lower,upper,result)
        case (3)
            call jf_skew_t_fit(data,lower,upper,result)
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
end program test_censored_beta_t
