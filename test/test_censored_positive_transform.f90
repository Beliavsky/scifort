! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_positive_transform
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        invgamma_censored_loglikelihood, invgamma_censored_nnlf, invgamma_fit, &
        invgamma_logpdf, invgamma_logcdf, invgamma_logsf, invgamma_ppf, &
        betaprime_censored_loglikelihood, betaprime_censored_nnlf, betaprime_fit, &
        betaprime_logpdf, betaprime_logcdf, betaprime_logsf, betaprime_ppf, &
        maxwell_censored_loglikelihood, maxwell_censored_nnlf, maxwell_fit, &
        maxwell_logpdf, maxwell_logcdf, maxwell_logsf, maxwell_ppf, &
        gamma_censored_loglikelihood, invweibull_censored_loglikelihood, &
        lomax_censored_loglikelihood, &
        fit_status_invalid_input, fit_status_no_finite_objective, fit_status_max_iter
    use scifort_math, only : quiet_nan, positive_infinity
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
    implicit none
    type(censored_data) :: data, single, empty, other
    type(fit_result) :: fit, old
    real(dp) :: inf, nan, actual, expected, x, c, p, previous, reference(3), lower(3), upper(3), scipy_parameter(3)
    integer :: family, i, j, checks, failures
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    reference = [-15.941918349997531_dp,-14.960817730539578_dp,-19.552136120161531_dp]
    scipy_parameter = [1.3391601562500006_dp,2.2651367187500027_dp,1.4474879497966218_dp]
    data = make_censored_data(uncensored=[0.5_dp,1.2_dp,2.5_dp],left=[0.3_dp], &
        right=[3.5_dp],interval=reshape([1.1_dp,4.0_dp,1.4_dp,4.4_dp],[2,2]))
    do family = 1, 3
        actual = likelihood(data,family,1.7_dp,scale=1.3_dp)
        call close(actual,reference(family),3.0e-13_dp,'mixed mpmath reference')
        select case (family)
        case (1)
            expected = invgamma_censored_nnlf(data,1.7_dp,scale=1.3_dp)
        case (2)
            expected = betaprime_censored_nnlf(data,1.7_dp,2.3_dp,scale=1.3_dp)
        case (3)
            expected = maxwell_censored_nnlf(data,scale=1.3_dp)
        end select
        call close(expected,-actual,0.0_dp,'negative likelihood')
        do j = 1, 3
            c = 0.4_dp
            if (j == 2) c = 1.0_dp
            if (j == 3) c = 3.0_dp
            previous = -inf
            do i = 1, 12
                x = real(i,dp)/3.0_dp
                single = make_censored_data(uncensored=[x])
                select case (family)
                case (1)
                    expected = invgamma_logpdf(x,c)
                case (2)
                    expected = betaprime_logpdf(x,c,2.3_dp)
                case (3)
                    expected = maxwell_logpdf(x)
                end select
                if (ieee_is_finite(expected)) then
                    call close(likelihood(single,family,c),expected,6.0e-14_dp,'scalar density')
                else
                    call check(likelihood(single,family,c) == expected,'scalar center density')
                end if
                single = make_censored_data(left=[x])
                select case (family)
                case (1)
                    expected = invgamma_logcdf(x,c)
                case (2)
                    expected = betaprime_logcdf(x,c,2.3_dp)
                case (3)
                    expected = maxwell_logcdf(x)
                end select
                actual = likelihood(single,family,c)
                call close(actual,expected,6.0e-14_dp,'scalar CDF')
                call check(actual >= previous .and. actual <= 0.0_dp,'CDF monotonicity and range')
                previous = actual
                single = make_censored_data(right=[x])
                select case (family)
                case (1)
                    expected = invgamma_logsf(x,c)
                case (2)
                    expected = betaprime_logsf(x,c,2.3_dp)
                case (3)
                    expected = maxwell_logsf(x)
                end select
                call close(likelihood(single,family,c),expected,6.0e-14_dp,'scalar SF')
            end do
            do i = 1, 9
                p = real(i,dp)/10.0_dp
                select case (family)
                case (1)
                    x = invgamma_ppf(p,c)
                case (2)
                    x = betaprime_ppf(p,c,2.3_dp)
                case (3)
                    x = maxwell_ppf(p)
                end select
                single = make_censored_data(left=[x])
                call close(exp(likelihood(single,family,c)),p,5.0e-14_dp,'CDF/PPF round trip')
            end do
        end do
        single = make_censored_data(uncensored=[0.5_dp,1.0_dp,2.0_dp])
        lower = [1.7_dp,0.0_dp,1.3_dp]
        upper = lower
        call fitting(single,family,lower,upper,fit)
        select case (family)
        case (1)
            call invgamma_fit(single%uncensored,lower,upper,old)
        case (2)
            call betaprime_fit(single%uncensored,[lower(1),2.3_dp,lower(2:3)], &
                [upper(1),2.3_dp,upper(2:3)],old)
        case (3)
            call maxwell_fit(single%uncensored,lower(2:3),upper(2:3),old)
        end select
        call check(fit%success .and. old%success,'fixed fits converge')
        call close(fit%nllf,old%nllf,2.0e-13_dp,'array fit compatibility')
        call close(fit%nllf,-likelihood(single,family,1.7_dp,scale=1.3_dp), &
            0.0_dp,'fit objective agreement')
        call close(likelihood(empty,family,1.0_dp),0.0_dp,0.0_dp,'empty valid sample')
        if (family /= 3) call check(ieee_is_nan(likelihood(empty,family,nan)),'NaN shape')
        if (family /= 3) call check(ieee_is_nan(likelihood(empty,family,inf)),'infinite shape')
        if (family /= 3) call check(ieee_is_nan(likelihood(empty,family,0.0_dp)),'zero shape')
        if (family /= 3) call check(ieee_is_nan(likelihood(empty,family,-1.0_dp)),'negative shape')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,loc=nan)),'NaN location')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,loc=inf)),'infinite location')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,scale=nan)),'NaN scale')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,scale=inf)),'infinite scale')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,scale=0.0_dp)),'zero scale')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,scale=-1.0_dp)),'negative scale')
        single = make_censored_data(uncensored=[nan])
        call check(ieee_is_nan(likelihood(single,family,1.0_dp)),'NaN exact value')
        single = make_censored_data(uncensored=[inf])
        call check(ieee_is_nan(likelihood(single,family,1.0_dp)),'infinite exact value')
        single = make_censored_data(interval=reshape([1.0_dp,0.0_dp],[1,2]))
        call check(ieee_is_nan(likelihood(single,family,1.0_dp)),'reversed interval')
        single = make_censored_data(interval=reshape([-inf,1.0_dp],[1,2]))
        other = make_censored_data(left=[1.0_dp])
        call close(likelihood(single,family,1.0_dp),likelihood(other,family,1.0_dp), &
            0.0_dp,'infinite lower interval endpoint')
        single = make_censored_data(interval=reshape([1.0_dp,inf],[1,2]))
        other = make_censored_data(right=[1.0_dp])
        call close(likelihood(single,family,1.0_dp),likelihood(other,family,1.0_dp), &
            0.0_dp,'infinite upper interval endpoint')
        single = make_censored_data(interval=reshape([1.0_dp,1.0_dp],[1,2]))
        other = make_censored_data(uncensored=[1.0_dp])
        call close(likelihood(single,family,1.0_dp),likelihood(other,family,1.0_dp), &
            0.0_dp,'equal endpoints represent exact data')
        lower = [0.3_dp,-0.1_dp,0.8_dp]
        upper = [3.0_dp,0.1_dp,2.0_dp]
        call fitting(data,family,lower,upper,fit)
        call check(fit%success,'all-free bounded fit')
        if (allocated(fit%params)) then
            if (family == 1) expected = likelihood(data,family,fit%params(1), &
                fit%params(2),fit%params(3))
            if (family == 2) expected = likelihood(data,family,fit%params(1), &
                fit%params(3),fit%params(4))
            if (family == 3) expected = likelihood(data,family,1.0_dp, &
                fit%params(1),fit%params(2))
            call close(fit%nllf,-expected,0.0_dp,'free fit objective')
        end if
        lower = [1.7_dp,0.0_dp,1.3_dp]
        upper = [1.7_dp,0.0_dp,3.0_dp]
        call fitting(data,family,lower,upper,fit)
        call check(fit%success,'partly fixed fit')
        lower = [0.01_dp,0.0_dp,1.3_dp]
        upper = [10.0_dp,0.0_dp,1.3_dp]
        call fitting(data,family,lower,upper,fit)
        call check(fit%success,'fixed support shape fit')
        call check(fit%nllf <= -reference(family),'shape fit improves objective')
        if (family < 3) then
            call close(fit%params(1),scipy_parameter(family),2.0e-4_dp,'SciPy shape fit')
        else
            call maxwell_fit(data,[0.0_dp,0.2_dp],[0.0_dp,3.0_dp],fit,tolerance=1.0e-7_dp)
            call check(fit%success,'Maxwell scale fit')
            call close(fit%params(2),scipy_parameter(3),2.0e-5_dp,'SciPy scale fit')
        end if
        call fitting(empty,family,lower,upper,fit)
        call check(fit%status == fit_status_invalid_input,'empty fit rejected')
        lower = [0.0_dp,0.0_dp,1.3_dp]
        upper = lower
        call fitting(data,family,lower,upper,fit)
        if (family /= 3) call check(.not. fit%success,'invalid fixed shape rejected')
    end do
    single = make_censored_data(left=[0.0_dp])
    do family = 1, 3
        call check(likelihood(single,family,1.0_dp) == -inf,'CDF at lower endpoint')
        other = make_censored_data(right=[-1.0_dp])
        call close(likelihood(other,family,1.0_dp),0.0_dp,0.0_dp,'SF below support')
        other = make_censored_data(interval=reshape([-1.0_dp,2.0_dp],[1,2]))
        empty = make_censored_data(left=[2.0_dp])
        call close(likelihood(other,family,1.0_dp),likelihood(empty,family,1.0_dp), &
            0.0_dp,'support clipping')
    end do
    single = make_censored_data(uncensored=[0.0_dp])
    call check(betaprime_censored_loglikelihood(single,0.4_dp,2.3_dp) == inf, &
        'beta prime singular endpoint')
    call check(betaprime_censored_loglikelihood(single,2.0_dp,2.3_dp) == -inf, &
        'beta prime zero endpoint')
    call close(betaprime_censored_loglikelihood(single,1.0_dp,2.3_dp),log(2.3_dp), &
        0.0_dp,'beta prime finite endpoint')
    call check(invgamma_censored_loglikelihood(single,0.4_dp) == -inf,'inverse gamma endpoint')
    call check(maxwell_censored_loglikelihood(single) == -inf,'Maxwell endpoint')
    empty = make_censored_data()
    call check(ieee_is_nan(betaprime_censored_loglikelihood(empty,1.0_dp,nan)), &
        'invalid second shape')
    call check(ieee_is_nan(betaprime_censored_loglikelihood(empty,huge(1.0_dp), &
        huge(1.0_dp))),'shape sum overflow rejected')
    call close(invgamma_censored_loglikelihood(data,1.0_dp,scale=1.3_dp), &
        invweibull_censored_loglikelihood(data,1.0_dp,scale=1.3_dp),3.0e-13_dp, &
        'inverse Weibull shape-one identity')
    call close(betaprime_censored_loglikelihood(data,1.0_dp,2.3_dp,scale=1.3_dp), &
        lomax_censored_loglikelihood(data,2.3_dp,scale=1.3_dp),3.0e-13_dp,'Lomax identity')
    single = make_censored_data(interval=reshape([1.0e-300_dp,2.0e-300_dp],[1,2]))
    actual = maxwell_censored_loglikelihood(single,scale=1.0e300_dp)
    expected = log(2.0_dp/acos(-1.0_dp))/2.0_dp+log(7.0_dp/3.0_dp) &
        -1800.0_dp*log(10.0_dp)
    call close(actual,expected,2.0e-12_dp,'Maxwell normalized interval underflow')
    single = make_censored_data(right=[1.0e300_dp])
    actual = invgamma_censored_loglikelihood(single,1.0_dp,scale=1.0e-300_dp)
    call close(actual,-600.0_dp*log(10.0_dp),3.0e-13_dp,'inverse gamma far tail')
    actual = betaprime_censored_loglikelihood(single,1.0_dp,2.3_dp,scale=1.0e-300_dp)
    call close(actual,-1380.0_dp*log(10.0_dp),2.0e-12_dp,'beta prime far tail')
    single = make_censored_data(interval=reshape([1.0_dp,nearest(1.0_dp,2.0_dp)],[1,2]))
    do family = 1, 3
        other = make_censored_data(uncensored=[1.0_dp])
        expected = likelihood(other,family,1.7_dp)+log(spacing(1.0_dp))
        call close(likelihood(single,family,1.7_dp),expected,3.0e-13_dp,'adjacent interval')
    end do
    single = make_censored_data(interval=reshape([1.0e-100_dp,1.0e100_dp],[1,2]))
    actual = invgamma_censored_loglikelihood(single,1.0_dp)
    call close(actual/(-1.0e-100_dp),1.0_dp,3.0e-13_dp,'nearly full inverse gamma interval')
    actual = betaprime_censored_loglikelihood(single,1.0_dp,1.0_dp)
    call close(actual/(-2.0e-100_dp),1.0_dp,3.0e-13_dp,'both beta prime omitted tails')
    single = make_censored_data(uncensored=[0.0_dp])
    call betaprime_fit(single,[0.4_dp,2.3_dp,0.0_dp,1.0_dp], &
        [0.4_dp,2.3_dp,0.0_dp,1.0_dp],fit)
    call check(fit%status == fit_status_no_finite_objective,'singular beta prime fixed fit')
    call check(ieee_is_nan(betaprime_censored_loglikelihood(empty,1.0_dp,0.0_dp)), &
        'zero second shape')
    call check(ieee_is_nan(betaprime_censored_loglikelihood(empty,1.0_dp,inf)), &
        'infinite second shape')
    single = make_censored_data(left=[1.25_dp])
    other = make_censored_data(left=[(1.25_dp/1.3_dp)**2/2.0_dp])
    call close(maxwell_censored_loglikelihood(single,scale=1.3_dp), &
        gamma_censored_loglikelihood(other,1.5_dp),3.0e-14_dp,'Maxwell gamma identity')
    single = make_censored_data(uncensored=[1.0e-300_dp])
    call check(betaprime_censored_loglikelihood(single,1.0e308_dp,1.0_dp,scale=1.0e300_dp) &
        == -inf,'unrepresentable beta prime log density')
    actual = betaprime_censored_loglikelihood(single,1.0e305_dp,1.7e308_dp,scale=1.0e300_dp)
    call close(actual/(-1.373112378206389e308_dp),1.0_dp,5.0e-13_dp, &
        'representable beta prime log density after normalization')
    single = make_censored_data(left=[1.0e-300_dp])
    actual = betaprime_censored_loglikelihood(single,1.0e305_dp,1.7e308_dp,scale=1.0e300_dp)
    call close(actual/(-1.373112378206389e308_dp),1.0_dp,5.0e-13_dp, &
        'representable beta prime log CDF after normalization')
    empty = make_censored_data()
    call close(betaprime_censored_loglikelihood(empty,nearest(0.0_dp,1.0_dp),huge(1.0_dp)), &
        0.0_dp,0.0_dp,'finite rounded shape sum with largest second shape')
    call close(betaprime_censored_loglikelihood(empty,huge(1.0_dp),nearest(0.0_dp,1.0_dp)), &
        0.0_dp,0.0_dp,'finite rounded shape sum with largest first shape')
    single = make_censored_data(uncensored=[0.0_dp])
    call close(betaprime_censored_loglikelihood(single,1.0_dp,huge(1.0_dp)), &
        log(huge(1.0_dp)),0.0_dp,'largest second shape endpoint density')
    write(*,'(a,i0,a)') 'Positive transform censored: ',checks,' checks'
    if (failures /= 0) error stop 'Positive transform censoring failures'
contains
    pure function likelihood(data, family, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        integer, intent(in) :: family !! 1 inverse gamma, 2 beta prime, 3 Maxwell
        real(dp), intent(in) :: c !! positive shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive scale, default one
        real(dp) :: value
        select case (family)
        case (1)
            value = invgamma_censored_loglikelihood(data,c,loc,scale)
        case (2)
            value = betaprime_censored_loglikelihood(data,c,2.3_dp,loc,scale)
        case (3)
            value = maxwell_censored_loglikelihood(data,loc,scale)
        end select
    end function likelihood
    subroutine fitting(data, family, lower, upper, result)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        integer, intent(in) :: family !! distribution identifier used by this test
        real(dp), intent(in) :: lower(3) !! shape, location, scale lower bounds
        real(dp), intent(in) :: upper(3) !! shape, location, scale upper bounds
        type(fit_result), intent(out) :: result !! bounded fit diagnostics
        select case (family)
        case (1)
            call invgamma_fit(data,lower,upper,result,max_iter=10000,tolerance=1.0e-7_dp)
        case (2)
            call betaprime_fit(data,[lower(1),2.3_dp,lower(2:3)], &
                [upper(1),2.3_dp,upper(2:3)],result,max_iter=10000,tolerance=1.0e-7_dp)
        case (3)
            call maxwell_fit(data,lower(2:3),upper(2:3),result,max_iter=10000,tolerance=1.0e-7_dp)
        end select
    end subroutine fitting
    subroutine check(condition,label)
        implicit none
        logical, intent(in) :: condition !! true when expected behavior holds
        character(*), intent(in) :: label !! failure diagnostic
        checks = checks+1
        if (condition) return
        failures = failures+1
        write(*,'(a)') 'FAIL: '//label
    end subroutine check
    subroutine close(value,reference,tolerance,label)
        implicit none
        real(dp), intent(in) :: value !! computed finite result
        real(dp), intent(in) :: reference !! expected finite value
        real(dp), intent(in) :: tolerance !! absolute tolerance for this regime
        character(*), intent(in) :: label !! failure diagnostic
        call check(ieee_is_finite(value),label//' finite')
        if (.not. ieee_is_finite(value)) return
        call check(abs(value-reference) <= tolerance,label)
        if (abs(value-reference) > tolerance) write(*,*) value,reference
    end subroutine close
end program test_censored_positive_transform
